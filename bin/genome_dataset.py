#!/usr/bin/env python3
"""Retrieve and audit explicitly accessioned NCBI genome assemblies."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import tempfile
from urllib.error import URLError
from urllib.request import Request, urlopen
import zipfile

FIELDS = ["cohort_id", "isolate_id", "pdt_accession", "biosample_accession",
          "assembly_accession", "organism", "assembly_source", "ast_available",
          "selection_reason"]
ACCESSION = re.compile(r"^(?:GCF|GCA)_\d{9}\.\d+$")
BIOSAMPLE = re.compile(r"^SAM[END]\d+$")
PDT = re.compile(r"^PDT\d+\.\d+$")
MAX_ARCHIVE = 100 * 1024 * 1024
API = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{}"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if reader.fieldnames != FIELDS:
            raise ValueError(f"manifest columns must be exactly: {' '.join(FIELDS)}")
        rows, seen_assemblies, seen_isolates = [], set(), set()
        for line, row in enumerate(reader, 2):
            if None in row or any(not (v or "").strip() for v in row.values()):
                raise ValueError(f"manifest line {line}: missing or extra field")
            row = {key: value.strip() for key, value in row.items()}
            if not ACCESSION.fullmatch(row["assembly_accession"]):
                raise ValueError(f"manifest line {line}: invalid Assembly accession")
            if not BIOSAMPLE.fullmatch(row["biosample_accession"]):
                raise ValueError(f"manifest line {line}: invalid BioSample accession")
            if not PDT.fullmatch(row["pdt_accession"]):
                raise ValueError(f"manifest line {line}: invalid PDT accession")
            if row["ast_available"] != "true":
                raise ValueError(f"manifest line {line}: cohort requires a linked AST observation")
            if row["organism"] != "Escherichia coli":
                raise ValueError(f"manifest line {line}: organism must be exact Escherichia coli")
            if row["assembly_accession"] in seen_assemblies:
                raise ValueError(f"duplicate assembly accession: {row['assembly_accession']}")
            if row["isolate_id"] in seen_isolates:
                raise ValueError(f"duplicate isolate mapping: {row['isolate_id']}")
            seen_assemblies.add(row["assembly_accession"])
            seen_isolates.add(row["isolate_id"])
            rows.append(row)
    if not rows:
        raise ValueError("manifest contains no validation isolates")
    return rows


def request_bytes(url, limit):
    if not url.startswith("https://"):
        raise ValueError("only HTTPS sources are accepted")
    req = Request(url, headers={"User-Agent": "AMRScan-genome-dataset/1.0"})
    with urlopen(req, timeout=120) as response:
        blocks, size = [], 0
        while True:
            block = response.read(1024 * 1024)
            if not block:
                break
            size += len(block)
            if size > limit:
                raise ValueError(f"download exceeds {limit} bytes")
            blocks.append(block)
        return b"".join(blocks), response.geturl()


def parse_fasta(data):
    text = data.decode("ascii")
    records, bases, header = 0, 0, False
    for line_number, line in enumerate(text.splitlines(), 1):
        if line.startswith(">"):
            if not line[1:].strip():
                raise ValueError(f"empty FASTA header at line {line_number}")
            records += 1
            header = True
        elif line.strip():
            if not header or not re.fullmatch(r"[A-Za-z*.-]+", line.strip()):
                raise ValueError(f"invalid FASTA sequence at line {line_number}")
            bases += len(line.strip())
    if records == 0 or bases == 0:
        raise ValueError("FASTA must contain a header and sequence")
    return records, bases


def package_fasta(archive, assembly):
    source = io.BytesIO(archive) if isinstance(archive, bytes) else archive
    with zipfile.ZipFile(source) as package:
        names = [n for n in package.namelist()
                 if n.startswith(f"ncbi_dataset/data/{assembly}/") and n.endswith("_genomic.fna")]
        if len(names) != 1:
            raise ValueError(f"NCBI package must contain exactly one FASTA for {assembly}; found {len(names)}")
        return package.read(names[0]), names[0]


def download(manifest, outdir):
    rows = read_manifest(manifest)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    provenance = outdir / "provenance.json"
    if provenance.exists():
        raise ValueError(f"refusing to replace existing provenance: {provenance}")
    entries = []
    for row in rows:
        assembly = row["assembly_accession"]
        fasta_path = outdir / f"{row['isolate_id']}_{assembly}.fna"
        if fasta_path.exists():
            raise ValueError(f"refusing to replace existing assembly file: {fasta_path}")
        url = API.format(assembly) + "/download?include_annotation_type=GENOME_FASTA,SEQUENCE_REPORT&filename=assembly.zip"
        try:
            archive, final_url = request_bytes(url, MAX_ARCHIVE)
            with tempfile.TemporaryDirectory(dir=outdir) as tmp:
                zip_path = Path(tmp) / "assembly.zip"
                zip_path.write_bytes(archive)
                fasta, member = package_fasta(zip_path, assembly)
                records, bases = parse_fasta(fasta)
                # The same official package carries the assembly report used to verify identity.
                with zipfile.ZipFile(zip_path) as package:
                    reports = [n for n in package.namelist() if n.endswith("assembly_data_report.jsonl")]
                    if len(reports) != 1:
                        raise ValueError("NCBI package lacks a unique assembly data report")
                    report_lines = package.read(reports[0]).decode("utf-8").splitlines()
                report = next((json.loads(line) for line in report_lines if line.strip()), {})
                assembly_info = report.get("assemblyInfo", {})
                biosample = report.get("biologicalMaterial", {}).get("biosample", {}).get("accession")
                organism = report.get("organism", {}).get("organismName")
                reported_accession = assembly_info.get("assemblyAccession")
                if reported_accession != assembly or biosample != row["biosample_accession"]:
                    raise ValueError(f"NCBI identity mismatch for {assembly}: accession={reported_accession}, BioSample={biosample}")
                if organism != row["organism"]:
                    raise ValueError(f"NCBI organism mismatch for {assembly}: {organism}")
                # Exclusive create plus atomic rename; never overwrite another file.
                with tempfile.NamedTemporaryFile(dir=outdir, prefix=".assembly-", delete=False) as tmp_fasta:
                    temp_path = Path(tmp_fasta.name)
                    tmp_fasta.write(fasta)
                if fasta_path.exists():
                    raise ValueError(f"assembly file appeared during download: {fasta_path}")
                temp_path.replace(fasta_path)
            entries.append({**row, "file": fasta_path.name, "sha256": sha256(fasta_path),
                            "source_url": final_url, "retrieved_at": datetime.now(timezone.utc).isoformat(),
                            "package_member": member, "sequence_records": records, "sequence_bases": bases})
        except Exception:
            if fasta_path.exists():
                fasta_path.unlink()
            raise
    payload = {"schema_version": "1.0", "manifest_sha256": sha256(manifest), "assemblies": entries}
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=outdir, prefix=".provenance-", delete=False) as stream:
        json.dump(payload, stream, indent=2, sort_keys=True)
        stream.write("\n")
        temp_manifest = Path(stream.name)
    if provenance.exists():
        temp_manifest.unlink()
        raise ValueError(f"refusing to replace existing provenance: {provenance}")
    temp_manifest.replace(provenance)


def validate(manifest, genomes):
    rows = read_manifest(manifest)
    genomes = Path(genomes)
    payload = json.loads((genomes / "provenance.json").read_text(encoding="utf-8"))
    if payload.get("manifest_sha256") != sha256(manifest):
        raise ValueError("cohort manifest hash differs from retrieval provenance")
    by_assembly = {item["assembly_accession"]: item for item in payload["assemblies"]}
    if len(by_assembly) != len(payload["assemblies"]):
        raise ValueError("duplicate assembly entries in provenance")
    if set(by_assembly) != {row["assembly_accession"] for row in rows}:
        raise ValueError("manifest and provenance assembly sets differ")
    for row in rows:
        item = by_assembly[row["assembly_accession"]]
        for field in ("biosample_accession", "organism", "isolate_id"):
            if item[field] != row[field]:
                raise ValueError(f"{field} mismatch for {row['assembly_accession']}")
        path = genomes / item["file"]
        if path.parent != genomes or not path.is_file() or sha256(path) != item["sha256"]:
            raise ValueError(f"missing file or SHA-256 mismatch for {row['assembly_accession']}")
        parse_fasta(path.read_bytes())
    return len(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    cmd = commands.add_parser("download")
    cmd.add_argument("--manifest", required=True, type=Path)
    cmd.add_argument("--outdir", required=True, type=Path)
    cmd = commands.add_parser("validate")
    cmd.add_argument("--manifest", required=True, type=Path)
    cmd.add_argument("--genomes", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.action == "download":
            download(args.manifest, args.outdir)
            print("download complete")
        else:
            print(f"validated assemblies: {validate(args.manifest, args.genomes)}")
    except (ValueError, OSError, KeyError, json.JSONDecodeError, UnicodeError, zipfile.BadZipFile, URLError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
