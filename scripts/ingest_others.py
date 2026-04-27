#!/usr/bin/env python3
"""
Parse AMRFinderPlus, ResFinder, MEGARes into uniform JSON.

Reads:
  ingest/sources/amrfinderplus/AMRProt.fa
  ingest/sources/amrfinderplus/AMR_CDS.fa
  ingest/sources/amrfinderplus/ReferenceGeneCatalog.txt
  ingest/sources/resfinder_db/*.fsa
  ingest/sources/megares_v3.fasta
  ingest/sources/megares_annotations_v3.csv

Writes:
  ingest/work/amrfinderplus.json
  ingest/work/resfinder.json
  ingest/work/megares.json
  ingest/work/all_proteins.fa
  ingest/work/all_nucleotide.fa
"""
from pathlib import Path
import json
import csv
import re
import sys

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "ingest" / "sources"
WORK = ROOT / "ingest" / "work"
WORK.mkdir(parents=True, exist_ok=True)


def parse_fasta(path):
    """Generic FASTA → list of (header_line_minus_gt, sequence)."""
    out = []
    cur_h, cur_s = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_h is not None:
                    out.append((cur_h, "".join(cur_s)))
                cur_h = line[1:]
                cur_s = []
            else:
                cur_s.append(line)
        if cur_h is not None:
            out.append((cur_h, "".join(cur_s)))
    return out


def ingest_amrfinderplus():
    """AMRFinderPlus uses pipe-delimited headers and a separate metadata TSV.
       AMRProt header format:
         >NCBI_acc|N|N|allele|gene_family||N|...|product_name
    """
    catalog = {}
    with (SRC / "amrfinderplus" / "ReferenceGeneCatalog.txt").open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            allele = row.get("allele", "").strip()
            if not allele:
                continue
            catalog[allele] = {
                "allele": allele,
                "gene_family": row.get("gene_family", "").strip(),
                "product_name": row.get("product_name", "").strip(),
                "class": row.get("class", "").strip(),
                "subclass": row.get("subclass", "").strip(),
                "type": row.get("type", "").strip(),
                "refseq_protein": row.get("refseq_protein_accession", "").strip(),
                "refseq_nucleotide": row.get("refseq_nucleotide_accession", "").strip(),
            }

    proteins = parse_fasta(SRC / "amrfinderplus" / "AMRProt.fa")
    nucleotides = parse_fasta(SRC / "amrfinderplus" / "AMR_CDS.fa")

    out = []
    for header, seq in proteins:
        parts = header.split("|")
        if len(parts) < 9:
            continue
        ncbi_acc = parts[0]
        allele = parts[3]
        meta = catalog.get(allele, {})
        out.append({
            "db": "AMRFinderPlus",
            "target_id": allele,
            "ncbi_protein": ncbi_acc,
            "level": "protein",
            "sequence": seq,
            "gene_family": meta.get("gene_family", ""),
            "drug_class": meta.get("class", ""),
            "subclass": meta.get("subclass", ""),
            "product_name": meta.get("product_name", ""),
        })

    for header, seq in nucleotides:
        parts = header.split("|")
        if len(parts) < 9:
            continue
        ncbi_acc = parts[0]
        allele = parts[3]
        meta = catalog.get(allele, {})
        out.append({
            "db": "AMRFinderPlus",
            "target_id": allele,
            "ncbi_nuccore": ncbi_acc,
            "level": "nucleotide",
            "sequence": seq,
            "gene_family": meta.get("gene_family", ""),
            "drug_class": meta.get("class", ""),
        })
    return out


def ingest_resfinder():
    """ResFinder FASTA headers: >gene_variant_accession (e.g. blaNDM-19_1_MF370080).
       File name = drug class.
    """
    out = []
    rfdir = SRC / "resfinder_db"
    for fsa in sorted(rfdir.glob("*.fsa")):
        if fsa.name == "all.fsa":
            continue
        drug_class = fsa.stem
        for header, seq in parse_fasta(fsa):
            # split: gene_variant_accession  →  parts may be e.g. ['blaNDM-19', '1', 'MF370080']
            parts = header.split("_")
            if len(parts) >= 3:
                accession = parts[-1]
                gene_id = "_".join(parts[:-1])  # gene + variant number
            else:
                accession = ""
                gene_id = header
            out.append({
                "db": "ResFinder",
                "target_id": header,        # full header is the unique ID
                "gene_id": gene_id,
                "ncbi_nuccore": accession,
                "level": "nucleotide",
                "sequence": seq,
                "drug_class": drug_class,
            })
    return out


def ingest_megares():
    """MEGARes header: >MEG_N|type|class|mechanism|group|requires_snp"""
    annotations = {}
    with (SRC / "megares_annotations_v3.csv").open() as f:
        r = csv.DictReader(f)
        for row in r:
            annotations[row["header"]] = row

    out = []
    for header, seq in parse_fasta(SRC / "megares_v3.fasta"):
        parts = header.split("|")
        meg_id = parts[0]
        meta = annotations.get(header, {})
        out.append({
            "db": "MEGARes",
            "target_id": meg_id,
            "full_header": header,
            "level": "nucleotide",
            "sequence": seq,
            "type": parts[1] if len(parts) > 1 else "",
            "drug_class": meta.get("class", parts[2] if len(parts) > 2 else ""),
            "mechanism": meta.get("mechanism", ""),
            "group": meta.get("group", parts[4] if len(parts) > 4 else ""),
            "requires_snp": "RequiresSNPConfirmation" in header,
        })
    return out


def main():
    print("Ingesting AMRFinderPlus...")
    afp = ingest_amrfinderplus()
    (WORK / "amrfinderplus.json").write_text(json.dumps(afp, indent=2))
    print(f"  {len(afp)} entries")

    print("Ingesting ResFinder...")
    rf = ingest_resfinder()
    (WORK / "resfinder.json").write_text(json.dumps(rf, indent=2))
    print(f"  {len(rf)} entries")

    print("Ingesting MEGARes...")
    mg = ingest_megares()
    (WORK / "megares.json").write_text(json.dumps(mg, indent=2))
    print(f"  {len(mg)} entries")

    # Combined FASTAs for MMseqs2 (CARD records added by ingest_card.py separately)
    print("Writing combined FASTAs...")
    with (WORK / "all_proteins.fa").open("w") as f:
        for e in afp:
            if e["level"] == "protein" and e["sequence"]:
                f.write(f">AMRFinderPlus__{e['target_id']}\n{e['sequence']}\n")

    with (WORK / "all_nucleotide.fa").open("w") as f:
        for e in afp:
            if e["level"] == "nucleotide" and e["sequence"]:
                f.write(f">AMRFinderPlus__{e['target_id']}\n{e['sequence']}\n")
        for e in rf:
            f.write(f">ResFinder__{e['target_id']}\n{e['sequence']}\n")
        for e in mg:
            f.write(f">MEGARes__{e['target_id']}\n{e['sequence']}\n")

    print("Done.")


if __name__ == "__main__":
    main()
