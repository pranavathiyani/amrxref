#!/usr/bin/env python3
"""
Parse CARD as the canonical anchor.

Reads:
  ingest/sources/card/aro_index.tsv
  ingest/sources/card/protein_fasta_protein_homolog_model.fasta
  ingest/sources/card/nucleotide_fasta_protein_homolog_model.fasta
  ingest/sources/card/nucleotide_fasta_rRNA_gene_variant_model.fasta

Writes:
  ingest/work/card_anchor.json   (one record per ARO)
  ingest/work/card_proteins.fa   (relabeled FASTA: >ARO:NNN)
  ingest/work/card_nucleotide.fa (relabeled FASTA: >ARO:NNN)
"""
from pathlib import Path
import json
import re
import csv
import sys

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "ingest" / "sources" / "card"
WORK = ROOT / "ingest" / "work"
WORK.mkdir(parents=True, exist_ok=True)


def parse_aro_index():
    """aro_index.tsv → {aro_id: {canonical_name, gene_family, drug_class, mechanism, ...}}"""
    out = {}
    with (SRC / "aro_index.tsv").open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            aro = row["ARO Accession"].strip()
            if not aro.startswith("ARO:"):
                continue
            out[aro] = {
                "aro_id": aro,
                "canonical_name": row["ARO Name"].strip(),
                "card_short_name": row["CARD Short Name"].strip(),
                "model_name": row["Model Name"].strip(),
                "gene_family": row["AMR Gene Family"].strip(),
                "drug_classes_raw": row["Drug Class"].strip(),
                "mechanism": row["Resistance Mechanism"].strip(),
                "ncbi_protein": row["Protein Accession"].strip() or None,
                "ncbi_nuccore": row["DNA Accession"].strip() or None,
            }
    return out


def parse_card_fasta(path):
    """CARD FASTA headers look like:  >gb|ACT97415.1|ARO:3002999|CblA-1 [...]
       Returns {aro_id: sequence}.
    """
    out = {}
    current_aro = None
    current_seq = []
    pattern = re.compile(r"\|(?:ARO:|ARO_)(\d+)\|")
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_aro and current_seq:
                    out[current_aro] = "".join(current_seq)
                m = pattern.search(line)
                current_aro = f"ARO:{m.group(1)}" if m else None
                current_seq = []
            else:
                if current_aro:
                    current_seq.append(line)
        if current_aro and current_seq:
            out[current_aro] = "".join(current_seq)
    return out


def main():
    print("Parsing CARD aro_index...")
    aro_meta = parse_aro_index()
    print(f"  {len(aro_meta)} ARO entries")

    print("Parsing protein homolog FASTA...")
    proteins = parse_card_fasta(SRC / "protein_fasta_protein_homolog_model.fasta")
    print(f"  {len(proteins)} protein sequences")

    print("Parsing nucleotide homolog FASTA...")
    nucleotides = parse_card_fasta(SRC / "nucleotide_fasta_protein_homolog_model.fasta")
    print(f"  {len(nucleotides)} nucleotide sequences")

    print("Parsing rRNA variant FASTA...")
    rrna = parse_card_fasta(SRC / "nucleotide_fasta_rRNA_gene_variant_model.fasta")
    print(f"  {len(rrna)} rRNA sequences")

    # Build records
    records = []
    for aro, meta in aro_meta.items():
        rec = dict(meta)
        rec["protein_seq"] = proteins.get(aro)
        rec["nucleotide_seq"] = nucleotides.get(aro) or rrna.get(aro)
        if aro in rrna:
            rec["variant_type"] = "rrna"
        elif rec["protein_seq"]:
            rec["variant_type"] = "protein_coding"
        else:
            rec["variant_type"] = "regulatory_region"  # fallback; refined later
        records.append(rec)

    out_path = WORK / "card_anchor.json"
    out_path.write_text(json.dumps(records, indent=2))
    print(f"Wrote {len(records)} records → {out_path}")

    # Relabeled FASTAs (header = ARO ID only) for clean MMseqs2 input
    with (WORK / "card_proteins.fa").open("w") as f:
        for r in records:
            if r.get("protein_seq"):
                f.write(f">{r['aro_id']}\n{r['protein_seq']}\n")

    with (WORK / "card_nucleotide.fa").open("w") as f:
        for r in records:
            if r.get("nucleotide_seq"):
                f.write(f">{r['aro_id']}\n{r['nucleotide_seq']}\n")

    print("Wrote relabeled FASTAs to ingest/work/")


if __name__ == "__main__":
    main()
