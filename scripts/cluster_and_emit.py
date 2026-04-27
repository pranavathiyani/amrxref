#!/usr/bin/env python3
"""
Cluster all sources at 100% identity and emit one bronze YAML per gene.

Reads:
  ingest/work/card_anchor.json
  ingest/work/card_proteins.fa
  ingest/work/card_nucleotide.fa
  ingest/work/all_proteins.fa
  ingest/work/all_nucleotide.fa

Calls:
  mmseqs easy-cluster (must be on PATH)

Writes:
  data/genes/ARO_*.yaml         (bronze records anchored to ARO)
  data/genes/AMRXREF_*.yaml     (bronze records without ARO anchor)
"""
from pathlib import Path
import json
import subprocess
import shutil
import sys
import re
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / "ingest" / "work"
GENES = ROOT / "data" / "genes"
GENES.mkdir(parents=True, exist_ok=True)

DATE = "2026-04-27"


def cat_fastas(out_path, *inputs):
    with open(out_path, "w") as out:
        for inp in inputs:
            if inp.exists():
                out.write(inp.read_text())


def run_mmseqs(input_fa, prefix):
    """Cluster at 100% identity. Returns path to TSV mapping rep -> member."""
    tmp = WORK / f"{prefix}_tmp"
    tmp.mkdir(exist_ok=True)
    out_prefix = WORK / prefix
    cmd = [
        "mmseqs", "easy-cluster",
        str(input_fa), str(out_prefix), str(tmp),
        "--min-seq-id", "1.0",
        "-c", "1.0",
        "--cov-mode", "0",
        "--threads", "4",
    ]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    return Path(f"{out_prefix}_cluster.tsv")


def parse_clusters(tsv):
    """MMseqs2 _cluster.tsv: <representative>\\t<member>"""
    clusters = defaultdict(set)
    with open(tsv) as f:
        for line in f:
            rep, mem = line.rstrip().split("\t")
            clusters[rep].add(mem)
    return clusters


def parse_member(m):
    """'ARO:3000873' or 'AMRFinderPlus__blaTEM-1' or 'ResFinder__blaTEM-1_1_AY...' """
    if m.startswith("ARO:"):
        return ("CARD", m)
    if "__" in m:
        db, tid = m.split("__", 1)
        return (db, tid)
    return ("Unknown", m)


def emit_yaml(record, path):
    """Write a YAML record manually (no yaml dependency required for output)."""
    lines = []
    lines.append(f'aro_id: "{record["aro_id"]}"')
    lines.append(f'canonical_name: {record["canonical_name"]}')
    lines.append(f'gene_family: {record["gene_family"]}')
    if record.get("preferred_label"):
        lines.append(f'preferred_label: "{record["preferred_label"]}"')
    lines.append(f'variant_type: {record["variant_type"]}')

    if record.get("drug_class_label"):
        lines.append("drug_classes:")
        lines.append(f'  - aro: "ARO:0000000"   # placeholder, refine via ARO drug branch')
        lines.append(f'    label: {record["drug_class_label"]}')

    if record.get("mechanism"):
        lines.append(f'mechanism: {record["mechanism"]}')
    if record.get("ncbi_protein"):
        lines.append("canonical_protein:")
        lines.append(f'  ncbi_protein: {record["ncbi_protein"]}')
        if record.get("length_aa"):
            lines.append(f'  length_aa: {record["length_aa"]}')
    if record.get("ncbi_nuccore"):
        lines.append("canonical_nucleotide:")
        lines.append(f'  ncbi_nuccore: {record["ncbi_nuccore"]}')
        if record.get("length_nt"):
            lines.append(f'  length_nt: {record["length_nt"]}')

    lines.append("")
    lines.append("mappings:")
    for m in record["mappings"]:
        lines.append(f'  - target_db: {m["target_db"]}')
        lines.append(f'    target_id: "{m["target_id"]}"')
        lines.append(f'    relation: "{m["relation"]}"')
        lines.append(f'    level: {m["level"]}')
        lines.append(f'    evidence:')
        lines.append(f'      type: {m["evidence_type"]}')
        lines.append(f'      value: "{m["evidence_value"]}"')
        lines.append(f'    quality: {m["quality"]}')
        lines.append(f'    date: "{m["date"]}"')
    path.write_text("\n".join(lines) + "\n")


def safe_filename(s):
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)


def main():
    if not shutil.which("mmseqs"):
        print("ERROR: mmseqs not on PATH. Install with: conda install -c bioconda mmseqs2", file=sys.stderr)
        sys.exit(1)

    # Load CARD anchor lookup
    card = {r["aro_id"]: r for r in json.loads((WORK / "card_anchor.json").read_text())}
    print(f"Loaded {len(card)} CARD anchor records")

    # Build merged FASTAs
    print("Merging FASTAs...")
    cat_fastas(WORK / "merged_proteins.fa",
               WORK / "card_proteins.fa", WORK / "all_proteins.fa")
    cat_fastas(WORK / "merged_nucleotide.fa",
               WORK / "card_nucleotide.fa", WORK / "all_nucleotide.fa")

    # Run MMseqs2 clustering
    print("Clustering proteins at 100% identity...")
    prot_tsv = run_mmseqs(WORK / "merged_proteins.fa", "cluster_prot")
    print("Clustering nucleotides at 100% identity...")
    nucl_tsv = run_mmseqs(WORK / "merged_nucleotide.fa", "cluster_nucl")

    prot_clusters = parse_clusters(prot_tsv)
    nucl_clusters = parse_clusters(nucl_tsv)
    print(f"Protein clusters: {len(prot_clusters)}")
    print(f"Nucleotide clusters: {len(nucl_clusters)}")

    # Build a unified ARO → mappings dict
    by_aro = defaultdict(lambda: {"protein": set(), "nucleotide": set()})

    for level, clusters in [("protein", prot_clusters), ("nucleotide", nucl_clusters)]:
        for rep, members in clusters.items():
            # find any ARO in this cluster
            aro_ids = [m for m in members if m.startswith("ARO:")]
            if not aro_ids:
                continue  # skip clusters with no CARD anchor for now
            anchor = aro_ids[0]
            for m in members:
                if m == anchor:
                    continue
                by_aro[anchor][level].add(m)

    # Emit one YAML per ARO with cross-references
    written = 0
    for aro, levels in by_aro.items():
        meta = card.get(aro)
        if not meta:
            continue
        mappings = [{
            "target_db": "CARD",
            "target_id": aro,
            "relation": "skos:exactMatch",
            "level": "protein" if meta["variant_type"] == "protein_coding" else "nucleotide",
            "evidence_type": "ontology_anchor",
            "evidence_value": "Canonical ARO term",
            "quality": "gold",
            "date": DATE,
        }]
        for level in ("protein", "nucleotide"):
            for m in sorted(levels[level]):
                db, tid = parse_member(m)
                if db == "CARD":
                    continue
                mappings.append({
                    "target_db": db,
                    "target_id": tid,
                    "relation": "skos:exactMatch",
                    "level": level,
                    "evidence_type": "sequence_identity_100",
                    "evidence_value": f"100% identity at {level} level via MMseqs2",
                    "quality": "bronze",
                    "date": DATE,
                })

        record = {
            "aro_id": aro,
            "canonical_name": meta["canonical_name"] or meta["card_short_name"] or aro,
            "gene_family": meta["gene_family"] or "unspecified",
            "variant_type": meta["variant_type"],
            "mechanism": meta["mechanism"].split(";")[0].strip().replace(" ", "_") if meta["mechanism"] else None,
            "drug_class_label": (meta["drug_classes_raw"].split(";")[0].strip()
                                 if meta["drug_classes_raw"] else None),
            "ncbi_protein": meta.get("ncbi_protein"),
            "ncbi_nuccore": meta.get("ncbi_nuccore"),
            "length_aa": len(meta["protein_seq"]) if meta.get("protein_seq") else None,
            "length_nt": len(meta["nucleotide_seq"]) if meta.get("nucleotide_seq") else None,
            "mappings": mappings,
        }

        fname = f"ARO_{aro.split(':')[1]}.yaml"
        emit_yaml(record, GENES / fname)
        written += 1

    print(f"Wrote {written} bronze YAML records to {GENES}")
    print("Note: clusters without CARD anchor were skipped in v0.1.")


if __name__ == "__main__":
    main()
