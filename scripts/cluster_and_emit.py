#!/usr/bin/env python3
"""Cluster sources by exact sequence identity (SHA256 hashing).

This replaces MMseqs2 because we only want strict 100% identity matches.
Hash-based equivalence is faster, correct by construction, and avoids
MMseqs2 index quirks for nucleotide databases.
"""
from pathlib import Path
import json
import hashlib
import sys
from collections import defaultdict

ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / "ingest" / "work"
GENES = ROOT / "data" / "genes"
GENES.mkdir(parents=True, exist_ok=True)

DATE = "2026-04-27"


def sha256_seq(seq):
    """Hash an upper-cased sequence (case-insensitive equivalence)."""
    return hashlib.sha256(seq.upper().strip().encode()).hexdigest()


def parse_fasta(path):
    """Yield (header, sequence) pairs."""
    if not path.exists():
        return
    cur_h, cur_s = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_h is not None:
                    yield cur_h, "".join(cur_s)
                cur_h = line[1:]
                cur_s = []
            else:
                cur_s.append(line)
        if cur_h is not None:
            yield cur_h, "".join(cur_s)


def parse_member(m):
    if m.startswith("ARO:"):
        return ("CARD", m)
    if "__" in m:
        db, tid = m.split("__", 1)
        return (db, tid)
    return ("Unknown", m)


def emit_yaml(record, path):
    lines = []
    lines.append(f'aro_id: "{record["aro_id"]}"')
    cn = record["canonical_name"].replace('"', "'")
    lines.append(f'canonical_name: "{cn}"')
    gf = record["gene_family"].replace('"', "'")
    lines.append(f'gene_family: "{gf}"')
    lines.append(f'variant_type: {record["variant_type"]}')
    if record.get("mechanism"):
        mc = record["mechanism"].replace('"', "'")
        lines.append(f'mechanism: "{mc}"')
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
        tid = m["target_id"].replace('"', "'")
        ev = m["evidence_value"].replace('"', "'")
        lines.append(f'  - target_db: {m["target_db"]}')
        lines.append(f'    target_id: "{tid}"')
        lines.append(f'    relation: "{m["relation"]}"')
        lines.append(f'    level: {m["level"]}')
        lines.append(f'    evidence:')
        lines.append(f'      type: {m["evidence_type"]}')
        lines.append(f'      value: "{ev}"')
        lines.append(f'    quality: {m["quality"]}')
        lines.append(f'    date: "{m["date"]}"')
    path.write_text("\n".join(lines) + "\n")


def cluster_by_hash(*fastas):
    """Build {sha256: set(headers)} across all input FASTAs."""
    clusters = defaultdict(set)
    total = 0
    for fa in fastas:
        for header, seq in parse_fasta(fa):
            if not seq:
                continue
            h = sha256_seq(seq)
            clusters[h].add(header)
            total += 1
    print(f"  {total} sequences -> {len(clusters)} unique sequence hashes")
    return clusters


def main():
    card = {r["aro_id"]: r for r in json.loads((WORK / "card_anchor.json").read_text())}
    print(f"Loaded {len(card)} CARD anchor records")

    print("\nHash-clustering proteins...")
    prot_clusters = cluster_by_hash(
        WORK / "card_proteins.fa",
        WORK / "all_proteins.fa",
    )

    print("\nHash-clustering nucleotides...")
    nucl_clusters = cluster_by_hash(
        WORK / "card_nucleotide.fa",
        WORK / "all_nucleotide.fa",
    )

    def count_multi_db(clusters):
        n = 0
        for members in clusters.values():
            dbs = {m.split("__")[0] if "__" in m else "CARD" for m in members}
            if len(dbs) > 1:
                n += 1
        return n

    print(f"\nProtein clusters spanning multiple DBs: {count_multi_db(prot_clusters)}")
    print(f"Nucleotide clusters spanning multiple DBs: {count_multi_db(nucl_clusters)}")

    by_aro = defaultdict(lambda: {"protein": set(), "nucleotide": set()})
    for level, clusters in [("protein", prot_clusters), ("nucleotide", nucl_clusters)]:
        for members in clusters.values():
            aro_ids = [m for m in members if m.startswith("ARO:")]
            if not aro_ids:
                continue
            anchor = aro_ids[0]
            for m in members:
                if m == anchor:
                    continue
                by_aro[anchor][level].add(m)

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
                    "evidence_value": f"SHA256 hash match at {level} level",
                    "quality": "bronze",
                    "date": DATE,
                })

        record = {
            "aro_id": aro,
            "canonical_name": meta["canonical_name"] or meta["card_short_name"] or aro,
            "gene_family": meta["gene_family"] or "unspecified",
            "variant_type": meta["variant_type"],
            "mechanism": meta["mechanism"].split(";")[0].strip().replace(" ", "_") if meta["mechanism"] else None,
            "ncbi_protein": meta.get("ncbi_protein"),
            "ncbi_nuccore": meta.get("ncbi_nuccore"),
            "length_aa": len(meta["protein_seq"]) if meta.get("protein_seq") else None,
            "length_nt": len(meta["nucleotide_seq"]) if meta.get("nucleotide_seq") else None,
            "mappings": mappings,
        }

        fname = f"ARO_{aro.split(':')[1]}.yaml"
        emit_yaml(record, GENES / fname)
        written += 1

    print(f"\nWrote {written} bronze YAML records to {GENES}")


if __name__ == "__main__":
    main()
