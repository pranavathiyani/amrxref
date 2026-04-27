#!/usr/bin/env python3
"""Generate coverage statistics from derived/crosswalk.json.

Writes docs/stats.json with computed numbers for the coverage page.
"""
import json
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).resolve().parent.parent
DERIVED = ROOT / "derived" / "crosswalk.json"
OUT = ROOT / "docs" / "stats.json"


def overall_quality(rec):
    qs = [m["quality"] for m in rec["mappings"]]
    if all(q == "gold" for q in qs):
        return "gold"
    if "bronze" in qs:
        return "bronze"
    return "silver"


def main():
    data = json.loads(DERIVED.read_text())
    total = len(data)

    # Coverage combos (which DB combinations exist)
    combos = Counter()
    for r in data:
        dbs = sorted(set(m["target_db"] for m in r["mappings"]))
        combos[" + ".join(dbs)] += 1

    # Variant type breakdown
    variants = Counter(r["variant_type"] for r in data)

    # Quality distribution
    qualities = Counter(overall_quality(r) for r in data)

    # Drug class distribution (CARD's gene_family is a proxy)
    families = Counter(r.get("gene_family", "unspecified") for r in data)

    # Per-database hit counts
    db_hits = Counter()
    for r in data:
        for m in r["mappings"]:
            db_hits[m["target_db"]] += 1

    # Records with each individual DB
    has_db = Counter()
    for r in data:
        seen = set(m["target_db"] for m in r["mappings"])
        for db in seen:
            has_db[db] += 1

    out = {
        "total_records": total,
        "coverage_combos": [
            {"combo": k, "count": v}
            for k, v in combos.most_common()
        ],
        "variant_types": dict(variants),
        "qualities": dict(qualities),
        "top_gene_families": [
            {"family": k, "count": v}
            for k, v in families.most_common(20)
        ],
        "total_mappings_per_db": dict(db_hits),
        "records_with_db": dict(has_db),
    }

    OUT.write_text(json.dumps(out, indent=2))
    print(f"Wrote {OUT}")

    # Console summary
    print(f"\nTotal records: {total}")
    print(f"Variant types: {dict(variants)}")
    print(f"Quality tiers: {dict(qualities)}")
    print(f"\nTop coverage combos:")
    for c in combos.most_common(5):
        print(f"  {c[1]:5d}  {c[0]}")
    print(f"\nRecords containing each DB:")
    for k, v in has_db.most_common():
        pct = 100 * v / total
        print(f"  {k:15s}  {v:5d}  ({pct:.1f}%)")


if __name__ == "__main__":
    main()
