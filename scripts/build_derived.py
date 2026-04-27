#!/usr/bin/env python3
"""Build derived artifacts from data/genes/*.yaml.

Outputs (in derived/):
  crosswalk.tsv     flat table for spreadsheet use
  crosswalk.json    array of full records
  crosswalk.jsonld  JSON-LD with SKOS context
  crosswalk.ttl     SKOS Turtle for RDF tooling
"""
from pathlib import Path
import json
import csv
import sys
import yaml
from jsonschema import validate, ValidationError

ROOT = Path(__file__).resolve().parent.parent
GENES_DIR = ROOT / "data" / "genes"
SCHEMA_PATH = ROOT / "data" / "schema" / "gene.schema.json"
OUT = ROOT / "derived"
W3ID = "https://w3id.org/amrxref"


# YAML loader that does NOT auto-convert ISO dates to date objects.
# This keeps "2026-04-27" as a string, which JSON Schema and JSON output expect.
class _NoDateLoader(yaml.SafeLoader):
    pass

_NoDateLoader.yaml_implicit_resolvers = {
    k: [(tag, regexp) for tag, regexp in v if tag != "tag:yaml.org,2002:timestamp"]
    for k, v in _NoDateLoader.yaml_implicit_resolvers.items()
}


def load_records():
    schema = json.loads(SCHEMA_PATH.read_text())
    records = []
    for yml in sorted(GENES_DIR.glob("*.yaml")):
        rec = yaml.load(yml.read_text(), Loader=_NoDateLoader)
        try:
            validate(rec, schema)
        except ValidationError as e:
            print(f"FAIL {yml.name}: {e.message}", file=sys.stderr)
            sys.exit(1)
        records.append(rec)
    return records


def overall_quality(rec):
    qs = [m["quality"] for m in rec["mappings"]]
    if all(q == "gold" for q in qs):
        return "gold"
    if "bronze" in qs:
        return "bronze"
    return "silver"


def write_tsv(records, path):
    fields = [
        "aro_id", "canonical_name", "gene_family", "variant_type",
        "amrfinderplus", "resfinder", "megares",
        "ncbi_protein", "uniprot", "drug_classes", "quality",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in records:
            row = {
                "aro_id": r["aro_id"],
                "canonical_name": r["canonical_name"],
                "gene_family": r["gene_family"],
                "variant_type": r["variant_type"],
                "drug_classes": ";".join(d["label"] for d in r.get("drug_classes", [])),
                "quality": overall_quality(r),
            }
            for m in r["mappings"]:
                db = m["target_db"].lower()
                if db == "amrfinderplus":
                    row["amrfinderplus"] = m["target_id"]
                elif db == "resfinder":
                    row["resfinder"] = m["target_id"]
                elif db == "megares":
                    row["megares"] = m["target_id"]
                elif db == "ncbi":
                    row["ncbi_protein"] = m["target_id"]
                elif db == "uniprot":
                    row["uniprot"] = m["target_id"]
            w.writerow(row)


def write_json(records, path):
    path.write_text(json.dumps(records, indent=2, default=str))


def write_jsonld(records, path):
    doc = {
        "@context": {
            "skos": "http://www.w3.org/2004/02/skos/core#",
            "aro": "http://purl.obolibrary.org/obo/ARO_",
            "amrxref": f"{W3ID}/",
            "canonical_name": "skos:prefLabel",
            "mappings": "skos:mappingRelation",
        },
        "@graph": records,
    }
    path.write_text(json.dumps(doc, indent=2, default=str))


def write_turtle(records, path):
    lines = [
        "@prefix skos: <http://www.w3.org/2004/02/skos/core#> .",
        "@prefix aro: <http://purl.obolibrary.org/obo/ARO_> .",
        f"@prefix amrxref: <{W3ID}/> .",
        "",
    ]
    for r in records:
        aro_local = r["aro_id"].split(":")[1]
        subj = f"amrxref:{r['canonical_name']}"
        lines.append(f"{subj} a skos:Concept ;")
        lines.append(f'    skos:prefLabel "{r["canonical_name"]}" ;')
        lines.append(f"    skos:exactMatch aro:{aro_local} ;")
        for m in r["mappings"]:
            if m["target_db"] == "CARD":
                continue
            tid = m["target_id"].replace('"', '\\"')
            lines.append(f'    {m["relation"]} "{m["target_db"]}:{tid}" ;')
        lines[-1] = lines[-1].rstrip(" ;") + " ."
        lines.append("")
    path.write_text("\n".join(lines))


def main():
    OUT.mkdir(exist_ok=True)
    records = load_records()
    print(f"Loaded {len(records)} validated records")
    write_tsv(records, OUT / "crosswalk.tsv")
    write_json(records, OUT / "crosswalk.json")
    write_jsonld(records, OUT / "crosswalk.jsonld")
    write_turtle(records, OUT / "crosswalk.ttl")
    print(f"Wrote derived artifacts to {OUT}")


if __name__ == "__main__":
    main()
