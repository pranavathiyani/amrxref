"""Command-line interface for AMRxref."""
import json
import sys
from pathlib import Path
import click

DERIVED = Path(__file__).resolve().parent.parent.parent / "derived" / "crosswalk.json"


def load_data():
    if not DERIVED.exists():
        click.echo(f"No data found at {DERIVED}. Run scripts/build_derived.py first.", err=True)
        sys.exit(1)
    return json.loads(DERIVED.read_text())


@click.group()
@click.version_option()
def main():
    """AMRxref: canonical cross-reference for AMR gene names."""


@main.command()
@click.argument("query")
def lookup(query):
    """Look up a gene by canonical name, ARO ID, or any database ID."""
    data = load_data()
    q = query.lower()
    matches = [
        r for r in data
        if q in r["canonical_name"].lower()
        or q in r["aro_id"].lower()
        or any(q in m["target_id"].lower() for m in r["mappings"])
    ]
    if not matches:
        click.echo(f"No match for {query!r}", err=True)
        sys.exit(1)
    for r in matches:
        click.echo(json.dumps(r, indent=2, default=str))


@main.command(name="list")
@click.option("--quality", type=click.Choice(["gold", "silver", "bronze"]))
def list_cmd(quality):
    """List all canonical gene names, optionally filtered by quality."""
    data = load_data()
    for r in data:
        qs = [m["quality"] for m in r["mappings"]]
        overall = "gold" if all(q == "gold" for q in qs) else "bronze" if "bronze" in qs else "silver"
        if quality and overall != quality:
            continue
        click.echo(f"{r['aro_id']}\t{r['canonical_name']}\t{overall}")


@main.command()
def stats():
    """Show summary statistics."""
    data = load_data()
    total = len(data)
    by_q = {"gold": 0, "silver": 0, "bronze": 0}
    by_v = {}
    for r in data:
        qs = [m["quality"] for m in r["mappings"]]
        overall = "gold" if all(q == "gold" for q in qs) else "bronze" if "bronze" in qs else "silver"
        by_q[overall] += 1
        v = r.get("variant_type", "unknown")
        by_v[v] = by_v.get(v, 0) + 1
    click.echo(f"Total records: {total}")
    click.echo("By quality:")
    for k, v in by_q.items():
        click.echo(f"  {k}: {v}")
    click.echo("By variant type:")
    for k, v in by_v.items():
        click.echo(f"  {k}: {v}")


if __name__ == "__main__":
    main()
