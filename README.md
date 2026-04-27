# AMRxref

**Canonical cross-reference for antimicrobial resistance gene names.**

AMRxref bridges the four major AMR gene databases — [CARD](https://card.mcmaster.ca/),
[AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/),
[ResFinder](https://cge.food.dtu.dk/services/ResFinder/), and
[MEGARes](https://www.meglab.org/megares/) — with curated SKOS mappings
anchored to the Antibiotic Resistance Ontology (ARO).

> Run three AMR detection tools on the same isolate and you get three different
> gene names for the same gene. AMRxref tells you which names refer to the same
> thing, with evidence and provenance.

## Why

Studies comparing AMRFinderPlus and ResFinder on identical isolates have found
gene-symbol disagreement on roughly 9% of calls. β-lactamases alone account for
much of this drift because the same allele often gets a fresh name in each new
paper. [hAMRonization](https://github.com/pha4ge/hAMRonization) standardises
the *report schema* across tools but explicitly does not standardise the
*gene names* themselves. AMRxref fills that last-mile gap.

## What you get

- One canonical name per gene, anchored to an ARO ID
- Cross-references to AMRFinderPlus, ResFinder, MEGARes, NCBI, UniProt
- Evidence for every mapping: sequence identity, length, literature
- Provenance: who curated, when, against which database version
- Quality tier: gold (reviewed), silver (auto, strong evidence), bronze (auto)
- Standards-compliant outputs: TSV, JSON, JSON-LD, SKOS Turtle

## Scope

AMRxref covers acquired AMR genes and resistance-conferring point mutations in
bacteria, with sources from CARD, AMRFinderPlus, ResFinder, and MEGARes. It
covers protein-coding genes, rRNA mutations, and promoter/regulatory variants.

Mycobacterium tuberculosis resistance is out of scope (different sources,
different conventions); a sibling project is planned.

## Quick start

\`\`\`bash
pip install amrxref
amrxref lookup blaTEM-1
\`\`\`

## FAIR

- **F** Persistent IRIs at https://w3id.org/amrxref/, Zenodo DOI per release
- **A** HTTPS, content negotiation (HTML, JSON, JSON-LD, Turtle, TSV)
- **I** Anchored to ARO; SKOS for cross-database mappings; aligned with PHA4GE/hAMRonization
- **R** CC0 for data, MIT for code; per-mapping provenance; SemVer

## Cite

See CITATION.cff. Each release also has a Zenodo DOI.

## Contribute

See CONTRIBUTING.md.

## Status

v0.0 scaffold. v0.1 target: automated bronze-tier records covering all genes
in the four source databases at pinned versions.
