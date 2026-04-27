# Contributing to AMRxref

Two ways to contribute:

1. **Add or fix a mapping.** Open an issue using the "New mapping" template,
   or send a PR adding/editing a YAML file in `data/genes/`. Every mapping
   needs evidence (sequence identity, literature, or curator endorsement).
2. **Improve code or docs.** Standard PR workflow.

## Mapping evidence tiers
- gold: manually reviewed, two or more independent evidence types
- silver: single strong evidence (100% sequence identity)
- bronze: automated only, awaits review

See `docs/methods.md` for full curation guidelines.
