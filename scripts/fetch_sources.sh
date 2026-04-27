#!/usr/bin/env bash
# Download AMR reference databases at pinned versions.
# Versions match data/sources.yaml.
set -euo pipefail

DEST="${1:-ingest/sources}"
mkdir -p "$DEST"
cd "$DEST"

echo "==> CARD (broadstreet-v4.0.1)"
if [ ! -f card.tar.bz2 ]; then
  wget -q --show-progress -O card.tar.bz2 \
    "https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2"
fi
mkdir -p card && tar -xjf card.tar.bz2 -C card

echo "==> AMRFinderPlus reference catalog"
if [ ! -d amrfinderplus ]; then
  mkdir -p amrfinderplus && cd amrfinderplus
  wget -q --show-progress \
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt.fa" \
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS" \
    "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt" || true
  cd ..
fi

echo "==> ResFinder DB"
if [ ! -d resfinder_db ]; then
  git clone --depth 1 https://bitbucket.org/genomicepidemiology/resfinder_db.git
fi

echo "==> MEGARes 3.0"
if [ ! -f megares_v3.fasta ]; then
  wget -q --show-progress -O megares_v3.fasta \
    "https://www.meglab.org/downloads/megares_v3.00/megares_database_v3.00.fasta" || \
  wget -q --show-progress -O megares_v3.fasta \
    "https://www.meglab.org/downloads/megares_v3.00/megares_drugs_database_v3.00.fasta" || \
  echo "WARN: MEGARes download failed, edit URL manually"
fi

echo
echo "Done. Sources at: $(pwd)"
ls -lh
