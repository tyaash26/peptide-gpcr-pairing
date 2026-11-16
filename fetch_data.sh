#!/bin/sh
# Download the GtoPdb bulk files used by build_pairs_v3.py.
#
# Note: guidetopharmacology.org/download.jsp is behind a login as of Nov 2026,
# but the static files under /DATA/ are served without one. The GPCR
# interaction file is lowercase (gpcr_interactions.csv); GPCR_interactions.csv
# returns 404.
#
# Data: IUPHAR/BPS Guide to PHARMACOLOGY (GtoPdb), ODbL / CC-BY-SA 4.0.
# Reference: Harding SD et al. (2026) Nucleic Acids Res 54(D1):D1446-D1456.

set -e
DIR="${1:-data/gtopdb}"
mkdir -p "$DIR"
cd "$DIR"

for f in gpcr_interactions.csv peptides.csv targets_and_families.csv ligand_id_mapping.csv; do
    echo "fetching $f"
    curl -fsSL -A "Mozilla/5.0" -O "https://www.guidetopharmacology.org/DATA/$f"
done

ls -lh
