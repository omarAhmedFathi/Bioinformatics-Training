#!/bin/bash
# ==========================================
# Genome-wide LD Decay Analysis
# Compute pairwise r² within 500kb windows across all chromosomes
# ==========================================

DATA="../data/Qatari156_filtered_pruned"
OUTDIR="../intermediate/ld_decay"
mkdir -p "$OUTDIR"

echo "=== Computing genome-wide pairwise LD (r² within 500kb windows) ==="

plink --bfile "$DATA" \
      --r2 \
      --ld-window-r2 0 \
      --ld-window 9999 \
      --ld-window-kb 500 \
      --out "$OUTDIR/genomewide_ld"

echo "=== LD computation done ==="
echo "Output file: $OUTDIR/genomewide_ld.ld"
wc -l "$OUTDIR/genomewide_ld.ld"
