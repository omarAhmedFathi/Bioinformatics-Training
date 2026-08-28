#!/bin/bash
# ==========================================
# Compute per-chromosome LD square matrices
# ==========================================

DATA="../data/Qatari156_filtered_pruned"
OUTDIR="../intermediate/ld_decay/per_chr"
mkdir -p "$OUTDIR"

for CHR in $(seq 1 22) 23; do
  echo "=== Chromosome $CHR ==="
  plink --bfile "$DATA" \
        --chr "$CHR" \
        --r2 square \
        --out "$OUTDIR/chr${CHR}_ld" \
        --silent
  echo "  Done. Matrix: $(wc -l < "$OUTDIR/chr${CHR}_ld.ld") SNPs"
done

echo "=== All chromosomes done ==="
