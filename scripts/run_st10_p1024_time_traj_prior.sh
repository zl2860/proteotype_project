#!/usr/bin/env zsh
set -euo pipefail

RUN="model_runs/st10_n131072_p1024_time_traj_prior_nope_mps_12ep"
CACHE="model_data/st10_core/st10_training_cache_n131072_p1024_time_traj_prior.npz"

mkdir -p "$RUN"

/usr/bin/python3 scripts/prepare_st10_training_cache.py \
  --raw-dir st10_core_rsid/raw_by_chr \
  --edges model_data/st10_core/st10_pqtl_edges.csv \
  --sequences model_data/traj/cancer_site_sequences.csv \
  --vocab model_data/traj/cancer_site_token_vocab.csv \
  --out "$CACHE" \
  --max-samples 131072 \
  --max-proteins 1024 \
  --max-seq-len 32 \
  --top-index-sites 20

/usr/bin/python3 scripts/train_st10_gip_trajectory_transformer.py \
  --cache "$CACHE" \
  --outdir "$RUN" \
  --device mps \
  --epochs 12 \
  --batch-size 128 \
  --d-model 64 \
  --decoder-layers 2 \
  --protein-encoder-layers 0 \
  --val-frac 0.1 \
  --seed 2026 \
  --time-loss-weight 0.5 \
  --event-loss-weight 0 \
  --site-loss-weight 0 \
  --site-binary-loss-weight 0 \
  --prior-strength 1.0
