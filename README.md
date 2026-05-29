# GIP Cancer Trajectory Modeling

Code for building genetically influenced proteotype (GIP) models from pQTL SNPs and evaluating cancer-site trajectory prediction.

This repository intentionally excludes raw genotype data, diagnosis event data, model caches, trained checkpoints, and literature supplements. Rebuild those locally from your authorized data sources.

## Main Workflow

1. Build pQTL edge tables from supplement sheets and matched rsIDs:

```bash
/usr/bin/python3 scripts/prepare_st10_pqtl_edges.py
```

2. Build cancer-site trajectory sequences from cleaned diagnosis events:

```bash
Rscript scripts/build_cancer_site_sequences.R
```

3. Prepare a training cache:

```bash
/usr/bin/python3 scripts/prepare_st10_training_cache.py \
  --max-samples 131072 \
  --max-proteins 1024 \
  --out model_data/st10_core/st10_training_cache_n131072_p1024_time_traj_prior.npz
```

4. Train a time-aware trajectory model:

```bash
/usr/bin/python3 scripts/train_st10_gip_trajectory_transformer.py \
  --cache model_data/st10_core/st10_training_cache_n131072_p1024_time_traj_prior.npz \
  --outdir model_runs/st10_n131072_p1024_time_traj_prior_nope_mps_12ep \
  --device mps \
  --epochs 12 \
  --batch-size 128 \
  --d-model 64 \
  --decoder-layers 2 \
  --protein-encoder-layers 0 \
  --time-loss-weight 0.5 \
  --event-loss-weight 0 \
  --site-loss-weight 0 \
  --site-binary-loss-weight 0 \
  --prior-strength 1.0
```

5. Evaluate strict and relaxed trajectory metrics:

```bash
/usr/bin/python3 scripts/evaluate_st10_endpoints.py \
  --cache model_data/st10_core/st10_training_cache_n131072_p1024_time_traj_prior.npz \
  --run-dir model_runs/st10_n131072_p1024_time_traj_prior_nope_mps_12ep
```

6. Infer proteotype clusters and trajectory associations:

```bash
/usr/bin/python3 scripts/infer_st10_proteotype_trajectory_links.py \
  --cache model_data/st10_core/st10_training_cache_n131072_p1024_time_traj_prior.npz \
  --run-dir model_runs/st10_n131072_p1024_time_traj_prior_nope_mps_12ep \
  --outdir model_runs/st10_n131072_p1024_time_traj_prior_nope_mps_12ep/proteotype_inference_k8_zscore
```

## Live Training Dashboard

```bash
/usr/bin/python3 scripts/serve_training_dashboard.py --host 127.0.0.1 --port 8765 --root .
```

Open <http://127.0.0.1:8765>. The dashboard scans `model_runs/*/history.jsonl` and refreshes every few seconds.

## Repository Contents

- `scripts/`: data preparation, model training, evaluation, dashboard, and proteotype inference.
- `diagnosis_traj/`: lightweight notes and ICD/site trajectory mapping documentation.
- `plink_variant_lists/`: lightweight extraction manifests and server-side helper scripts.

## Excluded Data

The following are ignored by git:

- `st10_core_rsid/`
- `model_data/`
- `model_runs/`
- raw diagnosis `.qs2` files
- literature PDFs, XLSX supplements, and other large/private artifacts

Keep those data files local or in controlled storage, not in GitHub.
