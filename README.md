# GIP Cancer Trajectory Modeling

Code for building genetically influenced proteotype (GIP) models from pQTL SNPs and evaluating cancer-site trajectory prediction.

This repository intentionally excludes raw genotype data, diagnosis event data, model caches, trained checkpoints, and literature supplements. Rebuild those locally from your authorized data sources.

## Main Workflow

1. Build pQTL edge tables from supplement sheets and matched rsIDs:

```bash
/usr/bin/python3 scripts/prepare_st10_pqtl_edges.py
```

2. Build Human Protein Atlas tissue priors for ST10 proteins:

```bash
/usr/bin/python3 scripts/build_hpa_protein_tissue_priors.py
```

3. Match SNP tissue-regulatory evidence from GTEx cis pQTL-eQTL colocalization:

```bash
/usr/bin/python3 scripts/build_snp_tissue_regulatory_priors.py
```

4. Build comprehensive pQTL feature tables:

```bash
/usr/bin/python3 scripts/build_comprehensive_pqtl_features.py
```

This creates edge-, SNP-, and protein-level feature tables under
`model_data/st10_core/features/`, including pQTL strength, consequence class,
fine-mapping/PIP, regulatory noncoding annotations, GTEx colocalized tissues,
pQTL-pQTL colocalization, PPI/receptor-ligand trans-network evidence,
covariate robustness, protein heritability, HPA tissue specificity, and
secretome/membrane labels.

5. Build cancer-site trajectory sequences from cleaned diagnosis events:

```bash
Rscript scripts/build_cancer_site_sequences.R
```

6. Prepare a training cache:

```bash
/usr/bin/python3 scripts/prepare_st10_training_cache.py \
  --max-samples 131072 \
  --max-proteins 1024 \
  --prior-version mode_tissue_reg \
  --out model_data/st10_core/st10_training_cache_n131072_p1024_time_traj_prior.npz
```

Use `--prior-version base` for the original beta/logp/cis/impact/panel prior, and
`--prior-version mode_tissue` for HPA protein tissue specificity without SNP-tissue
regulatory evidence. By default, the cache also includes standardized engineered
edge features from `model_data/st10_core/features/st10_edge_features.csv`; pass
`--no-edge-features` to disable them.

7. Train a time-aware trajectory model:

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

8. Evaluate strict and relaxed trajectory metrics:

```bash
/usr/bin/python3 scripts/evaluate_st10_endpoints.py \
  --cache model_data/st10_core/st10_training_cache_n131072_p1024_time_traj_prior.npz \
  --run-dir model_runs/st10_n131072_p1024_time_traj_prior_nope_mps_12ep
```

9. Infer proteotype clusters and trajectory associations:

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
