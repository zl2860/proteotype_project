#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.cluster import KMeans
from torch.utils.data import DataLoader, Subset


def load_trainer(path: Path):
    spec = importlib.util.spec_from_file_location("st10_trainer", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load trainer module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_model(trainer, z, ckpt, device: str):
    model_args = ckpt["args"]
    site_names = np.array(ckpt.get("site_names", z["site_names"] if "site_names" in z else []), dtype=object)
    time_labels = np.array(ckpt.get("time_bin_labels", z["time_bin_labels"] if "time_bin_labels" in z else []), dtype=object)
    encoder = trainer.ProteinAwareGIPEncoder(
        n_snps=z["X"].shape[1],
        n_proteins=len(z["protein_ids"]),
        edge_snp_idx=z["edge_snp_idx"],
        edge_protein_idx=z["edge_protein_idx"],
        edge_beta=z["edge_beta"],
        edge_logp=z["edge_logp"],
        edge_cis=z["edge_cis"],
        edge_prior=z["edge_prior"] if "edge_prior" in z else None,
        d_model=int(model_args["d_model"]),
        dropout=float(model_args["dropout"]),
        prior_strength=float(model_args.get("prior_strength", 1.0)),
        protein_encoder_layers=int(model_args.get("protein_encoder_layers", 1)),
    )
    model = trainer.GIPTrajectoryDecoder(
        encoder=encoder,
        vocab_size=int(z["vocab_size"][0]),
        d_model=int(model_args["d_model"]),
        n_layers=int(model_args["decoder_layers"]),
        dropout=float(model_args["dropout"]),
        n_site_classes=len(site_names) + 1,
        n_site_binary=len(site_names),
        n_time_bins=len(time_labels) or 1,
    )
    model.load_state_dict(ckpt["model_state_dict"])
    model.to(device)
    model.eval()
    return model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--vocab", type=Path, default=Path("model_data/traj/cancer_site_token_vocab.csv"))
    parser.add_argument("--edges", type=Path, default=Path("model_data/st10_core/st10_pqtl_edges.csv"))
    parser.add_argument("--trainer", type=Path, default=Path("scripts/train_st10_gip_trajectory_transformer.py"))
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--n-clusters", type=int, default=8)
    parser.add_argument("--max-samples", type=int, default=0)
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    trainer = load_trainer(args.trainer)
    z = np.load(args.cache, allow_pickle=True)
    ckpt = torch.load(args.run_dir / "model.pt", map_location="cpu", weights_only=False)
    model = build_model(trainer, z, ckpt, args.device)
    dataset = trainer.TrajDataset(args.cache)
    if args.max_samples and args.max_samples < len(dataset):
        dataset = Subset(dataset, list(range(args.max_samples)))
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False)

    gip_chunks, eid_chunks, y_chunks, time_chunks = [], [], [], []
    edge_dosage_sum = np.zeros(len(z["edge_snp_idx"]), dtype=np.float64)
    n_seen = 0
    with torch.no_grad():
        offset = 0
        for X, Y, time_bins, *_ in loader:
            X_dev = X.to(args.device)
            gip, _, _ = model.encoder(X_dev)
            gip_chunks.append(gip.cpu().numpy())
            batch_n = len(X)
            eid_chunks.append(z["eids"][offset : offset + batch_n])
            y_chunks.append(Y.numpy())
            time_chunks.append(time_bins.numpy())
            edge_dosage_sum += X[:, z["edge_snp_idx"]].numpy().sum(axis=0)
            n_seen += batch_n
            offset += batch_n

    gip = np.vstack(gip_chunks)
    eids = np.concatenate(eid_chunks)
    Y = np.vstack(y_chunks)
    time_bins = np.vstack(time_chunks)
    edge_dosage_mean = edge_dosage_sum / max(n_seen, 1)

    gip_scaled = (gip - gip.mean(axis=0, keepdims=True)) / np.maximum(gip.std(axis=0, keepdims=True), 1e-6)
    km = KMeans(n_clusters=args.n_clusters, random_state=args.seed, n_init=20)
    clusters = km.fit_predict(gip_scaled)
    pd.DataFrame({"eid": eids, "proteotype_cluster": clusters}).to_csv(
        args.outdir / "proteotype_clusters.csv", index=False
    )

    token_vocab = pd.read_csv(args.vocab).set_index("token_id")["token"].to_dict()
    time_labels = list(z["time_bin_labels"]) if "time_bin_labels" in z else []
    cluster_rows = []
    for c in range(args.n_clusters):
        mask = clusters == c
        cluster_rows.append(
            {
                "proteotype_cluster": c,
                "n": int(mask.sum()),
                "fraction": float(mask.mean()),
                "gip_norm_mean": float(np.linalg.norm(gip[mask], axis=1).mean()),
            }
        )
    pd.DataFrame(cluster_rows).to_csv(args.outdir / "proteotype_cluster_summary.csv", index=False)

    target_tokens = Y[:, 1:].reshape(-1)
    token_subject_clusters = np.repeat(clusters, Y.shape[1] - 1)
    token_keep = ~np.isin(target_tokens, [0, 1, 2])
    overall_token_counts = Counter(target_tokens[token_keep].tolist())
    traj_rows = []
    for c in range(args.n_clusters):
        c_tokens = target_tokens[token_keep & (token_subject_clusters == c)]
        c_counts = Counter(c_tokens.tolist())
        n_cluster = max(int((clusters == c).sum()), 1)
        n_total = max(len(clusters), 1)
        for token_id, count in c_counts.most_common(args.top_n):
            cluster_rate = count / n_cluster
            overall_rate = overall_token_counts[token_id] / n_total
            traj_rows.append(
                {
                    "proteotype_cluster": c,
                    "token_id": int(token_id),
                    "token": token_vocab.get(int(token_id), str(token_id)),
                    "cluster_count": int(count),
                    "cluster_rate_per_subject": cluster_rate,
                    "overall_rate_per_subject": overall_rate,
                    "log2_enrichment": float(np.log2((cluster_rate + 1e-6) / (overall_rate + 1e-6))),
                }
            )
    pd.DataFrame(traj_rows).sort_values(["proteotype_cluster", "log2_enrichment"], ascending=[True, False]).to_csv(
        args.outdir / "cluster_trajectory_token_enrichment.csv", index=False
    )

    target_times = time_bins[:, 1:].reshape(-1)
    time_subject_clusters = np.repeat(clusters, time_bins.shape[1] - 1)
    time_keep = target_times != 0
    overall_time_counts = Counter(target_times[time_keep].tolist())
    time_rows = []
    for c in range(args.n_clusters):
        c_times = target_times[time_keep & (time_subject_clusters == c)]
        c_counts = Counter(c_times.tolist())
        n_cluster = max(int((clusters == c).sum()), 1)
        n_total = max(len(clusters), 1)
        for time_id, count in c_counts.most_common():
            cluster_rate = count / n_cluster
            overall_rate = overall_time_counts[time_id] / n_total
            label = time_labels[int(time_id)] if int(time_id) < len(time_labels) else str(time_id)
            time_rows.append(
                {
                    "proteotype_cluster": c,
                    "time_bin": int(time_id),
                    "time_label": label,
                    "cluster_count": int(count),
                    "cluster_rate_per_subject": cluster_rate,
                    "overall_rate_per_subject": overall_rate,
                    "log2_enrichment": float(np.log2((cluster_rate + 1e-6) / (overall_rate + 1e-6))),
                }
            )
    pd.DataFrame(time_rows).sort_values(["proteotype_cluster", "log2_enrichment"], ascending=[True, False]).to_csv(
        args.outdir / "cluster_time_bin_enrichment.csv", index=False
    )

    edge_attention_sum = np.zeros((args.n_clusters, len(z["edge_snp_idx"])), dtype=np.float64)
    edge_dosage_sum_by_cluster = np.zeros_like(edge_attention_sum)
    cluster_counts = np.bincount(clusters, minlength=args.n_clusters).astype(np.float64)
    with torch.no_grad():
        offset = 0
        for X, *_ in loader:
            X_dev = X.to(args.device)
            _, attention, _ = model.encoder(X_dev)
            attention = attention.cpu().numpy()
            edge_dosage = X[:, z["edge_snp_idx"]].numpy()
            batch_clusters = clusters[offset : offset + len(X)]
            for c in range(args.n_clusters):
                mask = batch_clusters == c
                if mask.any():
                    edge_attention_sum[c] += attention[mask].sum(axis=0)
                    edge_dosage_sum_by_cluster[c] += edge_dosage[mask].sum(axis=0)
            offset += len(X)

    edge_table = pd.DataFrame(
        {
            "edge_index": np.arange(len(z["edge_snp_idx"])),
            "rsid": np.array(z["snp_rsids"], dtype=object)[z["edge_snp_idx"]],
            "protein_id": np.array(z["protein_ids"], dtype=object)[z["edge_protein_idx"]],
            "beta": z["edge_beta"],
            "log10p": z["edge_logp"],
            "cis": z["edge_cis"],
            "prior": z["edge_prior"] if "edge_prior" in z else np.zeros(len(z["edge_snp_idx"])),
        }
    )
    edge_rows = []
    for c in range(args.n_clusters):
        denom = max(cluster_counts[c], 1.0)
        mean_attention = edge_attention_sum[c] / denom
        mean_dosage = edge_dosage_sum_by_cluster[c] / denom
        dosage_delta = mean_dosage - edge_dosage_mean
        score = mean_attention * np.abs(dosage_delta) * (1.0 + np.abs(z["edge_beta"])) * (
            1.0 + np.log1p(z["edge_logp"])
        )
        top_idx = np.argsort(score)[::-1][: args.top_n]
        for edge_idx in top_idx:
            row = edge_table.iloc[int(edge_idx)].to_dict()
            row.update(
                {
                    "proteotype_cluster": c,
                    "cluster_mean_attention": float(mean_attention[edge_idx]),
                    "cluster_mean_dosage": float(mean_dosage[edge_idx]),
                    "overall_mean_dosage": float(edge_dosage_mean[edge_idx]),
                    "dosage_delta": float(dosage_delta[edge_idx]),
                    "cluster_edge_score": float(score[edge_idx]),
                }
            )
            edge_rows.append(row)
    pd.DataFrame(edge_rows).to_csv(args.outdir / "cluster_top_snp_protein_edges.csv", index=False)

    with (args.outdir / "run_metadata.json").open("w") as out:
        json.dump(
            {
                "cache": str(args.cache),
                "run_dir": str(args.run_dir),
                "n_samples": int(len(clusters)),
                "n_clusters": args.n_clusters,
                "top_n": args.top_n,
            },
            out,
            indent=2,
        )


if __name__ == "__main__":
    main()
