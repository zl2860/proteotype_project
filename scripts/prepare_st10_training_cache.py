#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


META_COLS = 6


def raw_snp_from_col(col: str) -> str:
    # PLINK --recode A uses rsid_allele.
    return col.rsplit("_", 1)[0]


def read_raw_header(path: Path) -> list[str]:
    with path.open() as fh:
        return fh.readline().rstrip("\n").split()


def time_bin(days: int | float) -> int:
    if pd.isna(days):
        return 0
    days = int(days)
    if days <= 0:
        return 1
    if days <= 30:
        return 2
    if days <= 180:
        return 3
    if days <= 365:
        return 4
    if days <= 365 * 2:
        return 5
    if days <= 365 * 5:
        return 6
    return 7


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=Path("st10_core_rsid/raw_by_chr"))
    parser.add_argument("--edges", type=Path, default=Path("model_data/st10_core/st10_pqtl_edges.csv"))
    parser.add_argument("--sequences", type=Path, default=Path("model_data/traj/cancer_site_sequences.csv"))
    parser.add_argument("--vocab", type=Path, default=Path("model_data/traj/cancer_site_token_vocab.csv"))
    parser.add_argument("--out", type=Path, default=Path("model_data/st10_core/st10_training_cache_smoke.npz"))
    parser.add_argument("--max-samples", type=int, default=4096)
    parser.add_argument("--max-proteins", type=int, default=384)
    parser.add_argument("--max-seq-len", type=int, default=32)
    parser.add_argument("--top-index-sites", type=int, default=20)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    edges = pd.read_csv(args.edges)
    vocab = pd.read_csv(args.vocab)
    seq_df = pd.read_csv(args.sequences)
    events_path = args.sequences.with_name("cancer_site_events.csv")
    events = pd.read_csv(events_path, parse_dates=["diagnosis_date"])
    seq_map = {
        int(r.eid): [int(x) for x in str(r.token_ids).split() if x.strip()]
        for r in seq_df.itertuples(index=False)
    }
    events = events.sort_values(["eid", "diagnosis_date", "event_index", "token"])
    event_map = {
        int(eid): grp[["token_id", "diagnosis_date"]].reset_index(drop=True)
        for eid, grp in events.groupby("eid", sort=False)
    }

    # Pick proteins with strongest available ST10 support.
    protein_rank = (
        edges.groupby("protein_id")
        .agg(n_edges=("rsid", "size"), max_log10p=("log10p", "max"), sum_abs_beta=("abs_beta", "sum"))
        .reset_index()
        .sort_values(["max_log10p", "n_edges", "sum_abs_beta"], ascending=False)
    )
    protein_ids = protein_rank.head(args.max_proteins)["protein_id"].tolist()
    protein_to_idx = {p: i for i, p in enumerate(protein_ids)}
    edges = edges[edges["protein_id"].isin(protein_to_idx)].copy()
    target_rsids = set(edges["rsid"].astype(str))

    raw_files = sorted(args.raw_dir.glob("st10_core_rsid.chr*.raw"), key=lambda p: int(p.stem.split("chr")[-1]))
    if not raw_files:
        raise FileNotFoundError(f"No raw files found in {args.raw_dir}")

    # Determine sample order from chr1.
    eids = []
    with raw_files[0].open() as fh:
        header = fh.readline().split()
        for i, line in enumerate(fh):
            if i >= args.max_samples:
                break
            parts = line.split(maxsplit=META_COLS)
            eids.append(int(parts[1]))
    n = len(eids)
    eid_to_row = {eid: i for i, eid in enumerate(eids)}

    snp_rsids: list[str] = []
    snp_chunks: list[np.ndarray] = []

    for raw in raw_files:
        header = read_raw_header(raw)
        raw_rsids = [raw_snp_from_col(c) for c in header[META_COLS:]]
        keep_idx = [i for i, rsid in enumerate(raw_rsids) if rsid in target_rsids]
        if not keep_idx:
            continue
        keep_cols = [META_COLS + i for i in keep_idx]
        keep_rsids = [raw_rsids[i] for i in keep_idx]
        mat = np.full((n, len(keep_cols)), np.nan, dtype=np.float32)
        with raw.open() as fh:
            fh.readline()
            for line_i, line in enumerate(fh):
                if line_i >= n:
                    break
                parts = line.rstrip("\n").split()
                for j, col_i in enumerate(keep_cols):
                    v = parts[col_i]
                    if v != "NA":
                        mat[line_i, j] = float(v)
        snp_rsids.extend(keep_rsids)
        snp_chunks.append(mat)

    if not snp_chunks:
        raise RuntimeError("No genotype columns matched pQTL edge rsIDs.")

    X = np.concatenate(snp_chunks, axis=1)
    # Mean-impute missing values for this first model.
    means = np.nanmean(X, axis=0)
    missing = np.isnan(X)
    X[missing] = np.take(means, np.where(missing)[1])
    X = X.astype(np.float32)

    snp_to_idx = {rsid: i for i, rsid in enumerate(snp_rsids)}
    edges = edges[edges["rsid"].isin(snp_to_idx)].copy()
    edge_snp_idx = edges["rsid"].map(snp_to_idx).astype("int64").to_numpy()
    edge_protein_idx = edges["protein_id"].map(protein_to_idx).astype("int64").to_numpy()
    edge_beta = edges["beta"].fillna(0).astype("float32").to_numpy()
    edge_logp = edges["log10p"].fillna(0).astype("float32").to_numpy()
    edge_cis = edges["cis_trans"].eq("cis").astype("float32").to_numpy()
    impact_boost = edges["impact"].fillna("").astype(str).str.upper().map(
        {"HIGH": 1.5, "MODERATE": 1.25, "LOW": 1.1, "MODIFIER": 1.0}
    ).fillna(1.0)
    panel_boost = edges["protein_panel"].fillna("").astype(str).str.contains("Oncology", case=False).map(
        {True: 1.15, False: 1.0}
    )
    prior = edges["beta"].abs().fillna(0) * np.sqrt(edges["log10p"].fillna(0).clip(lower=0) + 1.0)
    prior = prior * (1.0 + 0.5 * edges["cis_trans"].eq("cis").astype(float)) * impact_boost * panel_boost
    prior = prior.astype("float32").to_numpy()
    prior = np.log1p(prior / max(float(np.nanmedian(prior[prior > 0])) if np.any(prior > 0) else 1.0, 1e-6))

    pad_id, bos_id, eos_id = 0, 1, 2
    Y = np.full((n, args.max_seq_len), pad_id, dtype=np.int64)
    time_bins = np.zeros((n, args.max_seq_len), dtype=np.int64)
    y_len = np.zeros(n, dtype=np.int64)
    has_events = np.zeros(n, dtype=np.int64)
    for i, eid in enumerate(eids):
        event_df = event_map.get(eid)
        toks = seq_map.get(eid, [])
        if toks:
            has_events[i] = 1
        if event_df is not None:
            event_df = event_df.head(max(0, args.max_seq_len - 2))
            toks = event_df["token_id"].astype(int).tolist()
            dates = pd.to_datetime(event_df["diagnosis_date"]).tolist()
        else:
            dates = []
        seq = [bos_id] + toks + [eos_id]
        Y[i, : len(seq)] = seq
        y_len[i] = len(seq)
        if toks:
            time_bins[i, 1] = 1
            for j in range(1, len(toks)):
                time_bins[i, j + 1] = time_bin((dates[j] - dates[j - 1]).days)

    index_events = events[events["stage"].eq("index_cancer")].copy()
    site_rank = (
        index_events.groupby("cancer_site")["eid"]
        .nunique()
        .sort_values(ascending=False)
        .head(args.top_index_sites)
    )
    site_names = site_rank.index.astype(str).tolist()
    site_to_class = {site: i + 1 for i, site in enumerate(site_names)}
    first_index_site_class = np.zeros(n, dtype=np.int64)
    site_binary = np.zeros((n, len(site_names)), dtype=np.float32)
    if site_names:
        first_index = index_events[index_events["cancer_site"].isin(site_to_class)].sort_values(
            ["eid", "diagnosis_date", "cancer_site"]
        )
        first_index = first_index.drop_duplicates("eid", keep="first")
        first_site_map = dict(zip(first_index["eid"].astype(int), first_index["cancer_site"].astype(str)))
        any_site = index_events[index_events["cancer_site"].isin(site_to_class)]
        any_site_map = any_site.groupby("eid")["cancer_site"].agg(lambda s: sorted(set(map(str, s)))).to_dict()
        for i, eid in enumerate(eids):
            site = first_site_map.get(eid)
            if site is not None:
                first_index_site_class[i] = site_to_class[site]
            for site in any_site_map.get(eid, []):
                site_binary[i, site_to_class[site] - 1] = 1.0

    np.savez_compressed(
        args.out,
        X=X,
        eids=np.array(eids, dtype=np.int64),
        Y=Y,
        y_len=y_len,
        has_events=has_events,
        time_bins=time_bins,
        time_bin_labels=np.array(
            [
                "<PAD_OR_EOS>",
                "first_or_same_day",
                "1_30_days",
                "31_180_days",
                "181_365_days",
                "1_2_years",
                "2_5_years",
                "gt_5_years",
            ],
            dtype=object,
        ),
        snp_rsids=np.array(snp_rsids, dtype=object),
        protein_ids=np.array(protein_ids, dtype=object),
        edge_snp_idx=edge_snp_idx,
        edge_protein_idx=edge_protein_idx,
        edge_beta=edge_beta,
        edge_logp=edge_logp,
        edge_cis=edge_cis,
        edge_prior=prior.astype(np.float32),
        vocab_size=np.array([int(vocab["token_id"].max()) + 1], dtype=np.int64),
        max_seq_len=np.array([args.max_seq_len], dtype=np.int64),
        site_names=np.array(site_names, dtype=object),
        first_index_site_class=first_index_site_class,
        site_binary=site_binary,
    )

    summary = args.out.with_suffix(".summary.txt")
    with summary.open("w") as out:
        out.write(f"samples: {n}\n")
        out.write(f"samples_with_cancer_site_events: {int(has_events.sum())}\n")
        out.write(f"snps: {X.shape[1]}\n")
        out.write(f"proteins: {len(protein_ids)}\n")
        out.write(f"pqtl_edges: {len(edge_snp_idx)}\n")
        out.write(f"vocab_size: {int(vocab['token_id'].max()) + 1}\n")
        out.write(f"max_seq_len: {args.max_seq_len}\n")
        out.write(f"time_bins: 8\n")
        out.write(f"top_index_sites: {len(site_names)}\n")
        out.write(f"samples_with_top_index_site: {int((first_index_site_class > 0).sum())}\n")
    print(summary.read_text())


if __name__ == "__main__":
    main()
