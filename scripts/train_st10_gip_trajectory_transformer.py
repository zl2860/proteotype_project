#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import torch
from torch import nn
from torch.utils.data import DataLoader, Dataset, random_split


class TrajDataset(Dataset):
    def __init__(self, cache_path: Path):
        z = np.load(cache_path, allow_pickle=True)
        self.X = torch.tensor(z["X"], dtype=torch.float32)
        self.Y = torch.tensor(z["Y"], dtype=torch.long)
        self.time_bins = torch.tensor(
            z["time_bins"] if "time_bins" in z else np.zeros_like(z["Y"], dtype=np.int64),
            dtype=torch.long,
        )
        self.eids = z["eids"]
        self.has_events = torch.tensor(z["has_events"], dtype=torch.float32)
        self.first_index_site_class = torch.tensor(
            z["first_index_site_class"] if "first_index_site_class" in z else np.zeros(len(self.X), dtype=np.int64),
            dtype=torch.long,
        )
        self.site_binary = torch.tensor(
            z["site_binary"] if "site_binary" in z else np.zeros((len(self.X), 0), dtype=np.float32),
            dtype=torch.float32,
        )

    def __len__(self) -> int:
        return self.X.shape[0]

    def __getitem__(self, idx: int):
        return (
            self.X[idx],
            self.Y[idx],
            self.time_bins[idx],
            self.has_events[idx],
            self.first_index_site_class[idx],
            self.site_binary[idx],
        )


class ProteinAwareGIPEncoder(nn.Module):
    def __init__(
        self,
        n_snps: int,
        n_proteins: int,
        edge_snp_idx: np.ndarray,
        edge_protein_idx: np.ndarray,
        edge_beta: np.ndarray,
        edge_logp: np.ndarray,
        edge_cis: np.ndarray,
        edge_prior: np.ndarray | None,
        edge_features: np.ndarray | None,
        d_model: int,
        dropout: float,
        prior_strength: float,
        protein_encoder_layers: int,
    ):
        super().__init__()
        self.n_proteins = n_proteins
        self.d_model = d_model
        self.register_buffer("edge_snp_idx", torch.tensor(edge_snp_idx, dtype=torch.long))
        self.register_buffer("edge_protein_idx", torch.tensor(edge_protein_idx, dtype=torch.long))
        self.register_buffer("edge_beta", torch.tensor(edge_beta, dtype=torch.float32))
        self.register_buffer("edge_logp", torch.log1p(torch.tensor(edge_logp, dtype=torch.float32)))
        self.register_buffer("edge_cis", torch.tensor(edge_cis, dtype=torch.float32))
        if edge_prior is None:
            edge_prior = np.zeros_like(edge_beta, dtype=np.float32)
        self.register_buffer("edge_prior", torch.tensor(edge_prior, dtype=torch.float32))
        if edge_features is None:
            edge_features = np.zeros((len(edge_beta), 0), dtype=np.float32)
        self.register_buffer("edge_features", torch.tensor(edge_features, dtype=torch.float32))
        self.prior_strength = prior_strength

        self.snp_embedding = nn.Embedding(n_snps, d_model)
        self.protein_embedding = nn.Embedding(n_proteins, d_model)
        self.edge_mlp = nn.Sequential(
            nn.Linear(d_model + 4 + self.edge_features.shape[1], d_model),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(d_model, d_model),
        )
        self.attn_score = nn.Linear(d_model, 1)
        self.beta_bias = nn.Parameter(torch.tensor(0.1))
        self.logp_bias = nn.Parameter(torch.tensor(0.05))
        enc_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=4,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            batch_first=True,
            activation="gelu",
        )
        self.protein_encoder = (
            nn.TransformerEncoder(enc_layer, num_layers=protein_encoder_layers)
            if protein_encoder_layers > 0
            else nn.Identity()
        )
        self.out_norm = nn.LayerNorm(d_model)

        edge_to_protein = torch.zeros(len(edge_protein_idx), n_proteins, dtype=torch.float32)
        edge_to_protein[torch.arange(len(edge_protein_idx)), torch.tensor(edge_protein_idx, dtype=torch.long)] = 1.0
        self.register_buffer("edge_to_protein", edge_to_protein)

    def forward(self, dosage: torch.Tensor):
        bsz = dosage.shape[0]
        edge_dosage = dosage[:, self.edge_snp_idx]  # B x E
        snp_emb = self.snp_embedding(self.edge_snp_idx)  # E x D
        snp_emb = snp_emb.unsqueeze(0).expand(bsz, -1, -1)
        edge_feat = torch.stack(
            [
                edge_dosage,
                self.edge_beta.unsqueeze(0).expand(bsz, -1),
                self.edge_logp.unsqueeze(0).expand(bsz, -1),
                self.edge_cis.unsqueeze(0).expand(bsz, -1),
            ],
            dim=-1,
        )
        if self.edge_features.shape[1] > 0:
            edge_extra = self.edge_features.unsqueeze(0).expand(bsz, -1, -1)
            edge_feat = torch.cat([edge_feat, edge_extra], dim=-1)
        edge_hidden = self.edge_mlp(torch.cat([snp_emb, edge_feat], dim=-1))
        raw_scores = self.attn_score(edge_hidden).squeeze(-1)
        raw_scores = (
            raw_scores
            + self.beta_bias * self.edge_beta.abs()
            + self.logp_bias * self.edge_logp
            + self.prior_strength * self.edge_prior
        )

        exp_scores = torch.exp(raw_scores - raw_scores.max(dim=1, keepdim=True).values)
        protein_denoms = exp_scores @ self.edge_to_protein
        attention_by_edge = exp_scores / protein_denoms[:, self.edge_protein_idx].clamp_min(1e-8)
        weighted_edges = attention_by_edge.unsqueeze(-1) * edge_hidden
        proteins = torch.einsum("bed,ep->bpd", weighted_edges, self.edge_to_protein)
        proteins = proteins + self.protein_embedding.weight.unsqueeze(0)
        proteins = self.protein_encoder(proteins)
        gip = self.out_norm(proteins.mean(dim=1))
        return gip, attention_by_edge, proteins


class GIPTrajectoryDecoder(nn.Module):
    def __init__(
        self,
        encoder: ProteinAwareGIPEncoder,
        vocab_size: int,
        d_model: int,
        n_layers: int,
        dropout: float,
        n_site_classes: int = 1,
        n_site_binary: int = 0,
        n_time_bins: int = 1,
    ):
        super().__init__()
        self.encoder = encoder
        self.token_embedding = nn.Embedding(vocab_size, d_model, padding_idx=0)
        self.pos_embedding = nn.Embedding(256, d_model)
        dec_layer = nn.TransformerDecoderLayer(
            d_model=d_model,
            nhead=4,
            dim_feedforward=d_model * 4,
            dropout=dropout,
            batch_first=True,
            activation="gelu",
        )
        self.decoder = nn.TransformerDecoder(dec_layer, num_layers=n_layers)
        self.lm_head = nn.Linear(d_model, vocab_size)
        self.time_head = nn.Linear(d_model, n_time_bins)
        self.event_head = nn.Linear(d_model, 1)
        self.site_head = nn.Linear(d_model, n_site_classes)
        self.site_binary_head = nn.Linear(d_model, n_site_binary) if n_site_binary > 0 else None

    def forward(self, dosage: torch.Tensor, tokens_in: torch.Tensor):
        gip, attention, proteins = self.encoder(dosage)
        memory = gip.unsqueeze(1)
        pos = torch.arange(tokens_in.shape[1], device=tokens_in.device).unsqueeze(0)
        tgt = self.token_embedding(tokens_in) + self.pos_embedding(pos)
        causal = torch.triu(
            torch.full((tokens_in.shape[1], tokens_in.shape[1]), float("-inf"), device=tokens_in.device),
            diagonal=1,
        )
        decoded = self.decoder(tgt=tgt, memory=memory, tgt_mask=causal)
        logits = self.lm_head(decoded)
        time_logits = self.time_head(decoded)
        event_logit = self.event_head(gip).squeeze(-1)
        site_logits = self.site_head(gip)
        site_binary_logits = self.site_binary_head(gip) if self.site_binary_head is not None else None
        return logits, time_logits, event_logit, site_logits, site_binary_logits, attention


def train(args: argparse.Namespace) -> None:
    args.outdir.mkdir(parents=True, exist_ok=True)
    z = np.load(args.cache, allow_pickle=True)
    dataset = TrajDataset(args.cache)

    n = len(dataset)
    n_val = max(1, int(n * args.val_frac))
    n_train = n - n_val
    train_ds, val_ds = random_split(dataset, [n_train, n_val], generator=torch.Generator().manual_seed(args.seed))
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False)

    encoder = ProteinAwareGIPEncoder(
        n_snps=z["X"].shape[1],
        n_proteins=len(z["protein_ids"]),
        edge_snp_idx=z["edge_snp_idx"],
        edge_protein_idx=z["edge_protein_idx"],
        edge_beta=z["edge_beta"],
        edge_logp=z["edge_logp"],
        edge_cis=z["edge_cis"],
        edge_prior=z["edge_prior"] if "edge_prior" in z else None,
        edge_features=z["edge_features"] if args.use_edge_features and "edge_features" in z else None,
        d_model=args.d_model,
        dropout=args.dropout,
        prior_strength=args.prior_strength,
        protein_encoder_layers=args.protein_encoder_layers,
    )
    model = GIPTrajectoryDecoder(
        encoder=encoder,
        vocab_size=int(z["vocab_size"][0]),
        d_model=args.d_model,
        n_layers=args.decoder_layers,
        dropout=args.dropout,
        n_site_classes=len(z["site_names"]) + 1 if "site_names" in z else 1,
        n_site_binary=len(z["site_names"]) if "site_names" in z else 0,
        n_time_bins=len(z["time_bin_labels"]) if "time_bin_labels" in z else 1,
    )
    device = torch.device(args.device)
    model.to(device)

    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    ce = nn.CrossEntropyLoss(ignore_index=0)
    time_ce = nn.CrossEntropyLoss(ignore_index=0)
    bce = nn.BCEWithLogitsLoss()
    if args.label_weighting and "first_index_site_class" in z:
        site_counts = np.bincount(z["first_index_site_class"], minlength=len(z["site_names"]) + 1).astype(np.float32)
        if args.site_positive_only:
            positive_counts = site_counts.copy()
            positive_counts[0] = 0
            site_weights = positive_counts[1:].sum() / np.maximum(positive_counts[1:], 1.0)
            site_weights = site_weights / site_weights.mean()
        else:
            site_weights = site_counts.sum() / np.maximum(site_counts, 1.0)
            site_weights = site_weights / site_weights.mean()
        site_ce = nn.CrossEntropyLoss(weight=torch.tensor(site_weights, dtype=torch.float32, device=device))
    else:
        site_ce = nn.CrossEntropyLoss()
    if args.label_weighting and "site_binary" in z and z["site_binary"].shape[1] > 0:
        positives = z["site_binary"].sum(axis=0).astype(np.float32)
        negatives = z["site_binary"].shape[0] - positives
        site_pos_weight = torch.tensor(negatives / np.maximum(positives, 1.0), dtype=torch.float32, device=device)
        site_bce = nn.BCEWithLogitsLoss(pos_weight=site_pos_weight)
    else:
        site_bce = bce

    history = []
    history_path = args.outdir / "history.jsonl"
    history_path.write_text("")
    best_val_loss = math.inf
    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = total_lm = total_time = total_event = total_site = total_site_binary = 0.0
        for X, Y, time_bins, has_events, site_class, site_binary in train_loader:
            X, Y, has_events = X.to(device), Y.to(device), has_events.to(device)
            time_bins = time_bins.to(device)
            site_class, site_binary = site_class.to(device), site_binary.to(device)
            tokens_in = Y[:, :-1]
            targets = Y[:, 1:]
            time_targets = time_bins[:, 1:]
            logits, time_logits, event_logit, site_logits, site_binary_logits, _ = model(X, tokens_in)
            lm_loss = ce(logits.reshape(-1, logits.shape[-1]), targets.reshape(-1))
            time_loss = time_ce(time_logits.reshape(-1, time_logits.shape[-1]), time_targets.reshape(-1))
            event_loss = bce(event_logit, has_events)
            if args.site_positive_only:
                site_mask = site_class > 0
                site_loss = (
                    site_ce(site_logits[site_mask, 1:], site_class[site_mask] - 1)
                    if site_mask.any()
                    else site_logits.sum() * 0.0
                )
            else:
                site_loss = site_ce(site_logits, site_class)
            if site_binary_logits is not None and site_binary.shape[1] > 0:
                site_binary_loss = site_bce(site_binary_logits, site_binary)
            else:
                site_binary_loss = torch.tensor(0.0, device=device)
            loss = (
                lm_loss
                + args.time_loss_weight * time_loss
                + args.event_loss_weight * event_loss
                + args.site_loss_weight * site_loss
                + args.site_binary_loss_weight * site_binary_loss
            )
            opt.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), args.grad_clip)
            opt.step()
            total_loss += float(loss.detach())
            total_lm += float(lm_loss.detach())
            total_time += float(time_loss.detach())
            total_event += float(event_loss.detach())
            total_site += float(site_loss.detach())
            total_site_binary += float(site_binary_loss.detach())

        model.eval()
        val_loss = val_lm = val_time = val_event = val_site = val_site_binary = 0.0
        with torch.no_grad():
            for X, Y, time_bins, has_events, site_class, site_binary in val_loader:
                X, Y, has_events = X.to(device), Y.to(device), has_events.to(device)
                time_bins = time_bins.to(device)
                site_class, site_binary = site_class.to(device), site_binary.to(device)
                logits, time_logits, event_logit, site_logits, site_binary_logits, _ = model(X, Y[:, :-1])
                lm_loss = ce(logits.reshape(-1, logits.shape[-1]), Y[:, 1:].reshape(-1))
                time_loss = time_ce(time_logits.reshape(-1, time_logits.shape[-1]), time_bins[:, 1:].reshape(-1))
                event_loss = bce(event_logit, has_events)
                if args.site_positive_only:
                    site_mask = site_class > 0
                    site_loss = (
                        site_ce(site_logits[site_mask, 1:], site_class[site_mask] - 1)
                        if site_mask.any()
                        else site_logits.sum() * 0.0
                    )
                else:
                    site_loss = site_ce(site_logits, site_class)
                if site_binary_logits is not None and site_binary.shape[1] > 0:
                    site_binary_loss = site_bce(site_binary_logits, site_binary)
                else:
                    site_binary_loss = torch.tensor(0.0, device=device)
                loss = (
                    lm_loss
                    + args.time_loss_weight * time_loss
                    + args.event_loss_weight * event_loss
                    + args.site_loss_weight * site_loss
                    + args.site_binary_loss_weight * site_binary_loss
                )
                val_loss += float(loss)
                val_lm += float(lm_loss)
                val_time += float(time_loss)
                val_event += float(event_loss)
                val_site += float(site_loss)
                val_site_binary += float(site_binary_loss)

        row = {
            "epoch": epoch,
            "train_loss": total_loss / max(1, len(train_loader)),
            "train_lm_loss": total_lm / max(1, len(train_loader)),
            "train_time_loss": total_time / max(1, len(train_loader)),
            "train_event_loss": total_event / max(1, len(train_loader)),
            "train_site_loss": total_site / max(1, len(train_loader)),
            "train_site_binary_loss": total_site_binary / max(1, len(train_loader)),
            "val_loss": val_loss / max(1, len(val_loader)),
            "val_lm_loss": val_lm / max(1, len(val_loader)),
            "val_time_loss": val_time / max(1, len(val_loader)),
            "val_event_loss": val_event / max(1, len(val_loader)),
            "val_site_loss": val_site / max(1, len(val_loader)),
            "val_site_binary_loss": val_site_binary / max(1, len(val_loader)),
        }
        history.append(row)
        print(json.dumps(row, ensure_ascii=False), flush=True)
        with history_path.open("a") as out:
            out.write(json.dumps(row) + "\n")
        if row["val_loss"] < best_val_loss:
            best_val_loss = row["val_loss"]
            torch.save(
                {
                    "model_state_dict": model.state_dict(),
                    "args": vars(args),
                    "protein_ids": z["protein_ids"],
                    "snp_rsids": z["snp_rsids"],
                    "edge_feature_names": z["edge_feature_names"] if "edge_feature_names" in z else np.array([], dtype=object),
                    "site_names": z["site_names"] if "site_names" in z else np.array([], dtype=object),
                    "time_bin_labels": z["time_bin_labels"] if "time_bin_labels" in z else np.array([], dtype=object),
                    "best_epoch": epoch,
                    "best_val_loss": best_val_loss,
                },
                args.outdir / "model.best.pt",
            )

    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "args": vars(args),
            "protein_ids": z["protein_ids"],
            "snp_rsids": z["snp_rsids"],
            "edge_feature_names": z["edge_feature_names"] if "edge_feature_names" in z else np.array([], dtype=object),
            "site_names": z["site_names"] if "site_names" in z else np.array([], dtype=object),
            "time_bin_labels": z["time_bin_labels"] if "time_bin_labels" in z else np.array([], dtype=object),
        },
        args.outdir / "model.pt",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, default=Path("model_data/st10_core/st10_training_cache_smoke.npz"))
    parser.add_argument("--outdir", type=Path, default=Path("model_runs/st10_smoke"))
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--epochs", type=int, default=2)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--d-model", type=int, default=64)
    parser.add_argument("--decoder-layers", type=int, default=2)
    parser.add_argument("--protein-encoder-layers", type=int, default=1)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--lr", type=float, default=2e-4)
    parser.add_argument("--weight-decay", type=float, default=1e-4)
    parser.add_argument("--event-loss-weight", type=float, default=0.2)
    parser.add_argument("--time-loss-weight", type=float, default=0.5)
    parser.add_argument("--prior-strength", type=float, default=1.0)
    parser.add_argument("--no-edge-features", dest="use_edge_features", action="store_false")
    parser.set_defaults(use_edge_features=True)
    parser.add_argument("--site-loss-weight", type=float, default=0.5)
    parser.add_argument("--site-binary-loss-weight", type=float, default=0.5)
    parser.add_argument("--site-positive-only", action="store_true", default=True)
    parser.add_argument("--site-include-none", dest="site_positive_only", action="store_false")
    parser.add_argument("--no-label-weighting", dest="label_weighting", action="store_false")
    parser.set_defaults(label_weighting=True)
    parser.add_argument("--grad-clip", type=float, default=1.0)
    parser.add_argument("--val-frac", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=2026)
    train(parser.parse_args())


if __name__ == "__main__":
    main()
