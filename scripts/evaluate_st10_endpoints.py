#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import torch
from sklearn.metrics import accuracy_score, average_precision_score, balanced_accuracy_score, roc_auc_score
from torch.utils.data import DataLoader, random_split


STAGE_ORDER = {
    "risk_or_precancer_broad": 0,
    "precancer_strict": 1,
    "index_cancer": 2,
    "metastasis": 3,
}


def install_numpy_pickle_compat() -> None:
    if "numpy._core" in sys.modules:
        return
    import numpy.core as numpy_core

    sys.modules["numpy._core"] = numpy_core
    for name in ("multiarray", "numeric", "umath"):
        module_name = f"numpy.core.{name}"
        if module_name in sys.modules:
            sys.modules[f"numpy._core.{name}"] = sys.modules[module_name]


def load_trainer(path: Path):
    spec = importlib.util.spec_from_file_location("st10_trainer", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load trainer module from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def split_token(token: str) -> tuple[str, str]:
    if "__" not in token:
        return token, ""
    site, stage = token.rsplit("__", 1)
    return site, stage


def relaxed_token_match(true_token: str, pred_token: str) -> bool:
    true_site, true_stage = split_token(true_token)
    pred_site, pred_stage = split_token(pred_token)
    if true_token == pred_token:
        return True
    if true_site != pred_site:
        return False
    if not true_stage or not pred_stage:
        return False
    if true_stage in STAGE_ORDER and pred_stage in STAGE_ORDER:
        return abs(STAGE_ORDER[true_stage] - STAGE_ORDER[pred_stage]) <= 1
    return False


def precancer_compatible_match(true_token: str, pred_token: str) -> bool:
    true_site, true_stage = split_token(true_token)
    pred_site, pred_stage = split_token(pred_token)
    if true_token == pred_token:
        return True
    if true_site != pred_site:
        return False
    precancer_stages = {"risk_or_precancer_broad", "precancer_strict"}
    return true_stage in precancer_stages or pred_stage in precancer_stages


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--trainer", type=Path, default=Path("scripts/train_st10_gip_trajectory_transformer.py"))
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--val-frac", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()

    trainer = load_trainer(args.trainer)
    z = np.load(args.cache, allow_pickle=True)
    dataset = trainer.TrajDataset(args.cache)
    n_val = max(1, int(len(dataset) * args.val_frac))
    n_train = len(dataset) - n_val
    _, val_ds = random_split(dataset, [n_train, n_val], generator=torch.Generator().manual_seed(args.seed))
    loader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False)

    install_numpy_pickle_compat()
    ckpt = torch.load(args.run_dir / "model.pt", map_location="cpu", weights_only=False)
    model_args = ckpt["args"]
    site_names = np.array(ckpt.get("site_names", z["site_names"] if "site_names" in z else []), dtype=object)
    encoder = trainer.ProteinAwareGIPEncoder(
        n_snps=z["X"].shape[1],
        n_proteins=len(z["protein_ids"]),
        edge_snp_idx=z["edge_snp_idx"],
        edge_protein_idx=z["edge_protein_idx"],
        edge_beta=z["edge_beta"],
        edge_logp=z["edge_logp"],
        edge_cis=z["edge_cis"],
        edge_prior=z["edge_prior"] if "edge_prior" in z else None,
        edge_features=z["edge_features"] if model_args.get("use_edge_features", False) and "edge_features" in z else None,
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
        n_time_bins=len(ckpt.get("time_bin_labels", z["time_bin_labels"] if "time_bin_labels" in z else [])) or 1,
    )
    model.load_state_dict(ckpt["model_state_dict"])
    model.to(args.device)
    model.eval()

    event_y, event_p, site_y, site_pred, site_positive_pred, site_binary_y, site_binary_p = [], [], [], [], [], [], []
    token_y, token_pred, token_top3, time_y, time_pred = [], [], [], [], []
    token_vocab = {i: str(t) for i, t in enumerate(z["token_names"])} if "token_names" in z else {}
    if not token_vocab:
        vocab_path = Path("model_data/traj/cancer_site_token_vocab.csv")
        if vocab_path.exists():
            import pandas as pd

            token_vocab = pd.read_csv(vocab_path).set_index("token_id")["token"].astype(str).to_dict()
    with torch.no_grad():
        for X, Y, time_bins, has_events, site_class, site_binary in loader:
            X, Y = X.to(args.device), Y.to(args.device)
            logits, time_logits, event_logit, site_logits, site_binary_logits, _ = model(X, Y[:, :-1])
            target_tokens = Y[:, 1:].cpu().numpy()
            pred_tokens = logits.argmax(dim=-1).cpu().numpy()
            top3_tokens = torch.topk(logits, k=min(3, logits.shape[-1]), dim=-1).indices.cpu().numpy()
            token_mask = target_tokens != 0
            token_y.extend(target_tokens[token_mask].tolist())
            token_pred.extend(pred_tokens[token_mask].tolist())
            token_top3.extend(top3_tokens[token_mask].tolist())
            target_times = time_bins[:, 1:].numpy()
            pred_times = time_logits.argmax(dim=-1).cpu().numpy()
            time_mask = target_times != 0
            time_y.extend(target_times[time_mask].tolist())
            time_pred.extend(pred_times[time_mask].tolist())
            event_y.extend(has_events.numpy().tolist())
            event_p.extend(torch.sigmoid(event_logit).cpu().numpy().tolist())
            site_y.extend(site_class.numpy().tolist())
            site_pred.extend(site_logits.argmax(dim=1).cpu().numpy().tolist())
            site_positive_pred.extend((site_logits[:, 1:].argmax(dim=1) + 1).cpu().numpy().tolist())
            if site_binary_logits is not None:
                site_binary_y.append(site_binary.numpy())
                site_binary_p.append(torch.sigmoid(site_binary_logits).cpu().numpy())

    event_y = np.array(event_y)
    event_p = np.array(event_p)
    site_y = np.array(site_y)
    site_pred = np.array(site_pred)
    site_positive_pred = np.array(site_positive_pred)
    print(f"val_n\t{len(event_y)}")
    print(f"event_rate\t{event_y.mean():.6f}")
    if token_y:
        print(f"next_token_n\t{len(token_y)}")
        print(f"next_token_accuracy\t{accuracy_score(token_y, token_pred):.6f}")
        top3_acc = np.mean([y in preds for y, preds in zip(token_y, token_top3)])
        print(f"next_token_top3_accuracy\t{top3_acc:.6f}")
        token_majority = Counter(token_y).most_common(1)[0][1] / len(token_y)
        print(f"next_token_majority_baseline\t{token_majority:.6f}")
        true_text = [token_vocab.get(int(t), str(t)) for t in token_y]
        pred_text = [token_vocab.get(int(t), str(t)) for t in token_pred]
        true_site_stage = [split_token(t) for t in true_text]
        pred_site_stage = [split_token(t) for t in pred_text]
        site_acc = np.mean([a[0] == b[0] for a, b in zip(true_site_stage, pred_site_stage)])
        stage_acc = np.mean([a[1] == b[1] for a, b in zip(true_site_stage, pred_site_stage)])
        relaxed_acc = np.mean([relaxed_token_match(a, b) for a, b in zip(true_text, pred_text)])
        precancer_compatible_acc = np.mean(
            [precancer_compatible_match(a, b) for a, b in zip(true_text, pred_text)]
        )
        same_site_stage_adjacent = np.mean(
            [
                a[0] == b[0]
                and a[1] in STAGE_ORDER
                and b[1] in STAGE_ORDER
                and abs(STAGE_ORDER[a[1]] - STAGE_ORDER[b[1]]) <= 1
                for a, b in zip(true_site_stage, pred_site_stage)
            ]
        )
        print(f"next_site_accuracy\t{site_acc:.6f}")
        print(f"next_stage_accuracy\t{stage_acc:.6f}")
        print(f"next_token_relaxed_same_site_stage_adjacent_accuracy\t{relaxed_acc:.6f}")
        print(f"next_token_relaxed_same_site_precancer_compatible_accuracy\t{precancer_compatible_acc:.6f}")
        print(f"next_same_site_stage_adjacent_rate\t{same_site_stage_adjacent:.6f}")
    if time_y:
        print(f"time_bin_n\t{len(time_y)}")
        print(f"time_bin_accuracy\t{accuracy_score(time_y, time_pred):.6f}")
        time_majority = Counter(time_y).most_common(1)[0][1] / len(time_y)
        print(f"time_bin_majority_baseline\t{time_majority:.6f}")
        print(f"time_bin_balanced_accuracy\t{balanced_accuracy_score(time_y, time_pred):.6f}")
    print(f"event_auc\t{roc_auc_score(event_y, event_p):.6f}")
    print(f"event_accuracy\t{accuracy_score(event_y, event_p >= 0.5):.6f}")
    print(f"site_class_accuracy\t{accuracy_score(site_y, site_pred):.6f}")
    print(f"site_class_balanced_accuracy\t{balanced_accuracy_score(site_y, site_pred):.6f}")
    print(f"site_class_nonzero_rate\t{(site_y > 0).mean():.6f}")
    positive_mask = site_y > 0
    if positive_mask.any():
        print(f"site_class_positive_n\t{int(positive_mask.sum())}")
        print(
            "site_class_positive_accuracy\t"
            f"{accuracy_score(site_y[positive_mask], site_positive_pred[positive_mask]):.6f}"
        )
        print(
            "site_class_positive_balanced_accuracy\t"
            f"{balanced_accuracy_score(site_y[positive_mask], site_positive_pred[positive_mask]):.6f}"
        )

    if site_binary_y:
        y = np.vstack(site_binary_y)
        p = np.vstack(site_binary_p)
        print("site\tpositives\tauc\taverage_precision")
        for i, site in enumerate(site_names):
            positives = int(y[:, i].sum())
            if positives == 0 or positives == len(y):
                print(f"{site}\t{positives}\tNA\tNA")
                continue
            print(
                f"{site}\t{positives}\t{roc_auc_score(y[:, i], p[:, i]):.6f}\t"
                f"{average_precision_score(y[:, i], p[:, i]):.6f}"
            )


if __name__ == "__main__":
    main()
