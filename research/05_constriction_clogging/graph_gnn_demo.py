#!/usr/bin/env python3
"""Train the Flavor-B GNN to classify suspension configuration graphs.

The GNN classifies clog / no-clog from a configuration graph
(nodes = capsules, edges = near-contacts), using the edge features
(surface gap, relative speed, approach rate) that carry the arch signal:

    snapshot -> ParticleGraph -> PyG Data -> edge-aware message-passing GNN
             -> clog / no-clog

This trainer works on REAL data only. It requires the labeled graph dataset
produced by ``gen_gnn_data.py`` (which runs the clogging simulation across
apertures and extracts per-frame contact graphs). If that dataset is missing
or unusable it prints a warning and exits non-zero — it never falls back to
synthetic data.

Usage:
    # 1. generate the real dataset from simulations:
    python research/05_constriction_clogging/gen_gnn_data.py
    # 2. train the GNN on it:
    python research/05_constriction_clogging/graph_gnn_demo.py
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

_REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "python"))
for _d in ("build", "build_phase2", "build_phase1"):
    _c = _REPO / _d / "python"
    if _c.is_dir():
        sys.path.insert(0, str(_c))

from pysoftflow.ml import load_graph_dataset  # noqa: E402
from pysoftflow.ml.gnn import train_clog_gnn  # noqa: E402

REAL_DATASET = pathlib.Path(__file__).resolve().parent / "gnn_data" / "gnn_dataset.npz"

MIN_GRAPHS = 8  # need a meaningful train/val split


def _warn_and_exit(msg: str) -> None:
    print("\n" + "=" * 70)
    print("WARNING: cannot train the GNN — no usable REAL data.")
    print(msg)
    print("Generate it first:")
    print("    python research/05_constriction_clogging/gen_gnn_data.py")
    print("=" * 70)
    sys.exit(1)


def main() -> None:
    if not REAL_DATASET.exists():
        _warn_and_exit(f"Dataset not found: {REAL_DATASET}")

    graphs = [g for g in load_graph_dataset(str(REAL_DATASET)) if g.graph_label is not None]
    labels = {g.graph_label for g in graphs}

    if len(graphs) < MIN_GRAPHS:
        _warn_and_exit(f"Only {len(graphs)} labeled graphs in {REAL_DATASET} "
                       f"(need >= {MIN_GRAPHS}).")
    if len(labels) < 2:
        _warn_and_exit(
            f"Dataset has a single class {labels} — need both clog and no-clog "
            f"examples. Widen the aperture range in gen_gnn_data.py so some "
            f"runs clog and some pass.")

    n_clog = sum(1 for g in graphs if g.graph_label == 1.0)
    contacts = np.array([g.n_edges // 2 for g in graphs])
    print(f"REAL dataset: {REAL_DATASET}")
    print(f"  {len(graphs)} graphs ({n_clog} clog / {len(graphs) - n_clog} no-clog)")
    print(f"  contacts/graph: mean={contacts.mean():.1f}, max={int(contacts.max())}")
    print(f"  node features: {graphs[0].node_feature_names}")
    print(f"  edge features: {graphs[0].edge_feature_names}\n")

    print("Training edge-aware GNN (Flavor B) on real simulation data...")
    res = train_clog_gnn(graphs, epochs=120, hidden=48, n_layers=3, seed=0)

    print("\n  epoch   train_loss   val_acc")
    for ep, tl, va in res.history:
        print(f"  {ep:5d}   {tl:9.4f}   {va:6.3f}")
    print(f"\nFinal: train_acc={res.train_acc:.3f}  val_acc={res.val_acc:.3f}  "
          f"val_balanced_acc={res.val_balanced_acc:.3f}  val_AUC={res.val_auc:.3f}")
    print("(balanced_acc and AUC are the meaningful metrics under class imbalance)")


if __name__ == "__main__":
    main()
