#!/usr/bin/env python3
"""End-to-end results + publication figures for the clogging / GNN study.

Consumes the REAL data produced by ``gen_gnn_data.py`` (per-aperture
trajectories + pooled graph dataset) and produces:

  fig1_phase.png        clog physics: escaped fraction & neck contact number
                        Z vs aperture D/d, with the clog/no-clog transition.
  fig2_gnn_training.png pooled GNN training curve + ROC (AUC) on a held-out
                        frame split.
  fig3_loro.png         leave-one-run-out per-fold accuracy (honest
                        generalization), coloured by class, with the mean.
  fig4_graphs.png       what the GNN sees: a clog vs a no-clog configuration
                        graph (capsule nodes + contact edges).

Also prints a results summary. Run AFTER gen_gnn_data.py.

Usage:
    python research/05_constriction_clogging/make_figures.py
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

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pysoftflow.ml import load_graph_dataset  # noqa: E402
from pysoftflow.ml.gnn import fit_and_score, graphs_to_pyg, train_clog_gnn  # noqa: E402
import torch  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "gnn_data"
FIGS = HERE / "figures"
NX, NY = 220.0, 48.0


# --------------------------------------------------------------------------
# Load per-run groups (each cell dir = one aperture run)
# --------------------------------------------------------------------------

def load_runs() -> dict[str, dict]:
    """Return {cell_name: {aperture, label, graphs}} from gnn_data/*/graphs.npz."""
    runs = {}
    for npz in sorted(DATA.glob("*/graphs.npz")):
        gs = [g for g in load_graph_dataset(str(npz)) if g.graph_label is not None]
        if not gs:
            continue
        # cell name like ap2p75_s0 -> aperture 2.75
        ap_str = npz.parent.name.split("_")[0].replace("ap", "").replace("p", ".")
        runs[npz.parent.name] = {
            "aperture": float(ap_str),
            "label": gs[0].graph_label,
            "graphs": gs,
        }
    return runs


def neck_contacts(graphs, neck_x=162.5, band=15.0) -> float:
    """Mean contact degree of nodes within the neck band, over the run's frames."""
    vals = []
    for g in graphs:
        xi = g.node_feature_names.index("x_norm")
        di = g.node_feature_names.index("contact_degree")
        x = g.node_features[:, xi] * NX
        m = np.abs(x - neck_x) <= band
        if m.any():
            vals.append(float(g.node_features[m, di].mean()))
    return float(np.mean(vals)) if vals else 0.0


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

def fig_phase(runs: dict, path: pathlib.Path) -> None:
    """Z_neck and clog label vs aperture — the clog->pass transition."""
    aps = sorted({r["aperture"] for r in runs.values()})
    z_by_ap, clog_by_ap = {}, {}
    for ap in aps:
        rs = [r for r in runs.values() if r["aperture"] == ap]
        z_by_ap[ap] = np.mean([neck_contacts(r["graphs"]) for r in rs])
        clog_by_ap[ap] = np.mean([r["label"] for r in rs])  # fraction of seeds that clog

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(aps, [z_by_ap[a] for a in aps], "o-", color="#1f77b4", label="neck contact number $Z$")
    ax.set_xlabel("aperture ratio  $D/d$")
    ax.set_ylabel("neck contact number  $Z$", color="#1f77b4")
    ax.tick_params(axis="y", labelcolor="#1f77b4")
    ax2 = ax.twinx()
    ax2.step(aps, [clog_by_ap[a] for a in aps], where="mid", color="#d62728",
             label="clog fraction", lw=2)
    ax2.set_ylabel("clog fraction (per-seed)", color="#d62728")
    ax2.tick_params(axis="y", labelcolor="#d62728")
    ax2.set_ylim(-0.05, 1.05)
    # Mark the transition.
    ax.axvspan(2.75, 3.0, color="grey", alpha=0.15)
    ax.text(2.875, ax.get_ylim()[1] * 0.95, "clog→pass", ha="center",
            va="top", fontsize=9, color="grey")
    ax.set_title("Clogging transition vs constriction aperture (LBM–IBM simulations)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_gnn_training(graphs, path: pathlib.Path) -> dict:
    """Pooled training curve + ROC. Returns metrics dict."""
    res = train_clog_gnn(graphs, epochs=120, hidden=48, n_layers=3, seed=0)
    # ROC on a held-out split (reuse train_clog_gnn's val via fresh scores).
    from pysoftflow.ml.gnn import _scores  # noqa
    from torch_geometric.loader import DataLoader
    data = graphs_to_pyg(graphs)
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(data))
    n_val = max(1, int(0.25 * len(data)))
    val = [data[i] for i in idx[:n_val]]
    y, s = _scores(res.model, DataLoader(val, batch_size=16))

    # ROC curve (no sklearn).
    order = np.argsort(-s)
    y_sorted = y[order]
    P, N = (y == 1).sum(), (y == 0).sum()
    tpr = np.cumsum(y_sorted == 1) / max(P, 1)
    fpr = np.cumsum(y_sorted == 0) / max(N, 1)
    tpr = np.concatenate([[0], tpr]); fpr = np.concatenate([[0], fpr])

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
    epochs = [h[0] for h in res.history]
    a1.plot(epochs, [h[1] for h in res.history], "o-", label="train loss")
    a1.plot(epochs, [h[2] for h in res.history], "s-", label="val accuracy")
    a1.set_xlabel("epoch"); a1.set_ylabel("value"); a1.legend()
    a1.set_title("GNN training (pooled frames)")
    a2.plot(fpr, tpr, color="#2ca02c", lw=2, label=f"ROC (AUC={res.val_auc:.2f})")
    a2.plot([0, 1], [0, 1], "k--", lw=1)
    a2.set_xlabel("false positive rate"); a2.set_ylabel("true positive rate")
    a2.legend(); a2.set_title("Clog classifier ROC (held-out frames)")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return {"train_acc": res.train_acc, "val_acc": res.val_acc,
            "val_balanced_acc": res.val_balanced_acc, "val_auc": res.val_auc}


def fig_loro(runs: dict, path: pathlib.Path) -> float:
    """Leave-one-run-out per-fold accuracy bar chart. Returns mean."""
    names = list(runs)
    accs, colors, labels = [], [], []
    for held in names:
        train = [g for n in names if n != held for g in runs[n]["graphs"]]
        if len({g.graph_label for g in train}) < 2:
            continue
        m = fit_and_score(train, runs[held]["graphs"], epochs=120, seed=0)
        accs.append(m["acc"])
        colors.append("#d62728" if runs[held]["label"] == 1.0 else "#1f77b4")
        labels.append(f"{runs[held]['aperture']:.2f}")
    mean = float(np.mean(accs)) if accs else float("nan")

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.bar(range(len(accs)), accs, color=colors)
    ax.axhline(mean, color="k", ls="--", label=f"mean = {mean:.3f}")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("fraction of unseen-run frames classified correctly")
    ax.set_xlabel("held-out run (aperture $D/d$)")
    ax.set_ylim(0, 1.05)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#d62728", label="clog run"),
                       Patch(color="#1f77b4", label="no-clog run"),
                       plt.Line2D([], [], color="k", ls="--", label=f"mean={mean:.3f}")])
    ax.set_title("Leave-one-run-out generalization of the clog GNN")
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)
    return mean


def fig_graphs(runs: dict, path: pathlib.Path) -> None:
    """Streamwise contact-degree profile: where contacts concentrate.

    A whole-domain contact graph looks alike for clog and no-clog (both are
    dense suspensions); the discriminating quantity is *where* the contacts sit
    relative to the neck. We bin each capsule's contact degree by its streamwise
    position and average over all frames of the clogging vs the flowing runs.
    The clogging ensemble shows a pronounced accumulation just upstream of the
    neck (the arch); the flowing ensemble is flatter and shifted downstream.
    """
    bins = np.linspace(0, NX, 23)
    centres = 0.5 * (bins[:-1] + bins[1:])

    def profile(graphs):
        deg_sum = np.zeros(len(centres))
        cnt = np.zeros(len(centres))
        for g in graphs:
            xi = g.node_feature_names.index("x_norm")
            di = g.node_feature_names.index("contact_degree")
            x = g.node_features[:, xi] * NX
            d = g.node_features[:, di]
            idx = np.clip(np.digitize(x, bins) - 1, 0, len(centres) - 1)
            for b, dd in zip(idx, d):
                deg_sum[b] += dd
                cnt[b] += 1
        with np.errstate(invalid="ignore"):
            return np.where(cnt > 0, deg_sum / cnt, np.nan)

    clog_graphs = [g for r in runs.values() if r["label"] == 1.0 for g in r["graphs"]]
    flow_graphs = [g for r in runs.values() if r["label"] == 0.0 for g in r["graphs"]]

    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(centres, profile(clog_graphs), "o-", color="#d62728",
            label=r"clogging runs ($D/d \leq 2.75$)")
    ax.plot(centres, profile(flow_graphs), "s-", color="#1f77b4",
            label=r"flowing runs ($D/d \geq 3.0$)")
    ax.axvspan(150, 175, color="grey", alpha=0.2, label="constriction")
    ax.set_xlabel("streamwise position  $x$ (LU)")
    ax.set_ylabel("mean contact degree per capsule")
    ax.set_title("Where contacts concentrate: clogging vs flowing ensembles")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main() -> None:
    if not (DATA / "gnn_dataset.npz").exists():
        print(f"No dataset at {DATA}. Run gen_gnn_data.py first.")
        sys.exit(1)
    FIGS.mkdir(exist_ok=True)
    torch.manual_seed(0)

    runs = load_runs()
    pooled = [g for r in runs.values() for g in r["graphs"]]
    n_clog = sum(1 for g in pooled if g.graph_label == 1.0)
    print(f"Loaded {len(runs)} runs, {len(pooled)} graphs "
          f"({n_clog} clog / {len(pooled) - n_clog} no-clog)\n")

    print("[1/4] phase diagram ...")
    fig_phase(runs, FIGS / "fig1_phase.png")
    print("[2/4] GNN training + ROC ...")
    m = fig_gnn_training(pooled, FIGS / "fig2_gnn_training.png")
    print("[3/4] leave-one-run-out ...")
    loro = fig_loro(runs, FIGS / "fig3_loro.png")
    print("[4/4] configuration graphs ...")
    fig_graphs(runs, FIGS / "fig4_graphs.png")

    print("\n================ RESULTS SUMMARY ================")
    print(f"Runs / graphs        : {len(runs)} runs, {len(pooled)} graphs")
    print(f"GNN pooled val acc   : {m['val_acc']:.3f}")
    print(f"GNN balanced acc     : {m['val_balanced_acc']:.3f}")
    print(f"GNN ROC AUC          : {m['val_auc']:.3f}")
    print(f"LORO generalization  : {loro:.3f}  (honest, no frame leakage)")
    print(f"Figures              : {FIGS}/fig1..4.png")
    print("=================================================")


if __name__ == "__main__":
    main()
