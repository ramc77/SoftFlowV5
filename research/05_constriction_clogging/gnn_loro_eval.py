#!/usr/bin/env python3
"""Leave-one-run-out (LORO) evaluation of the Flavor-B clog GNN.

The pooled training metrics use a random split over *frames*, but frames from
the same simulation run are correlated, so that split leaks information. The
honest generalization test holds out an entire run (one aperture) as
validation and trains on the others — the GNN must classify a configuration
whose frames it never saw.

Each per-run dataset is read from ``gnn_data/<cell>/graphs.npz`` (written by
``gen_gnn_data.py``). Reports per-fold balanced accuracy + AUC and the mean.

Usage:
    python research/05_constriction_clogging/gnn_loro_eval.py
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
from pysoftflow.ml.gnn import fit_and_score  # noqa: E402

DATA_DIR = pathlib.Path(__file__).resolve().parent / "gnn_data"


def main() -> None:
    cells = sorted(p for p in DATA_DIR.glob("*/graphs.npz"))
    if len(cells) < 3:
        print(f"Need >= 3 per-run datasets in {DATA_DIR}; found {len(cells)}. "
              f"Run gen_gnn_data.py first.")
        sys.exit(1)

    # Load each run as a group; its label is the (single) clog flag of the run.
    groups = {}
    for npz in cells:
        gs = [g for g in load_graph_dataset(str(npz)) if g.graph_label is not None]
        if gs:
            groups[npz.parent.name] = gs
    labels = {name: gs[0].graph_label for name, gs in groups.items()}
    n_clog = sum(1 for v in labels.values() if v == 1.0)
    print(f"{len(groups)} runs ({n_clog} clog / {len(groups) - n_clog} no-clog):")
    for name, lab in labels.items():
        print(f"    {name:10s} label={'clog' if lab == 1.0 else 'no-clog':7s} "
              f"({len(groups[name])} graphs)")

    print("\nLeave-one-run-out folds (each held-out run is single-class, so the"
          "\nmetric is the fraction of that unseen run's frames classified right):")
    accs, correct_by_class = [], {0.0: [], 1.0: []}
    for held in groups:
        train = [g for name, gs in groups.items() if name != held for g in gs]
        if len({g.graph_label for g in train}) < 2:
            print(f"    hold {held:10s}: skipped (training would be single-class)")
            continue
        m = fit_and_score(train, groups[held], epochs=120, seed=0)
        accs.append(m["acc"])
        correct_by_class[labels[held]].append(m["acc"])
        tag = "clog" if labels[held] == 1.0 else "no-clog"
        print(f"    hold {held:10s} ({tag:7s}): {m['acc']:.3f} of frames correct")

    if accs:
        print(f"\nMean LORO accuracy: {np.mean(accs):.3f} over {len(accs)} folds")
        for lab, name in ((1.0, "clog"), (0.0, "no-clog")):
            if correct_by_class[lab]:
                print(f"    {name:7s} runs: mean {np.mean(correct_by_class[lab]):.3f}")
        print("\nThis is the honest generalization number (no frame leakage): can"
              "\nthe GNN label a configuration from a run it never trained on?")


if __name__ == "__main__":
    main()
