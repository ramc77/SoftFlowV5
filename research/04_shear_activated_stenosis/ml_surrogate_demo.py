#!/usr/bin/env python3
"""End-to-end demo of the reduced-order ML layer on the stenosis sweep.

Runs three things on the existing 4x4 ``sweep_results.npz``, using only
NumPy/SciPy (no PyTorch / scikit-learn):

    (b) Reduced-order surrogate  — fit a Gaussian process to
        (severity, log10 gamma_th) -> eta_deposit, with leave-one-out CV.
    (c) Explainable ML           — ARD feature-importance for the two inputs.
    (a) Inverse design           — find the (severity, gamma_th) that maximises
        predicted deposition; suggest the next most-informative simulation.

Usage
-----
    python research/04_shear_activated_stenosis/ml_surrogate_demo.py
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

# Make pysoftflow importable without an editable install.
_REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "python"))

from pysoftflow.ml import (  # noqa: E402
    GaussianProcessSurrogate,
    load_sweep,
    optimise_target,
    suggest_next,
)

NPZ = pathlib.Path(__file__).resolve().parent / "sweep_results.npz"


def main() -> None:
    print(f"Loading sweep: {NPZ}")
    ds = load_sweep(NPZ, features=("severity", "gamma_th"), target="eta_deposit")
    print(f"  {ds.n_samples} cells, features={ds.feature_names}, target={ds.target_name}\n")

    gp = GaussianProcessSurrogate(n_restarts=12, random_state=0).fit(ds.X, ds.y)
    print("(b) Surrogate fit")
    print(f"    log marginal likelihood : {gp.log_marginal_likelihood:8.3f}")
    print(f"    ARD length-scales        : {np.array2string(gp.lengthscales, precision=3)}")
    print(f"    signal var / noise var   : {gp.signal_var:.3f} / {gp.noise_var:.4f}")

    loo = gp.leave_one_out()
    max_z = np.nanmax(np.abs(loo.standardised_residuals))
    print("\n    Leave-one-out cross-validation")
    print(f"      RMSE = {loo.rmse:.4f}   MAE = {loo.mae:.4f}   R^2 = {loo.r2:.3f}")
    print(f"      max |standardised residual| = {max_z:.2f}")

    print("\n(c) Feature importance (ARD relevance, sums to 1)")
    for name, rel in gp.feature_importance(ds.feature_names).items():
        bar = "#" * int(round(rel * 40))
        print(f"      {name:>16s} : {rel:5.3f} {bar}")

    print("\n(a) Inverse design — maximise predicted eta_deposit")
    res = optimise_target(gp, ds.bounds(), maximise=True)
    phys = ds.invert_features(res.x)
    print(f"      optimum (feature space) : {np.array2string(res.x, precision=3)}")
    print(f"      optimum (physical)      : severity = {phys[0]:.1f} %, "
          f"gamma_th = {phys[1]:.2e} LU")
    print(f"      predicted eta_deposit   : {res.predicted_value:.3f} "
          f"+/- {res.predicted_std:.3f}")
    print("      (compare: measured sweet-spot peak eta_deposit = 0.772 "
          "at sigma=25%, gamma_th=2e-3)")

    print("\n      Next simulation to run (active learning)")
    nxt = suggest_next(gp, ds.bounds(), mode="explore")
    nxt_phys = ds.invert_features(nxt["point"])
    print(f"        {nxt['mode']}")
    print(f"        suggest severity = {nxt_phys[0]:.1f} %, "
          f"gamma_th = {nxt_phys[1]:.2e} LU  (pred std = {nxt['score']:.3f})")


if __name__ == "__main__":
    main()
