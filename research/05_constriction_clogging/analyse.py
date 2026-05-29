"""
Project 5 analysis — clog detection, jamming labels, and GNN graph dataset.
===========================================================================
Reads a ``particle_data.csv`` trajectory and computes:

  • First-passage escapee curve N(t): capsules that have crossed the neck.
  • Clog flag: True if the suspension stops draining with capsules still
    trapped upstream (the Bielinski 2021 "clog" state).
  • Evacuation time: when 90 % of capsules have passed (NaN if clogged).
  • Neck contact number Z: mean contact degree of capsules in the neck band
    (the arch signature), via pysoftflow.ml.snapshot_to_graph.

It then builds a GNN graph dataset (one labelled graph per frame, label =
run-level clog flag) with pysoftflow.ml.build_graph_dataset, ready for the
Flavor-B graph neural network.

Run on a real trajectory:
    python research/05_constriction_clogging/analyse.py \
        --csv vtk_clog/particle_data.csv --neck-x 185

Verify the logic without the engine (synthetic clog + pass trajectories):
    python research/05_constriction_clogging/analyse.py --self-test
"""

from __future__ import annotations

import argparse
import pathlib
import sys
from dataclasses import dataclass

import numpy as np

script_dir = pathlib.Path(__file__).resolve().parent
project_dir = script_dir.parents[1]
sys.path.insert(0, str(project_dir / "python"))

from pysoftflow.ml import (  # noqa: E402
    build_graph_dataset,
    load_particle_csv,
    save_graph_dataset,
    snapshot_to_graph,
)

# Geometry defaults (match run.py).
NX, NY = 220.0, 48.0
R_CAPSULE = 4.0
NECK_BAND = 15.0   # half-width of the x-band around the neck used for Z
PERIODIC_X = True   # run.py uses a periodic channel + body force


@dataclass
class ClogMetrics:
    """Summary of one trajectory's clogging behaviour."""

    n_cells: int
    escaped_fraction: float
    clogged: bool
    evacuation_time: float
    neck_contact_number: float
    n_frames: int


def clog_metrics(
    snaps,
    neck_x: float,
    *,
    escape_threshold: float = 0.9,
    stall_window: int = 5,
) -> ClogMetrics:
    """Compute clog metrics from a list of SimulationSnapshots.

    Parameters
    ----------
    snaps:
        Frames in time order (from :func:`load_particle_csv`).
    neck_x:
        Streamwise position of the constriction centre.
    escape_threshold:
        Fraction of capsules that must pass for the run to count as drained.
    stall_window:
        A run is "clogged" if no new capsule escapes over the last
        ``stall_window`` frames while capsules remain trapped upstream.

    Returns
    -------
    ClogMetrics
    """
    if not snaps:
        return ClogMetrics(0, 0.0, False, float("nan"), 0.0, 0)

    n_cells = snaps[0].n_particles
    escaped_ids: set[int] = set()
    escaped_count = np.zeros(len(snaps), dtype=int)
    evac_time = float("nan")

    for f, snap in enumerate(snaps):
        xs = snap.positions[:, 0]
        for i in range(snap.n_particles):
            if xs[i] > neck_x:
                escaped_ids.add(i)
        escaped_count[f] = len(escaped_ids)
        if np.isnan(evac_time) and len(escaped_ids) >= escape_threshold * n_cells:
            evac_time = float(snap.time)

    escaped_fraction = escaped_count[-1] / n_cells if n_cells else 0.0
    # Clogged: drained below threshold AND flux stalled at the end.
    recent = escaped_count[-stall_window:] if len(escaped_count) >= stall_window else escaped_count
    flux_stalled = recent[-1] == recent[0]
    clogged = (escaped_fraction < escape_threshold) and flux_stalled and (escaped_fraction < 1.0)

    z_neck = _neck_contact_number(snaps, neck_x)
    return ClogMetrics(
        n_cells=n_cells,
        escaped_fraction=float(escaped_fraction),
        clogged=bool(clogged),
        evacuation_time=evac_time,
        neck_contact_number=z_neck,
        n_frames=len(snaps),
    )


def _neck_contact_number(snaps, neck_x: float) -> float:
    """Mean contact degree of capsules within the neck band, over all frames."""
    degs: list[float] = []
    for snap in snaps:
        g = snapshot_to_graph(snap, contact_cutoff=1.0)
        deg_col = g.node_feature_names.index("contact_degree")
        xnorm_col = g.node_feature_names.index("x_norm")
        x = g.node_features[:, xnorm_col] * snap.Lx
        in_band = np.abs(x - neck_x) <= NECK_BAND
        if np.any(in_band):
            degs.append(float(np.mean(g.node_features[in_band, deg_col])))
    return float(np.mean(degs)) if degs else 0.0


def analyse_trajectory(csv_path: str, neck_x: float, *, save_graphs: str | None = None) -> ClogMetrics:
    """Load a trajectory CSV, compute clog metrics, optionally save a GNN dataset."""
    snaps = load_particle_csv(
        csv_path, radius=R_CAPSULE, box=(NX, NY), periodic_x=PERIODIC_X, frame_stride=1
    )
    metrics = clog_metrics(snaps, neck_x)
    if save_graphs:
        graphs = build_graph_dataset(
            csv_path, radius=R_CAPSULE, box=(NX, NY), periodic_x=PERIODIC_X,
            contact_cutoff=1.0, frame_stride=5,
            label=float(metrics.clogged),
        )
        save_graph_dataset(graphs, save_graphs)
    return metrics


# =================================================================
# Self-test: synthesise a clogging run and a passing run (no engine)
# =================================================================

def _write_synthetic_csv(path: pathlib.Path, *, clog: bool, n=12, frames=40) -> None:
    """Write a synthetic particle_data.csv for one of two regimes."""
    neck_x = 185.0
    rng = np.random.default_rng(0)
    y0 = np.linspace(10, NY - 10, n)
    rows = []
    for f in range(frames):
        t = f * 200.0
        for cid in range(n):
            if clog:
                # Capsules crowd up to the neck and stop (arch).
                x = min(40 + 3.0 * f, neck_x - 6.0 - 0.3 * cid)
            else:
                # Capsules stream steadily through the neck.
                x = 40 + 6.0 * f + 4.0 * cid
            x = float(np.clip(x, 5, NX - 5))
            y = float(np.clip(y0[cid] + rng.normal(0, 0.2), 2, NY - 2))
            rows.append((f * 200, t, cid, 0, x, y, 0.001, 0.0))
    header = "timestep,time,capsule_id,type,cx,cy,vx,vy"
    arr = "\n".join(f"{a},{b:.6e},{c},{d},{e:.6e},{g:.6e},{h:.6e},{i:.6e}"
                    for (a, b, c, d, e, g, h, i) in rows)
    path.write_text(header + "\n" + arr + "\n")


def self_test() -> int:
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        clog_csv = tmp / "clog.csv"
        pass_csv = tmp / "pass.csv"
        _write_synthetic_csv(clog_csv, clog=True)
        _write_synthetic_csv(pass_csv, clog=False)

        m_clog = analyse_trajectory(str(clog_csv), neck_x=185.0)
        m_pass = analyse_trajectory(str(pass_csv), neck_x=185.0)

    print("Self-test (synthetic trajectories, no engine required)")
    print(f"  CLOG  run : escaped={m_clog.escaped_fraction:.2f}  clogged={m_clog.clogged}  "
          f"Z_neck={m_clog.neck_contact_number:.2f}  evac_t={m_clog.evacuation_time}")
    print(f"  PASS  run : escaped={m_pass.escaped_fraction:.2f}  clogged={m_pass.clogged}  "
          f"Z_neck={m_pass.neck_contact_number:.2f}  evac_t={m_pass.evacuation_time}")

    ok = m_clog.clogged and (not m_pass.clogged) and (m_pass.escaped_fraction > 0.9)
    print("  RESULT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", type=str, help="Path to particle_data.csv")
    p.add_argument("--neck-x", type=float, default=162.5, help="Neck centre x (LU)")
    p.add_argument("--save-graphs", type=str, default=None,
                   help="Optional path to save the GNN graph dataset (.npz)")
    p.add_argument("--self-test", action="store_true",
                   help="Run on synthetic trajectories (no engine needed)")
    args = p.parse_args(argv if argv is not None else sys.argv[1:])

    if args.self_test:
        return self_test()
    if not args.csv:
        p.error("provide --csv PATH (or use --self-test)")

    m = analyse_trajectory(args.csv, args.neck_x, save_graphs=args.save_graphs)
    print(f"n_cells            = {m.n_cells}")
    print(f"escaped_fraction   = {m.escaped_fraction:.3f}")
    print(f"clogged            = {m.clogged}")
    print(f"evacuation_time    = {m.evacuation_time}")
    print(f"neck_contact_number= {m.neck_contact_number:.3f}")
    print(f"n_frames           = {m.n_frames}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
