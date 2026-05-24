"""
Parameter sweep over (severity, γ̇_th) for Project 4
====================================================
4 × 4 = 16 cells. Each cell runs the shear-activated stenosis
simulation at one (severity, threshold) pair and records the
deposition efficiency, off-target loss, release fraction, and
delivery timing.

Resume-safe (cell_result.json cache). --force overrides cache.

Usage
-----
    python sweep.py --smoke           # 2×2 × 5000 steps (~10 min)
    python sweep.py                    # 4×4 × 30000 steps (~1.5 h)
    python sweep.py --force            # ignore cache
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import pathlib
import sys
import time
from dataclasses import dataclass
from typing import Optional

import numpy as np


# -- Path setup --
script_dir   = pathlib.Path(__file__).resolve().parent
project_dir  = script_dir.parents[1]
sys.path.insert(0, str(project_dir / "python"))
for d in ("build", "build_phase1", "build_phase2"):
    cand = project_dir / d / "python"
    if cand.is_dir():
        sys.path.insert(0, str(cand))

from run import (                                                      # noqa: E402
    build_sim, save_history,
    N_CARRIER, PAYLOAD_PER_CARRIER, OUT_EVERY,
)


# =================================================================
# Sweep grid
# =================================================================

SEVERITIES_PROD  = (0, 25, 50, 75)
SEVERITIES_SMOKE = (0, 50)

GAMMA_TH_PROD  = (5e-4, 1e-3, 2e-3, 5e-3)
GAMMA_TH_SMOKE = (1e-3, 2e-3)


# =================================================================
# Cell result
# =================================================================

CELL_RESULT_FILE = "cell_result.json"


@dataclass
class CellResult:
    """Per-cell summary. All metrics are monotone-growing time
    integrals — no snapshot fragility."""
    severity:           int           # %
    gamma_th:           float
    eta_deposit:        float         # target_absorbed / released
    eta_offtarget:      float         # off_target_absorbed / released
    eta_remaining:      float         # 1 - eta_deposit - eta_offtarget
    release_fraction:   float         # released / loaded
    target_absorbed:    float
    off_target_absorbed: float
    payload_released:   float
    payload_remaining:  float
    t_50_deposit:       float         # step at which target_absorbed reaches 50% of final
    walltime_s:         float


# =================================================================
# Cache helpers
# =================================================================

_LAST_N_STEPS_USED: int = 0


def _save_cell_result(result: CellResult, cell_dir: pathlib.Path) -> None:
    payload = dataclasses.asdict(result)
    for k, v in list(payload.items()):
        if isinstance(v, float) and (v != v):       # NaN
            payload[k] = None
    payload["_format_version"] = 1
    payload["_n_steps"] = int(_LAST_N_STEPS_USED)
    (cell_dir / CELL_RESULT_FILE).write_text(json.dumps(payload, indent=2))


def _load_cell_result(cell_dir: pathlib.Path,
                       expected_n_steps: int) -> Optional[CellResult]:
    f = cell_dir / CELL_RESULT_FILE
    if not f.exists():
        return None
    try:
        data = json.loads(f.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    if not isinstance(data, dict) or data.get("_format_version") != 1:
        return None
    if data.get("_n_steps") != expected_n_steps:
        return None
    data.pop("_format_version", None)
    data.pop("_n_steps", None)
    for k, v in list(data.items()):
        if v is None:
            data[k] = float("nan")
    try:
        return CellResult(**data)
    except (TypeError, ValueError):
        return None


# =================================================================
# Reduce simulation → CellResult
# =================================================================

def _reduce_to_result(*, severity: int, gamma_th: float,
                       run, n_steps: int,
                       walltime_s: float) -> CellResult:
    nan = float("nan")
    if not run.history:
        return CellResult(
            severity=severity, gamma_th=gamma_th,
            eta_deposit=nan, eta_offtarget=nan, eta_remaining=nan,
            release_fraction=0.0,
            target_absorbed=0.0, off_target_absorbed=0.0,
            payload_released=0.0,
            payload_remaining=float(N_CARRIER * PAYLOAD_PER_CARRIER),
            t_50_deposit=nan,
            walltime_s=walltime_s,
        )

    last     = run.history[-1]
    released = float(last.total_released)
    target   = float(last.target_absorbed)
    off      = float(last.off_target_absorbed)
    remain   = float(last.total_M_p)

    if released > 1e-12:
        eta_d = target / released
        eta_o = off / released
        eta_r = max(0.0, 1.0 - eta_d - eta_o)
    else:
        eta_d = eta_o = eta_r = nan

    release_fraction = released / float(N_CARRIER * PAYLOAD_PER_CARRIER)

    # t_50_deposit: first step at which target_absorbed reaches 50% of final
    if target > 1e-12:
        half = 0.5 * target
        t_50 = next(
            (float(r.step) for r in run.history if r.target_absorbed >= half),
            nan)
    else:
        t_50 = nan

    return CellResult(
        severity=severity, gamma_th=gamma_th,
        eta_deposit=eta_d, eta_offtarget=eta_o, eta_remaining=eta_r,
        release_fraction=release_fraction,
        target_absorbed=target, off_target_absorbed=off,
        payload_released=released, payload_remaining=remain,
        t_50_deposit=t_50,
        walltime_s=walltime_s,
    )


# =================================================================
# One simulation cell
# =================================================================

def run_one_cell(*, severity: int, gamma_th: float,
                  output_dir: str, n_steps: int) -> CellResult:
    t_start = time.time()
    sim, run = build_sim(severity=severity, gamma_th=gamma_th,
                          output_dir=os.path.join(output_dir, "vtk_stenosis"))

    sim.thermo(interval=max(n_steps // 10, 200))
    sim.warmup(steps=500, ramp_steps=1500)
    sim.run(n_steps)

    save_history(run, out_dir=output_dir,
                  severity=severity, gamma_th=gamma_th)
    result = _reduce_to_result(severity=severity, gamma_th=gamma_th,
                                run=run, n_steps=n_steps,
                                walltime_s=time.time() - t_start)
    _save_cell_result(result, pathlib.Path(output_dir))
    return result


# =================================================================
# Sweep driver
# =================================================================

def _grid_pairs(severities, thresholds):
    return [(s, t) for s in severities for t in thresholds]


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--smoke", action="store_true",
                    help="2×2 grid × 5000 steps (~10 min)")
    p.add_argument("--steps", type=int, default=None,
                    help="override per-cell N_STEPS (smoke=5k, prod=30k)")
    p.add_argument("--force", action="store_true",
                    help="ignore cached cell_result.json files")
    return p.parse_args(argv)


def _cell_dirname(severity: int, gamma_th: float, thresholds) -> str:
    """Deterministic short directory name for a (severity, γ̇_th) cell."""
    th_idx = thresholds.index(gamma_th) + 1
    return f"sev{severity:02d}_th{th_idx}"


def main(argv=None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    if args.smoke:
        sev_vals, th_vals = SEVERITIES_SMOKE, GAMMA_TH_SMOKE
        n_steps_default = 5_000
    else:
        sev_vals, th_vals = SEVERITIES_PROD, GAMMA_TH_PROD
        n_steps_default = 200_000
    n_steps = args.steps if args.steps is not None else n_steps_default

    pairs = _grid_pairs(sev_vals, th_vals)
    sweep_root = script_dir / "sweep_out"
    sweep_root.mkdir(exist_ok=True)

    global _LAST_N_STEPS_USED
    _LAST_N_STEPS_USED = int(n_steps)

    cached: list[tuple[int, CellResult]] = []
    to_run: list[tuple[int, int, float, pathlib.Path]] = []
    for k, (sev, th) in enumerate(pairs):
        cell_dir = sweep_root / _cell_dirname(sev, th, th_vals)
        cell_dir.mkdir(exist_ok=True)
        if args.force:
            cached_r = None
        else:
            cached_r = _load_cell_result(cell_dir, expected_n_steps=n_steps)
        if cached_r is not None:
            cached.append((k, cached_r))
        else:
            to_run.append((k, sev, th, cell_dir))

    print(f"\nSweep: {len(pairs)} cells, n_steps={n_steps}")
    print(f"  severities:  {sev_vals} %")
    print(f"  γ̇_th:        {th_vals}")
    print(f"  resume:      {len(cached)} cached, {len(to_run)} to run"
           + ("  (--force)" if args.force else ""))
    eta_min = len(to_run) * (3 if args.smoke else 60)
    print(f"  ETA:         ~{eta_min} min\n")

    if cached:
        print("Reusing cached cells:")
        for k, r in cached:
            print(f"  [{k + 1:2d}/{len(pairs)}] sev={r.severity:>3d}%  "
                   f"γ̇_th={r.gamma_th:.1e}  η_dep={r.eta_deposit:.3f}   (cache)")
        print()

    results: list[CellResult] = [r for _, r in cached]
    for (k, sev, th, cell_dir) in to_run:
        print(f"[{k + 1:2d}/{len(pairs)}] sev={sev:>3d}%  γ̇_th={th:.1e}  "
               f"→ {cell_dir.name}", flush=True)
        try:
            res = run_one_cell(severity=sev, gamma_th=th,
                                output_dir=str(cell_dir),
                                n_steps=n_steps)
            results.append(res)
            print(f"     η_dep={res.eta_deposit:5.3f}  "
                   f"η_off={res.eta_offtarget:5.3f}  "
                   f"rel_frac={res.release_fraction:5.3f}  "
                   f"wall {res.walltime_s:.0f}s")
        except KeyboardInterrupt:
            print(f"\n  ⚠ interrupted by user; {len(results)} cells "
                   f"complete. Re-launch to resume.")
            break
        except Exception as e:
            print(f"     FAILED: {e}")
            continue

    if not results:
        print("No successful cells. Bailing.")
        return 1

    field_names = list(CellResult.__dataclass_fields__.keys())
    cols = {name: np.asarray([getattr(r, name) for r in results])
             for name in field_names}
    np.savez(script_dir / "sweep_results.npz",
              **cols,
              severity_values=np.asarray(sev_vals),
              gamma_th_values=np.asarray(th_vals),
              n_steps=n_steps,
              n_cells=len(results))

    print(f"\nWrote {script_dir / 'sweep_results.npz'} ({len(results)} cells)")
    # Top-3 by deposition efficiency
    best = sorted([r for r in results if not np.isnan(r.eta_deposit)],
                   key=lambda r: -r.eta_deposit)[:3]
    print(f"\nTop-3 cells by η_deposit:")
    for r in best:
        print(f"  sev={r.severity:>3d}%  γ̇_th={r.gamma_th:.1e}  "
               f"η_deposit = {r.eta_deposit:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
