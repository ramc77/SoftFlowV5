"""Project 4 — status reporter (read-only, instant).

    python status.py            # production grid
    python status.py --smoke    # smoke grid
"""
from __future__ import annotations
import argparse, json, pathlib, sys
script_dir = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(script_dir))
from sweep import (SEVERITIES_PROD, SEVERITIES_SMOKE,                  # noqa: E402
                    GAMMA_TH_PROD, GAMMA_TH_SMOKE, _cell_dirname)
SWEEP_ROOT = script_dir / "sweep_out"
def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--smoke", action="store_true")
    return p.parse_args(argv)
def main(argv=None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    sev_vals = SEVERITIES_SMOKE if args.smoke else SEVERITIES_PROD
    th_vals  = GAMMA_TH_SMOKE if args.smoke else GAMMA_TH_PROD
    label    = "SMOKE" if args.smoke else "PRODUCTION"
    pairs = [(s, t) for s in sev_vals for t in th_vals]
    if not SWEEP_ROOT.exists():
        print(f"sweep_out/ not found. {len(pairs)} cells pending.")
        return 0
    n_done = n_pending = n_partial = 0
    print(f"\nProject 4 stenosis — status ({label}, {len(pairs)} cells)")
    print("=" * 96)
    print(f"  {'cell':<14s} {'status':<12s} {'η_dep':>8s} {'η_off':>8s} "
           f"{'rel_frac':>9s} {'walltime':>9s}")
    print("-" * 96)
    for s, th in pairs:
        cd = SWEEP_ROOT / _cell_dirname(s, th, th_vals)
        f = cd / "cell_result.json"
        if f.exists():
            try:
                d = json.loads(f.read_text()); n_done += 1
                ed = d.get("eta_deposit"); eo = d.get("eta_offtarget")
                rf = d.get("release_fraction"); wt = d.get("walltime_s")
                ns = d.get("_n_steps", "?")
                print(f"  {cd.name:<14s} {'✓ done':<12s} "
                       f"{('—' if ed is None else f'{ed:.3f}'):>8s} "
                       f"{('—' if eo is None else f'{eo:.3f}'):>8s} "
                       f"{('—' if rf is None else f'{rf:.3f}'):>9s} "
                       f"{('—' if wt is None else f'{wt:.0f}s'):>9s}"
                       f"   (n_steps={ns})")
            except (json.JSONDecodeError, OSError):
                n_partial += 1; print(f"  {cd.name:<14s} ⚠ corrupt")
        elif (cd / "vtk_stenosis").exists():
            n_partial += 1; print(f"  {cd.name:<14s} ⚠ partial")
        else:
            n_pending += 1; print(f"  {cd.name:<14s} · pending")
    print("=" * 96)
    print(f"  Total: {len(pairs)}   Done: {n_done}   "
           f"Pending: {n_pending}   Partial: {n_partial}")
    if n_pending + n_partial > 0:
        eta = (n_pending + n_partial) * (3 if args.smoke else 6)
        print(f"  → Resume:  python sweep.py"
               f"{' --smoke' if args.smoke else ''}   (~{eta} min)")
    else:
        print("  → All done. Run sweep_analyse.py for figures.")
    return 0
if __name__ == "__main__":
    sys.exit(main())
