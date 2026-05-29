#!/usr/bin/env python3
"""Demo: extract a graph dataset from a stenosis trajectory (Flavor-B data layer).

Builds one configuration graph per saved frame from a cell's
``particle_data.csv`` — nodes = carriers, edges = near-contacts — and reports
the node/edge feature schema and basic graph statistics. This is the
dependency-free data layer a graph neural network would consume; no PyTorch is
needed to produce or inspect it.

Usage
-----
    python research/04_shear_activated_stenosis/graph_extraction_demo.py
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

_REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "python"))

from pysoftflow.ml import build_graph_dataset, save_graph_dataset  # noqa: E402

# Geometry / particle parameters from research/04.../run.py
NX, NY = 400.0, 80.0
R_CARRIER = 2.0
CSV = (
    pathlib.Path(__file__).resolve().parent
    / "sweep_out/sev50_th2/vtk_stenosis/particle_data.csv"
)


def main() -> None:
    print(f"Trajectory: {CSV}")
    graphs = build_graph_dataset(
        str(CSV),
        radius=R_CARRIER,
        box=(NX, NY),
        periodic_x=True,
        contact_cutoff=1.0,
        frame_stride=50,   # subsample frames for the demo
    )
    print(f"Built {len(graphs)} graphs (one per sampled frame)\n")

    g0 = graphs[0]
    print("Node feature schema:", g0.node_feature_names)
    print("Edge feature schema:", g0.edge_feature_names)
    print(f"Per-graph nodes    : {g0.n_nodes}")

    degrees = np.array([g.mean_degree for g in graphs])
    n_edges = np.array([g.n_edges for g in graphs])
    contacts = n_edges // 2
    print("\nContact statistics across frames")
    print(f"  mean degree      : {degrees.mean():.3f}  (min {degrees.min():.2f}, "
          f"max {degrees.max():.2f})")
    print(f"  undirected edges : mean {contacts.mean():.1f}, max {int(contacts.max())}")
    busiest = int(np.argmax(n_edges))
    print(f"  busiest frame    : step {graphs[busiest].step} "
          f"with {int(contacts[busiest])} contacts")

    # Show a sample node-feature row and an edge-feature row.
    print("\nExample node feature (capsule 0, first frame):")
    for name, val in zip(g0.node_feature_names, g0.node_features[0]):
        print(f"    {name:>16s} = {val:+.4f}")
    gb = graphs[busiest]
    if gb.n_edges > 0:
        print(f"\nExample edge feature (busiest frame, edge 0):")
        for name, val in zip(gb.edge_feature_names, gb.edge_features[0]):
            print(f"    {name:>16s} = {val:+.4f}")

    out = pathlib.Path(__file__).resolve().parent / "graph_dataset_demo.npz"
    save_graph_dataset(graphs, str(out))
    size_kb = out.stat().st_size / 1024
    print(f"\nSaved graph dataset -> {out.name}  ({size_kb:.1f} KB, {len(graphs)} graphs)")
    print("PyG bridge: call ParticleGraph.to_pyg() once torch_geometric is installed.")


if __name__ == "__main__":
    main()
