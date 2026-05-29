#!/usr/bin/env python3
"""Composite the four ParaView snapshots into one 2x2 publication figure.

Inputs (in paper/figures/, supplied from ParaView):
  fig5_ss1.png         clog,   initial      fig5_ss1-final.png  clog,   final
  fig6_ss_initial.png  no-clog, initial     fig6_ss_final.png   no-clog, final

Each frame shows the fluid speed field (background) with the capsules coloured
by their speed (blue slow -> red fast), i.e. fluid + particle tracking. The
script auto-crops the uniform border of each frame and arranges them as
rows = {clog, no-clog}, columns = {initial, final}, with panel labels.

Output: paper/figures/fig_snapshots.png
"""

from __future__ import annotations

import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

FIGDIR = pathlib.Path(__file__).resolve().parent / "paper" / "figures"

# Common panel size (pixels) so every subplot is identical.
TARGET_H, TARGET_W = 210, 940


def crop_to_channel(img: np.ndarray, sat_tol: float = 0.18, pad: int = 3) -> np.ndarray:
    """Crop tightly to the coloured fluid channel.

    The channel (blue->red speed field with capsules) is highly saturated,
    while both the white and the grey ParaView backgrounds are unsaturated.
    Masking on colour saturation therefore removes the background uniformly
    regardless of its shade. The largest connected saturated component is the
    channel; smaller ones (e.g. the orientation-axis widget) are discarded.
    """
    rgb = img[..., :3]
    sat = mcolors.rgb_to_hsv(rgb)[..., 1]
    mask = sat > sat_tol
    if not mask.any():
        return img
    labels, n = ndimage.label(mask)
    if n > 1:
        sizes = np.bincount(labels.ravel())
        sizes[0] = 0
        mask = labels == int(sizes.argmax())
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, img.shape[0])
    c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, img.shape[1])
    return rgb[r0:r1, c0:c1]


def resize(img: np.ndarray, h: int = TARGET_H, w: int = TARGET_W) -> np.ndarray:
    """Resize an RGB image to a fixed (h, w) so all panels are equal in size."""
    zh, zw = h / img.shape[0], w / img.shape[1]
    return np.clip(ndimage.zoom(img, (zh, zw, 1), order=1), 0.0, 1.0)


def main() -> None:
    panels = [
        ("fig5_ss1.png", "(a) clog — initial"),
        ("fig5_ss1-final.png", "(b) clog — final"),
        ("fig6_ss_initial.png", "(c) no-clog — initial"),
        ("fig6_ss_final.png", "(d) no-clog — final"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 4.2))
    for ax, (fname, label) in zip(axes.ravel(), panels):
        img = resize(crop_to_channel(plt.imread(str(FIGDIR / fname))))
        ax.imshow(img)
        ax.set_title(label, fontsize=11, loc="left")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    fig.suptitle("Fluid speed field with capsules coloured by speed "
                 "(blue = slow, red = fast); flow is left to right",
                 fontsize=10, y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = FIGDIR / "fig_snapshots.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
