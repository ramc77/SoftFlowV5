# DLD deformability-threshold manuscript

Submission-ready draft for **Microfluidics and Nanofluidics**
(Springer).

## Files

| File | What it is |
|---|---|
| `main.tex`           | The article (article class, ~5000 words) |
| `references.bib`     | BibTeX (DLD + LBM/IBM + capsule mechanics) |
| `figures/*.png`      | Headline + supplementary figures from `sweep_analyse.py` |
| `Makefile`           | One-command build (`make`) |
| `README.md`          | This file |

## Build (local)

Prerequisites: a TeX distribution that ships `pdflatex` and `bibtex`
(TeX Live, MacTeX, or MiKTeX).

```bash
make            # full build (pdflatex → bibtex → pdflatex × 2)
make clean      # remove aux files
make distclean  # also remove main.pdf
```

Or manually:

```bash
pdflatex main.tex
bibtex   main
pdflatex main.tex
pdflatex main.tex
```

The output is `main.pdf`.

## Switching to the Springer Nature submission template

For the actual submission, replace the first non-comment line of
`main.tex`:

```latex
\documentclass[11pt,a4paper]{article}
```

with the Springer Nature template:

```latex
\documentclass[sn-mathphys-num]{sn-jnl}
```

after downloading `sn-jnl.cls` + `sn-mathphys.bst` from the
publisher's LaTeX support page:
https://www.springernature.com/gp/authors/campaigns/latex-author-support

Comment out the `\linenumbers` and `\onehalfspacing` lines (the
Springer template handles those for you).

## Submission checklist

- [ ] Title, author, affiliation, ORCID
- [ ] Abstract under 300 words
- [ ] 5-6 keywords
- [ ] Body length 5000–7000 words (currently ~5200)
- [ ] All figures referenced and present in `figures/`
- [ ] All BibTeX entries have DOIs
- [ ] Compliance / conflict-of-interest declaration filled
- [ ] Data availability statement (already includes the GitHub URL)
- [ ] Cover letter to editor (draft separately)
- [ ] Suggest 3–5 reviewers in the submission portal

## Figure inventory

| Fig | File | Purpose |
|---|---|---|
| 1 | `fig_sweep_dtheta_corrected.png`        | Δθ_corr heat-map (main result) |
| 2 | `fig_sweep_dtheta_corrected_vs_Ca.png`  | Ca-collapse + threshold (headline) |
| S1 | `fig_sweep_dtheta.png`                  | Raw Δθ heat-map (before correction) |
| S2 | `fig_sweep_dtheta_vs_Ca.png`            | Raw Δθ vs Ca |
| S3 | `fig_sweep_lane_order.png`              | Lane order (negative control) |

## Notes / TODOs before submission

1. **Replicate seeds for the headline row.** Currently a single seed
   per cell. Run 3 additional seeds for the `G_s_soft = 0.030` row
   (5 cells × 3 seeds ≈ 1.5 h compute) and add ± σ error bars to
   Fig. 2. This is the single biggest credibility improvement
   reviewers will ask for.

2. **Replace `henon2017deterministic` BibTeX entry.** The cite-key
   refers to the published 2017 *Biomicrofluidics* paper; the exact
   title/page may need verification at proof stage.

3. **Cover-letter draft.** Three sentences: novel finding (threshold),
   methodological contribution (matched-membrane control + diagonal
   subtraction), open-source release.

4. **Confirm the Springer Nature template version** matches the
   journal's current submission system before swapping.
