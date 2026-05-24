# Stenosis manuscript — build instructions

Generic LaTeX format (article class) — submission-portable.

## Files

| File | What it is |
|---|---|
| `main.tex`           | Full manuscript (~6 000 words, 6 sections, 7 figures) |
| `references.bib`     | BibTeX (15 entries — all Crossref-verified) |
| `figures/*.png`      | All 7 figures from the analysis pipeline |
| `Makefile`           | One-command build (`make`) |
| `README.md`          | This file |

## Build (local)

Prerequisites: any TeX distribution with `pdflatex` and `bibtex`
(TeX Live, MacTeX, MiKTeX).

```bash
make            # full build (pdflatex → bibtex → pdflatex × 2)
make clean      # remove aux files
make distclean  # also remove main.pdf
```

Output: `main.pdf` — submit-ready PDF in generic article format.

## Switching to a specific journal template

The manuscript body is fully portable. Only the *header block*
(documentclass, title macros, author/affiliation block, bibliography
style) needs to change for each journal.

### Journal of Controlled Release (Elsevier, recommended)

Download `elsarticle.cls` from
<https://www.elsevier.com/researcher/author/policies-and-guidelines/latex-instructions>
and replace the first ~50 lines of `main.tex` with:

```latex
\documentclass[review,12pt]{elsarticle}
\usepackage{lineno}
\usepackage{amsmath,amssymb,graphicx,booktabs,hyperref}
\graphicspath{{figures/}}
\bibliographystyle{elsarticle-num}
\journal{Journal of Controlled Release}

\title{...}
\author[1]{Ram Chand\corref{cor1}}
\cortext[cor1]{Corresponding author: ram.chand@bnbwu.edu.pk}
\affiliation[1]{organization={Department of Natural Sciences,
   Begum Nusrat Bhutto Women University},
   city={Sukkur, Sindh}, country={Pakistan}}
```

### Lab on a Chip (RSC)

Download `rsc.cls`. Replace header with:

```latex
\documentclass[journal=loc]{rsc}
\graphicspath{{figures/}}
\bibliographystyle{rsc}
```

### Microfluidics and Nanofluidics (Springer)

Download `sn-jnl.cls` from
<https://www.springernature.com/gp/authors/campaigns/latex-author-support>:

```latex
\documentclass[sn-mathphys-num]{sn-jnl}
\graphicspath{{figures/}}
\bibliographystyle{spbasic}
```

## Submission checklist

- [ ] Title, author, affiliation, ORCID
- [ ] Abstract ≤ 300 words (currently ~280)
- [ ] 5–7 keywords (currently 7)
- [ ] Body length 5 000–7 000 words (currently ~6 200)
- [ ] All 7 figures present in `figures/` and referenced in body
- [ ] All 15 BibTeX entries have DOIs (verified via Crossref)
- [ ] Conflict-of-interest declaration filled
- [ ] Data-availability statement filled (references the project's
      `run_manifest.json` per-cell provenance pattern)
- [ ] Cover letter (draft separately — pitch in section above)
- [ ] Suggest 3–5 reviewers in submission portal

## Figure inventory

| Fig | File | Purpose |
|---|---|---|
| 1 | `fig_schematic.png`         | System schematic + zone annotations |
| 2 | `fig_snapshots.png`         | 3-panel spatial trajectory of sweet-spot cell |
| 3 | `fig_eta_deposit.png`       | **HEADLINE** — phase diagram of η_deposit |
| 4 | `fig_sweet_spot.png`        | η_deposit vs threshold, one curve per severity |
| 5 | `fig_trajectory.png`        | Single-cell deposition cascade time series |
| 6 (supp) | `fig_eta_offtarget.png`     | Off-target loss diagnostic |
| 7 (supp) | `fig_release_fraction.png`  | Release activity diagnostic |

## Notes / TODOs before submission

1. **Cover letter** (~150 words) — drafted in earlier conversation;
   pitch the Goldilocks finding as the headline novelty.

2. **Replicate seeds for sweet-spot cells** (optional, ~16 h compute) —
   would tighten error-bar reporting if a reviewer asks. The Korin
   2012 *Science* paper used single-seed simulations, so this is
   not blocking for a first submission.

3. **Optional ParaView snapshots with drug-field overlay** — current
   `fig_snapshots.png` shows carrier positions only. Adding a drug
   concentration heatmap on top would strengthen Figure 2 visually.
   ParaView instructions in the project README.

4. **3D extension** — flagged as future work in §5. Reviewers may
   request, but it is a substantial separate project. The 2D
   Goldilocks finding is publishable on its own.
