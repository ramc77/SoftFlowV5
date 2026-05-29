# Manuscript — clogging prediction with a graph neural network

Target journal: **Physica A: Statistical Mechanics and its Applications**
(Elsevier). The clogging/jamming transition + contact-network + machine-
learning framing fits Physica A's scope (statistical mechanics of complex,
disordered and athermal systems).

## Files

| File | What it is |
|---|---|
| `main.tex`        | Full manuscript (~10 pp, 4 figures, physics + GNN equations) |
| `references.bib`  | 15 references; journal DOIs verified via Crossref |
| `figures/`        | fig1 phase transition, fig2 GNN training+ROC, fig3 LORO, fig4 graphs |
| `Makefile`        | `make` builds the PDF |

## Build

```bash
make            # pdflatex -> bibtex -> pdflatex x2  -> main.pdf
make clean      # remove aux files
```
Builds with any TeX distribution (generic `article` class; falls back to
`plainnat` if `elsarticle-num.bst` is not installed).

## Switching to the Elsevier (Physica A) template for submission

Download `elsarticle.cls` and `elsarticle-num.bst` from
<https://www.elsevier.com/researcher/author/policies-and-guidelines/latex-instructions>
and replace the header block of `main.tex` (documentclass + title/author)
with:

```latex
\documentclass[preprint,12pt]{elsarticle}
\usepackage{amsmath,amssymb,graphicx,booktabs,hyperref}
\graphicspath{{figures/}}
\journal{Physica A}
\bibliographystyle{elsarticle-num}

\begin{frontmatter}
\title{Predicting clogging from the contact network: a graph neural network
       for deformable-particle suspensions at a microfluidic constriction}
\author{Ram Chand}
\affiliation{organization={Department of Natural Sciences,
   The Begum Nusrat Bhutto Women University}, city={Sukkur, Sindh},
   country={Pakistan}}
\begin{abstract} ... \end{abstract}
\begin{keyword} clogging \sep jamming \sep contact network \sep
   graph neural network \sep lattice Boltzmann \sep immersed boundary \end{keyword}
\end{frontmatter}
```
The manuscript body (Sections 1–5) is portable and needs no changes.

## Headline result

A 2-D LBM–IBM model of soft capsules at a constriction reproduces the
arching threshold near `D/d ≈ 3`; an edge-aware GNN classifies clog vs flow
from the instantaneous contact graph with ROC AUC `0.93` (held-out frames)
and a leave-one-run-out generalisation accuracy of `0.87` to unseen
apertures, failing only at the single aperture on the transition.
```
