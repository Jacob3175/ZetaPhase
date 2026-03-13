# ZetaPhase

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18994912.svg)](https://doi.org/10.5281/zenodo.18994912)

**A logarithmic phase mapping of physical constants using the imaginary part of the first Riemann zeta zero.**

Author: Jacob Powell  
ORCID: https://orcid.org/0009-0003-5641-4413

---

## Core Result

- Cosmological hierarchy: ln(R / ℓ_P) ≈ 141.46
- First zeta zero scaling: 10 × γ₁ ≈ 141.34725 (difference < 0.08 %)
- 39 standard constants (PDG, NIST, CODATA) map to a unit circle in log-space.

Clear clustering appears near simple fractional phases:

Dimensionless couplings → 0  
Particle masses → 1/3  
Thermal / neutrino scales → 1/2  
Gravitational couplings → 2/3  

Monte Carlo null test (100,000 trials): p ≲ 10⁻⁵  
Optimal scaling confirmed at C ≈ 10.

---

## Mathematical Mapping

Physical constants are mapped to a phase coordinate using

φ(Q) = ln(Q / Q_P) / (10γ₁) mod 1

where

Q = physical constant  
Q_P = corresponding Planck-scale normalization  
γ₁ ≈ 14.134725 = first non-trivial zero of the Riemann zeta function

This converts quantities spanning many orders of magnitude into a circular coordinate on the interval [0,1).

---

## Files in this repository

data/constants.csv  
Exact dataset of the 39 physical constants used in the analysis.

analysis/zeta_phase_analysis.py  
Full reproduction script using mpmath, numpy, and matplotlib.

paper/zeta_phase_paper.tex  
LaTeX source of the research paper.

plots/clusters.png  
Circular phase mapping visualization.

monte_carlo_results.txt  
Monte Carlo statistical summary.

---

## How to reproduce the analysis

Install dependencies

pip install mpmath numpy matplotlib

Run the analysis

python analysis/zeta_phase_analysis.py

This regenerates the clustering figure and statistical results.

---

## Citation

Powell, Jacob (2026)

ZetaPhase: A Log-Phase Mapping of Physical Constants Using the First Riemann Zeta Zero

DOI  
https://doi.org/10.5281/zenodo.18994912

Repository  
https://github.com/Jacob3175

---

## Notes

This repository contains an exploratory numerical study connecting logarithmic physical scale hierarchies with the first non-trivial zero of the Riemann zeta function.

The results are presented as a phenomenological observation rather than a derived physical theory.

Further testing using larger datasets and blind prediction methods is encouraged.
