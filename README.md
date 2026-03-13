# ZetaPhase

**A logarithmic phase mapping of physical constants using the imaginary part of the first Riemann zeta zero.**

### Core Result
- Cosmological hierarchy: `ln(R / ℓ_P) ≈ 141.46`  
- First zeta zero scaling: `10 × γ₁ ≈ 141.34725` (difference < 0.08 %)  
- 39 standard constants (PDG, NIST, CODATA) map to a unit circle in log-space  
- Clear clustering:  
  - Dimensionless couplings → **0**  
  - Particle masses → **1/3**  
  - Thermal/neutrino scales → **1/2**  
  - Gravitational couplings → **2/3**  

Monte Carlo null test (100,000 trials): **p ≲ 10^{-5}**.  
Optimal scaling confirmed at C ≈ 10.

### Files in this repo
- `data/constants.csv` — your exact 39-constant dataset  
- `analysis/zeta_phase_analysis.py` — full reproduction script (mpmath + numpy + matplotlib)  
- `paper/zeta_phase_paper.tex` — the final LaTeX (compiles to the PDF)  
- `plots/clusters.png` — high-resolution circular phase plot  
- `monte_carlo_results.txt` — statistical summary

### How to reproduce
```bash
pip install mpmath numpy matplotlib
python analysis/zeta_phase_analysis.py
