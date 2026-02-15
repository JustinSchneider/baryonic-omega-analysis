# Baryonic Omega Analysis

Upgrading the empirical omega velocity correction model (Flynn & Cannaliato 2025) from point-mass approximations to full baryonic mass decomposition using the Corbelli method.

## Overview

This project fits the omega parameter to galaxy rotation curves using decomposed baryonic velocity components (gas, disk, bulge) from the SPARC database, rather than relying on Keplerian point-mass approximations.

**Core model (Flynn & Cannaliato 2025 — Linear):**

$$V_{model}(R) = V_{bary}(R) + \omega \cdot R$$

**Extended model (Schneider 2026 — Rational Taper):**

$$V_{model}(R) = V_{bary}(R) + \frac{\omega \cdot R}{1 + R / R_t}$$

The taper introduces a transition radius $R_t$ where the linear correction saturates, producing flat rotation at large radii. We find that $R_t \approx k \cdot R_d$, where $R_d$ is the disk scale length and $k \approx 2.8$ is a candidate universal coupling constant.

## Key Results

### 1. Mechanism: Kinematic, Not Dynamic

BIC model comparison on the M33 calibration target strongly favors **linear velocity addition** over **quadrature force addition** ($\Delta\text{BIC} = -4559$). The omega correction acts as a velocity boost — not a dark matter-like potential. This rejects the standard DM halo addition mechanism.

### 2. The Tapered Model

The pure linear model diverges at large radii. The rational taper saturates the correction at a characteristic radius $R_t$, reducing M33 RMSE by **69%** (31.1 → 9.5 km/s) and yielding $\chi^2_\nu = 4.6$ vs. 72.9 for the linear form.

![M33 Linear vs Tapered](results/figures/M33_linear_reanalysis.png)

_M33 — the calibration target. The Flynn-Cannaliato (2025) linear model (red dashed, $\omega = 6.97$) diverges beyond ~8 kpc while the Schneider (2026) tapered model (orange, $\omega = 42.97$, $R_t = 2.0$ kpc) tracks the flat observed rotation curve. RMSE drops from 31.1 to 9.5 km/s ($\Delta$BIC = 3896)._

![NGC 3198 Gallery](results/figures/gallery_NGC3198.png)

_NGC 3198 — an intermediate spiral showing the same pattern across galaxy types. The linear model overshoots at large radii while the tapered model tracks the flat rotation curve with RMSE = 6.4 km/s._

### 3. Universal Coupling Constant

Batch analysis of 118 SPARC galaxies reveals a candidate scaling law: the transition radius $R_t$ is proportional to the disk scale length $R_d$ with a median coupling factor **$k = 2.81$** (full sample) or **$k = 2.55$** among the 89 galaxies with $\chi^2_\nu < 5$.

![Coupling Constant Distribution](results/figures/k_distribution.png)

_Distribution of the coupling constant $k = R_t / R_d$ across 118 SPARC galaxies. The primary peak near $k \approx 2$ and the secondary pile-up at $k = 20$ (the parameter bound) reveal two distinct populations: galaxies where the taper is well-constrained (84%) and those that exhibit extended linear-like rise before tapering (16%)._

### 4. Two Populations: Interior vs. Boundary Solutions

The $k$ distribution is bimodal. Galaxies splitting into two populations:

- **Interior solutions** ($k < 20$, N=99, 84%): The taper is well-constrained. These are predominantly lower-luminosity, lower surface brightness systems.
- **Boundary solutions** ($k = 20$, N=19, 16%): The optimizer hits the parameter bound — these galaxies prefer the pure linear model. They are systematically brighter, more massive, and have higher surface brightness.

![Population Split](results/figures/split_populations.png)

_Boxplots comparing luminosity, central surface brightness, and flat velocity between the two populations. Boundary-solution galaxies (linear-preferred) are systematically more luminous and denser — a physically meaningful distinction, not random fitting noise._

### 5. Robustness

Mass-to-light ratio sensitivity testing ($\Upsilon_d \in \{0.3, 0.5, 0.8\}$) shows that $\omega$ is most stable for gas-dominated (LSB) systems (28% variation for DDO 161) and most sensitive for disk-dominated (HSB) systems (134% for NGC 2841). The model's strength lies in the LSB regime where baryonic uncertainties are smallest.

### 6. Model Gallery

Head-to-head comparison of the Flynn & Cannaliato (2025) linear model vs. the Schneider (2026) tapered model across diverse galaxy types:

| Galaxy   | $\Sigma_0$ ($L_\odot$/pc$^2$) | Prediction    | Preferred (BIC) | Verdict   |
| -------- | ----------------------------: | ------------- | --------------- | --------- |
| DDO 154  |                            62 | Tapered (LSB) | Tapered         | Correct   |
| DDO 161  |                            59 | Tapered (LSB) | Tapered         | Correct   |
| NGC 0300 |                           152 | Tapered (LSB) | Tapered         | Correct   |
| NGC 3198 |                           618 | Transition    | Tapered         | —         |
| NGC 2841 |                          2260 | Linear (HSB)  | Tapered         | Incorrect |
| NGC 7331 |                          1583 | Linear (HSB)  | Linear          | Correct   |
| M33      |                             — | —             | Tapered         | —         |

The surface brightness predictor achieves **80% accuracy** (4/5 testable cases). The NGC 2841 failure suggests the HSB threshold needs refinement or that the tapered model is more broadly applicable than the population split implies.

## Documentation

- [**RESULTS.md**](docs/RESULTS.md) — Full Phase I & II results with figures and analysis
- [**METHODOLOGY.md**](docs/METHODOLOGY.md) — Detailed methods, equations, and fitting procedures
- [**CLAUDE.md**](CLAUDE.md) — Project plan, phases, and developer guidelines

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Initialize database
python src/database.py --init

# Ingest M33 from Corbelli 2014 Table 1
python src/ingest.py --m33

# Ingest all SPARC galaxies from MRT file
python src/ingest.py --mrt data/raw/MassModels_Lelli2016c.mrt --metadata data/raw/SPARC_Lelli2016c.mrt

# Run omega fit on a single galaxy
python src/fit.py --galaxy M33 --plot
```

## Notebooks

| #   | Notebook                                                                   | Description                                                |
| --- | -------------------------------------------------------------------------- | ---------------------------------------------------------- |
| 01  | [M33 Calibration](notebooks/01_m33_calibration.ipynb)                      | Pipeline validation against Corbelli 2014 data             |
| 02  | [Linear vs Quadrature](notebooks/02_linear_vs_quadrature_comparison.ipynb) | Mechanism test: kinematic boost vs. force addition         |
| 03  | [Tapered Models](notebooks/03_tapered_linear_model.ipynb)                  | Rational taper and tanh taper on M33                       |
| 04  | [SPARC Batch](notebooks/04_sparc_batch_analysis.ipynb)                     | 118-galaxy batch fit with $k \cdot R_d$ parameterization   |
| 05  | [Phase II Analysis](notebooks/05_phase2_density_coupling.ipynb)            | Population split, density coupling, $\Upsilon$ sensitivity |
| 06  | [Model Gallery](notebooks/06_model_gallery.ipynb)                          | Head-to-head Linear vs. Tapered across galaxy types        |

## Project Structure

- `src/` — Core Python modules (database, physics, ingestion, fitting)
- `tests/` — Pytest test suite
- `notebooks/` — Analysis and visualization notebooks
- `data/raw/` — Original SPARC data files
- `data/extracted/` — Data extracted from published papers (e.g., Corbelli 2014 Table 1)
- `data/processed/` — SQLite database
- `results/figures/` — Publication-quality plots
- `results/tables/` — Summary statistics and fit results (CSV)
- `docs/` — Methodology, results, and internal documentation

## Data Sources

- **SPARC Database** (Lelli, McGaugh, & Schombert 2016): 175 disk galaxies with Spitzer photometry and rotation curves. http://astroweb.cwru.edu/SPARC/
- **M33 Calibration** (Corbelli et al. 2014): Surface density profiles from Table 1, converted to velocity components via Casertano (1983) thin-disk method.

## References

1. Flynn, D. C. & Cannaliato, J. (2025). "A New Empirical Fit to Galaxy Rotation Curves."
2. Corbelli, E. & Salucci, P. (2000). MNRAS, 311, 441. "The Extended Rotation Curve and the Dark Matter Halo of M33."
3. Corbelli, E., et al. (2014). A&A, 572, A23. "Dynamical signatures of a $\Lambda$CDM-halo and the distribution of the baryons in M33."
4. Casertano, S. (1983). MNRAS, 203, 735. "Rotation curve of the edge-on spiral galaxy NGC 5907."
5. Kass, R. E. & Raftery, A. E. (1995). JASA, 90, 773. "Bayes Factors."
6. Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). AJ, 152, 157. "SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves."
