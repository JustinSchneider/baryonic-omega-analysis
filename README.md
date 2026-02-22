# Baryonic Coupling in Galaxy Rotation Curves

[![Status: Under Review](https://img.shields.io/badge/Status-Under_Review-yellow.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This repository contains the data, fitting pipeline, and full results table for the manuscript:
**"A Baryonically-Coupled Rational Taper Model for Galaxy Rotation Curves: Evidence from the Full SPARC Catalog"** (Schneider, 2026; _Awaiting submission_).

## Key Finding

![Dynamic Coupling: Asymptotic Velocity vs. Baryonic Mass](results/figures/pub_fig2_dynamic_coupling.png)

The Rational Taper model ($V_{model} = V_{bary} + \omega R / (1 + R/R_t)$) saturates at an asymptotic velocity $V_{sat} = \omega \cdot R_t$. Across 127 well-fit SPARC galaxies, $V_{sat}$ scales with total baryonic mass as a power law with slope $0.221 \pm 0.026$ — consistent with the expected $V \propto M^{0.25}$ Baryonic Tully-Fisher Relation ($R^2 = 0.364$, $p = 6.2 \times 10^{-14}$). This demonstrates that the model's parameters are not arbitrary curve-fitting artifacts but are physically coupled to global galactic mass distributions.

## Overview

This project explores a phenomenological extension to the empirical velocity correction model proposed by Flynn & Cannaliato (2025). By upgrading from point-mass approximations to full baryonic mass decompositions using the SPARC database, we test two kinematic models:

**1. Linear Model (Flynn & Cannaliato 2025):**
$$V_{model}(R) = V_{bary}(R) + \omega R$$

**2. Rational Taper Model (Schneider 2026):**
$$V_{model}(R) = V_{bary}(R) + \frac{\omega R}{1 + R / R_t}$$

The Tapered model introduces a transition radius $R_t$ where the linear correction saturates to $V_{sat} = \omega \cdot R_t$. Across 171 quality-controlled SPARC galaxies, BIC selects the Tapered model in 74.3% of cases, and the saturation scale couples to the baryonic disk: $R_t \approx 2.4 R_d$.

## Gemini Gem

An LLM has been configured to [engage with this research on Google Gemini](https://gemini.google.com/gem/1XhA_9T5KTLUSwiFIHFfAHLfPLNAtQLt5). This chatbot is capable of answering many questions about the analysis performed in this work, though as is always the case, responses should be verified via the data directly, available in this repository.

## Reproducing the Manuscript Figures

Reviewers and readers can reproduce the exact figures and statistical analyses found in the manuscript using the provided Jupyter Notebooks:

## Notebooks

| #   | Notebook                                                                   | Description                                                                                                                          |
| --- | -------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| 01  | [M33 Calibration](notebooks/01_m33_calibration.ipynb)                      | Pipeline validation against Corbelli 2014 data                                                                                       |
| 02  | [Linear vs Quadrature](notebooks/02_linear_vs_quadrature_comparison.ipynb) | Mechanism test: kinematic boost vs. force addition                                                                                   |
| 03  | [Tapered Models](notebooks/03_tapered_linear_model.ipynb)                  | Rational taper and tanh taper on M33                                                                                                 |
| 04  | [SPARC Batch](notebooks/04_sparc_batch_analysis.ipynb)                     | 118-galaxy batch fit with $k \cdot R_d$ parameterization                                                                             |
| 05  | [Phase II Analysis](notebooks/05_phase2_density_coupling.ipynb)            | Population split, density coupling, $\Upsilon$ sensitivity                                                                           |
| 06  | [Model Gallery](notebooks/06_model_gallery.ipynb)                          | Head-to-head Linear vs. Tapered across galaxy types                                                                                  |
| 07  | [Full Catalog Analysis](notebooks/07_full_catalog_analysis.ipynb)          | Phase III: BIC selection, $R_t$–$R_d$ scaling, $\Sigma_0$ regime test on all 175 SPARC galaxies                                      |
| 08  | [Full Gallery](notebooks/08_full_gallery.ipynb)                            | Rotation-curve gallery for all 171 quality-controlled galaxies (29 pages, sorted by $\Delta$BIC)                                     |
| 09  | [Publication Figures](notebooks/09_publication_figures.ipynb)              | Generates the four manuscript-ready figures: M33 calibration, BTFR recovery, data truncation artifact, and $\Delta$ BIC distribution |
| 10  | [Surface Brightness vs. $\omega$](notebooks/10_sigma0_vs_omega_correlation.ipynb) | Pearson/Spearman correlation of $\omega$ against central surface brightness $\Sigma_0$; establishes that model preference is $\Sigma_0$-independent while coupling amplitude is weakly modulated by baryonic surface density |

## Quick Start

The pipeline is written in Python and uses SQLite for data management.

```bash
# Clone the repository
git clone https://github.com/JustinSchneider/baryonic-omega-analysis.git
cd baryonic-omega-analysis

# Install dependencies
pip install -r requirements.txt

# Initialize database and ingest SPARC/M33 data
python src/database.py --init
python src/ingest.py --m33
python src/ingest.py --mrt data/raw/MassModels_Lelli2016c.mrt --metadata data/raw/SPARC_Lelli2016c.mrt

# Run the comparative fit on a single galaxy (e.g., NGC 3198)
python src/fit.py --galaxy NGC3198 --plot
```

## Repository Structure

- `src/` — Core Python modules (database, physics, ingestion, fitting pipeline)
- `notebooks/` — Analysis and visualization notebooks mapped to manuscript figures
- `data/` — Raw SPARC data (LMS16) and the compiled SQLite database
- `results/figures/` — Publication-quality plots generated by the notebooks
- `results/tables/` — Full fit parameters and summary statistics (CSV)
- `docs/` — Internal documentation and mathematical methodology

## Data Sources

- **SPARC Database** (Lelli, McGaugh, & Schombert 2016): 175 disk galaxies with Spitzer photometry and rotation curves. http://astroweb.cwru.edu/SPARC/
- **M33 Calibration** (Corbelli et al. 2014): Surface density profiles from Table 1, converted to velocity components via Casertano (1983) thin-disk method.

## Citation

(Citation details will be updated upon publication. For preprint inquiries, please reference the GitHub URL).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
