# Baryonic Omega Analysis

Upgrading the empirical omega velocity correction model (Flynn & Cannaliato 2025) from point-mass approximations to full baryonic mass decomposition using the Corbelli method.

## Overview

This project fits the omega parameter to galaxy rotation curves using decomposed baryonic velocity components (gas, disk, bulge) from the SPARC database, rather than relying on Keplerian point-mass approximations.

**Core model (Linear):**

$$V_{model}(R) = V_{bary}(R) + \omega \cdot R$$

where $V_{bary} = \sqrt{|V_{gas}| \cdot V_{gas} + \Upsilon_d \cdot |V_{disk}| \cdot V_{disk} + \Upsilon_b \cdot |V_{bulge}| \cdot V_{bulge}}$

**Extended model (Rational Taper):**

$$V_{model}(R) = V_{bary}(R) + \frac{\omega \cdot R}{1 + R / R_t}$$

The taper introduces a transition radius $R_t$ where the linear correction saturates, producing flat rotation at large radii.

## Phase I Results

Phase I analysis (118 SPARC galaxies) established:

- **Linear > Quadrature:** BIC strongly favors a kinematic (linear velocity addition) over a dynamic (quadrature force addition) model, with $\Delta\text{BIC} = 4559$ on the M33 calibration target.
- **Tapered model:** A rational taper reduces M33 residual RMSE by 69% (31.1 km/s to 9.5 km/s) by saturating the linear correction at a characteristic radius.
- **Candidate scaling law:** The transition radius $R_t \approx k \cdot R_d$ with median $k \approx 2.5$ among well-fit galaxies, suggesting the correction scale is tied to the baryonic disk.

See [docs/PHASE_I_RESULTS.md](docs/PHASE_I_RESULTS.md) for full analysis and [docs/METHODOLOGY.md](docs/METHODOLOGY.md) for detailed methods.

### Key Figures

![M33 Tapered Linear Models](results/figures/M33_tapered_models.png)

*Comparison of pure linear, rational taper, and tanh taper fits to M33 (Corbelli 2014 data). The rational taper (green) best reproduces the flat outer rotation curve (RMSE = 9.5 km/s), while the pure linear model (red) diverges. Lower panel shows residuals.*

![Coupling Constant Universality](results/figures/k_universality.png)

*The coupling constant $k = R_t / R_d$ across 118 SPARC galaxies, plotted against flat velocity, Hubble type, and disk scale length. The lack of strong systematic trends is consistent with a universal value ($k \approx 2.5\text{–}2.8$), though 19% of fits hit parameter bounds.*

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

- Casertano, S. (1983). MNRAS, 203, 735. *Rotation curve of the edge-on spiral galaxy NGC 5907.*
- Corbelli, E. & Salucci, P. (2000). MNRAS, 311, 441. *The Extended Rotation Curve and the Dark Matter Halo of M33.*
- Corbelli, E., et al. (2014). A&A, 572, A23. *Dynamical signatures of a LCDM-halo and the distribution of the baryons in M33.*
- Flynn, D. C. & Cannaliato, J. (2025). *A New Empirical Fit to Galaxy Rotation Curves.*
- Lelli, F., McGaugh, S. S., & Schombert, J. M. (2016). AJ, 152, 157. *SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves.*
