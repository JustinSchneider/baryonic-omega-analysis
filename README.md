# Baryonic Omega Analysis

Upgrading the empirical omega velocity correction model (Flynn & Cannaliato 2025) from point-mass approximations to full baryonic mass decomposition using the Corbelli method.

## Overview

This project fits the omega parameter to galaxy rotation curves using decomposed baryonic velocity components (gas, disk, bulge) from the SPARC database, rather than relying on Keplerian point-mass approximations.

**Fitting equation:**

```
V_model(R) = V_bary(R) + omega * R
```

where `V_bary = sqrt(V_gas^2 + Upsilon_disk * V_disk^2 + Upsilon_bulge * V_bulge^2)`

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Initialize database
python src/database.py --init

# Ingest a SPARC galaxy
python src/ingest.py --file data/raw/NGC5055_rotmod.dat

# Run omega fit
python src/fit.py --galaxy NGC5055 --plot
```

## Project Structure

- `src/` — Core Python modules (database, physics, ingestion, fitting)
- `tests/` — Pytest test suite
- `notebooks/` — Analysis and visualization notebooks
- `data/raw/` — Original SPARC data files
- `data/processed/` — SQLite database
- `results/` — Figures and summary tables

## References

- Corbelli & Salucci (2000). *The Extended Rotation Curve and the Dark Matter Halo of M33.*
- Lelli, McGaugh, & Schombert (2016). *SPARC: Mass Models for 175 Disk Galaxies.*
- Flynn & Cannaliato (2025). *A New Empirical Fit to Galaxy Rotation Curves.*
