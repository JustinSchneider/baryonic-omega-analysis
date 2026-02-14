# M33 Calibration Data — Status & Options

**Date:** 2026-02-14
**Status:** BLOCKED — No machine-readable M33 rotation curve data available

## The Problem

M33 (NGC 598) is our Phase I calibration target. We need to replicate Corbelli & Salucci (2000)'s mass model before trusting our pipeline on other galaxies. However:

1. **M33 is not in the SPARC catalog.** The 175-galaxy SPARC database (Lelli et al. 2016) does not include M33/NGC0598. Confirmed by searching the ingested data.

2. **Corbelli & Salucci (2000) did not publish machine-readable data.** The rotation curve (Figure 2b) and mass decomposition (Figure 6) exist only as plots in the paper. There is no supplementary data table, no CDS/VizieR archive entry.

3. **Corbelli et al. (2014) — "Dynamical signatures of a ΛCDM-halo and the distribution of the baryons in M33" (A&A 572, A23)** — is referenced in the Setup Guide as having tabulated data. However:
   - The VizieR catalog `J/A+A/572/A23` does not exist (returns "not found").
   - No machine-readable table was located via VizieR search, web search, or the A&A article page.

4. **Digitizing from figures is not reproducible.** Different team members would extract different values, and the uncertainty is unquantifiable.

## What We Know from Corbelli & Salucci (2000)

- **Distance:** D = 0.7 Mpc (1 arcmin = 0.2 kpc, or 5 arcmin = 1 kpc)
- **Disk scale length:** R_d = 1.2 kpc
- **Best-fit M_d/L_B:** 0.8 ± 0.2
- **V_gas:** Peaks at ~40 km/s near R ~ 8 kpc, declines to ~30 km/s at R ~ 16 kpc
- **V_lum (stars+gas):** Peaks at ~80 km/s near R ~ 3 kpc (Figure 5)
- **V_obs:** Rises to ~120 km/s at R ~ 3 kpc, continues rising to ~130 km/s at R ~ 16 kpc
- **Binned rotation curve:** 20 data points from ~0.5 to 16 kpc (Figure 2b), errors ~2-6%
- **Bulge:** Negligible — "the central bulge is very small and can be completely neglected" (Section 2a)
- **β_d(1):** 0.42 (fraction of disk velocity to total at R = R_opt)

## Options for the Team

### Option A: Contact the Authors Directly
Email Edvige Corbelli (Osservatorio di Arcetri, Florence) requesting the tabulated rotation curve and component velocities. This is standard practice in astronomy and would give us the definitive dataset.

**Pros:** Authoritative, reproducible, may include unpublished updates.
**Cons:** Response time uncertain.

### Option B: Use Corbelli (2003) / Corbelli (2014) Data
Corbelli published updated M33 rotation curves in later papers:
- Corbelli (2003) MNRAS 342, 199 — "ISM and star formation in M33"
- Corbelli et al. (2014) A&A 572, A23 — "Dynamical signatures of a ΛCDM-halo"

These may have tabulated data that we haven't been able to locate on VizieR. Check the actual journal PDF supplementary materials, or contact CDS directly.

### Option C: Use a SPARC Galaxy as Interim Calibration Target
Pick a well-studied SPARC galaxy that has been independently validated in the literature, and use it as our Phase I benchmark instead of (or alongside) M33. Good candidates:
- **NGC2403** — Late-type spiral similar to M33, extensively modeled in the literature
- **NGC3198** — Classic rotation curve benchmark, well-decomposed
- **DDO154** — Gas-dominated dwarf, simple decomposition

**Pros:** Data already ingested and ready. Immediate progress.
**Cons:** Not the specific Corbelli benchmark specified in the project plan.

### Option D: Use THINGS Survey Data for M33
The HI Nearby Galaxy Survey (THINGS) includes M33 and provides high-resolution HI data. The derived rotation curves from de Blok et al. (2008) are available as FITS tables. This would require additional parsing but gives us a peer-reviewed, machine-readable M33 dataset.

**Reference:** de Blok et al. (2008), AJ, 136, 2648 — "High-Resolution Rotation Curves and Galaxy Mass Models from THINGS"

**Update (2026-02-14):** Investigation revealed that M33 is likely **NOT** in the THINGS mass model sample. de Blok et al. (2008) published rotation curve decompositions for only 19 of 34 THINGS galaxies, and M33 is not among them — the paper actually cites Corbelli for M33 data. However, THINGS remains valuable for adding 19 additional calibration targets beyond SPARC.

## Recommendation

~~Pursue **Option A** (author contact) and **Option C** (interim SPARC calibration) in parallel.~~

**Updated approach (2026-02-14):** We have implemented:
1. **VizieR query pipeline** (`src/vizier_query.py`) that attempts to fetch Corbelli (2000) and Corbelli (2014) tabulated data from `J/MNRAS/311/441` and `J/A+A/572/A23`
2. **Manual fallback** (`create_m33_manual_data()` in `src/ingest.py`) with ~19 data points transcribed from Corbelli & Salucci (2000) paper text/figures — sufficient for pipeline validation
3. **THINGS ingestion** (`ingest_things_catalog()`) to add ~19 additional galaxies from de Blok et al. (2008) via VizieR (`J/AJ/136/2648`)
4. **Data provenance** — `data_source` column added to `galaxies` table to track origin (SPARC, Corbelli2000_VizieR, Corbelli2000_manual, THINGS)

Run `python src/ingest.py --m33` to ingest M33, or `python src/ingest.py --things` for the THINGS catalog.

## Action Items

- [ ] David/Jim: Reach out to Corbelli for tabulated M33 data (still valuable for authoritative values)
- [x] Justin: Implement VizieR query + manual fallback pipeline for M33
- [x] Justin: Implement THINGS survey ingestion pipeline
- [x] All: Review Option D — **finding: M33 not in THINGS mass models, but 19 other galaxies available**
- [ ] Justin: Run M33 calibration notebook end-to-end once VizieR queries are tested
- [ ] All: If Corbelli VizieR data is available, compare against manual transcription and update
