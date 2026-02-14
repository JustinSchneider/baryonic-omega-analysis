# Data Sources

## SPARC (Spitzer Photometry & Accurate Rotation Curves)

- **URL:** http://astroweb.cwru.edu/SPARC/
- **Reference:** Lelli, McGaugh, & Schombert (2016) - "SPARC: Mass Models for 175 Disk Galaxies"
- **Format:** Per-galaxy `_rotmod.dat` files with columns: Rad, Vobs, errV, Vgas, Vdisk, Vbul, SBdisk, SBbul
- **Units:** kpc (radius), km/s (velocities), L_sun/pc^2 (surface brightness)
- **Note:** V_disk and V_bul are provided at mass-to-light ratio Upsilon = 1 M_sun/L_sun

## M33 / NGC 598 (Calibration Target)

- **Primary:** SPARC catalog (if included as NGC0598)
- **Secondary:** Corbelli et al. (2014) via VizieR catalog J/MNRAS/442/2883
- **Reference:** Corbelli & Salucci (2000) - "The Extended Rotation Curve and the Dark Matter Halo of M33"
