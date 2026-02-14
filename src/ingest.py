"""Parsers and ingestion pipeline for galaxy rotation curve data.

Handles:
  - Individual SPARC _rotmod.dat files (one galaxy per file)
  - Combined MRT (Machine-Readable Table) files from Lelli et al. (2016)
    containing all 175 galaxies in a single fixed-width file.
  - Corbelli & Salucci (2000) / Corbelli et al. (2014) M33 data from VizieR
  - THINGS survey rotation curves from de Blok et al. (2008) via VizieR
"""

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from src.database import get_session, init_db, insert_galaxy, insert_radial_profiles
from src.physics import HELIUM_FACTOR, circular_velocity_thin_disk, compute_v_bary
from src.utils import get_project_root, setup_logger

logger = setup_logger(__name__)

SPARC_COLUMNS = ["Rad", "Vobs", "errV", "Vgas", "Vdisk", "Vbul", "SBdisk", "SBbul"]


def parse_sparc_rotmod(filepath: str | Path) -> pd.DataFrame:
    """Parse a single SPARC _rotmod.dat file.

    Skips comment lines starting with '!'.

    Args:
        filepath: Path to the .dat file.

    Returns:
        DataFrame with columns: Rad, Vobs, errV, Vgas, Vdisk, Vbul, SBdisk, SBbul.

    Raises:
        FileNotFoundError: If filepath does not exist.
        ValueError: If file has unexpected number of columns.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"SPARC file not found: {filepath}")

    df = pd.read_csv(
        filepath,
        sep=r"\s+",
        comment="!",
        names=SPARC_COLUMNS,
        header=None,
    )

    if len(df.columns) != 8:
        raise ValueError(
            f"Expected 8 columns, got {len(df.columns)} in {filepath.name}"
        )

    # Ensure all columns are float
    for col in SPARC_COLUMNS:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    logger.info("Parsed %s: %d data points", filepath.name, len(df))
    return df


def extract_galaxy_name_from_filename(filepath: str | Path) -> str:
    """Extract galaxy name from SPARC filename convention.

    Examples:
        'NGC5055_rotmod.dat' -> 'NGC5055'
        'UGC02885_rotmod.dat' -> 'UGC02885'
    """
    name = Path(filepath).stem  # e.g., 'NGC5055_rotmod'
    # Remove the '_rotmod' suffix if present
    if "_rotmod" in name:
        name = name.split("_rotmod")[0]
    return name


def ingest_sparc_file(
    filepath: str | Path,
    galaxy_name: Optional[str] = None,
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    db_path: Optional[str] = None,
) -> str:
    """Full pipeline: parse SPARC file -> compute V_bary -> insert into DB.

    Args:
        filepath: Path to SPARC _rotmod.dat file.
        galaxy_name: Override for galaxy identifier. If None, extracted from filename.
        upsilon_disk: Mass-to-light ratio for disk. Default 0.5.
        upsilon_bulge: Mass-to-light ratio for bulge. Default 0.7.
        db_path: Optional database path override.

    Returns:
        galaxy_id that was inserted.
    """
    filepath = Path(filepath)

    # Parse the raw data
    raw_df = parse_sparc_rotmod(filepath)

    # Determine galaxy name
    galaxy_id = galaxy_name or extract_galaxy_name_from_filename(filepath)

    # Compute baryonic velocity
    v_bary = compute_v_bary(
        raw_df["Vgas"].values,
        raw_df["Vdisk"].values,
        raw_df["Vbul"].values,
        upsilon_disk=upsilon_disk,
        upsilon_bulge=upsilon_bulge,
    )

    # Build the profiles DataFrame with DB column names
    profiles_df = pd.DataFrame(
        {
            "radius_kpc": raw_df["Rad"],
            "v_obs": raw_df["Vobs"],
            "v_err": raw_df["errV"],
            "v_gas": raw_df["Vgas"],
            "v_disk": raw_df["Vdisk"],
            "v_bulge": raw_df["Vbul"],
            "v_baryon_total": v_bary,
        }
    )

    # Initialize DB and insert
    engine = init_db(db_path)
    session = get_session(engine)
    try:
        insert_galaxy(session, galaxy_id)
        n_rows = insert_radial_profiles(session, galaxy_id, profiles_df)
        logger.info(
            "Ingested %s: %d profiles (Upsilon_d=%.2f, Upsilon_b=%.2f)",
            galaxy_id,
            n_rows,
            upsilon_disk,
            upsilon_bulge,
        )
    finally:
        session.close()

    return galaxy_id


def ingest_sparc_directory(
    dirpath: str | Path,
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    db_path: Optional[str] = None,
) -> list[str]:
    """Ingest all _rotmod.dat files in a directory.

    Args:
        dirpath: Directory containing SPARC .dat files.
        upsilon_disk: Mass-to-light ratio for disk.
        upsilon_bulge: Mass-to-light ratio for bulge.
        db_path: Optional database path override.

    Returns:
        List of galaxy_ids that were ingested.
    """
    dirpath = Path(dirpath)
    dat_files = sorted(dirpath.glob("*_rotmod.dat"))

    if not dat_files:
        logger.warning("No _rotmod.dat files found in %s", dirpath)
        return []

    logger.info("Found %d SPARC files in %s", len(dat_files), dirpath)

    ingested = []
    for f in dat_files:
        try:
            gid = ingest_sparc_file(
                f,
                upsilon_disk=upsilon_disk,
                upsilon_bulge=upsilon_bulge,
                db_path=db_path,
            )
            ingested.append(gid)
        except Exception as e:
            logger.error("Failed to ingest %s: %s", f.name, e)

    logger.info("Successfully ingested %d / %d galaxies", len(ingested), len(dat_files))
    return ingested


def parse_massmodels_mrt(filepath: str | Path) -> dict[str, pd.DataFrame]:
    """Parse the combined MassModels MRT file (Lelli et al. 2016, Table 2).

    This file contains rotation curve data for all 175 SPARC galaxies in a
    single fixed-width file. Format per the MRT header:
        Col  1-11: Galaxy ID (A11)
        Col 13-18: Distance (F6.2, Mpc)
        Col 20-25: Radius (F6.2, kpc)
        Col 27-32: Vobs (F6.2, km/s)
        Col 34-38: errV (F5.2, km/s)
        Col 40-45: Vgas (F6.2, km/s)
        Col 47-52: Vdisk (F6.2, km/s)
        Col 54-59: Vbul (F6.2, km/s)
        Col 61-67: SBdisk (F7.2, solLum/pc2)
        Col 69-76: SBbul (F8.2, solLum/pc2)

    Args:
        filepath: Path to MassModels_Lelli2016c.mrt file.

    Returns:
        Dict mapping galaxy_id -> DataFrame with standard SPARC columns.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"MRT file not found: {filepath}")

    # Read all lines, skip the header (lines before the data start)
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Find where data starts (after the last '---' separator line)
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("---"):
            data_start = i + 1

    # Parse fixed-width data
    records = []
    for line in lines[data_start:]:
        if len(line.strip()) == 0:
            continue
        try:
            galaxy_id = line[0:11].strip()
            # distance = float(line[12:18])  # Not needed for profiles
            rad = float(line[19:25])
            vobs = float(line[26:32])
            errv = float(line[33:38])
            vgas = float(line[39:45])
            vdisk = float(line[46:52])
            vbul = float(line[53:59])
            sbdisk = float(line[60:67])
            sbbul = float(line[68:76])
            records.append({
                "galaxy_id": galaxy_id,
                "Rad": rad,
                "Vobs": vobs,
                "errV": errv,
                "Vgas": vgas,
                "Vdisk": vdisk,
                "Vbul": vbul,
                "SBdisk": sbdisk,
                "SBbul": sbbul,
            })
        except (ValueError, IndexError) as e:
            logger.warning("Skipping malformed line %d: %s", data_start + 1, e)

    all_df = pd.DataFrame(records)
    logger.info(
        "Parsed MRT file: %d data points across %d galaxies",
        len(all_df),
        all_df["galaxy_id"].nunique(),
    )

    # Split into per-galaxy DataFrames
    result = {}
    for gid, group in all_df.groupby("galaxy_id", sort=False):
        result[gid] = group[SPARC_COLUMNS].reset_index(drop=True)

    return result


def parse_sparc_metadata_mrt(filepath: str | Path) -> pd.DataFrame:
    """Parse the SPARC galaxy metadata MRT file (Lelli et al. 2016, Table 1).

    Whitespace-separated with 19 fields per line:
        [0] Galaxy, [1] T (Hubble type), [2] D (Mpc), [3] e_D, [4] f_D,
        [5] Inc (deg), [6] e_Inc, [7] L[3.6] (10^9 Lsun), [8] e_L,
        [9] Reff (kpc), [10] SBeff, [11] Rdisk (kpc), [12] SBdisk,
        [13] MHI (10^9 Msun), [14] RHI (kpc), [15] Vflat (km/s),
        [16] e_Vflat, [17] Q (quality flag), [18] References

    Returns:
        DataFrame with columns: galaxy_id, hubble_type, distance_mpc,
        inclination, luminosity_band_36, quality_flag, v_flat.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Metadata MRT file not found: {filepath}")

    with open(filepath, "r") as f:
        lines = f.readlines()

    # Find data start (after last '---' separator)
    data_start = 0
    for i, line in enumerate(lines):
        if line.strip().startswith("---"):
            data_start = i + 1

    records = []
    for line in lines[data_start:]:
        if len(line.strip()) == 0:
            continue
        try:
            fields = line.split()
            if len(fields) < 18:
                continue
            records.append({
                "galaxy_id": fields[0],
                "hubble_type": int(fields[1]),
                "distance_mpc": float(fields[2]),
                "inclination": float(fields[5]),
                "luminosity_band_36": float(fields[7]),  # 10^9 Lsun
                "quality_flag": int(fields[17]),
                "v_flat": float(fields[15]),
            })
        except (ValueError, IndexError) as e:
            logger.warning("Skipping metadata line: %s", e)

    df = pd.DataFrame(records)
    logger.info("Parsed metadata for %d galaxies", len(df))
    return df


def ingest_massmodels_mrt(
    massmodels_path: str | Path,
    metadata_path: Optional[str | Path] = None,
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    db_path: Optional[str] = None,
) -> list[str]:
    """Ingest all galaxies from the combined MassModels MRT file.

    Optionally also ingests metadata from the SPARC Table 1 MRT file.

    Args:
        massmodels_path: Path to MassModels_Lelli2016c.mrt.
        metadata_path: Optional path to SPARC_Lelli2016c.mrt for metadata.
        upsilon_disk: Mass-to-light ratio for disk. Default 0.5.
        upsilon_bulge: Mass-to-light ratio for bulge. Default 0.7.
        db_path: Optional database path override.

    Returns:
        List of galaxy_ids that were ingested.
    """
    # Parse rotation curve data
    galaxies_data = parse_massmodels_mrt(massmodels_path)

    # Parse metadata if provided
    metadata_df = None
    if metadata_path:
        metadata_df = parse_sparc_metadata_mrt(metadata_path)
        metadata_df = metadata_df.set_index("galaxy_id")

    # Initialize DB
    engine = init_db(db_path)
    session = get_session(engine)

    ingested = []
    try:
        for galaxy_id, raw_df in galaxies_data.items():
            # Compute baryonic velocity
            v_bary = compute_v_bary(
                raw_df["Vgas"].values,
                raw_df["Vdisk"].values,
                raw_df["Vbul"].values,
                upsilon_disk=upsilon_disk,
                upsilon_bulge=upsilon_bulge,
            )

            # Build profiles DataFrame
            profiles_df = pd.DataFrame({
                "radius_kpc": raw_df["Rad"],
                "v_obs": raw_df["Vobs"],
                "v_err": raw_df["errV"],
                "v_gas": raw_df["Vgas"],
                "v_disk": raw_df["Vdisk"],
                "v_bulge": raw_df["Vbul"],
                "v_baryon_total": v_bary,
            })

            # Build galaxy metadata kwargs
            meta_kwargs = {}
            if metadata_df is not None and galaxy_id in metadata_df.index:
                row = metadata_df.loc[galaxy_id]
                meta_kwargs = {
                    "distance_mpc": float(row["distance_mpc"]),
                    "inclination": float(row["inclination"]),
                    "luminosity_band_36": float(row["luminosity_band_36"]),
                    "quality_flag": int(row["quality_flag"]),
                }

            # Insert into DB
            insert_galaxy(session, galaxy_id, **meta_kwargs)
            insert_radial_profiles(session, galaxy_id, profiles_df)
            ingested.append(galaxy_id)

        logger.info(
            "Ingested %d galaxies from MRT (Upsilon_d=%.2f, Upsilon_b=%.2f)",
            len(ingested),
            upsilon_disk,
            upsilon_bulge,
        )
    finally:
        session.close()

    return ingested


# ---------------------------------------------------------------------------
# Corbelli M33 (Track 1)
# ---------------------------------------------------------------------------

# Flexible column mapping for VizieR tables — maps possible VizieR column
# names to our SPARC convention.  Extended as needed after inspecting actual
# catalog structure.
_VIZIER_COL_MAP = {
    # Radius
    "Rmaj": "Rad", "R": "Rad", "Rkpc": "Rad", "Rad": "Rad", "r": "Rad",
    # Observed velocity
    "Vrot": "Vobs", "Vobs": "Vobs", "V": "Vobs", "Vc": "Vobs",
    # Error
    "e_Vrot": "errV", "e_Vobs": "errV", "errV": "errV", "e_V": "errV",
    "e_Vc": "errV",
    # Gas
    "Vgas": "Vgas", "VHI": "Vgas", "Vg": "Vgas",
    # Disk / stellar
    "Vdisk": "Vdisk", "Vd": "Vdisk", "Vstar": "Vdisk", "Vlum": "Vdisk",
    "Vs": "Vdisk",
    # Bulge
    "Vbul": "Vbul", "Vbulge": "Vbul", "Vb": "Vbul",
}


def _map_vizier_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename VizieR columns to SPARC convention using flexible mapping."""
    renamed = {}
    for viz_col, sparc_col in _VIZIER_COL_MAP.items():
        if viz_col in df.columns and sparc_col not in renamed.values():
            renamed[viz_col] = sparc_col
    return df.rename(columns=renamed)


def parse_corbelli_vizier(df_vizier: pd.DataFrame) -> pd.DataFrame:
    """Parse a VizieR table from Corbelli et al. into SPARC-compatible format.

    Applies flexible column name mapping and fills missing columns
    (M33 has negligible bulge per Corbelli 2000).

    Args:
        df_vizier: Raw DataFrame from a VizieR query.

    Returns:
        DataFrame with SPARC-compatible columns.

    Raises:
        ValueError: If required columns (Rad, Vobs, errV) cannot be mapped.
    """
    df = _map_vizier_columns(df_vizier.copy())

    # M33 has no bulge — fill if missing
    if "Vbul" not in df.columns:
        df["Vbul"] = 0.0

    # Surface brightness columns aren't critical for fitting
    if "SBdisk" not in df.columns:
        df["SBdisk"] = 0.0
    if "SBbul" not in df.columns:
        df["SBbul"] = 0.0

    # Fill Vgas/Vdisk with 0 and warn if missing
    for col in ("Vgas", "Vdisk"):
        if col not in df.columns:
            logger.warning(
                "Column %s not found in VizieR data — setting to 0.0", col
            )
            df[col] = 0.0

    # Validate required columns
    for col in ("Rad", "Vobs", "errV"):
        if col not in df.columns:
            raise ValueError(
                f"Could not map required column '{col}' from VizieR data. "
                f"Available columns after mapping: {list(df.columns)}"
            )

    return df[SPARC_COLUMNS]


# M33 physical parameters from Corbelli et al. (2014)
M33_DISTANCE_MPC = 0.84
M33_INCLINATION_DEG = 52.0


def load_m33_corbelli2014_data(
    csv_path: Optional[str | Path] = None,
) -> pd.DataFrame:
    """Load M33 data from Corbelli et al. (2014) Table 1 and compute velocities.

    Reads the extracted CSV containing rotation curve, HI gas surface density,
    and stellar mass surface density. Converts surface densities to velocity
    contributions using the thin-disk gravitational calculation (Casertano 1983).

    Args:
        csv_path: Path to the extracted CSV file. Defaults to
            data/extracted/corbelli2014_table1.csv relative to project root.

    Returns:
        DataFrame with SPARC-format columns (Rad, Vobs, errV, Vgas, Vdisk,
        Vbul, SBdisk, SBbul) with 57 radial bins from R=0.24 to R=22.72 kpc.

    Raises:
        FileNotFoundError: If the CSV file is not found.
    """
    if csv_path is None:
        csv_path = get_project_root() / "data" / "extracted" / "corbelli2014_table1.csv"
    csv_path = Path(csv_path)

    if not csv_path.exists():
        raise FileNotFoundError(
            f"Corbelli 2014 Table 1 CSV not found: {csv_path}. "
            "Expected at data/extracted/corbelli2014_table1.csv"
        )

    df = pd.read_csv(csv_path, comment="#")
    logger.info(
        "Loaded Corbelli 2014 Table 1: %d data points (R = %.2f – %.2f kpc)",
        len(df), df["R_kpc"].min(), df["R_kpc"].max(),
    )

    # Convert surface densities to velocity contributions
    radii = df["R_kpc"].values
    sigma_hi = df["Sigma_HI_Msun_pc2"].values
    sigma_star = df["Sigma_star_Msun_pc2"].values

    v_gas = circular_velocity_thin_disk(
        radii, radii, sigma_hi, helium_factor=HELIUM_FACTOR,
    )
    v_disk = circular_velocity_thin_disk(
        radii, radii, sigma_star, helium_factor=1.0,
    )

    logger.info(
        "Computed velocity components: V_gas peak=%.1f km/s, V_disk peak=%.1f km/s",
        np.max(v_gas), np.max(v_disk),
    )

    n = len(df)
    return pd.DataFrame({
        "Rad": radii,
        "Vobs": df["V_r_kms"].values,
        "errV": df["sigma_V_kms"].values,
        "Vgas": v_gas,
        "Vdisk": v_disk,
        "Vbul": np.zeros(n),
        "SBdisk": np.zeros(n),
        "SBbul": np.zeros(n),
    })


def create_m33_manual_data() -> pd.DataFrame:
    """Create approximate M33 rotation curve from Corbelli & Salucci (2000).

    LEGACY: This function provides approximate, manually-transcribed data
    suitable for pipeline validation only. Prefer load_m33_corbelli2014_data()
    for actual analysis.

    Key parameters from the paper:
      - Distance: D = 0.7 Mpc (5 arcmin = 1 kpc)
      - Disk scale length: R_d = 1.2 kpc
      - Best-fit M_d/L_B: 0.8 ± 0.2
      - Bulge: negligible

    Returns:
        DataFrame with SPARC-format columns (~20 data points, 0.5-16 kpc).
    """
    logger.warning(
        "Using manually transcribed M33 data from Corbelli & Salucci (2000). "
        "Values are approximate — suitable for pipeline validation only."
    )
    data = {
        #       Rad   Vobs  errV  Vgas  Vdisk  Vbul  SBdisk  SBbul
        # From Figure 2b (binned rotation curve) and Section 4 (gas/stellar)
        "Rad":    [ 0.5,  1.0,  1.5,  2.0,  2.5,  3.0,  4.0,  5.0,
                    6.0,  7.0,  8.0,  9.0, 10.0, 11.0, 12.0, 13.0,
                   14.0, 15.0, 16.0],
        "Vobs":   [  25,   50,   75,   95,  107,  115,  118,  115,
                    113,  115,  120,  122,  125,  128,  130,  128,
                    127,  126,  125],
        "errV":   [ 5.0,  4.0,  3.5,  3.0,  3.0,  3.0,  3.0,  3.5,
                    3.5,  3.5,  4.0,  4.0,  5.0,  5.0,  6.0,  6.0,
                    7.0,  7.0,  8.0],
        "Vgas":   [ 5.0, 10.0, 16.0, 22.0, 28.0, 33.0, 38.0, 40.0,
                   40.0, 39.0, 37.0, 34.0, 30.0, 26.0, 22.0, 18.0,
                   15.0, 12.0, 10.0],
        "Vdisk":  [20.0, 40.0, 55.0, 65.0, 72.0, 76.0, 78.0, 76.0,
                   72.0, 68.0, 64.0, 60.0, 56.0, 52.0, 48.0, 44.0,
                   40.0, 37.0, 34.0],
        "Vbul":   [0.0] * 19,
        "SBdisk": [0.0] * 19,
        "SBbul":  [0.0] * 19,
    }
    return pd.DataFrame(data)


def ingest_m33(
    use_vizier: bool = False,
    use_legacy_manual: bool = False,
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    db_path: Optional[str] = None,
) -> str:
    """Ingest M33 rotation curve data from Corbelli sources.

    Default: loads Corbelli et al. (2014) Table 1 from the extracted CSV,
    converting surface densities to velocity contributions via thin-disk
    gravitational calculation.

    Args:
        use_vizier: If True, try VizieR queries first (Corbelli 2000/2014),
            falling back to the CSV data.
        use_legacy_manual: If True, use the old approximate manually-
            transcribed data from Corbelli & Salucci (2000).
        upsilon_disk: Mass-to-light ratio for disk.
        upsilon_bulge: Mass-to-light ratio for bulge.
        db_path: Optional database path override.

    Returns:
        Galaxy ID ("M33").
    """
    raw_df = None
    data_source = None
    distance_mpc = M33_DISTANCE_MPC
    inclination = M33_INCLINATION_DEG

    if use_legacy_manual:
        raw_df = create_m33_manual_data()
        data_source = "Corbelli2000_manual"
        distance_mpc = 0.7
        inclination = 54.0
    elif use_vizier:
        from src.vizier_query import query_corbelli_2000, query_corbelli_2014

        # Try Corbelli & Salucci (2000) first
        df_viz = query_corbelli_2000()
        if df_viz is not None:
            try:
                raw_df = parse_corbelli_vizier(df_viz)
                data_source = "Corbelli2000_VizieR"
                logger.info("Successfully parsed Corbelli (2000) from VizieR")
            except ValueError as e:
                logger.warning("Corbelli (2000) parse failed: %s", e)

        # Try Corbelli et al. (2014) as backup
        if raw_df is None:
            df_viz = query_corbelli_2014()
            if df_viz is not None:
                try:
                    raw_df = parse_corbelli_vizier(df_viz)
                    data_source = "Corbelli2014_VizieR"
                    logger.info("Successfully parsed Corbelli (2014) from VizieR")
                except ValueError as e:
                    logger.warning("Corbelli (2014) parse failed: %s", e)

    # Default: Corbelli 2014 Table 1 from extracted CSV
    if raw_df is None:
        raw_df = load_m33_corbelli2014_data()
        data_source = "Corbelli2014_Table1"

    # Compute V_bary
    v_bary = compute_v_bary(
        raw_df["Vgas"].values,
        raw_df["Vdisk"].values,
        raw_df["Vbul"].values,
        upsilon_disk=upsilon_disk,
        upsilon_bulge=upsilon_bulge,
    )

    profiles_df = pd.DataFrame({
        "radius_kpc": raw_df["Rad"],
        "v_obs": raw_df["Vobs"],
        "v_err": raw_df["errV"],
        "v_gas": raw_df["Vgas"],
        "v_disk": raw_df["Vdisk"],
        "v_bulge": raw_df["Vbul"],
        "v_baryon_total": v_bary,
    })

    # Insert into DB
    engine = init_db(db_path)
    session = get_session(engine)
    try:
        insert_galaxy(
            session, "M33",
            distance_mpc=distance_mpc,
            inclination=inclination,
            data_source=data_source,
        )
        n_rows = insert_radial_profiles(session, "M33", profiles_df)
        logger.info(
            "Ingested M33: %d profiles from %s (Upsilon_d=%.2f, Upsilon_b=%.2f)",
            n_rows, data_source, upsilon_disk, upsilon_bulge,
        )
    finally:
        session.close()

    return "M33"


# ---------------------------------------------------------------------------
# THINGS Survey (Track 2)
# ---------------------------------------------------------------------------

def parse_things_galaxy(df_vizier: pd.DataFrame, galaxy_name: str) -> pd.DataFrame:
    """Parse a single THINGS galaxy's VizieR data into SPARC format.

    Applies flexible column name mapping. Most THINGS galaxies are
    late-type spirals with negligible bulge.

    Args:
        df_vizier: Raw DataFrame for one galaxy from VizieR.
        galaxy_name: Galaxy identifier (for error messages).

    Returns:
        DataFrame with SPARC-compatible columns.

    Raises:
        ValueError: If required columns cannot be mapped.
    """
    df = _map_vizier_columns(df_vizier.copy())

    # Most THINGS galaxies have no bulge data
    if "Vbul" not in df.columns:
        df["Vbul"] = 0.0

    # Surface brightness not provided by THINGS
    df["SBdisk"] = 0.0
    df["SBbul"] = 0.0

    # Fill missing gas/disk with 0 and warn
    for col in ("Vgas", "Vdisk"):
        if col not in df.columns:
            logger.warning(
                "THINGS galaxy %s: column %s not found — setting to 0.0",
                galaxy_name, col,
            )
            df[col] = 0.0

    # Validate required columns
    missing = [c for c in ("Rad", "Vobs", "errV") if c not in df.columns]
    if missing:
        raise ValueError(
            f"THINGS galaxy {galaxy_name}: cannot map required columns {missing}. "
            f"Available after mapping: {list(df.columns)}"
        )

    return df[SPARC_COLUMNS]


def ingest_things_catalog(
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    skip_existing: bool = True,
    db_path: Optional[str] = None,
) -> list[str]:
    """Ingest THINGS survey rotation curves from VizieR.

    Queries the de Blok et al. (2008) catalog (J/AJ/136/2648) and ingests
    all available galaxy mass models.

    Args:
        upsilon_disk: Mass-to-light ratio for disk.
        upsilon_bulge: Mass-to-light ratio for bulge.
        skip_existing: If True, skip galaxies already in the database.
        db_path: Optional database path override.

    Returns:
        List of galaxy_ids successfully ingested.
    """
    from src.vizier_query import query_things_rotation_curves

    galaxies_data = query_things_rotation_curves()
    if galaxies_data is None:
        logger.error("THINGS catalog not available from VizieR")
        return []

    engine = init_db(db_path)
    session = get_session(engine)

    ingested = []
    skipped = []
    failed = []
    try:
        for galaxy_name, df_viz in galaxies_data.items():
            galaxy_id = galaxy_name

            # Check for existing galaxy if skip_existing
            if skip_existing:
                from src.database import Galaxy
                existing = (
                    session.query(Galaxy)
                    .filter(Galaxy.galaxy_id == galaxy_id)
                    .first()
                )
                if existing is not None:
                    skipped.append(galaxy_id)
                    continue

            try:
                raw_df = parse_things_galaxy(df_viz, galaxy_name)

                v_bary = compute_v_bary(
                    raw_df["Vgas"].values,
                    raw_df["Vdisk"].values,
                    raw_df["Vbul"].values,
                    upsilon_disk=upsilon_disk,
                    upsilon_bulge=upsilon_bulge,
                )

                profiles_df = pd.DataFrame({
                    "radius_kpc": raw_df["Rad"],
                    "v_obs": raw_df["Vobs"],
                    "v_err": raw_df["errV"],
                    "v_gas": raw_df["Vgas"],
                    "v_disk": raw_df["Vdisk"],
                    "v_bulge": raw_df["Vbul"],
                    "v_baryon_total": v_bary,
                })

                insert_galaxy(session, galaxy_id, data_source="THINGS")
                insert_radial_profiles(session, galaxy_id, profiles_df)
                ingested.append(galaxy_id)

            except Exception as e:
                logger.error("Failed to ingest THINGS galaxy %s: %s", galaxy_name, e)
                failed.append(galaxy_name)

        logger.info(
            "THINGS ingestion complete: %d ingested, %d skipped (existing), %d failed",
            len(ingested), len(skipped), len(failed),
        )
    finally:
        session.close()

    return ingested


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ingest galaxy rotation curve data")
    # SPARC sources
    parser.add_argument("--file", type=str, help="Path to single _rotmod.dat file")
    parser.add_argument("--dir", type=str, help="Path to directory of _rotmod.dat files")
    parser.add_argument("--mrt", type=str, help="Path to MassModels MRT file (all galaxies)")
    parser.add_argument("--metadata", type=str, help="Path to SPARC metadata MRT file")
    parser.add_argument("--name", type=str, help="Galaxy name override (single file only)")
    # Corbelli M33
    parser.add_argument("--m33", action="store_true",
                        help="Ingest M33 from Corbelli 2014 Table 1 (default)")
    parser.add_argument("--m33-vizier", action="store_true",
                        help="Ingest M33 trying VizieR first, falling back to Table 1")
    parser.add_argument("--m33-manual", action="store_true",
                        help="Ingest M33 using legacy manually transcribed data")
    # THINGS survey
    parser.add_argument("--things", action="store_true",
                        help="Ingest THINGS survey galaxies from VizieR")
    # Shared options
    parser.add_argument("--upsilon-disk", type=float, default=0.5)
    parser.add_argument("--upsilon-bulge", type=float, default=0.7)
    args = parser.parse_args()

    if args.m33 or args.m33_vizier or args.m33_manual:
        ingest_m33(
            use_vizier=args.m33_vizier,
            use_legacy_manual=args.m33_manual,
            upsilon_disk=args.upsilon_disk,
            upsilon_bulge=args.upsilon_bulge,
        )
    elif args.things:
        ingest_things_catalog(
            upsilon_disk=args.upsilon_disk,
            upsilon_bulge=args.upsilon_bulge,
        )
    elif args.mrt:
        ingest_massmodels_mrt(
            args.mrt,
            metadata_path=args.metadata,
            upsilon_disk=args.upsilon_disk,
            upsilon_bulge=args.upsilon_bulge,
        )
    elif args.file:
        ingest_sparc_file(
            args.file,
            galaxy_name=args.name,
            upsilon_disk=args.upsilon_disk,
            upsilon_bulge=args.upsilon_bulge,
        )
    elif args.dir:
        ingest_sparc_directory(
            args.dir,
            upsilon_disk=args.upsilon_disk,
            upsilon_bulge=args.upsilon_bulge,
        )
    else:
        parser.print_help()
