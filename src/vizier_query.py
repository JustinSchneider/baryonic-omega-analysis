"""VizieR query utilities for fetching rotation curve data from external surveys.

Provides functions to query specific catalogs:
  - Corbelli & Salucci (2000) M33 rotation curve
  - Corbelli et al. (2014) updated M33 mass model
  - de Blok et al. (2008) THINGS survey rotation curves
"""

from concurrent.futures import ThreadPoolExecutor, TimeoutError as FuturesTimeoutError
from typing import Optional

import pandas as pd
from astroquery.vizier import Vizier

from src.utils import setup_logger

logger = setup_logger(__name__)


DEFAULT_TIMEOUT_S = 15


def search_vizier_catalog(
    catalog_id: str, table_index: int = 0, timeout: int = DEFAULT_TIMEOUT_S
) -> Optional[pd.DataFrame]:
    """Query a VizieR catalog by ID and return as a DataFrame.

    Args:
        catalog_id: VizieR catalog identifier (e.g. "J/MNRAS/311/441").
        table_index: Which table to return if catalog has multiple (default 0).
        timeout: Request timeout in seconds (default 15).

    Returns:
        DataFrame of the requested table, or None if catalog not found.
    """
    v = Vizier(row_limit=-1, timeout=timeout)
    try:
        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(v.get_catalogs, catalog_id)
            tables = future.result(timeout=timeout)
    except FuturesTimeoutError:
        logger.error("VizieR query timed out for %s after %ds", catalog_id, timeout)
        return None
    except Exception as e:
        logger.error("VizieR query failed for %s: %s", catalog_id, e)
        return None

    if tables is None or len(tables) == 0:
        logger.warning("Catalog %s not found in VizieR", catalog_id)
        return None

    if table_index >= len(tables):
        logger.error(
            "Requested table_index %d but catalog %s only has %d table(s)",
            table_index,
            catalog_id,
            len(tables),
        )
        return None

    logger.info(
        "VizieR catalog %s: %d table(s), columns in table[%d]: %s",
        catalog_id,
        len(tables),
        table_index,
        list(tables[table_index].colnames),
    )

    return tables[table_index].to_pandas()


def list_vizier_tables(catalog_id: str) -> Optional[list[str]]:
    """List all table names/descriptions in a VizieR catalog.

    Useful for discovering which table contains rotation curve data
    when a catalog has multiple tables.

    Returns:
        List of table descriptions, or None if catalog not found.
    """
    v = Vizier(row_limit=1, timeout=DEFAULT_TIMEOUT_S)
    try:
        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(v.get_catalogs, catalog_id)
            tables = future.result(timeout=DEFAULT_TIMEOUT_S)
    except FuturesTimeoutError:
        logger.error("VizieR query timed out for %s after %ds", catalog_id, DEFAULT_TIMEOUT_S)
        return None
    except Exception as e:
        logger.error("VizieR query failed for %s: %s", catalog_id, e)
        return None

    if tables is None or len(tables) == 0:
        return None

    descriptions = []
    for i, table in enumerate(tables):
        desc = table.meta.get("description", f"Table {i}")
        descriptions.append(f"[{i}] {desc} — columns: {list(table.colnames)}")
    return descriptions


def query_corbelli_2000() -> Optional[pd.DataFrame]:
    """Query Corbelli & Salucci (2000) M33 data from VizieR.

    Catalog: J/MNRAS/311/441
    Paper: "The Extended Rotation Curve and the Dark Matter Halo of M33"

    Returns:
        DataFrame with M33 rotation curve data, or None if not available.
    """
    logger.info("Querying VizieR for Corbelli & Salucci (2000) — J/MNRAS/311/441")
    return search_vizier_catalog("J/MNRAS/311/441")


def query_corbelli_2014() -> Optional[pd.DataFrame]:
    """Query Corbelli et al. (2014) M33 data from VizieR.

    Catalog: J/A+A/572/A23
    Paper: "Dynamical signatures of a ΛCDM-halo and the distribution
           of the baryons in M33"

    Returns:
        DataFrame with M33 rotation curve data, or None if not available.
    """
    logger.info("Querying VizieR for Corbelli et al. (2014) — J/A+A/572/A23")
    return search_vizier_catalog("J/A+A/572/A23")


def query_things_rotation_curves() -> Optional[dict[str, pd.DataFrame]]:
    """Query THINGS survey rotation curves from VizieR.

    Catalog: J/AJ/136/2648
    Paper: de Blok et al. (2008) "High-Resolution Rotation Curves and
           Galaxy Mass Models from THINGS"

    Returns:
        Dict mapping galaxy name -> DataFrame of rotation curve data,
        or None if catalog not available.
    """
    logger.info("Querying VizieR for THINGS (de Blok et al. 2008) — J/AJ/136/2648")

    # THINGS catalog may have multiple tables; list them first to find
    # the rotation curve table
    table_info = list_vizier_tables("J/AJ/136/2648")
    if table_info is None:
        logger.warning("THINGS catalog J/AJ/136/2648 not found in VizieR")
        return None

    for desc in table_info:
        logger.info("  THINGS table: %s", desc)

    # Try the first table; may need adjustment after inspecting actual structure
    df = search_vizier_catalog("J/AJ/136/2648", table_index=0)
    if df is None:
        return None

    # Identify the galaxy name column — common names: "Galaxy", "Name", "Gal"
    name_col = None
    for candidate in ["Galaxy", "Name", "Gal", "galaxy", "name"]:
        if candidate in df.columns:
            name_col = candidate
            break

    if name_col is None:
        logger.error(
            "Cannot identify galaxy name column in THINGS data. "
            "Available columns: %s",
            list(df.columns),
        )
        return None

    result = {}
    for gal_name, group in df.groupby(name_col, sort=False):
        result[str(gal_name).strip()] = group.reset_index(drop=True)

    logger.info("Retrieved %d THINGS galaxies: %s", len(result), list(result.keys()))
    return result
