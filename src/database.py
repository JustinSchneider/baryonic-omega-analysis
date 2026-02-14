"""SQLite database schema, connection management, and query helpers.

Uses SQLAlchemy 2.0 ORM with DeclarativeBase for the three-table schema:
  - Galaxies: metadata (distance, inclination, luminosity, quality flag)
  - RadialProfiles: per-radius velocity components from SPARC data
  - OmegaFits: omega fitting results with full diagnostics
"""

import argparse
from datetime import datetime
from typing import Optional

import pandas as pd
from sqlalchemy import (
    Float,
    Integer,
    String,
    DateTime,
    Boolean,
    ForeignKey,
    create_engine,
    delete,
)
from sqlalchemy.orm import (
    DeclarativeBase,
    Mapped,
    Session,
    mapped_column,
    relationship,
    sessionmaker,
)

from src.utils import get_db_path, setup_logger

logger = setup_logger(__name__)


# ---------------------------------------------------------------------------
# ORM Models
# ---------------------------------------------------------------------------

class Base(DeclarativeBase):
    pass


class Galaxy(Base):
    __tablename__ = "galaxies"

    galaxy_id: Mapped[str] = mapped_column(String, primary_key=True)
    distance_mpc: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    inclination: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    position_angle: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    luminosity_band_36: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    quality_flag: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    data_source: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    profiles: Mapped[list["RadialProfile"]] = relationship(
        back_populates="galaxy", cascade="all, delete-orphan"
    )
    fits: Mapped[list["OmegaFit"]] = relationship(
        back_populates="galaxy", cascade="all, delete-orphan"
    )


class RadialProfile(Base):
    __tablename__ = "radial_profiles"

    profile_id: Mapped[int] = mapped_column(
        Integer, primary_key=True, autoincrement=True
    )
    galaxy_id: Mapped[str] = mapped_column(ForeignKey("galaxies.galaxy_id"))
    radius_kpc: Mapped[float] = mapped_column(Float)
    v_obs: Mapped[float] = mapped_column(Float)
    v_err: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    v_gas: Mapped[float] = mapped_column(Float)
    v_disk: Mapped[float] = mapped_column(Float)
    v_bulge: Mapped[float] = mapped_column(Float)
    v_baryon_total: Mapped[Optional[float]] = mapped_column(Float, nullable=True)

    galaxy: Mapped["Galaxy"] = relationship(back_populates="profiles")


class OmegaFit(Base):
    __tablename__ = "omega_fits"

    fit_id: Mapped[int] = mapped_column(
        Integer, primary_key=True, autoincrement=True
    )
    galaxy_id: Mapped[str] = mapped_column(ForeignKey("galaxies.galaxy_id"))
    omega_value: Mapped[float] = mapped_column(Float)
    omega_uncertainty: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    chi_squared: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    reduced_chi_squared: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    residuals_rmse: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    method_version: Mapped[str] = mapped_column(String, default="v1_fixed_ML")
    upsilon_disk: Mapped[float] = mapped_column(Float, default=0.5)
    upsilon_bulge: Mapped[float] = mapped_column(Float, default=0.7)
    n_points: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    converged: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)
    flag_v_obs_lt_v_bary: Mapped[Optional[bool]] = mapped_column(
        Boolean, nullable=True
    )
    timestamp: Mapped[datetime] = mapped_column(
        DateTime, default=datetime.utcnow
    )

    galaxy: Mapped["Galaxy"] = relationship(back_populates="fits")


# ---------------------------------------------------------------------------
# Engine / Session helpers
# ---------------------------------------------------------------------------

def get_engine(db_path: str | None = None):
    """Create a SQLAlchemy engine for the galaxy_dynamics database."""
    path = db_path or str(get_db_path())
    engine = create_engine(f"sqlite:///{path}", echo=False)
    return engine


def init_db(db_path: str | None = None):
    """Create all tables if they don't exist. Safe to call repeatedly."""
    engine = get_engine(db_path)
    Base.metadata.create_all(engine)
    logger.info("Database initialized at %s", db_path or get_db_path())
    return engine


def get_session(engine=None) -> Session:
    """Create and return a new SQLAlchemy Session."""
    if engine is None:
        engine = get_engine()
    factory = sessionmaker(bind=engine)
    return factory()


# ---------------------------------------------------------------------------
# Insert / Query helpers
# ---------------------------------------------------------------------------

def insert_galaxy(session: Session, galaxy_id: str, **metadata) -> Galaxy:
    """Insert or update a galaxy record using merge (upsert semantics)."""
    galaxy = Galaxy(galaxy_id=galaxy_id, **metadata)
    galaxy = session.merge(galaxy)
    session.commit()
    logger.info("Upserted galaxy: %s", galaxy_id)
    return galaxy


def insert_radial_profiles(
    session: Session, galaxy_id: str, df: pd.DataFrame
) -> int:
    """Bulk-insert radial profile rows from a DataFrame.

    Clears existing profiles for this galaxy_id first (replace semantics).

    Args:
        session: Active SQLAlchemy session.
        galaxy_id: Galaxy identifier.
        df: DataFrame with columns: radius_kpc, v_obs, v_err, v_gas,
            v_disk, v_bulge, and optionally v_baryon_total.

    Returns:
        Number of rows inserted.
    """
    # Delete existing profiles for this galaxy
    session.execute(
        delete(RadialProfile).where(RadialProfile.galaxy_id == galaxy_id)
    )

    rows = []
    for _, row in df.iterrows():
        profile = RadialProfile(
            galaxy_id=galaxy_id,
            radius_kpc=float(row["radius_kpc"]),
            v_obs=float(row["v_obs"]),
            v_err=float(row.get("v_err", 0.0)) if pd.notna(row.get("v_err")) else None,
            v_gas=float(row["v_gas"]),
            v_disk=float(row["v_disk"]),
            v_bulge=float(row["v_bulge"]),
            v_baryon_total=(
                float(row["v_baryon_total"])
                if "v_baryon_total" in row and pd.notna(row.get("v_baryon_total"))
                else None
            ),
        )
        rows.append(profile)

    session.add_all(rows)
    session.commit()
    logger.info("Inserted %d radial profiles for %s", len(rows), galaxy_id)
    return len(rows)


def query_profiles_as_dataframe(session: Session, galaxy_id: str) -> pd.DataFrame:
    """Retrieve all radial profiles for a galaxy as a Pandas DataFrame."""
    profiles = (
        session.query(RadialProfile)
        .filter(RadialProfile.galaxy_id == galaxy_id)
        .order_by(RadialProfile.radius_kpc)
        .all()
    )
    if not profiles:
        logger.warning("No profiles found for galaxy: %s", galaxy_id)
        return pd.DataFrame()

    data = [
        {
            "radius_kpc": p.radius_kpc,
            "v_obs": p.v_obs,
            "v_err": p.v_err,
            "v_gas": p.v_gas,
            "v_disk": p.v_disk,
            "v_bulge": p.v_bulge,
            "v_baryon_total": p.v_baryon_total,
        }
        for p in profiles
    ]
    return pd.DataFrame(data)


def insert_omega_fit(session: Session, fit_result: dict) -> OmegaFit:
    """Insert a new fit result into the OmegaFits table.

    Args:
        fit_result: dict with keys matching OmegaFit column names.
    """
    fit = OmegaFit(**fit_result)
    session.add(fit)
    session.commit()
    logger.info(
        "Stored omega fit for %s: omega=%.4f +/- %.4f (chi2=%.2f)",
        fit.galaxy_id,
        fit.omega_value,
        fit.omega_uncertainty or 0.0,
        fit.chi_squared or 0.0,
    )
    return fit


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Galaxy dynamics database management")
    parser.add_argument(
        "--init", action="store_true", help="Initialize database schema"
    )
    args = parser.parse_args()
    if args.init:
        init_db()
