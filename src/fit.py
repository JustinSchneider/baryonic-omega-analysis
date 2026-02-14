"""CLI fitting script: run omega fits on galaxies and generate plots.

Usage:
    python src/fit.py --galaxy NGC5055 --plot
    python src/fit.py --all --quality 1 --plot
"""

import argparse
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from src.database import (
    Galaxy,
    get_engine,
    get_session,
    init_db,
    insert_omega_fit,
    query_profiles_as_dataframe,
)
from src.physics import OmegaFitResult, compute_v_bary, fit_omega
from src.utils import get_project_root, setup_logger

logger = setup_logger(__name__)


def run_fit_for_galaxy(
    galaxy_id: str,
    method_version: str = "v1_fixed_ML",
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    db_path: Optional[str] = None,
    plot: bool = False,
    output_dir: Optional[str] = None,
) -> Optional[OmegaFitResult]:
    """Execute the full fitting pipeline for a single galaxy.

    Steps:
        1. Query RadialProfiles from database.
        2. Compute V_bary from component velocities.
        3. Fit omega using physics.fit_omega().
        4. Store result in OmegaFits table.
        5. Optionally generate decomposition plot.

    Returns:
        OmegaFitResult, or None if no data found.
    """
    engine = init_db(db_path)
    session = get_session(engine)

    try:
        df = query_profiles_as_dataframe(session, galaxy_id)
        if df.empty:
            logger.error("No data found for galaxy: %s", galaxy_id)
            return None

        # Extract arrays
        radius = df["radius_kpc"].values
        v_obs = df["v_obs"].values
        v_err = df["v_err"].values
        v_gas = df["v_gas"].values
        v_disk = df["v_disk"].values
        v_bulge = df["v_bulge"].values

        # Compute baryonic velocity
        v_bary = compute_v_bary(
            v_gas, v_disk, v_bulge,
            upsilon_disk=upsilon_disk,
            upsilon_bulge=upsilon_bulge,
        )

        # Fit omega
        result = fit_omega(
            radius, v_obs, v_err, v_bary,
            galaxy_id=galaxy_id,
            method_version=method_version,
            upsilon_disk=upsilon_disk,
            upsilon_bulge=upsilon_bulge,
        )

        # Store in database
        insert_omega_fit(session, result.to_dict())

        # Generate plot if requested
        if plot:
            out_dir = output_dir or str(get_project_root() / "results" / "figures")
            generate_decomposition_plot(
                galaxy_id=galaxy_id,
                radius=radius,
                v_obs=v_obs,
                v_err=v_err,
                v_gas=v_gas,
                v_disk=v_disk,
                v_bulge=v_bulge,
                v_bary=v_bary,
                v_model=result.v_model,
                omega_value=result.omega_value,
                chi_squared=result.reduced_chi_squared,
                upsilon_disk=upsilon_disk,
                upsilon_bulge=upsilon_bulge,
                output_path=str(Path(out_dir) / f"{galaxy_id}_decomposition.png"),
            )

        return result

    finally:
        session.close()


def generate_decomposition_plot(
    galaxy_id: str,
    radius: np.ndarray,
    v_obs: np.ndarray,
    v_err: np.ndarray,
    v_gas: np.ndarray,
    v_disk: np.ndarray,
    v_bulge: np.ndarray,
    v_bary: np.ndarray,
    v_model: np.ndarray,
    omega_value: float,
    chi_squared: float,
    upsilon_disk: float = 0.5,
    upsilon_bulge: float = 0.7,
    output_path: Optional[str] = None,
):
    """Generate a SPARC-style rotation curve decomposition plot.

    Visual style follows Lelli et al. (2016) Figure 5 conventions.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Observed data with error bars
    ax.errorbar(
        radius, v_obs, yerr=v_err,
        fmt="ko", markersize=4, capsize=2, label=r"$V_\mathrm{obs}$",
    )

    # Scaled component velocities for display
    # Show the absolute magnitude of each component contribution
    v_gas_display = np.sqrt(np.abs(v_gas) * np.abs(v_gas))  # |V_gas|
    v_disk_display = np.sqrt(upsilon_disk) * np.abs(v_disk)
    v_bulge_display = np.sqrt(upsilon_bulge) * np.abs(v_bulge)

    ax.plot(radius, v_gas_display, "g--", linewidth=1.5, label=r"$V_\mathrm{gas}$")
    ax.plot(radius, v_disk_display, "r--", linewidth=1.5, label=r"$V_\mathrm{disk}$")
    if np.any(v_bulge != 0):
        ax.plot(
            radius, v_bulge_display, "m-.",
            linewidth=1.5, label=r"$V_\mathrm{bulge}$",
        )

    # Total baryonic
    ax.plot(radius, v_bary, "b-", linewidth=2, label=r"$V_\mathrm{bary}$")

    # Model with omega correction
    ax.plot(
        radius, v_model, color="orange", linewidth=2.5,
        label=rf"$V_\mathrm{{bary}} + \omega R$",
    )

    # Annotation
    ax.text(
        0.95, 0.05,
        (
            rf"$\omega = {omega_value:.3f}$ km/s/kpc"
            "\n"
            rf"$\chi^2_\nu = {chi_squared:.2f}$"
            "\n"
            rf"$\Upsilon_d = {upsilon_disk:.1f}$, $\Upsilon_b = {upsilon_bulge:.1f}$"
        ),
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8),
    )

    ax.set_xlabel("Radius (kpc)", fontsize=12)
    ax.set_ylabel("Velocity (km/s)", fontsize=12)
    ax.set_title(f"{galaxy_id} — Rotation Curve Decomposition", fontsize=14)
    ax.legend(loc="upper left", fontsize=10)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        logger.info("Plot saved: %s", output_path)

    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Fit omega parameter to galaxy rotation curves"
    )
    parser.add_argument("--galaxy", type=str, help="Single galaxy ID to fit")
    parser.add_argument("--all", action="store_true", help="Fit all galaxies in DB")
    parser.add_argument(
        "--quality", type=int, default=None,
        help="Only fit galaxies with this quality flag (1, 2, or 3)",
    )
    parser.add_argument(
        "--method", type=str, default="v1_fixed_ML",
        help="Method version string for tracking",
    )
    parser.add_argument("--upsilon-disk", type=float, default=0.5)
    parser.add_argument("--upsilon-bulge", type=float, default=0.7)
    parser.add_argument("--plot", action="store_true", help="Generate decomposition plots")
    parser.add_argument("--output-dir", type=str, default=None)
    args = parser.parse_args()

    if args.galaxy:
        result = run_fit_for_galaxy(
            args.galaxy,
            method_version=args.method,
            upsilon_disk=args.upsilon_disk,
            upsilon_bulge=args.upsilon_bulge,
            plot=args.plot,
            output_dir=args.output_dir,
        )
        if result and result.converged:
            print(
                f"{result.galaxy_id}: omega={result.omega_value:.4f} "
                f"+/- {result.omega_uncertainty:.4f} km/s/kpc  "
                f"chi2_r={result.reduced_chi_squared:.2f}"
            )

    elif args.all:
        engine = init_db()
        session = get_session(engine)
        query = session.query(Galaxy.galaxy_id)
        if args.quality is not None:
            query = query.filter(Galaxy.quality_flag == args.quality)
        galaxy_ids = [row[0] for row in query.all()]
        session.close()

        logger.info("Fitting %d galaxies", len(galaxy_ids))
        converged_count = 0
        for gid in galaxy_ids:
            result = run_fit_for_galaxy(
                gid,
                method_version=args.method,
                upsilon_disk=args.upsilon_disk,
                upsilon_bulge=args.upsilon_bulge,
                plot=args.plot,
                output_dir=args.output_dir,
            )
            if result and result.converged:
                converged_count += 1

        logger.info(
            "Converged: %d / %d (%.1f%%)",
            converged_count,
            len(galaxy_ids),
            100 * converged_count / max(len(galaxy_ids), 1),
        )

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
