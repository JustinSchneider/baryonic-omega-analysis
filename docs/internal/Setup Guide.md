This is an ambitious and scientifically grounded upgrade to the EPS research project. Moving from Flynn & Cannaliato's "point-mass" approximation to Corbelli & Salucci's "baryonic mass decomposition" will significantly strengthen your results and likely reduce the scatter you see in your correlations.

Here is your roadmap to building a Corbelli-grade analysis pipeline.

### 1. Data Collection: The "Corbelli Standard" for M33 and Others

Corbelli’s method is superior because she does not treat the galaxy as a point mass. She decomposes the velocity into its physical components: **Gas**, **Disk**, and **Bulge**. To match this for M33 and others, you need specific datasets, not just total velocity.

#### **For M33 (The Gold Standard)**

You need the specific radial profiles used in Corbelli et al. (2014) or (2000).

- **Source:** The **CDS (Centre de Données astronomiques de Strasbourg)** via the **VizieR** service.
- **Action:** You can query VizieR programmatically in Python using `astroquery`.
- **Data Needed:**
    - `R` (Radius)
    - `V_obs` (Observed Velocity)
    - `V_gas` (Contribution from HI + H2 gas)
    - `V_stars` (Contribution from stellar disk, derived from photometry)

#### **For Other Galaxies (Generalizing)**

You cannot manually do this for 100+ galaxies. You need a database that has _already_ done the photometry-to-mass conversion.

- **Primary Source:** **SPARC (Spitzer Photometry & Accurate Rotation Curves)**.
    - This is the dataset used by Lelli, McGaugh, & Schombert. It is the modern standard for this work.
    - **Why:** For 175 galaxies, they provide the pre-calculated `V_gas`, `V_disk`, and `V_bulge`. This saves you years of reducing raw telescope data.
- **Secondary Source:** **THINGS (The HI Nearby Galaxy Survey)**. High-resolution data, but you often have to do the mass modeling yourself. Stick to SPARC for Phase II.
    

### 2. Database Design (SQLite)

A SQLite database is perfect here. It’s portable (single file), fast, and works with Unity (C#) and Python.

**Schema Design**

You need three core tables: `Galaxies` (metadata), `RadialProfiles` (the raw curves), and `OmegaResults` (your calculated fits).

SQL

```
-- Table 1: Galaxy Metadata (Static properties)
CREATE TABLE Galaxies (
    galaxy_id TEXT PRIMARY KEY,    -- e.g., 'NGC5055'
    distance_mpc REAL,
    inclination REAL,
    position_angle REAL,
    luminosity_band_36 REAL,       -- 3.6 micron luminosity (stellar mass proxy)
    quality_flag INTEGER           -- 1=High (Corbelli/SPARC standard), 0=Low
);

-- Table 2: The Physical Data (1-to-Many relationship with Galaxies)
CREATE TABLE RadialProfiles (
    profile_id INTEGER PRIMARY KEY AUTOINCREMENT,
    galaxy_id TEXT,
    radius_kpc REAL,
    v_obs REAL,           -- Observed velocity (km/s)
    v_err REAL,           -- Error in observation
    v_gas REAL,           -- Velocity contribution from Gas
    v_disk REAL,          -- Velocity contribution from Disk
    v_bulge REAL,         -- Velocity contribution from Bulge
    v_baryon_total REAL,  -- Sqrt(V_gas^2 + V_disk^2 + V_bulge^2)
    FOREIGN KEY(galaxy_id) REFERENCES Galaxies(galaxy_id)
);

-- Table 3: Your Analysis Results
CREATE TABLE OmegaFits (
    fit_id INTEGER PRIMARY KEY AUTOINCREMENT,
    galaxy_id TEXT,
    omega_value REAL,     -- The best-fit omega (km/s/kpc)
    chi_squared REAL,     -- Goodness of fit
    residuals_rmse REAL,  -- How much velocity is still unexplained?
    method_version TEXT,  -- e.g., 'Flynn_Phase2_Baryonic'
    timestamp DATETIME DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY(galaxy_id) REFERENCES Galaxies(galaxy_id)
);
```

### 3. Updating the Fitting Process (The Math)

This is the most critical scientific update.

**The "Old" Flynn Method (Phase I)**

$$V_{obs} = V_{Kepler} + R \cdot \omega$$

- _Flaw:_ It calculated $\omega$ using only the last data point and assumed the galaxy acted like a point mass ($V_{Kepler} \propto 1/\sqrt{R}$) at the edge.

**The "Corbelli" Method (Phase II)**

You must fit $\omega$ to the **entire curve**, using the **baryonic potential**.

1. **Calculate Baryonic Velocity ($V_{bary}$):**
    $$V_{bary}(R) = \sqrt{V_{gas}(R)^2 + V_{disk}(R)^2 \cdot \Upsilon_{disk} + V_{bulge}(R)^2 \cdot \Upsilon_{bulge}}$$
    - Note: $\Upsilon$ (Upsilon) is the Mass-to-Light ratio. SPARC provides reasonable defaults, but you can leave them as free parameters or fix them like Corbelli did.
2. **The New Fitting Equation:**
    Instead of a dark matter halo, you represent the "missing physics" as your $\omega$ term.
    $$V_{model}(R) = V_{bary}(R) + R \cdot \omega$$
3. **The Algorithm:**
    Instead of solving for $\omega$ at the last point, you use **Non-Linear Least Squares** (e.g., `scipy.optimize.curve_fit`) to find the single $\omega$ value that minimizes the error across all radii $R$.
    This makes your $\omega$ statistically robust and less sensitive to outliers at the galaxy's edge.

### 4. Implementation: Python Script

Here is a Python script to fetch M33 data from SPARC (or a local file), store it in your SQLite DB, and perform the new "Corbelli-style" fit.

**Prerequisites:**

`pip install pandas numpy scipy sqlalchemy astroquery`

Python

```
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sqlalchemy import create_engine, text
import matplotlib.pyplot as plt

# 1. SETUP DATABASE
engine = create_engine('sqlite:///galaxy_dynamics.db')

def init_db():
    with engine.connect() as conn:
        conn.execute(text("""
            CREATE TABLE IF NOT EXISTS RadialProfiles (
                galaxy_id TEXT, radius_kpc REAL, v_obs REAL, v_err REAL, 
                v_gas REAL, v_disk REAL, v_bulge REAL
            )
        """))
        conn.execute(text("""
            CREATE TABLE IF NOT EXISTS OmegaFits (
                galaxy_id TEXT, omega_value REAL, chi_squared REAL
            )
        """))

# 2. LOAD DATA (Example using SPARC format for M33 or similar)
# Note: In reality, download the SPARC .dat files and parse them.
def load_galaxy_data(filepath, galaxy_name):
    # SPARC files usually have columns: Rad, Vobs, Vgas, Vdisk, Vbulge
    df = pd.read_csv(filepath, sep='\s+', names=['Rad', 'Vobs', 'Verr', 'Vgas', 'Vdisk', 'Vbulge', 'SBdisk', 'SBbulge'], skiprows=1)
    
    # Save to SQLite
    df['galaxy_id'] = galaxy_name
    df[['galaxy_id', 'Rad', 'Vobs', 'Verr', 'Vgas', 'Vdisk', 'Vbulge']].to_sql(
        'RadialProfiles', engine, if_exists='append', index=False
    )
    return df

# 3. THE NEW "CORBELLI-STYLE" OMEGA FIT
def fit_omega_baryonic(df):
    """
    Fits V_obs = V_bary + R * omega
    where V_bary is the vector sum of gas, disk, and bulge components.
    """
    # Calculate Baryonic Velocity (No Dark Matter)
    # Note: absolute values needed because V^2 contributions can be negative in some raw data conventions if V < 0
    v_gas_sq = np.abs(df['Vgas']) * df['Vgas']
    v_disk_sq = np.abs(df['Vdisk']) * df['Vdisk']
    v_bulge_sq = np.abs(df['Vbulge']) * df['Vbulge']
    
    # Total Baryonic Velocity (V_bary)
    # We assume Mass-to-Light ratio (Upsilon) is 1.0 for simplicity here, 
    # but Corbelli treats this as a variable.
    v_bary = np.sqrt(np.abs(v_gas_sq + v_disk_sq + v_bulge_sq))
    
    # Define the model function
    def omega_model(r, omega):
        # The hypothesis: Observed = Baryonic + Omega * R
        return v_bary + (r * omega)
    
    # Perform the fit
    # We weight by 1/Verr to trust precise measurements more
    popt, pcov = curve_fit(
        omega_model, 
        df['Rad'], 
        df['Vobs'], 
        sigma=df['Verr'],
        absolute_sigma=True
    )
    
    omega_best = popt[0]
    perr = np.sqrt(np.diag(pcov))[0]
    
    # Calculate Chi-Squared
    residuals = df['Vobs'] - omega_model(df['Rad'], omega_best)
    chi2 = np.sum((residuals / df['Verr'])**2)
    
    return omega_best, perr, chi2, v_bary

# --- USAGE EXAMPLE ---
# init_db()
# # You would point this to a real SPARC file
# # df = load_galaxy_data('NGC5055_rotmod.dat', 'NGC5055') 
# # omega, err, chi2, v_bary = fit_omega_baryonic(df)
# # print(f"New Corbelli-Style Omega: {omega:.4f} +/- {err:.4f}")
```

### Next Steps

1. **Download the SPARC dataset:** It contains `_rotmod.dat` files for 175 galaxies which exactly fit the format needed for the script above.
2. **Run M33 First:** Verify you can replicate a curve that looks like Corbelli’s. Note that Corbelli finds $\omega$ isn't needed if you add a Dark Matter halo; your goal is to show that $\omega$ _replaces_ the need for that halo.
3. **Batch Process:** Once M33 works, loop through all 175 SPARC galaxies and populate your `OmegaFits` table. This will give you a massive, high-quality dataset for your Phase II analysis.