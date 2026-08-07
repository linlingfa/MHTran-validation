"""
Stage 3: Integrated Multi-Hazard Validation

Constructs a county-level composite MHTran multi-hazard exposure
score and compares it with overall EAGLE-I outage frequency.

Validation is performed using:
1. Spearman rank correlation
2. Top-5% hotspot overlap
3. Geographically Weighted Regression (GWR)

For each hazard:
- Substation exposure is aggregated to the county level using
  the maximum valid exposure value.
- County exposure is converted to a percentile rank from 0 to 1.

The composite multi-hazard score is calculated as the mean
percentile rank across hazards available for each county.
Missing hazard values remain unavailable rather than being
assigned zero, and hazard_count records the number of hazards
contributing to each county's score.

GWR evaluates spatial variation in the relationship between
the composite exposure score and log-transformed outage frequency.
"""

import numpy as np
import pandas as pd
import geopandas as gpd

from pathlib import Path
from functools import reduce
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from mgwr.sel_bw import Sel_BW
from mgwr.gwr import GWR
from joblib import parallel_config


# ============================================================
# 1. Configuration
# ============================================================

DATA_DIR = Path("Data")
MODEL_DIR = DATA_DIR / "model_output" / "data" / "hotspots"
EAGLE_DIR = DATA_DIR / "eaglei_outages"

COUNTY_URL = (
    "https://www2.census.gov/geo/tiger/GENZ2022/shp/"
    "cb_2022_us_county_500k.zip"
)

# Define the nine MHTran hazards and their exposure metrics
HAZARDS = {
    "flood": "flood_depth_m",
    "fzg": "R_50RP",
    "hail": "hail_rate",
    "landslide": "ls_susc",
    "lightning": "ltng_rate",
    "seismic": "pga_475",
    "tornado": "annual_hit_rate",
    "wildfire": "whp_cnt",
    "wind": "pe_5pct_maxwind",
}


# ============================================================
# 2. Load U.S. county boundaries
# ============================================================

counties = gpd.read_file(
    COUNTY_URL
)[["GEOID", "geometry"]]

# Standardize county GEOIDs to five digits
counties["GEOID"] = (
    counties["GEOID"]
    .astype(str)
    .str.zfill(5)
)


# ============================================================
# 3. Create county-level exposure for each hazard
# ============================================================

def create_county_exposure(
    hazard_name,
    exposure_metric,
):
    """
    Load one MHTran hazard layer, assign substations to counties,
    and calculate the maximum valid exposure value per county.
    """

    hazard_path = (
        MODEL_DIR
        / f"subs_{hazard_name}.parquet"
    )

    if not hazard_path.exists():
        raise FileNotFoundError(
            f"Hazard file not found: {hazard_path}"
        )

    subs = pd.read_parquet(
        hazard_path
    )

    # Convert substation coordinates to geographic points
    gdf = gpd.GeoDataFrame(
        subs,
        geometry=gpd.points_from_xy(
            subs["lon"],
            subs["lat"],
        ),
        crs="EPSG:4326",
    )

    # Match the county CRS
    gdf = gdf.to_crs(
        counties.crs
    )

    # Remove any pre-existing GEOID so county IDs are
    # assigned consistently from substation coordinates
    gdf = gdf.drop(
        columns=["GEOID"],
        errors="ignore",
    )

    # Spatially assign each substation to a county
    gdf = gpd.sjoin(
        gdf,
        counties,
        how="inner",
        predicate="intersects",
    ).drop(
        columns=["index_right"],
        errors="ignore",
    )

    gdf["GEOID"] = (
        gdf["GEOID"]
        .astype(str)
        .str.zfill(5)
    )

    # Convert the selected exposure metric to numeric
    gdf[exposure_metric] = pd.to_numeric(
        gdf[exposure_metric],
        errors="coerce",
    )

    # Remove substations without valid exposure
    gdf = gdf.dropna(
        subset=[exposure_metric]
    )

    # Use maximum valid substation exposure as county exposure
    county_exposure = (
        gdf
        .groupby("GEOID")[exposure_metric]
        .max()
        .reset_index()
        .rename(
            columns={
                exposure_metric:
                f"{hazard_name}_exposure"
            }
        )
    )

    return county_exposure


# Generate county-level exposure table for each hazard
exposure_tables = []

for hazard, metric in HAZARDS.items():

    exposure = create_county_exposure(
        hazard,
        metric,
    )

    exposure_tables.append(
        exposure
    )


# ============================================================
# 4. Construct composite multi-hazard exposure score
# ============================================================

# Outer-join all hazard layers so counties are retained even
# when one or more hazard exposure values are unavailable
multi_hazard = reduce(
    lambda left, right: pd.merge(
        left,
        right,
        on="GEOID",
        how="outer",
    ),
    exposure_tables,
)

hazard_columns = [
    column
    for column in multi_hazard.columns
    if column.endswith("_exposure")
]

# Convert each hazard exposure layer to a percentile rank
# from 0 to 1 relative to counties with data for that hazard
for column in hazard_columns:

    multi_hazard[
        f"{column}_norm"
    ] = (
        multi_hazard[column]
        .rank(pct=True)
    )

norm_columns = [
    column
    for column in multi_hazard.columns
    if column.endswith("_norm")
]

# Calculate the composite score as the mean percentile rank
# across hazards available for each county.
# Missing hazard values remain NaN and are skipped by mean().
multi_hazard["multi_hazard_score"] = (
    multi_hazard[norm_columns]
    .mean(axis=1)
)

# Record the number of hazard layers contributing to each score
multi_hazard["hazard_count"] = (
    multi_hazard[norm_columns]
    .notna()
    .sum(axis=1)
)


# ============================================================
# 5. Load and calculate EAGLE-I outage frequency
# ============================================================

eagle_files = sorted(
    EAGLE_DIR.glob(
        "eaglei_outages_*.csv"
    )
)

if not eagle_files:
    raise FileNotFoundError(
        f"No EAGLE-I files found in {EAGLE_DIR}"
    )

eagle_list = []

for file in eagle_files:

    year = int(
        file.stem.split("_")[-1]
    )

    # EAGLE-I changed outage-value column beginning in 2024
    outage_col = (
        "sum"
        if year <= 2023
        else "customers_out"
    )

    df = pd.read_csv(
        file,
        usecols=[
            "fips_code",
            outage_col,
        ],
        dtype={
            "fips_code": "string"
        },
    )

    df[outage_col] = pd.to_numeric(
        df[outage_col],
        errors="coerce",
    ).fillna(0)

    # Retain only observations with positive recorded outages
    df = df[
        df[outage_col] > 0
    ].copy()

    eagle_list.append(df)

eagle = pd.concat(
    eagle_list,
    ignore_index=True,
)

eagle["GEOID"] = (
    eagle["fips_code"]
    .astype("string")
    .str.zfill(5)
)

# Count positive EAGLE-I observations in each county
outage_frequency = (
    eagle
    .groupby("GEOID")
    .size()
    .reset_index(
        name="outage_frequency"
    )
)

outage_frequency["GEOID"] = (
    outage_frequency["GEOID"]
    .astype(str)
    .str.zfill(5)
)


# ============================================================
# 6. Merge composite exposure with outage frequency
# ============================================================

final_validation = multi_hazard.merge(
    outage_frequency,
    on="GEOID",
    how="left",
)

# Counties with modeled exposure but no positive EAGLE-I
# observations are assigned an outage frequency of zero
final_validation["outage_frequency"] = (
    final_validation[
        "outage_frequency"
    ]
    .fillna(0)
)

# ============================================================
# 7. National Spearman rank correlation
# ============================================================

rho, p_value = spearmanr(
    final_validation[
        "multi_hazard_score"
    ],
    final_validation[
        "outage_frequency"
    ],
    nan_policy="omit",
)

valid_spearman = final_validation[
    [
        "multi_hazard_score",
        "outage_frequency",
    ]
].dropna()

print("\nSPEARMAN CORRELATION RESULTS")
print(f"Spearman rho = {rho:.3f}")
print(f"p-value = {p_value:.3e}")
print(f"N = {len(valid_spearman)}")


# ============================================================
# 8. Top-5% hotspot overlap
# ============================================================

# Define the top 5% of validation counties as hotspots
top_5_percent = int(
    len(final_validation) * 0.05
)

# Modeled hotspots: counties with highest composite exposure
model_hotspots = set(
    final_validation
    .nlargest(
        top_5_percent,
        "multi_hazard_score",
    )["GEOID"]
)

# Observed hotspots: counties with highest outage frequency
observed_hotspots = set(
    final_validation
    .nlargest(
        top_5_percent,
        "outage_frequency",
    )["GEOID"]
)

overlap = (
    model_hotspots
    & observed_hotspots
)

print("\nHOTSPOT OVERLAP RESULTS")
print(
    "Counties in each hotspot set:",
    top_5_percent,
)
print(
    "Overlap counties:",
    len(overlap),
)
print(
    f"Percentage hotspot overlap: "
    f"{100 * len(overlap) / len(model_hotspots):.1f}%"
)


# ============================================================
# 9. Prepare variables for GWR
# ============================================================

# Attach county geometry to the validation dataset
gwr_data = counties.merge(
    final_validation,
    on="GEOID",
    how="inner",
).copy()

# Retain counties with valid composite exposure and outage data
gwr_data = gwr_data.dropna(
    subset=[
        "multi_hazard_score",
        "outage_frequency",
    ]
)

# Log-transform outage frequency to reduce right-skewness
# and the influence of counties with extremely high counts
gwr_data["log_outage"] = np.log1p(
    gwr_data["outage_frequency"]
)

# Project counties to a projected CRS suitable for U.S.
# distance-based spatial analysis
gwr_data = gwr_data.to_crs(
    "EPSG:5070"
)

# Represent each county by its polygon centroid
centroids = (
    gwr_data.geometry.centroid
)

coords = np.column_stack(
    [
        centroids.x,
        centroids.y,
    ]
)

# Dependent variable: log-transformed outage frequency
y = gwr_data[
    ["log_outage"]
].values

# Independent variable: standardized composite exposure score
X = StandardScaler().fit_transform(
    gwr_data[
        ["multi_hazard_score"]
    ]
)


# ============================================================
# 10. Fit Geographically Weighted Regression
# ============================================================

# Select an adaptive bandwidth using AICc
selector = Sel_BW(
    coords,
    y,
    X,
)

# Restrict parallel processing for reproducibility and
# compatibility across execution environments
with parallel_config(n_jobs=1):

    bw = selector.search(
        criterion="AICc"
    )

print(
    "\nSelected GWR bandwidth:",
    bw,
)

# Fit GWR using the selected adaptive bandwidth.
# mgwr automatically includes the intercept.
gwr_model = GWR(
    coords,
    y,
    X,
    bw,
)

with parallel_config(n_jobs=1):

    gwr_results = (
        gwr_model.fit()
    )


# ============================================================
# 11. Extract spatial GWR results
# ============================================================

# Local coefficient for standardized multi-hazard exposure
# Column 0 is the intercept; column 1 is exposure
gwr_data["local_coefficient"] = (
    gwr_results.params[:, 1]
)

# Local model explanatory power
gwr_data["local_R2"] = (
    gwr_results.localR2
)

# Standard error of each local exposure coefficient
gwr_data["coefficient_SE"] = (
    gwr_results.bse[:, 1]
)

print("\nGWR RESULTS")
print(
    gwr_results.summary()
)

print(
    "Local coefficient range:",
    f"{gwr_data['local_coefficient'].min():.3f}",
    "to",
    f"{gwr_data['local_coefficient'].max():.3f}",
)