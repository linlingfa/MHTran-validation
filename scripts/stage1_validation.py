"""
Stage 1: Individual Hazard Exposure Validation

Compares county-level MHTran hazard exposure with EAGLE-I outage
frequency and severity using:
1. Spearman rank correlation
2. Top-5% hotspot overlap

County exposure is represented by the maximum valid substation
exposure value within each county.
"""

import pandas as pd
import geopandas as gpd
from pathlib import Path
from scipy.stats import spearmanr

DATA_DIR = Path("Data")
MODEL_DIR = DATA_DIR / "model_output" / "data" / "hotspots"
EAGLE_DIR = DATA_DIR / "eaglei_outages"

COUNTY_URL = (
    "https://www2.census.gov/geo/tiger/GENZ2022/shp/"
    "cb_2022_us_county_500k.zip"
)

# ============================================================
# 1. Configuration
# ============================================================

# Select hazard to validate
hazard_name = "flood"

HAZARDS = {
    "flood":      {"file": "subs_flood.parquet",      "metric": "flood_depth_m"},
    "fzg":        {"file": "subs_fzg.parquet",        "metric": "R_50RP"},
    "hail":       {"file": "subs_hail.parquet",       "metric": "hail_rate"},
    "landslide":  {"file": "subs_landslide.parquet",  "metric": "ls_susc"},
    "lightning":  {"file": "subs_lightning.parquet",  "metric": "ltng_rate"},
    "seismic":    {"file": "subs_seismic.parquet",    "metric": "pga_475"},
    "tornado":    {"file": "subs_tornado.parquet",    "metric": "annual_hit_rate"},
    "wildfire":   {"file": "subs_wildfire.parquet",   "metric": "whp_cnt"},
    "wind":       {"file": "subs_wind.parquet",       "metric": "pe_5pct_maxwind"},
}

if hazard_name not in HAZARDS:
    raise ValueError(
        f"Invalid hazard '{hazard_name}'. "
        f"Choose from: {list(HAZARDS.keys())}"
    )

hazard_file = HAZARDS[hazard_name]["file"]
hazard_metric = HAZARDS[hazard_name]["metric"]

# ============================================================
# 2. Load and aggregate MHTran exposure
# ============================================================

#reads the parquet file for that specific subs hazard
hazard_path = MODEL_DIR / hazard_file

if not hazard_path.exists():
    raise FileNotFoundError(f"Hazard file not found: {hazard_path}")

subs = pd.read_parquet(hazard_path)

# Load county boundaries
counties = gpd.read_file(COUNTY_URL)[["GEOID", "geometry"]]

counties["GEOID"] = counties["GEOID"].astype(str).str.zfill(5)

# Convert substations to GeoDataFrame using coordinates
gdf = gpd.GeoDataFrame(
    subs,
    geometry=gpd.points_from_xy(
        subs["lon"],
        subs["lat"]
    ),
    crs="EPSG:4326"
)

# Match coordinate CRS
counties = counties.to_crs(gdf.crs)

gdf = gdf.drop(columns=["GEOID"], errors="ignore")

# Assign each substation a county GEOID
gdf = gpd.sjoin(
    gdf,
    counties,
    how="inner",
    predicate="intersects"
).drop(columns=["index_right"], errors="ignore")

gdf[hazard_metric] = pd.to_numeric(
    gdf[hazard_metric],
    errors="coerce"
)

gdf = gdf.dropna(subset=[hazard_metric])

# Aggregate hazard exposure by county
county_exposure = (
    gdf.groupby("GEOID")[hazard_metric]
    .max()
    .reset_index()
    .rename(columns={hazard_metric: "hazard_exposure"})
)

# County-level modeled exposure used for validation
validation_model = county_exposure

files = sorted(EAGLE_DIR.glob("eaglei_outages_*.csv"))

if not files:
    raise FileNotFoundError(
        f"No EAGLE-I files found in {EAGLE_DIR}"
    )

# ============================================================
# 3. Calculate EAGLE-I outage frequency
# ============================================================

eagle_list_freq = []

for file in files:
    year = int(file.stem.split("_")[-1])
    col = "sum" if year <= 2023 else "customers_out"

    df = pd.read_csv(
        file,
        usecols=["fips_code", col],
        dtype={"fips_code": "string"}
    )

    df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[df[col] > 0]

    eagle_list_freq.append(df)

eagle_freq = pd.concat(eagle_list_freq, ignore_index=True)
eagle_freq["GEOID"] = eagle_freq["fips_code"].str.zfill(5)

outage_county_freq = (
    eagle_freq.groupby("GEOID")
    .size()
    .reset_index(name="outage_frequency")
)

final_freq = validation_model.merge(
    outage_county_freq,
    on="GEOID",
    how="left"
)

final_freq["outage_frequency"] = final_freq["outage_frequency"].fillna(0)

# ============================================================
# 4. Calculate EAGLE-I outage severity
# ============================================================

eagle_list_sev = []

for file in files:
    year = int(file.stem.split("_")[-1])
    col = "sum" if year <= 2023 else "customers_out"

    df = pd.read_csv(
        file,
        usecols=["fips_code", col],
        dtype={"fips_code": "string"}
    ).rename(columns={col: "outage_severity"})

    eagle_list_sev.append(df)

eagle_sev = pd.concat(eagle_list_sev, ignore_index=True)

eagle_sev["GEOID"] = eagle_sev["fips_code"].str.zfill(5)
eagle_sev["outage_severity"] = pd.to_numeric(
    eagle_sev["outage_severity"],
    errors="coerce"
).fillna(0)

outage_county_sev = (
    eagle_sev.groupby("GEOID", as_index=False)["outage_severity"]
    .sum()
)

final_sev = validation_model.merge(
    outage_county_sev,
    on="GEOID",
    how="left"
)

final_sev["outage_severity"] = final_sev["outage_severity"].fillna(0)

# ============================================================
# 5. Spearman rank correlation
# ============================================================

rho_freq, p_freq = spearmanr(
    final_freq["hazard_exposure"],
    final_freq["outage_frequency"]
)

rho_sev, p_sev = spearmanr(
    final_sev["hazard_exposure"],
    final_sev["outage_severity"]
)

print("\nSPEARMAN CORRELATION RESULTS")
print(
    f"Frequency: rho = {rho_freq:.3f}, "
    f"p = {p_freq:.3e}"
)
print(
    f"Severity:  rho = {rho_sev:.3f}, "
    f"p = {p_sev:.3e}"
)

# ============================================================
# 6. Top-5% hotspot overlap
# ============================================================

n_hot_freq = max(1, round(len(final_freq) * 0.05))
n_hot_sev = max(1, round(len(final_sev) * 0.05))

model_hot_freq = set(
    final_freq.sort_values(
        ["hazard_exposure", "GEOID"],
        ascending=[False, True]
    )
    .head(n_hot_freq)["GEOID"]
)

model_hot_sev = set(
    final_sev.sort_values(
        ["hazard_exposure", "GEOID"],
        ascending=[False, True]
    )
    .head(n_hot_sev)["GEOID"]
)

real_hot_freq = set(
    final_freq.sort_values(
        ["outage_frequency", "GEOID"],
        ascending=[False, True]
    )
    .head(n_hot_freq)["GEOID"]
)

real_hot_sev = set(
    final_sev.sort_values(
        ["outage_severity", "GEOID"],
        ascending=[False, True]
    )
    .head(n_hot_sev)["GEOID"]
)

overlap_freq = model_hot_freq & real_hot_freq
overlap_sev = model_hot_sev & real_hot_sev

print()

print("\nFREQUENCY HOTSPOT RESULTS")
print("Counties in each hotspot set:", n_hot_freq)
print("Overlap counties:", len(overlap_freq))
print(
    f"Percentage hotspot overlap: "
    f"{100 * len(overlap_freq) / n_hot_freq:.1f}%"
)

print("\nSEVERITY HOTSPOT RESULTS")
print("Counties in each hotspot set:", n_hot_sev)
print("Overlap counties:", len(overlap_sev))
print(
    f"Percentage hotspot overlap: "
    f"{100 * len(overlap_sev) / n_hot_sev:.1f}%"
)