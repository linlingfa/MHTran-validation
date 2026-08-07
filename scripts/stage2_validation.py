"""
Stage 2: Hazard-Attributed Outage Validation

Evaluates hazards with sufficient and suitable NOAA Storm Events
coverage by comparing county-level MHTran exposure with
hazard-attributed EAGLE-I outage frequency.

Hazards retained for Stage 2:
- Flood
- Hail
- Landslide
- Lightning
- Tornado

Validation is performed using:
1. Spearman rank correlation
2. Top-5% hotspot overlap

Both maximum and mean county-level exposure are evaluated.

An EAGLE-I outage observation is attributed to a NOAA event when:
1. Both occur in the same county.
2. The outage timestamp occurs between the NOAA event start time
   and the event end time plus the hazard-specific attribution window.

Duplicate EAGLE-I observations matched to overlapping NOAA events
are removed before county-level outage frequency is calculated.
"""

import pandas as pd
import geopandas as gpd
from pathlib import Path
from scipy.stats import spearmanr


# ============================================================
# 1. Configuration
# ============================================================

DATA_DIR = Path("Data")
MODEL_DIR = DATA_DIR / "model_output" / "data" / "hotspots"
EAGLE_DIR = DATA_DIR / "eaglei_outages"
NOAA_DIR = DATA_DIR / "NOAA storms data"

COUNTY_URL = (
    "https://www2.census.gov/geo/tiger/GENZ2022/shp/"
    "cb_2022_us_county_500k.zip"
)

# Select hazard to validate
hazard_name = "flood"

HAZARDS = {
    "flood": {
        "file": "subs_flood.parquet",
        "metric": "flood_depth_m",
        "noaa_events": ["Flood", "Flash Flood", "Heavy Rain"],
        "time_window_hours": 24,
    },
    "hail": {
        "file": "subs_hail.parquet",
        "metric": "hail_rate",
        "noaa_events": ["Hail"],
        "time_window_hours": 24,
    },
    "landslide": {
        "file": "subs_landslide.parquet",
        "metric": "ls_susc",
        "noaa_events": ["Debris Flow"],
        "time_window_hours": 48,
    },
    "lightning": {
        "file": "subs_lightning.parquet",
        "metric": "ltng_rate",
        "noaa_events": ["Lightning"],
        "time_window_hours": 24,
    },
    "tornado": {
        "file": "subs_tornado.parquet",
        "metric": "annual_hit_rate",
        "noaa_events": ["Tornado"],
        "time_window_hours": 24,
    },
}

hazard_name = hazard_name.lower()

if hazard_name not in HAZARDS:
    raise ValueError(
        f"Invalid hazard '{hazard_name}'. "
        f"Choose from: {list(HAZARDS.keys())}"
    )

config = HAZARDS[hazard_name]

hazard_file = config["file"]
hazard_metric = config["metric"]
hazards_list = config["noaa_events"]
time_window_hours = config["time_window_hours"]

if hazards_list is None:
    raise ValueError(
        f"{hazard_name.title()} has no corresponding NOAA Storm Events "
        "category for Stage 2."
    )


# ============================================================
# 2. Load and aggregate MHTran exposure
# ============================================================

hazard_path = MODEL_DIR / hazard_file

if not hazard_path.exists():
    raise FileNotFoundError(
        f"Hazard file not found: {hazard_path}"
    )

subs = pd.read_parquet(hazard_path)

counties = gpd.read_file(
    COUNTY_URL
)[["GEOID", "geometry"]]

counties["GEOID"] = (
    counties["GEOID"]
    .astype(str)
    .str.zfill(5)
)

# Convert substations to geographic points
gdf = gpd.GeoDataFrame(
    subs,
    geometry=gpd.points_from_xy(
        subs["lon"],
        subs["lat"]
    ),
    crs="EPSG:4326"
)

counties = counties.to_crs(gdf.crs)

# Remove any pre-existing GEOID so all county IDs
# are assigned consistently from coordinates
gdf = gdf.drop(
    columns=["GEOID"],
    errors="ignore"
)

# Assign each substation to a county
gdf = gpd.sjoin(
    gdf,
    counties,
    how="inner",
    predicate="intersects"
).drop(
    columns=["index_right"],
    errors="ignore"
)

# Convert exposure metric to numeric
gdf[hazard_metric] = pd.to_numeric(
    gdf[hazard_metric],
    errors="coerce"
)

# Remove substations without valid exposure
gdf = gdf.dropna(
    subset=[hazard_metric]
)

# Calculate maximum and mean county exposure
county_exposure = (
    gdf.groupby("GEOID")[hazard_metric]
    .agg(["max", "mean"])
    .reset_index()
    .rename(
        columns={
            "max": "max_exposure",
            "mean": "mean_exposure",
        }
    )
)

# ============================================================
# 3. Load EAGLE-I outage observations
# ============================================================

eagle_files = sorted(
    EAGLE_DIR.glob("eaglei_outages_*.csv")
)

if not eagle_files:
    raise FileNotFoundError(
        f"No EAGLE-I files found in {EAGLE_DIR}"
    )

eagle_list = []

for file in eagle_files:

    year = int(file.stem.split("_")[-1])

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
            "run_start_time",
            outage_col,
        ],
        dtype={"fips_code": "string"},
    )

    df[outage_col] = pd.to_numeric(
        df[outage_col],
        errors="coerce"
    ).fillna(0)

    # Keep only observations with positive outages
    df = df[
        df[outage_col] > 0
    ].copy()

    df = df.rename(
        columns={
            outage_col: "outage_value"
        }
    )

    eagle_list.append(df)

eagle = pd.concat(
    eagle_list,
    ignore_index=True
)

eagle["run_start_time"] = pd.to_datetime(
    eagle["run_start_time"],
    errors="coerce"
)

eagle["GEOID"] = (
    eagle["fips_code"]
    .astype("string")
    .str.zfill(5)
)

eagle = eagle.dropna(
    subset=[
        "run_start_time",
        "GEOID",
    ]
)

eagle = eagle.sort_values(
    ["GEOID", "run_start_time"]
).reset_index(drop=True)


# ============================================================
# 4. Load and prepare NOAA Storm Events
# ============================================================

noaa_files = sorted(
    NOAA_DIR.glob("StormEvents_details*.csv")
)

if not noaa_files:
    raise FileNotFoundError(
        f"No NOAA Storm Events files found in {NOAA_DIR}"
    )

noaa_list = []

for file in noaa_files:

    df = pd.read_csv(
        file,
        usecols=[
            "BEGIN_DATE_TIME",
            "END_DATE_TIME",
            "STATE",
            "CZ_TYPE",
            "CZ_FIPS",
            "EVENT_TYPE",
        ],
    )

    noaa_list.append(df)

noaa = pd.concat(
    noaa_list,
    ignore_index=True
)

noaa["BEGIN_DATE_TIME"] = pd.to_datetime(
    noaa["BEGIN_DATE_TIME"],
    format="%d-%b-%y %H:%M:%S",
    errors="coerce"
)

noaa["END_DATE_TIME"] = pd.to_datetime(
    noaa["END_DATE_TIME"],
    format="%d-%b-%y %H:%M:%S",
    errors="coerce"
)

state_lookup = {
    "ALABAMA": "01",
    "ALASKA": "02",
    "ARIZONA": "04",
    "ARKANSAS": "05",
    "CALIFORNIA": "06",
    "COLORADO": "08",
    "CONNECTICUT": "09",
    "DELAWARE": "10",
    "DISTRICT OF COLUMBIA": "11",
    "FLORIDA": "12",
    "GEORGIA": "13",
    "HAWAII": "15",
    "IDAHO": "16",
    "ILLINOIS": "17",
    "INDIANA": "18",
    "IOWA": "19",
    "KANSAS": "20",
    "KENTUCKY": "21",
    "LOUISIANA": "22",
    "MAINE": "23",
    "MARYLAND": "24",
    "MASSACHUSETTS": "25",
    "MICHIGAN": "26",
    "MINNESOTA": "27",
    "MISSISSIPPI": "28",
    "MISSOURI": "29",
    "MONTANA": "30",
    "NEBRASKA": "31",
    "NEVADA": "32",
    "NEW HAMPSHIRE": "33",
    "NEW JERSEY": "34",
    "NEW MEXICO": "35",
    "NEW YORK": "36",
    "NORTH CAROLINA": "37",
    "NORTH DAKOTA": "38",
    "OHIO": "39",
    "OKLAHOMA": "40",
    "OREGON": "41",
    "PENNSYLVANIA": "42",
    "RHODE ISLAND": "44",
    "SOUTH CAROLINA": "45",
    "SOUTH DAKOTA": "46",
    "TENNESSEE": "47",
    "TEXAS": "48",
    "UTAH": "49",
    "VERMONT": "50",
    "VIRGINIA": "51",
    "WASHINGTON": "53",
    "WEST VIRGINIA": "54",
    "WISCONSIN": "55",
    "WYOMING": "56",
}

# Convert state names to state FIPS codes
noaa["STATEFP"] = noaa["STATE"].map(
    state_lookup
)

noaa = noaa.dropna(
    subset=[
        "STATEFP",
        "BEGIN_DATE_TIME",
        "END_DATE_TIME",
    ]
)

# Retain only county-based NOAA events
noaa = noaa[
    noaa["CZ_TYPE"] == "C"
].copy()

# Format county FIPS
noaa["CZ_FIPS"] = (
    pd.to_numeric(
        noaa["CZ_FIPS"],
        errors="coerce"
    )
    .astype("Int64")
    .astype("string")
    .str.zfill(3)
)

noaa = noaa.dropna(
    subset=["CZ_FIPS"]
)

# Construct five-digit county GEOID
noaa["GEOID"] = (
    noaa["STATEFP"]
    + noaa["CZ_FIPS"]
)

# Retain NOAA categories corresponding to selected hazard
noaa_hazard = noaa[
    noaa["EVENT_TYPE"].isin(
        hazards_list
    )
].copy()

noaa_hazard = noaa_hazard.rename(
    columns={
        "BEGIN_DATE_TIME": "event_start",
        "END_DATE_TIME": "event_end",
    }
)

noaa_hazard = noaa_hazard.sort_values(
    ["GEOID", "event_start"]
).reset_index(drop=True)

print("\nNOAA EVENT COVERAGE")
print(noaa_hazard["EVENT_TYPE"].value_counts())
print(f"Total qualifying NOAA events: {len(noaa_hazard)}")


# ============================================================
# 5. Attribute EAGLE-I outages to NOAA events
# ============================================================

TIME_WINDOW = pd.Timedelta(
    hours=time_window_hours
)

hazard_matches = []

for county, outages in eagle.groupby("GEOID"):

    hazards = noaa_hazard[
        noaa_hazard["GEOID"] == county
    ]

    for _, hazard in hazards.iterrows():

        matches = outages[
            (
                outages["run_start_time"]
                >= hazard["event_start"]
            )
            &
            (
                outages["run_start_time"]
                <= hazard["event_end"] + TIME_WINDOW
            )
        ]

        if not matches.empty:
            hazard_matches.append(matches)

hazard_outages = (
    pd.concat(
        hazard_matches,
        ignore_index=True
    )
    if hazard_matches
    else eagle.iloc[0:0].copy()
)

# Prevent one EAGLE-I observation from being counted multiple
# times when it matches overlapping NOAA events
hazard_outages = hazard_outages.drop_duplicates(
    subset=[
        "GEOID",
        "run_start_time",
    ]
)

# Count hazard-attributed outage observations per county
hazard_outage_count = (
    hazard_outages
    .groupby("GEOID")
    .size()
    .reset_index(
        name="hazard_outage_frequency"
    )
)

# ============================================================
# 6. Merge modeled exposure with attributed outages
# ============================================================

validation = county_exposure.merge(
    hazard_outage_count,
    on="GEOID",
    how="left"
)

# Modeled counties with no attributed outages receive zero
validation["hazard_outage_frequency"] = (
    validation["hazard_outage_frequency"]
    .fillna(0)
)

# ============================================================
# 7. Spearman rank correlation
# ============================================================

rho_max, p_max = spearmanr(
    validation["max_exposure"],
    validation["hazard_outage_frequency"],
    nan_policy="omit"
)

rho_mean, p_mean = spearmanr(
    validation["mean_exposure"],
    validation["hazard_outage_frequency"],
    nan_policy="omit"
)

valid_max = validation[
    [
        "max_exposure",
        "hazard_outage_frequency",
    ]
].dropna()

valid_mean = validation[
    [
        "mean_exposure",
        "hazard_outage_frequency",
    ]
].dropna()

print("\nSPEARMAN CORRELATION RESULTS")

print(
    f"Maximum exposure: "
    f"rho = {rho_max:.3f}, "
    f"p = {p_max:.3e}, "
    f"N = {len(valid_max)}"
)

print(
    f"Mean exposure:    "
    f"rho = {rho_mean:.3f}, "
    f"p = {p_mean:.3e}, "
    f"N = {len(valid_mean)}"
)


# ============================================================
# 8. Top-5% hotspot overlap
# ============================================================

n_hot = max(
    1,
    round(len(validation) * 0.05)
)

model_hot_max = set(
    validation.sort_values(
        ["max_exposure", "GEOID"],
        ascending=[False, True]
    )
    .head(n_hot)["GEOID"]
)

model_hot_mean = set(
    validation.sort_values(
        ["mean_exposure", "GEOID"],
        ascending=[False, True]
    )
    .head(n_hot)["GEOID"]
)

observed_hot = set(
    validation.sort_values(
        ["hazard_outage_frequency", "GEOID"],
        ascending=[False, True]
    )
    .head(n_hot)["GEOID"]
)

overlap_max = (
    model_hot_max
    & observed_hot
)

overlap_mean = (
    model_hot_mean
    & observed_hot
)

print("\nHOTSPOT OVERLAP RESULTS")
print(
    "Counties in each hotspot set:",
    n_hot
)

print(
    f"Maximum exposure overlap: "
    f"{len(overlap_max)} counties "
    f"({100 * len(overlap_max) / n_hot:.1f}%)"
)

print(
    f"Mean exposure overlap: "
    f"{len(overlap_mean)} counties "
    f"({100 * len(overlap_mean) / n_hot:.1f}%)"
)