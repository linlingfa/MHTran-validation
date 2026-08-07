# Data

The raw datasets used in this project are not included in the repository.

This project uses MHTran hazard exposure data, EAGLE-I outage observations, NOAA Storm Events, and U.S. Census Bureau county boundaries.

## Expected Local Data Structure

The validation scripts expect the following local directory structure:

```text
Data/
├── model_output/
│   └── data/
│       └── hotspots/
│           ├── subs_flood.parquet
│           ├── subs_fzg.parquet
│           ├── subs_hail.parquet
│           ├── subs_landslide.parquet
│           ├── subs_lightning.parquet
│           ├── subs_seismic.parquet
│           ├── subs_tornado.parquet
│           ├── subs_wildfire.parquet
│           └── subs_wind.parquet
│
├── eaglei_outages/
│   ├── eaglei_outages_2014.csv
│   ├── eaglei_outages_2015.csv
│   ├── ...
│   └── eaglei_outages_2025.csv
│
└── NOAA storms data/
    └── StormEvents_details*.csv
```

## MHTran Hazard Exposure Data

The MHTran data contain substation-level natural-hazard exposure metrics.

The following files and metrics are used:

| Hazard | File | Exposure Metric |
|---|---|---|
| Flood | `subs_flood.parquet` | `flood_depth_m` |
| Freezing rain and wind gust (FZG) | `subs_fzg.parquet` | `R_50RP` |
| Hail | `subs_hail.parquet` | `hail_rate` |
| Landslide | `subs_landslide.parquet` | `ls_susc` |
| Lightning | `subs_lightning.parquet` | `ltng_rate` |
| Seismic | `subs_seismic.parquet` | `pga_475` |
| Tornado | `subs_tornado.parquet` | `annual_hit_rate` |
| Wildfire | `subs_wildfire.parquet` | `whp_cnt` |
| Wind | `subs_wind.parquet` | `pe_5pct_maxwind` |

Substation latitude and longitude coordinates are used to assign each substation to a U.S. county.

## EAGLE-I Outage Data

EAGLE-I data are used to calculate county-level outage frequency and severity.

The analysis uses annual files from 2014 through 2025.

For outage values:

- Files through 2023 use the `sum` column.
- Files from 2024 onward use the `customers_out` column.

County identifiers are obtained from the `fips_code` field and standardized to five-digit GEOIDs.

For outage frequency, only observations with positive outage values are counted.

For outage severity, recorded outage values are summed by county.

## NOAA Storm Events

NOAA Storm Events are used in Stage 2 to associate EAGLE-I outage observations with corresponding natural-hazard events.

Only county-level NOAA events are used.

The retained Stage 2 hazard mappings are:

| MHTran Hazard | NOAA Event Types | Attribution Window |
|---|---|---|
| Flood | Flood, Flash Flood, Heavy Rain | 24 hours |
| Hail | Hail | 24 hours |
| Landslide | Debris Flow | 48 hours |
| Lightning | Lightning | 24 hours |
| Tornado | Tornado | 24 hours |

An EAGLE-I observation is associated with a NOAA event when both occur in the same county and the outage timestamp falls between the NOAA event start time and the event end time plus the attribution window.

Duplicate outage observations matched to overlapping qualifying NOAA events are removed before county-level outage frequency is calculated.

## U.S. Census Bureau County Boundaries

County boundaries are loaded directly by the validation scripts from the U.S. Census Bureau 2022 TIGER/Line cartographic boundary dataset.

The county `GEOID` field is used as the common geographic identifier across datasets.

## Data Availability

Raw datasets are excluded from this repository.

Users who wish to reproduce the analysis should obtain the required datasets separately and place them in the directory structure shown above.

The repository `.gitignore` should exclude the local `Data/` directory so that raw or restricted data are not accidentally committed.