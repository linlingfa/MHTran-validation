# MHTran Power Grid Outage Validation

This repository contains a three-stage validation of MHTran natural-hazard exposure metrics against historical U.S. power outage observations from EAGLE-I and NOAA Storm Events.

The analysis evaluates whether modeled natural-hazard exposure corresponds with observed county-level outage patterns across the United States.

## Research Questions

1. To what extent do individual MHTran hazard exposure metrics correspond with overall observed county-level power outage patterns?
2. To what extent do individual MHTran hazard exposure metrics correspond with historical outage patterns associated with the same hazards?
3. To what extent does an integrated MHTran multi-hazard score align with observed county-level power outage frequency, including spatial variation across the contiguous United States?

## Validation Framework

### Stage 1: Individual Hazard Exposure Validation

Each MHTran hazard is evaluated independently against overall EAGLE-I outage frequency and severity.

County exposure is represented by the maximum valid substation exposure value within each county.

Validation metrics:

- Spearman rank correlation
- Top-5% hotspot overlap

Hazards evaluated:

- Flood
- Freezing rain and wind gust (FZG)
- Hail
- Landslide
- Lightning
- Seismic
- Tornado
- Wildfire
- Wind

### Stage 2: Hazard-Attributed Outage Validation

EAGLE-I observations are associated with corresponding NOAA Storm Events using county location and hazard-specific temporal windows.

Both maximum and mean county exposure are evaluated.

Hazards retained for Stage 2:

- Flood
- Hail
- Landslide
- Lightning
- Tornado

Validation metrics:

- Spearman rank correlation
- Top-5% hotspot overlap

An EAGLE-I outage observation is attributed to a NOAA event when both occur in the same county and the outage timestamp falls between the NOAA event start time and the event end time plus the hazard-specific attribution window.

Duplicate EAGLE-I observations matched to overlapping NOAA events are removed before county-level outage frequency is calculated.

### Stage 3: Integrated Multi-Hazard Validation

County exposure for all nine hazards is converted to percentile ranks ranging from 0 to 1.

The composite multi-hazard exposure score is calculated as the mean percentile rank across hazards available for each county.

Missing hazard values remain unavailable rather than being assigned zero, and the number of hazards contributing to each county score is retained as `hazard_count`.

Validation methods:

- Spearman rank correlation
- Top-5% hotspot overlap
- Geographically Weighted Regression (GWR)

The GWR analysis evaluates spatial variation in the relationship between the composite multi-hazard exposure score and log-transformed outage frequency.

## Repository Structure

```text
MHTran-validation/
├── README.md
├── requirements.txt
├── .gitignore
│
├── scripts/
│   ├── stage1_validation.py
│   ├── stage2_validation.py
│   └── stage3_validation.py
│
├── data/
│   └── README.md
│
└── results/
    ├── figures/
    └── tables/
```

## Data

Raw datasets are not included in this repository.

The analysis uses:

- MHTran substation-level hazard exposure data
- EAGLE-I county-level outage observations
- NOAA Storm Events
- U.S. Census Bureau county boundaries

See `data/README.md` for the expected local data structure and filenames.

## Installation

Install required packages with:

```bash
pip install -r requirements.txt
```

## Running the Analysis

Run each stage from the repository root:

```bash
python scripts/stage1_validation.py
python scripts/stage2_validation.py
python scripts/stage3_validation.py
```

For Stages 1 and 2, change the `hazard_name` variable near the top of the script to select the hazard being evaluated.

### Execution Note

For reproducibility, the validation scripts should be run as standard Python files from the terminal rather than through the VS Code debugger.

Stage 3 in particular uses GWR bandwidth selection and parallel processing, which may encounter issues when executed through the debugger.