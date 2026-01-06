# Lower Guyandotte Watershed Water Quality Analysis

Comprehensive Python analysis tools for water quality monitoring data from the Lower Guyandotte River watershed (HUC 05070102) in West Virginia.

## Overview

This repository contains tools for analyzing TMDL (Total Maximum Daily Load) compliance and water quality trends in the Lower Guyandotte watershed, aligned with **West Virginia 47CSR2** water quality standards.

## Files

| File | Description |
|------|-------------|
| `water_quality_analysis.py` | Version 1.0 - Basic analysis script |
| `water_quality_analysis_v2.py` | Version 2.0 - Enhanced analysis with WV-specific standards |
| `DATASET_ENHANCEMENT_RECOMMENDATIONS.md` | Recommendations for improving the monitoring dataset |

## Key Features (v2.0)

- **WV 47CSR2 Standards** - Correct state-specific criteria (Fe=1.5 mg/L, pH 6.0-9.0, Al Dissolved=0.75 mg/L)
- **Ambient vs Point Source Separation** - Filters discharge monitoring from ambient stream data
- **Fecal Coliform Analysis** - Proper geometric mean + single sample max (200/400 CFU/100mL)
- **Kaplan-Meier for Non-Detects** - Improved handling of left-censored data
- **Assessment Unit Analysis** - Per-AU impairment determination for 303(d) listing
- **HUC12 Hotspot Prioritization** - Ranks subwatersheds by impairment severity
- **Mining Impact Signatures** - Detects AMD and surface mining indicators

## TMDL-Listed Parameters

| Parameter | WV Standard | Key Finding |
|-----------|-------------|-------------|
| Fecal Coliform | 200 geomean / 400 max | **30% exceedance rate** - primary impairment |
| Specific Conductance | 500 µS/cm benchmark | **+96% increase** since 1970s - mining impact |
| Fe, Total | 1.5 mg/L | **-76% decrease** - improving |
| Al, Dissolved | 0.75 mg/L | Generally compliant |
| pH | 6.0 - 9.0 | Generally compliant |

## Requirements

```bash
pip install pandas numpy scipy matplotlib seaborn
```

## Usage

```bash
# Place CSV data file in same directory
python water_quality_analysis_v2.py
```

Output is saved to `wq_analysis_output_v2/`:
- `summary_statistics.csv`
- `exceedance_analysis.csv`
- `trend_analysis.csv`
- `fecal_coliform_analysis.csv`
- `assessment_unit_analysis.csv`
- `huc12_hotspots.csv`
- `analysis_report.txt`
- `plots/` - Visualization PNGs

## Data Source

Data from West Virginia DEP ambient water quality monitoring program for the Lower Guyandotte watershed.

## License

MIT
