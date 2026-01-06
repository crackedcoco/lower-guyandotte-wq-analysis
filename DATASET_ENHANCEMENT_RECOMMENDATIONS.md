# Dataset Enhancement Recommendations
## Lower Guyandotte Water Quality Monitoring Program

This document outlines recommended improvements to enhance the analytical power of the water quality dataset for TMDL assessment and watershed management.

---

## 1. CRITICAL MISSING COLUMNS

### 1.1 Flow/Discharge Data
**Current Gap:** No streamflow data linked to samples
**Impact:** Cannot calculate loads, normalize concentrations, or perform flow-adjusted trend analysis

| Recommended Column | Type | Source | Purpose |
|-------------------|------|--------|---------|
| `DISCHARGE_CFS` | float | USGS gages / field measurement | Load calculations, flow-adjusted trends |
| `GAGE_ID` | string | USGS NWIS | Link to continuous flow data |
| `FLOW_CONDITION` | enum | Field observation | Baseflow/stormflow/snowmelt classification |
| `DAYS_SINCE_RAIN` | int | Weather data | Antecedent conditions |

**Why Critical:**
```
Load (lb/day) = Concentration (mg/L) × Flow (CFS) × 5.39

Without flow, you cannot:
- Calculate TMDL allocations
- Distinguish dilution from source reduction
- Normalize for wet vs. dry conditions
- Perform proper trend analysis (Weighted Regressions on Time, Discharge, Season - WRTDS)
```

### 1.2 Biological Assessment Data
**Current Gap:** No macroinvertebrate or fish community data
**Impact:** Cannot assess biological impairment directly (current data uses chemistry as surrogate)

| Recommended Column | Type | Source | Purpose |
|-------------------|------|--------|---------|
| `WVSCI_SCORE` | float | WV Stream Condition Index | Biological impairment determination |
| `GENUS_EPT_RICHNESS` | int | Benthic sampling | Sensitive taxa indicator |
| `HILSENHOFF_INDEX` | float | Biotic index | Organic pollution indicator |
| `FISH_IBI` | float | Fish Index of Biotic Integrity | Aquatic life use assessment |
| `HABITAT_SCORE` | float | Rapid Bioassessment Protocol | Physical habitat quality |

**Why Critical:**
- WV uses WVSCI < 68 as biological impairment threshold
- Many streams are "biologically impaired" without chemistry exceedances
- Biological data integrates stressor effects over time
- Required for weight-of-evidence impairment determinations

### 1.3 Land Use / Source Attribution
**Current Gap:** No linkage to pollution sources
**Impact:** Cannot attribute impairments to specific sources for TMDL allocations

| Recommended Column | Type | Source | Purpose |
|-------------------|------|--------|---------|
| `UPSTREAM_MINING_ACRES` | float | WVDEP permits | Surface mining impact quantification |
| `UPSTREAM_NPDES_COUNT` | int | ECHO/ICIS | Point source inventory |
| `NEAREST_NPDES_ID` | string | WVDEP | Link to discharge monitoring reports |
| `SEPTIC_DENSITY_PER_MI2` | float | Census/WVDHHR | Fecal coliform source assessment |
| `AGRICULTURAL_PCT` | float | NLCD | Nutrient/bacteria source |
| `IMPERVIOUS_PCT` | float | NLCD | Urban stormwater indicator |
| `FOREST_PCT` | float | NLCD | Reference condition indicator |

### 1.4 Sample Context
**Current Gap:** Limited information about sample representativeness
**Impact:** Difficulty interpreting anomalous values

| Recommended Column | Type | Source | Purpose |
|-------------------|------|--------|---------|
| `WEATHER_24HR` | enum | Weather data | Rain/dry classification |
| `PRECIP_INCHES_48HR` | float | PRISM/weather station | Storm event identification |
| `SAMPLE_PURPOSE` | enum | Sampler | Routine/storm/complaint/special study |
| `UPSTREAM_DISTURBANCE` | text | Field notes | Recent mining/construction/spills |
| `VISUAL_ASSESSMENT` | text | Field notes | Odor, color, algae, deposits |

---

## 2. SENSOR DATA INTEGRATION

### 2.1 Continuous Water Quality Sensors
**Benefit:** High-frequency data captures dynamics missed by grab samples

| Parameter | Sensor Type | Deployment | Value Added |
|-----------|-------------|------------|-------------|
| **Specific Conductance** | YSI/Hydrolab sonde | Continuous at key sites | Real-time mining impact detection, diel patterns |
| **pH** | Multi-parameter sonde | Continuous | Diel variation (photosynthesis), event detection |
| **Dissolved Oxygen** | Optical DO sensor | Continuous | Diel minimum (most critical), stratification |
| **Temperature** | Thermistor | Continuous | Thermal regime, climate change tracking |
| **Turbidity** | Optical turbidity | Continuous | Sediment transport, storm response, proxy for TSS |
| **Nitrate** | UV-Vis sensor (SUNA/S::CAN) | Continuous | Nutrient dynamics, storm loading |

**Data Columns Needed for Sensor Integration:**
```
SENSOR_ID          - Unique sensor identifier
MEASUREMENT_TYPE   - "continuous" vs "grab"
AGGREGATION_PERIOD - 15min, hourly, daily
QC_FLAG_SENSOR     - Sensor-specific QA codes
CALIBRATION_DATE   - Last calibration
SENSOR_DEPTH_M     - Deployment depth
```

### 2.2 Remote Sensing Integration
**Benefit:** Watershed-scale context, areas without monitoring

| Data Source | Parameters | Resolution | Application |
|-------------|------------|------------|-------------|
| **Landsat/Sentinel-2** | Chlorophyll-a, turbidity, TSS | 10-30m, 5-16 days | Reservoir monitoring, algae detection |
| **MODIS** | Surface temperature | 250m, daily | Thermal impairment screening |
| **Planet/SkySat** | Surface mining extent | 3m, daily | Real-time disturbance tracking |
| **LiDAR** | Riparian buffer width, channel morphology | 1m | Habitat assessment, erosion risk |

### 2.3 Weather Station Network
**Benefit:** Essential for flow estimation, stormwater analysis

| Parameter | Source | Integration |
|-----------|--------|-------------|
| Precipitation | CoCoRaHS, PRISM, NWS | Link each sample to nearest station |
| Air temperature | NOAA ISD | Seasonal analysis, ice-free periods |
| Evapotranspiration | gridMET | Water balance calculations |

---

## 3. STRUCTURAL IMPROVEMENTS

### 3.1 Separate Tables (Normalized Database)
Current: All data in one flat CSV
Recommended: Relational database structure

```
STATIONS table:
  - STATION_ID (PK)
  - WATERBODY_NAME
  - AU_ID
  - HUC12
  - LAT_DD, LON_DD
  - STATION_TYPE (ambient/effluent/reference)
  - DRAINAGE_AREA_MI2
  - ELEVATION_FT

SAMPLES table:
  - SAMPLE_ID (PK)
  - STATION_ID (FK)
  - SAMPLE_DATETIME
  - SAMPLE_PURPOSE
  - FLOW_CFS
  - WEATHER_CONDITION
  - SAMPLER_ID

RESULTS table:
  - RESULT_ID (PK)
  - SAMPLE_ID (FK)
  - PARAMETER
  - FRACTION
  - VALUE
  - UNIT
  - MDL
  - QL
  - QC_FLAG
  - LAB_ID
  - ANALYSIS_DATE

STANDARDS table:
  - PARAMETER
  - STANDARD_VALUE
  - STANDARD_TYPE
  - DESIGNATED_USE
  - REGULATORY_CITATION
  - EFFECTIVE_DATE
```

### 3.2 Method Tracking Improvements
**Current Gap:** ANL_METHOD is inconsistent, makes comparability difficult

| Recommended Column | Purpose |
|-------------------|---------|
| `METHOD_ID` | Standardized EPA method number |
| `METHOD_VERSION` | Track method updates over time |
| `LAB_ID` | Identify lab-specific biases |
| `HOLDING_TIME_EXCEEDED` | Flag samples with QC issues |
| `SAMPLE_CONTAINER` | Glass vs plastic affects metals |
| `PRESERVATION` | Acid/unpreserved/chilled |

### 3.3 Impairment Status Tracking
**Current Gap:** No link to regulatory decisions

| Recommended Column | Purpose |
|-------------------|---------|
| `IR_CYCLE` | Which Integrated Report (2016, 2018, 2020, etc.) |
| `LISTING_STATUS` | 303(d) listed, TMDL complete, delisted |
| `IMPAIRMENT_CAUSE` | ATTAINS cause code |
| `IMPAIRMENT_SOURCE` | ATTAINS source code |
| `TMDL_ID` | Link to approved TMDL document |

---

## 4. ADDITIONAL PARAMETERS TO MONITOR

### 4.1 Currently Under-Sampled (Add to Routine)
| Parameter | Current n | Recommended | Reason |
|-----------|-----------|-------------|--------|
| `E. coli` | 105 | All bacteria sites | Modern indicator replacing fecal coliform |
| `Chloride` | 2,483 | All SC sites | Distinguish mining types (surface vs deep) |
| `Bromide` | 1,432 | Mining-impacted sites | Marcellus shale wastewater tracer |
| `Selenium, Dissolved` | 62 | All Se sites | Bioavailable fraction |
| `Hardness` | 3,220 | All metals sites | Required for hardness-dependent criteria |

### 4.2 Emerging Contaminants
| Parameter | Why Important | Frequency |
|-----------|---------------|-----------|
| Microplastics | Emerging concern, unknown extent | Annual screening |
| PFAS | Industrial/landfill contamination | Targeted sites |
| Pharmaceuticals | WWTP indicator | WWTP-impacted streams |
| Coal combustion residuals | Ash pond leaching | Near power plants |

### 4.3 Nutrient Expansion
| Parameter | Current Status | Recommendation |
|-----------|----------------|----------------|
| Total N | 1,382 samples | Expand - needed for nutrient TMDLs |
| Total P | 1,490 samples | Expand - limiting nutrient |
| Chlorophyll-a | 457 samples | Expand to all wadeable streams |
| Orthophosphate | 343 samples | Add for bioavailable P |

---

## 5. QUALITY ASSURANCE ENHANCEMENTS

### 5.1 Required QA Columns
```
FIELD_DUPLICATE     - Was a field duplicate collected?
LAB_DUPLICATE       - Lab duplicate result
MATRIX_SPIKE_RECOVERY - MS/MSD recovery percent
BLANK_RESULT        - Associated blank value
RPD_DUPLICATE       - Relative percent difference
QC_BATCH_ID         - Link to lab QC batch
```

### 5.2 Data Validation Flags
Standardize FLAG_CODE to EPA-compatible codes:
```
U  = Non-detect at MDL
J  = Estimated (detected below QL)
R  = Rejected (QC failure)
B  = Blank contamination
H  = Holding time exceeded
E  = Estimated (other reason)
*  = Laboratory verified
```

---

## 6. IMPLEMENTATION PRIORITY

### Phase 1 (Immediate - High Impact)
1. Add `FLOW_CFS` or link to USGS gages
2. Add `IS_AMBIENT` flag (already implemented in v2 script)
3. Standardize `FLAG_CODE` vocabulary
4. Add `DRAINAGE_AREA_MI2` to stations

### Phase 2 (Near-term - TMDL Support)
1. Integrate biological data (WVSCI scores)
2. Add land use percentages per station
3. Link to NPDES permit database
4. Add precipitation data

### Phase 3 (Long-term - Program Enhancement)
1. Deploy continuous SC sensors at key sites
2. Integrate remote sensing for mining extent
3. Add emerging contaminants screening
4. Implement relational database structure

---

## 7. DATA SOURCES FOR ENHANCEMENT

| Data Type | Source | API/Access |
|-----------|--------|------------|
| Streamflow | USGS NWIS | `dataretrieval` Python package |
| Weather | PRISM / NOAA | `prism` R package, NOAA API |
| Land cover | NLCD | `rasterstats` + MRLC download |
| Mining permits | WVDEP | Public records request |
| NPDES data | EPA ECHO | ECHO API |
| Biological | WVDEP WAB | Data sharing agreement |
| Remote sensing | Google Earth Engine | GEE Python API |

---

## SUMMARY

The current dataset is valuable but limited by:
1. **No flow data** - Cannot calculate loads or normalize trends
2. **No biological data** - Must use chemistry as surrogate
3. **No source linkage** - Cannot attribute impairments
4. **Mixed ambient/effluent** - Requires careful filtering (now addressed in v2)
5. **Inconsistent QC flags** - Complicates data quality assessment

Addressing these gaps would transform the dataset from a **compliance monitoring record** to a **comprehensive watershed assessment tool** capable of supporting:
- Defensible TMDL calculations
- Source identification and allocation
- Trend analysis with proper controls
- Biological condition assessment
- Real-time impairment detection
- Climate change impact tracking
