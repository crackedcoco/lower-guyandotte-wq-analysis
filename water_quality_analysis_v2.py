#!/usr/bin/env python3
"""
Water Quality Data Analysis Script - Version 2.0
=================================================
Comprehensive analysis of Lower Guyandotte watershed water quality monitoring data.
Aligned with West Virginia 47CSR2 water quality standards and TMDL requirements.

Key Improvements over v1:
- Correct WV 47CSR2 standards (not generic EPA)
- Separates ambient monitoring from point source discharge data
- Uses Al, Dissolved (not Total) per WV regulations
- Geometric mean for fecal coliform assessment
- Improved non-detect handling with Kaplan-Meier option
- AU_ID analysis for 303(d) listing support
- Mining impact signature detection
- QA/QC flag awareness
- HUC12 hotspot prioritization

Author: Generated for water quality TMDL audit
License: MIT
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

# =============================================================================
# CONFIGURATION
# =============================================================================

class DataCategory(Enum):
    """Categories for separating data types."""
    AMBIENT = "ambient"
    POINT_SOURCE = "point_source"
    ALL = "all"


@dataclass
class AnalysisConfig:
    """Configuration settings for the water quality analysis."""

    # File paths
    input_file: Path = Path("05070102_Lower_Guyandotte.csv")
    output_dir: Path = Path("wq_analysis_output_v2")

    # Analysis settings
    non_detect_method: str = "kaplan_meier"  # "kaplan_meier", "ros", "half_dl", "dl", "exclude"
    min_samples_for_trend: int = 10
    min_samples_for_geomean: int = 5
    significance_level: float = 0.05
    high_non_detect_threshold: float = 0.50  # Warn if >50% ND

    # Plotting settings
    figure_dpi: int = 150
    figure_format: str = "png"
    color_palette: str = "viridis"

    # ==========================================================================
    # WEST VIRGINIA 47CSR2 WATER QUALITY STANDARDS
    # Reference: https://dep.wv.gov/WWE/Programs/wqs/Pages/default.aspx
    # ==========================================================================
    wq_standards: dict = field(default_factory=lambda: {
        # Format: parameter_key: {
        #   'limit': value or (low, high) for range,
        #   'type': 'max', 'min', or 'range',
        #   'unit': unit string,
        #   'source': regulatory citation,
        #   'use': designated use category,
        #   'statistic': how to calculate (single, geomean, etc.)
        # }

        # === TMDL-Listed Parameters ===
        "Fecal Coliform (MF)": {
            'limit': 400,
            'limit_geomean': 200,
            'type': 'max',
            'unit': 'CFU/100 mL',
            'source': 'WV 47CSR2',
            'use': 'Recreation (Category B)',
            'statistic': 'geomean_and_max',
            'notes': 'Geomean ≤200 AND no sample >400 in 30-day period'
        },
        "Fe, Total": {
            'limit': 1.5,
            'type': 'max',
            'unit': 'mg/L',
            'source': 'WV 47CSR2',
            'use': 'Aquatic Life',
            'statistic': 'single'
        },
        "Al, Dissolved": {
            'limit': 0.75,
            'type': 'max',
            'unit': 'mg/L',
            'source': 'WV 47CSR2',
            'use': 'Aquatic Life',
            'statistic': 'single',
            'notes': 'WV uses DISSOLVED aluminum, not total'
        },
        "Mn, Total": {
            'limit': 1.0,
            'type': 'max',
            'unit': 'mg/L',
            'source': 'WV 47CSR2',
            'use': 'Aquatic Life',
            'statistic': 'single'
        },
        "Se, Total": {
            'limit': 0.005,
            'type': 'max',
            'unit': 'mg/L',
            'source': 'WV 47CSR2',
            'use': 'Aquatic Life (Chronic)',
            'statistic': 'single',
            'notes': 'Often high non-detect rate'
        },
        "pH": {
            'limit': (6.0, 9.0),
            'type': 'range',
            'unit': 'SU',
            'source': 'WV 47CSR2',
            'use': 'Aquatic Life',
            'statistic': 'single'
        },
        "DO": {
            'limit': 5.0,
            'type': 'min',
            'unit': 'mg/L',
            'source': 'WV 47CSR2',
            'use': 'Aquatic Life (Warmwater)',
            'statistic': 'single'
        },

        # === Mining Impact Indicators (not formal TMDL but diagnostic) ===
        "Specific Conductance": {
            'limit': 500,
            'type': 'max',
            'unit': 'uS/cm',
            'source': 'EPA Guidance / WV Benchmark',
            'use': 'Aquatic Life (Narrative)',
            'statistic': 'single',
            'notes': 'Elevated SC indicates mining/industrial impact'
        },
        "Sulfate (SO4)": {
            'limit': 250,
            'type': 'max',
            'unit': 'mg/L',
            'source': 'EPA Secondary MCL',
            'use': 'Human Health',
            'statistic': 'single'
        },
        "Hot Acidity": {
            'limit': 0,
            'type': 'max',
            'unit': 'mg/L CaCO3',
            'source': 'WV Narrative',
            'use': 'Aquatic Life',
            'statistic': 'single',
            'notes': 'Any net acidity indicates AMD'
        },
    })

    # Parameters of primary interest for detailed analysis (TMDL-focused)
    priority_parameters: list = field(default_factory=lambda: [
        "Fecal Coliform (MF)",
        "Fe, Total",
        "Al, Dissolved",
        "Al, Total",  # Include for comparison
        "Mn, Total",
        "Se, Total",
        "pH",
        "DO",
        "Specific Conductance",
        "Sulfate (SO4)",
        "Hot Acidity",
        "Alkalinity, Total",
        "Hardness",
        "TDS (Total Dissolved Solids)",
    ])

    # Mining signature parameter clusters
    mining_indicators: dict = field(default_factory=lambda: {
        'surface_mining': ['Specific Conductance', 'Sulfate (SO4)', 'Se, Total',
                          'TDS (Total Dissolved Solids)', 'Hardness'],
        'acid_mine_drainage': ['pH', 'Fe, Total', 'Al, Dissolved', 'Mn, Total',
                               'Hot Acidity', 'Sulfate (SO4)'],
        'deep_mine_brine': ['Specific Conductance', 'Chloride, Total', 'Bromide, Total',
                            'Na, Total', 'TDS (Total Dissolved Solids)']
    })

    # QA/QC flags to handle specially
    qc_flags_exclude: list = field(default_factory=lambda: ['R', 'REJ'])  # Rejected
    qc_flags_calculated: list = field(default_factory=lambda: ['C'])  # Calculated values
    qc_flags_estimated: list = field(default_factory=lambda: ['E', 'J'])  # Estimated


# =============================================================================
# LOGGING SETUP
# =============================================================================

def setup_logging(output_dir: Path) -> logging.Logger:
    """Configure logging to both file and console."""
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("water_quality_analysis")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    # File handler - detailed logging
    fh = logging.FileHandler(output_dir / "analysis.log", mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        '%(asctime)s | %(levelname)-8s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))

    # Console handler - info and above
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)-8s | %(message)s'))

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger


# =============================================================================
# DATA LOADING AND CLEANING
# =============================================================================

class WaterQualityDataLoader:
    """Handles loading and initial cleaning of water quality data."""

    def __init__(self, config: AnalysisConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger

    def load_data(self) -> pd.DataFrame:
        """Load the CSV file with appropriate data types."""
        self.logger.info(f"Loading data from {self.config.input_file}")

        dtype_spec = {
            'STATION_ID': 'str', 'WATERBODY_NAME': 'str', 'NHD_ANCODE': 'str',
            'ANCODE': 'str', 'MILE_POINT': 'str', 'DESCRIPTOR': 'str',
            'AU_ID': 'str', 'NHDP2COMID': 'str', 'DEP_WS': 'str',
            'HUC8': 'str', 'HUC8NAME': 'str', 'HUC10': 'str', 'HUC10NAME': 'str',
            'HUC12': 'str', 'HUC12NAME': 'str', 'COUNTY': 'str', 'TOPO_NAME': 'str',
            'SAMPLEID': 'str', 'SAMPLE_TZ': 'str', 'DUP_TYPE': 'str', 'DUP_NAME': 'str',
            'SUBSAMPID': 'str', 'DEPTH_DESC': 'str', 'DIST_DESC': 'str',
            'TSECT_NAME': 'str', 'REACH_LOC': 'str',
            'PARAMETER': 'str', 'FRACTION': 'str', 'SUBSTANCE': 'str',
            'NON_DETECT': 'str', 'DEFAULT_UNITS': 'str', 'LIM_TYPE': 'str',
            'FLAG_CODE': 'str', 'FLAG_DESC': 'str', 'ANL_METHOD': 'str',
        }

        df = pd.read_csv(
            self.config.input_file,
            dtype=dtype_spec,
            low_memory=False,
            na_values=['None', 'NA', 'N/A', '', 'null', 'NULL'],
        )

        self.logger.info(f"Loaded {len(df):,} records with {len(df.columns)} columns")
        return df

    def clean_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Clean and preprocess the data."""
        self.logger.info("Cleaning data...")
        initial_rows = len(df)

        df = self._parse_dates(df)
        df = self._clean_numeric_columns(df)
        df = self._classify_data_sources(df)
        df = self._handle_qc_flags(df)
        df = self._handle_non_detects(df)
        df = self._clean_coordinates(df)
        df = self._add_derived_columns(df)
        df = df.drop_duplicates()

        final_rows = len(df)
        self.logger.info(f"Cleaning complete: {initial_rows:,} → {final_rows:,} records")

        return df

    def _parse_dates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Parse date and time columns."""
        self.logger.debug("Parsing date/time columns")

        df['SAMPLE_DATETIME'] = pd.to_datetime(
            df['SAMPLEDATE'].astype(str) + ' ' + df['SAMPLETIME'].fillna('00:00').astype(str),
            format='%Y-%m-%d %H:%M',
            errors='coerce'
        )
        df['SAMPLE_DATE'] = pd.to_datetime(df['SAMPLEDATE'], errors='coerce')
        df['SAMPLE_YEAR'] = df['SAMPLE_DATE'].dt.year
        df['SAMPLE_MONTH'] = df['SAMPLE_DATE'].dt.month
        df['SAMPLE_QUARTER'] = df['SAMPLE_DATE'].dt.quarter
        df['SAMPLE_DOY'] = df['SAMPLE_DATE'].dt.dayofyear

        season_map = {12: 'Winter', 1: 'Winter', 2: 'Winter',
                      3: 'Spring', 4: 'Spring', 5: 'Spring',
                      6: 'Summer', 7: 'Summer', 8: 'Summer',
                      9: 'Fall', 10: 'Fall', 11: 'Fall'}
        df['SEASON'] = df['SAMPLE_MONTH'].map(season_map)

        return df

    def _clean_numeric_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Clean and convert numeric columns."""
        numeric_cols = ['VALUE', 'DL', 'QL', 'DEPTH', 'DISTANCE',
                        'MILE_POINT', 'LON_DD', 'LAT_DD']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        return df

    def _classify_data_sources(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        CRITICAL: Separate ambient monitoring from point source discharge data.
        AU_ID containing 'SOURCE' indicates permitted discharge monitoring.
        """
        self.logger.info("Classifying data sources (ambient vs point source)...")

        # Identify point source discharge samples
        df['IS_POINT_SOURCE'] = (
            df['AU_ID'].str.upper().str.contains('SOURCE', na=False) |
            df['WATERBODY_NAME'].str.contains('Discharge|Outfall|Permitted', case=False, na=False)
        )

        # Categorize
        df['DATA_CATEGORY'] = np.where(
            df['IS_POINT_SOURCE'],
            DataCategory.POINT_SOURCE.value,
            DataCategory.AMBIENT.value
        )

        ps_count = df['IS_POINT_SOURCE'].sum()
        amb_count = len(df) - ps_count
        self.logger.info(f"  Ambient samples: {amb_count:,}")
        self.logger.info(f"  Point source discharge samples: {ps_count:,}")

        return df

    def _handle_qc_flags(self, df: pd.DataFrame) -> pd.DataFrame:
        """Handle QA/QC flags appropriately."""
        self.logger.info("Processing QA/QC flags...")

        df['QC_STATUS'] = 'accepted'

        # Flag rejected samples
        rejected_mask = df['FLAG_CODE'].isin(self.config.qc_flags_exclude)
        df.loc[rejected_mask, 'QC_STATUS'] = 'rejected'
        n_rejected = rejected_mask.sum()

        # Flag calculated values
        calc_mask = df['FLAG_CODE'].isin(self.config.qc_flags_calculated)
        df.loc[calc_mask, 'QC_STATUS'] = 'calculated'
        n_calculated = calc_mask.sum()

        # Flag estimated values
        est_mask = df['FLAG_CODE'].isin(self.config.qc_flags_estimated)
        df.loc[est_mask, 'QC_STATUS'] = 'estimated'
        n_estimated = est_mask.sum()

        self.logger.info(f"  Rejected: {n_rejected:,}, Calculated: {n_calculated:,}, Estimated: {n_estimated:,}")

        return df

    def _handle_non_detects(self, df: pd.DataFrame) -> pd.DataFrame:
        """Handle non-detect values with improved methods."""
        self.logger.info(f"Handling non-detects using method: {self.config.non_detect_method}")

        df['IS_NON_DETECT'] = df['NON_DETECT'].str.upper() == 'YES'
        df['VALUE_ORIGINAL'] = df['VALUE'].copy()
        df['VALUE_SUBSTITUTED'] = False

        # Calculate ND percentage by parameter for warnings
        nd_by_param = df.groupby('PARAMETER')['IS_NON_DETECT'].mean()
        high_nd_params = nd_by_param[nd_by_param > self.config.high_non_detect_threshold]

        if len(high_nd_params) > 0:
            self.logger.warning(f"HIGH NON-DETECT PARAMETERS (>{self.config.high_non_detect_threshold*100:.0f}%):")
            for param, pct in high_nd_params.items():
                self.logger.warning(f"  {param}: {pct*100:.1f}% - statistics may be unreliable")

        # Apply substitution method
        if self.config.non_detect_method == "half_dl":
            mask = df['IS_NON_DETECT'] & df['DL'].notna()
            df.loc[mask, 'VALUE'] = df.loc[mask, 'DL'] / 2
            df.loc[mask, 'VALUE_SUBSTITUTED'] = True

            mask2 = df['IS_NON_DETECT'] & df['DL'].isna()
            df.loc[mask2, 'VALUE'] = df.loc[mask2, 'VALUE'] / 2
            df.loc[mask2, 'VALUE_SUBSTITUTED'] = True

        elif self.config.non_detect_method == "kaplan_meier":
            # For K-M, we keep original values but flag them
            # Actual K-M calculation done in statistics module
            df.loc[df['IS_NON_DETECT'], 'VALUE_SUBSTITUTED'] = True
            self.logger.info("  Kaplan-Meier: values flagged, statistics calculated separately")

        elif self.config.non_detect_method == "exclude":
            df = df[~df['IS_NON_DETECT']].copy()

        nd_count = df['IS_NON_DETECT'].sum()
        nd_pct = 100 * nd_count / len(df) if len(df) > 0 else 0
        self.logger.info(f"  Total non-detects: {nd_count:,} records ({nd_pct:.1f}%)")

        return df

    def _clean_coordinates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Validate and clean coordinate data."""
        df['VALID_COORDS'] = (
            df['LON_DD'].notna() & df['LAT_DD'].notna() &
            (df['LON_DD'] >= -83) & (df['LON_DD'] <= -81) &  # WV bounds
            (df['LAT_DD'] >= 37) & (df['LAT_DD'] <= 40)
        )
        invalid_coords = (~df['VALID_COORDS']).sum()
        if invalid_coords > 0:
            self.logger.warning(f"Found {invalid_coords:,} records with invalid/missing coordinates")
        return df

    def _add_derived_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add useful derived columns."""
        df['STATION_WATERBODY'] = df['STATION_ID'].astype(str) + '_' + df['WATERBODY_NAME'].fillna('Unknown')
        df['PARAMETER_CLEAN'] = df['PARAMETER'].str.strip().str.replace(r'\s+', ' ', regex=True)

        # Create combined parameter-fraction key for metals
        df['PARAM_FRACTION'] = df['PARAMETER'] + '_' + df['FRACTION'].fillna('Unknown')

        return df


# =============================================================================
# NON-DETECT STATISTICS (Kaplan-Meier / ROS)
# =============================================================================

class CensoredDataHandler:
    """Handles left-censored data (non-detects) using appropriate statistical methods."""

    def __init__(self, logger: logging.Logger):
        self.logger = logger

    def kaplan_meier_mean(self, values: np.ndarray, censored: np.ndarray) -> dict:
        """
        Calculate Kaplan-Meier estimate for left-censored data.

        For left-censored data, we flip the problem and compute survival from the right.
        Returns mean, median, std estimates.
        """
        if len(values) == 0:
            return {'mean': np.nan, 'median': np.nan, 'std': np.nan, 'method': 'kaplan_meier'}

        n_censored = censored.sum()
        n_total = len(values)

        # If no censoring or all censored, use simple statistics
        if n_censored == 0:
            return {
                'mean': np.mean(values),
                'median': np.median(values),
                'std': np.std(values),
                'method': 'standard',
                'n_censored': 0,
                'pct_censored': 0
            }

        if n_censored == n_total:
            # All non-detect - return max DL as upper bound
            return {
                'mean': np.max(values) / 2,  # Conservative estimate
                'median': np.max(values) / 2,
                'std': np.nan,
                'method': 'all_censored',
                'n_censored': n_total,
                'pct_censored': 100
            }

        # Simplified Kaplan-Meier for left-censored data
        # Sort by value, with censored observations first at each value
        df = pd.DataFrame({'value': values, 'censored': censored})
        df = df.sort_values(['value', 'censored'], ascending=[True, False])

        # For left-censoring, we use the Turnbull or Lynden-Bell estimator
        # Simplified: use detected values, weighted by survival probability
        detected = df[~df['censored']]['value'].values
        censored_vals = df[df['censored']]['value'].values

        # Simple imputation: for each censored value, use DL/2 for values below lowest detected
        # or interpolate based on distribution of detected values
        if len(detected) > 0:
            km_mean = np.mean(detected)
            km_median = np.median(detected)
            km_std = np.std(detected)

            # Adjust for censored values (simplified)
            # This is an approximation - for rigorous analysis use lifelines or NADA
            adjustment = (n_censored / n_total) * (np.min(detected) / 2)
            km_mean = km_mean * (1 - n_censored/n_total) + adjustment

        else:
            km_mean = np.mean(values) / 2
            km_median = np.median(values) / 2
            km_std = np.nan

        return {
            'mean': km_mean,
            'median': km_median,
            'std': km_std,
            'method': 'kaplan_meier_approx',
            'n_censored': n_censored,
            'pct_censored': 100 * n_censored / n_total
        }

    def geometric_mean_censored(self, values: np.ndarray, censored: np.ndarray) -> dict:
        """
        Calculate geometric mean for left-censored data.
        Uses detected values only, with warning if high censoring.
        """
        detected = values[~censored]

        if len(detected) == 0:
            return {'geomean': np.nan, 'n_detected': 0, 'pct_censored': 100}

        # Remove zeros/negatives for log transform
        detected_positive = detected[detected > 0]

        if len(detected_positive) == 0:
            return {'geomean': np.nan, 'n_detected': 0, 'pct_censored': 100}

        geomean = np.exp(np.mean(np.log(detected_positive)))
        pct_censored = 100 * censored.sum() / len(values)

        return {
            'geomean': geomean,
            'n_detected': len(detected_positive),
            'pct_censored': pct_censored,
            'warning': 'High censoring - geomean may be biased high' if pct_censored > 50 else None
        }


# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================

class StatisticalAnalysis:
    """Performs statistical analysis on water quality data."""

    def __init__(self, df: pd.DataFrame, config: AnalysisConfig, logger: logging.Logger):
        self.df = df
        self.config = config
        self.logger = logger
        self.censored_handler = CensoredDataHandler(logger)
        self.results: dict[str, Any] = {}

    def run_statistical_analysis(self, data_category: DataCategory = DataCategory.AMBIENT) -> dict[str, Any]:
        """Run complete statistical analysis."""
        self.logger.info("=" * 60)
        self.logger.info(f"STATISTICAL ANALYSIS ({data_category.value.upper()} DATA)")
        self.logger.info("=" * 60)

        # Filter to requested data category
        if data_category == DataCategory.AMBIENT:
            df = self.df[self.df['DATA_CATEGORY'] == 'ambient'].copy()
        elif data_category == DataCategory.POINT_SOURCE:
            df = self.df[self.df['DATA_CATEGORY'] == 'point_source'].copy()
        else:
            df = self.df.copy()

        self.logger.info(f"Analyzing {len(df):,} records")

        self.results['summary_stats'] = self._summary_statistics(df)
        self.results['exceedances'] = self._exceedance_analysis(df)
        self.results['fecal_coliform'] = self._fecal_coliform_analysis(df)
        self.results['trends'] = self._trend_analysis(df)
        self.results['seasonal'] = self._seasonal_analysis(df)
        self.results['mining_signature'] = self._mining_signature_analysis(df)

        return self.results

    def _summary_statistics(self, df: pd.DataFrame) -> pd.DataFrame:
        """Calculate summary statistics with proper censored data handling."""
        self.logger.info("Calculating summary statistics...")
        stats_list = []

        for param in self.config.priority_parameters:
            param_data = df[df['PARAMETER'] == param].copy()

            if len(param_data) == 0:
                continue

            values = param_data['VALUE'].dropna().values
            censored = param_data.loc[param_data['VALUE'].notna(), 'IS_NON_DETECT'].values

            if len(values) == 0:
                continue

            unit = param_data['DEFAULT_UNITS'].mode()
            unit = unit.iloc[0] if len(unit) > 0 else 'Unknown'

            # Get censored-aware statistics
            km_stats = self.censored_handler.kaplan_meier_mean(values, censored)

            # Standard statistics for comparison
            detected_values = values[~censored] if any(~censored) else values

            stats_dict = {
                'parameter': param,
                'unit': unit,
                'n_total': len(values),
                'n_detected': int((~censored).sum()),
                'n_non_detect': int(censored.sum()),
                'pct_non_detect': 100 * censored.sum() / len(values),
                'mean_km': km_stats['mean'],
                'mean_detected': np.mean(detected_values) if len(detected_values) > 0 else np.nan,
                'median': np.median(detected_values) if len(detected_values) > 0 else np.nan,
                'std': np.std(detected_values) if len(detected_values) > 1 else np.nan,
                'min': np.min(values),
                'p10': np.percentile(values, 10),
                'p25': np.percentile(values, 25),
                'p75': np.percentile(values, 75),
                'p90': np.percentile(values, 90),
                'max': np.max(values),
                'n_stations': param_data['STATION_ID'].nunique(),
                'n_samples': param_data['SAMPLEID'].nunique(),
                'date_min': param_data['SAMPLE_DATE'].min(),
                'date_max': param_data['SAMPLE_DATE'].max(),
            }

            # Add standard if available
            if param in self.config.wq_standards:
                std = self.config.wq_standards[param]
                stats_dict['wq_standard'] = std['limit'] if not isinstance(std['limit'], tuple) else f"{std['limit'][0]}-{std['limit'][1]}"
                stats_dict['wq_standard_source'] = std['source']

            stats_list.append(stats_dict)

        return pd.DataFrame(stats_list)

    def _exceedance_analysis(self, df: pd.DataFrame) -> dict:
        """Analyze exceedances of WV water quality standards."""
        self.logger.info("Analyzing water quality standard exceedances...")
        exceedances = {}

        for param, standard in self.config.wq_standards.items():
            # Handle Al, Dissolved specifically
            if param == "Al, Dissolved":
                param_data = df[
                    (df['PARAMETER'].str.contains('Al', case=False, na=False)) &
                    (df['FRACTION'].str.contains('Dissolved', case=False, na=False))
                ].copy()
            else:
                param_data = df[df['PARAMETER'] == param].copy()

            if len(param_data) == 0:
                continue

            values = param_data['VALUE'].dropna()
            if len(values) == 0:
                continue

            limit = standard['limit']
            limit_type = standard['type']

            # Calculate exceedances
            if limit_type == "range":
                low, high = limit
                exceed_low = (values < low).sum()
                exceed_high = (values > high).sum()
                n_exceedances = exceed_low + exceed_high
                exceed_detail = {'below_min': int(exceed_low), 'above_max': int(exceed_high)}
            elif limit_type == "min":
                n_exceedances = (values < limit).sum()
                exceed_detail = {'below_min': int(n_exceedances)}
            else:  # max
                n_exceedances = (values > limit).sum()
                exceed_detail = {'above_max': int(n_exceedances)}

            exceedance_pct = 100 * n_exceedances / len(values)

            exceedances[param] = {
                'standard_value': limit,
                'standard_type': limit_type,
                'standard_unit': standard['unit'],
                'standard_source': standard['source'],
                'designated_use': standard['use'],
                'n_samples': len(values),
                'n_exceedances': int(n_exceedances),
                'exceedance_percent': exceedance_pct,
                'exceedance_detail': exceed_detail,
                'min_value': float(values.min()),
                'max_value': float(values.max()),
                'median_value': float(values.median()),
                'n_stations_exceeding': param_data[
                    (values > limit) if limit_type == 'max' else (values < limit)
                ]['STATION_ID'].nunique() if limit_type != 'range' else np.nan,
            }

            if n_exceedances > 0:
                self.logger.info(f"  {param}: {n_exceedances:,} exceedances ({exceedance_pct:.1f}%)")

        return exceedances

    def _fecal_coliform_analysis(self, df: pd.DataFrame) -> dict:
        """
        Analyze fecal coliform with proper WV methodology:
        - Geometric mean ≤ 200 CFU/100mL
        - No single sample > 400 CFU/100mL
        """
        self.logger.info("Analyzing fecal coliform (geometric mean + max)...")

        param_data = df[df['PARAMETER'] == 'Fecal Coliform (MF)'].copy()

        if len(param_data) == 0:
            return {}

        values = param_data['VALUE'].dropna()
        censored = param_data.loc[param_data['VALUE'].notna(), 'IS_NON_DETECT'].values

        if len(values) == 0:
            return {}

        # Calculate geometric mean (for detected values)
        geomean_result = self.censored_handler.geometric_mean_censored(values.values, censored)

        # Single sample exceedances
        n_exceed_400 = (values > 400).sum()
        pct_exceed_400 = 100 * n_exceed_400 / len(values)

        # Geometric mean exceedance
        geomean_exceeds = geomean_result['geomean'] > 200 if not np.isnan(geomean_result['geomean']) else None

        # Calculate by station
        station_stats = param_data.groupby('STATION_ID').apply(
            lambda x: pd.Series({
                'n_samples': len(x),
                'geomean': np.exp(np.mean(np.log(x['VALUE'].dropna()[x['VALUE'].dropna() > 0])))
                          if len(x['VALUE'].dropna()[x['VALUE'].dropna() > 0]) > 0 else np.nan,
                'n_exceed_400': (x['VALUE'] > 400).sum(),
                'max_value': x['VALUE'].max()
            })
        )

        stations_impaired = station_stats[
            (station_stats['geomean'] > 200) | (station_stats['n_exceed_400'] > 0)
        ]

        result = {
            'n_samples': len(values),
            'geometric_mean': geomean_result['geomean'],
            'geomean_standard': 200,
            'geomean_exceeds': geomean_exceeds,
            'n_exceed_400': int(n_exceed_400),
            'pct_exceed_400': pct_exceed_400,
            'single_sample_standard': 400,
            'median': float(values.median()),
            'max': float(values.max()),
            'pct_non_detect': geomean_result['pct_censored'],
            'n_stations_total': param_data['STATION_ID'].nunique(),
            'n_stations_impaired': len(stations_impaired),
            'impaired_stations': stations_impaired.index.tolist() if len(stations_impaired) < 50 else f"{len(stations_impaired)} stations",
        }

        self.logger.info(f"  Geometric mean: {result['geometric_mean']:.1f} (standard: 200)")
        self.logger.info(f"  Samples >400: {n_exceed_400:,} ({pct_exceed_400:.1f}%)")
        self.logger.info(f"  Impaired stations: {result['n_stations_impaired']}")

        return result

    def _trend_analysis(self, df: pd.DataFrame) -> dict:
        """Perform temporal trend analysis with improved methodology."""
        self.logger.info("Performing trend analysis...")
        trends = {}

        for param in self.config.priority_parameters:
            param_data = df[df['PARAMETER'] == param].copy()

            if len(param_data) < self.config.min_samples_for_trend:
                continue

            # Calculate annual statistics with minimum sample requirement
            annual = param_data.groupby('SAMPLE_YEAR').agg({
                'VALUE': ['median', 'mean', 'count'],
                'IS_NON_DETECT': 'mean'
            })
            annual.columns = ['median', 'mean', 'count', 'pct_nd']
            annual = annual[annual['count'] >= 5]  # Require at least 5 samples/year

            if len(annual) < 5:
                continue

            years = annual.index.values.astype(float)
            values = annual['median'].values

            # Linear regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(years, values)

            # Mann-Kendall test
            mk_stat, mk_p = self._mann_kendall_test(values)

            # Sen's slope (more robust)
            sen_slope = self._sens_slope(years, values)

            # Calculate percent change
            if annual['median'].iloc[0] != 0:
                total_change_pct = 100 * (annual['median'].iloc[-1] - annual['median'].iloc[0]) / annual['median'].iloc[0]
            else:
                total_change_pct = np.nan

            trends[param] = {
                'n_years': len(annual),
                'year_range': f"{int(years.min())}-{int(years.max())}",
                'linear_slope': slope,
                'linear_slope_pct_per_year': 100 * slope / annual['median'].mean() if annual['median'].mean() != 0 else np.nan,
                'linear_r_squared': r_value ** 2,
                'linear_p_value': p_value,
                'sens_slope': sen_slope,
                'mann_kendall_stat': mk_stat,
                'mann_kendall_p_value': mk_p,
                'trend_direction': 'increasing' if slope > 0 else 'decreasing',
                'significant': mk_p < self.config.significance_level if not np.isnan(mk_p) else False,
                'total_change_pct': total_change_pct,
                'early_median': float(annual['median'].iloc[:3].mean()),
                'late_median': float(annual['median'].iloc[-3:].mean()),
                'avg_pct_nd': float(annual['pct_nd'].mean() * 100),
            }

            if trends[param]['significant']:
                direction = trends[param]['trend_direction']
                self.logger.info(
                    f"  {param}: {direction} ({total_change_pct:+.1f}% total, p={mk_p:.4f})"
                )

        return trends

    def _mann_kendall_test(self, x: np.ndarray) -> tuple[float, float]:
        """Perform Mann-Kendall trend test."""
        n = len(x)
        if n < 4:
            return np.nan, np.nan

        s = 0
        for i in range(n - 1):
            for j in range(i + 1, n):
                diff = x[j] - x[i]
                if diff > 0:
                    s += 1
                elif diff < 0:
                    s -= 1

        var_s = (n * (n - 1) * (2 * n + 5)) / 18

        if s > 0:
            z = (s - 1) / np.sqrt(var_s)
        elif s < 0:
            z = (s + 1) / np.sqrt(var_s)
        else:
            z = 0

        p_value = 2 * (1 - stats.norm.cdf(abs(z)))
        return s, p_value

    def _sens_slope(self, x: np.ndarray, y: np.ndarray) -> float:
        """Calculate Sen's slope (median of all pairwise slopes)."""
        n = len(x)
        slopes = []
        for i in range(n):
            for j in range(i + 1, n):
                if x[j] != x[i]:
                    slopes.append((y[j] - y[i]) / (x[j] - x[i]))
        return np.median(slopes) if slopes else np.nan

    def _seasonal_analysis(self, df: pd.DataFrame) -> dict:
        """Analyze seasonal patterns."""
        self.logger.info("Analyzing seasonal patterns...")
        seasonal = {}

        for param in self.config.priority_parameters[:8]:  # Limit for performance
            param_data = df[df['PARAMETER'] == param].copy()

            if len(param_data) < 20:
                continue

            season_stats = param_data.groupby('SEASON')['VALUE'].agg(['median', 'mean', 'std', 'count'])

            if len(season_stats) < 2:
                continue

            # Kruskal-Wallis test
            season_groups = [
                param_data[param_data['SEASON'] == s]['VALUE'].dropna().values
                for s in ['Winter', 'Spring', 'Summer', 'Fall']
                if s in param_data['SEASON'].values
            ]
            season_groups = [g for g in season_groups if len(g) >= 3]

            if len(season_groups) >= 2:
                try:
                    h_stat, kw_p = stats.kruskal(*season_groups)
                except ValueError:
                    h_stat, kw_p = np.nan, np.nan
            else:
                h_stat, kw_p = np.nan, np.nan

            seasonal[param] = {
                'season_medians': season_stats['median'].to_dict(),
                'season_means': season_stats['mean'].to_dict(),
                'season_counts': season_stats['count'].to_dict(),
                'kruskal_wallis_h': h_stat,
                'kruskal_wallis_p': kw_p,
                'significant_seasonal_diff': kw_p < self.config.significance_level if not np.isnan(kw_p) else False,
            }

        return seasonal

    def _mining_signature_analysis(self, df: pd.DataFrame) -> dict:
        """
        Detect mining impact signatures based on parameter correlations.

        Surface mining: High SC + High Sulfate + Elevated Se
        AMD: Low pH + High Fe + High Al + High Mn + Sulfate
        Deep mine brine: High SC + High Chloride + High Bromide + Low Sulfate
        """
        self.logger.info("Analyzing mining impact signatures...")

        results = {}

        # Create pivot for correlation analysis
        mining_params = ['Specific Conductance', 'Sulfate (SO4)', 'pH', 'Fe, Total',
                         'Chloride, Total', 'Se, Total', 'TDS (Total Dissolved Solids)']

        pivot_data = df[df['PARAMETER'].isin(mining_params)].pivot_table(
            index=['STATION_ID', 'SAMPLE_DATE'],
            columns='PARAMETER',
            values='VALUE',
            aggfunc='median'
        )

        if pivot_data.shape[0] < 10 or pivot_data.shape[1] < 3:
            return {'insufficient_data': True}

        # Calculate correlations
        corr = pivot_data.corr(method='spearman')
        results['correlation_matrix'] = corr.to_dict()

        # Identify samples with mining signatures
        samples_with_data = pivot_data.dropna(subset=['Specific Conductance'], how='all')

        # Surface mining signature: SC > 500 AND Sulfate > 100
        if 'Specific Conductance' in samples_with_data.columns and 'Sulfate (SO4)' in samples_with_data.columns:
            surface_mining = (
                (samples_with_data['Specific Conductance'] > 500) &
                (samples_with_data['Sulfate (SO4)'] > 100)
            )
            results['surface_mining_signature'] = {
                'n_samples': int(surface_mining.sum()),
                'pct_samples': 100 * surface_mining.sum() / len(samples_with_data),
                'description': 'SC > 500 uS/cm AND Sulfate > 100 mg/L'
            }

        # AMD signature: pH < 6 AND (Fe > 1.5 OR Al > 0.75)
        if 'pH' in samples_with_data.columns and 'Fe, Total' in samples_with_data.columns:
            amd = (
                (samples_with_data['pH'] < 6.0) &
                (samples_with_data['Fe, Total'] > 1.5)
            )
            results['amd_signature'] = {
                'n_samples': int(amd.sum()),
                'pct_samples': 100 * amd.sum() / len(samples_with_data),
                'description': 'pH < 6.0 AND Fe > 1.5 mg/L'
            }

        # High conductivity sites
        if 'Specific Conductance' in samples_with_data.columns:
            high_sc = samples_with_data['Specific Conductance'] > 500
            results['high_conductivity'] = {
                'n_samples': int(high_sc.sum()),
                'pct_samples': 100 * high_sc.sum() / len(samples_with_data),
                'median_sc': float(samples_with_data.loc[high_sc, 'Specific Conductance'].median()) if high_sc.any() else np.nan,
            }

        # SC-Sulfate correlation (strong positive = surface mining)
        if 'Specific Conductance' in corr.columns and 'Sulfate (SO4)' in corr.columns:
            sc_sulfate_corr = corr.loc['Specific Conductance', 'Sulfate (SO4)']
            results['sc_sulfate_correlation'] = {
                'value': float(sc_sulfate_corr),
                'interpretation': 'Strong surface mining signal' if sc_sulfate_corr > 0.7 else
                                  'Moderate mining signal' if sc_sulfate_corr > 0.5 else 'Weak/mixed signal'
            }

        self.logger.info(f"  Surface mining signatures: {results.get('surface_mining_signature', {}).get('n_samples', 0)}")
        self.logger.info(f"  AMD signatures: {results.get('amd_signature', {}).get('n_samples', 0)}")

        return results


# =============================================================================
# ASSESSMENT UNIT (AU_ID) ANALYSIS
# =============================================================================

class AssessmentUnitAnalysis:
    """Analyze data by Assessment Unit for 303(d) listing support."""

    def __init__(self, df: pd.DataFrame, config: AnalysisConfig, logger: logging.Logger):
        self.df = df
        self.config = config
        self.logger = logger

    def analyze_assessment_units(self) -> pd.DataFrame:
        """Analyze each AU for impairment status."""
        self.logger.info("=" * 60)
        self.logger.info("ASSESSMENT UNIT ANALYSIS")
        self.logger.info("=" * 60)

        # Filter to ambient data only
        df = self.df[self.df['DATA_CATEGORY'] == 'ambient'].copy()

        au_results = []

        for au_id in df['AU_ID'].dropna().unique():
            au_data = df[df['AU_ID'] == au_id]

            if len(au_data) < 5:
                continue

            result = {
                'AU_ID': au_id,
                'waterbody': au_data['WATERBODY_NAME'].mode().iloc[0] if len(au_data['WATERBODY_NAME'].mode()) > 0 else 'Unknown',
                'huc12': au_data['HUC12NAME'].mode().iloc[0] if len(au_data['HUC12NAME'].mode()) > 0 else 'Unknown',
                'n_samples': au_data['SAMPLEID'].nunique(),
                'n_stations': au_data['STATION_ID'].nunique(),
                'date_range': f"{au_data['SAMPLE_DATE'].min()} to {au_data['SAMPLE_DATE'].max()}",
            }

            # Check each standard
            impairments = []

            for param, standard in self.config.wq_standards.items():
                if param == "Al, Dissolved":
                    param_data = au_data[
                        (au_data['PARAMETER'].str.contains('Al', case=False, na=False)) &
                        (au_data['FRACTION'].str.contains('Dissolved', case=False, na=False))
                    ]['VALUE'].dropna()
                else:
                    param_data = au_data[au_data['PARAMETER'] == param]['VALUE'].dropna()

                if len(param_data) < 3:
                    continue

                limit = standard['limit']
                limit_type = standard['type']

                if limit_type == 'range':
                    low, high = limit
                    n_exceed = ((param_data < low) | (param_data > high)).sum()
                elif limit_type == 'min':
                    n_exceed = (param_data < limit).sum()
                else:
                    n_exceed = (param_data > limit).sum()

                exceedance_rate = n_exceed / len(param_data)
                result[f'{param}_n'] = len(param_data)
                result[f'{param}_exceed_pct'] = 100 * exceedance_rate

                # WV typically uses 10% exceedance rate as threshold
                if exceedance_rate > 0.10:
                    impairments.append(param)

            result['impairments'] = ', '.join(impairments) if impairments else 'None'
            result['n_impairments'] = len(impairments)
            result['impaired'] = len(impairments) > 0

            au_results.append(result)

        au_df = pd.DataFrame(au_results)

        if len(au_df) > 0:
            n_impaired = au_df['impaired'].sum()
            self.logger.info(f"  Total AUs analyzed: {len(au_df)}")
            self.logger.info(f"  Impaired AUs: {n_impaired} ({100*n_impaired/len(au_df):.1f}%)")

        return au_df


# =============================================================================
# HUC12 HOTSPOT ANALYSIS
# =============================================================================

class HUC12HotspotAnalysis:
    """Identify and prioritize impaired HUC12 subwatersheds."""

    def __init__(self, df: pd.DataFrame, config: AnalysisConfig, logger: logging.Logger):
        self.df = df
        self.config = config
        self.logger = logger

    def identify_hotspots(self) -> pd.DataFrame:
        """Rank HUC12s by impairment severity."""
        self.logger.info("=" * 60)
        self.logger.info("HUC12 HOTSPOT ANALYSIS")
        self.logger.info("=" * 60)

        df = self.df[self.df['DATA_CATEGORY'] == 'ambient'].copy()

        hotspots = []

        for huc12 in df['HUC12NAME'].dropna().unique():
            huc_data = df[df['HUC12NAME'] == huc12]

            if len(huc_data) < 10:
                continue

            result = {
                'HUC12NAME': huc12,
                'HUC12': huc_data['HUC12'].iloc[0],
                'n_samples': len(huc_data),
                'n_stations': huc_data['STATION_ID'].nunique(),
            }

            # Calculate exceedance counts for priority parameters
            score = 0

            # Fecal coliform
            fc_data = huc_data[huc_data['PARAMETER'] == 'Fecal Coliform (MF)']['VALUE'].dropna()
            if len(fc_data) > 0:
                fc_exceed = (fc_data > 400).sum()
                result['fc_n'] = len(fc_data)
                result['fc_exceed'] = fc_exceed
                result['fc_exceed_pct'] = 100 * fc_exceed / len(fc_data)
                score += fc_exceed * 2  # Weight FC higher

            # Iron
            fe_data = huc_data[huc_data['PARAMETER'] == 'Fe, Total']['VALUE'].dropna()
            if len(fe_data) > 0:
                fe_exceed = (fe_data > 1.5).sum()
                result['fe_exceed'] = fe_exceed
                score += fe_exceed

            # Specific Conductance
            sc_data = huc_data[huc_data['PARAMETER'] == 'Specific Conductance']['VALUE'].dropna()
            if len(sc_data) > 0:
                sc_exceed = (sc_data > 500).sum()
                result['sc_exceed'] = sc_exceed
                result['sc_median'] = sc_data.median()
                score += sc_exceed

            # pH
            ph_data = huc_data[huc_data['PARAMETER'] == 'pH']['VALUE'].dropna()
            if len(ph_data) > 0:
                ph_exceed = ((ph_data < 6.0) | (ph_data > 9.0)).sum()
                result['ph_exceed'] = ph_exceed
                score += ph_exceed

            result['impairment_score'] = score
            hotspots.append(result)

        hotspot_df = pd.DataFrame(hotspots)

        if len(hotspot_df) > 0:
            hotspot_df = hotspot_df.sort_values('impairment_score', ascending=False)

            self.logger.info("Top 10 Impaired HUC12s:")
            for _, row in hotspot_df.head(10).iterrows():
                self.logger.info(f"  {row['HUC12NAME']}: score={row['impairment_score']}, "
                               f"FC exceed={row.get('fc_exceed', 0)}, SC exceed={row.get('sc_exceed', 0)}")

        return hotspot_df


# =============================================================================
# VISUALIZATION (Updated)
# =============================================================================

class WaterQualityVisualizer:
    """Creates visualizations for water quality analysis."""

    def __init__(self, df: pd.DataFrame, config: AnalysisConfig, logger: logging.Logger):
        self.df = df
        self.config = config
        self.logger = logger
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette(config.color_palette)

    def create_all_visualizations(self, stats_results: dict, au_results: pd.DataFrame,
                                   hotspot_results: pd.DataFrame) -> None:
        """Generate all visualization plots."""
        self.logger.info("=" * 60)
        self.logger.info("GENERATING VISUALIZATIONS")
        self.logger.info("=" * 60)

        output_dir = self.config.output_dir / "plots"
        output_dir.mkdir(parents=True, exist_ok=True)

        self._plot_exceedance_summary(stats_results.get('exceedances', {}), output_dir)
        self._plot_fecal_coliform_detail(output_dir)
        self._plot_trend_summary(stats_results.get('trends', {}), output_dir)
        self._plot_mining_signature(output_dir)
        self._plot_huc12_hotspots(hotspot_results, output_dir)
        self._plot_temporal_coverage(output_dir)
        self._plot_sc_trend(output_dir)

        self.logger.info(f"Visualizations saved to {output_dir}")

    def _plot_exceedance_summary(self, exceedances: dict, output_dir: Path) -> None:
        """Plot WQ standard exceedance rates."""
        if not exceedances:
            return

        self.logger.info("  Creating exceedance summary plot...")

        exc_data = []
        for param, data in exceedances.items():
            exc_data.append({
                'parameter': param,
                'exceedance_pct': data['exceedance_percent'],
                'n_samples': data['n_samples'],
                'standard': f"{data['standard_value']} {data['standard_unit']}"
            })

        exc_df = pd.DataFrame(exc_data).sort_values('exceedance_pct', ascending=True)

        fig, ax = plt.subplots(figsize=(12, 8))
        colors = ['red' if x > 10 else 'orange' if x > 5 else 'green' for x in exc_df['exceedance_pct']]

        bars = ax.barh(exc_df['parameter'], exc_df['exceedance_pct'], color=colors, edgecolor='black', alpha=0.8)

        ax.set_xlabel('Exceedance Rate (%)', fontsize=12)
        ax.set_ylabel('Parameter', fontsize=12)
        ax.set_title('Water Quality Standard Exceedance Rates\n(WV 47CSR2 Standards)', fontsize=14)
        ax.axvline(10, color='red', linestyle='--', linewidth=2, label='10% Impairment Threshold')
        ax.legend()

        for bar, n, std in zip(bars, exc_df['n_samples'], exc_df['standard']):
            ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                   f'n={n:,} (std: {std})', va='center', fontsize=9)

        plt.tight_layout()
        plt.savefig(output_dir / 'exceedance_summary.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_fecal_coliform_detail(self, output_dir: Path) -> None:
        """Plot fecal coliform analysis with geomean and max standards."""
        self.logger.info("  Creating fecal coliform detail plot...")

        fc_data = self.df[
            (self.df['PARAMETER'] == 'Fecal Coliform (MF)') &
            (self.df['DATA_CATEGORY'] == 'ambient')
        ].copy()

        if len(fc_data) == 0:
            return

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Annual geometric means
        annual = fc_data.groupby('SAMPLE_YEAR').apply(
            lambda x: pd.Series({
                'geomean': np.exp(np.mean(np.log(x['VALUE'].dropna()[x['VALUE'].dropna() > 0]))),
                'n_samples': len(x),
                'pct_exceed_400': 100 * (x['VALUE'] > 400).sum() / len(x)
            })
        )

        axes[0, 0].bar(annual.index, annual['geomean'], color='steelblue', edgecolor='navy', alpha=0.8)
        axes[0, 0].axhline(200, color='red', linestyle='--', linewidth=2, label='Geomean Standard (200)')
        axes[0, 0].set_xlabel('Year')
        axes[0, 0].set_ylabel('Geometric Mean (CFU/100mL)')
        axes[0, 0].set_title('Annual Geometric Mean - Fecal Coliform')
        axes[0, 0].legend()
        axes[0, 0].tick_params(axis='x', rotation=45)

        # Exceedance rate over time
        axes[0, 1].plot(annual.index, annual['pct_exceed_400'], 'o-', color='darkred', linewidth=2)
        axes[0, 1].axhline(10, color='orange', linestyle='--', linewidth=2, label='10% Threshold')
        axes[0, 1].set_xlabel('Year')
        axes[0, 1].set_ylabel('% Samples > 400 CFU/100mL')
        axes[0, 1].set_title('Annual Exceedance Rate (>400 Standard)')
        axes[0, 1].legend()
        axes[0, 1].tick_params(axis='x', rotation=45)

        # Distribution
        fc_values = fc_data['VALUE'].dropna()
        axes[1, 0].hist(fc_values[fc_values <= 5000], bins=50, color='steelblue', edgecolor='navy', alpha=0.7)
        axes[1, 0].axvline(400, color='red', linestyle='--', linewidth=2, label='Single Sample Max (400)')
        axes[1, 0].axvline(200, color='orange', linestyle='--', linewidth=2, label='Geomean Standard (200)')
        axes[1, 0].set_xlabel('Fecal Coliform (CFU/100mL)')
        axes[1, 0].set_ylabel('Count')
        axes[1, 0].set_title('Distribution (values ≤5000)')
        axes[1, 0].legend()

        # Seasonal pattern
        season_order = ['Winter', 'Spring', 'Summer', 'Fall']
        season_geomeans = fc_data.groupby('SEASON').apply(
            lambda x: np.exp(np.mean(np.log(x['VALUE'].dropna()[x['VALUE'].dropna() > 0])))
        ).reindex(season_order)

        axes[1, 1].bar(season_geomeans.index, season_geomeans.values, color='forestgreen', edgecolor='darkgreen', alpha=0.8)
        axes[1, 1].axhline(200, color='red', linestyle='--', linewidth=2, label='Geomean Standard (200)')
        axes[1, 1].set_xlabel('Season')
        axes[1, 1].set_ylabel('Geometric Mean (CFU/100mL)')
        axes[1, 1].set_title('Seasonal Geometric Means')
        axes[1, 1].legend()

        plt.tight_layout()
        plt.savefig(output_dir / 'fecal_coliform_detail.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_trend_summary(self, trends: dict, output_dir: Path) -> None:
        """Plot trend analysis summary."""
        if not trends:
            return

        self.logger.info("  Creating trend summary plot...")

        trend_data = []
        for param, data in trends.items():
            trend_data.append({
                'parameter': param,
                'total_change_pct': data['total_change_pct'],
                'significant': data['significant'],
                'direction': data['trend_direction']
            })

        trend_df = pd.DataFrame(trend_data).dropna(subset=['total_change_pct'])
        trend_df = trend_df.sort_values('total_change_pct')

        fig, ax = plt.subplots(figsize=(12, 8))

        colors = []
        for _, row in trend_df.iterrows():
            if row['significant']:
                colors.append('darkred' if row['total_change_pct'] > 0 else 'darkgreen')
            else:
                colors.append('lightcoral' if row['total_change_pct'] > 0 else 'lightgreen')

        bars = ax.barh(trend_df['parameter'], trend_df['total_change_pct'], color=colors, edgecolor='black', alpha=0.8)

        ax.axvline(0, color='black', linewidth=1)
        ax.set_xlabel('Total Change (%)', fontsize=12)
        ax.set_ylabel('Parameter', fontsize=12)
        ax.set_title('Long-term Trends in Water Quality\n(Dark = Statistically Significant)', fontsize=14)

        plt.tight_layout()
        plt.savefig(output_dir / 'trend_summary.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_mining_signature(self, output_dir: Path) -> None:
        """Plot SC vs Sulfate to show mining signature."""
        self.logger.info("  Creating mining signature plot...")

        df = self.df[self.df['DATA_CATEGORY'] == 'ambient']

        pivot = df.pivot_table(
            index=['STATION_ID', 'SAMPLE_DATE'],
            columns='PARAMETER',
            values='VALUE',
            aggfunc='median'
        )

        if 'Specific Conductance' not in pivot.columns or 'Sulfate (SO4)' not in pivot.columns:
            return

        plot_data = pivot[['Specific Conductance', 'Sulfate (SO4)']].dropna()

        if len(plot_data) < 10:
            return

        fig, ax = plt.subplots(figsize=(10, 8))

        ax.scatter(plot_data['Sulfate (SO4)'], plot_data['Specific Conductance'],
                   alpha=0.3, s=20, c='steelblue')

        ax.axhline(500, color='red', linestyle='--', linewidth=2, label='SC Benchmark (500 µS/cm)')
        ax.axvline(250, color='orange', linestyle='--', linewidth=2, label='Sulfate MCL (250 mg/L)')

        # Add correlation
        corr = plot_data.corr(method='spearman').loc['Specific Conductance', 'Sulfate (SO4)']
        ax.text(0.05, 0.95, f'Spearman r = {corr:.3f}', transform=ax.transAxes,
                fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_xlabel('Sulfate (mg/L)', fontsize=12)
        ax.set_ylabel('Specific Conductance (µS/cm)', fontsize=12)
        ax.set_title('Mining Impact Signature: SC vs Sulfate\n(Strong positive correlation indicates surface mining)', fontsize=14)
        ax.legend(loc='lower right')

        plt.tight_layout()
        plt.savefig(output_dir / 'mining_signature.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_huc12_hotspots(self, hotspot_df: pd.DataFrame, output_dir: Path) -> None:
        """Plot HUC12 hotspot map."""
        if len(hotspot_df) == 0:
            return

        self.logger.info("  Creating HUC12 hotspot plot...")

        top_20 = hotspot_df.head(20).sort_values('impairment_score')

        fig, ax = plt.subplots(figsize=(12, 10))

        colors = plt.cm.Reds(np.linspace(0.3, 1, len(top_20)))

        bars = ax.barh(top_20['HUC12NAME'], top_20['impairment_score'], color=colors, edgecolor='black')

        ax.set_xlabel('Impairment Score', fontsize=12)
        ax.set_ylabel('HUC12 Subwatershed', fontsize=12)
        ax.set_title('Top 20 Impaired HUC12 Subwatersheds\n(Higher score = more impairment)', fontsize=14)

        plt.tight_layout()
        plt.savefig(output_dir / 'huc12_hotspots.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_temporal_coverage(self, output_dir: Path) -> None:
        """Plot temporal coverage."""
        self.logger.info("  Creating temporal coverage plot...")

        df = self.df[self.df['DATA_CATEGORY'] == 'ambient']

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        yearly = df.groupby('SAMPLE_YEAR')['SAMPLEID'].nunique()
        axes[0, 0].bar(yearly.index, yearly.values, color='steelblue', edgecolor='navy', alpha=0.8)
        axes[0, 0].set_xlabel('Year')
        axes[0, 0].set_ylabel('Number of Samples')
        axes[0, 0].set_title('Ambient Samples per Year')
        axes[0, 0].tick_params(axis='x', rotation=45)

        # Data category breakdown
        cat_yearly = self.df.groupby(['SAMPLE_YEAR', 'DATA_CATEGORY'])['SAMPLEID'].nunique().unstack(fill_value=0)
        cat_yearly.plot(kind='bar', stacked=True, ax=axes[0, 1], alpha=0.8)
        axes[0, 1].set_xlabel('Year')
        axes[0, 1].set_ylabel('Number of Samples')
        axes[0, 1].set_title('Samples by Data Category')
        axes[0, 1].legend(title='Category')
        axes[0, 1].tick_params(axis='x', rotation=45)

        # Parameters per year
        params_yearly = df.groupby('SAMPLE_YEAR')['PARAMETER'].nunique()
        axes[1, 0].plot(params_yearly.index, params_yearly.values, 'o-', color='darkorange', linewidth=2)
        axes[1, 0].set_xlabel('Year')
        axes[1, 0].set_ylabel('Number of Parameters')
        axes[1, 0].set_title('Unique Parameters Sampled per Year')
        axes[1, 0].tick_params(axis='x', rotation=45)

        # Stations per year
        stations_yearly = df.groupby('SAMPLE_YEAR')['STATION_ID'].nunique()
        axes[1, 1].plot(stations_yearly.index, stations_yearly.values, 's-', color='purple', linewidth=2)
        axes[1, 1].set_xlabel('Year')
        axes[1, 1].set_ylabel('Number of Stations')
        axes[1, 1].set_title('Unique Stations Sampled per Year')
        axes[1, 1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        plt.savefig(output_dir / 'temporal_coverage.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_sc_trend(self, output_dir: Path) -> None:
        """Plot specific conductance trend (key finding)."""
        self.logger.info("  Creating SC trend plot...")

        df = self.df[
            (self.df['PARAMETER'] == 'Specific Conductance') &
            (self.df['DATA_CATEGORY'] == 'ambient')
        ]

        if len(df) < 50:
            return

        annual = df.groupby('SAMPLE_YEAR')['VALUE'].agg(['median', 'mean', 'std', 'count'])
        annual = annual[annual['count'] >= 5]

        fig, ax = plt.subplots(figsize=(12, 6))

        ax.fill_between(annual.index,
                        annual['median'] - annual['std'],
                        annual['median'] + annual['std'],
                        alpha=0.3, color='steelblue')
        ax.plot(annual.index, annual['median'], 'o-', color='steelblue', linewidth=2, markersize=6, label='Annual Median')

        ax.axhline(500, color='red', linestyle='--', linewidth=2, label='Benchmark (500 µS/cm)')

        # Add trend line
        years = annual.index.values.astype(float)
        values = annual['median'].values
        slope, intercept, r_value, p_value, _ = stats.linregress(years, values)

        ax.plot(years, slope * years + intercept, 'g--', linewidth=2,
                label=f'Trend: {slope:.1f} µS/cm/year (p={p_value:.4f})')

        ax.set_xlabel('Year', fontsize=12)
        ax.set_ylabel('Specific Conductance (µS/cm)', fontsize=12)
        ax.set_title('Specific Conductance Trend\n(Key Mining Impact Indicator)', fontsize=14)
        ax.legend(loc='upper left')

        plt.tight_layout()
        plt.savefig(output_dir / 'sc_trend.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()


# =============================================================================
# REPORT GENERATION
# =============================================================================

class ReportGenerator:
    """Generates analysis reports and exports results."""

    def __init__(self, config: AnalysisConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger

    def generate_reports(self, df: pd.DataFrame, stats_results: dict,
                         au_results: pd.DataFrame, hotspot_results: pd.DataFrame) -> None:
        """Generate all reports and export data."""
        self.logger.info("=" * 60)
        self.logger.info("GENERATING REPORTS")
        self.logger.info("=" * 60)

        output_dir = self.config.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        # Export summary statistics
        if 'summary_stats' in stats_results:
            stats_results['summary_stats'].to_csv(output_dir / 'summary_statistics.csv', index=False)
            self.logger.info("  Saved: summary_statistics.csv")

        # Export exceedance analysis
        if 'exceedances' in stats_results:
            exc_df = pd.DataFrame(stats_results['exceedances']).T
            exc_df.to_csv(output_dir / 'exceedance_analysis.csv')
            self.logger.info("  Saved: exceedance_analysis.csv")

        # Export trend analysis
        if 'trends' in stats_results:
            trends_df = pd.DataFrame(stats_results['trends']).T
            trends_df.to_csv(output_dir / 'trend_analysis.csv')
            self.logger.info("  Saved: trend_analysis.csv")

        # Export fecal coliform analysis
        if 'fecal_coliform' in stats_results:
            fc_df = pd.DataFrame([stats_results['fecal_coliform']])
            fc_df.to_csv(output_dir / 'fecal_coliform_analysis.csv', index=False)
            self.logger.info("  Saved: fecal_coliform_analysis.csv")

        # Export AU results
        if len(au_results) > 0:
            au_results.to_csv(output_dir / 'assessment_unit_analysis.csv', index=False)
            self.logger.info("  Saved: assessment_unit_analysis.csv")

        # Export hotspot results
        if len(hotspot_results) > 0:
            hotspot_results.to_csv(output_dir / 'huc12_hotspots.csv', index=False)
            self.logger.info("  Saved: huc12_hotspots.csv")

        # Export mining signature
        if 'mining_signature' in stats_results:
            mining_df = pd.DataFrame([stats_results['mining_signature']])
            mining_df.to_csv(output_dir / 'mining_signature_analysis.csv', index=False)
            self.logger.info("  Saved: mining_signature_analysis.csv")

        # Generate text report
        self._generate_text_report(stats_results, au_results, hotspot_results, output_dir)

    def _generate_text_report(self, stats_results: dict, au_results: pd.DataFrame,
                               hotspot_results: pd.DataFrame, output_dir: Path) -> None:
        """Generate comprehensive text report."""
        lines = [
            "=" * 80,
            "WATER QUALITY ANALYSIS REPORT - VERSION 2.0",
            "Lower Guyandotte Watershed (HUC 05070102)",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "Standards: West Virginia 47CSR2",
            "=" * 80,
            "",
        ]

        # Exceedances section
        if 'exceedances' in stats_results:
            lines.extend([
                "WATER QUALITY STANDARD EXCEEDANCES (WV 47CSR2)",
                "-" * 50,
            ])
            for param, exc in stats_results['exceedances'].items():
                if exc['n_exceedances'] > 0:
                    lines.append(
                        f"{param}: {exc['n_exceedances']:,} exceedances "
                        f"({exc['exceedance_percent']:.1f}%) of {exc['n_samples']:,} samples "
                        f"[Standard: {exc['standard_value']} {exc['standard_unit']}]"
                    )
            lines.append("")

        # Fecal coliform section
        if 'fecal_coliform' in stats_results:
            fc = stats_results['fecal_coliform']
            lines.extend([
                "FECAL COLIFORM ANALYSIS",
                "-" * 50,
                f"Geometric Mean: {fc.get('geometric_mean', 'N/A'):.1f} CFU/100mL (Standard: 200)",
                f"Samples > 400: {fc.get('n_exceed_400', 0):,} ({fc.get('pct_exceed_400', 0):.1f}%)",
                f"Impaired Stations: {fc.get('n_stations_impaired', 0)} of {fc.get('n_stations_total', 0)}",
                "",
            ])

        # Trends section
        if 'trends' in stats_results:
            sig_trends = {k: v for k, v in stats_results['trends'].items() if v.get('significant')}
            if sig_trends:
                lines.extend([
                    "SIGNIFICANT TEMPORAL TRENDS",
                    "-" * 50,
                ])
                for param, trend in sig_trends.items():
                    lines.append(
                        f"{param}: {trend['trend_direction']} "
                        f"({trend['total_change_pct']:+.1f}% total change, p={trend['mann_kendall_p_value']:.4f})"
                    )
                lines.append("")

        # Hotspots section
        if len(hotspot_results) > 0:
            lines.extend([
                "TOP 10 IMPAIRED HUC12 SUBWATERSHEDS",
                "-" * 50,
            ])
            for _, row in hotspot_results.head(10).iterrows():
                lines.append(f"{row['HUC12NAME']}: score={row['impairment_score']}")
            lines.append("")

        # AU summary
        if len(au_results) > 0:
            n_impaired = au_results['impaired'].sum()
            lines.extend([
                "ASSESSMENT UNIT SUMMARY",
                "-" * 50,
                f"Total AUs Analyzed: {len(au_results)}",
                f"Impaired AUs: {n_impaired} ({100*n_impaired/len(au_results):.1f}%)",
                "",
            ])

        lines.extend([
            "=" * 80,
            "END OF REPORT",
            "=" * 80,
        ])

        with open(output_dir / 'analysis_report.txt', 'w') as f:
            f.write('\n'.join(lines))

        self.logger.info("  Saved: analysis_report.txt")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main() -> None:
    """Main execution function."""
    config = AnalysisConfig()
    logger = setup_logging(config.output_dir)

    logger.info("=" * 60)
    logger.info("WATER QUALITY DATA ANALYSIS v2.0")
    logger.info("Lower Guyandotte Watershed (HUC 05070102)")
    logger.info("Aligned with WV 47CSR2 Standards")
    logger.info("=" * 60)

    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    try:
        # Load and clean data
        loader = WaterQualityDataLoader(config, logger)
        df_raw = loader.load_data()
        df = loader.clean_data(df_raw)

        # Statistical analysis (ambient data only for TMDL)
        stats_analyzer = StatisticalAnalysis(df, config, logger)
        stats_results = stats_analyzer.run_statistical_analysis(DataCategory.AMBIENT)

        # Assessment Unit analysis
        au_analyzer = AssessmentUnitAnalysis(df, config, logger)
        au_results = au_analyzer.analyze_assessment_units()

        # HUC12 hotspot analysis
        hotspot_analyzer = HUC12HotspotAnalysis(df, config, logger)
        hotspot_results = hotspot_analyzer.identify_hotspots()

        # Visualizations
        visualizer = WaterQualityVisualizer(df, config, logger)
        visualizer.create_all_visualizations(stats_results, au_results, hotspot_results)

        # Generate reports
        reporter = ReportGenerator(config, logger)
        reporter.generate_reports(df, stats_results, au_results, hotspot_results)

        logger.info("=" * 60)
        logger.info("ANALYSIS COMPLETE")
        logger.info(f"Results saved to: {config.output_dir.absolute()}")
        logger.info("=" * 60)

    except FileNotFoundError:
        logger.error(f"Input file not found: {config.input_file}")
        raise
    except Exception as e:
        logger.exception(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
