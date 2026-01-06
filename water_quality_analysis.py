#!/usr/bin/env python3
"""
Water Quality Data Analysis Script
===================================
Comprehensive analysis of Lower Guyandotte watershed water quality monitoring data.

This script performs:
- Data loading, validation, and cleaning
- Exploratory data analysis
- Statistical summaries with proper handling of non-detects
- Temporal and spatial trend analysis
- Water quality standards exceedance analysis
- Visualization generation
- Export of processed results

Author: Generated for water quality analysis
License: MIT
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

# =============================================================================
# CONFIGURATION
# =============================================================================

@dataclass
class AnalysisConfig:
    """Configuration settings for the water quality analysis."""

    # File paths
    input_file: Path = Path("05070102_Lower_Guyandotte.csv")
    output_dir: Path = Path("wq_analysis_output")

    # Analysis settings
    non_detect_substitution: str = "half_dl"  # "half_dl", "dl", "zero", "exclude"
    min_samples_for_trend: int = 10
    significance_level: float = 0.05

    # Plotting settings
    figure_dpi: int = 150
    figure_format: str = "png"
    color_palette: str = "viridis"

    # Water quality standards (EPA recommended, adjust as needed)
    # Format: {parameter_pattern: (limit, comparison, unit, standard_name)}
    wq_standards: dict = field(default_factory=lambda: {
        "pH": (6.5, 9.0, "range", "pH units", "EPA Freshwater"),
        "DO": (5.0, "min", "mg/L", "EPA Freshwater"),
        "Temperature": (32.2, "max", "Degrees C", "EPA Warmwater"),
        "Fecal Coliform": (400, "max", "CFU/100 mL", "EPA Recreation"),
        "Fe, Total": (1.0, "max", "mg/L", "EPA Secondary"),
        "Al, Total": (0.087, "max", "mg/L", "EPA Freshwater Chronic"),
        "Cu, Total": (0.013, "max", "mg/L", "EPA Freshwater Chronic"),
        "Zn, Total": (0.120, "max", "mg/L", "EPA Freshwater Chronic"),
        "Ammonia-N": (1.9, "max", "mg/L", "EPA Freshwater Acute"),
        "Specific Conductance": (500, "max", "uS or umhos/cm", "General Guidance"),
    })

    # Parameters of primary interest for detailed analysis
    priority_parameters: list = field(default_factory=lambda: [
        "pH", "DO", "Temperature", "Specific Conductance",
        "Fecal Coliform", "Fe, Total", "Al, Total",
        "Hardness", "Alkalinity, Total", "Sulfate (SO4)"
    ])


# =============================================================================
# LOGGING SETUP
# =============================================================================

def setup_logging(output_dir: Path) -> logging.Logger:
    """Configure logging to both file and console."""
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("water_quality_analysis")
    logger.setLevel(logging.DEBUG)

    # Clear existing handlers
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

        # Define dtypes for efficiency and correctness
        dtype_spec = {
            'STATION_ID': 'str',
            'WATERBODY_NAME': 'str',
            'NHD_ANCODE': 'str',
            'ANCODE': 'str',
            'MILE_POINT': 'str',
            'DESCRIPTOR': 'str',
            'AU_ID': 'str',
            'NHDP2COMID': 'str',
            'DEP_WS': 'str',
            'HUC8': 'str',
            'HUC8NAME': 'str',
            'HUC10': 'str',
            'HUC10NAME': 'str',
            'HUC12': 'str',
            'HUC12NAME': 'str',
            'COUNTY': 'str',
            'TOPO_NAME': 'str',
            'SAMPLEID': 'str',
            'SAMPLE_TZ': 'str',
            'DUP_TYPE': 'str',
            'DUP_NAME': 'str',
            'SUBSAMPID': 'str',
            'DEPTH_DESC': 'str',
            'DIST_DESC': 'str',
            'TSECT_NAME': 'str',
            'REACH_LOC': 'str',
            'PARAMETER': 'str',
            'FRACTION': 'str',
            'SUBSTANCE': 'str',
            'NON_DETECT': 'str',
            'DEFAULT_UNITS': 'str',
            'LIM_TYPE': 'str',
            'FLAG_CODE': 'str',
            'FLAG_DESC': 'str',
            'ANL_METHOD': 'str',
        }

        # Load data
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

        # Parse dates
        df = self._parse_dates(df)

        # Clean numeric columns
        df = self._clean_numeric_columns(df)

        # Handle non-detects
        df = self._handle_non_detects(df)

        # Clean coordinates
        df = self._clean_coordinates(df)

        # Remove complete duplicates
        df = df.drop_duplicates()

        # Add derived columns
        df = self._add_derived_columns(df)

        final_rows = len(df)
        self.logger.info(f"Cleaning complete: {initial_rows:,} → {final_rows:,} records "
                        f"({initial_rows - final_rows:,} removed)")

        return df

    def _parse_dates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Parse date and time columns."""
        self.logger.debug("Parsing date/time columns")

        # Combine date and time for sample datetime
        df['SAMPLE_DATETIME'] = pd.to_datetime(
            df['SAMPLEDATE'].astype(str) + ' ' + df['SAMPLETIME'].fillna('00:00').astype(str),
            format='%Y-%m-%d %H:%M',
            errors='coerce'
        )

        # Also parse just the date for grouping
        df['SAMPLE_DATE'] = pd.to_datetime(df['SAMPLEDATE'], errors='coerce')

        # Extract useful date components
        df['SAMPLE_YEAR'] = df['SAMPLE_DATE'].dt.year
        df['SAMPLE_MONTH'] = df['SAMPLE_DATE'].dt.month
        df['SAMPLE_QUARTER'] = df['SAMPLE_DATE'].dt.quarter
        df['SAMPLE_DOY'] = df['SAMPLE_DATE'].dt.dayofyear

        # Season mapping
        season_map = {12: 'Winter', 1: 'Winter', 2: 'Winter',
                      3: 'Spring', 4: 'Spring', 5: 'Spring',
                      6: 'Summer', 7: 'Summer', 8: 'Summer',
                      9: 'Fall', 10: 'Fall', 11: 'Fall'}
        df['SEASON'] = df['SAMPLE_MONTH'].map(season_map)

        return df

    def _clean_numeric_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Clean and convert numeric columns."""
        self.logger.debug("Cleaning numeric columns")

        numeric_cols = ['VALUE', 'DL', 'QL', 'DEPTH', 'DISTANCE',
                        'MILE_POINT', 'LON_DD', 'LAT_DD']

        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        return df

    def _handle_non_detects(self, df: pd.DataFrame) -> pd.DataFrame:
        """Handle non-detect values according to configuration."""
        self.logger.debug(f"Handling non-detects using method: {self.config.non_detect_substitution}")

        # Create flag for non-detects
        df['IS_NON_DETECT'] = df['NON_DETECT'].str.upper() == 'YES'

        # Store original value
        df['VALUE_ORIGINAL'] = df['VALUE'].copy()

        # Apply substitution method
        if self.config.non_detect_substitution == "half_dl":
            # Use half the detection limit
            mask = df['IS_NON_DETECT'] & df['DL'].notna()
            df.loc[mask, 'VALUE'] = df.loc[mask, 'DL'] / 2
            # If DL not available, use half the reported value
            mask2 = df['IS_NON_DETECT'] & df['DL'].isna()
            df.loc[mask2, 'VALUE'] = df.loc[mask2, 'VALUE'] / 2

        elif self.config.non_detect_substitution == "dl":
            # Use the detection limit as-is (already in VALUE typically)
            pass

        elif self.config.non_detect_substitution == "zero":
            df.loc[df['IS_NON_DETECT'], 'VALUE'] = 0

        elif self.config.non_detect_substitution == "exclude":
            df = df[~df['IS_NON_DETECT']].copy()

        nd_count = df['IS_NON_DETECT'].sum()
        nd_pct = 100 * nd_count / len(df) if len(df) > 0 else 0
        self.logger.info(f"Non-detects: {nd_count:,} records ({nd_pct:.1f}%)")

        return df

    def _clean_coordinates(self, df: pd.DataFrame) -> pd.DataFrame:
        """Validate and clean coordinate data."""
        self.logger.debug("Cleaning coordinate data")

        # Flag invalid coordinates
        df['VALID_COORDS'] = (
            df['LON_DD'].notna() &
            df['LAT_DD'].notna() &
            (df['LON_DD'] >= -180) & (df['LON_DD'] <= 180) &
            (df['LAT_DD'] >= -90) & (df['LAT_DD'] <= 90)
        )

        invalid_coords = (~df['VALID_COORDS']).sum()
        if invalid_coords > 0:
            self.logger.warning(f"Found {invalid_coords:,} records with invalid coordinates")

        return df

    def _add_derived_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add useful derived columns."""
        self.logger.debug("Adding derived columns")

        # Create station-waterbody identifier
        df['STATION_WATERBODY'] = df['STATION_ID'].astype(str) + '_' + df['WATERBODY_NAME'].fillna('Unknown')

        # Standardize parameter names for easier grouping
        df['PARAMETER_CLEAN'] = df['PARAMETER'].str.strip().str.replace(r'\s+', ' ', regex=True)

        return df


# =============================================================================
# EXPLORATORY DATA ANALYSIS
# =============================================================================

class ExploratoryAnalysis:
    """Performs exploratory data analysis on water quality data."""

    def __init__(self, df: pd.DataFrame, config: AnalysisConfig, logger: logging.Logger):
        self.df = df
        self.config = config
        self.logger = logger
        self.results: dict[str, Any] = {}

    def run_full_eda(self) -> dict[str, Any]:
        """Run complete exploratory data analysis."""
        self.logger.info("=" * 60)
        self.logger.info("EXPLORATORY DATA ANALYSIS")
        self.logger.info("=" * 60)

        self.results['overview'] = self._data_overview()
        self.results['temporal'] = self._temporal_analysis()
        self.results['spatial'] = self._spatial_analysis()
        self.results['parameters'] = self._parameter_analysis()
        self.results['quality'] = self._data_quality_assessment()

        return self.results

    def _data_overview(self) -> dict:
        """Generate high-level data overview."""
        self.logger.info("Generating data overview...")

        overview = {
            'total_records': len(self.df),
            'unique_stations': self.df['STATION_ID'].nunique(),
            'unique_waterbodies': self.df['WATERBODY_NAME'].nunique(),
            'unique_parameters': self.df['PARAMETER'].nunique(),
            'unique_samples': self.df['SAMPLEID'].nunique(),
            'date_range': {
                'min': self.df['SAMPLE_DATE'].min(),
                'max': self.df['SAMPLE_DATE'].max(),
                'span_years': (self.df['SAMPLE_DATE'].max() - self.df['SAMPLE_DATE'].min()).days / 365.25
            },
            'huc8_count': self.df['HUC8'].nunique(),
            'huc10_count': self.df['HUC10'].nunique(),
            'huc12_count': self.df['HUC12'].nunique(),
            'counties': self.df['COUNTY'].dropna().unique().tolist(),
        }

        self.logger.info(f"  Total records: {overview['total_records']:,}")
        self.logger.info(f"  Unique stations: {overview['unique_stations']:,}")
        self.logger.info(f"  Unique parameters: {overview['unique_parameters']:,}")
        self.logger.info(f"  Date range: {overview['date_range']['min']} to {overview['date_range']['max']}")

        return overview

    def _temporal_analysis(self) -> dict:
        """Analyze temporal patterns in the data."""
        self.logger.info("Analyzing temporal patterns...")

        temporal = {}

        # Samples per year
        temporal['samples_by_year'] = (
            self.df.groupby('SAMPLE_YEAR')['SAMPLEID']
            .nunique()
            .to_dict()
        )

        # Samples per month (across all years)
        temporal['samples_by_month'] = (
            self.df.groupby('SAMPLE_MONTH')['SAMPLEID']
            .nunique()
            .to_dict()
        )

        # Samples per season
        temporal['samples_by_season'] = (
            self.df.groupby('SEASON')['SAMPLEID']
            .nunique()
            .to_dict()
        )

        # Sampling frequency statistics
        station_dates = self.df.groupby('STATION_ID')['SAMPLE_DATE'].agg(['min', 'max', 'nunique'])
        temporal['station_sampling_stats'] = {
            'mean_visits_per_station': station_dates['nunique'].mean(),
            'median_visits_per_station': station_dates['nunique'].median(),
            'max_visits_per_station': station_dates['nunique'].max(),
            'stations_single_visit': (station_dates['nunique'] == 1).sum(),
        }

        return temporal

    def _spatial_analysis(self) -> dict:
        """Analyze spatial distribution of data."""
        self.logger.info("Analyzing spatial distribution...")

        spatial = {}

        # Records per HUC12
        spatial['records_by_huc12'] = (
            self.df.groupby(['HUC12', 'HUC12NAME'])
            .agg({
                'SAMPLEID': 'nunique',
                'STATION_ID': 'nunique',
                'PARAMETER': 'nunique'
            })
            .rename(columns={
                'SAMPLEID': 'n_samples',
                'STATION_ID': 'n_stations',
                'PARAMETER': 'n_parameters'
            })
            .reset_index()
            .to_dict('records')
        )

        # Records per waterbody
        spatial['top_waterbodies'] = (
            self.df.groupby('WATERBODY_NAME')
            .agg({
                'SAMPLEID': 'nunique',
                'STATION_ID': 'nunique',
            })
            .rename(columns={'SAMPLEID': 'n_samples', 'STATION_ID': 'n_stations'})
            .sort_values('n_samples', ascending=False)
            .head(20)
            .to_dict()
        )

        # Coordinate bounding box
        valid_coords = self.df[self.df['VALID_COORDS']]
        if len(valid_coords) > 0:
            spatial['bounding_box'] = {
                'lon_min': valid_coords['LON_DD'].min(),
                'lon_max': valid_coords['LON_DD'].max(),
                'lat_min': valid_coords['LAT_DD'].min(),
                'lat_max': valid_coords['LAT_DD'].max(),
            }

        return spatial

    def _parameter_analysis(self) -> dict:
        """Analyze parameter coverage and statistics."""
        self.logger.info("Analyzing parameters...")

        params = {}

        # Parameter frequency
        param_counts = self.df['PARAMETER'].value_counts()
        params['parameter_counts'] = param_counts.to_dict()
        params['top_20_parameters'] = param_counts.head(20).to_dict()

        # Parameter statistics
        param_stats = (
            self.df.groupby('PARAMETER')
            .agg({
                'VALUE': ['count', 'mean', 'std', 'min', 'median', 'max'],
                'IS_NON_DETECT': 'sum',
                'STATION_ID': 'nunique',
                'SAMPLE_DATE': ['min', 'max']
            })
        )
        param_stats.columns = ['_'.join(col).strip() for col in param_stats.columns]
        param_stats['non_detect_pct'] = 100 * param_stats['IS_NON_DETECT_sum'] / param_stats['VALUE_count']
        params['parameter_statistics'] = param_stats.to_dict()

        # Unit consistency check
        param_units = self.df.groupby('PARAMETER')['DEFAULT_UNITS'].nunique()
        inconsistent_units = param_units[param_units > 1]
        if len(inconsistent_units) > 0:
            params['inconsistent_units'] = inconsistent_units.to_dict()
            self.logger.warning(f"Found {len(inconsistent_units)} parameters with inconsistent units")

        return params

    def _data_quality_assessment(self) -> dict:
        """Assess data quality issues."""
        self.logger.info("Assessing data quality...")

        quality = {}

        # Missing values by column
        missing = self.df.isnull().sum()
        missing_pct = 100 * missing / len(self.df)
        quality['missing_values'] = {
            col: {'count': int(missing[col]), 'percent': float(missing_pct[col])}
            for col in missing.index if missing[col] > 0
        }

        # Non-detect summary
        quality['non_detect_summary'] = {
            'total': int(self.df['IS_NON_DETECT'].sum()),
            'percent': float(100 * self.df['IS_NON_DETECT'].mean()),
            'by_parameter': (
                self.df.groupby('PARAMETER')['IS_NON_DETECT']
                .agg(['sum', 'mean'])
                .rename(columns={'sum': 'count', 'mean': 'proportion'})
                .sort_values('proportion', ascending=False)
                .head(20)
                .to_dict()
            )
        }

        # Potential outliers (values > 3 std from mean, by parameter)
        quality['potential_outliers'] = {}
        for param in self.config.priority_parameters:
            param_data = self.df[self.df['PARAMETER'] == param]['VALUE'].dropna()
            if len(param_data) > 10:
                mean, std = param_data.mean(), param_data.std()
                if std > 0:
                    outliers = param_data[np.abs(param_data - mean) > 3 * std]
                    if len(outliers) > 0:
                        quality['potential_outliers'][param] = {
                            'count': len(outliers),
                            'values': outliers.tolist()[:10]  # Limit to first 10
                        }

        # Duplicate check
        dup_cols = ['STATION_ID', 'SAMPLE_DATETIME', 'PARAMETER', 'VALUE']
        available_cols = [c for c in dup_cols if c in self.df.columns]
        duplicates = self.df.duplicated(subset=available_cols, keep=False).sum()
        quality['potential_duplicates'] = int(duplicates)

        return quality


# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================

class StatisticalAnalysis:
    """Performs statistical analysis on water quality data."""

    def __init__(self, df: pd.DataFrame, config: AnalysisConfig, logger: logging.Logger):
        self.df = df
        self.config = config
        self.logger = logger
        self.results: dict[str, Any] = {}

    def run_statistical_analysis(self) -> dict[str, Any]:
        """Run complete statistical analysis."""
        self.logger.info("=" * 60)
        self.logger.info("STATISTICAL ANALYSIS")
        self.logger.info("=" * 60)

        self.results['summary_stats'] = self._summary_statistics()
        self.results['trends'] = self._trend_analysis()
        self.results['seasonal'] = self._seasonal_analysis()
        self.results['correlations'] = self._correlation_analysis()
        self.results['exceedances'] = self._exceedance_analysis()

        return self.results

    def _summary_statistics(self) -> pd.DataFrame:
        """Calculate summary statistics for priority parameters."""
        self.logger.info("Calculating summary statistics...")

        stats_list = []

        for param in self.config.priority_parameters:
            param_data = self.df[self.df['PARAMETER'] == param].copy()

            if len(param_data) == 0:
                continue

            values = param_data['VALUE'].dropna()

            if len(values) == 0:
                continue

            # Get the most common unit
            unit = param_data['DEFAULT_UNITS'].mode()
            unit = unit.iloc[0] if len(unit) > 0 else 'Unknown'

            stats_dict = {
                'parameter': param,
                'unit': unit,
                'n': len(values),
                'n_non_detect': param_data['IS_NON_DETECT'].sum(),
                'pct_non_detect': 100 * param_data['IS_NON_DETECT'].mean(),
                'mean': values.mean(),
                'std': values.std(),
                'min': values.min(),
                'p10': values.quantile(0.10),
                'p25': values.quantile(0.25),
                'median': values.median(),
                'p75': values.quantile(0.75),
                'p90': values.quantile(0.90),
                'max': values.max(),
                'skewness': stats.skew(values) if len(values) > 2 else np.nan,
                'kurtosis': stats.kurtosis(values) if len(values) > 3 else np.nan,
                'n_stations': param_data['STATION_ID'].nunique(),
                'n_samples': param_data['SAMPLEID'].nunique(),
            }

            stats_list.append(stats_dict)

        return pd.DataFrame(stats_list)

    def _trend_analysis(self) -> dict:
        """Perform temporal trend analysis using Mann-Kendall test."""
        self.logger.info("Performing trend analysis...")

        trends = {}

        for param in self.config.priority_parameters:
            param_data = self.df[self.df['PARAMETER'] == param].copy()

            if len(param_data) < self.config.min_samples_for_trend:
                continue

            # Aggregate by year (annual median)
            annual = param_data.groupby('SAMPLE_YEAR')['VALUE'].median().dropna()

            if len(annual) < 5:
                continue

            # Simple linear regression for trend
            years = annual.index.values
            values = annual.values

            slope, intercept, r_value, p_value, std_err = stats.linregress(years, values)

            # Mann-Kendall test (simplified implementation)
            mk_stat, mk_p = self._mann_kendall_test(values)

            trends[param] = {
                'n_years': len(annual),
                'linear_slope': slope,
                'linear_r_squared': r_value ** 2,
                'linear_p_value': p_value,
                'mann_kendall_stat': mk_stat,
                'mann_kendall_p_value': mk_p,
                'trend_direction': 'increasing' if slope > 0 else 'decreasing',
                'significant': p_value < self.config.significance_level,
                'annual_values': annual.to_dict(),
            }

            if trends[param]['significant']:
                direction = trends[param]['trend_direction']
                self.logger.info(f"  {param}: Significant {direction} trend (p={p_value:.4f})")

        return trends

    def _mann_kendall_test(self, x: np.ndarray) -> tuple[float, float]:
        """Perform Mann-Kendall trend test."""
        n = len(x)

        if n < 4:
            return np.nan, np.nan

        # Calculate S statistic
        s = 0
        for i in range(n - 1):
            for j in range(i + 1, n):
                diff = x[j] - x[i]
                if diff > 0:
                    s += 1
                elif diff < 0:
                    s -= 1

        # Calculate variance
        var_s = (n * (n - 1) * (2 * n + 5)) / 18

        # Calculate Z statistic
        if s > 0:
            z = (s - 1) / np.sqrt(var_s)
        elif s < 0:
            z = (s + 1) / np.sqrt(var_s)
        else:
            z = 0

        # Calculate p-value (two-tailed)
        p_value = 2 * (1 - stats.norm.cdf(abs(z)))

        return s, p_value

    def _seasonal_analysis(self) -> dict:
        """Analyze seasonal patterns in water quality."""
        self.logger.info("Analyzing seasonal patterns...")

        seasonal = {}

        for param in self.config.priority_parameters:
            param_data = self.df[self.df['PARAMETER'] == param].copy()

            if len(param_data) < 20:
                continue

            # Calculate seasonal medians
            season_stats = param_data.groupby('SEASON')['VALUE'].agg(['median', 'mean', 'std', 'count'])

            if len(season_stats) < 2:
                continue

            # Kruskal-Wallis test for seasonal differences
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

    def _correlation_analysis(self) -> dict:
        """Analyze correlations between parameters."""
        self.logger.info("Analyzing parameter correlations...")

        # Create wide-format data for correlation
        priority_data = self.df[self.df['PARAMETER'].isin(self.config.priority_parameters)].copy()

        # Pivot to get parameters as columns
        pivot = priority_data.pivot_table(
            index=['STATION_ID', 'SAMPLE_DATE'],
            columns='PARAMETER',
            values='VALUE',
            aggfunc='median'
        )

        if pivot.shape[1] < 2:
            return {}

        # Calculate correlation matrix
        corr_matrix = pivot.corr(method='spearman')

        # Find significant correlations
        significant_corrs = []
        for i, param1 in enumerate(corr_matrix.columns):
            for j, param2 in enumerate(corr_matrix.columns):
                if i < j:  # Upper triangle only
                    r = corr_matrix.loc[param1, param2]
                    if abs(r) > 0.5 and not np.isnan(r):
                        significant_corrs.append({
                            'param1': param1,
                            'param2': param2,
                            'correlation': r,
                            'strength': 'strong' if abs(r) > 0.7 else 'moderate'
                        })

        return {
            'correlation_matrix': corr_matrix.to_dict(),
            'significant_correlations': significant_corrs,
        }

    def _exceedance_analysis(self) -> dict:
        """Analyze exceedances of water quality standards."""
        self.logger.info("Analyzing water quality standard exceedances...")

        exceedances = {}

        for param, standard in self.config.wq_standards.items():
            # Find matching parameter in data
            param_data = self.df[self.df['PARAMETER'].str.contains(param, case=False, na=False)].copy()

            if len(param_data) == 0:
                continue

            values = param_data['VALUE'].dropna()

            if len(values) == 0:
                continue

            # Handle different standard types
            if standard[1] == "range":
                # pH-style range check
                low_limit, high_limit = standard[0], standard[1] if isinstance(standard[1], (int, float)) else 9.0
                if param == "pH":
                    low_limit, high_limit = 6.5, 9.0
                exceed_low = (values < low_limit).sum()
                exceed_high = (values > high_limit).sum()
                n_exceedances = exceed_low + exceed_high
            elif standard[1] == "min":
                # Minimum value (e.g., DO)
                n_exceedances = (values < standard[0]).sum()
            else:  # "max"
                # Maximum value
                n_exceedances = (values > standard[0]).sum()

            exceedance_pct = 100 * n_exceedances / len(values)

            exceedances[param] = {
                'standard_value': standard[0],
                'standard_type': standard[1],
                'standard_unit': standard[2],
                'standard_source': standard[3],
                'n_samples': len(values),
                'n_exceedances': int(n_exceedances),
                'exceedance_percent': exceedance_pct,
                'min_value': values.min(),
                'max_value': values.max(),
                'median_value': values.median(),
            }

            if n_exceedances > 0:
                self.logger.info(f"  {param}: {n_exceedances:,} exceedances ({exceedance_pct:.1f}%)")

        return exceedances


# =============================================================================
# VISUALIZATION
# =============================================================================

class WaterQualityVisualizer:
    """Creates visualizations for water quality analysis."""

    def __init__(self, df: pd.DataFrame, config: AnalysisConfig, logger: logging.Logger):
        self.df = df
        self.config = config
        self.logger = logger

        # Set style
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette(config.color_palette)

    def create_all_visualizations(self) -> None:
        """Generate all visualization plots."""
        self.logger.info("=" * 60)
        self.logger.info("GENERATING VISUALIZATIONS")
        self.logger.info("=" * 60)

        output_dir = self.config.output_dir / "plots"
        output_dir.mkdir(parents=True, exist_ok=True)

        self._plot_temporal_coverage(output_dir)
        self._plot_parameter_distributions(output_dir)
        self._plot_time_series(output_dir)
        self._plot_seasonal_boxplots(output_dir)
        self._plot_station_heatmap(output_dir)
        self._plot_correlation_heatmap(output_dir)
        self._plot_exceedance_summary(output_dir)
        self._plot_spatial_distribution(output_dir)

        self.logger.info(f"Visualizations saved to {output_dir}")

    def _plot_temporal_coverage(self, output_dir: Path) -> None:
        """Plot temporal coverage of sampling."""
        self.logger.info("  Creating temporal coverage plot...")

        fig, axes = plt.subplots(2, 2, figsize=(14, 10))

        # Samples per year
        yearly = self.df.groupby('SAMPLE_YEAR')['SAMPLEID'].nunique()
        axes[0, 0].bar(yearly.index, yearly.values, color='steelblue', edgecolor='navy', alpha=0.8)
        axes[0, 0].set_xlabel('Year')
        axes[0, 0].set_ylabel('Number of Samples')
        axes[0, 0].set_title('Samples per Year')
        axes[0, 0].tick_params(axis='x', rotation=45)

        # Samples per month
        monthly = self.df.groupby('SAMPLE_MONTH')['SAMPLEID'].nunique()
        axes[0, 1].bar(monthly.index, monthly.values, color='forestgreen', edgecolor='darkgreen', alpha=0.8)
        axes[0, 1].set_xlabel('Month')
        axes[0, 1].set_ylabel('Number of Samples')
        axes[0, 1].set_title('Samples per Month (All Years)')
        axes[0, 1].set_xticks(range(1, 13))
        axes[0, 1].set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                                     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

        # Parameters per year
        params_yearly = self.df.groupby('SAMPLE_YEAR')['PARAMETER'].nunique()
        axes[1, 0].plot(params_yearly.index, params_yearly.values, 'o-', color='darkorange', linewidth=2)
        axes[1, 0].set_xlabel('Year')
        axes[1, 0].set_ylabel('Number of Parameters')
        axes[1, 0].set_title('Unique Parameters Sampled per Year')
        axes[1, 0].tick_params(axis='x', rotation=45)

        # Stations per year
        stations_yearly = self.df.groupby('SAMPLE_YEAR')['STATION_ID'].nunique()
        axes[1, 1].plot(stations_yearly.index, stations_yearly.values, 's-', color='purple', linewidth=2)
        axes[1, 1].set_xlabel('Year')
        axes[1, 1].set_ylabel('Number of Stations')
        axes[1, 1].set_title('Unique Stations Sampled per Year')
        axes[1, 1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        plt.savefig(output_dir / 'temporal_coverage.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_parameter_distributions(self, output_dir: Path) -> None:
        """Plot distributions of priority parameters."""
        self.logger.info("  Creating parameter distribution plots...")

        n_params = len(self.config.priority_parameters)
        n_cols = 3
        n_rows = (n_params + n_cols - 1) // n_cols

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4 * n_rows))
        axes = axes.flatten()

        for i, param in enumerate(self.config.priority_parameters):
            param_data = self.df[self.df['PARAMETER'] == param]['VALUE'].dropna()

            if len(param_data) > 0:
                # Use log scale for highly skewed data
                if param_data.max() / (param_data.median() + 0.001) > 100:
                    param_data = param_data[param_data > 0]
                    axes[i].hist(param_data, bins=50, color='steelblue', edgecolor='navy', alpha=0.7, log=True)
                    axes[i].set_ylabel('Count (log scale)')
                else:
                    axes[i].hist(param_data, bins=50, color='steelblue', edgecolor='navy', alpha=0.7)
                    axes[i].set_ylabel('Count')

                # Add median line
                med = param_data.median()
                axes[i].axvline(med, color='red', linestyle='--', linewidth=2, label=f'Median: {med:.2f}')
                axes[i].legend()

            axes[i].set_title(param)
            axes[i].set_xlabel('Value')

        # Hide unused subplots
        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        plt.tight_layout()
        plt.savefig(output_dir / 'parameter_distributions.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_time_series(self, output_dir: Path) -> None:
        """Plot time series for priority parameters."""
        self.logger.info("  Creating time series plots...")

        n_params = min(6, len(self.config.priority_parameters))
        fig, axes = plt.subplots(n_params, 1, figsize=(14, 3 * n_params), sharex=True)

        if n_params == 1:
            axes = [axes]

        for i, param in enumerate(self.config.priority_parameters[:n_params]):
            param_data = self.df[self.df['PARAMETER'] == param].copy()

            if len(param_data) > 0:
                # Plot annual medians with confidence intervals
                annual = param_data.groupby('SAMPLE_YEAR')['VALUE'].agg(['median', 'mean', 'std', 'count'])
                annual = annual[annual['count'] >= 3]  # Only years with 3+ samples

                if len(annual) > 0:
                    axes[i].fill_between(
                        annual.index,
                        annual['median'] - annual['std'],
                        annual['median'] + annual['std'],
                        alpha=0.3, color='steelblue'
                    )
                    axes[i].plot(annual.index, annual['median'], 'o-', color='steelblue',
                                linewidth=2, markersize=6, label='Annual Median')

                    # Add standard if available
                    if param in self.config.wq_standards:
                        std_val = self.config.wq_standards[param][0]
                        axes[i].axhline(std_val, color='red', linestyle='--',
                                       linewidth=2, label=f'Standard: {std_val}')

                    axes[i].legend(loc='upper right')

            axes[i].set_ylabel(param)
            axes[i].set_title(f'{param} Over Time')

        axes[-1].set_xlabel('Year')
        plt.tight_layout()
        plt.savefig(output_dir / 'time_series.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_seasonal_boxplots(self, output_dir: Path) -> None:
        """Plot seasonal boxplots for priority parameters."""
        self.logger.info("  Creating seasonal boxplots...")

        n_params = min(6, len(self.config.priority_parameters))
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()

        season_order = ['Winter', 'Spring', 'Summer', 'Fall']

        for i, param in enumerate(self.config.priority_parameters[:n_params]):
            param_data = self.df[self.df['PARAMETER'] == param].copy()

            if len(param_data) > 0:
                sns.boxplot(
                    data=param_data,
                    x='SEASON',
                    y='VALUE',
                    order=season_order,
                    ax=axes[i],
                    palette='Set2'
                )

            axes[i].set_title(param)
            axes[i].set_xlabel('Season')
            axes[i].set_ylabel('Value')

        # Hide unused subplots
        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)

        plt.tight_layout()
        plt.savefig(output_dir / 'seasonal_boxplots.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_station_heatmap(self, output_dir: Path) -> None:
        """Plot heatmap of parameter values by station."""
        self.logger.info("  Creating station heatmap...")

        # Get top stations by sample count
        top_stations = self.df['STATION_ID'].value_counts().head(30).index

        # Filter to top stations and priority parameters
        subset = self.df[
            (self.df['STATION_ID'].isin(top_stations)) &
            (self.df['PARAMETER'].isin(self.config.priority_parameters))
        ]

        if len(subset) == 0:
            return

        # Pivot to create heatmap data
        pivot = subset.pivot_table(
            index='STATION_ID',
            columns='PARAMETER',
            values='VALUE',
            aggfunc='median'
        )

        if pivot.shape[0] < 2 or pivot.shape[1] < 2:
            return

        # Normalize by column for better visualization
        pivot_norm = (pivot - pivot.min()) / (pivot.max() - pivot.min())

        fig, ax = plt.subplots(figsize=(12, max(8, len(pivot) * 0.4)))
        sns.heatmap(
            pivot_norm,
            cmap='RdYlBu_r',
            center=0.5,
            annot=False,
            ax=ax,
            cbar_kws={'label': 'Normalized Value'}
        )
        ax.set_title('Parameter Values by Station (Normalized)')
        ax.set_xlabel('Parameter')
        ax.set_ylabel('Station ID')

        plt.tight_layout()
        plt.savefig(output_dir / 'station_heatmap.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_correlation_heatmap(self, output_dir: Path) -> None:
        """Plot correlation heatmap between parameters."""
        self.logger.info("  Creating correlation heatmap...")

        # Create pivot table
        priority_data = self.df[self.df['PARAMETER'].isin(self.config.priority_parameters)]

        pivot = priority_data.pivot_table(
            index=['STATION_ID', 'SAMPLE_DATE'],
            columns='PARAMETER',
            values='VALUE',
            aggfunc='median'
        )

        if pivot.shape[1] < 2:
            return

        # Calculate correlation
        corr = pivot.corr(method='spearman')

        # Create mask for upper triangle
        mask = np.triu(np.ones_like(corr, dtype=bool))

        fig, ax = plt.subplots(figsize=(10, 8))
        sns.heatmap(
            corr,
            mask=mask,
            cmap='RdBu_r',
            center=0,
            annot=True,
            fmt='.2f',
            square=True,
            ax=ax,
            cbar_kws={'label': 'Spearman Correlation'}
        )
        ax.set_title('Parameter Correlations (Spearman)')

        plt.tight_layout()
        plt.savefig(output_dir / 'correlation_heatmap.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_exceedance_summary(self, output_dir: Path) -> None:
        """Plot summary of water quality standard exceedances."""
        self.logger.info("  Creating exceedance summary plot...")

        exceedance_data = []

        for param, standard in self.config.wq_standards.items():
            param_data = self.df[self.df['PARAMETER'].str.contains(param, case=False, na=False)]['VALUE'].dropna()

            if len(param_data) == 0:
                continue

            if standard[1] == "min":
                n_exceed = (param_data < standard[0]).sum()
            elif standard[1] == "max":
                n_exceed = (param_data > standard[0]).sum()
            else:  # range
                n_exceed = ((param_data < 6.5) | (param_data > 9.0)).sum()

            exceedance_data.append({
                'parameter': param,
                'exceedance_pct': 100 * n_exceed / len(param_data),
                'n_samples': len(param_data)
            })

        if not exceedance_data:
            return

        exc_df = pd.DataFrame(exceedance_data).sort_values('exceedance_pct', ascending=True)

        fig, ax = plt.subplots(figsize=(10, 6))
        colors = ['red' if x > 10 else 'orange' if x > 5 else 'green' for x in exc_df['exceedance_pct']]
        bars = ax.barh(exc_df['parameter'], exc_df['exceedance_pct'], color=colors, edgecolor='black', alpha=0.8)

        ax.set_xlabel('Exceedance Rate (%)')
        ax.set_ylabel('Parameter')
        ax.set_title('Water Quality Standard Exceedance Rates')
        ax.axvline(10, color='red', linestyle='--', linewidth=1, label='10% threshold')
        ax.legend()

        # Add sample count annotations
        for i, (bar, n) in enumerate(zip(bars, exc_df['n_samples'])):
            ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                   f'n={n:,}', va='center', fontsize=9)

        plt.tight_layout()
        plt.savefig(output_dir / 'exceedance_summary.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()

    def _plot_spatial_distribution(self, output_dir: Path) -> None:
        """Plot spatial distribution of sampling stations."""
        self.logger.info("  Creating spatial distribution plot...")

        # Get unique stations with coordinates
        stations = self.df[self.df['VALID_COORDS']].groupby('STATION_ID').agg({
            'LON_DD': 'first',
            'LAT_DD': 'first',
            'WATERBODY_NAME': 'first',
            'SAMPLEID': 'nunique'
        }).reset_index()

        if len(stations) == 0:
            return

        fig, ax = plt.subplots(figsize=(12, 10))

        scatter = ax.scatter(
            stations['LON_DD'],
            stations['LAT_DD'],
            c=stations['SAMPLEID'],
            s=50 + stations['SAMPLEID'] * 2,
            cmap='viridis',
            alpha=0.7,
            edgecolors='black',
            linewidths=0.5
        )

        plt.colorbar(scatter, ax=ax, label='Number of Samples')

        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_title('Sampling Station Locations\n(Size/Color = Number of Samples)')
        ax.set_aspect('equal')

        plt.tight_layout()
        plt.savefig(output_dir / 'spatial_distribution.png', dpi=self.config.figure_dpi, bbox_inches='tight')
        plt.close()


# =============================================================================
# REPORT GENERATION
# =============================================================================

class ReportGenerator:
    """Generates analysis reports and exports results."""

    def __init__(self, config: AnalysisConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger

    def generate_reports(
        self,
        df: pd.DataFrame,
        eda_results: dict,
        stats_results: dict
    ) -> None:
        """Generate all reports and export data."""
        self.logger.info("=" * 60)
        self.logger.info("GENERATING REPORTS")
        self.logger.info("=" * 60)

        output_dir = self.config.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)

        # Export summary statistics
        if 'summary_stats' in stats_results and isinstance(stats_results['summary_stats'], pd.DataFrame):
            stats_results['summary_stats'].to_csv(output_dir / 'summary_statistics.csv', index=False)
            self.logger.info("  Saved: summary_statistics.csv")

        # Export exceedance analysis
        if 'exceedances' in stats_results:
            exc_df = pd.DataFrame(stats_results['exceedances']).T
            exc_df.to_csv(output_dir / 'exceedance_analysis.csv')
            self.logger.info("  Saved: exceedance_analysis.csv")

        # Export trend analysis
        if 'trends' in stats_results:
            trends_summary = []
            for param, trend in stats_results['trends'].items():
                trends_summary.append({
                    'parameter': param,
                    'n_years': trend['n_years'],
                    'slope': trend['linear_slope'],
                    'r_squared': trend['linear_r_squared'],
                    'p_value': trend['linear_p_value'],
                    'direction': trend['trend_direction'],
                    'significant': trend['significant']
                })
            if trends_summary:
                pd.DataFrame(trends_summary).to_csv(output_dir / 'trend_analysis.csv', index=False)
                self.logger.info("  Saved: trend_analysis.csv")

        # Export cleaned data sample
        df.head(10000).to_csv(output_dir / 'cleaned_data_sample.csv', index=False)
        self.logger.info("  Saved: cleaned_data_sample.csv (first 10,000 rows)")

        # Generate text summary report
        self._generate_text_report(eda_results, stats_results, output_dir)

    def _generate_text_report(
        self,
        eda_results: dict,
        stats_results: dict,
        output_dir: Path
    ) -> None:
        """Generate a text summary report."""
        report_lines = [
            "=" * 80,
            "WATER QUALITY ANALYSIS REPORT",
            f"Lower Guyandotte Watershed (HUC 05070102)",
            f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            "=" * 80,
            "",
            "DATA OVERVIEW",
            "-" * 40,
        ]

        if 'overview' in eda_results:
            ov = eda_results['overview']
            report_lines.extend([
                f"Total Records: {ov['total_records']:,}",
                f"Unique Stations: {ov['unique_stations']:,}",
                f"Unique Waterbodies: {ov['unique_waterbodies']:,}",
                f"Unique Parameters: {ov['unique_parameters']:,}",
                f"Date Range: {ov['date_range']['min']} to {ov['date_range']['max']}",
                f"Time Span: {ov['date_range']['span_years']:.1f} years",
                "",
            ])

        # Add exceedance summary
        if 'exceedances' in stats_results:
            report_lines.extend([
                "WATER QUALITY STANDARD EXCEEDANCES",
                "-" * 40,
            ])
            for param, exc in stats_results['exceedances'].items():
                if exc['n_exceedances'] > 0:
                    report_lines.append(
                        f"{param}: {exc['n_exceedances']:,} exceedances "
                        f"({exc['exceedance_percent']:.1f}%) of {exc['n_samples']:,} samples"
                    )
            report_lines.append("")

        # Add trend summary
        if 'trends' in stats_results:
            sig_trends = {k: v for k, v in stats_results['trends'].items() if v['significant']}
            if sig_trends:
                report_lines.extend([
                    "SIGNIFICANT TEMPORAL TRENDS",
                    "-" * 40,
                ])
                for param, trend in sig_trends.items():
                    report_lines.append(
                        f"{param}: {trend['trend_direction']} "
                        f"(slope={trend['linear_slope']:.4f}, p={trend['linear_p_value']:.4f})"
                    )
                report_lines.append("")

        report_lines.extend([
            "=" * 80,
            "END OF REPORT",
            "=" * 80,
        ])

        report_path = output_dir / 'analysis_report.txt'
        with open(report_path, 'w') as f:
            f.write('\n'.join(report_lines))

        self.logger.info(f"  Saved: analysis_report.txt")


# =============================================================================
# MAIN EXECUTION
# =============================================================================

def main() -> None:
    """Main execution function."""
    # Initialize configuration
    config = AnalysisConfig()

    # Setup logging
    logger = setup_logging(config.output_dir)

    logger.info("=" * 60)
    logger.info("WATER QUALITY DATA ANALYSIS")
    logger.info("Lower Guyandotte Watershed (HUC 05070102)")
    logger.info("=" * 60)

    # Suppress warnings for cleaner output
    warnings.filterwarnings('ignore', category=FutureWarning)
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    try:
        # Load and clean data
        loader = WaterQualityDataLoader(config, logger)
        df_raw = loader.load_data()
        df = loader.clean_data(df_raw)

        # Exploratory data analysis
        eda = ExploratoryAnalysis(df, config, logger)
        eda_results = eda.run_full_eda()

        # Statistical analysis
        stats_analyzer = StatisticalAnalysis(df, config, logger)
        stats_results = stats_analyzer.run_statistical_analysis()

        # Visualizations
        visualizer = WaterQualityVisualizer(df, config, logger)
        visualizer.create_all_visualizations()

        # Generate reports
        reporter = ReportGenerator(config, logger)
        reporter.generate_reports(df, eda_results, stats_results)

        logger.info("=" * 60)
        logger.info("ANALYSIS COMPLETE")
        logger.info(f"Results saved to: {config.output_dir.absolute()}")
        logger.info("=" * 60)

    except FileNotFoundError:
        logger.error(f"Input file not found: {config.input_file}")
        logger.error("Please ensure the CSV file is in the current directory.")
        raise
    except Exception as e:
        logger.exception(f"Analysis failed: {e}")
        raise


if __name__ == "__main__":
    main()
