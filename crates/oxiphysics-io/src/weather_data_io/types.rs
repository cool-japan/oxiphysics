//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::should_implement_trait)]
#[allow(unused_imports)]
use super::functions::*;
use std::collections::HashMap;

/// A coordinate variable in a NetCDF-CF dataset.
#[derive(Debug, Clone)]
pub struct CfCoordinate {
    /// Variable name.
    pub name: String,
    /// Units string (e.g. `"degrees_north"`, `"Pa"`, `"days since 1850-01-01"`).
    pub units: String,
    /// Axis designation: `"X"`, `"Y"`, `"Z"`, or `"T"`.
    pub axis: String,
    /// Coordinate values.
    pub values: Vec<f64>,
    /// Whether the coordinate is monotonically increasing.
    pub positive: Option<String>,
}
/// Climate time series with anomaly calculation and trend detection.
#[derive(Debug, Clone)]
pub struct ClimateTimeSeries {
    /// Station or grid-point identifier.
    pub id: String,
    /// Temporal resolution.
    pub resolution: TimeResolution,
    /// Records sorted chronologically.
    pub records: Vec<ClimateRecord>,
    /// Baseline period start year for anomaly calculation.
    pub baseline_start: i32,
    /// Baseline period end year for anomaly calculation.
    pub baseline_end: i32,
}
impl ClimateTimeSeries {
    /// Create a new empty time series.
    pub fn new(id: &str, resolution: TimeResolution) -> Self {
        Self {
            id: id.to_string(),
            resolution,
            records: Vec::new(),
            baseline_start: 1961,
            baseline_end: 1990,
        }
    }
    /// Set the baseline period for anomaly calculations.
    pub fn set_baseline(&mut self, start: i32, end: i32) {
        self.baseline_start = start;
        self.baseline_end = end;
    }
    /// Add a record.
    pub fn add_record(&mut self, rec: ClimateRecord) {
        self.records.push(rec);
    }
    /// Mean temperature across all records (Celsius).
    pub fn mean_temperature(&self) -> Option<f64> {
        if self.records.is_empty() {
            return None;
        }
        let sum: f64 = self.records.iter().map(|r| r.temperature).sum();
        Some(sum / self.records.len() as f64)
    }
    /// Total precipitation across all records (mm).
    pub fn total_precipitation(&self) -> f64 {
        self.records.iter().map(|r| r.precipitation).sum()
    }
    /// Compute temperature anomalies relative to the baseline period.
    ///
    /// Returns a vector of anomaly values (one per record).
    pub fn temperature_anomalies(&self) -> Vec<f64> {
        let baseline_records: Vec<&ClimateRecord> = self
            .records
            .iter()
            .filter(|r| r.year >= self.baseline_start && r.year <= self.baseline_end)
            .collect();
        if baseline_records.is_empty() {
            return vec![0.0; self.records.len()];
        }
        let baseline_mean: f64 = baseline_records.iter().map(|r| r.temperature).sum::<f64>()
            / baseline_records.len() as f64;
        self.records
            .iter()
            .map(|r| r.temperature - baseline_mean)
            .collect()
    }
    /// Detect a linear temperature trend using least-squares regression.
    ///
    /// Returns `(slope, intercept)` where slope is degrees per record step.
    pub fn temperature_trend(&self) -> Option<(f64, f64)> {
        let n = self.records.len();
        if n < 2 {
            return None;
        }
        let nf = n as f64;
        let mut sx = 0.0_f64;
        let mut sy = 0.0_f64;
        let mut sxx = 0.0_f64;
        let mut sxy = 0.0_f64;
        for (i, r) in self.records.iter().enumerate() {
            let x = i as f64;
            let y = r.temperature;
            sx += x;
            sy += y;
            sxx += x * x;
            sxy += x * y;
        }
        let denom = nf * sxx - sx * sx;
        if denom.abs() < 1e-15 {
            return None;
        }
        let slope = (nf * sxy - sx * sy) / denom;
        let intercept = (sy - slope * sx) / nf;
        Some((slope, intercept))
    }
    /// Monthly mean temperatures (index 0 = January).
    ///
    /// Only meaningful when resolution is [`TimeResolution::Daily`] or
    /// [`TimeResolution::Monthly`].
    pub fn monthly_means(&self) -> [f64; 12] {
        let mut sums = [0.0f64; 12];
        let mut counts = [0u32; 12];
        for r in &self.records {
            if r.month >= 1 && r.month <= 12 {
                let idx = (r.month - 1) as usize;
                sums[idx] += r.temperature;
                counts[idx] += 1;
            }
        }
        let mut means = [0.0f64; 12];
        for i in 0..12 {
            if counts[i] > 0 {
                means[i] = sums[i] / counts[i] as f64;
            }
        }
        means
    }
}
/// Section type identifier for GRIB2 section parsing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GribSectionType {
    /// Section 0 — Indicator.
    Indicator,
    /// Section 1 — Identification.
    Identification,
    /// Section 2 — Local use (optional).
    LocalUse,
    /// Section 3 — Grid definition.
    GridDefinitionSection,
    /// Section 4 — Product definition.
    ProductDefinition,
    /// Section 5 — Data representation.
    DataRepresentationSection,
    /// Section 6 — Bitmap.
    Bitmap,
    /// Section 7 — Data.
    Data,
    /// Section 8 — End.
    End,
}
/// Decoded METAR report.
#[derive(Debug, Clone)]
pub struct MetarReport {
    /// ICAO airport identifier.
    pub station_id: String,
    /// Observation day of month.
    pub day: u8,
    /// Observation hour (UTC).
    pub hour: u8,
    /// Observation minute (UTC).
    pub minute: u8,
    /// Wind information.
    pub wind: MetarWind,
    /// Visibility in metres.
    pub visibility: f64,
    /// Temperature in degrees Celsius.
    pub temperature: f64,
    /// Dew-point temperature in degrees Celsius.
    pub dewpoint: f64,
    /// QNH (altimeter setting) in hPa.
    pub altimeter: f64,
    /// Raw METAR string.
    pub raw: String,
}
/// Wind rose: directional wind frequency distribution.
#[derive(Debug, Clone)]
pub struct WindRose {
    /// Number of direction bins (typically 16 or 36).
    pub direction_bins: usize,
    /// Speed bin edges in m/s (e.g. `[0, 2, 5, 10, 15, 20]`).
    pub speed_edges: Vec<f64>,
    /// Frequency matrix: `frequencies[dir_bin][speed_bin]`.
    /// Each value is a count.
    pub frequencies: Vec<Vec<u32>>,
    /// Total observation count (including calms).
    pub total_count: u32,
    /// Number of calm observations (wind speed below threshold).
    pub calm_count: u32,
    /// Calm threshold in m/s.
    pub calm_threshold: f64,
}
impl WindRose {
    /// Create a new wind rose with the specified binning.
    pub fn new(direction_bins: usize, speed_edges: Vec<f64>, calm_threshold: f64) -> Self {
        let n_speed = if speed_edges.len() > 1 {
            speed_edges.len() - 1
        } else {
            1
        };
        let frequencies = vec![vec![0u32; n_speed]; direction_bins];
        Self {
            direction_bins,
            speed_edges,
            frequencies,
            total_count: 0,
            calm_count: 0,
            calm_threshold,
        }
    }
    /// Add a wind observation.
    pub fn add_observation(&mut self, direction_deg: f64, speed: f64) {
        self.total_count += 1;
        if speed < self.calm_threshold {
            self.calm_count += 1;
            return;
        }
        let bin_width = 360.0 / self.direction_bins as f64;
        let mut dir = direction_deg % 360.0;
        if dir < 0.0 {
            dir += 360.0;
        }
        let dir_bin = ((dir + bin_width / 2.0) / bin_width) as usize % self.direction_bins;
        let speed_bin = self.speed_bin_index(speed);
        if dir_bin < self.frequencies.len() && speed_bin < self.frequencies[dir_bin].len() {
            self.frequencies[dir_bin][speed_bin] += 1;
        }
    }
    /// Find the speed bin index for a given speed.
    fn speed_bin_index(&self, speed: f64) -> usize {
        for i in 1..self.speed_edges.len() {
            if speed < self.speed_edges[i] {
                return i - 1;
            }
        }
        if self.speed_edges.len() > 1 {
            self.speed_edges.len() - 2
        } else {
            0
        }
    }
    /// Total non-calm observation count.
    pub fn non_calm_count(&self) -> u32 {
        self.total_count - self.calm_count
    }
    /// Frequency of a given direction bin as a fraction of total count.
    pub fn direction_frequency(&self, dir_bin: usize) -> f64 {
        if self.total_count == 0 || dir_bin >= self.direction_bins {
            return 0.0;
        }
        let bin_total: u32 = self.frequencies[dir_bin].iter().sum();
        bin_total as f64 / self.total_count as f64
    }
    /// Sum of all frequencies (should equal `total_count - calm_count`).
    pub fn frequency_sum(&self) -> u32 {
        self.frequencies.iter().flat_map(|row| row.iter()).sum()
    }
    /// Fit a Weibull distribution to the overall speed data.
    ///
    /// Returns `(shape_k, scale_c)` estimated from mean and variance.
    pub fn weibull_fit(&self, observations: &[(f64, f64)]) -> (f64, f64) {
        let speeds: Vec<f64> = observations
            .iter()
            .map(|&(_, s)| s)
            .filter(|&s| s >= self.calm_threshold)
            .collect();
        if speeds.is_empty() {
            return (1.0, 1.0);
        }
        let n = speeds.len() as f64;
        let mean: f64 = speeds.iter().sum::<f64>() / n;
        let variance: f64 = speeds.iter().map(|&s| (s - mean).powi(2)).sum::<f64>() / n;
        if variance < 1e-15 || mean < 1e-15 {
            return (1.0, mean.max(1.0));
        }
        let cv = variance.sqrt() / mean;
        let k = (1.0 / cv).powf(1.086);
        let k = k.clamp(0.5, 10.0);
        let c = mean / gamma_approx(1.0 + 1.0 / k);
        (k, c)
    }
    /// Mean direction (vector-average) in degrees.
    pub fn mean_direction(&self) -> f64 {
        let bin_width = 360.0 / self.direction_bins as f64;
        let mut sin_sum = 0.0_f64;
        let mut cos_sum = 0.0_f64;
        for (i, row) in self.frequencies.iter().enumerate() {
            let angle = (i as f64 * bin_width).to_radians();
            let count: u32 = row.iter().sum();
            sin_sum += angle.sin() * count as f64;
            cos_sum += angle.cos() * count as f64;
        }
        let mean_rad = sin_sum.atan2(cos_sum);
        let mut deg = mean_rad.to_degrees();
        if deg < 0.0 {
            deg += 360.0;
        }
        deg
    }
}
/// Intensity-Duration-Frequency (IDF) curve entry.
#[derive(Debug, Clone, Copy)]
pub struct IdfEntry {
    /// Duration in minutes.
    pub duration_min: f64,
    /// Return period in years.
    pub return_period_yr: f64,
    /// Intensity in mm/h.
    pub intensity: f64,
}
/// A single surface weather observation.
#[derive(Debug, Clone)]
pub struct SynopObservation {
    /// Observation time (seconds since epoch).
    pub time: i64,
    /// Temperature in Kelvin.
    pub temperature: f64,
    /// Dew-point temperature in Kelvin.
    pub dewpoint: f64,
    /// Station-level pressure in Pa.
    pub pressure: f64,
    /// Wind direction in degrees (0–360, from which the wind blows).
    pub wind_direction: f64,
    /// Wind speed in m/s.
    pub wind_speed: f64,
    /// Present weather code (WMO table 4677).
    pub weather_code: u8,
    /// Total cloud cover in oktas (0–8).
    pub cloud_cover: u8,
    /// Visibility in metres.
    pub visibility: f64,
    /// Precipitation accumulated since last observation (mm).
    pub precipitation: f64,
}
/// Temporal resolution for a climate time series.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeResolution {
    /// Daily values.
    Daily,
    /// Monthly values.
    Monthly,
    /// Annual values.
    Annual,
}
/// CF cell method describing how data was derived over a cell.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CellMethod {
    /// Point value (instantaneous).
    Point,
    /// Mean over the cell.
    Mean,
    /// Sum accumulated over the cell.
    Sum,
    /// Maximum within the cell.
    Maximum,
    /// Minimum within the cell.
    Minimum,
    /// Standard deviation over the cell.
    StandardDeviation,
}
impl CellMethod {
    /// Parse from a CF cell_methods attribute fragment.
    pub fn from_str(s: &str) -> Self {
        match s.trim() {
            "point" => Self::Point,
            "mean" => Self::Mean,
            "sum" => Self::Sum,
            "maximum" => Self::Maximum,
            "minimum" => Self::Minimum,
            "standard_deviation" => Self::StandardDeviation,
            _ => Self::Mean,
        }
    }
}
/// Parameters of a Generalised Extreme Value (GEV) distribution.
#[derive(Debug, Clone, Copy)]
pub struct GevParameters {
    /// Location parameter (mu).
    pub location: f64,
    /// Scale parameter (sigma > 0).
    pub scale: f64,
    /// Shape parameter (xi).
    pub shape: f64,
}
impl GevParameters {
    /// Create new GEV parameters.
    pub fn new(location: f64, scale: f64, shape: f64) -> Self {
        Self {
            location,
            scale: scale.max(1e-10),
            shape,
        }
    }
    /// Quantile (return level) for a given exceedance probability.
    ///
    /// `p` is the exceedance probability (e.g. 0.01 for 1% annual chance).
    pub fn quantile(&self, p: f64) -> f64 {
        let y = -(-(1.0_f64 - p).ln()).ln();
        if self.shape.abs() < 1e-10 {
            self.location + self.scale * y
        } else {
            self.location
                + self.scale / self.shape * ((-(1.0_f64 - p).ln()).powf(-self.shape) - 1.0)
        }
    }
    /// Return level for a given return period in years.
    pub fn return_level(&self, return_period: f64) -> f64 {
        let p = 1.0 / return_period;
        self.quantile(p)
    }
    /// Probability density at value `x`.
    pub fn pdf(&self, x: f64) -> f64 {
        let z = (x - self.location) / self.scale;
        if self.shape.abs() < 1e-10 {
            let e = (-z).exp();
            (1.0 / self.scale) * e * (-e).exp()
        } else {
            let t = 1.0 + self.shape * z;
            if t <= 0.0 {
                return 0.0;
            }
            let tp = t.powf(-1.0 / self.shape - 1.0);
            (1.0 / self.scale) * tp * (-t.powf(-1.0 / self.shape)).exp()
        }
    }
}
/// Writer for weather data in CSV and gridded formats.
#[derive(Debug, Clone)]
pub struct WeatherWriter {
    /// Whether to include a header row in CSV output.
    pub include_header: bool,
    /// Separator character for CSV.
    pub separator: char,
    /// Missing value representation string.
    pub missing_value: String,
    /// Number of decimal places for floating-point values.
    pub precision: usize,
}
impl WeatherWriter {
    /// Create a new writer with default settings.
    pub fn new() -> Self {
        Self::default()
    }
    /// Write station observations to CSV-format string.
    pub fn write_station_csv(&self, station: &WeatherStation) -> String {
        let mut out = String::new();
        if self.include_header {
            out.push_str(&format!(
                "# Station: {} ({}){}\n",
                station.name, station.id, ""
            ));
            out.push_str(&format!(
                "# Location: {:.4} N{}{:.4} E{}elevation {:.1} m\n",
                station.latitude,
                self.separator,
                station.longitude,
                self.separator,
                station.elevation
            ));
            let fields = [
                "time",
                "temperature_K",
                "dewpoint_K",
                "pressure_Pa",
                "wind_dir_deg",
                "wind_speed_ms",
                "weather_code",
                "cloud_oktas",
                "visibility_m",
                "precip_mm",
            ];
            out.push_str(&fields.join(&self.separator.to_string()));
            out.push('\n');
        }
        for obs in &station.observations {
            let line = format!(
                "{1}{0}{2:.prec$}{0}{3:.prec$}{0}{4:.prec$}{0}{5:.prec$}{0}{6:.prec$}{0}{7}{0}{8}{0}{9:.prec$}{0}{10:.prec$}\n",
                self.separator,
                obs.time,
                obs.temperature,
                obs.dewpoint,
                obs.pressure,
                obs.wind_direction,
                obs.wind_speed,
                obs.weather_code,
                obs.cloud_cover,
                obs.visibility,
                obs.precipitation,
                prec = self.precision,
            );
            out.push_str(&line);
        }
        out
    }
    /// Write a climate time series to CSV-format string.
    pub fn write_climate_csv(&self, ts: &ClimateTimeSeries) -> String {
        let mut out = String::new();
        if self.include_header {
            out.push_str(&format!("# Climate series: {}\n", ts.id));
            let fields = ["year", "month", "day", "temperature_C", "precip_mm"];
            out.push_str(&fields.join(&self.separator.to_string()));
            out.push('\n');
        }
        for r in &ts.records {
            let line = format!(
                "{1}{0}{2}{0}{3}{0}{4:.prec$}{0}{5:.prec$}\n",
                self.separator,
                r.year,
                r.month,
                r.day,
                r.temperature,
                r.precipitation,
                prec = self.precision,
            );
            out.push_str(&line);
        }
        out
    }
    /// Write gridded data in a simple NetCDF-like text format.
    ///
    /// Format: one header line per dimension, then row-major data values.
    pub fn write_gridded(&self, msg: &GribMessage) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "# GRIB2 discipline={} category={} parameter={}\n",
            msg.discipline, msg.parameter_category, msg.parameter_number
        ));
        out.push_str(&format!("# description: {}\n", msg.parameter_description()));
        match &msg.grid_definition {
            GridDefinition::LatLon {
                ni,
                nj,
                lat_first,
                lat_last,
                lon_first,
                lon_last,
                ..
            } => {
                out.push_str(&format!(
                    "# grid: lat_lon ni={} nj={} lat=[{},{}] lon=[{},{}]\n",
                    ni, nj, lat_first, lat_last, lon_first, lon_last
                ));
            }
            GridDefinition::LambertConformal { nx, ny, .. } => {
                out.push_str(&format!("# grid: lambert nx={} ny={}\n", nx, ny));
            }
            GridDefinition::PolarStereographic { nx, ny, .. } => {
                out.push_str(&format!("# grid: polar_stereo nx={} ny={}\n", nx, ny));
            }
        }
        for (i, &val) in msg.data.iter().enumerate() {
            if i > 0 {
                out.push(self.separator);
            }
            out.push_str(&format!("{:.prec$}", val, prec = self.precision));
        }
        out.push('\n');
        out
    }
    /// Write wind rose statistics to a CSV-format string.
    pub fn write_wind_rose_csv(&self, wr: &WindRose) -> String {
        let mut out = String::new();
        if self.include_header {
            out.push_str("# Wind Rose\n");
            out.push_str(&format!(
                "# direction_bins={} calm_count={} total_count={}\n",
                wr.direction_bins, wr.calm_count, wr.total_count
            ));
            out.push_str("direction_deg");
            for i in 0..wr.speed_edges.len().saturating_sub(1) {
                out.push(self.separator);
                out.push_str(&format!(
                    "{:.1}-{:.1}_ms",
                    wr.speed_edges[i],
                    wr.speed_edges[i + 1]
                ));
            }
            out.push('\n');
        }
        let bin_width = 360.0 / wr.direction_bins as f64;
        for (i, row) in wr.frequencies.iter().enumerate() {
            let dir = i as f64 * bin_width;
            out.push_str(&format!("{:.1}", dir));
            for &count in row {
                out.push(self.separator);
                out.push_str(&count.to_string());
            }
            out.push('\n');
        }
        out
    }
    /// Format a floating-point value, using the missing value string for NaN.
    pub fn format_value(&self, val: f64) -> String {
        if val.is_nan() {
            self.missing_value.clone()
        } else {
            format!("{:.prec$}", val, prec = self.precision)
        }
    }
}
/// A CF data variable.
#[derive(Debug, Clone)]
pub struct CfVariable {
    /// Variable name.
    pub name: String,
    /// CF standard name.
    pub standard_name: CfStandardName,
    /// Units.
    pub units: String,
    /// Cell method.
    pub cell_method: CellMethod,
    /// Dimension names (e.g. `["time", "lat", "lon"]`).
    pub dimensions: Vec<String>,
    /// Fill value for missing data.
    pub fill_value: f64,
    /// Flattened data array (row-major).
    pub data: Vec<f64>,
}
/// A climate time-series record.
#[derive(Debug, Clone)]
pub struct ClimateRecord {
    /// Year.
    pub year: i32,
    /// Month (1–12), or 0 for annual.
    pub month: u8,
    /// Day of month (1–31), or 0 for monthly/annual.
    pub day: u8,
    /// Temperature in degrees Celsius.
    pub temperature: f64,
    /// Precipitation in mm.
    pub precipitation: f64,
}
/// Precipitation analysis: IDF curves and GEV distribution fitting.
#[derive(Debug, Clone)]
pub struct PrecipitationAnalysis {
    /// Annual maximum series for each duration (key = duration in minutes).
    pub annual_maxima: HashMap<u32, Vec<f64>>,
    /// Fitted GEV parameters per duration.
    pub gev_fits: HashMap<u32, GevParameters>,
    /// Computed IDF entries.
    pub idf_entries: Vec<IdfEntry>,
}
impl PrecipitationAnalysis {
    /// Create a new empty analysis.
    pub fn new() -> Self {
        Self {
            annual_maxima: HashMap::new(),
            gev_fits: HashMap::new(),
            idf_entries: Vec::new(),
        }
    }
    /// Add an annual maximum value for a given duration.
    pub fn add_annual_max(&mut self, duration_min: u32, value: f64) {
        self.annual_maxima
            .entry(duration_min)
            .or_default()
            .push(value);
    }
    /// Fit GEV distribution to annual maxima using the method of L-moments.
    pub fn fit_gev(&mut self, duration_min: u32) -> Option<GevParameters> {
        let data = self.annual_maxima.get(&duration_min)?;
        if data.len() < 3 {
            return None;
        }
        let params = fit_gev_lmoments(data);
        self.gev_fits.insert(duration_min, params);
        Some(params)
    }
    /// Fit GEV for all durations.
    pub fn fit_all(&mut self) {
        let durations: Vec<u32> = self.annual_maxima.keys().copied().collect();
        for d in durations {
            self.fit_gev(d);
        }
    }
    /// Compute IDF entries for given return periods.
    pub fn compute_idf(&mut self, return_periods: &[f64]) {
        self.idf_entries.clear();
        let mut durations: Vec<u32> = self.gev_fits.keys().copied().collect();
        durations.sort();
        for &dur in &durations {
            if let Some(params) = self.gev_fits.get(&dur) {
                for &rp in return_periods {
                    let depth_mm = params.return_level(rp);
                    let intensity = depth_mm / (dur as f64 / 60.0);
                    self.idf_entries.push(IdfEntry {
                        duration_min: dur as f64,
                        return_period_yr: rp,
                        intensity,
                    });
                }
            }
        }
    }
    /// Check IDF curve monotonicity: for a fixed return period,
    /// intensity should decrease as duration increases.
    pub fn check_idf_monotonicity(&self, return_period: f64) -> bool {
        let mut entries: Vec<&IdfEntry> = self
            .idf_entries
            .iter()
            .filter(|e| (e.return_period_yr - return_period).abs() < 1e-6)
            .collect();
        entries.sort_by(|a, b| {
            a.duration_min
                .partial_cmp(&b.duration_min)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for w in entries.windows(2) {
            if w[1].intensity > w[0].intensity + 1e-10 {
                return false;
            }
        }
        true
    }
}
/// NetCDF-CF dataset reader (in-memory representation).
#[derive(Debug, Clone)]
pub struct NetCdfCfReader {
    /// Global attributes.
    pub global_attributes: HashMap<String, String>,
    /// Coordinate variables.
    pub coordinates: Vec<CfCoordinate>,
    /// Data variables.
    pub variables: Vec<CfVariable>,
    /// Time axis (if present).
    pub time_axis: Option<CfTimeAxis>,
}
impl NetCdfCfReader {
    /// Create an empty reader.
    pub fn new() -> Self {
        Self {
            global_attributes: HashMap::new(),
            coordinates: Vec::new(),
            variables: Vec::new(),
            time_axis: None,
        }
    }
    /// Set a global attribute.
    pub fn set_attribute(&mut self, key: &str, value: &str) {
        self.global_attributes
            .insert(key.to_string(), value.to_string());
    }
    /// Add a coordinate variable.
    pub fn add_coordinate(&mut self, coord: CfCoordinate) {
        if coord.axis == "T" {
            self.time_axis = Some(CfTimeAxis {
                units: coord.units.clone(),
                calendar: "standard".to_string(),
                values: coord.values.clone(),
            });
        }
        self.coordinates.push(coord);
    }
    /// Add a data variable.
    pub fn add_variable(&mut self, var: CfVariable) {
        self.variables.push(var);
    }
    /// Find a variable by its CF standard name.
    pub fn find_by_standard_name(&self, name: &CfStandardName) -> Option<&CfVariable> {
        self.variables.iter().find(|v| &v.standard_name == name)
    }
    /// Find a coordinate by axis designation.
    pub fn find_coordinate_by_axis(&self, axis: &str) -> Option<&CfCoordinate> {
        self.coordinates.iter().find(|c| c.axis == axis)
    }
    /// List all variable names.
    pub fn variable_names(&self) -> Vec<&str> {
        self.variables.iter().map(|v| v.name.as_str()).collect()
    }
}
/// A parsed GRIB2 section.
#[derive(Debug, Clone)]
pub struct GribSection {
    /// Section type.
    pub section_type: GribSectionType,
    /// Total byte length of the section.
    pub length: u32,
    /// Raw section payload (excluding the length and section-number bytes).
    pub payload: Vec<u8>,
}
/// A fixed weather station.
#[derive(Debug, Clone)]
pub struct WeatherStation {
    /// Station identifier.
    pub id: String,
    /// Station name.
    pub name: String,
    /// Latitude in degrees (north positive).
    pub latitude: f64,
    /// Longitude in degrees (east positive).
    pub longitude: f64,
    /// Elevation above mean sea level in metres.
    pub elevation: f64,
    /// Chronological list of SYNOP observations.
    pub observations: Vec<SynopObservation>,
}
impl WeatherStation {
    /// Create a new station.
    pub fn new(id: &str, name: &str, lat: f64, lon: f64, elev: f64) -> Self {
        Self {
            id: id.to_string(),
            name: name.to_string(),
            latitude: lat,
            longitude: lon,
            elevation: elev,
            observations: Vec::new(),
        }
    }
    /// Append an observation.
    pub fn add_observation(&mut self, obs: SynopObservation) {
        self.observations.push(obs);
    }
    /// Mean temperature (K) across all observations.
    pub fn mean_temperature(&self) -> Option<f64> {
        if self.observations.is_empty() {
            return None;
        }
        let sum: f64 = self.observations.iter().map(|o| o.temperature).sum();
        Some(sum / self.observations.len() as f64)
    }
    /// Maximum wind speed observed (m/s).
    pub fn max_wind_speed(&self) -> Option<f64> {
        self.observations
            .iter()
            .map(|o| o.wind_speed)
            .fold(None, |acc, v| {
                Some(match acc {
                    Some(mx) => {
                        if v > mx {
                            v
                        } else {
                            mx
                        }
                    }
                    None => v,
                })
            })
    }
    /// Total precipitation accumulated over all observations (mm).
    pub fn total_precipitation(&self) -> f64 {
        self.observations.iter().map(|o| o.precipitation).sum()
    }
}
/// Data representation method inside a GRIB2 data section.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataRepresentation {
    /// Simple packing (template 5.0).
    SimplePacking,
    /// Complex packing (template 5.2).
    ComplexPacking,
    /// JPEG-2000 compression (template 5.40).
    Jpeg2000,
    /// PNG compression (template 5.41).
    Png,
    /// Run-length packing for pre-defined bitmaps.
    RunLength,
}
/// CF-convention standard name for a coordinate or data variable.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum CfStandardName {
    /// Air temperature.
    AirTemperature,
    /// Air pressure.
    AirPressure,
    /// Specific humidity.
    SpecificHumidity,
    /// Precipitation flux.
    PrecipitationFlux,
    /// Wind speed.
    WindSpeed,
    /// Geopotential height.
    GeopotentialHeight,
    /// Sea surface temperature.
    SeaSurfaceTemperature,
    /// Custom standard name.
    Custom(String),
}
impl CfStandardName {
    /// Parse a CF standard_name attribute string.
    pub fn from_str(s: &str) -> Self {
        match s {
            "air_temperature" => Self::AirTemperature,
            "air_pressure" => Self::AirPressure,
            "specific_humidity" => Self::SpecificHumidity,
            "precipitation_flux" => Self::PrecipitationFlux,
            "wind_speed" => Self::WindSpeed,
            "geopotential_height" => Self::GeopotentialHeight,
            "sea_surface_temperature" => Self::SeaSurfaceTemperature,
            other => Self::Custom(other.to_string()),
        }
    }
    /// String representation matching the CF table.
    pub fn as_str(&self) -> &str {
        match self {
            Self::AirTemperature => "air_temperature",
            Self::AirPressure => "air_pressure",
            Self::SpecificHumidity => "specific_humidity",
            Self::PrecipitationFlux => "precipitation_flux",
            Self::WindSpeed => "wind_speed",
            Self::GeopotentialHeight => "geopotential_height",
            Self::SeaSurfaceTemperature => "sea_surface_temperature",
            Self::Custom(s) => s.as_str(),
        }
    }
}
/// A single GRIB2 message.
#[derive(Debug, Clone)]
pub struct GribMessage {
    /// WMO discipline code (0 = meteorological, 1 = hydrological, …).
    pub discipline: u8,
    /// Parameter category within the discipline.
    pub parameter_category: u8,
    /// Parameter number within the category.
    pub parameter_number: u8,
    /// Grid definition describing the spatial layout.
    pub grid_definition: GridDefinition,
    /// Data representation method.
    pub data_representation: DataRepresentation,
    /// Reference time (seconds since 1970-01-01T00:00 UTC).
    pub reference_time: i64,
    /// Forecast time offset in seconds.
    pub forecast_time: i64,
    /// Generating centre ID.
    pub centre_id: u16,
    /// Decoded data values (row-major order).
    pub data: Vec<f64>,
    /// Optional bitmap: `true` = valid, `false` = missing.
    pub bitmap: Option<Vec<bool>>,
}
impl GribMessage {
    /// Create a new GRIB message with the given parameters.
    pub fn new(
        discipline: u8,
        parameter_category: u8,
        parameter_number: u8,
        grid_definition: GridDefinition,
        data_representation: DataRepresentation,
    ) -> Self {
        Self {
            discipline,
            parameter_category,
            parameter_number,
            grid_definition,
            data_representation,
            reference_time: 0,
            forecast_time: 0,
            centre_id: 0,
            data: Vec::new(),
            bitmap: None,
        }
    }
    /// Total number of grid points implied by the grid definition.
    pub fn grid_point_count(&self) -> usize {
        match &self.grid_definition {
            GridDefinition::LatLon { ni, nj, .. } => (*ni as usize) * (*nj as usize),
            GridDefinition::LambertConformal { nx, ny, .. } => (*nx as usize) * (*ny as usize),
            GridDefinition::PolarStereographic { nx, ny, .. } => (*nx as usize) * (*ny as usize),
        }
    }
    /// Return non-missing data values (applying the bitmap mask).
    pub fn valid_data(&self) -> Vec<f64> {
        match &self.bitmap {
            Some(bm) => self
                .data
                .iter()
                .zip(bm.iter())
                .filter_map(|(&v, &valid)| if valid { Some(v) } else { None })
                .collect(),
            None => self.data.clone(),
        }
    }
    /// WMO parameter short description.
    pub fn parameter_description(&self) -> &'static str {
        match (
            self.discipline,
            self.parameter_category,
            self.parameter_number,
        ) {
            (0, 0, 0) => "Temperature",
            (0, 0, 6) => "Dew point temperature",
            (0, 1, 1) => "Relative humidity",
            (0, 1, 8) => "Total precipitation",
            (0, 2, 2) => "U-component of wind",
            (0, 2, 3) => "V-component of wind",
            (0, 3, 0) => "Pressure",
            (0, 3, 5) => "Geopotential height",
            (0, 6, 1) => "Total cloud cover",
            (1, 0, 0) => "Flash flood guidance",
            (2, 0, 0) => "Land cover",
            _ => "Unknown parameter",
        }
    }
}
/// GRIB2 grid definition template.
#[derive(Debug, Clone, PartialEq)]
pub enum GridDefinition {
    /// Regular latitude–longitude grid.
    LatLon {
        /// Number of points along the latitude axis.
        ni: u32,
        /// Number of points along the longitude axis.
        nj: u32,
        /// First latitude in degrees (south is negative).
        lat_first: f64,
        /// Last latitude in degrees.
        lat_last: f64,
        /// First longitude in degrees (west is negative).
        lon_first: f64,
        /// Last longitude in degrees.
        lon_last: f64,
        /// Latitude increment in degrees.
        di: f64,
        /// Longitude increment in degrees.
        dj: f64,
    },
    /// Lambert conformal conic projection.
    LambertConformal {
        /// Number of x-direction points.
        nx: u32,
        /// Number of y-direction points.
        ny: u32,
        /// Reference latitude in degrees.
        lat_ref: f64,
        /// Reference longitude in degrees.
        lon_ref: f64,
        /// First standard parallel in degrees.
        latin1: f64,
        /// Second standard parallel in degrees.
        latin2: f64,
        /// Grid spacing in the x direction (metres).
        dx: f64,
        /// Grid spacing in the y direction (metres).
        dy: f64,
    },
    /// Polar stereographic projection.
    PolarStereographic {
        /// Number of x-direction points.
        nx: u32,
        /// Number of y-direction points.
        ny: u32,
        /// Latitude of first grid point.
        lat_first: f64,
        /// Longitude of first grid point.
        lon_first: f64,
        /// Grid spacing in metres.
        dx: f64,
        /// Grid spacing in metres.
        dy: f64,
    },
}
/// Reader that parses GRIB2 byte streams into messages.
#[derive(Debug)]
pub struct GribReader {
    /// Raw byte buffer.
    pub(super) buf: Vec<u8>,
    /// Current read position.
    pub(super) pos: usize,
    /// Parsed sections for the current message.
    pub(super) sections: Vec<GribSection>,
}
impl GribReader {
    /// Create a reader from a byte slice.
    pub fn new(data: &[u8]) -> Self {
        Self {
            buf: data.to_vec(),
            pos: 0,
            sections: Vec::new(),
        }
    }
    /// Try to read the next GRIB2 message from the stream.
    pub fn next_message(&mut self) -> Option<GribMessage> {
        self.sections.clear();
        let magic = b"GRIB";
        while self.pos + 4 <= self.buf.len() {
            if &self.buf[self.pos..self.pos + 4] == magic {
                break;
            }
            self.pos += 1;
        }
        if self.pos + 16 > self.buf.len() {
            return None;
        }
        let discipline = self.buf[self.pos + 6];
        let edition = self.buf[self.pos + 7];
        if edition != 2 {
            return None;
        }
        let total_len = self.read_u64_at(self.pos + 8) as usize;
        if self.pos + total_len > self.buf.len() {
            return None;
        }
        self.sections.push(GribSection {
            section_type: GribSectionType::Indicator,
            length: 16,
            payload: self.buf[self.pos..self.pos + 16].to_vec(),
        });
        let msg_end = self.pos + total_len;
        let mut cur = self.pos + 16;
        while cur + 4 <= msg_end {
            if &self.buf[cur..cur + 4] == b"7777" {
                self.sections.push(GribSection {
                    section_type: GribSectionType::End,
                    length: 4,
                    payload: b"7777".to_vec(),
                });
                let _ = cur + 4;
                break;
            }
            let sec_len = self.read_u32_at(cur) as usize;
            if sec_len < 5 || cur + sec_len > msg_end {
                break;
            }
            let sec_num = self.buf[cur + 4];
            let sec_type = match sec_num {
                1 => GribSectionType::Identification,
                2 => GribSectionType::LocalUse,
                3 => GribSectionType::GridDefinitionSection,
                4 => GribSectionType::ProductDefinition,
                5 => GribSectionType::DataRepresentationSection,
                6 => GribSectionType::Bitmap,
                7 => GribSectionType::Data,
                _ => GribSectionType::LocalUse,
            };
            self.sections.push(GribSection {
                section_type: sec_type,
                length: sec_len as u32,
                payload: self.buf[cur + 5..cur + sec_len].to_vec(),
            });
            cur += sec_len;
        }
        self.pos = msg_end;
        let msg = self.assemble_message(discipline);
        Some(msg)
    }
    /// Assemble a message from parsed sections.
    fn assemble_message(&self, discipline: u8) -> GribMessage {
        let mut param_cat: u8 = 0;
        let mut param_num: u8 = 0;
        let mut centre: u16 = 0;
        let mut ref_time: i64 = 0;
        let mut forecast: i64 = 0;
        let mut grid_def = GridDefinition::LatLon {
            ni: 1,
            nj: 1,
            lat_first: 0.0,
            lat_last: 0.0,
            lon_first: 0.0,
            lon_last: 0.0,
            di: 1.0,
            dj: 1.0,
        };
        let mut data_rep = DataRepresentation::SimplePacking;
        let mut data_vals: Vec<f64> = Vec::new();
        let mut bitmap_vals: Option<Vec<bool>> = None;
        for sec in &self.sections {
            match sec.section_type {
                GribSectionType::Identification if sec.payload.len() >= 12 => {
                    centre = u16::from_be_bytes([sec.payload[0], sec.payload[1]]);
                    let year = u16::from_be_bytes([sec.payload[7], sec.payload[8]]) as i64;
                    let month = sec.payload[9] as i64;
                    let day = sec.payload[10] as i64;
                    let hour = sec.payload[11] as i64;
                    ref_time = ((year - 1970) * 365 * 86400)
                        + ((month - 1) * 30 * 86400)
                        + ((day - 1) * 86400)
                        + (hour * 3600);
                }
                GribSectionType::GridDefinitionSection if sec.payload.len() >= 25 => {
                    let template = u16::from_be_bytes([sec.payload[7], sec.payload[8]]);
                    if template == 0 {
                        let ni = self.read_payload_u32(&sec.payload, 9);
                        let nj = self.read_payload_u32(&sec.payload, 13);
                        grid_def = GridDefinition::LatLon {
                            ni,
                            nj,
                            lat_first: 0.0,
                            lat_last: 0.0,
                            lon_first: 0.0,
                            lon_last: 0.0,
                            di: 1.0,
                            dj: 1.0,
                        };
                    }
                }
                GribSectionType::ProductDefinition => {
                    if sec.payload.len() >= 4 {
                        param_cat = sec.payload[2];
                        param_num = sec.payload[3];
                    }
                    if sec.payload.len() >= 19 {
                        forecast = sec.payload[18] as i64 * 3600;
                    }
                }
                GribSectionType::DataRepresentationSection if sec.payload.len() >= 4 => {
                    let template = u16::from_be_bytes([sec.payload[2], sec.payload[3]]);
                    data_rep = match template {
                        0 => DataRepresentation::SimplePacking,
                        2 => DataRepresentation::ComplexPacking,
                        40 => DataRepresentation::Jpeg2000,
                        41 => DataRepresentation::Png,
                        _ => DataRepresentation::SimplePacking,
                    };
                }
                GribSectionType::Bitmap if !sec.payload.is_empty() => {
                    let indicator = sec.payload[0];
                    if indicator == 0 && sec.payload.len() > 1 {
                        let mut bm = Vec::new();
                        for &byte in &sec.payload[1..] {
                            for bit in (0..8).rev() {
                                bm.push((byte >> bit) & 1 != 0);
                            }
                        }
                        bitmap_vals = Some(bm);
                    }
                }
                GribSectionType::Data => {
                    let chunk_size = 4;
                    let count = sec.payload.len() / chunk_size;
                    data_vals.reserve(count);
                    for i in 0..count {
                        let off = i * chunk_size;
                        if off + chunk_size <= sec.payload.len() {
                            let bytes = [
                                sec.payload[off],
                                sec.payload[off + 1],
                                sec.payload[off + 2],
                                sec.payload[off + 3],
                            ];
                            data_vals.push(f32::from_be_bytes(bytes) as f64);
                        }
                    }
                }
                _ => {}
            }
        }
        GribMessage {
            discipline,
            parameter_category: param_cat,
            parameter_number: param_num,
            grid_definition: grid_def,
            data_representation: data_rep,
            reference_time: ref_time,
            forecast_time: forecast,
            centre_id: centre,
            data: data_vals,
            bitmap: bitmap_vals,
        }
    }
    /// Read a big-endian u64 at an absolute position.
    fn read_u64_at(&self, pos: usize) -> u64 {
        let mut bytes = [0u8; 8];
        bytes.copy_from_slice(&self.buf[pos..pos + 8]);
        u64::from_be_bytes(bytes)
    }
    /// Read a big-endian u32 at an absolute position.
    fn read_u32_at(&self, pos: usize) -> u32 {
        let mut bytes = [0u8; 4];
        bytes.copy_from_slice(&self.buf[pos..pos + 4]);
        u32::from_be_bytes(bytes)
    }
    /// Read a big-endian u32 from a payload slice at a relative offset.
    fn read_payload_u32(&self, payload: &[u8], off: usize) -> u32 {
        if off + 4 <= payload.len() {
            u32::from_be_bytes([
                payload[off],
                payload[off + 1],
                payload[off + 2],
                payload[off + 3],
            ])
        } else {
            0
        }
    }
    /// Return the parsed sections from the most recent `next_message` call.
    pub fn sections(&self) -> &[GribSection] {
        &self.sections
    }
}
/// Wind information decoded from a METAR string.
#[derive(Debug, Clone, PartialEq)]
pub struct MetarWind {
    /// Direction in degrees true north (0 = calm/VRB).
    pub direction: u16,
    /// Sustained speed in knots.
    pub speed: u16,
    /// Gust speed in knots (0 if no gust).
    pub gust: u16,
    /// Whether direction is variable.
    pub variable: bool,
}
/// CF time axis utilities.
#[derive(Debug, Clone)]
pub struct CfTimeAxis {
    /// Reference date as "days since YYYY-MM-DD".
    pub units: String,
    /// Calendar type.
    pub calendar: String,
    /// Time coordinate values in the native unit.
    pub values: Vec<f64>,
}
impl CfTimeAxis {
    /// Create a new CF time axis.
    pub fn new(units: &str, calendar: &str, values: Vec<f64>) -> Self {
        Self {
            units: units.to_string(),
            calendar: calendar.to_string(),
            values,
        }
    }
    /// Parse the base date from the units string.
    ///
    /// Expects format `"`unit` since YYYY-MM-DD"`.
    pub fn base_date(&self) -> Option<(i32, u32, u32)> {
        let since_pos = self.units.find("since ")?;
        let date_str = self.units[since_pos + 6..].trim();
        let parts: Vec<&str> = date_str.split('-').collect();
        if parts.len() >= 3 {
            let year: i32 = parts[0].parse().ok()?;
            let month: u32 = parts[1].parse().ok()?;
            let day: u32 = parts[2].split_whitespace().next()?.parse().ok()?;
            Some((year, month, day))
        } else {
            None
        }
    }
    /// Return the time unit multiplier (seconds per unit).
    pub fn seconds_per_unit(&self) -> f64 {
        if self.units.starts_with("seconds") {
            1.0
        } else if self.units.starts_with("minutes") {
            60.0
        } else if self.units.starts_with("hours") {
            3600.0
        } else if self.units.starts_with("days") {
            86400.0
        } else {
            1.0
        }
    }
    /// Convert time values to seconds since the base date.
    pub fn to_seconds(&self) -> Vec<f64> {
        let mult = self.seconds_per_unit();
        self.values.iter().map(|&v| v * mult).collect()
    }
    /// Number of time steps.
    pub fn len(&self) -> usize {
        self.values.len()
    }
    /// Whether the axis is empty.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}
