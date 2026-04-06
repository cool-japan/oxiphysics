//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::should_implement_trait)]
use super::functions::detect_delimiter;
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;

/// A single record (row) in a CSV file.
#[allow(dead_code)]
pub struct CsvRecord {
    /// The fields (columns) in this record.
    pub fields: Vec<String>,
}
/// An in-memory representation of a CSV file with headers.
#[allow(dead_code)]
pub struct CsvFile {
    /// Column header names.
    pub headers: Vec<String>,
    /// All data records.
    pub records: Vec<CsvRecord>,
}
impl CsvFile {
    /// Create a new `CsvFile` with the given headers and no records.
    pub fn new(headers: Vec<String>) -> Self {
        CsvFile {
            headers,
            records: Vec::new(),
        }
    }
    /// Append a record from string fields.
    pub fn add_record(&mut self, fields: Vec<String>) {
        self.records.push(CsvRecord { fields });
    }
    /// Append a record from `f64` values (formatted with full precision).
    pub fn add_record_f64(&mut self, values: &[f64]) {
        let fields = values.iter().map(|v| format!("{}", v)).collect();
        self.records.push(CsvRecord { fields });
    }
    /// Return the number of data records (excluding header).
    pub fn record_count(&self) -> usize {
        self.records.len()
    }
    /// Return the number of columns (header count).
    pub fn column_count(&self) -> usize {
        self.headers.len()
    }
    /// Extract all values from a column by index as `f64`.
    pub fn get_column_f64(&self, col_idx: usize) -> Result<Vec<f64>, String> {
        if col_idx >= self.headers.len() {
            return Err(format!("Column index {} out of range", col_idx));
        }
        let mut out = Vec::with_capacity(self.records.len());
        for (row, rec) in self.records.iter().enumerate() {
            let s = rec
                .fields
                .get(col_idx)
                .ok_or_else(|| format!("Row {} has no field at column {}", row, col_idx))?;
            let v: f64 = s
                .trim()
                .parse()
                .map_err(|e| format!("Row {}, col {}: parse error: {}", row, col_idx, e))?;
            out.push(v);
        }
        Ok(out)
    }
    /// Return the column index of a header name, or `None` if not found.
    pub fn get_column_by_name(&self, name: &str) -> Option<usize> {
        self.headers.iter().position(|h| h == name)
    }
    /// Serialize the CSV file to a `String` (header + records, comma-separated).
    #[allow(clippy::inherent_to_string)]
    pub fn to_string(&self) -> String {
        self.to_string_with_delimiter(',')
    }
    /// Serialize the CSV file using a custom delimiter character.
    #[allow(dead_code)]
    pub fn to_string_with_delimiter(&self, delim: char) -> String {
        let mut out = String::new();
        let d = delim.to_string();
        out.push_str(&self.headers.join(&d));
        out.push('\n');
        for rec in &self.records {
            out.push_str(&rec.fields.join(&d));
            out.push('\n');
        }
        out
    }
    /// Parse a CSV string (first line = header, comma delimiter).
    pub fn from_str(s: &str) -> Result<Self, String> {
        Self::from_str_with_delimiter(s, ',')
    }
    /// Parse a CSV string with a specified delimiter.
    #[allow(dead_code)]
    pub fn from_str_with_delimiter(s: &str, delim: char) -> Result<Self, String> {
        let mut lines = s.lines();
        let header_line = lines.next().ok_or("Empty CSV input")?;
        let headers: Vec<String> = header_line
            .split(delim)
            .map(|f| f.trim().to_string())
            .collect();
        if headers.is_empty() || headers.iter().all(|h| h.is_empty()) {
            return Err("No headers found".to_string());
        }
        let mut records = Vec::new();
        for line in lines {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<String> = line.split(delim).map(|f| f.trim().to_string()).collect();
            records.push(CsvRecord { fields });
        }
        Ok(CsvFile { headers, records })
    }
    /// Return a new `CsvFile` containing only rows where `pred(value)` is true
    /// for the value in `col_idx`.  Rows where parsing fails are excluded.
    pub fn filter_rows(&self, col_idx: usize, pred: impl Fn(f64) -> bool) -> CsvFile {
        let mut out = CsvFile::new(self.headers.clone());
        for rec in &self.records {
            if let Some(s) = rec.fields.get(col_idx)
                && let Ok(v) = s.trim().parse::<f64>()
                && pred(v)
            {
                out.records.push(CsvRecord {
                    fields: rec.fields.clone(),
                });
            }
        }
        out
    }
    /// Infer the type of a column (Integer, Float, or Text).
    #[allow(dead_code)]
    pub fn infer_column_type(&self, col_idx: usize) -> ColumnType {
        if col_idx >= self.headers.len() {
            return ColumnType::Text;
        }
        let mut all_int = true;
        let mut all_float = true;
        let mut any_value = false;
        for rec in &self.records {
            if let Some(s) = rec.fields.get(col_idx) {
                let s = s.trim();
                if s.is_empty() {
                    continue;
                }
                any_value = true;
                if s.parse::<i64>().is_err() {
                    all_int = false;
                }
                if s.parse::<f64>().is_err() {
                    all_float = false;
                }
            }
        }
        if !any_value {
            return ColumnType::Text;
        }
        if all_int {
            ColumnType::Integer
        } else if all_float {
            ColumnType::Float
        } else {
            ColumnType::Text
        }
    }
    /// Return a new `CsvFile` containing only the specified columns (by index).
    #[allow(dead_code)]
    pub fn select_columns(&self, col_indices: &[usize]) -> CsvFile {
        let headers: Vec<String> = col_indices
            .iter()
            .filter_map(|&i| self.headers.get(i).cloned())
            .collect();
        let mut out = CsvFile::new(headers);
        for rec in &self.records {
            let fields: Vec<String> = col_indices
                .iter()
                .map(|&i| rec.fields.get(i).cloned().unwrap_or_default())
                .collect();
            out.records.push(CsvRecord { fields });
        }
        out
    }
    /// Return a new `CsvFile` containing only columns whose headers match
    /// the given names (preserving the order of `names`).
    #[allow(dead_code)]
    pub fn select_columns_by_name(&self, names: &[&str]) -> CsvFile {
        let indices: Vec<usize> = names
            .iter()
            .filter_map(|n| self.get_column_by_name(n))
            .collect();
        self.select_columns(&indices)
    }
    /// Normalize headers: lowercase, replace spaces/special chars with underscores,
    /// strip leading/trailing whitespace.
    #[allow(dead_code)]
    pub fn normalize_headers(&mut self) {
        for h in &mut self.headers {
            let normalized: String = h
                .trim()
                .to_lowercase()
                .chars()
                .map(|c| {
                    if c.is_alphanumeric() || c == '_' {
                        c
                    } else {
                        '_'
                    }
                })
                .collect();
            *h = normalized;
        }
    }
    /// Compute statistics (min, max, mean, sum, count) for a numeric column.
    /// Returns `None` if no numeric values are found.
    #[allow(dead_code)]
    pub fn column_stats(&self, col_idx: usize) -> Option<ColumnStats> {
        let values = self.get_column_f64(col_idx).ok()?;
        if values.is_empty() {
            return None;
        }
        let mut min = f64::INFINITY;
        let mut max = f64::NEG_INFINITY;
        let mut sum = 0.0;
        for &v in &values {
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
            sum += v;
        }
        let count = values.len();
        Some(ColumnStats {
            min,
            max,
            mean: sum / count as f64,
            count,
            sum,
        })
    }
    /// Compute statistics for all columns that are numeric.
    /// Returns a vector of `(column_name, ColumnStats)`.
    #[allow(dead_code)]
    pub fn all_column_stats(&self) -> Vec<(String, ColumnStats)> {
        let mut result = Vec::new();
        for i in 0..self.headers.len() {
            if let Some(stats) = self.column_stats(i) {
                result.push((self.headers[i].clone(), stats));
            }
        }
        result
    }
    /// Extract all values from a column as strings.
    #[allow(dead_code)]
    pub fn get_column_strings(&self, col_idx: usize) -> Result<Vec<String>, String> {
        if col_idx >= self.headers.len() {
            return Err(format!("Column index {} out of range", col_idx));
        }
        let mut out = Vec::with_capacity(self.records.len());
        for (row, rec) in self.records.iter().enumerate() {
            let s = rec
                .fields
                .get(col_idx)
                .ok_or_else(|| format!("Row {} has no field at column {}", row, col_idx))?;
            out.push(s.trim().to_string());
        }
        Ok(out)
    }
    /// Extract all values from a column as `i64`.
    #[allow(dead_code)]
    pub fn get_column_i64(&self, col_idx: usize) -> Result<Vec<i64>, String> {
        if col_idx >= self.headers.len() {
            return Err(format!("Column index {} out of range", col_idx));
        }
        let mut out = Vec::with_capacity(self.records.len());
        for (row, rec) in self.records.iter().enumerate() {
            let s = rec
                .fields
                .get(col_idx)
                .ok_or_else(|| format!("Row {} has no field at column {}", row, col_idx))?;
            let v: i64 = s
                .trim()
                .parse()
                .map_err(|e| format!("Row {}, col {}: parse error: {}", row, col_idx, e))?;
            out.push(v);
        }
        Ok(out)
    }
    /// Sort rows by a column (ascending, numeric).
    #[allow(dead_code)]
    pub fn sort_by_column(&mut self, col_idx: usize) {
        self.records.sort_by(|a, b| {
            let va = a
                .fields
                .get(col_idx)
                .and_then(|s| s.trim().parse::<f64>().ok())
                .unwrap_or(f64::NAN);
            let vb = b
                .fields
                .get(col_idx)
                .and_then(|s| s.trim().parse::<f64>().ok())
                .unwrap_or(f64::NAN);
            va.partial_cmp(&vb).unwrap_or(std::cmp::Ordering::Equal)
        });
    }
}
/// A streaming CSV writer that builds output row-by-row.
#[allow(dead_code)]
pub struct CsvWriter {
    pub(super) headers: Vec<String>,
    pub(super) delimiter: char,
    pub(super) lines: Vec<String>,
}
impl CsvWriter {
    /// Create a writer with the given column headers and delimiter.
    pub fn new(headers: Vec<String>, delimiter: char) -> Self {
        let header_line = headers.join(&delimiter.to_string());
        Self {
            headers,
            delimiter,
            lines: vec![header_line],
        }
    }
    /// Append a row of string values.  Extra values are truncated; missing
    /// values are filled with empty strings.
    pub fn write_row(&mut self, values: &[&str]) {
        let n = self.headers.len();
        let row: Vec<&str> = (0..n)
            .map(|i| values.get(i).copied().unwrap_or(""))
            .collect();
        self.lines.push(row.join(&self.delimiter.to_string()));
    }
    /// Append a row of f64 values formatted to 6 significant figures.
    pub fn write_row_f64(&mut self, values: &[f64]) {
        let strs: Vec<String> = values.iter().map(|v| format!("{v:.6}")).collect();
        let refs: Vec<&str> = strs.iter().map(String::as_str).collect();
        self.write_row(&refs);
    }
    /// Return the accumulated CSV as a string.
    pub fn finish(self) -> String {
        self.lines.join("\n")
    }
    /// Number of data rows written (not counting the header).
    pub fn row_count(&self) -> usize {
        self.lines.len().saturating_sub(1)
    }
}
/// Aggregation operations for a numeric CSV column.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AggOp {
    /// Sum of all values.
    Sum,
    /// Arithmetic mean.
    Mean,
    /// Minimum value.
    Min,
    /// Maximum value.
    Max,
    /// Standard deviation (population).
    Std,
    /// Count of rows.
    Count,
}
/// A time-series CSV file: a [`CsvFile`] with an explicit time column.
///
/// The time column is identified by name (default `"time"`).
#[allow(dead_code)]
pub struct TimeSeriesCsv {
    /// Underlying CSV data.
    pub csv: CsvFile,
    /// Name of the time column.
    pub time_column: String,
}
impl TimeSeriesCsv {
    /// Wrap an existing [`CsvFile`], treating column `time_column` as time.
    pub fn new(csv: CsvFile, time_column: &str) -> Self {
        Self {
            csv,
            time_column: time_column.to_owned(),
        }
    }
    /// Parse a CSV string and treat `time_column` as the time axis.
    pub fn from_str(s: &str, time_column: &str) -> Result<Self, String> {
        let csv = CsvFile::from_str(s)?;
        Ok(Self::new(csv, time_column))
    }
    /// Extract the time column as `Vec`f64`.  Returns `None` if the column
    /// doesn't exist or contains non-numeric values.
    pub fn times(&self) -> Option<Vec<f64>> {
        let idx = self
            .csv
            .headers
            .iter()
            .position(|h| h == &self.time_column)?;
        self.csv.get_column_f64(idx).ok()
    }
    /// Extract a named data column as `Vec`f64`.
    pub fn column_f64(&self, name: &str) -> Option<Vec<f64>> {
        let idx = self.csv.headers.iter().position(|h| h == name)?;
        self.csv.get_column_f64(idx).ok()
    }
    /// Number of time steps (rows) in the series.
    pub fn n_steps(&self) -> usize {
        self.csv.records.len()
    }
    /// Duration = max(time) - min(time).  Returns `0.0` on empty or single-row.
    pub fn duration(&self) -> f64 {
        let ts = self.times().unwrap_or_default();
        if ts.len() < 2 {
            return 0.0;
        }
        let min = ts.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = ts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        max - min
    }
}
/// A schema definition for a CSV file: each column has a name and expected type.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CsvSchema {
    /// Ordered column definitions.
    pub columns: Vec<(String, ColumnType)>,
}
impl CsvSchema {
    /// Create a schema from paired `(name, type)` entries.
    pub fn new(columns: Vec<(String, ColumnType)>) -> Self {
        Self { columns }
    }
    /// Number of columns in this schema.
    pub fn len(&self) -> usize {
        self.columns.len()
    }
    /// Returns `true` if the schema has no columns.
    pub fn is_empty(&self) -> bool {
        self.columns.is_empty()
    }
    /// Validate that a [`CsvFile`] matches this schema.
    ///
    /// Checks column count and, for numeric types, that every value in the
    /// column parses successfully.  Returns a list of human-readable error
    /// messages (empty if valid).
    pub fn validate(&self, csv: &CsvFile) -> Vec<String> {
        let mut errors = Vec::new();
        if csv.headers.len() != self.columns.len() {
            errors.push(format!(
                "column count mismatch: schema has {} columns, file has {}",
                self.columns.len(),
                csv.headers.len()
            ));
            return errors;
        }
        for (col_idx, (name, expected_type)) in self.columns.iter().enumerate() {
            if csv.headers[col_idx] != *name {
                errors.push(format!(
                    "column {} name mismatch: expected '{}', got '{}'",
                    col_idx, name, csv.headers[col_idx]
                ));
            }
            for (row_idx, record) in csv.records.iter().enumerate() {
                if col_idx >= record.fields.len() {
                    errors.push(format!("row {} column {} missing", row_idx, col_idx));
                    continue;
                }
                let v = &record.fields[col_idx];
                match expected_type {
                    ColumnType::Integer => {
                        if v.parse::<i64>().is_err() {
                            errors.push(format!(
                                "row {} column '{}': expected Integer, got '{}'",
                                row_idx, name, v
                            ));
                        }
                    }
                    ColumnType::Float => {
                        if v.parse::<f64>().is_err() {
                            errors.push(format!(
                                "row {} column '{}': expected Float, got '{}'",
                                row_idx, name, v
                            ));
                        }
                    }
                    ColumnType::Text => {}
                }
            }
        }
        errors
    }
}
/// A lazy CSV line iterator that yields one record at a time without
/// loading the entire file into memory.
///
/// The iterator owns the input string and yields `Vec`String` for each
/// non-header row.
pub struct LazyCsvIter<'a> {
    pub(super) lines: std::str::Lines<'a>,
    pub(super) delimiter: char,
    /// The parsed headers from the first line.
    pub headers: Vec<String>,
}
impl<'a> LazyCsvIter<'a> {
    /// Create a lazy iterator over `input`.  The first line is consumed
    /// immediately as the header row.
    pub fn new(input: &'a str, delimiter: char) -> Self {
        let mut lines = input.lines();
        let headers = lines
            .next()
            .map(|h| {
                h.split(delimiter)
                    .map(str::trim)
                    .map(String::from)
                    .collect()
            })
            .unwrap_or_default();
        Self {
            lines,
            delimiter,
            headers,
        }
    }
}
/// Summary of validation errors found in a CSV file.
#[allow(dead_code)]
#[derive(Debug, Default)]
pub struct CsvValidationReport {
    /// All error messages.
    pub errors: Vec<String>,
}
impl CsvValidationReport {
    /// Returns `true` if no errors were found.
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }
    /// Number of errors.
    pub fn error_count(&self) -> usize {
        self.errors.len()
    }
}
/// Inferred column types.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub enum ColumnType {
    /// All values parse as integers.
    Integer,
    /// All values parse as floating-point numbers.
    Float,
    /// Fallback: treated as text.
    Text,
}
/// Typed column data in a [`CsvDataFrame`].
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum CsvColumnData {
    /// Column of 64-bit integers.
    Integer(Vec<i64>),
    /// Column of 64-bit floats.
    Float(Vec<f64>),
    /// Column of strings.
    Text(Vec<String>),
}
#[allow(dead_code)]
impl CsvColumnData {
    /// Return the length of this column.
    pub fn len(&self) -> usize {
        match self {
            CsvColumnData::Integer(v) => v.len(),
            CsvColumnData::Float(v) => v.len(),
            CsvColumnData::Text(v) => v.len(),
        }
    }
    /// Returns `true` if the column has no elements.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
    /// Inferred type as [`ColumnType`].
    pub fn column_type(&self) -> ColumnType {
        match self {
            CsvColumnData::Integer(_) => ColumnType::Integer,
            CsvColumnData::Float(_) => ColumnType::Float,
            CsvColumnData::Text(_) => ColumnType::Text,
        }
    }
}
/// A DataFrame-like structure with named, typed columns.
///
/// Each column carries its inferred type. Construct via
/// [`CsvDataFrame::from_csv`] or [`CsvDataFrame::from_str`].
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CsvDataFrame {
    /// Column names, in order.
    pub column_names: Vec<String>,
    /// Typed column data (parallel to `column_names`).
    pub columns: Vec<CsvColumnData>,
}
#[allow(dead_code)]
impl CsvDataFrame {
    /// Build a `CsvDataFrame` from an existing [`CsvFile`] with auto type inference.
    pub fn from_csv(csv: &CsvFile) -> Self {
        let mut column_names = csv.headers.clone();
        let mut columns: Vec<CsvColumnData> = Vec::with_capacity(csv.headers.len());
        for col_idx in 0..csv.headers.len() {
            let col_type = csv.infer_column_type(col_idx);
            let col_data = match col_type {
                ColumnType::Integer => {
                    let vals: Vec<i64> = csv
                        .records
                        .iter()
                        .map(|r| {
                            r.fields
                                .get(col_idx)
                                .and_then(|s| s.trim().parse::<i64>().ok())
                                .unwrap_or(0)
                        })
                        .collect();
                    CsvColumnData::Integer(vals)
                }
                ColumnType::Float => {
                    let vals: Vec<f64> = csv
                        .records
                        .iter()
                        .map(|r| {
                            r.fields
                                .get(col_idx)
                                .and_then(|s| s.trim().parse::<f64>().ok())
                                .unwrap_or(f64::NAN)
                        })
                        .collect();
                    CsvColumnData::Float(vals)
                }
                ColumnType::Text => {
                    let vals: Vec<String> = csv
                        .records
                        .iter()
                        .map(|r| {
                            r.fields
                                .get(col_idx)
                                .map(|s| s.trim().to_string())
                                .unwrap_or_default()
                        })
                        .collect();
                    CsvColumnData::Text(vals)
                }
            };
            columns.push(col_data);
        }
        for (i, name) in column_names.iter_mut().enumerate() {
            if name.is_empty() {
                *name = format!("col_{}", i);
            }
        }
        CsvDataFrame {
            column_names,
            columns,
        }
    }
    /// Parse a CSV string and build a `CsvDataFrame` with type inference.
    ///
    /// Returns an error message with the problematic line number if parsing fails.
    pub fn from_str(s: &str) -> std::result::Result<Self, String> {
        let csv = CsvFile::from_str(s).map_err(|e| format!("line 1: {e}"))?;
        Ok(Self::from_csv(&csv))
    }
    /// Parse a delimiter-separated string and build a `CsvDataFrame`.
    pub fn from_str_with_delimiter(s: &str, delim: char) -> std::result::Result<Self, String> {
        let csv = CsvFile::from_str_with_delimiter(s, delim).map_err(|e| format!("line 1: {e}"))?;
        Ok(Self::from_csv(&csv))
    }
    /// Number of rows.
    pub fn n_rows(&self) -> usize {
        self.columns.first().map(|c| c.len()).unwrap_or(0)
    }
    /// Number of columns.
    pub fn n_cols(&self) -> usize {
        self.columns.len()
    }
    /// Return the index of a column by name, or `None`.
    pub fn column_index(&self, name: &str) -> Option<usize> {
        self.column_names.iter().position(|n| n == name)
    }
    /// Return a reference to the typed data of a column by index.
    pub fn column(&self, idx: usize) -> Option<&CsvColumnData> {
        self.columns.get(idx)
    }
    /// Return a reference to the typed data of a column by name.
    pub fn column_by_name(&self, name: &str) -> Option<&CsvColumnData> {
        let idx = self.column_index(name)?;
        self.column(idx)
    }
    /// Extract a float column by name.  Returns `None` if not found or wrong type.
    pub fn float_column(&self, name: &str) -> Option<&Vec<f64>> {
        match self.column_by_name(name)? {
            CsvColumnData::Float(v) => Some(v),
            _ => None,
        }
    }
    /// Extract an integer column by name.
    pub fn integer_column(&self, name: &str) -> Option<&Vec<i64>> {
        match self.column_by_name(name)? {
            CsvColumnData::Integer(v) => Some(v),
            _ => None,
        }
    }
    /// Extract a text column by name.
    pub fn text_column(&self, name: &str) -> Option<&Vec<String>> {
        match self.column_by_name(name)? {
            CsvColumnData::Text(v) => Some(v),
            _ => None,
        }
    }
    /// Serialize the DataFrame back to a CSV string (all values as strings).
    pub fn to_csv_string(&self) -> String {
        let mut out = self.column_names.join(",");
        out.push('\n');
        let n_rows = self.n_rows();
        for row in 0..n_rows {
            let fields: Vec<String> = self
                .columns
                .iter()
                .map(|col| match col {
                    CsvColumnData::Integer(v) => {
                        v.get(row).map(|x| x.to_string()).unwrap_or_default()
                    }
                    CsvColumnData::Float(v) => {
                        v.get(row).map(|x| format!("{}", x)).unwrap_or_default()
                    }
                    CsvColumnData::Text(v) => v.get(row).cloned().unwrap_or_default(),
                })
                .collect();
            out.push_str(&fields.join(","));
            out.push('\n');
        }
        out
    }
}
/// A line-by-line streaming CSV reader backed by an in-memory string.
///
/// Unlike [`LazyCsvIter`], this reader tracks the current row number and
/// supports peeking at headers before iterating data rows.
#[allow(dead_code)]
pub struct StreamingCsvReader<'a> {
    /// Delimiter detected or supplied.
    pub delimiter: char,
    /// Parsed column headers.
    pub headers: Vec<String>,
    pub(super) lines: std::str::Lines<'a>,
    /// Current row (0-based, not counting header).
    pub(super) row: usize,
}
#[allow(dead_code)]
impl<'a> StreamingCsvReader<'a> {
    /// Create a reader with explicit delimiter.
    pub fn new(input: &'a str, delimiter: char) -> Self {
        let mut lines = input.lines();
        let headers = lines
            .next()
            .map(|h| {
                h.split(delimiter)
                    .map(str::trim)
                    .map(String::from)
                    .collect()
            })
            .unwrap_or_default();
        Self {
            delimiter,
            headers,
            lines,
            row: 0,
        }
    }
    /// Create a reader with auto-detected delimiter.
    pub fn auto(input: &'a str) -> Self {
        let delim = detect_delimiter(input);
        Self::new(input, delim)
    }
    /// Number of columns (from header).
    pub fn n_cols(&self) -> usize {
        self.headers.len()
    }
    /// Current row number (rows consumed so far).
    pub fn current_row(&self) -> usize {
        self.row
    }
    /// Read the next data row as a `Vec`String`, or `None` at EOF.
    pub fn next_row(&mut self) -> Option<Vec<String>> {
        loop {
            let line = self.lines.next()?;
            if line.trim().is_empty() {
                continue;
            }
            self.row += 1;
            return Some(
                line.split(self.delimiter)
                    .map(str::trim)
                    .map(String::from)
                    .collect(),
            );
        }
    }
    /// Read all remaining rows into a [`CsvFile`].
    pub fn collect_all(mut self) -> CsvFile {
        let mut file = CsvFile::new(self.headers.clone());
        while let Some(fields) = self.next_row() {
            file.add_record(fields);
        }
        file
    }
}
/// Statistics for a numeric column.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ColumnStats {
    /// Minimum value.
    pub min: f64,
    /// Maximum value.
    pub max: f64,
    /// Arithmetic mean.
    pub mean: f64,
    /// Number of numeric values used for the statistics.
    pub count: usize,
    /// Sum of all numeric values.
    pub sum: f64,
}
/// A single frame in a 3D coordinate trajectory CSV.
///
/// One frame per "block", where each block has optional comment/header
/// lines followed by rows of `x,y,z` coordinates (one atom per row).
/// Blocks are separated by a blank line.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct TrajectoryFrame {
    /// Optional frame title/comment.
    pub title: String,
    /// 3D positions for each atom: `[x, y, z]` in simulation units.
    pub positions: Vec<[f64; 3]>,
}
#[allow(dead_code)]
impl TrajectoryFrame {
    /// Create an empty frame.
    pub fn new() -> Self {
        Self::default()
    }
    /// Number of atoms in this frame.
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }
}
