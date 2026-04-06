//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::should_implement_trait)]
#[allow(unused_imports)]
use super::functions::*;
use crate::{Error, Result};
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Configuration for [`ConfigurableCsvWriter`].
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CsvWriterConfig {
    /// Field delimiter (default `,`).
    pub delimiter: char,
    /// Line ending (default `\n`).
    pub line_ending: String,
    /// Whether to quote all fields (default `false`).
    pub quote_all: bool,
    /// Float precision (default 6 decimal places).
    pub precision: usize,
}
/// Full RFC-4180-compliant CSV parser.
///
/// Handles:
/// - Quoted fields containing delimiters, newlines and escaped quotes (`""`)
/// - Backslash escape sequences inside quoted fields (`\"`, `\\`, `\n`, `\t`, `\r`)
/// - Multi-line records (newlines inside quoted fields)
/// - Comment lines (optional, configurable prefix)
/// - Configurable delimiter
///
/// # Example
///
/// ```no_run
/// use oxiphysics_io::csv::CsvParser;
///
/// let data = "name,notes\nAlice,\"line1\nline2\"\nBob,\"he said \"\"hi\"\"\"\n";
/// let parser = CsvParser::new(data, ',');
/// let records = parser.parse_all().unwrap();
/// assert_eq!(records[0].get(0), "name");
/// assert_eq!(records[1].get(1), "line1\nline2");
/// assert_eq!(records[2].get(1), "he said \"hi\"");
/// ```
#[allow(dead_code)]
pub struct CsvParser<'a> {
    pub(super) input: &'a str,
    pub(super) delimiter: char,
    pub(super) comment_prefix: Option<char>,
}
#[allow(dead_code)]
impl<'a> CsvParser<'a> {
    /// Create a new parser for `input` with the given `delimiter`.
    pub fn new(input: &'a str, delimiter: char) -> Self {
        Self {
            input,
            delimiter,
            comment_prefix: None,
        }
    }
    /// Enable comment-line skipping with the given prefix character (e.g. `'#'`).
    pub fn with_comment_prefix(mut self, prefix: char) -> Self {
        self.comment_prefix = Some(prefix);
        self
    }
    /// Parse all records from the input.
    ///
    /// Returns a `Vec`CsvRecord` where the first record is usually the header.
    pub fn parse_all(self) -> std::result::Result<Vec<CsvRecord>, Error> {
        let mut records = Vec::new();
        let mut chars = self.input.chars().peekable();
        'outer: loop {
            if chars.peek().is_none() {
                break;
            }
            let mut fields: Vec<String> = Vec::new();
            let mut field = String::new();
            let mut in_quotes = false;
            loop {
                match chars.next() {
                    None => {
                        if in_quotes {
                            return Err(Error::Parse(
                                "unterminated quoted field at EOF".to_string(),
                            ));
                        }
                        fields.push(field);
                        break;
                    }
                    Some('"') if !in_quotes => {
                        in_quotes = true;
                    }
                    Some('"') if in_quotes => {
                        if chars.peek() == Some(&'"') {
                            chars.next();
                            field.push('"');
                        } else {
                            in_quotes = false;
                        }
                    }
                    Some('\\') if in_quotes => match chars.next() {
                        Some('n') => field.push('\n'),
                        Some('t') => field.push('\t'),
                        Some('r') => field.push('\r'),
                        Some('"') => field.push('"'),
                        Some('\\') => field.push('\\'),
                        Some(c) => {
                            field.push('\\');
                            field.push(c);
                        }
                        None => {
                            return Err(Error::Parse("trailing backslash at EOF".to_string()));
                        }
                    },
                    Some('\r') if !in_quotes => {
                        if chars.peek() == Some(&'\n') {
                            chars.next();
                        }
                        fields.push(field);
                        break;
                    }
                    Some('\n') if !in_quotes => {
                        fields.push(field);
                        break;
                    }
                    Some(c) if c == self.delimiter && !in_quotes => {
                        fields.push(field.clone());
                        field = String::new();
                    }
                    Some(c) => {
                        field.push(c);
                    }
                    #[allow(unreachable_patterns)]
                    _ => {}
                }
            }
            if let Some(prefix) = self.comment_prefix
                && fields
                    .first()
                    .map(|f| f.trim_start().starts_with(prefix))
                    .unwrap_or(false)
            {
                continue 'outer;
            }
            if fields.len() == 1 && fields[0].trim().is_empty() {
                continue 'outer;
            }
            records.push(CsvRecord { fields });
        }
        Ok(records)
    }
}
/// A parsed CSV record (one logical row, possibly spanning multiple physical lines).
#[derive(Debug, Clone, PartialEq)]
pub struct CsvRecord {
    /// Raw field strings (unquoted, unescape-handled).
    pub fields: Vec<String>,
}
impl CsvRecord {
    /// Number of fields in this record.
    pub fn len(&self) -> usize {
        self.fields.len()
    }
    /// Returns `true` if this record has no fields.
    pub fn is_empty(&self) -> bool {
        self.fields.is_empty()
    }
    /// Get a field by index, returning an empty string for out-of-range indices.
    pub fn get(&self, index: usize) -> &str {
        self.fields.get(index).map(|s| s.as_str()).unwrap_or("")
    }
}
/// Aggregation function for pivot table values.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PivotAgg {
    /// Sum of values.
    Sum,
    /// Arithmetic mean.
    Mean,
    /// Count of non-empty values.
    Count,
    /// Minimum value.
    Min,
    /// Maximum value.
    Max,
}
/// An in-memory table of string-valued cells.
///
/// Columns are named; rows are `Vec`String`. Provides the foundation for
/// type inference, merge/join, pivot, and diff operations.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CsvTable {
    /// Column headers.
    pub headers: Vec<String>,
    /// Data rows (each row has `headers.len()` fields, padded with `""` if short).
    pub rows: Vec<Vec<String>>,
}
#[allow(dead_code)]
impl CsvTable {
    /// Create an empty table with the given headers.
    pub fn new(headers: Vec<String>) -> Self {
        Self {
            headers,
            rows: Vec::new(),
        }
    }
    /// Parse a CSV string into a `CsvTable` using the given delimiter.
    pub fn from_str(data: &str, delimiter: char) -> std::result::Result<Self, Error> {
        let parser = CsvParser::new(data, delimiter).with_comment_prefix('#');
        let mut records = parser.parse_all()?;
        if records.is_empty() {
            return Err(Error::Parse("CSV table is empty".to_string()));
        }
        let header_rec = records.remove(0);
        let headers: Vec<String> = header_rec
            .fields
            .iter()
            .map(|s| s.trim().to_string())
            .collect();
        let ncols = headers.len();
        let mut rows = Vec::new();
        for rec in records {
            let mut row: Vec<String> = rec
                .fields
                .into_iter()
                .map(|s| s.trim().to_string())
                .collect();
            while row.len() < ncols {
                row.push(String::new());
            }
            rows.push(row);
        }
        Ok(Self { headers, rows })
    }
    /// Serialize back to a CSV string with the given delimiter.
    pub fn to_csv_string(&self, delimiter: char) -> String {
        let mut out = String::new();
        out.push_str(&self.headers.join(&delimiter.to_string()));
        out.push('\n');
        for row in &self.rows {
            let line: Vec<String> = row.iter().map(|f| quote_field(f, delimiter)).collect();
            out.push_str(&line.join(&delimiter.to_string()));
            out.push('\n');
        }
        out
    }
    /// Look up a column index by name.
    pub fn column_index(&self, name: &str) -> std::result::Result<usize, Error> {
        self.headers
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| Error::Parse(format!("column '{}' not found", name)))
    }
    /// Get all values from a named column as strings.
    pub fn column_values(&self, name: &str) -> std::result::Result<Vec<&str>, Error> {
        let idx = self.column_index(name)?;
        Ok(self.rows.iter().map(|r| r[idx].as_str()).collect())
    }
    /// Get all values from a named column parsed as `f64`.
    pub fn column_f64(&self, name: &str) -> std::result::Result<Vec<f64>, Error> {
        let idx = self.column_index(name)?;
        self.rows
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let s = r[idx].trim();
                if s.is_empty() {
                    Ok(f64::NAN)
                } else {
                    s.parse::<f64>().map_err(|_| {
                        Error::Parse(format!("row {}: cannot parse '{}' as f64", i, s))
                    })
                }
            })
            .collect()
    }
    /// Number of data rows (excluding header).
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }
    /// Number of columns.
    pub fn col_count(&self) -> usize {
        self.headers.len()
    }
}
/// In-memory CSV writer with configurable delimiter and floating-point precision.
#[allow(dead_code)]
pub struct InMemoryCsvWriter {
    pub(super) columns: Vec<String>,
    pub(super) delimiter: char,
    pub(super) precision: usize,
    pub(super) rows: Vec<Vec<f64>>,
}
#[allow(dead_code)]
impl InMemoryCsvWriter {
    /// Create a new `InMemoryCsvWriter` with the given column names and delimiter.
    pub fn new(columns: &[&str], delimiter: char) -> Self {
        Self {
            columns: columns.iter().map(|s| s.to_string()).collect(),
            delimiter,
            precision: 6,
            rows: Vec::new(),
        }
    }
    /// Set the floating-point precision used when formatting values.
    pub fn with_precision(mut self, precision: usize) -> Self {
        self.precision = precision;
        self
    }
    /// Format a single row of values as a delimited string (no newline).
    pub fn write_row(&self, values: &[f64]) -> std::result::Result<String, Error> {
        if values.len() != self.columns.len() {
            return Err(Error::Parse(format!(
                "expected {} values, got {}",
                self.columns.len(),
                values.len()
            )));
        }
        let parts: Vec<String> = values
            .iter()
            .map(|v| format!("{:.prec$}", v, prec = self.precision))
            .collect();
        Ok(parts.join(&self.delimiter.to_string()))
    }
    /// Return the header line (no newline).
    pub fn write_header(&self) -> String {
        self.columns.join(&self.delimiter.to_string())
    }
    /// Accumulate a row.
    pub fn add_row(&mut self, values: Vec<f64>) {
        self.rows.push(values);
    }
    /// Render the entire accumulated dataset (header + all rows) as a String.
    pub fn write_all(&self, rows: &[Vec<f64>]) -> String {
        let mut out = self.write_header();
        out.push('\n');
        for row in rows {
            if let Ok(line) = self.write_row(row) {
                out.push_str(&line);
                out.push('\n');
            }
        }
        out
    }
}
/// Reader for CSV files with numeric data.
pub struct CsvReader;
impl CsvReader {
    /// Read a CSV file, returning `(headers, rows)`.
    ///
    /// Assumes the first line is a header and all subsequent lines are numeric.
    pub fn read(path: &str) -> Result<(Vec<String>, Vec<Vec<f64>>)> {
        let file = File::open(Path::new(path))?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();
        let header_line = lines
            .next()
            .ok_or_else(|| Error::Parse("empty CSV file".to_string()))??;
        let headers: Vec<String> = header_line
            .split(',')
            .map(|s| s.trim().to_string())
            .collect();
        let mut rows = Vec::new();
        for line in lines {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }
            let row: std::result::Result<Vec<f64>, _> = trimmed
                .split(',')
                .map(|s| s.trim().parse::<f64>())
                .collect();
            let row = row.map_err(|e| Error::Parse(e.to_string()))?;
            rows.push(row);
        }
        Ok((headers, rows))
    }
}
/// Inferred column type after scanning all values in a column.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ColumnType {
    /// All non-empty values parse as integers.
    Integer,
    /// All non-empty values parse as floats (but not all as integers).
    Float,
    /// All non-empty values are "true"/"false"/"yes"/"no"/"1"/"0".
    Boolean,
    /// At least one value cannot be parsed as a number or boolean.
    Text,
    /// The column contains only empty values.
    Empty,
}
/// A row-at-a-time streaming CSV parser for large files.
///
/// Reads the file line by line with minimal memory allocation. Quoted
/// multi-line fields are assembled across physical lines before yielding a
/// logical record.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_io::csv::CsvStreamParser;
///
/// let mut parser = CsvStreamParser::open("/tmp/big.csv", ',').unwrap();
/// let headers = parser.headers().to_vec();
/// while let Some(record) = parser.next_record().unwrap() {
///     // process record.fields
///     let _ = record.len();
/// }
/// ```
#[allow(dead_code)]
pub struct CsvStreamParser {
    pub(super) reader: BufReader<File>,
    pub(super) delimiter: char,
    pub(super) headers: Vec<String>,
    /// Buffer for partial multi-line records across physical lines.
    pub(super) pending: String,
    pub(super) open_quotes: bool,
}
#[allow(dead_code)]
impl CsvStreamParser {
    /// Open a CSV file for streaming parsing.
    pub fn open(path: &str, delimiter: char) -> std::result::Result<Self, Error> {
        let file = File::open(Path::new(path))?;
        let mut reader = BufReader::new(file);
        let mut header_line = String::new();
        reader.read_line(&mut header_line)?;
        let headers = split_csv_line(
            header_line.trim_end_matches('\n').trim_end_matches('\r'),
            delimiter,
        )
        .into_iter()
        .map(|s| s.trim().to_string())
        .collect();
        Ok(Self {
            reader,
            delimiter,
            headers,
            pending: String::new(),
            open_quotes: false,
        })
    }
    /// Return the column headers read from the first row.
    pub fn headers(&self) -> &[String] {
        &self.headers
    }
    /// Read and return the next logical CSV record.
    ///
    /// Returns `Ok(None)` at EOF, `Err(_)` on I/O or parse errors.
    pub fn next_record(&mut self) -> std::result::Result<Option<CsvRecord>, Error> {
        loop {
            let mut line = String::new();
            let bytes_read = self.reader.read_line(&mut line)?;
            if bytes_read == 0 {
                if self.pending.is_empty() {
                    return Ok(None);
                }
                if self.open_quotes {
                    return Err(Error::Parse("unterminated quoted field at EOF".to_string()));
                }
                let record = self.flush_pending()?;
                return Ok(Some(record));
            }
            for ch in line.chars() {
                if ch == '"' {
                    self.open_quotes = !self.open_quotes;
                }
            }
            self.pending.push_str(&line);
            if !self.open_quotes {
                let record = self.flush_pending()?;
                if record.fields.len() == 1 && record.fields[0].trim().is_empty() {
                    continue;
                }
                return Ok(Some(record));
            }
        }
    }
    fn flush_pending(&mut self) -> std::result::Result<CsvRecord, Error> {
        let line = std::mem::take(&mut self.pending);
        let trimmed = line.trim_end_matches('\n').trim_end_matches('\r');
        let fields = split_csv_line(trimmed, self.delimiter);
        Ok(CsvRecord {
            fields: fields.into_iter().map(|s| s.trim().to_string()).collect(),
        })
    }
}
/// In-memory CSV reader that parses a complete CSV string.
///
/// Supports:
/// - configurable delimiter
/// - comment lines starting with `#`
/// - empty fields (treated as `f64::NAN`)
/// - quoted strings (double-quoted fields with escaped commas)
#[allow(dead_code)]
pub struct InMemoryCsvReader {
    pub(super) headers: Vec<String>,
    pub(super) rows: Vec<Vec<Option<f64>>>,
}
#[allow(dead_code)]
impl InMemoryCsvReader {
    /// Parse a CSV string with the default comma delimiter.
    pub fn from_str(data: &str) -> std::result::Result<Self, Error> {
        Self::parse_with_delimiter(data, ',')
    }
    /// Parse a CSV string with a custom delimiter.
    pub fn parse_with_delimiter(data: &str, delim: char) -> std::result::Result<Self, Error> {
        let mut non_empty_lines: Vec<&str> = data
            .lines()
            .filter(|l| {
                let t = l.trim();
                !t.is_empty() && !t.starts_with('#')
            })
            .collect();
        if non_empty_lines.is_empty() {
            return Err(Error::Parse("CSV input is empty".to_string()));
        }
        let header_line = non_empty_lines.remove(0);
        let headers: Vec<String> = split_csv_line(header_line, delim)
            .into_iter()
            .map(|s| s.trim().to_string())
            .collect();
        if headers.is_empty() {
            return Err(Error::Parse("no headers found".to_string()));
        }
        let mut rows: Vec<Vec<Option<f64>>> = Vec::new();
        for line in &non_empty_lines {
            let fields = split_csv_line(line, delim);
            let parsed: Vec<Option<f64>> = fields
                .iter()
                .map(|f| {
                    let t = f.trim();
                    if t.is_empty() {
                        None
                    } else {
                        t.parse::<f64>().ok()
                    }
                })
                .collect();
            rows.push(parsed);
        }
        Ok(Self { headers, rows })
    }
    /// Return the values in the named column as `f64`.
    ///
    /// Missing (`None`) values are replaced with `f64::NAN`.
    pub fn get_column_f64(&self, name: &str) -> std::result::Result<Vec<f64>, Error> {
        let idx = self
            .headers
            .iter()
            .position(|h| h == name)
            .ok_or_else(|| Error::Parse(format!("column '{}' not found", name)))?;
        let col: Vec<f64> = self
            .rows
            .iter()
            .map(|row| row.get(idx).copied().flatten().unwrap_or(f64::NAN))
            .collect();
        Ok(col)
    }
    /// Return the number of data rows (excluding header).
    pub fn get_row_count(&self) -> usize {
        self.rows.len()
    }
    /// Return the column headers.
    pub fn headers(&self) -> &[String] {
        &self.headers
    }
    /// Compute (min, max, mean, std) for a named column, ignoring NaN values.
    pub fn column_stats(&self, name: &str) -> std::result::Result<(f64, f64, f64, f64), Error> {
        let col = self.get_column_f64(name)?;
        let valid: Vec<f64> = col.into_iter().filter(|v| !v.is_nan()).collect();
        if valid.is_empty() {
            return Err(Error::Parse(format!(
                "column '{}' has no valid numeric data",
                name
            )));
        }
        let n = valid.len() as f64;
        let min = valid.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = valid.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mean = valid.iter().sum::<f64>() / n;
        let variance = valid.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
        let std = variance.sqrt();
        Ok((min, max, mean, std))
    }
}
/// A CSV reader that exposes per-column typed access after type inference.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_io::csv::TypedCsvReader;
///
/// let data = "id,x,active\n1,3.14,true\n2,2.71,false\n";
/// let reader = TypedCsvReader::from_str(data).unwrap();
/// let ids = reader.column_as_i64("id").unwrap();
/// assert_eq!(ids, vec![1, 2]);
/// let xs = reader.column_as_f64("x").unwrap();
/// assert!((xs[0] - 3.14).abs() < 1e-10);
/// let active = reader.column_as_bool("active").unwrap();
/// assert_eq!(active, vec![true, false]);
/// ```
#[allow(dead_code)]
pub struct TypedCsvReader {
    pub(super) table: CsvTable,
}
#[allow(dead_code)]
impl TypedCsvReader {
    /// Parse a CSV string with the default comma delimiter.
    pub fn from_str(data: &str) -> std::result::Result<Self, Error> {
        let table = CsvTable::from_str(data, ',')?;
        Ok(Self { table })
    }
    /// Parse a CSV string with a custom delimiter.
    pub fn with_delimiter(data: &str, delimiter: char) -> std::result::Result<Self, Error> {
        let table = CsvTable::from_str(data, delimiter)?;
        Ok(Self { table })
    }
    /// Return the inferred type of a named column.
    pub fn column_type(&self, name: &str) -> std::result::Result<ColumnType, Error> {
        let idx = self.table.column_index(name)?;
        let values: Vec<&str> = self.table.rows.iter().map(|r| r[idx].as_str()).collect();
        Ok(infer_column_type(&values))
    }
    /// Return a named column parsed as `i64`.
    pub fn column_as_i64(&self, name: &str) -> std::result::Result<Vec<i64>, Error> {
        let idx = self.table.column_index(name)?;
        self.table
            .rows
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let s = r[idx].trim();
                s.parse::<i64>()
                    .map_err(|_| Error::Parse(format!("row {}: cannot parse '{}' as i64", i, s)))
            })
            .collect()
    }
    /// Return a named column parsed as `f64`.
    pub fn column_as_f64(&self, name: &str) -> std::result::Result<Vec<f64>, Error> {
        self.table.column_f64(name)
    }
    /// Return a named column parsed as `bool`.
    ///
    /// Accepts "true"/"1"/"yes" → `true`; "false"/"0"/"no" → `false`.
    pub fn column_as_bool(&self, name: &str) -> std::result::Result<Vec<bool>, Error> {
        let idx = self.table.column_index(name)?;
        self.table
            .rows
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let s = r[idx].trim().to_lowercase();
                match s.as_str() {
                    "true" | "1" | "yes" => Ok(true),
                    "false" | "0" | "no" => Ok(false),
                    other => Err(Error::Parse(format!(
                        "row {}: cannot parse '{}' as bool",
                        i, other
                    ))),
                }
            })
            .collect()
    }
    /// Return the underlying table reference.
    pub fn table(&self) -> &CsvTable {
        &self.table
    }
    /// Return the column headers.
    pub fn headers(&self) -> &[String] {
        &self.table.headers
    }
    /// Number of data rows.
    pub fn row_count(&self) -> usize {
        self.table.row_count()
    }
}
/// The result of comparing two CSV datasets.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CsvDiff {
    /// Rows present in `left` but not in `right` (by key column value).
    pub removed: Vec<Vec<String>>,
    /// Rows present in `right` but not in `left`.
    pub added: Vec<Vec<String>>,
    /// Rows where the key exists in both but non-key columns differ.
    pub changed: Vec<CsvChangedRow>,
}
/// Append-mode CSV writer for time series data.
pub struct CsvWriter {
    pub(super) writer: BufWriter<File>,
}
impl CsvWriter {
    /// Create a new CSV file with the given headers.
    pub fn new(path: &str, headers: &[&str]) -> Result<Self> {
        let file = File::create(Path::new(path))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "{}", headers.join(","))?;
        writer.flush()?;
        let file = OpenOptions::new().append(true).open(Path::new(path))?;
        let writer = BufWriter::new(file);
        Ok(Self { writer })
    }
    /// Append one row of values.
    pub fn write_row(&mut self, values: &[f64]) -> Result<()> {
        let row: Vec<String> = values.iter().map(|v| v.to_string()).collect();
        writeln!(self.writer, "{}", row.join(","))?;
        self.writer.flush()?;
        Ok(())
    }
}
/// A single changed row in a [`CsvDiff`].
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CsvChangedRow {
    /// The key value identifying this row.
    pub key: String,
    /// The row as it appeared in the `left` table.
    pub before: Vec<String>,
    /// The row as it appeared in the `right` table.
    pub after: Vec<String>,
}
/// CSV writer with configurable delimiter, quoting, line endings and precision.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_io::csv::{ConfigurableCsvWriter, CsvWriterConfig};
///
/// let cfg = CsvWriterConfig { delimiter: ';', precision: 3, ..Default::default() };
/// let mut w = ConfigurableCsvWriter::new(cfg);
/// w.write_header(&["x", "y"]);
/// w.write_f64_row(&[1.0, 2.5]);
/// let out = w.finish();
/// assert!(out.starts_with("x;y"));
/// assert!(out.contains("1.000;2.500"));
/// ```
#[allow(dead_code)]
pub struct ConfigurableCsvWriter {
    pub(super) config: CsvWriterConfig,
    pub(super) buffer: String,
}
#[allow(dead_code)]
impl ConfigurableCsvWriter {
    /// Create a new writer with the given configuration.
    pub fn new(config: CsvWriterConfig) -> Self {
        Self {
            config,
            buffer: String::new(),
        }
    }
    /// Write a header row from string slices.
    pub fn write_header(&mut self, headers: &[&str]) {
        let line = if self.config.quote_all {
            headers
                .iter()
                .map(|h| format!("\"{}\"", h.replace('"', "\"\"")))
                .collect::<Vec<_>>()
                .join(&self.config.delimiter.to_string())
        } else {
            headers
                .iter()
                .map(|h| quote_field(h, self.config.delimiter))
                .collect::<Vec<_>>()
                .join(&self.config.delimiter.to_string())
        };
        self.buffer.push_str(&line);
        self.buffer.push_str(&self.config.line_ending);
    }
    /// Write a row of `f64` values.
    pub fn write_f64_row(&mut self, values: &[f64]) {
        let prec = self.config.precision;
        let line: Vec<String> = values
            .iter()
            .map(|v| format!("{:.prec$}", v, prec = prec))
            .collect();
        self.buffer
            .push_str(&line.join(&self.config.delimiter.to_string()));
        self.buffer.push_str(&self.config.line_ending);
    }
    /// Write a row of string values.
    pub fn write_str_row(&mut self, values: &[&str]) {
        let delim = self.config.delimiter;
        let line: Vec<String> = values
            .iter()
            .map(|v| {
                if self.config.quote_all {
                    format!("\"{}\"", v.replace('"', "\"\""))
                } else {
                    quote_field(v, delim)
                }
            })
            .collect();
        self.buffer
            .push_str(&line.join(&self.config.delimiter.to_string()));
        self.buffer.push_str(&self.config.line_ending);
    }
    /// Consume the writer and return the accumulated CSV string.
    pub fn finish(self) -> String {
        self.buffer
    }
}
