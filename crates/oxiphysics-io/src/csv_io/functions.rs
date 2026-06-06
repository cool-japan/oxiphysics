//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{AggOp, CsvFile, CsvRecord, CsvSchema, CsvValidationReport, TrajectoryFrame};

/// Build a two-column time-series CSV string from parallel slices.
pub fn write_timeseries(
    times: &[f64],
    values: &[f64],
    header_time: &str,
    header_val: &str,
) -> String {
    let mut csv = CsvFile::new(vec![header_time.to_string(), header_val.to_string()]);
    for (&t, &v) in times.iter().zip(values.iter()) {
        csv.add_record_f64(&[t, v]);
    }
    csv.to_string()
}
/// Extract a named column from a `CsvFile` as `f64`.
pub fn parse_column_f64(csv: &CsvFile, name: &str) -> Result<Vec<f64>, String> {
    let idx = csv
        .get_column_by_name(name)
        .ok_or_else(|| format!("Column '{}' not found", name))?;
    csv.get_column_f64(idx)
}
/// Auto-detect the delimiter used in a CSV string.
///
/// Tests comma, tab, semicolon, and pipe. Returns the delimiter that produces
/// the most consistent column count across the first few lines.
pub fn detect_delimiter(s: &str) -> char {
    let candidates = [',', '\t', ';', '|'];
    let lines: Vec<&str> = s.lines().take(10).collect();
    if lines.is_empty() {
        return ',';
    }
    let mut best_delim = ',';
    let mut best_score: usize = 0;
    for &delim in &candidates {
        let counts: Vec<usize> = lines.iter().map(|l| l.split(delim).count()).collect();
        if counts.is_empty() {
            continue;
        }
        let first = counts[0];
        if first < 2 {
            continue;
        }
        let consistent = counts.iter().filter(|&&c| c == first).count();
        let score = consistent * first;
        if score > best_score {
            best_score = score;
            best_delim = delim;
        }
    }
    best_delim
}
/// Parse a CSV string with auto-detected delimiter.
pub fn parse_auto(s: &str) -> Result<CsvFile, String> {
    let delim = detect_delimiter(s);
    CsvFile::from_str_with_delimiter(s, delim)
}
/// Read a CSV in chunks: returns an iterator of `CsvFile` objects, each
/// containing at most `chunk_size` records.
pub fn read_chunked(s: &str, chunk_size: usize) -> Vec<CsvFile> {
    let full = match s.parse::<CsvFile>() {
        Ok(f) => f,
        Err(_) => return Vec::new(),
    };
    if chunk_size == 0 {
        return vec![full];
    }
    let mut chunks = Vec::new();
    let mut start = 0;
    while start < full.records.len() {
        let end = (start + chunk_size).min(full.records.len());
        let mut chunk = CsvFile::new(full.headers.clone());
        for i in start..end {
            chunk.records.push(CsvRecord {
                fields: full.records[i].fields.clone(),
            });
        }
        chunks.push(chunk);
        start = end;
    }
    chunks
}
/// Normalize a header string: lowercase, replace non-alphanumeric with underscores.
pub fn normalize_header(s: &str) -> String {
    s.trim()
        .to_lowercase()
        .chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}
/// Aggregate a numeric column of a [`CsvFile`] using the given [`AggOp`].
///
/// Returns `None` if the column index is out of range or no numeric values
/// are found.
pub fn aggregate_column(csv: &CsvFile, col: usize, op: AggOp) -> Option<f64> {
    let vals = csv.get_column_f64(col).ok()?;
    if vals.is_empty() {
        return None;
    }
    let n = vals.len() as f64;
    Some(match op {
        AggOp::Sum => vals.iter().sum(),
        AggOp::Mean => vals.iter().sum::<f64>() / n,
        AggOp::Min => vals.iter().cloned().fold(f64::INFINITY, f64::min),
        AggOp::Max => vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        AggOp::Std => {
            let mean = vals.iter().sum::<f64>() / n;
            let var = vals.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / n;
            var.sqrt()
        }
        AggOp::Count => vals.len() as f64,
    })
}
/// Validate a [`CsvFile`] against a [`CsvSchema`] and return a report.
pub fn validate_csv(csv: &CsvFile, schema: &CsvSchema) -> CsvValidationReport {
    CsvValidationReport {
        errors: schema.validate(csv),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::csv_io::types::*;
    #[test]
    fn test_new_empty() {
        let csv = CsvFile::new(vec!["x".into(), "y".into()]);
        assert_eq!(csv.column_count(), 2);
        assert_eq!(csv.record_count(), 0);
    }
    #[test]
    fn test_add_record_string() {
        let mut csv = CsvFile::new(vec!["a".into(), "b".into()]);
        csv.add_record(vec!["1".into(), "2".into()]);
        assert_eq!(csv.record_count(), 1);
        assert_eq!(csv.records[0].fields[0], "1");
    }
    #[test]
    fn test_add_record_f64() {
        let mut csv = CsvFile::new(vec!["t".into(), "v".into()]);
        csv.add_record_f64(&[0.0, 9.81]);
        assert_eq!(csv.record_count(), 1);
        assert_eq!(csv.records[0].fields[1], "9.81");
    }
    #[test]
    fn test_get_column_f64_valid() {
        let mut csv = CsvFile::new(vec!["x".into()]);
        csv.add_record_f64(&[1.5]);
        csv.add_record_f64(&[3.0]);
        let col = csv.get_column_f64(0).unwrap();
        assert_eq!(col, vec![1.5, 3.0]);
    }
    #[test]
    fn test_get_column_f64_out_of_range() {
        let csv = CsvFile::new(vec!["x".into()]);
        assert!(csv.get_column_f64(5).is_err());
    }
    #[test]
    fn test_get_column_f64_parse_error() {
        let mut csv = CsvFile::new(vec!["x".into()]);
        csv.add_record(vec!["not_a_number".into()]);
        assert!(csv.get_column_f64(0).is_err());
    }
    #[test]
    fn test_get_column_by_name_found() {
        let csv = CsvFile::new(vec!["time".into(), "energy".into()]);
        assert_eq!(csv.get_column_by_name("energy"), Some(1));
    }
    #[test]
    fn test_get_column_by_name_missing() {
        let csv = CsvFile::new(vec!["time".into()]);
        assert!(csv.get_column_by_name("missing").is_none());
    }
    #[test]
    fn test_to_string_roundtrip() {
        let mut csv = CsvFile::new(vec!["t".into(), "x".into()]);
        csv.add_record_f64(&[0.0, 1.0]);
        csv.add_record_f64(&[1.0, 2.0]);
        let s = csv.to_string();
        let parsed = &s.parse::<CsvFile>().unwrap();
        assert_eq!(parsed.headers, vec!["t", "x"]);
        assert_eq!(parsed.record_count(), 2);
    }
    #[test]
    fn test_from_str_with_spaces() {
        let s = "time , value\n0.0 , 1.0\n1.0 , 2.0\n";
        let csv = s.parse::<CsvFile>().unwrap();
        assert_eq!(csv.headers[0], "time");
        assert_eq!(csv.record_count(), 2);
    }
    #[test]
    fn test_from_str_empty_lines_ignored() {
        let s = "a,b\n1,2\n\n3,4\n";
        let csv = s.parse::<CsvFile>().unwrap();
        assert_eq!(csv.record_count(), 2);
    }
    #[test]
    fn test_from_str_empty_input() {
        assert!("".parse::<CsvFile>().is_err());
    }
    #[test]
    fn test_filter_rows_positive() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[-1.0]);
        csv.add_record_f64(&[2.0]);
        csv.add_record_f64(&[3.0]);
        let filtered = csv.filter_rows(0, |v| v > 0.0);
        assert_eq!(filtered.record_count(), 2);
    }
    #[test]
    fn test_filter_rows_none_match() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[1.0]);
        let filtered = csv.filter_rows(0, |v| v > 100.0);
        assert_eq!(filtered.record_count(), 0);
    }
    #[test]
    fn test_filter_rows_preserves_headers() {
        let mut csv = CsvFile::new(vec!["t".into(), "x".into()]);
        csv.add_record_f64(&[0.0, 1.0]);
        let filtered = csv.filter_rows(0, |_| true);
        assert_eq!(filtered.headers, vec!["t", "x"]);
    }
    #[test]
    fn test_write_timeseries_format() {
        let s = write_timeseries(&[0.0, 1.0], &[10.0, 20.0], "time", "pos");
        assert!(s.starts_with("time,pos\n"));
        assert!(s.contains("0,10") || s.contains("0.0") || s.contains("10"));
    }
    #[test]
    fn test_parse_column_f64_by_name() {
        let mut csv = CsvFile::new(vec!["t".into(), "e".into()]);
        csv.add_record_f64(&[0.0, 5.0]);
        csv.add_record_f64(&[1.0, 6.0]);
        let col = parse_column_f64(&csv, "e").unwrap();
        assert_eq!(col, vec![5.0, 6.0]);
    }
    #[test]
    fn test_parse_column_f64_missing_name() {
        let csv = CsvFile::new(vec!["t".into()]);
        assert!(parse_column_f64(&csv, "nope").is_err());
    }
    #[test]
    fn test_multiple_columns_round_trip() {
        let mut csv = CsvFile::new(vec!["x".into(), "y".into(), "z".into()]);
        csv.add_record_f64(&[1.0, 2.0, 3.0]);
        let s = csv.to_string();
        let parsed = &s.parse::<CsvFile>().unwrap();
        let x = parsed.get_column_f64(0).unwrap();
        let y = parsed.get_column_f64(1).unwrap();
        let z = parsed.get_column_f64(2).unwrap();
        assert!((x[0] - 1.0).abs() < 1e-12);
        assert!((y[0] - 2.0).abs() < 1e-12);
        assert!((z[0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_infer_column_type_integer() {
        let mut csv = CsvFile::new(vec!["a".into()]);
        csv.add_record(vec!["1".into()]);
        csv.add_record(vec!["2".into()]);
        csv.add_record(vec!["-10".into()]);
        assert_eq!(csv.infer_column_type(0), ColumnType::Integer);
    }
    #[test]
    fn test_infer_column_type_float() {
        let mut csv = CsvFile::new(vec!["a".into()]);
        csv.add_record(vec!["1.5".into()]);
        csv.add_record(vec!["2.7".into()]);
        assert_eq!(csv.infer_column_type(0), ColumnType::Float);
    }
    #[test]
    fn test_infer_column_type_text() {
        let mut csv = CsvFile::new(vec!["a".into()]);
        csv.add_record(vec!["hello".into()]);
        csv.add_record(vec!["world".into()]);
        assert_eq!(csv.infer_column_type(0), ColumnType::Text);
    }
    #[test]
    fn test_infer_column_type_mixed_int_float() {
        let mut csv = CsvFile::new(vec!["a".into()]);
        csv.add_record(vec!["1".into()]);
        csv.add_record(vec!["2.5".into()]);
        assert_eq!(csv.infer_column_type(0), ColumnType::Float);
    }
    #[test]
    fn test_infer_column_type_empty() {
        let csv = CsvFile::new(vec!["a".into()]);
        assert_eq!(csv.infer_column_type(0), ColumnType::Text);
    }
    #[test]
    fn test_infer_column_type_out_of_range() {
        let csv = CsvFile::new(vec!["a".into()]);
        assert_eq!(csv.infer_column_type(99), ColumnType::Text);
    }
    #[test]
    fn test_select_columns_by_index() {
        let mut csv = CsvFile::new(vec!["a".into(), "b".into(), "c".into()]);
        csv.add_record(vec!["1".into(), "2".into(), "3".into()]);
        let subset = csv.select_columns(&[0, 2]);
        assert_eq!(subset.headers, vec!["a", "c"]);
        assert_eq!(subset.records[0].fields, vec!["1", "3"]);
    }
    #[test]
    fn test_select_columns_by_name() {
        let mut csv = CsvFile::new(vec!["time".into(), "x".into(), "y".into()]);
        csv.add_record(vec!["0".into(), "1.0".into(), "2.0".into()]);
        let subset = csv.select_columns_by_name(&["y", "time"]);
        assert_eq!(subset.headers, vec!["y", "time"]);
        assert_eq!(subset.records[0].fields, vec!["2.0", "0"]);
    }
    #[test]
    fn test_select_columns_missing_name_ignored() {
        let mut csv = CsvFile::new(vec!["a".into(), "b".into()]);
        csv.add_record(vec!["1".into(), "2".into()]);
        let subset = csv.select_columns_by_name(&["a", "missing"]);
        assert_eq!(subset.headers, vec!["a"]);
    }
    #[test]
    fn test_normalize_headers() {
        let mut csv = CsvFile::new(vec![
            " Time Step ".into(),
            "X Position".into(),
            "energy (J)".into(),
        ]);
        csv.normalize_headers();
        assert_eq!(csv.headers[0], "time_step");
        assert_eq!(csv.headers[1], "x_position");
        assert_eq!(csv.headers[2], "energy__j_");
    }
    #[test]
    fn test_column_stats_basic() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[1.0]);
        csv.add_record_f64(&[3.0]);
        csv.add_record_f64(&[5.0]);
        let stats = csv.column_stats(0).unwrap();
        assert!((stats.min - 1.0).abs() < 1e-12);
        assert!((stats.max - 5.0).abs() < 1e-12);
        assert!((stats.mean - 3.0).abs() < 1e-12);
        assert_eq!(stats.count, 3);
        assert!((stats.sum - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_column_stats_single_value() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[42.0]);
        let stats = csv.column_stats(0).unwrap();
        assert!((stats.min - 42.0).abs() < 1e-12);
        assert!((stats.max - 42.0).abs() < 1e-12);
        assert!((stats.mean - 42.0).abs() < 1e-12);
    }
    #[test]
    fn test_column_stats_negative_values() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[-10.0]);
        csv.add_record_f64(&[-5.0]);
        csv.add_record_f64(&[0.0]);
        let stats = csv.column_stats(0).unwrap();
        assert!((stats.min - (-10.0)).abs() < 1e-12);
        assert!((stats.max - 0.0).abs() < 1e-12);
        assert!((stats.mean - (-5.0)).abs() < 1e-12);
    }
    #[test]
    fn test_column_stats_text_column_returns_none() {
        let mut csv = CsvFile::new(vec!["name".into()]);
        csv.add_record(vec!["alice".into()]);
        assert!(csv.column_stats(0).is_none());
    }
    #[test]
    fn test_all_column_stats() {
        let mut csv = CsvFile::new(vec!["x".into(), "label".into(), "y".into()]);
        csv.add_record(vec!["1.0".into(), "a".into(), "10.0".into()]);
        csv.add_record(vec!["2.0".into(), "b".into(), "20.0".into()]);
        let stats = csv.all_column_stats();
        assert_eq!(stats.len(), 2);
        assert_eq!(stats[0].0, "x");
        assert_eq!(stats[1].0, "y");
    }
    #[test]
    fn test_detect_delimiter_comma() {
        let s = "a,b,c\n1,2,3\n4,5,6\n";
        assert_eq!(detect_delimiter(s), ',');
    }
    #[test]
    fn test_detect_delimiter_tab() {
        let s = "a\tb\tc\n1\t2\t3\n4\t5\t6\n";
        assert_eq!(detect_delimiter(s), '\t');
    }
    #[test]
    fn test_detect_delimiter_semicolon() {
        let s = "a;b;c\n1;2;3\n4;5;6\n";
        assert_eq!(detect_delimiter(s), ';');
    }
    #[test]
    fn test_detect_delimiter_pipe() {
        let s = "a|b|c\n1|2|3\n";
        assert_eq!(detect_delimiter(s), '|');
    }
    #[test]
    fn test_parse_auto_tab() {
        let s = "time\tvalue\n0.0\t1.0\n1.0\t2.0\n";
        let csv = parse_auto(s).unwrap();
        assert_eq!(csv.headers, vec!["time", "value"]);
        assert_eq!(csv.record_count(), 2);
    }
    #[test]
    fn test_parse_auto_semicolon() {
        let s = "x;y\n1;2\n3;4\n";
        let csv = parse_auto(s).unwrap();
        assert_eq!(csv.headers, vec!["x", "y"]);
        let col = csv.get_column_f64(1).unwrap();
        assert_eq!(col, vec![2.0, 4.0]);
    }
    #[test]
    fn test_read_chunked_basic() {
        let s = "x\n1\n2\n3\n4\n5\n";
        let chunks = read_chunked(s, 2);
        assert_eq!(chunks.len(), 3);
        assert_eq!(chunks[0].record_count(), 2);
        assert_eq!(chunks[1].record_count(), 2);
        assert_eq!(chunks[2].record_count(), 1);
    }
    #[test]
    fn test_read_chunked_exact_multiple() {
        let s = "x\n1\n2\n3\n4\n";
        let chunks = read_chunked(s, 2);
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].record_count(), 2);
        assert_eq!(chunks[1].record_count(), 2);
    }
    #[test]
    fn test_read_chunked_larger_than_data() {
        let s = "x\n1\n2\n";
        let chunks = read_chunked(s, 100);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].record_count(), 2);
    }
    #[test]
    fn test_read_chunked_preserves_headers() {
        let s = "a,b\n1,2\n3,4\n5,6\n";
        let chunks = read_chunked(s, 2);
        for chunk in &chunks {
            assert_eq!(chunk.headers, vec!["a", "b"]);
        }
    }
    #[test]
    fn test_to_string_with_delimiter() {
        let mut csv = CsvFile::new(vec!["a".into(), "b".into()]);
        csv.add_record(vec!["1".into(), "2".into()]);
        let s = csv.to_string_with_delimiter(';');
        assert!(s.starts_with("a;b\n"));
        assert!(s.contains("1;2"));
    }
    #[test]
    fn test_from_str_with_delimiter() {
        let s = "x;y\n1;2\n3;4\n";
        let csv = CsvFile::from_str_with_delimiter(s, ';').unwrap();
        assert_eq!(csv.headers, vec!["x", "y"]);
        assert_eq!(csv.record_count(), 2);
    }
    #[test]
    fn test_get_column_i64() {
        let mut csv = CsvFile::new(vec!["n".into()]);
        csv.add_record(vec!["42".into()]);
        csv.add_record(vec!["-7".into()]);
        let col = csv.get_column_i64(0).unwrap();
        assert_eq!(col, vec![42, -7]);
    }
    #[test]
    fn test_get_column_i64_parse_error() {
        let mut csv = CsvFile::new(vec!["n".into()]);
        csv.add_record(vec!["1.5".into()]);
        assert!(csv.get_column_i64(0).is_err());
    }
    #[test]
    fn test_get_column_strings() {
        let mut csv = CsvFile::new(vec!["name".into()]);
        csv.add_record(vec!["alice".into()]);
        csv.add_record(vec!["bob".into()]);
        let col = csv.get_column_strings(0).unwrap();
        assert_eq!(col, vec!["alice", "bob"]);
    }
    #[test]
    fn test_sort_by_column() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[3.0]);
        csv.add_record_f64(&[1.0]);
        csv.add_record_f64(&[2.0]);
        csv.sort_by_column(0);
        let col = csv.get_column_f64(0).unwrap();
        assert_eq!(col, vec![1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_sort_by_column_already_sorted() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[1.0]);
        csv.add_record_f64(&[2.0]);
        csv.sort_by_column(0);
        let col = csv.get_column_f64(0).unwrap();
        assert_eq!(col, vec![1.0, 2.0]);
    }
    #[test]
    fn test_normalize_header_fn() {
        assert_eq!(normalize_header(" Time Step "), "time_step");
        assert_eq!(normalize_header("X(m/s)"), "x_m_s_");
        assert_eq!(normalize_header("abc_def"), "abc_def");
    }
    #[test]
    fn test_delimiter_roundtrip_semicolon() {
        let mut csv = CsvFile::new(vec!["a".into(), "b".into()]);
        csv.add_record_f64(&[1.0, 2.0]);
        let s = csv.to_string_with_delimiter(';');
        let parsed = CsvFile::from_str_with_delimiter(&s, ';').unwrap();
        assert_eq!(parsed.headers, vec!["a", "b"]);
        let col = parsed.get_column_f64(0).unwrap();
        assert!((col[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_detect_delimiter_empty() {
        assert_eq!(detect_delimiter(""), ',');
    }
    #[test]
    fn test_chunked_empty_input() {
        let chunks = read_chunked("", 5);
        assert!(chunks.is_empty());
    }
    #[test]
    fn test_select_columns_empty_indices() {
        let mut csv = CsvFile::new(vec!["a".into(), "b".into()]);
        csv.add_record(vec!["1".into(), "2".into()]);
        let subset = csv.select_columns(&[]);
        assert!(subset.headers.is_empty());
        assert_eq!(subset.records[0].fields.len(), 0);
    }
    #[test]
    fn schema_validate_ok() {
        let schema = CsvSchema::new(vec![
            ("x".into(), ColumnType::Float),
            ("label".into(), ColumnType::Text),
        ]);
        let mut csv = CsvFile::new(vec!["x".into(), "label".into()]);
        csv.add_record(vec!["3.14".into(), "hello".into()]);
        let errors = schema.validate(&csv);
        assert!(errors.is_empty(), "unexpected errors: {errors:?}");
    }
    #[test]
    fn schema_validate_type_mismatch() {
        let schema = CsvSchema::new(vec![("x".into(), ColumnType::Integer)]);
        let mut csv = CsvFile::new(vec!["x".into()]);
        csv.add_record(vec!["not_an_int".into()]);
        let errors = schema.validate(&csv);
        assert!(!errors.is_empty(), "should report type error");
    }
    #[test]
    fn schema_validate_column_count_mismatch() {
        let schema = CsvSchema::new(vec![
            ("a".into(), ColumnType::Float),
            ("b".into(), ColumnType::Float),
        ]);
        let csv = CsvFile::new(vec!["a".into()]);
        let errors = schema.validate(&csv);
        assert!(!errors.is_empty());
    }
    #[test]
    fn schema_validate_name_mismatch() {
        let schema = CsvSchema::new(vec![("expected".into(), ColumnType::Text)]);
        let csv = CsvFile::new(vec!["actual".into()]);
        let errors = schema.validate(&csv);
        assert!(!errors.is_empty());
    }
    #[test]
    fn time_series_times_extracted() {
        let input = "time,temp\n0.0,300.0\n1.0,301.0\n2.0,302.0\n";
        let ts = TimeSeriesCsv::parse(input, "time").unwrap();
        let times = ts.times().unwrap();
        assert_eq!(times, vec![0.0, 1.0, 2.0]);
    }
    #[test]
    fn time_series_duration() {
        let input = "time,v\n1.0,0.0\n3.0,1.0\n5.0,2.0\n";
        let ts = TimeSeriesCsv::parse(input, "time").unwrap();
        assert!((ts.duration() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn time_series_n_steps() {
        let input = "time,v\n0.0,1.0\n0.5,2.0\n";
        let ts = TimeSeriesCsv::parse(input, "time").unwrap();
        assert_eq!(ts.n_steps(), 2);
    }
    #[test]
    fn time_series_missing_column() {
        let input = "x,y\n1.0,2.0\n";
        let ts = TimeSeriesCsv::parse(input, "time").unwrap();
        assert!(ts.times().is_none());
    }
    fn sample_csv() -> CsvFile {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[1.0]);
        csv.add_record_f64(&[2.0]);
        csv.add_record_f64(&[3.0]);
        csv.add_record_f64(&[4.0]);
        csv
    }
    #[test]
    fn aggregate_sum() {
        let csv = sample_csv();
        let s = aggregate_column(&csv, 0, AggOp::Sum).unwrap();
        assert!((s - 10.0).abs() < 1e-10);
    }
    #[test]
    fn aggregate_mean() {
        let csv = sample_csv();
        let m = aggregate_column(&csv, 0, AggOp::Mean).unwrap();
        assert!((m - 2.5).abs() < 1e-10);
    }
    #[test]
    fn aggregate_min_max() {
        let csv = sample_csv();
        assert!((aggregate_column(&csv, 0, AggOp::Min).unwrap() - 1.0).abs() < 1e-10);
        assert!((aggregate_column(&csv, 0, AggOp::Max).unwrap() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn aggregate_count() {
        let csv = sample_csv();
        assert!((aggregate_column(&csv, 0, AggOp::Count).unwrap() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn aggregate_std() {
        let csv = sample_csv();
        let std = aggregate_column(&csv, 0, AggOp::Std).unwrap();
        assert!((std - 1.25f64.sqrt()).abs() < 1e-9);
    }
    #[test]
    fn aggregate_out_of_range() {
        let csv = sample_csv();
        assert!(aggregate_column(&csv, 99, AggOp::Sum).is_none());
    }
    #[test]
    fn csv_writer_basic() {
        let mut w = CsvWriter::new(vec!["x".into(), "y".into()], ',');
        w.write_row(&["1", "2"]);
        w.write_row(&["3", "4"]);
        let s = w.finish();
        assert!(s.starts_with("x,y\n"));
        assert!(s.contains("1,2"));
    }
    #[test]
    fn csv_writer_f64() {
        let mut w = CsvWriter::new(vec!["val".into()], ',');
        w.write_row_f64(&[2.54321]);
        let s = w.finish();
        assert!(s.contains("2.543210"));
    }
    #[test]
    fn csv_writer_row_count() {
        let mut w = CsvWriter::new(vec!["a".into()], ',');
        for _ in 0..5 {
            w.write_row(&["x"]);
        }
        assert_eq!(w.row_count(), 5);
    }
    #[test]
    fn csv_writer_semicolon_delimiter() {
        let mut w = CsvWriter::new(vec!["a".into(), "b".into()], ';');
        w.write_row(&["1", "2"]);
        let s = w.finish();
        assert!(s.contains("a;b"));
        assert!(s.contains("1;2"));
    }
    #[test]
    fn validation_report_valid() {
        let schema = CsvSchema::new(vec![("x".into(), ColumnType::Float)]);
        let mut csv = CsvFile::new(vec!["x".into()]);
        csv.add_record_f64(&[1.0]);
        let report = validate_csv(&csv, &schema);
        assert!(report.is_valid());
        assert_eq!(report.error_count(), 0);
    }
    #[test]
    fn validation_report_invalid() {
        let schema = CsvSchema::new(vec![("x".into(), ColumnType::Integer)]);
        let mut csv = CsvFile::new(vec!["x".into()]);
        csv.add_record(vec!["hello".into()]);
        let report = validate_csv(&csv, &schema);
        assert!(!report.is_valid());
        assert!(report.error_count() > 0);
    }
    #[test]
    fn lazy_iter_yields_rows() {
        let input = "a,b,c\n1,2,3\n4,5,6\n";
        let mut iter = LazyCsvIter::new(input, ',');
        assert_eq!(iter.headers, vec!["a", "b", "c"]);
        let r1 = iter.next().unwrap();
        assert_eq!(r1, vec!["1", "2", "3"]);
        let r2 = iter.next().unwrap();
        assert_eq!(r2, vec!["4", "5", "6"]);
        assert!(iter.next().is_none());
    }
    #[test]
    fn lazy_iter_empty_input() {
        let mut iter = LazyCsvIter::new("", ',');
        assert!(iter.headers.is_empty());
        assert!(iter.next().is_none());
    }
    #[test]
    fn lazy_iter_header_only() {
        let mut iter = LazyCsvIter::new("x,y\n", ',');
        assert_eq!(iter.headers, vec!["x", "y"]);
        assert!(iter.next().is_none());
    }
    #[test]
    fn lazy_iter_semicolon_delimiter() {
        let input = "a;b\n10;20\n";
        let mut iter = LazyCsvIter::new(input, ';');
        assert_eq!(iter.headers, vec!["a", "b"]);
        let row = iter.next().unwrap();
        assert_eq!(row, vec!["10", "20"]);
    }
}
/// Write a sequence of trajectory frames to a CSV string.
///
/// Format: each frame is written as:
/// ```text
/// # frame_title
/// x,y,z
/// x,y,z
/// ...
/// <blank line>
/// ```
pub fn write_trajectory_csv(frames: &[TrajectoryFrame]) -> String {
    let mut out = String::new();
    for frame in frames {
        if !frame.title.is_empty() {
            out.push_str(&format!("# {}\n", frame.title));
        }
        for pos in &frame.positions {
            out.push_str(&format!("{},{},{}\n", pos[0], pos[1], pos[2]));
        }
        out.push('\n');
    }
    out
}
/// Read trajectory frames from a CSV string.
///
/// Each block of non-blank lines is one frame; lines starting with `#` are
/// treated as the frame title.  Each data line must have exactly 3
/// comma-separated values (`x,y,z`). Returns an error string with the line
/// number on failure.
pub fn read_trajectory_csv(s: &str) -> std::result::Result<Vec<TrajectoryFrame>, String> {
    let mut frames: Vec<TrajectoryFrame> = Vec::new();
    let mut current = TrajectoryFrame::new();
    let mut in_frame = false;
    for (line_no, raw_line) in s.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() {
            if in_frame {
                frames.push(std::mem::take(&mut current));
                in_frame = false;
            }
            continue;
        }
        if line.starts_with('#') {
            current.title = line.trim_start_matches('#').trim().to_string();
            in_frame = true;
            continue;
        }
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() != 3 {
            return Err(format!(
                "line {}: expected 3 comma-separated values, got {}",
                line_no + 1,
                parts.len()
            ));
        }
        let x: f64 = parts[0]
            .trim()
            .parse()
            .map_err(|e| format!("line {}: x parse error: {}", line_no + 1, e))?;
        let y: f64 = parts[1]
            .trim()
            .parse()
            .map_err(|e| format!("line {}: y parse error: {}", line_no + 1, e))?;
        let z: f64 = parts[2]
            .trim()
            .parse()
            .map_err(|e| format!("line {}: z parse error: {}", line_no + 1, e))?;
        current.positions.push([x, y, z]);
        in_frame = true;
    }
    if in_frame {
        frames.push(current);
    }
    Ok(frames)
}
#[cfg(test)]
mod tests_dataframe {
    use super::*;
    use crate::csv_io::types::*;
    #[test]
    fn dataframe_from_csv_types() {
        let csv_str = "id,x,label\n1,3.14,hello\n2,2.71,world\n";
        let df = csv_str.parse::<CsvDataFrame>().unwrap();
        assert_eq!(df.n_cols(), 3);
        assert_eq!(df.n_rows(), 2);
        assert_eq!(df.column(0).unwrap().column_type(), ColumnType::Integer);
        assert_eq!(df.column(1).unwrap().column_type(), ColumnType::Float);
        assert_eq!(df.column(2).unwrap().column_type(), ColumnType::Text);
    }
    #[test]
    fn dataframe_float_column_by_name() {
        let csv_str = "t,v\n0.0,1.5\n1.0,2.5\n";
        let df = csv_str.parse::<CsvDataFrame>().unwrap();
        let v = df.float_column("v").unwrap();
        assert_eq!(v.len(), 2);
        assert!((v[0] - 1.5).abs() < 1e-12);
        assert!((v[1] - 2.5).abs() < 1e-12);
    }
    #[test]
    fn dataframe_integer_column_by_name() {
        let csv_str = "n,label\n10,a\n20,b\n";
        let df = csv_str.parse::<CsvDataFrame>().unwrap();
        let n = df.integer_column("n").unwrap();
        assert_eq!(n, &vec![10_i64, 20_i64]);
    }
    #[test]
    fn dataframe_text_column_by_name() {
        let csv_str = "name,val\nalice,1.0\nbob,2.0\n";
        let df = csv_str.parse::<CsvDataFrame>().unwrap();
        let names = df.text_column("name").unwrap();
        assert_eq!(names, &vec!["alice".to_string(), "bob".to_string()]);
    }
    #[test]
    fn dataframe_column_index_missing() {
        let csv_str = "a,b\n1,2\n";
        let df = csv_str.parse::<CsvDataFrame>().unwrap();
        assert!(df.column_index("nope").is_none());
    }
    #[test]
    fn dataframe_to_csv_string_roundtrip() {
        let csv_str = "x,y\n1.5,2.5\n3.5,4.5\n";
        let df = csv_str.parse::<CsvDataFrame>().unwrap();
        let out = df.to_csv_string();
        assert!(out.contains("x,y"));
        let df2 = out.parse::<CsvDataFrame>().unwrap();
        let x = df2.float_column("x").unwrap();
        assert!((x[0] - 1.5).abs() < 1e-12);
        assert!((x[1] - 3.5).abs() < 1e-12);
    }
    #[test]
    fn dataframe_empty_input() {
        assert!("".parse::<CsvDataFrame>().is_err());
    }
    #[test]
    fn dataframe_n_rows_n_cols() {
        let csv_str = "a,b,c\n1,2,3\n4,5,6\n7,8,9\n";
        let df = csv_str.parse::<CsvDataFrame>().unwrap();
        assert_eq!(df.n_rows(), 3);
        assert_eq!(df.n_cols(), 3);
    }
    #[test]
    fn dataframe_column_by_name_returns_none_for_missing() {
        let csv_str = "x\n1.0\n";
        let df = csv_str.parse::<CsvDataFrame>().unwrap();
        assert!(df.column_by_name("missing").is_none());
    }
    #[test]
    fn dataframe_with_delimiter() {
        let csv_str = "x;y\n1.0;2.0\n3.0;4.0\n";
        let df = CsvDataFrame::from_str_with_delimiter(csv_str, ';').unwrap();
        assert_eq!(df.n_cols(), 2);
        let y = df.float_column("y").unwrap();
        assert!((y[0] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn trajectory_write_read_roundtrip() {
        let frames = vec![
            TrajectoryFrame {
                title: "frame 0".to_string(),
                positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            },
            TrajectoryFrame {
                title: "frame 1".to_string(),
                positions: vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
            },
        ];
        let csv = write_trajectory_csv(&frames);
        let parsed = read_trajectory_csv(&csv).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].title, "frame 0");
        assert_eq!(parsed[0].n_atoms(), 2);
        assert!((parsed[0].positions[0][0] - 1.0).abs() < 1e-12);
        assert!((parsed[1].positions[1][2] - 0.6).abs() < 1e-12);
    }
    #[test]
    fn trajectory_single_frame_no_title() {
        let csv = "1.0,2.0,3.0\n4.0,5.0,6.0\n";
        let frames = read_trajectory_csv(csv).unwrap();
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].n_atoms(), 2);
        assert!(frames[0].title.is_empty());
    }
    #[test]
    fn trajectory_empty_input() {
        let frames = read_trajectory_csv("").unwrap();
        assert!(frames.is_empty());
    }
    #[test]
    fn trajectory_error_on_bad_value() {
        let csv = "1.0,not_a_float,3.0\n";
        assert!(read_trajectory_csv(csv).is_err());
    }
    #[test]
    fn trajectory_error_on_wrong_column_count() {
        let csv = "1.0,2.0\n";
        assert!(read_trajectory_csv(csv).is_err());
    }
    #[test]
    fn trajectory_multiple_frames_position_accuracy() {
        let frames = vec![TrajectoryFrame {
            title: String::new(),
            positions: vec![[0.123456789, -0.987654321, 1.111111111]],
        }];
        let csv = write_trajectory_csv(&frames);
        let parsed = read_trajectory_csv(&csv).unwrap();
        assert!((parsed[0].positions[0][0] - 0.123456789).abs() < 1e-9);
        assert!((parsed[0].positions[0][1] - (-0.987654321)).abs() < 1e-9);
        assert!((parsed[0].positions[0][2] - 1.111111111).abs() < 1e-9);
    }
    #[test]
    fn trajectory_frame_n_atoms() {
        let f = TrajectoryFrame {
            title: "t".into(),
            positions: vec![[0.0; 3]; 5],
        };
        assert_eq!(f.n_atoms(), 5);
    }
    #[test]
    fn csv_column_data_len_and_type() {
        let col = CsvColumnData::Float(vec![1.0, 2.0, 3.0]);
        assert_eq!(col.len(), 3);
        assert!(!col.is_empty());
        assert_eq!(col.column_type(), ColumnType::Float);
        let int_col = CsvColumnData::Integer(vec![10, 20]);
        assert_eq!(int_col.len(), 2);
        assert_eq!(int_col.column_type(), ColumnType::Integer);
        let text_col = CsvColumnData::Text(vec!["a".into()]);
        assert_eq!(text_col.column_type(), ColumnType::Text);
        assert!(!text_col.is_empty());
    }
}
/// Merge two [`CsvFile`] objects that share the same column schema.
///
/// Rows from `other` are appended after the rows of `base`.
/// Returns `Err` if the header lists differ.
pub fn merge_csv_files(base: &CsvFile, other: &CsvFile) -> Result<CsvFile, String> {
    if base.headers != other.headers {
        return Err(format!(
            "header mismatch: {:?} vs {:?}",
            base.headers, other.headers
        ));
    }
    let mut result = CsvFile::new(base.headers.clone());
    for rec in &base.records {
        result.records.push(CsvRecord {
            fields: rec.fields.clone(),
        });
    }
    for rec in &other.records {
        result.records.push(CsvRecord {
            fields: rec.fields.clone(),
        });
    }
    Ok(result)
}
/// Transpose a [`CsvFile`]: rows become columns and columns become rows.
///
/// The resulting file has column headers `"col_0"`, `"col_1"`, …
/// (original row indices), and each new row corresponds to an original column.
pub fn transpose_csv(csv: &CsvFile) -> CsvFile {
    let n_rows = csv.records.len();
    let n_cols = csv.headers.len();
    if n_rows == 0 || n_cols == 0 {
        return CsvFile::new(vec![]);
    }
    let new_headers: Vec<String> = (0..n_rows).map(|i| format!("col_{}", i)).collect();
    let mut result = CsvFile::new(new_headers);
    for col in 0..n_cols {
        let fields: Vec<String> = (0..n_rows)
            .map(|row| {
                csv.records[row]
                    .fields
                    .get(col)
                    .cloned()
                    .unwrap_or_default()
            })
            .collect();
        result.records.push(CsvRecord { fields });
    }
    result
}
/// Inner-join two [`CsvFile`]s on a shared key column.
///
/// Rows are matched when the key column value is equal (string comparison).
/// The output contains all columns from `left` followed by the non-key columns
/// of `right`.
pub fn inner_join_csv(left: &CsvFile, right: &CsvFile, key: &str) -> Result<CsvFile, String> {
    let left_key_idx = left
        .get_column_by_name(key)
        .ok_or_else(|| format!("key '{}' not in left file", key))?;
    let right_key_idx = right
        .get_column_by_name(key)
        .ok_or_else(|| format!("key '{}' not in right file", key))?;
    let mut headers = left.headers.clone();
    for (i, h) in right.headers.iter().enumerate() {
        if i != right_key_idx {
            headers.push(h.clone());
        }
    }
    let mut result = CsvFile::new(headers);
    for l_rec in &left.records {
        let l_key = l_rec.fields.get(left_key_idx).cloned().unwrap_or_default();
        for r_rec in &right.records {
            let r_key = r_rec.fields.get(right_key_idx).cloned().unwrap_or_default();
            if l_key == r_key {
                let mut fields = l_rec.fields.clone();
                for (i, f) in r_rec.fields.iter().enumerate() {
                    if i != right_key_idx {
                        fields.push(f.clone());
                    }
                }
                result.records.push(CsvRecord { fields });
            }
        }
    }
    Ok(result)
}
/// Compute the per-cell numeric difference between two [`CsvFile`]s.
///
/// Both files must have the same dimensions and numeric-only columns.
/// Returns `Err` on dimension mismatch or parse failure.
pub fn diff_csv(a: &CsvFile, b: &CsvFile) -> Result<CsvFile, String> {
    if a.headers != b.headers {
        return Err("header mismatch".to_string());
    }
    if a.records.len() != b.records.len() {
        return Err(format!(
            "row count mismatch: {} vs {}",
            a.records.len(),
            b.records.len()
        ));
    }
    let mut result = CsvFile::new(a.headers.clone());
    for (i, (ar, br)) in a.records.iter().zip(b.records.iter()).enumerate() {
        let mut fields = Vec::with_capacity(a.headers.len());
        for (j, (af, bf)) in ar.fields.iter().zip(br.fields.iter()).enumerate() {
            let av: f64 = af
                .parse()
                .map_err(|_| format!("row {} col {}: not numeric", i, j))?;
            let bv: f64 = bf
                .parse()
                .map_err(|_| format!("row {} col {}: not numeric", i, j))?;
            fields.push(format!("{}", av - bv));
        }
        result.records.push(CsvRecord { fields });
    }
    Ok(result)
}
/// Detect whether the first row of a CSV string is likely a header row.
///
/// A row is considered a header if it contains at least one non-numeric,
/// non-empty token (i.e., any field that cannot be parsed as `f64`).
pub fn has_header(s: &str, delim: char) -> bool {
    let first = s.lines().next().unwrap_or("");
    first.split(delim).any(|f| {
        let t = f.trim();
        !t.is_empty() && t.parse::<f64>().is_err()
    })
}
/// Infer headers for a headerless CSV.
///
/// Returns `"col_0"`, `"col_1"`, … based on the number of fields in the
/// first data row.
pub fn infer_headers(s: &str, delim: char) -> Vec<String> {
    let n = s
        .lines()
        .next()
        .map(|l| l.split(delim).count())
        .unwrap_or(0);
    (0..n).map(|i| format!("col_{}", i)).collect()
}
/// Parse a CSV string that may or may not have a header row.
///
/// If `auto_header` is `true` (default behaviour), header presence is
/// inferred using [`has_header`].  If no header is detected, synthetic
/// column names `col_0`, `col_1`, … are used.
pub fn parse_smart(s: &str) -> Result<CsvFile, String> {
    let delim = detect_delimiter(s);
    if has_header(s, delim) {
        CsvFile::from_str_with_delimiter(s, delim)
    } else {
        let headers = infer_headers(s, delim);
        let mut result = CsvFile::new(headers);
        for line in s.lines() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<String> = line.split(delim).map(str::trim).map(String::from).collect();
            result.add_record(fields);
        }
        Ok(result)
    }
}
/// Infer whether a column can be interpreted as boolean.
///
/// Recognised truthy values: `"true"`, `"1"`, `"yes"`, `"on"` (case-insensitive).
/// Recognised falsy values: `"false"`, `"0"`, `"no"`, `"off"`.
/// Any other value → not boolean.
pub fn is_boolean_column(csv: &CsvFile, col_idx: usize) -> bool {
    let truthy = ["true", "1", "yes", "on"];
    let falsy = ["false", "0", "no", "off"];
    if col_idx >= csv.headers.len() {
        return false;
    }
    for rec in &csv.records {
        if let Some(v) = rec.fields.get(col_idx) {
            let lower = v.trim().to_lowercase();
            if !truthy.contains(&lower.as_str()) && !falsy.contains(&lower.as_str()) {
                return false;
            }
        }
    }
    true
}
/// Parse a boolean column as `Vec`bool`.
///
/// Returns `Err` if any value is not recognisable as boolean.
pub fn get_column_bool(csv: &CsvFile, col_idx: usize) -> Result<Vec<bool>, String> {
    let truthy = ["true", "1", "yes", "on"];
    let falsy = ["false", "0", "no", "off"];
    if col_idx >= csv.headers.len() {
        return Err(format!("column index {} out of range", col_idx));
    }
    let mut out = Vec::with_capacity(csv.records.len());
    for (row, rec) in csv.records.iter().enumerate() {
        let raw = rec
            .fields
            .get(col_idx)
            .ok_or_else(|| format!("row {} has no field at col {}", row, col_idx))?;
        let lower = raw.trim().to_lowercase();
        if truthy.contains(&lower.as_str()) {
            out.push(true);
        } else if falsy.contains(&lower.as_str()) {
            out.push(false);
        } else {
            return Err(format!("row {}: '{}' is not a boolean value", row, raw));
        }
    }
    Ok(out)
}
/// Sample every `n`-th row from a [`CsvFile`].
///
/// The header is preserved. `stride = 1` returns all rows unchanged.
pub fn sample_every_nth(csv: &CsvFile, stride: usize) -> CsvFile {
    if stride == 0 {
        return CsvFile::new(csv.headers.clone());
    }
    let mut result = CsvFile::new(csv.headers.clone());
    for (i, rec) in csv.records.iter().enumerate() {
        if i % stride == 0 {
            result.records.push(CsvRecord {
                fields: rec.fields.clone(),
            });
        }
    }
    result
}
/// Deduplicate rows in a [`CsvFile`] based on a key column (keep first occurrence).
pub fn dedup_by_column(csv: &CsvFile, col_idx: usize) -> CsvFile {
    use std::collections::HashSet;
    let mut seen: HashSet<String> = HashSet::new();
    let mut result = CsvFile::new(csv.headers.clone());
    for rec in &csv.records {
        let key = rec.fields.get(col_idx).cloned().unwrap_or_default();
        if seen.insert(key) {
            result.records.push(CsvRecord {
                fields: rec.fields.clone(),
            });
        }
    }
    result
}
#[cfg(test)]
mod tests_csv_new {
    use super::*;
    use crate::csv_io::types::*;
    #[test]
    fn streaming_reader_headers() {
        let input = "time,x,y\n0.0,1.0,2.0\n1.0,3.0,4.0\n";
        let r = StreamingCsvReader::new(input, ',');
        assert_eq!(r.headers, vec!["time", "x", "y"]);
        assert_eq!(r.n_cols(), 3);
    }
    #[test]
    fn streaming_reader_next_row() {
        let input = "a,b\n1,2\n3,4\n";
        let mut r = StreamingCsvReader::new(input, ',');
        let row = r.next_row().unwrap();
        assert_eq!(row, vec!["1", "2"]);
        assert_eq!(r.current_row(), 1);
        let row2 = r.next_row().unwrap();
        assert_eq!(row2, vec!["3", "4"]);
    }
    #[test]
    fn streaming_reader_eof() {
        let input = "a\n1\n";
        let mut r = StreamingCsvReader::new(input, ',');
        r.next_row();
        assert!(r.next_row().is_none());
    }
    #[test]
    fn streaming_reader_skips_blank_lines() {
        let input = "a\n1\n\n2\n";
        let mut r = StreamingCsvReader::new(input, ',');
        r.next_row();
        let row = r.next_row().unwrap();
        assert_eq!(row, vec!["2"]);
    }
    #[test]
    fn streaming_reader_collect_all() {
        let input = "x,y\n1,2\n3,4\n5,6\n";
        let r = StreamingCsvReader::new(input, ',');
        let csv = r.collect_all();
        assert_eq!(csv.record_count(), 3);
        assert_eq!(csv.headers, vec!["x", "y"]);
    }
    #[test]
    fn streaming_reader_auto_delimiter() {
        let input = "a\tb\tc\n1\t2\t3\n";
        let r = StreamingCsvReader::auto(input);
        assert_eq!(r.delimiter, '\t');
        assert_eq!(r.headers, vec!["a", "b", "c"]);
    }
    #[test]
    fn merge_csv_files_basic() {
        let mut a = CsvFile::new(vec!["x".into()]);
        a.add_record_f64(&[1.0]);
        let mut b = CsvFile::new(vec!["x".into()]);
        b.add_record_f64(&[2.0]);
        b.add_record_f64(&[3.0]);
        let merged = merge_csv_files(&a, &b).unwrap();
        assert_eq!(merged.record_count(), 3);
    }
    #[test]
    fn merge_csv_files_header_mismatch() {
        let a = CsvFile::new(vec!["x".into()]);
        let b = CsvFile::new(vec!["y".into()]);
        assert!(merge_csv_files(&a, &b).is_err());
    }
    #[test]
    fn transpose_basic() {
        let mut csv = CsvFile::new(vec!["a".into(), "b".into()]);
        csv.add_record(vec!["1".into(), "2".into()]);
        csv.add_record(vec!["3".into(), "4".into()]);
        let t = transpose_csv(&csv);
        assert_eq!(t.headers.len(), 2);
        assert_eq!(t.record_count(), 2);
        assert_eq!(t.records[0].fields, vec!["1", "3"]);
        assert_eq!(t.records[1].fields, vec!["2", "4"]);
    }
    #[test]
    fn transpose_empty() {
        let csv = CsvFile::new(vec![]);
        let t = transpose_csv(&csv);
        assert!(t.headers.is_empty());
    }
    #[test]
    fn inner_join_basic() {
        let mut left = CsvFile::new(vec!["id".into(), "name".into()]);
        left.add_record(vec!["1".into(), "alice".into()]);
        left.add_record(vec!["2".into(), "bob".into()]);
        let mut right = CsvFile::new(vec!["id".into(), "score".into()]);
        right.add_record(vec!["1".into(), "90".into()]);
        right.add_record(vec!["3".into(), "80".into()]);
        let joined = inner_join_csv(&left, &right, "id").unwrap();
        assert_eq!(joined.record_count(), 1);
        assert_eq!(joined.headers, vec!["id", "name", "score"]);
        assert_eq!(joined.records[0].fields[1], "alice");
        assert_eq!(joined.records[0].fields[2], "90");
    }
    #[test]
    fn inner_join_missing_key() {
        let left = CsvFile::new(vec!["a".into()]);
        let right = CsvFile::new(vec!["b".into()]);
        assert!(inner_join_csv(&left, &right, "id").is_err());
    }
    #[test]
    fn diff_csv_basic() {
        let mut a = CsvFile::new(vec!["v".into()]);
        a.add_record(vec!["5.0".into()]);
        a.add_record(vec!["3.0".into()]);
        let mut b = CsvFile::new(vec!["v".into()]);
        b.add_record(vec!["1.0".into()]);
        b.add_record(vec!["1.0".into()]);
        let d = diff_csv(&a, &b).unwrap();
        let vals = d.get_column_f64(0).unwrap();
        assert!((vals[0] - 4.0).abs() < 1e-12);
        assert!((vals[1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn diff_csv_row_count_mismatch() {
        let mut a = CsvFile::new(vec!["v".into()]);
        a.add_record(vec!["1.0".into()]);
        let b = CsvFile::new(vec!["v".into()]);
        assert!(diff_csv(&a, &b).is_err());
    }
    #[test]
    fn has_header_true() {
        assert!(has_header("time,x,y\n0,1,2\n", ','));
    }
    #[test]
    fn has_header_false_all_numbers() {
        assert!(!has_header("0,1,2\n3,4,5\n", ','));
    }
    #[test]
    fn infer_headers_count() {
        let headers = infer_headers("1,2,3,4\n", ',');
        assert_eq!(headers, vec!["col_0", "col_1", "col_2", "col_3"]);
    }
    #[test]
    fn parse_smart_with_header() {
        let s = "a,b\n1,2\n3,4\n";
        let csv = parse_smart(s).unwrap();
        assert_eq!(csv.headers, vec!["a", "b"]);
        assert_eq!(csv.record_count(), 2);
    }
    #[test]
    fn parse_smart_without_header() {
        let s = "1,2\n3,4\n";
        let csv = parse_smart(s).unwrap();
        assert_eq!(csv.headers, vec!["col_0", "col_1"]);
        assert_eq!(csv.record_count(), 2);
    }
    #[test]
    fn boolean_column_detection() {
        let mut csv = CsvFile::new(vec!["flag".into()]);
        csv.add_record(vec!["true".into()]);
        csv.add_record(vec!["false".into()]);
        csv.add_record(vec!["yes".into()]);
        assert!(is_boolean_column(&csv, 0));
    }
    #[test]
    fn boolean_column_rejection() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record(vec!["true".into()]);
        csv.add_record(vec!["maybe".into()]);
        assert!(!is_boolean_column(&csv, 0));
    }
    #[test]
    fn get_column_bool_values() {
        let mut csv = CsvFile::new(vec!["b".into()]);
        csv.add_record(vec!["1".into()]);
        csv.add_record(vec!["0".into()]);
        csv.add_record(vec!["yes".into()]);
        csv.add_record(vec!["no".into()]);
        let vals = get_column_bool(&csv, 0).unwrap();
        assert_eq!(vals, vec![true, false, true, false]);
    }
    #[test]
    fn get_column_bool_error_on_bad() {
        let mut csv = CsvFile::new(vec!["b".into()]);
        csv.add_record(vec!["maybe".into()]);
        assert!(get_column_bool(&csv, 0).is_err());
    }
    #[test]
    fn sample_every_nth_basic() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        for i in 0..10_usize {
            csv.add_record(vec![i.to_string()]);
        }
        let sampled = sample_every_nth(&csv, 3);
        assert_eq!(sampled.record_count(), 4);
    }
    #[test]
    fn sample_every_nth_stride_one() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[1.0]);
        csv.add_record_f64(&[2.0]);
        let s = sample_every_nth(&csv, 1);
        assert_eq!(s.record_count(), 2);
    }
    #[test]
    fn sample_every_nth_zero_stride() {
        let mut csv = CsvFile::new(vec!["v".into()]);
        csv.add_record_f64(&[1.0]);
        let s = sample_every_nth(&csv, 0);
        assert_eq!(s.record_count(), 0);
    }
    #[test]
    fn dedup_by_column_basic() {
        let mut csv = CsvFile::new(vec!["id".into(), "val".into()]);
        csv.add_record(vec!["1".into(), "a".into()]);
        csv.add_record(vec!["2".into(), "b".into()]);
        csv.add_record(vec!["1".into(), "c".into()]);
        let deduped = dedup_by_column(&csv, 0);
        assert_eq!(deduped.record_count(), 2);
        assert_eq!(deduped.records[0].fields[1], "a");
    }
    #[test]
    fn dedup_by_column_all_unique() {
        let mut csv = CsvFile::new(vec!["id".into()]);
        for i in 0..5_usize {
            csv.add_record(vec![i.to_string()]);
        }
        let d = dedup_by_column(&csv, 0);
        assert_eq!(d.record_count(), 5);
    }
    #[test]
    fn lazy_iter_yields_correct_rows() {
        let input = "x,y\n1,2\n3,4\n5,6\n";
        let mut it = LazyCsvIter::new(input, ',');
        assert_eq!(it.headers, vec!["x", "y"]);
        let r1 = it.next().unwrap();
        assert_eq!(r1, vec!["1", "2"]);
        let r2 = it.next().unwrap();
        assert_eq!(r2, vec!["3", "4"]);
        let r3 = it.next().unwrap();
        assert_eq!(r3, vec!["5", "6"]);
        assert!(it.next().is_none());
    }
    #[test]
    fn lazy_iter_tab_delimiter() {
        let input = "a\tb\n10\t20\n";
        let mut it = LazyCsvIter::new(input, '\t');
        assert_eq!(it.headers, vec!["a", "b"]);
        let row = it.next().unwrap();
        assert_eq!(row[0], "10");
        assert_eq!(row[1], "20");
    }
    #[test]
    fn validation_report_is_valid() {
        let schema = CsvSchema::new(vec![("x".into(), ColumnType::Float)]);
        let mut csv = CsvFile::new(vec!["x".into()]);
        csv.add_record(vec!["3.14".into()]);
        let report = validate_csv(&csv, &schema);
        assert!(report.is_valid());
        assert_eq!(report.error_count(), 0);
    }
    #[test]
    fn validation_report_has_errors() {
        let schema = CsvSchema::new(vec![("x".into(), ColumnType::Integer)]);
        let mut csv = CsvFile::new(vec!["x".into()]);
        csv.add_record(vec!["not_int".into()]);
        let report = validate_csv(&csv, &schema);
        assert!(!report.is_valid());
        assert!(report.error_count() > 0);
    }
    #[test]
    fn normalize_header_special_chars() {
        assert_eq!(normalize_header("E (J/mol)"), "e__j_mol_");
    }
    #[test]
    fn normalize_header_already_clean() {
        assert_eq!(normalize_header("velocity_x"), "velocity_x");
    }
}
