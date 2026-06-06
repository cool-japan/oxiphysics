//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::Error;
use std::collections::HashMap;

use super::types::{ColumnType, CsvChangedRow, CsvDiff, CsvTable, PivotAgg};

/// Split a CSV line by the given delimiter, respecting double-quoted fields.
pub(super) fn split_csv_line(line: &str, delim: char) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch == '"' {
            if in_quotes {
                if chars.peek() == Some(&'"') {
                    chars.next();
                    current.push('"');
                } else {
                    in_quotes = false;
                }
            } else {
                in_quotes = true;
            }
        } else if ch == delim && !in_quotes {
            fields.push(current.clone());
            current.clear();
        } else {
            current.push(ch);
        }
    }
    fields.push(current);
    fields
}
/// Wrap a field in double-quotes if it contains the delimiter, a quote, or a newline.
pub(super) fn quote_field(field: &str, delimiter: char) -> String {
    if field.contains(delimiter)
        || field.contains('"')
        || field.contains('\n')
        || field.contains('\r')
    {
        let escaped = field.replace('"', "\"\"");
        format!("\"{}\"", escaped)
    } else {
        field.to_string()
    }
}
/// Infer the column type for a slice of string values.
///
/// Empty strings are treated as missing values and do not affect type inference.
pub fn infer_column_type(values: &[&str]) -> ColumnType {
    let non_empty: Vec<&str> = values
        .iter()
        .copied()
        .filter(|s| !s.trim().is_empty())
        .collect();
    if non_empty.is_empty() {
        return ColumnType::Empty;
    }
    let all_int = non_empty.iter().all(|s| s.trim().parse::<i64>().is_ok());
    if all_int {
        return ColumnType::Integer;
    }
    let all_float = non_empty.iter().all(|s| s.trim().parse::<f64>().is_ok());
    if all_float {
        return ColumnType::Float;
    }
    let bool_values = ["true", "false", "yes", "no", "1", "0"];
    let all_bool = non_empty
        .iter()
        .all(|s| bool_values.contains(&s.trim().to_lowercase().as_str()));
    if all_bool {
        return ColumnType::Boolean;
    }
    ColumnType::Text
}
/// Infer types for all columns in a [`CsvTable`].
///
/// Returns a `Vec<(column_name, ColumnType)>` in column order.
pub fn infer_table_types(table: &CsvTable) -> Vec<(String, ColumnType)> {
    table
        .headers
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let values: Vec<&str> = table.rows.iter().map(|r| r[i].as_str()).collect();
            let typ = infer_column_type(&values);
            (name.clone(), typ)
        })
        .collect()
}
/// Append all rows of `right` to `left`, returning a new table.
///
/// Both tables must have identical headers (in the same order).
///
/// # Example
///
/// ```no_run
/// use oxiphysics_io::csv::{CsvTable, csv_merge};
///
/// let a = CsvTable::parse("x,y\n1,2\n", ',').unwrap();
/// let b = CsvTable::parse("x,y\n3,4\n", ',').unwrap();
/// let merged = csv_merge(&a, &b).unwrap();
/// assert_eq!(merged.row_count(), 2);
/// ```
pub fn csv_merge(left: &CsvTable, right: &CsvTable) -> std::result::Result<CsvTable, Error> {
    if left.headers != right.headers {
        return Err(Error::Parse(format!(
            "csv_merge: header mismatch: {:?} vs {:?}",
            left.headers, right.headers
        )));
    }
    let mut result = left.clone();
    result.rows.extend(right.rows.iter().cloned());
    Ok(result)
}
/// Perform an inner join of `left` and `right` on the named key column.
///
/// The key column must exist in both tables. Duplicate key values in `right`
/// are matched to every occurrence in `left` (standard SQL inner join semantics).
/// Columns from the right table (except the key) are appended to the result.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_io::csv::{CsvTable, csv_join};
///
/// let left  = CsvTable::parse("id,x\n1,10\n2,20\n3,30\n", ',').unwrap();
/// let right = CsvTable::parse("id,y\n1,100\n3,300\n", ',').unwrap();
/// let joined = csv_join(&left, &right, "id").unwrap();
/// assert_eq!(joined.row_count(), 2);
/// let x_col = joined.column_f64("x").unwrap();
/// assert!((x_col[0] - 10.0).abs() < 1e-10);
/// assert!((x_col[1] - 30.0).abs() < 1e-10);
/// ```
pub fn csv_join(
    left: &CsvTable,
    right: &CsvTable,
    key: &str,
) -> std::result::Result<CsvTable, Error> {
    let left_key_idx = left.column_index(key)?;
    let right_key_idx = right.column_index(key)?;
    let mut right_map: HashMap<String, Vec<Vec<String>>> = HashMap::new();
    for row in &right.rows {
        right_map
            .entry(row[right_key_idx].clone())
            .or_default()
            .push(row.clone());
    }
    let right_extra_headers: Vec<String> = right
        .headers
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != right_key_idx)
        .map(|(_, h)| h.clone())
        .collect();
    let right_extra_indices: Vec<usize> = right
        .headers
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != right_key_idx)
        .map(|(i, _)| i)
        .collect();
    let mut result_headers = left.headers.clone();
    result_headers.extend(right_extra_headers);
    let mut result_rows = Vec::new();
    for left_row in &left.rows {
        let key_val = &left_row[left_key_idx];
        if let Some(right_rows) = right_map.get(key_val) {
            for right_row in right_rows {
                let mut merged = left_row.clone();
                for &ri in &right_extra_indices {
                    merged.push(right_row[ri].clone());
                }
                result_rows.push(merged);
            }
        }
    }
    Ok(CsvTable {
        headers: result_headers,
        rows: result_rows,
    })
}
/// Create a pivot table from `table`.
///
/// - `row_col`: column whose distinct values become pivot rows (index).
/// - `col_col`: column whose distinct values become pivot columns.
/// - `val_col`: numeric column to aggregate.
/// - `agg`: aggregation function.
///
/// Returns a new [`CsvTable`] where the first column is the row index
/// and subsequent columns are the pivot column values.
///
/// Missing combinations produce an empty cell.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_io::csv::{CsvTable, csv_pivot, PivotAgg};
///
/// let data = "region,product,sales\nNorth,A,10\nNorth,B,20\nSouth,A,30\nSouth,B,40\n";
/// let table = CsvTable::parse(data, ',').unwrap();
/// let pivot = csv_pivot(&table, "region", "product", "sales", PivotAgg::Sum).unwrap();
/// assert!(pivot.headers.contains(&"A".to_string()));
/// assert!(pivot.headers.contains(&"B".to_string()));
/// ```
pub fn csv_pivot(
    table: &CsvTable,
    row_col: &str,
    col_col: &str,
    val_col: &str,
    agg: PivotAgg,
) -> std::result::Result<CsvTable, Error> {
    let row_idx = table.column_index(row_col)?;
    let col_idx = table.column_index(col_col)?;
    let val_idx = table.column_index(val_col)?;
    let mut row_keys: Vec<String> = Vec::new();
    let mut col_keys: Vec<String> = Vec::new();
    for row in &table.rows {
        let rk = row[row_idx].clone();
        if !row_keys.contains(&rk) {
            row_keys.push(rk);
        }
        let ck = row[col_idx].clone();
        if !col_keys.contains(&ck) {
            col_keys.push(ck);
        }
    }
    let mut buckets: HashMap<(String, String), Vec<f64>> = HashMap::new();
    for row in &table.rows {
        let rk = row[row_idx].clone();
        let ck = row[col_idx].clone();
        let vs = row[val_idx].trim();
        if let Ok(v) = vs.parse::<f64>() {
            buckets.entry((rk, ck)).or_default().push(v);
        }
    }
    let mut headers = vec![row_col.to_string()];
    headers.extend(col_keys.iter().cloned());
    let rows: Vec<Vec<String>> = row_keys
        .iter()
        .map(|rk| {
            let mut row_out = vec![rk.clone()];
            for ck in &col_keys {
                let cell = if let Some(vals) = buckets.get(&(rk.clone(), ck.clone())) {
                    let agg_val = match agg {
                        PivotAgg::Sum => vals.iter().sum::<f64>(),
                        PivotAgg::Mean => vals.iter().sum::<f64>() / vals.len() as f64,
                        PivotAgg::Count => vals.len() as f64,
                        PivotAgg::Min => vals.iter().cloned().fold(f64::INFINITY, f64::min),
                        PivotAgg::Max => vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
                    };
                    format!("{}", agg_val)
                } else {
                    String::new()
                };
                row_out.push(cell);
            }
            row_out
        })
        .collect();
    Ok(CsvTable { headers, rows })
}
/// Compute the diff between two CSV tables keyed by `key_col`.
///
/// Rows are matched by the value of `key_col`. For each key:
/// - If it exists only in `left` → `removed`.
/// - If it exists only in `right` → `added`.
/// - If it exists in both but the rows differ → `changed`.
///
/// Both tables must have the same headers.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_io::csv::{CsvTable, csv_diff};
///
/// let a = CsvTable::parse("id,v\n1,10\n2,20\n3,30\n", ',').unwrap();
/// let b = CsvTable::parse("id,v\n1,10\n2,99\n4,40\n", ',').unwrap();
/// let diff = csv_diff(&a, &b, "id").unwrap();
/// assert_eq!(diff.removed.len(), 1); // id=3
/// assert_eq!(diff.added.len(), 1);   // id=4
/// assert_eq!(diff.changed.len(), 1); // id=2 changed v
/// ```
pub fn csv_diff(
    left: &CsvTable,
    right: &CsvTable,
    key_col: &str,
) -> std::result::Result<CsvDiff, Error> {
    if left.headers != right.headers {
        return Err(Error::Parse(format!(
            "csv_diff: header mismatch: {:?} vs {:?}",
            left.headers, right.headers
        )));
    }
    let key_idx = left.column_index(key_col)?;
    let left_map: HashMap<String, &Vec<String>> =
        left.rows.iter().map(|r| (r[key_idx].clone(), r)).collect();
    let right_map: HashMap<String, &Vec<String>> =
        right.rows.iter().map(|r| (r[key_idx].clone(), r)).collect();
    let mut removed = Vec::new();
    let mut changed = Vec::new();
    for (key, left_row) in &left_map {
        match right_map.get(key) {
            None => removed.push((*left_row).clone()),
            Some(right_row) => {
                if left_row != right_row {
                    changed.push(CsvChangedRow {
                        key: key.clone(),
                        before: (*left_row).clone(),
                        after: (*right_row).clone(),
                    });
                }
            }
        }
    }
    let added: Vec<Vec<String>> = right_map
        .iter()
        .filter(|(key, _)| !left_map.contains_key(*key))
        .map(|(_, row)| (*row).clone())
        .collect();
    Ok(CsvDiff {
        removed,
        added,
        changed,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::CsvReader;
    use crate::CsvWriter;
    use crate::csv::ConfigurableCsvWriter;
    use crate::csv::CsvParser;
    use crate::csv::CsvRecord;
    use crate::csv::CsvStreamParser;
    use crate::csv::CsvWriterConfig;
    use crate::csv::InMemoryCsvReader;
    use crate::csv::InMemoryCsvWriter;
    use crate::csv::TypedCsvReader;
    use crate::csv::types::*;
    use std::fs::File;
    use std::io::Write;
    #[test]
    fn test_csv_write_and_read_roundtrip() {
        let path = std::env::temp_dir().join("oxiphy_test.csv");
        {
            let mut w = CsvWriter::new(
                path.to_str().unwrap_or(""),
                &["time", "energy", "temperature"],
            )
            .unwrap();
            w.write_row(&[0.0, 100.0, 300.0]).unwrap();
            w.write_row(&[1.0, 99.5, 299.8]).unwrap();
            w.write_row(&[2.0, 99.0, 299.5]).unwrap();
        }
        let (headers, rows) = CsvReader::read(path.to_str().unwrap_or("")).unwrap();
        assert_eq!(headers, vec!["time", "energy", "temperature"]);
        assert_eq!(rows.len(), 3);
        assert!((rows[0][0] - 0.0).abs() < 1e-10);
        assert!((rows[2][2] - 299.5).abs() < 1e-10);
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_csv_write_read() {
        let path = std::env::temp_dir().join("oxiphy_test_positions.csv");
        let particle_positions: Vec<[f64; 3]> = vec![
            [1.0, 2.0, 3.0],
            [4.5, 5.5, 6.5],
            [-1.0, -2.0, -3.0],
            [0.0, 0.0, 0.0],
        ];
        {
            let mut w = CsvWriter::new(path.to_str().unwrap_or(""), &["x", "y", "z"]).unwrap();
            for pos in &particle_positions {
                w.write_row(pos).unwrap();
            }
        }
        let (headers, rows) = CsvReader::read(path.to_str().unwrap_or("")).unwrap();
        assert_eq!(headers, vec!["x", "y", "z"]);
        assert_eq!(rows.len(), 4, "expected 4 rows");
        for (i, expected) in particle_positions.iter().enumerate() {
            assert!(
                (rows[i][0] - expected[0]).abs() < 1e-10,
                "row {} x mismatch: {} vs {}",
                i,
                rows[i][0],
                expected[0]
            );
            assert!(
                (rows[i][1] - expected[1]).abs() < 1e-10,
                "row {} y mismatch: {} vs {}",
                i,
                rows[i][1],
                expected[1]
            );
            assert!(
                (rows[i][2] - expected[2]).abs() < 1e-10,
                "row {} z mismatch: {} vs {}",
                i,
                rows[i][2],
                expected[2]
            );
        }
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_csv_with_empty_data() {
        let path = std::env::temp_dir().join("oxiphy_test_empty.csv");
        {
            let _w = CsvWriter::new(path.to_str().unwrap_or(""), &["x", "y"]).unwrap();
        }
        let (headers, rows) = CsvReader::read(path.to_str().unwrap_or("")).unwrap();
        assert_eq!(headers, vec!["x", "y"]);
        assert!(rows.is_empty());
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_in_memory_writer_header() {
        let w = InMemoryCsvWriter::new(&["x", "y", "z"], ',');
        assert_eq!(w.write_header(), "x,y,z");
    }
    #[test]
    fn test_in_memory_writer_tab_delimiter() {
        let w = InMemoryCsvWriter::new(&["a", "b"], '\t');
        assert_eq!(w.write_header(), "a\tb");
        let row = w.write_row(&[1.0, 2.0]).unwrap();
        assert!(row.contains('\t'));
    }
    #[test]
    fn test_in_memory_writer_row_wrong_len() {
        let w = InMemoryCsvWriter::new(&["x", "y"], ',');
        assert!(w.write_row(&[1.0]).is_err());
    }
    #[test]
    fn test_in_memory_writer_write_all() {
        let w = InMemoryCsvWriter::new(&["t", "v"], ',').with_precision(2);
        let rows = vec![vec![0.0, 1.0], vec![1.0, 2.5]];
        let out = w.write_all(&rows);
        assert!(out.starts_with("t,v\n"));
        assert!(out.contains("0.00"));
        assert!(out.contains("2.50"));
    }
    #[test]
    fn test_in_memory_writer_precision() {
        let w = InMemoryCsvWriter::new(&["v"], ',').with_precision(3);
        let row = w.write_row(&[2.54321]).unwrap();
        assert_eq!(row, "2.543");
    }
    #[test]
    fn test_in_memory_reader_basic() {
        let data = "x,y,z\n1.0,2.0,3.0\n4.0,5.0,6.0\n";
        let reader = data.parse::<InMemoryCsvReader>().unwrap();
        assert_eq!(reader.get_row_count(), 2);
        let x = reader.get_column_f64("x").unwrap();
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_in_memory_reader_missing_values() {
        let data = "a,b\n1.0,\n,3.0\n";
        let reader = data.parse::<InMemoryCsvReader>().unwrap();
        let a = reader.get_column_f64("a").unwrap();
        let b = reader.get_column_f64("b").unwrap();
        assert!((a[0] - 1.0).abs() < 1e-10);
        assert!(a[1].is_nan());
        assert!(b[0].is_nan());
        assert!((b[1] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_in_memory_reader_comment_lines() {
        let data = "# This is a comment\nx,y\n# Another comment\n1.0,2.0\n3.0,4.0\n";
        let reader = data.parse::<InMemoryCsvReader>().unwrap();
        assert_eq!(reader.get_row_count(), 2);
        assert_eq!(reader.headers(), &["x", "y"]);
    }
    #[test]
    fn test_in_memory_reader_semicolon_delimiter() {
        let data = "a;b;c\n1.0;2.0;3.0\n";
        let reader = InMemoryCsvReader::parse_with_delimiter(data, ';').unwrap();
        let b = reader.get_column_f64("b").unwrap();
        assert!((b[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_in_memory_reader_empty_input() {
        assert!("".parse::<InMemoryCsvReader>().is_err());
        assert!("# only comments\n".parse::<InMemoryCsvReader>().is_err());
    }
    #[test]
    fn test_in_memory_reader_column_not_found() {
        let data = "x,y\n1.0,2.0\n";
        let reader = data.parse::<InMemoryCsvReader>().unwrap();
        assert!(reader.get_column_f64("z").is_err());
    }
    #[test]
    fn test_column_stats() {
        let data = "v\n1.0\n2.0\n3.0\n4.0\n5.0\n";
        let reader = data.parse::<InMemoryCsvReader>().unwrap();
        let (min, max, mean, std) = reader.column_stats("v").unwrap();
        assert!((min - 1.0).abs() < 1e-10);
        assert!((max - 5.0).abs() < 1e-10);
        assert!((mean - 3.0).abs() < 1e-10);
        assert!((std - 2.0_f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_column_stats_with_missing() {
        let data = "v\n2.0\n\n4.0\n";
        let reader = data.parse::<InMemoryCsvReader>().unwrap();
        let (min, max, mean, _std) = reader.column_stats("v").unwrap();
        assert!((min - 2.0).abs() < 1e-10);
        assert!((max - 4.0).abs() < 1e-10);
        assert!((mean - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_write_read_roundtrip_in_memory() {
        let writer = InMemoryCsvWriter::new(&["time", "energy"], ',').with_precision(4);
        let rows = vec![vec![0.0, 100.0], vec![0.5, 99.5], vec![1.0, 99.0]];
        let csv_str = writer.write_all(&rows);
        let reader = csv_str.parse::<InMemoryCsvReader>().unwrap();
        assert_eq!(reader.get_row_count(), 3);
        let time = reader.get_column_f64("time").unwrap();
        let energy = reader.get_column_f64("energy").unwrap();
        assert!((time[1] - 0.5).abs() < 1e-3);
        assert!((energy[2] - 99.0).abs() < 1e-3);
    }
    #[test]
    fn test_quoted_fields() {
        let data = "name,value\n\"hello,world\",42.0\n";
        let reader = data.parse::<InMemoryCsvReader>().unwrap();
        assert_eq!(reader.get_row_count(), 1);
        let value = reader.get_column_f64("value").unwrap();
        assert!((value[0] - 42.0).abs() < 1e-10);
    }
    #[test]
    fn test_csv_parser_basic() {
        let data = "name,age\nAlice,30\nBob,25\n";
        let parser = CsvParser::new(data, ',');
        let records = parser.parse_all().unwrap();
        assert_eq!(records.len(), 3);
        assert_eq!(records[0].get(0), "name");
        assert_eq!(records[1].get(0), "Alice");
        assert_eq!(records[2].get(1), "25");
    }
    #[test]
    fn test_csv_parser_quoted_with_delimiter() {
        let data = "a,b\n\"x,y\",2\n";
        let parser = CsvParser::new(data, ',');
        let records = parser.parse_all().unwrap();
        assert_eq!(records[1].get(0), "x,y");
        assert_eq!(records[1].get(1), "2");
    }
    #[test]
    fn test_csv_parser_escaped_quote() {
        let data = "q\n\"he said \"\"hi\"\"\"\n";
        let parser = CsvParser::new(data, ',');
        let records = parser.parse_all().unwrap();
        assert_eq!(records[1].get(0), "he said \"hi\"");
    }
    #[test]
    fn test_csv_parser_backslash_escape() {
        let data = "s\n\"line1\\nline2\"\n";
        let parser = CsvParser::new(data, ',');
        let records = parser.parse_all().unwrap();
        assert_eq!(records[1].get(0), "line1\nline2");
    }
    #[test]
    fn test_csv_parser_comment_skip() {
        let data = "# header comment\nx,y\n# data comment\n1,2\n";
        let parser = CsvParser::new(data, ',').with_comment_prefix('#');
        let records = parser.parse_all().unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].get(0), "x");
    }
    #[test]
    fn test_csv_parser_empty_fields() {
        let data = "a,b,c\n1,,3\n";
        let parser = CsvParser::new(data, ',');
        let records = parser.parse_all().unwrap();
        assert_eq!(records[1].get(1), "");
        assert_eq!(records[1].get(2), "3");
    }
    #[test]
    fn test_csv_table_from_str() {
        let data = "x,y\n1,2\n3,4\n";
        let table = CsvTable::parse(data, ',').unwrap();
        assert_eq!(table.row_count(), 2);
        assert_eq!(table.col_count(), 2);
    }
    #[test]
    fn test_csv_table_column_f64() {
        let data = "x,y\n1.5,2.5\n3.5,4.5\n";
        let table = CsvTable::parse(data, ',').unwrap();
        let x = table.column_f64("x").unwrap();
        assert!((x[0] - 1.5).abs() < 1e-10);
        assert!((x[1] - 3.5).abs() < 1e-10);
    }
    #[test]
    fn test_csv_table_to_csv_string() {
        let data = "a,b\n1,2\n3,4\n";
        let table = CsvTable::parse(data, ',').unwrap();
        let out = table.to_csv_string(',');
        assert!(out.starts_with("a,b\n"));
        assert!(out.contains("1,2"));
    }
    #[test]
    fn test_csv_table_quote_field_with_delimiter() {
        let data = "name,val\n\"x,y\",5\n";
        let table = CsvTable::parse(data, ',').unwrap();
        let names = table.column_values("name").unwrap();
        assert_eq!(names[0], "x,y");
    }
    #[test]
    fn test_configurable_writer_semicolon() {
        let cfg = CsvWriterConfig {
            delimiter: ';',
            precision: 2,
            ..Default::default()
        };
        let mut w = ConfigurableCsvWriter::new(cfg);
        w.write_header(&["x", "y"]);
        w.write_f64_row(&[1.0, 2.5]);
        let out = w.finish();
        assert!(out.starts_with("x;y\n"));
        assert!(out.contains("1.00;2.50"));
    }
    #[test]
    fn test_configurable_writer_quote_all() {
        let cfg = CsvWriterConfig {
            quote_all: true,
            ..Default::default()
        };
        let mut w = ConfigurableCsvWriter::new(cfg);
        w.write_str_row(&["hello", "world"]);
        let out = w.finish();
        assert!(out.contains("\"hello\""));
        assert!(out.contains("\"world\""));
    }
    #[test]
    fn test_configurable_writer_crlf() {
        let cfg = CsvWriterConfig {
            line_ending: "\r\n".to_string(),
            ..Default::default()
        };
        let mut w = ConfigurableCsvWriter::new(cfg);
        w.write_header(&["v"]);
        let out = w.finish();
        assert!(out.ends_with("\r\n"));
    }
    #[test]
    fn test_infer_integer() {
        assert_eq!(infer_column_type(&["1", "2", "3"]), ColumnType::Integer);
        assert_eq!(infer_column_type(&["-5", "0", "100"]), ColumnType::Integer);
    }
    #[test]
    fn test_infer_float() {
        assert_eq!(
            infer_column_type(&["1.0", "2.5", "3.14"]),
            ColumnType::Float
        );
        assert_eq!(infer_column_type(&["1e10", "-0.5"]), ColumnType::Float);
    }
    #[test]
    fn test_infer_bool() {
        assert_eq!(infer_column_type(&["true", "false"]), ColumnType::Boolean);
        assert_eq!(infer_column_type(&["yes", "no"]), ColumnType::Boolean);
        assert_eq!(infer_column_type(&["1", "0"]), ColumnType::Integer);
    }
    #[test]
    fn test_infer_text() {
        assert_eq!(infer_column_type(&["Alice", "Bob"]), ColumnType::Text);
        assert_eq!(infer_column_type(&["1.0", "abc"]), ColumnType::Text);
    }
    #[test]
    fn test_infer_empty() {
        assert_eq!(infer_column_type(&["", ""]), ColumnType::Empty);
        assert_eq!(infer_column_type(&[]), ColumnType::Empty);
    }
    #[test]
    fn test_infer_with_missing() {
        assert_eq!(infer_column_type(&["1.5", "", "3.5"]), ColumnType::Float);
    }
    #[test]
    fn test_typed_reader_i64() {
        let data = "id,val\n1,10\n2,20\n";
        let reader = data.parse::<TypedCsvReader>().unwrap();
        let ids = reader.column_as_i64("id").unwrap();
        assert_eq!(ids, vec![1, 2]);
    }
    #[test]
    fn test_typed_reader_f64() {
        let data = "x\n1.5\n2.5\n";
        let reader = data.parse::<TypedCsvReader>().unwrap();
        let x = reader.column_as_f64("x").unwrap();
        assert!((x[0] - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_typed_reader_bool() {
        let data = "active\ntrue\nfalse\nyes\nno\n";
        let reader = data.parse::<TypedCsvReader>().unwrap();
        let active = reader.column_as_bool("active").unwrap();
        assert_eq!(active, vec![true, false, true, false]);
    }
    #[test]
    fn test_typed_reader_column_type() {
        let data = "id,name,score\n1,Alice,9.5\n2,Bob,8.0\n";
        let reader = data.parse::<TypedCsvReader>().unwrap();
        assert_eq!(reader.column_type("id").unwrap(), ColumnType::Integer);
        assert_eq!(reader.column_type("name").unwrap(), ColumnType::Text);
        assert_eq!(reader.column_type("score").unwrap(), ColumnType::Float);
    }
    #[test]
    fn test_typed_reader_headers_and_count() {
        let data = "a,b\n1,2\n3,4\n5,6\n";
        let reader = data.parse::<TypedCsvReader>().unwrap();
        assert_eq!(reader.row_count(), 3);
        assert_eq!(reader.headers(), &["a", "b"]);
    }
    #[test]
    fn test_csv_merge_basic() {
        let a = CsvTable::parse("x,y\n1,2\n", ',').unwrap();
        let b = CsvTable::parse("x,y\n3,4\n", ',').unwrap();
        let merged = csv_merge(&a, &b).unwrap();
        assert_eq!(merged.row_count(), 2);
        let x = merged.column_f64("x").unwrap();
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_csv_merge_header_mismatch() {
        let a = CsvTable::parse("x,y\n1,2\n", ',').unwrap();
        let b = CsvTable::parse("x,z\n3,4\n", ',').unwrap();
        assert!(csv_merge(&a, &b).is_err());
    }
    #[test]
    fn test_csv_merge_empty_right() {
        let a = CsvTable::parse("x\n1\n2\n", ',').unwrap();
        let b = CsvTable::new(vec!["x".to_string()]);
        let merged = csv_merge(&a, &b).unwrap();
        assert_eq!(merged.row_count(), 2);
    }
    #[test]
    fn test_csv_join_basic() {
        let left = CsvTable::parse("id,x\n1,10\n2,20\n3,30\n", ',').unwrap();
        let right = CsvTable::parse("id,y\n1,100\n3,300\n", ',').unwrap();
        let joined = csv_join(&left, &right, "id").unwrap();
        assert_eq!(joined.row_count(), 2);
        let y = joined.column_f64("y").unwrap();
        assert!((y[0] - 100.0).abs() < 1e-10);
        assert!((y[1] - 300.0).abs() < 1e-10);
    }
    #[test]
    fn test_csv_join_no_matches() {
        let left = CsvTable::parse("id,x\n1,10\n", ',').unwrap();
        let right = CsvTable::parse("id,y\n9,99\n", ',').unwrap();
        let joined = csv_join(&left, &right, "id").unwrap();
        assert_eq!(joined.row_count(), 0);
    }
    #[test]
    fn test_csv_join_key_not_found() {
        let left = CsvTable::parse("id,x\n1,10\n", ',').unwrap();
        let right = CsvTable::parse("id,y\n1,100\n", ',').unwrap();
        assert!(csv_join(&left, &right, "missing").is_err());
    }
    #[test]
    fn test_csv_pivot_sum() {
        let data = "region,product,sales\nNorth,A,10\nNorth,B,20\nSouth,A,30\nSouth,B,40\n";
        let table = CsvTable::parse(data, ',').unwrap();
        let pivot = csv_pivot(&table, "region", "product", "sales", PivotAgg::Sum).unwrap();
        assert_eq!(pivot.row_count(), 2);
        assert!(pivot.headers.contains(&"A".to_string()));
        assert!(pivot.headers.contains(&"B".to_string()));
    }
    #[test]
    fn test_csv_pivot_count() {
        let data = "cat,sub,v\nA,X,1\nA,X,2\nA,Y,3\nB,X,4\n";
        let table = CsvTable::parse(data, ',').unwrap();
        let pivot = csv_pivot(&table, "cat", "sub", "v", PivotAgg::Count).unwrap();
        let x_col = pivot.column_values("X").unwrap();
        assert_eq!(x_col[0], "2");
    }
    #[test]
    fn test_csv_pivot_mean() {
        let data = "g,c,v\nA,X,10\nA,X,20\n";
        let table = CsvTable::parse(data, ',').unwrap();
        let pivot = csv_pivot(&table, "g", "c", "v", PivotAgg::Mean).unwrap();
        let x_col = pivot.column_values("X").unwrap();
        let mean: f64 = x_col[0].parse().unwrap();
        assert!((mean - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_csv_diff_basic() {
        let a = CsvTable::parse("id,v\n1,10\n2,20\n3,30\n", ',').unwrap();
        let b = CsvTable::parse("id,v\n1,10\n2,99\n4,40\n", ',').unwrap();
        let diff = csv_diff(&a, &b, "id").unwrap();
        assert_eq!(diff.removed.len(), 1);
        assert_eq!(diff.added.len(), 1);
        assert_eq!(diff.changed.len(), 1);
        assert_eq!(diff.changed[0].key, "2");
    }
    #[test]
    fn test_csv_diff_no_changes() {
        let a = CsvTable::parse("id,v\n1,10\n2,20\n", ',').unwrap();
        let b = CsvTable::parse("id,v\n1,10\n2,20\n", ',').unwrap();
        let diff = csv_diff(&a, &b, "id").unwrap();
        assert_eq!(diff.removed.len(), 0);
        assert_eq!(diff.added.len(), 0);
        assert_eq!(diff.changed.len(), 0);
    }
    #[test]
    fn test_csv_diff_header_mismatch() {
        let a = CsvTable::parse("id,v\n1,10\n", ',').unwrap();
        let b = CsvTable::parse("id,w\n1,10\n", ',').unwrap();
        assert!(csv_diff(&a, &b, "id").is_err());
    }
    #[test]
    fn test_csv_diff_all_added() {
        let a = CsvTable::new(vec!["id".to_string(), "v".to_string()]);
        let b = CsvTable::parse("id,v\n1,10\n", ',').unwrap();
        let diff = csv_diff(&a, &b, "id").unwrap();
        assert_eq!(diff.added.len(), 1);
        assert_eq!(diff.removed.len(), 0);
    }
    #[test]
    fn test_csv_stream_parser_file() {
        let path = std::env::temp_dir().join("oxiphy_stream_test.csv");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, "x,y,z").unwrap();
            writeln!(f, "1,2,3").unwrap();
            writeln!(f, "4,5,6").unwrap();
            writeln!(f, "7,8,9").unwrap();
        }
        let mut parser = CsvStreamParser::open(path.to_str().unwrap_or(""), ',').unwrap();
        assert_eq!(parser.headers(), &["x", "y", "z"]);
        let mut count = 0;
        while let Some(rec) = parser.next_record().unwrap() {
            count += 1;
            assert_eq!(rec.len(), 3);
        }
        assert_eq!(count, 3);
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_csv_stream_parser_large() {
        let path = std::env::temp_dir().join("oxiphy_stream_large.csv");
        let nrows = 500_usize;
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, "i,v").unwrap();
            for i in 0..nrows {
                writeln!(f, "{},{}", i, i as f64 * 1.5).unwrap();
            }
        }
        let mut parser = CsvStreamParser::open(path.to_str().unwrap_or(""), ',').unwrap();
        let mut count = 0usize;
        while let Some(_rec) = parser.next_record().unwrap() {
            count += 1;
        }
        assert_eq!(count, nrows);
        std::fs::remove_file(&path).ok();
    }
    #[test]
    fn test_csv_record_accessors() {
        let rec = CsvRecord {
            fields: vec!["a".to_string(), "b".to_string(), "c".to_string()],
        };
        assert_eq!(rec.len(), 3);
        assert!(!rec.is_empty());
        assert_eq!(rec.get(0), "a");
        assert_eq!(rec.get(10), "");
    }
    #[test]
    fn test_quote_field_no_quote_needed() {
        assert_eq!(quote_field("hello", ','), "hello");
    }
    #[test]
    fn test_quote_field_with_delimiter() {
        assert_eq!(quote_field("a,b", ','), "\"a,b\"");
    }
    #[test]
    fn test_quote_field_with_embedded_quote() {
        assert_eq!(quote_field("say \"hi\"", ','), "\"say \"\"hi\"\"\"");
    }
    #[test]
    fn test_infer_table_types_mixed() {
        let data = "id,name,score\n1,Alice,9.5\n2,Bob,8.0\n";
        let table = CsvTable::parse(data, ',').unwrap();
        let types = infer_table_types(&table);
        assert_eq!(types[0].1, ColumnType::Integer);
        assert_eq!(types[1].1, ColumnType::Text);
        assert_eq!(types[2].1, ColumnType::Float);
    }
}
