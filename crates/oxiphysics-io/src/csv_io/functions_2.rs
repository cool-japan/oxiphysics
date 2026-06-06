//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests_csv_extra {
    use crate::csv_io::*;
    #[test]
    fn csv_file_new_empty() {
        let f = CsvFile::new(vec!["a".into(), "b".into()]);
        assert_eq!(f.record_count(), 0);
        assert_eq!(f.column_count(), 2);
    }
    #[test]
    fn csv_file_add_record_f64_roundtrip() {
        let mut f = CsvFile::new(vec!["x".into(), "y".into()]);
        f.add_record_f64(&[1.5, 2.5]);
        let vals = f.get_column_f64(0).unwrap();
        assert!((vals[0] - 1.5).abs() < 1e-12);
    }
    #[test]
    fn csv_file_get_column_out_of_range() {
        let f = CsvFile::new(vec!["x".into()]);
        assert!(f.get_column_f64(5).is_err());
    }
    #[test]
    fn csv_file_get_column_by_name_found() {
        let f = CsvFile::new(vec!["time".into(), "energy".into()]);
        assert_eq!(f.get_column_by_name("energy"), Some(1));
    }
    #[test]
    fn csv_file_get_column_by_name_missing() {
        let f = CsvFile::new(vec!["x".into()]);
        assert_eq!(f.get_column_by_name("z"), None);
    }
    #[test]
    fn csv_file_to_string_with_delimiter_tab() {
        let mut f = CsvFile::new(vec!["a".into(), "b".into()]);
        f.add_record(vec!["1".into(), "2".into()]);
        let s = f.to_string_with_delimiter('\t');
        assert!(s.contains("a\tb"));
        assert!(s.contains("1\t2"));
    }
    #[test]
    fn csv_file_from_str_with_delimiter_semicolon() {
        let data = "x;y\n1;2\n3;4\n";
        let f = CsvFile::from_str_with_delimiter(data, ';').unwrap();
        assert_eq!(f.column_count(), 2);
        assert_eq!(f.record_count(), 2);
    }
    #[test]
    fn csv_file_roundtrip_comma() {
        let mut f = CsvFile::new(vec!["t".into(), "v".into()]);
        f.add_record_f64(&[0.0, 1.0]);
        f.add_record_f64(&[1.0, 2.0]);
        let s = f.to_string();
        let f2 = &s.parse::<CsvFile>().unwrap();
        assert_eq!(f2.record_count(), 2);
        let v = f2.get_column_f64(1).unwrap();
        assert!((v[1] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn csv_file_filter_rows_positive() {
        let mut f = CsvFile::new(vec!["x".into()]);
        f.add_record(vec!["-1".into()]);
        f.add_record(vec!["2".into()]);
        f.add_record(vec!["3".into()]);
        let filtered = f.filter_rows(0, |v| v > 0.0);
        assert_eq!(filtered.record_count(), 2);
    }
    #[test]
    fn csv_file_select_columns_by_name() {
        let mut f = CsvFile::new(vec!["t".into(), "x".into(), "y".into()]);
        f.add_record(vec!["0".into(), "1".into(), "2".into()]);
        let sub = f.select_columns_by_name(&["x", "y"]);
        assert_eq!(sub.column_count(), 2);
        assert_eq!(sub.headers[0], "x");
    }
    #[test]
    fn csv_file_normalize_headers_lowercases() {
        let mut f = CsvFile::new(vec!["Time (s)".into(), "Force (N)".into()]);
        f.normalize_headers();
        assert_eq!(f.headers[0], "time__s_");
        assert_eq!(f.headers[1], "force__n_");
    }
    #[test]
    fn csv_file_column_stats_basic() {
        let mut f = CsvFile::new(vec!["v".into()]);
        for v in [1.0, 2.0, 3.0, 4.0] {
            f.add_record_f64(&[v]);
        }
        let stats = f.column_stats(0).unwrap();
        assert!((stats.min - 1.0).abs() < 1e-12);
        assert!((stats.max - 4.0).abs() < 1e-12);
        assert!((stats.mean - 2.5).abs() < 1e-12);
        assert_eq!(stats.count, 4);
    }
    #[test]
    fn csv_file_column_stats_empty_column() {
        let f = CsvFile::new(vec!["v".into()]);
        assert!(f.column_stats(0).is_none());
    }
    #[test]
    fn csv_file_get_column_strings() {
        let mut f = CsvFile::new(vec!["name".into()]);
        f.add_record(vec!["Alice".into()]);
        f.add_record(vec!["Bob".into()]);
        let names = f.get_column_strings(0).unwrap();
        assert_eq!(names, vec!["Alice", "Bob"]);
    }
    #[test]
    fn csv_file_get_column_i64() {
        let mut f = CsvFile::new(vec!["n".into()]);
        f.add_record(vec!["42".into()]);
        f.add_record(vec!["100".into()]);
        let ns = f.get_column_i64(0).unwrap();
        assert_eq!(ns, vec![42_i64, 100]);
    }
    #[test]
    fn csv_file_sort_by_column() {
        let mut f = CsvFile::new(vec!["v".into()]);
        f.add_record(vec!["3".into()]);
        f.add_record(vec!["1".into()]);
        f.add_record(vec!["2".into()]);
        f.sort_by_column(0);
        let vals = f.get_column_f64(0).unwrap();
        assert_eq!(vals, vec![1.0, 2.0, 3.0]);
    }
    #[test]
    fn detect_delimiter_comma() {
        assert_eq!(detect_delimiter("a,b,c\n1,2,3"), ',');
    }
    #[test]
    fn detect_delimiter_tab() {
        assert_eq!(detect_delimiter("a\tb\tc\n1\t2\t3"), '\t');
    }
    #[test]
    fn detect_delimiter_semicolon() {
        assert_eq!(detect_delimiter("a;b;c\n1;2;3"), ';');
    }
    #[test]
    fn parse_auto_with_comma() {
        let data = "x,y\n1,2\n3,4\n";
        let f = parse_auto(data).unwrap();
        assert_eq!(f.record_count(), 2);
    }
    #[test]
    fn read_chunked_basic() {
        let data = "x,y\n1,2\n3,4\n5,6\n7,8\n";
        let chunks = read_chunked(data, 2);
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].record_count(), 2);
    }
    #[test]
    fn read_chunked_zero_size_returns_full() {
        let data = "x\n1\n2\n";
        let chunks = read_chunked(data, 0);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0].record_count(), 2);
    }
    #[test]
    fn csv_schema_len_matches_columns() {
        let schema = CsvSchema::new(vec![
            ("x".into(), ColumnType::Float),
            ("n".into(), ColumnType::Integer),
        ]);
        assert_eq!(schema.len(), 2);
        assert!(!schema.is_empty());
    }
    #[test]
    fn csv_schema_is_empty_when_no_cols() {
        let schema = CsvSchema::new(vec![]);
        assert!(schema.is_empty());
    }
    #[test]
    fn aggregate_column_sum() {
        let mut f = CsvFile::new(vec!["v".into()]);
        for v in [1.0, 2.0, 3.0] {
            f.add_record_f64(&[v]);
        }
        let sum = aggregate_column(&f, 0, AggOp::Sum).unwrap();
        assert!((sum - 6.0).abs() < 1e-12);
    }
    #[test]
    fn aggregate_column_mean() {
        let mut f = CsvFile::new(vec!["v".into()]);
        for v in [10.0, 20.0, 30.0] {
            f.add_record_f64(&[v]);
        }
        let mean = aggregate_column(&f, 0, AggOp::Mean).unwrap();
        assert!((mean - 20.0).abs() < 1e-12);
    }
    #[test]
    fn aggregate_column_min_max() {
        let mut f = CsvFile::new(vec!["v".into()]);
        for v in [5.0, 2.0, 8.0, 1.0] {
            f.add_record_f64(&[v]);
        }
        let min = aggregate_column(&f, 0, AggOp::Min).unwrap();
        let max = aggregate_column(&f, 0, AggOp::Max).unwrap();
        assert!((min - 1.0).abs() < 1e-12);
        assert!((max - 8.0).abs() < 1e-12);
    }
    #[test]
    fn aggregate_column_count() {
        let mut f = CsvFile::new(vec!["v".into()]);
        for v in [1.0, 2.0, 3.0, 4.0, 5.0] {
            f.add_record_f64(&[v]);
        }
        let count = aggregate_column(&f, 0, AggOp::Count).unwrap();
        assert!((count - 5.0).abs() < 1e-12);
    }
    #[test]
    fn aggregate_column_std_nonneg() {
        let mut f = CsvFile::new(vec!["v".into()]);
        for v in [1.0, 2.0, 3.0] {
            f.add_record_f64(&[v]);
        }
        let std_val = aggregate_column(&f, 0, AggOp::Std).unwrap();
        assert!(std_val >= 0.0);
    }
    #[test]
    fn timeseries_csv_n_steps() {
        let data = "t,v\n0.0,1.0\n1.0,2.0\n2.0,3.0\n";
        let ts = TimeSeriesCsv::parse(data, "t").unwrap();
        assert_eq!(ts.n_steps(), 3);
    }
    #[test]
    fn timeseries_csv_times() {
        let data = "t,v\n0.0,10.0\n0.5,20.0\n1.0,30.0\n";
        let ts = TimeSeriesCsv::parse(data, "t").unwrap();
        let times = ts.times().unwrap();
        assert_eq!(times.len(), 3);
        assert!((times[1] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn timeseries_csv_duration() {
        let data = "t,x\n0.0,1.0\n2.0,3.0\n";
        let ts = TimeSeriesCsv::parse(data, "t").unwrap();
        assert!((ts.duration() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn timeseries_csv_duration_single_row() {
        let data = "t,x\n1.0,2.0\n";
        let ts = TimeSeriesCsv::parse(data, "t").unwrap();
        assert_eq!(ts.duration(), 0.0);
    }
    #[test]
    fn timeseries_csv_column_f64() {
        let data = "t,v\n0.0,5.0\n1.0,10.0\n";
        let ts = TimeSeriesCsv::parse(data, "t").unwrap();
        let v = ts.column_f64("v").unwrap();
        assert!((v[0] - 5.0).abs() < 1e-12);
        assert!((v[1] - 10.0).abs() < 1e-12);
    }
    #[test]
    fn csv_writer_basic_output() {
        let mut w = CsvWriter::new(vec!["x".into(), "y".into()], ',');
        w.write_row(&["1", "2"]);
        w.write_row(&["3", "4"]);
        let s = w.finish();
        assert!(s.contains("x,y"));
        assert!(s.contains("1,2"));
        assert!(s.contains("3,4"));
    }
    #[test]
    fn csv_writer_row_count() {
        let mut w = CsvWriter::new(vec!["v".into()], ',');
        w.write_row_f64(&[1.0]);
        w.write_row_f64(&[2.0]);
        w.write_row_f64(&[3.0]);
        assert_eq!(w.row_count(), 3);
    }
    #[test]
    fn csv_writer_tab_delimiter() {
        let mut w = CsvWriter::new(vec!["a".into(), "b".into()], '\t');
        w.write_row(&["10", "20"]);
        let s = w.finish();
        assert!(s.contains("a\tb"));
        assert!(s.contains("10\t20"));
    }
    #[test]
    fn dataframe_from_csv_types() {
        let data = "n,x,name\n1,1.5,Alice\n2,2.5,Bob\n";
        let csv = data.parse::<CsvFile>().unwrap();
        let df = CsvDataFrame::from_csv(&csv);
        assert_eq!(df.n_cols(), 3);
        assert_eq!(df.n_rows(), 2);
        assert_eq!(
            df.column_by_name("n").unwrap().column_type(),
            ColumnType::Integer
        );
        assert_eq!(
            df.column_by_name("x").unwrap().column_type(),
            ColumnType::Float
        );
        assert_eq!(
            df.column_by_name("name").unwrap().column_type(),
            ColumnType::Text
        );
    }
    #[test]
    fn dataframe_float_column() {
        let data = "v\n1.1\n2.2\n3.3\n";
        let csv = data.parse::<CsvFile>().unwrap();
        let df = CsvDataFrame::from_csv(&csv);
        let col = df.float_column("v").unwrap();
        assert_eq!(col.len(), 3);
        assert!((col[0] - 1.1).abs() < 1e-9);
    }
    #[test]
    fn dataframe_integer_column() {
        let data = "n\n10\n20\n30\n";
        let csv = data.parse::<CsvFile>().unwrap();
        let df = CsvDataFrame::from_csv(&csv);
        let col = df.integer_column("n").unwrap();
        assert_eq!(col, &vec![10_i64, 20, 30]);
    }
    #[test]
    fn dataframe_text_column() {
        let data = "label\nalpha\nbeta\n";
        let csv = data.parse::<CsvFile>().unwrap();
        let df = CsvDataFrame::from_csv(&csv);
        let col = df.text_column("label").unwrap();
        assert_eq!(col[0], "alpha");
    }
    #[test]
    fn dataframe_column_index_by_name() {
        let data = "a,b,c\n1,2,3\n";
        let csv = data.parse::<CsvFile>().unwrap();
        let df = CsvDataFrame::from_csv(&csv);
        assert_eq!(df.column_index("b"), Some(1));
        assert_eq!(df.column_index("z"), None);
    }
    #[test]
    fn all_column_stats_mixed() {
        let data = "n,label\n1,foo\n2,bar\n3,baz\n";
        let csv = data.parse::<CsvFile>().unwrap();
        let stats = csv.all_column_stats();
        assert_eq!(stats.len(), 1);
        assert_eq!(stats[0].0, "n");
    }
    #[test]
    fn csv_column_data_len_empty() {
        let col = CsvColumnData::Float(vec![]);
        assert_eq!(col.len(), 0);
        assert!(col.is_empty());
    }
    #[test]
    fn csv_column_data_type_integer() {
        let col = CsvColumnData::Integer(vec![1, 2, 3]);
        assert_eq!(col.column_type(), ColumnType::Integer);
        assert_eq!(col.len(), 3);
    }
}
