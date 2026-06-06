//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod tests {

    use crate::database_io::CacheEntry;
    use crate::database_io::DataCatalog;
    use crate::database_io::DatabaseQuery;
    use crate::database_io::DatabaseSerializer;
    use crate::database_io::DbRow;
    use crate::database_io::DbValue;
    use crate::database_io::ExportFilter;
    use crate::database_io::ExportFormat;
    use crate::database_io::ExportPipeline;
    use crate::database_io::MaterialDatabase;
    use crate::database_io::ResultCache;
    use crate::database_io::SimulationDatabase;
    use crate::database_io::SimulationMetadata;
    use crate::database_io::SimulationRecord;
    use crate::database_io::SimulationRecordDatabase;
    use crate::database_io::SnapshotCatalog;
    use crate::database_io::SnapshotEntry;
    use crate::database_io::TimeSeriesStore;

    #[test]
    fn test_dbvalue_as_f64_float() {
        let v = DbValue::Float(3.125);
        assert!((v.as_f64().unwrap() - 3.125).abs() < 1e-10);
    }
    #[test]
    fn test_dbvalue_as_f64_int() {
        let v = DbValue::Int(7);
        assert!((v.as_f64().unwrap() - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_dbvalue_as_f64_text_is_none() {
        let v = DbValue::Text("hello".into());
        assert!(v.as_f64().is_none());
    }
    #[test]
    fn test_dbvalue_as_str() {
        let v = DbValue::Text("world".into());
        assert_eq!(v.as_str(), Some("world"));
    }
    #[test]
    fn test_dbvalue_from_conversions() {
        let _: DbValue = 3.125f64.into();
        let _: DbValue = 42i64.into();
        let _: DbValue = "hello".into();
        let _: DbValue = true.into();
    }
    #[test]
    fn test_dbrow_set_get() {
        let mut row = DbRow::new();
        row.set("x", 1.5f64);
        assert_eq!(row.get("x"), Some(&DbValue::Float(1.5)));
    }
    #[test]
    fn test_dbrow_missing_key() {
        let row = DbRow::new();
        assert!(row.get("missing").is_none());
    }
    #[test]
    fn test_simdb_insert_count() {
        let mut db = SimulationDatabase::new("test");
        let mut row = DbRow::new();
        row.set("step", 1i64);
        db.insert(row);
        assert_eq!(db.row_count(), 1);
    }
    #[test]
    fn test_simdb_query_eq() {
        let mut db = SimulationDatabase::new("test");
        for i in 0i64..5 {
            let mut row = DbRow::new();
            row.set("id", i);
            db.insert(row);
        }
        let results = db.query_eq("id", &DbValue::Int(3));
        assert_eq!(results.len(), 1);
    }
    #[test]
    fn test_simdb_query_range() {
        let mut db = SimulationDatabase::new("test");
        for i in 0..10 {
            let mut row = DbRow::new();
            row.set("t", i as f64);
            db.insert(row);
        }
        let results = db.query_range("t", 3.0, 6.0);
        assert_eq!(results.len(), 4);
    }
    #[test]
    fn test_simdb_delete_eq() {
        let mut db = SimulationDatabase::new("test");
        for i in 0i64..5 {
            let mut row = DbRow::new();
            row.set("id", i);
            db.insert(row);
        }
        let deleted = db.delete_eq("id", &DbValue::Int(2));
        assert_eq!(deleted, 1);
        assert_eq!(db.row_count(), 4);
    }
    #[test]
    fn test_simdb_clear() {
        let mut db = SimulationDatabase::new("test");
        for _ in 0..3 {
            db.insert(DbRow::new());
        }
        db.clear();
        assert_eq!(db.row_count(), 0);
    }
    #[test]
    fn test_simdb_column_mean() {
        let mut db = SimulationDatabase::new("test");
        for i in 0..4 {
            let mut row = DbRow::new();
            row.set("v", i as f64);
            db.insert(row);
        }
        let mean = db.column_mean("v").unwrap();
        assert!((mean - 1.5).abs() < 1e-10, "mean = {mean}");
    }
    #[test]
    fn test_simdb_column_mean_empty() {
        let db = SimulationDatabase::new("test");
        assert!(db.column_mean("v").is_none());
    }
    #[test]
    fn test_simdb_project() {
        let mut db = SimulationDatabase::new("test");
        let mut row = DbRow::new();
        row.set("x", 1.0f64);
        row.set("y", 2.0f64);
        row.set("z", 3.0f64);
        db.insert(row);
        let proj = db.project(&["x", "z"]);
        assert_eq!(proj.len(), 1);
        assert!(proj[0].contains_key("x"));
        assert!(!proj[0].contains_key("y"));
    }
    #[test]
    fn test_ts_append_len() {
        let mut ts = TimeSeriesStore::new("pressure");
        ts.append(0.0, 100.0);
        ts.append(0.1, 101.0);
        assert_eq!(ts.len(), 2);
    }
    #[test]
    fn test_ts_range_query() {
        let mut ts = TimeSeriesStore::new("temp");
        for i in 0..10 {
            ts.append(i as f64 * 0.1, i as f64);
        }
        let results = ts.range_query(0.2, 0.5);
        assert_eq!(results.len(), 4);
    }
    #[test]
    fn test_ts_decimate() {
        let mut ts = TimeSeriesStore::new("v");
        for i in 0..10 {
            ts.append(i as f64, i as f64);
        }
        let dec = ts.decimate(2);
        assert_eq!(dec.len(), 5);
    }
    #[test]
    fn test_ts_decimate_zero() {
        let mut ts = TimeSeriesStore::new("v");
        ts.append(0.0, 1.0);
        let dec = ts.decimate(0);
        assert!(dec.is_empty());
    }
    #[test]
    fn test_ts_mean_in_range() {
        let mut ts = TimeSeriesStore::new("x");
        ts.append(0.0, 0.0);
        ts.append(1.0, 2.0);
        ts.append(2.0, 4.0);
        let mean = ts.mean_in_range(0.0, 2.0).unwrap();
        assert!((mean - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_ts_max_min() {
        let mut ts = TimeSeriesStore::new("y");
        ts.append(0.0, 3.0);
        ts.append(1.0, 7.0);
        ts.append(2.0, 1.0);
        assert!((ts.max_value().unwrap() - 7.0).abs() < 1e-10);
        assert!((ts.min_value().unwrap() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_ts_interpolate_midpoint() {
        let mut ts = TimeSeriesStore::new("f");
        ts.append(0.0, 0.0);
        ts.append(1.0, 10.0);
        let v = ts.interpolate(0.5).unwrap();
        assert!((v - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_ts_interpolate_empty() {
        let ts = TimeSeriesStore::new("f");
        assert!(ts.interpolate(0.5).is_none());
    }
    #[test]
    fn test_ts_to_csv() {
        let mut ts = TimeSeriesStore::new("p");
        ts.append(0.0, 1.0);
        let csv = ts.to_csv();
        assert!(csv.starts_with("time,value\n"));
        assert!(csv.contains("0"));
    }
    #[test]
    fn test_matdb_defaults_non_empty() {
        let db = MaterialDatabase::with_defaults();
        assert!(!db.is_empty());
    }
    #[test]
    fn test_matdb_lookup_steel() {
        let db = MaterialDatabase::with_defaults();
        let rec = db.lookup("steel_1020").unwrap();
        assert!((rec.density - 7850.0).abs() < 1.0);
    }
    #[test]
    fn test_matdb_lookup_missing() {
        let db = MaterialDatabase::with_defaults();
        assert!(db.lookup("unobtanium").is_none());
    }
    #[test]
    fn test_matdb_fuzzy_search() {
        let db = MaterialDatabase::with_defaults();
        let results = db.fuzzy_search("metal");
        assert!(!results.is_empty());
    }
    #[test]
    fn test_matdb_fuzzy_search_no_match() {
        let db = MaterialDatabase::with_defaults();
        let results = db.fuzzy_search("unobtanium_xyz");
        assert!(results.is_empty());
    }
    #[test]
    fn test_matdb_interpolate() {
        let db = MaterialDatabase::with_defaults();
        let blended = db.interpolate("steel_1020", "aluminium_6061", 0.5).unwrap();
        assert!(blended.density > 2700.0 && blended.density < 7850.0);
    }
    #[test]
    fn test_matdb_interpolate_t0() {
        let db = MaterialDatabase::with_defaults();
        let rec = db.interpolate("steel_1020", "copper", 0.0).unwrap();
        assert!((rec.density - 7850.0).abs() < 1.0);
    }
    #[test]
    fn test_matdb_names() {
        let db = MaterialDatabase::with_defaults();
        let names = db.names();
        assert!(names.contains(&"steel_1020"));
    }
    #[test]
    fn test_cache_put_get() {
        let mut cache = ResultCache::new(4);
        cache.put(CacheEntry::new("k1", vec![1.0, 2.0], 0.0, 0));
        let entry = cache.get("k1").unwrap();
        assert_eq!(entry.data, vec![1.0, 2.0]);
    }
    #[test]
    fn test_cache_eviction() {
        let mut cache = ResultCache::new(2);
        cache.put(CacheEntry::new("k1", vec![], 0.0, 0));
        cache.put(CacheEntry::new("k2", vec![], 0.0, 0));
        cache.put(CacheEntry::new("k3", vec![], 0.0, 0));
        assert!(cache.get("k1").is_none());
        assert!(cache.get("k2").is_some() || cache.get("k3").is_some());
    }
    #[test]
    fn test_cache_invalidate_all() {
        let mut cache = ResultCache::new(4);
        cache.put(CacheEntry::new("k1", vec![], 0.0, 0));
        cache.invalidate_all();
        assert!(cache.get("k1").is_none());
    }
    #[test]
    fn test_cache_remove() {
        let mut cache = ResultCache::new(4);
        cache.put(CacheEntry::new("k1", vec![], 0.0, 0));
        cache.remove("k1");
        assert!(cache.get("k1").is_none());
    }
    #[test]
    fn test_cache_evict_stale() {
        let mut cache = ResultCache::new(10);
        cache.put(CacheEntry::new("old", vec![], 0.0, 0));
        cache.invalidate_all();
        cache.put(CacheEntry::new("fresh", vec![], 1.0, 1));
        cache.evict_stale();
        assert_eq!(cache.len(), 1);
    }
    #[test]
    fn test_cache_lru_order() {
        let mut cache = ResultCache::new(3);
        cache.put(CacheEntry::new("a", vec![], 0.0, 0));
        cache.put(CacheEntry::new("b", vec![], 0.0, 0));
        cache.put(CacheEntry::new("c", vec![], 0.0, 0));
        let _ = cache.get("a");
        cache.put(CacheEntry::new("d", vec![], 0.0, 0));
        assert!(cache.get("a").is_some());
        assert!(cache.get("b").is_none());
    }
    #[test]
    fn test_catalog_register_lookup() {
        let mut cat = DataCatalog::new();
        let meta = SimulationMetadata::new("sim001", "Test run", 0.0);
        cat.register(meta);
        assert!(cat.lookup("sim001").is_some());
    }
    #[test]
    fn test_catalog_remove() {
        let mut cat = DataCatalog::new();
        cat.register(SimulationMetadata::new("s1", "desc", 0.0));
        let removed = cat.remove("s1");
        assert!(removed);
        assert!(cat.lookup("s1").is_none());
    }
    #[test]
    fn test_catalog_search_description() {
        let mut cat = DataCatalog::new();
        cat.register(SimulationMetadata::new("s1", "turbulence simulation", 0.0));
        cat.register(SimulationMetadata::new("s2", "heat transfer", 0.0));
        let results = cat.search_description("turbulence");
        assert_eq!(results.len(), 1);
    }
    #[test]
    fn test_catalog_metadata_params() {
        let mut cat = DataCatalog::new();
        cat.register(SimulationMetadata::new("s1", "run", 0.0));
        let meta = cat.lookup_mut("s1").unwrap();
        meta.set_param("dt", "0.01");
        meta.add_artifact(
            std::env::temp_dir()
                .join("result.vtk")
                .to_str()
                .unwrap_or(""),
        );
        let meta2 = cat.lookup("s1").unwrap();
        assert_eq!(meta2.parameters.get("dt"), Some(&"0.01".to_string()));
        assert_eq!(meta2.artifacts.len(), 1);
    }
    #[test]
    fn test_catalog_all_ids() {
        let mut cat = DataCatalog::new();
        cat.register(SimulationMetadata::new("a", "", 0.0));
        cat.register(SimulationMetadata::new("b", "", 0.0));
        let ids = cat.all_ids();
        assert_eq!(ids.len(), 2);
    }
    fn make_test_db() -> SimulationDatabase {
        let mut db = SimulationDatabase::new("runs");
        for i in 0..3 {
            let mut row = DbRow::new();
            row.set("step", i as i64);
            row.set("energy", (i as f64) * 10.0);
            row.set("label", format!("step_{i}").as_str());
            db.insert(row);
        }
        db
    }
    #[test]
    fn test_export_csv_basic() {
        let db = make_test_db();
        let pipe = ExportPipeline::new(ExportFormat::Csv);
        let out = pipe.export(&db);
        assert!(out.contains("step"));
        assert!(out.contains("energy"));
    }
    #[test]
    fn test_export_json_basic() {
        let db = make_test_db();
        let pipe = ExportPipeline::new(ExportFormat::Json);
        let out = pipe.export(&db);
        assert!(out.starts_with('['));
        assert!(out.ends_with(']'));
    }
    #[test]
    fn test_export_hdf5text_basic() {
        let db = make_test_db();
        let pipe = ExportPipeline::new(ExportFormat::Hdf5Text);
        let out = pipe.export(&db);
        assert!(out.contains("# HDF5-like"));
        assert!(out.contains("runs"));
    }
    #[test]
    fn test_export_with_range_filter() {
        let db = make_test_db();
        let mut filter = ExportFilter::new();
        filter.min_value = Some(("energy".into(), 5.0));
        let pipe = ExportPipeline::new(ExportFormat::Csv).with_filter(filter);
        let out = pipe.export(&db);
        let data_lines = out.lines().skip(1).count();
        assert_eq!(data_lines, 2);
    }
    #[test]
    fn test_export_with_column_projection() {
        let db = make_test_db();
        let mut filter = ExportFilter::new();
        filter.columns = vec!["step".into()];
        let pipe = ExportPipeline::new(ExportFormat::Csv).with_filter(filter);
        let out = pipe.export(&db);
        assert!(out.contains("step"));
        assert!(!out.contains("energy"));
    }
    #[test]
    fn test_export_empty_db() {
        let db = SimulationDatabase::new("empty");
        let pipe = ExportPipeline::new(ExportFormat::Csv);
        let out = pipe.export(&db);
        assert!(out.is_empty() || !out.contains('\n').eq(&false));
    }
    #[test]
    fn test_export_filter_max_value() {
        let db = make_test_db();
        let mut filter = ExportFilter::new();
        filter.max_value = Some(("energy".into(), 5.0));
        let pipe = ExportPipeline::new(ExportFormat::Json).with_filter(filter);
        let out = pipe.export(&db);
        let obj_count = out.matches('{').count();
        assert_eq!(obj_count, 1);
    }
    #[test]
    fn test_sim_record_new() {
        let rec = SimulationRecord::new("run_001", 1234.0_f64);
        assert_eq!(rec.name, "run_001");
        assert!((rec.timestamp - 1234.0_f64).abs() < 1e-10);
        assert!(rec.params.is_empty());
        assert!(rec.output_path.is_none());
    }
    #[test]
    fn test_sim_record_set_get_param() {
        let mut rec = SimulationRecord::new("r", 0.0_f64);
        rec.set_param("dt", "0.01");
        rec.set_param("steps", "1000");
        assert_eq!(rec.get_param("dt"), Some("0.01"));
        assert_eq!(rec.get_param("steps"), Some("1000"));
        assert!(rec.get_param("missing").is_none());
    }
    #[test]
    fn test_sim_record_set_output() {
        let path = std::env::temp_dir().join("result.vtk");
        let path_str = path.to_str().unwrap_or("").to_string();
        let mut rec = SimulationRecord::new("r", 0.0_f64);
        rec.set_output(&path_str);
        assert_eq!(rec.output_path, Some(path_str));
    }
    #[test]
    fn test_query_time_range_pass() {
        let rec = SimulationRecord::new("x", 5.0_f64);
        let q = DatabaseQuery::new().with_time_range(0.0_f64, 10.0_f64);
        assert!(q.matches(&rec));
    }
    #[test]
    fn test_query_time_range_fail_before() {
        let rec = SimulationRecord::new("x", -1.0_f64);
        let q = DatabaseQuery::new().with_time_range(0.0_f64, 10.0_f64);
        assert!(!q.matches(&rec));
    }
    #[test]
    fn test_query_time_range_fail_after() {
        let rec = SimulationRecord::new("x", 100.0_f64);
        let q = DatabaseQuery::new().with_time_range(0.0_f64, 10.0_f64);
        assert!(!q.matches(&rec));
    }
    #[test]
    fn test_query_name_prefix_pass() {
        let rec = SimulationRecord::new("run_001", 0.0_f64);
        let q = DatabaseQuery::new().with_name_prefix("run_");
        assert!(q.matches(&rec));
    }
    #[test]
    fn test_query_name_prefix_fail() {
        let rec = SimulationRecord::new("sim_001", 0.0_f64);
        let q = DatabaseQuery::new().with_name_prefix("run_");
        assert!(!q.matches(&rec));
    }
    #[test]
    fn test_query_param_filter_pass() {
        let mut rec = SimulationRecord::new("r", 0.0_f64);
        rec.set_param("solver", "cg");
        let q = DatabaseQuery::new().with_param("solver", "cg");
        assert!(q.matches(&rec));
    }
    #[test]
    fn test_query_param_filter_fail() {
        let mut rec = SimulationRecord::new("r", 0.0_f64);
        rec.set_param("solver", "direct");
        let q = DatabaseQuery::new().with_param("solver", "cg");
        assert!(!q.matches(&rec));
    }
    #[test]
    fn test_query_empty_passes_all() {
        let rec = SimulationRecord::new("anything", 999.0_f64);
        let q = DatabaseQuery::new();
        assert!(q.matches(&rec));
    }
    #[test]
    fn test_record_db_insert_get() {
        let mut db = SimulationRecordDatabase::new();
        db.insert(SimulationRecord::new("r1", 1.0_f64));
        assert!(db.get("r1").is_some());
        assert_eq!(db.count(), 1);
    }
    #[test]
    fn test_record_db_delete() {
        let mut db = SimulationRecordDatabase::new();
        db.insert(SimulationRecord::new("r1", 0.0_f64));
        let removed = db.delete("r1");
        assert!(removed);
        assert!(db.get("r1").is_none());
        assert_eq!(db.count(), 0);
    }
    #[test]
    fn test_record_db_delete_missing() {
        let mut db = SimulationRecordDatabase::new();
        assert!(!db.delete("nonexistent"));
    }
    #[test]
    fn test_record_db_update() {
        let mut db = SimulationRecordDatabase::new();
        db.insert(SimulationRecord::new("r1", 0.0_f64));
        db.get_mut("r1").unwrap().set_param("dt", "0.01");
        assert_eq!(db.get("r1").unwrap().get_param("dt"), Some("0.01"));
    }
    #[test]
    fn test_record_db_query() {
        let mut db = SimulationRecordDatabase::new();
        for i in 0..5 {
            let mut rec = SimulationRecord::new(format!("run_{i}"), i as f64);
            rec.set_param("type", "fluid");
            db.insert(rec);
        }
        let q = DatabaseQuery::new()
            .with_time_range(1.0_f64, 3.0_f64)
            .with_param("type", "fluid");
        let results = db.query(&q);
        assert_eq!(results.len(), 3);
    }
    #[test]
    fn test_record_db_names() {
        let mut db = SimulationRecordDatabase::new();
        db.insert(SimulationRecord::new("a", 0.0_f64));
        db.insert(SimulationRecord::new("b", 0.0_f64));
        let names = db.names();
        assert_eq!(names.len(), 2);
    }
    #[test]
    fn test_record_db_clear() {
        let mut db = SimulationRecordDatabase::new();
        db.insert(SimulationRecord::new("x", 0.0_f64));
        db.clear();
        assert!(db.is_empty());
    }
    #[test]
    fn test_serializer_round_trip_name_timestamp() {
        let ser = DatabaseSerializer::new();
        let rec = SimulationRecord::new("run_test", 42.5_f64);
        let s = ser.serialize(&rec);
        let rec2 = ser.deserialize(&s).expect("deserialize failed");
        assert_eq!(rec2.name, "run_test");
        assert!((rec2.timestamp - 42.5_f64).abs() < 1e-6_f64);
    }
    #[test]
    fn test_serializer_serialize_contains_name() {
        let ser = DatabaseSerializer::new();
        let rec = SimulationRecord::new("my_run", 0.0_f64);
        let s = ser.serialize(&rec);
        assert!(s.contains("my_run"));
        assert!(s.contains("timestamp"));
    }
    #[test]
    fn test_serializer_serialize_all() {
        let ser = DatabaseSerializer::new();
        let r1 = SimulationRecord::new("a", 0.0_f64);
        let r2 = SimulationRecord::new("b", 1.0_f64);
        let out = ser.serialize_all(&[&r1, &r2]);
        assert!(out.starts_with('['));
        assert!(out.ends_with(']'));
    }
    #[test]
    fn test_serializer_output_path() {
        let ser = DatabaseSerializer::new();
        let mut rec = SimulationRecord::new("r", 0.0_f64);
        rec.set_output("/data/result.vtu");
        let s = ser.serialize(&rec);
        assert!(s.contains("/data/result.vtu"));
    }
    #[test]
    fn test_snapshot_catalog_register_get() {
        let mut cat = SnapshotCatalog::new();
        let path = std::env::temp_dir().join("f001.vtk");
        let entry = SnapshotEntry::new("frame_001", 0.1_f64, path.to_str().unwrap_or(""));
        cat.register(entry);
        assert!(cat.get("frame_001").is_some());
        assert_eq!(cat.len(), 1);
    }
    #[test]
    fn test_snapshot_catalog_auto_id() {
        let mut cat = SnapshotCatalog::new();
        let path = std::env::temp_dir().join("snap.vtk");
        let id = cat.register_auto(0.5_f64, path.to_str().unwrap_or(""));
        assert!(id.starts_with("snap_"));
        assert_eq!(cat.len(), 1);
    }
    #[test]
    fn test_snapshot_catalog_remove() {
        let mut cat = SnapshotCatalog::new();
        cat.register(SnapshotEntry::new("s1", 0.0_f64, ""));
        assert!(cat.remove("s1"));
        assert!(cat.is_empty());
        assert!(!cat.remove("s1"));
    }
    #[test]
    fn test_snapshot_catalog_range_query() {
        let mut cat = SnapshotCatalog::new();
        for i in 0..10 {
            cat.register(SnapshotEntry::new(format!("s{i}"), i as f64 * 0.1_f64, ""));
        }
        let results = cat.range_query(0.2_f64, 0.5_f64);
        assert_eq!(results.len(), 4);
    }
    #[test]
    fn test_snapshot_catalog_query_by_tag() {
        let mut cat = SnapshotCatalog::new();
        let mut s1 = SnapshotEntry::new("s1", 0.0_f64, "");
        s1.add_tag("checkpoint");
        let mut s2 = SnapshotEntry::new("s2", 1.0_f64, "");
        s2.add_tag("debug");
        cat.register(s1);
        cat.register(s2);
        let cps = cat.query_by_tag("checkpoint");
        assert_eq!(cps.len(), 1);
    }
    #[test]
    fn test_snapshot_catalog_mark_loaded() {
        let mut cat = SnapshotCatalog::new();
        cat.register(SnapshotEntry::new("s1", 0.0_f64, ""));
        cat.register(SnapshotEntry::new("s2", 5.0_f64, ""));
        cat.mark_range_loaded(0.0_f64, 3.0_f64);
        assert!(cat.get("s1").unwrap().loaded);
        assert!(!cat.get("s2").unwrap().loaded);
    }
    #[test]
    fn test_snapshot_catalog_unloaded_ids() {
        let mut cat = SnapshotCatalog::new();
        cat.register(SnapshotEntry::new("s1", 0.0_f64, ""));
        cat.register(SnapshotEntry::new("s2", 1.0_f64, ""));
        cat.get_mut("s1").unwrap().mark_loaded();
        let unloaded = cat.unloaded_ids();
        assert_eq!(unloaded.len(), 1);
        assert_eq!(unloaded[0], "s2");
    }
    #[test]
    fn test_snapshot_catalog_total_file_size() {
        let mut cat = SnapshotCatalog::new();
        let mut s1 = SnapshotEntry::new("s1", 0.0_f64, "");
        s1.file_size = 1024;
        let mut s2 = SnapshotEntry::new("s2", 1.0_f64, "");
        s2.file_size = 2048;
        cat.register(s1);
        cat.register(s2);
        assert_eq!(cat.total_file_size(), 3072);
    }
    #[test]
    fn test_snapshot_catalog_sorted_by_time() {
        let mut cat = SnapshotCatalog::new();
        cat.register(SnapshotEntry::new("c", 3.0_f64, ""));
        cat.register(SnapshotEntry::new("a", 1.0_f64, ""));
        cat.register(SnapshotEntry::new("b", 2.0_f64, ""));
        let sorted = cat.sorted_by_time();
        assert_eq!(sorted.len(), 3);
        assert!((sorted[0].sim_time - 1.0_f64).abs() < 1e-10);
        assert!((sorted[2].sim_time - 3.0_f64).abs() < 1e-10);
    }
}
