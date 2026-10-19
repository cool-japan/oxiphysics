// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the material JSON writer and reader: Unicode text, structural
//! characters in names and keys, exact number round trips, non-finite values,
//! deterministic output and the rejection of malformed records.

use super::{JsonMaterialDb, MaterialProperty, MaterialRecord};

fn write(record: &MaterialRecord) -> String {
    JsonMaterialDb::write_json_material(record)
}

fn parse(json: &str) -> Option<MaterialRecord> {
    JsonMaterialDb::parse_json_material(json)
}

fn round_trip(record: &MaterialRecord) -> MaterialRecord {
    let json = write(record);
    parse(&json).unwrap_or_else(|| panic!("the record did not parse back:\n{json}"))
}

fn assert_same_record(parsed: &MaterialRecord, original: &MaterialRecord) {
    assert_eq!(parsed.name, original.name);
    assert_eq!(parsed.category, original.category);
    assert_eq!(parsed.source, original.source);
    assert_eq!(parsed.properties, original.properties);
}

#[test]
fn non_ascii_names_and_property_keys_round_trip_without_panic() {
    let mut record = MaterialRecord::new("Ti-6Al-4V (α+β)", "Métal>鉄鋼");
    record.source = "Källa: 日本工業規格 😀".to_string();
    record.set_scalar("ρ (kg/m³)", 4430.0);
    record.set_range("σ_y (MPa) — 降伏", 800.0, 1000.0);
    record.properties.insert(
        "κ[W/m·K]".to_string(),
        MaterialProperty::Matrix(vec![1.5, 2.5, 3.5]),
    );
    let parsed = round_trip(&record);
    assert_same_record(&parsed, &record);
    assert_eq!(parsed.get_value("ρ (kg/m³)"), Some(4430.0));
}

#[test]
fn structural_characters_in_text_do_not_confuse_the_parser() {
    let mut record = MaterialRecord::new("a}\"properties\": {", "x{y}z");
    record.source = "\"type\": \"scalar\"".to_string();
    record.set_scalar("}{", 1.0);
    record.set_scalar("\"quoted\\key\"", 2.0);
    // Keys equal to the words the format uses.
    record.set_scalar("scalar", 3.0);
    record.set_scalar("range", 4.0);
    record.set_scalar("matrix", 5.0);
    record
        .properties
        .insert("type".to_string(), MaterialProperty::Range(6.0, 7.0));
    record.set_scalar("value", 8.0);
    assert_same_record(&round_trip(&record), &record);
}

#[test]
fn control_characters_and_unicode_edge_cases_round_trip() {
    let mut text: String = (0u8..0x20).map(char::from).collect();
    text.push_str("\"\\/\u{7f}\u{80}\u{2028}\u{2029}\u{fffd}\u{ffff}\u{10000}\u{10ffff}é😀");
    let mut record = MaterialRecord::new(&text, &text);
    record.source = text.clone();
    record.set_scalar(&text, 1.0);
    let json = write(&record);
    assert!(
        !json.chars().any(|c| (c as u32) < 0x20 && c != '\n'),
        "control characters are escaped in the output"
    );
    assert_same_record(&parse(&json).unwrap(), &record);
}

#[test]
fn numbers_round_trip_bit_exactly() {
    let values = [
        0.0,
        -0.0,
        1.0,
        7850.0,
        7850.123456789,
        1.23e-7,
        1e-300,
        5e-324,
        f64::MIN_POSITIVE,
        f64::MAX,
        -2.5e10,
        0.1 + 0.2,
        std::f64::consts::PI,
    ];
    let mut record = MaterialRecord::new("exact", "Test");
    for (i, v) in values.iter().enumerate() {
        record.set_scalar(&format!("s{i}"), *v);
    }
    record.set_range("r", 1.23e-7, 4.56e8);
    record
        .properties
        .insert("m".to_string(), MaterialProperty::Matrix(values.to_vec()));
    let parsed = round_trip(&record);
    for (i, v) in values.iter().enumerate() {
        match parsed.properties.get(&format!("s{i}")) {
            Some(MaterialProperty::Scalar(got)) => {
                assert_eq!(got.to_bits(), v.to_bits(), "scalar {i}: {got} vs {v}");
            }
            other => panic!("scalar {i} came back as {other:?}"),
        }
    }
    match parsed.properties.get("r") {
        Some(MaterialProperty::Range(lo, hi)) => {
            assert_eq!(lo.to_bits(), 1.23e-7_f64.to_bits());
            assert_eq!(hi.to_bits(), 4.56e8_f64.to_bits());
        }
        other => panic!("range came back as {other:?}"),
    }
    match parsed.properties.get("m") {
        Some(MaterialProperty::Matrix(got)) => {
            assert_eq!(got.len(), values.len());
            for (g, v) in got.iter().zip(values) {
                assert_eq!(g.to_bits(), v.to_bits());
            }
        }
        other => panic!("matrix came back as {other:?}"),
    }
}

#[test]
fn small_values_are_not_rounded_to_six_decimals() {
    let mut record = MaterialRecord::new("alloy", "Metal");
    record.set_scalar("thermal_expansion", 1.234567e-5);
    let parsed = round_trip(&record);
    assert_eq!(parsed.get_value("thermal_expansion"), Some(1.234567e-5));
}

#[test]
fn non_finite_values_round_trip_and_the_text_stays_valid_json() {
    let mut record = MaterialRecord::new("limits", "Test");
    record.set_scalar("nan", f64::NAN);
    record.set_range("open", 0.0, f64::INFINITY);
    record.set_range("below", f64::NEG_INFINITY, 1.0);
    record.properties.insert(
        "m".to_string(),
        MaterialProperty::Matrix(vec![f64::NEG_INFINITY, f64::NAN, 1.0, f64::INFINITY]),
    );
    let json = write(&record);
    assert!(
        serde_json::from_str::<serde_json::Value>(&json).is_ok(),
        "NaN and infinities must not appear as bare tokens:\n{json}"
    );
    let parsed = parse(&json).unwrap();
    match parsed.properties.get("nan") {
        Some(MaterialProperty::Scalar(v)) => assert!(v.is_nan()),
        other => panic!("{other:?}"),
    }
    assert_eq!(
        parsed.properties.get("open"),
        Some(&MaterialProperty::Range(0.0, f64::INFINITY))
    );
    assert_eq!(
        parsed.properties.get("below"),
        Some(&MaterialProperty::Range(f64::NEG_INFINITY, 1.0))
    );
    match parsed.properties.get("m") {
        Some(MaterialProperty::Matrix(v)) => {
            assert_eq!(v.len(), 4);
            assert_eq!(v[0], f64::NEG_INFINITY);
            assert!(v[1].is_nan());
            assert_eq!(v[2], 1.0);
            assert_eq!(v[3], f64::INFINITY);
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn writer_output_has_the_documented_layout() {
    let mut record = MaterialRecord::new("Steel", "Metal");
    record.source = "Matweb".to_string();
    record.set_scalar("density", 7850.0);
    record.set_range("yield", 250.0, 350.0);
    record
        .properties
        .insert("k".to_string(), MaterialProperty::Matrix(vec![1.0, 2.5]));
    let expected = "{\n  \"name\": \"Steel\",\n  \"category\": \"Metal\",\n  \"source\": \"Matweb\",\n  \
         \"properties\": {\n    \"density\": {\"type\": \"scalar\", \"value\": 7850.0},\n    \
         \"k\": {\"type\": \"matrix\", \"values\": [1.0, 2.5]},\n    \
         \"yield\": {\"type\": \"range\", \"lo\": 250.0, \"hi\": 350.0}\n  }\n}";
    assert_eq!(write(&record), expected);
    let empty = MaterialRecord::new("Void", "");
    assert_eq!(
        write(&empty),
        "{\n  \"name\": \"Void\",\n  \"category\": \"\",\n  \"source\": \"\",\n  \
         \"properties\": {\n  }\n}"
    );
    assert_same_record(&round_trip(&empty), &empty);
}

#[test]
fn output_does_not_depend_on_hash_map_order() {
    let build = |reverse: bool| {
        let mut record = MaterialRecord::new("many", "Test");
        let mut keys: Vec<usize> = (0..40).collect();
        if reverse {
            keys.reverse();
        }
        for k in keys {
            record.set_scalar(&format!("property_{k:02}"), k as f64);
        }
        record
    };
    let first = write(&build(false));
    assert_eq!(first, write(&build(true)));
    assert_eq!(first, write(&build(false)));
    // Keys are listed in sorted order.
    let positions: Vec<usize> = (0..40)
        .map(|k| first.find(&format!("\"property_{k:02}\"")).unwrap())
        .collect();
    assert!(positions.windows(2).all(|w| w[0] < w[1]));
}

#[test]
fn optional_members_default_to_empty_and_unknown_members_are_ignored() {
    let parsed = parse(r#"{"name": "x", "properties": {}}"#).unwrap();
    assert_eq!(
        (
            parsed.name.as_str(),
            parsed.category.as_str(),
            parsed.source.as_str()
        ),
        ("x", "", "")
    );
    assert!(parsed.properties.is_empty());
    let parsed =
        parse(r#"{"extra": [1, 2, {"a": null}], "name": "x", "properties": {}, "id": 7}"#).unwrap();
    assert_eq!(parsed.name, "x");
}

#[test]
fn any_json_layout_is_accepted() {
    let json = r#"{"properties":{"density":{"value":7850,"type":"scalar"},
        "yield":{"hi":3.5E2,"lo":2.5e+2,"type":"range"},
        "m":{"values":[1,2,3],"type":"matrix"}},"source":"S","category":"C","name":"N"}"#;
    let parsed = parse(json).unwrap();
    assert_eq!(
        (
            parsed.name.as_str(),
            parsed.category.as_str(),
            parsed.source.as_str()
        ),
        ("N", "C", "S")
    );
    assert_eq!(
        parsed.properties.get("density"),
        Some(&MaterialProperty::Scalar(7850.0))
    );
    assert_eq!(
        parsed.properties.get("yield"),
        Some(&MaterialProperty::Range(250.0, 350.0))
    );
    assert_eq!(
        parsed.properties.get("m"),
        Some(&MaterialProperty::Matrix(vec![1.0, 2.0, 3.0]))
    );
}

#[test]
fn property_type_comes_from_the_type_member_only() {
    let json = r#"{"name": "x", "properties": {
        "a": {"type": "scalar", "value": 1, "note": "range matrix", "lo": 5, "hi": 6}}}"#;
    let parsed = parse(json).unwrap();
    assert_eq!(
        parsed.properties.get("a"),
        Some(&MaterialProperty::Scalar(1.0))
    );
}

#[test]
fn json_escapes_in_the_input_are_decoded() {
    let json = "{\"name\": \"a\\\"b\\\\c\\/d\\n\\t\\u00e5\\u9444\", \"properties\": {}}";
    assert_eq!(parse(json).unwrap().name, "a\"b\\c/d\n\tå鑄");
    // Surrogate pair for U+1F600.
    let pair = "{\"name\": \"\\ud83d\\ude00!\", \"properties\": {}}";
    assert_eq!(parse(pair).unwrap().name, "\u{1F600}!");
}

#[test]
fn invalid_escapes_and_unterminated_strings_are_rejected() {
    for bad_name in [
        r"\q", r"\u12", r"\u12zz", r"\ud83d", r"\ud83dx", r"\ude00", r"\ud83dA",
    ] {
        let json = format!("{{\"name\": \"{bad_name}\", \"properties\": {{}}}}");
        assert!(parse(&json).is_none(), "{json}");
    }
    for truncated in [
        "{\"name\": \"abc\\",
        "{\"name\": \"abc",
        "{\"name\": \"",
        "{\"name\": \"鋼\\",
        "{\"name\":",
        "{\"name\": \"Ab\\",
    ] {
        assert!(parse(truncated).is_none(), "{truncated}");
    }
}

#[test]
fn a_record_cut_short_anywhere_is_rejected_and_never_panics() {
    let mut record = MaterialRecord::new("Ståhl 鋼 (α+β)", "ρ>σ");
    record.source = "Källa 😀".to_string();
    record.set_scalar("ρ (kg/m³)", 7850.0);
    record.set_range("σ_y", 250.0, 350.0);
    record
        .properties
        .insert("κ".to_string(), MaterialProperty::Matrix(vec![1.0, 2.0]));
    let json = write(&record);
    for (end, _) in json.char_indices() {
        assert!(
            parse(&json[..end]).is_none(),
            "the first {end} bytes parsed as a record"
        );
    }
    assert_same_record(&parse(&json).unwrap(), &record);
}

#[test]
fn malformed_records_are_rejected_instead_of_half_parsed() {
    let bad = [
        "",
        "null",
        "[]",
        "42",
        "\"name\"",
        "{}",
        r#"{"properties": {}}"#,
        r#"{"name": 5, "properties": {}}"#,
        r#"{"name": "x"}"#,
        r#"{"name": "x", "properties": []}"#,
        r#"{"name": "x", "properties": {"a": 5}}"#,
        r#"{"name": "x", "properties": {"a": {}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": 1}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "tensor", "value": 1}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "scalar"}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "scalar", "value": "abc"}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "scalar", "value": null}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "scalar", "value": true}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "range", "lo": 1}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "range", "lo": 1, "hi": "oops"}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "matrix"}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "matrix", "values": 3}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "matrix", "values": [1, "x"]}}}"#,
        r#"{"name": "x", "properties": {"a": {"type": "matrix", "values": [1, null]}}}"#,
        r#"{"name": "x", "category": 5, "properties": {}}"#,
        r#"{"name": "x", "source": null, "properties": {}}"#,
        r#"{"name": "x", "properties": {}} trailing"#,
        r#"{"name": "x", "properties": {}}{"name": "y", "properties": {}}"#,
        // One bad property poisons the record: nothing is silently dropped.
        r#"{"name": "x", "properties": {"ok": {"type": "scalar", "value": 1},
            "bad": {"type": "scalar", "value": "n/a"}}}"#,
    ];
    for json in bad {
        assert!(parse(json).is_none(), "{json}");
    }
}

#[test]
fn deeply_nested_input_is_rejected_without_overflowing_the_stack() {
    let depth = 100_000;
    assert!(parse(&"[".repeat(depth)).is_none());
    let nested = format!(
        "{{\"name\": \"x\", \"properties\": {{\"a\": {{\"type\": \"scalar\", \"value\": {}1{}}}}}}}",
        "[".repeat(depth),
        "]".repeat(depth)
    );
    assert!(parse(&nested).is_none());
}
