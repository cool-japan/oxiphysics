// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! JSON restart encoding.
//!
//! One JSON object per file; one member per line, arrays on a single line:
//!
//! ```text
//! {
//!   "format": "oxiphysics-restart",
//!   "format_version": 2,
//!   "version": "1.0",
//!   "timestamp": 1234567890,
//!   "step": 42,
//!   "time": 4.2,
//!   "crate": "oxiphysics-md",
//!   "description": "unit test",
//!   "n_particles": 3,
//!   "positions": [[0.0,0.0,0.0],[1.0,0.5,0.0],[2.0,1.0,0.0]],
//!   "velocities": [...], "forces": [...], "masses": [...], "types": [...],
//!   "box_matrix": [[10.0,0.0,0.0],[0.0,10.0,0.0],[0.0,0.0,10.0]],
//!   "extra_scalars": [{"name": "charge", "values": [0.1, -0.2]}],
//!   "extra_vectors": [{"name": "spin", "values": [[0.0,0.0,1.0]]}]
//! }
//! ```
//!
//! Floating-point numbers are written with the shortest decimal form that
//! parses back to the identical `f64`, and the reader is built with serde_json's
//! `float_roundtrip` feature, so values survive bit-exactly (including `-0.0`
//! and subnormals). JSON has no spelling for NaN and infinities, so those are
//! written as the strings `"NaN"`, `"Infinity"` and `"-Infinity"` and accepted
//! in the same form on read. Strings are escaped by serde_json.
//!
//! On read, `forces`, `velocities`, `masses`, `types`, `box_matrix` and the
//! extras may be absent (older files omit some of them); unknown members are
//! ignored.

use serde_json::{Map, Value};

use super::{FORMAT_VERSION, RestartData, RestartMetadata};

/// Value of the `"format"` member.
const FORMAT_TAG: &str = "oxiphysics-restart";

// ── Encoder ───────────────────────────────────────────────────────────────────

/// Render an `f64` as a JSON token.
fn json_f64(v: f64) -> String {
    if v.is_nan() {
        "\"NaN\"".to_string()
    } else if v == f64::INFINITY {
        "\"Infinity\"".to_string()
    } else if v == f64::NEG_INFINITY {
        "\"-Infinity\"".to_string()
    } else {
        // `{:?}` is the shortest round-trip form and is valid JSON for every
        // finite value (never `+`, a leading `.`, or a trailing `.`).
        format!("{v:?}")
    }
}

fn json_vec3(v: &[f64; 3]) -> String {
    format!("[{},{},{}]", json_f64(v[0]), json_f64(v[1]), json_f64(v[2]))
}

fn json_string(s: &str) -> Result<String, String> {
    serde_json::to_string(s).map_err(|e| format!("cannot encode string as JSON: {e}"))
}

fn json_array<T>(items: &[T], render: impl Fn(&T) -> String) -> String {
    let body: Vec<String> = items.iter().map(render).collect();
    format!("[{}]", body.join(","))
}

fn json_named_arrays<T>(
    arrays: &[(String, Vec<T>)],
    render: impl Fn(&T) -> String,
) -> Result<String, String> {
    let mut entries = Vec::with_capacity(arrays.len());
    for (name, values) in arrays {
        entries.push(format!(
            "{{\"name\": {}, \"values\": {}}}",
            json_string(name)?,
            json_array(values, &render)
        ));
    }
    Ok(format!("[{}]", entries.join(", ")))
}

/// Encode `data` as a JSON document. Fails if [`RestartData::validate`] fails.
pub(super) fn encode(data: &RestartData) -> Result<String, String> {
    data.validate()?;
    let meta = &data.metadata;
    let members: Vec<(&str, String)> = vec![
        ("format", json_string(FORMAT_TAG)?),
        ("format_version", FORMAT_VERSION.to_string()),
        ("version", json_string(&meta.version)?),
        ("timestamp", meta.timestamp.to_string()),
        ("step", meta.step.to_string()),
        ("time", json_f64(meta.time)),
        ("crate", json_string(&meta.crate_name)?),
        ("description", json_string(&meta.description)?),
        ("n_particles", data.n_particles().to_string()),
        ("positions", json_array(&data.positions, json_vec3)),
        ("velocities", json_array(&data.velocities, json_vec3)),
        ("forces", json_array(&data.forces, json_vec3)),
        ("masses", json_array(&data.masses, |m| json_f64(*m))),
        ("types", json_array(&data.types, |t| t.to_string())),
        ("box_matrix", json_array(&data.box_matrix, json_vec3)),
        (
            "extra_scalars",
            json_named_arrays(&data.extra_scalars, |v| json_f64(*v))?,
        ),
        (
            "extra_vectors",
            json_named_arrays(&data.extra_vectors, json_vec3)?,
        ),
    ];
    let mut out = String::from("{\n");
    for (i, (key, value)) in members.iter().enumerate() {
        let comma = if i + 1 < members.len() { "," } else { "" };
        out.push_str(&format!("  \"{key}\": {value}{comma}\n"));
    }
    out.push_str("}\n");
    Ok(out)
}

// ── Decoder ───────────────────────────────────────────────────────────────────

fn value_f64(value: &Value, what: &str) -> Result<f64, String> {
    match value {
        Value::Number(n) => n
            .as_f64()
            .ok_or_else(|| format!("restart JSON: {what} is not a finite number")),
        Value::String(s) => match s.as_str() {
            "NaN" => Ok(f64::NAN),
            "Infinity" | "inf" => Ok(f64::INFINITY),
            "-Infinity" | "-inf" => Ok(f64::NEG_INFINITY),
            other => Err(format!(
                "restart JSON: {what} is the string `{other}`, expected a number, \"NaN\", \"Infinity\" or \"-Infinity\""
            )),
        },
        other => Err(format!(
            "restart JSON: {what} must be a number, found {}",
            json_kind(other)
        )),
    }
}

fn json_kind(value: &Value) -> &'static str {
    match value {
        Value::Null => "null",
        Value::Bool(_) => "a boolean",
        Value::Number(_) => "a number",
        Value::String(_) => "a string",
        Value::Array(_) => "an array",
        Value::Object(_) => "an object",
    }
}

fn value_vec3(value: &Value, what: &str) -> Result<[f64; 3], String> {
    let items = value.as_array().ok_or_else(|| {
        format!(
            "restart JSON: {what} must be an array of 3 numbers, found {}",
            json_kind(value)
        )
    })?;
    match items.as_slice() {
        [x, y, z] => Ok([
            value_f64(x, &format!("{what}[0]"))?,
            value_f64(y, &format!("{what}[1]"))?,
            value_f64(z, &format!("{what}[2]"))?,
        ]),
        _ => Err(format!(
            "restart JSON: {what} must have exactly 3 entries, found {}",
            items.len()
        )),
    }
}

fn member_array<'a>(
    obj: &'a Map<String, Value>,
    key: &str,
) -> Result<Option<&'a Vec<Value>>, String> {
    match obj.get(key) {
        None => Ok(None),
        Some(Value::Array(items)) => Ok(Some(items)),
        Some(other) => Err(format!(
            "restart JSON: `{key}` must be an array, found {}",
            json_kind(other)
        )),
    }
}

fn member_vec3_array(obj: &Map<String, Value>, key: &str) -> Result<Vec<[f64; 3]>, String> {
    match member_array(obj, key)? {
        None => Ok(Vec::new()),
        Some(items) => items
            .iter()
            .enumerate()
            .map(|(i, v)| value_vec3(v, &format!("{key}[{i}]")))
            .collect(),
    }
}

fn member_f64_array(obj: &Map<String, Value>, key: &str) -> Result<Vec<f64>, String> {
    match member_array(obj, key)? {
        None => Ok(Vec::new()),
        Some(items) => items
            .iter()
            .enumerate()
            .map(|(i, v)| value_f64(v, &format!("{key}[{i}]")))
            .collect(),
    }
}

fn member_u32_array(obj: &Map<String, Value>, key: &str) -> Result<Vec<u32>, String> {
    match member_array(obj, key)? {
        None => Ok(Vec::new()),
        Some(items) => items
            .iter()
            .enumerate()
            .map(|(i, v)| {
                v.as_u64()
                    .and_then(|n| u32::try_from(n).ok())
                    .ok_or_else(|| {
                        format!("restart JSON: {key}[{i}] must be an unsigned 32-bit integer")
                    })
            })
            .collect(),
    }
}

fn member_string(obj: &Map<String, Value>, key: &str) -> Result<Option<String>, String> {
    match obj.get(key) {
        None => Ok(None),
        Some(Value::String(s)) => Ok(Some(s.clone())),
        Some(other) => Err(format!(
            "restart JSON: `{key}` must be a string, found {}",
            json_kind(other)
        )),
    }
}

fn member_u64(obj: &Map<String, Value>, key: &str) -> Result<Option<u64>, String> {
    match obj.get(key) {
        None => Ok(None),
        Some(v) => v
            .as_u64()
            .map(Some)
            .ok_or_else(|| format!("restart JSON: `{key}` must be an unsigned integer, found {v}")),
    }
}

fn member_f64(obj: &Map<String, Value>, key: &str) -> Result<Option<f64>, String> {
    match obj.get(key) {
        None => Ok(None),
        Some(v) => value_f64(v, &format!("`{key}`")).map(Some),
    }
}

fn member_named_arrays<T>(
    obj: &Map<String, Value>,
    key: &str,
    convert: impl Fn(&Value, &str) -> Result<T, String>,
) -> Result<Vec<(String, Vec<T>)>, String> {
    let Some(entries) = member_array(obj, key)? else {
        return Ok(Vec::new());
    };
    let mut out = Vec::with_capacity(entries.len());
    for (i, entry) in entries.iter().enumerate() {
        let entry = entry.as_object().ok_or_else(|| {
            format!(
                "restart JSON: {key}[{i}] must be an object with `name` and `values`, found {}",
                json_kind(entry)
            )
        })?;
        let name = match entry.get("name") {
            Some(Value::String(s)) => s.clone(),
            _ => {
                return Err(format!(
                    "restart JSON: {key}[{i}] needs a string `name` member"
                ));
            }
        };
        let values = match entry.get("values") {
            Some(Value::Array(items)) => items
                .iter()
                .enumerate()
                .map(|(j, v)| convert(v, &format!("{key}[{i}].values[{j}]")))
                .collect::<Result<Vec<T>, String>>()?,
            _ => {
                return Err(format!(
                    "restart JSON: {key}[{i}] needs a `values` array member"
                ));
            }
        };
        out.push((name, values));
    }
    Ok(out)
}

fn member_box(obj: &Map<String, Value>) -> Result<[[f64; 3]; 3], String> {
    let rows = member_vec3_array(obj, "box_matrix")?;
    match rows.as_slice() {
        [] => Ok([[0.0; 3]; 3]),
        [a, b, c] => Ok([*a, *b, *c]),
        other => Err(format!(
            "restart JSON: `box_matrix` must have exactly 3 rows, found {}",
            other.len()
        )),
    }
}

/// Decode a JSON restart document.
pub(super) fn decode(text: &str) -> Result<RestartData, String> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let root: Value = serde_json::from_str(text).map_err(|e| format!("restart JSON: {e}"))?;
    let obj = root.as_object().ok_or_else(|| {
        format!(
            "restart JSON: the top-level value must be an object, found {}",
            json_kind(&root)
        )
    })?;

    if let Some(tag) = member_string(obj, "format")?
        && tag != FORMAT_TAG
    {
        return Err(format!(
            "restart JSON: unexpected `format` \"{tag}\" (expected \"{FORMAT_TAG}\")"
        ));
    }
    if let Some(version) = member_u64(obj, "format_version")?
        && (version == 0 || version > u64::from(FORMAT_VERSION))
    {
        return Err(format!(
            "restart JSON: unsupported format_version {version} (this reader supports 1..={FORMAT_VERSION})"
        ));
    }

    let metadata = RestartMetadata {
        version: member_string(obj, "version")?.unwrap_or_default(),
        timestamp: member_u64(obj, "timestamp")?.unwrap_or(0),
        step: member_u64(obj, "step")?.unwrap_or(0),
        time: member_f64(obj, "time")?.unwrap_or(0.0),
        crate_name: member_string(obj, "crate")?.unwrap_or_default(),
        description: member_string(obj, "description")?.unwrap_or_default(),
    };

    let positions = member_vec3_array(obj, "positions")?;
    let velocities = member_vec3_array(obj, "velocities")?;
    let forces = member_vec3_array(obj, "forces")?;
    let masses = member_f64_array(obj, "masses")?;
    let types = member_u32_array(obj, "types")?;
    let box_matrix = member_box(obj)?;
    let extra_scalars = member_named_arrays(obj, "extra_scalars", value_f64)?;
    let extra_vectors = member_named_arrays(obj, "extra_vectors", value_vec3)?;

    if let Some(declared) = member_u64(obj, "n_particles")?
        && declared != positions.len() as u64
    {
        return Err(format!(
            "restart JSON: `n_particles` is {declared} but `positions` holds {} entries",
            positions.len()
        ));
    }

    let data = RestartData {
        metadata,
        positions,
        velocities,
        forces,
        masses,
        types,
        box_matrix,
        extra_scalars,
        extra_vectors,
    };
    data.validate().map_err(|e| format!("restart JSON: {e}"))?;
    Ok(data)
}
