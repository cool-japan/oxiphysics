#![allow(clippy::manual_strip, clippy::should_implement_trait)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simple hand-rolled JSON serialization for physics data (no serde dependency).

#![allow(dead_code)]

/// A JSON value that covers all standard JSON types.
#[derive(Debug, Clone, PartialEq)]
pub enum JsonValue {
    /// JSON `null`.
    Null,
    /// JSON boolean.
    Bool(bool),
    /// JSON number (stored as f64).
    Number(f64),
    /// JSON string.
    Str(String),
    /// JSON array.
    Array(Vec<JsonValue>),
    /// JSON object as ordered key-value pairs.
    Object(Vec<(String, JsonValue)>),
}

impl JsonValue {
    /// Serialize this value to a JSON string.
    pub fn to_json_string(&self) -> String {
        match self {
            JsonValue::Null => "null".to_string(),
            JsonValue::Bool(b) => b.to_string(),
            JsonValue::Number(n) => {
                if n.is_nan() {
                    "null".to_string()
                } else if n.is_infinite() {
                    if *n > 0.0 {
                        "1e308".to_string()
                    } else {
                        "-1e308".to_string()
                    }
                } else if n.fract() == 0.0 && n.abs() < 1e15 {
                    format!("{}", *n as i64)
                } else {
                    format!("{}", n)
                }
            }
            JsonValue::Str(s) => {
                let escaped = s
                    .replace('\\', "\\\\")
                    .replace('"', "\\\"")
                    .replace('\n', "\\n")
                    .replace('\r', "\\r")
                    .replace('\t', "\\t");
                format!("\"{}\"", escaped)
            }
            JsonValue::Array(arr) => {
                let parts: Vec<String> = arr.iter().map(|v| v.to_json_string()).collect();
                format!("[{}]", parts.join(","))
            }
            JsonValue::Object(pairs) => {
                let parts: Vec<String> = pairs
                    .iter()
                    .map(|(k, v)| format!("\"{}\":{}", k, v.to_json_string()))
                    .collect();
                format!("{{{}}}", parts.join(","))
            }
        }
    }

    /// Serialize to a pretty-printed JSON string with indentation.
    pub fn to_json_pretty(&self, indent: usize) -> String {
        self.to_json_pretty_inner(indent, 0)
    }

    fn to_json_pretty_inner(&self, indent: usize, depth: usize) -> String {
        let pad = " ".repeat(indent * depth);
        let pad_inner = " ".repeat(indent * (depth + 1));
        match self {
            JsonValue::Null | JsonValue::Bool(_) | JsonValue::Number(_) | JsonValue::Str(_) => {
                self.to_json_string()
            }
            JsonValue::Array(arr) => {
                if arr.is_empty() {
                    return "[]".to_string();
                }
                let parts: Vec<String> = arr
                    .iter()
                    .map(|v| format!("{}{}", pad_inner, v.to_json_pretty_inner(indent, depth + 1)))
                    .collect();
                format!("[\n{}\n{}]", parts.join(",\n"), pad)
            }
            JsonValue::Object(pairs) => {
                if pairs.is_empty() {
                    return "{}".to_string();
                }
                let parts: Vec<String> = pairs
                    .iter()
                    .map(|(k, v)| {
                        format!(
                            "{}\"{}\": {}",
                            pad_inner,
                            k,
                            v.to_json_pretty_inner(indent, depth + 1)
                        )
                    })
                    .collect();
                format!("{{\n{}\n{}}}", parts.join(",\n"), pad)
            }
        }
    }

    /// Parse a JSON string into a `JsonValue`.
    pub fn from_str(s: &str) -> Result<Self, String> {
        let s = s.trim();
        let (val, rest) = parse_value(s)?;
        if rest.trim().is_empty() {
            Ok(val)
        } else {
            Err(format!("Unexpected trailing input: {:?}", rest.trim()))
        }
    }

    /// If this is a `Number`, return its `f64` value.
    pub fn as_f64(&self) -> Option<f64> {
        if let JsonValue::Number(n) = self {
            Some(*n)
        } else {
            None
        }
    }

    /// If this is a `Bool`, return its value.
    pub fn as_bool(&self) -> Option<bool> {
        if let JsonValue::Bool(b) = self {
            Some(*b)
        } else {
            None
        }
    }

    /// If this is a `Str`, return a reference to the inner string.
    pub fn as_str(&self) -> Option<&str> {
        if let JsonValue::Str(s) = self {
            Some(s.as_str())
        } else {
            None
        }
    }

    /// If this is an `Array`, return a reference to the inner `Vec`.
    pub fn as_array(&self) -> Option<&Vec<JsonValue>> {
        if let JsonValue::Array(arr) = self {
            Some(arr)
        } else {
            None
        }
    }

    /// If this is an `Object`, return a reference to the inner pairs.
    pub fn as_object(&self) -> Option<&Vec<(String, JsonValue)>> {
        if let JsonValue::Object(pairs) = self {
            Some(pairs)
        } else {
            None
        }
    }

    /// Look up a key in an `Object`.  Returns `None` if not an object or key absent.
    pub fn get(&self, key: &str) -> Option<&JsonValue> {
        if let JsonValue::Object(pairs) = self {
            pairs.iter().find(|(k, _)| k == key).map(|(_, v)| v)
        } else {
            None
        }
    }

    /// Index into an `Array`.  Returns `None` if not an array or index out of bounds.
    pub fn get_index(&self, i: usize) -> Option<&JsonValue> {
        if let JsonValue::Array(arr) = self {
            arr.get(i)
        } else {
            None
        }
    }

    /// Navigate a dot-separated path like "a.b.c" through nested objects.
    pub fn get_path(&self, path: &str) -> Option<&JsonValue> {
        let mut current = self;
        for key in path.split('.') {
            current = current.get(key)?;
        }
        Some(current)
    }

    /// Check if this value is null.
    pub fn is_null(&self) -> bool {
        matches!(self, JsonValue::Null)
    }

    /// Return the number of elements in an array or object, or 0 otherwise.
    pub fn len(&self) -> usize {
        match self {
            JsonValue::Array(arr) => arr.len(),
            JsonValue::Object(pairs) => pairs.len(),
            _ => 0,
        }
    }

    /// Return true if this value is an empty array, empty object, or null.
    pub fn is_empty(&self) -> bool {
        match self {
            JsonValue::Null => true,
            JsonValue::Array(arr) => arr.is_empty(),
            JsonValue::Object(pairs) => pairs.is_empty(),
            _ => false,
        }
    }

    /// Merge two JSON objects. Keys from `other` override keys in `self`.
    /// Non-object values return `other`.
    pub fn merge(&self, other: &JsonValue) -> JsonValue {
        match (self, other) {
            (JsonValue::Object(base), JsonValue::Object(overlay)) => {
                let mut merged = base.clone();
                for (key, val) in overlay {
                    if let Some(pos) = merged.iter().position(|(k, _)| k == key) {
                        merged[pos] = (key.clone(), merged[pos].1.merge(val));
                    } else {
                        merged.push((key.clone(), val.clone()));
                    }
                }
                JsonValue::Object(merged)
            }
            _ => other.clone(),
        }
    }

    /// Return all keys if this is an object.
    pub fn keys(&self) -> Vec<&str> {
        if let JsonValue::Object(pairs) = self {
            pairs.iter().map(|(k, _)| k.as_str()).collect()
        } else {
            Vec::new()
        }
    }

    /// Return all values if this is an object.
    pub fn values(&self) -> Vec<&JsonValue> {
        if let JsonValue::Object(pairs) = self {
            pairs.iter().map(|(_, v)| v).collect()
        } else {
            Vec::new()
        }
    }
}

// ─── JSON Schema Validation ──────────────────────────────────────────────────

/// Simple JSON schema for validation.
#[derive(Debug, Clone)]
pub enum JsonSchema {
    /// Any value is accepted.
    Any,
    /// Must be null.
    Null,
    /// Must be a boolean.
    Bool,
    /// Must be a number.
    Number,
    /// Must be a string.
    StringType,
    /// Must be an array where each element matches the inner schema.
    ArrayOf(Box<JsonSchema>),
    /// Must be an object with the specified required fields.
    ObjectWith {
        /// Required fields: (key, schema).
        required: Vec<(String, JsonSchema)>,
        /// Whether additional fields are allowed.
        allow_extra: bool,
    },
}

impl JsonSchema {
    /// Validate a JSON value against this schema. Returns Ok(()) or an error message.
    pub fn validate(&self, value: &JsonValue) -> Result<(), String> {
        match self {
            JsonSchema::Any => Ok(()),
            JsonSchema::Null => {
                if value.is_null() {
                    Ok(())
                } else {
                    Err("expected null".to_string())
                }
            }
            JsonSchema::Bool => {
                if value.as_bool().is_some() {
                    Ok(())
                } else {
                    Err("expected boolean".to_string())
                }
            }
            JsonSchema::Number => {
                if value.as_f64().is_some() {
                    Ok(())
                } else {
                    Err("expected number".to_string())
                }
            }
            JsonSchema::StringType => {
                if value.as_str().is_some() {
                    Ok(())
                } else {
                    Err("expected string".to_string())
                }
            }
            JsonSchema::ArrayOf(inner) => {
                let arr = value.as_array().ok_or("expected array")?;
                for (i, elem) in arr.iter().enumerate() {
                    inner
                        .validate(elem)
                        .map_err(|e| format!("array[{i}]: {e}"))?;
                }
                Ok(())
            }
            JsonSchema::ObjectWith {
                required,
                allow_extra,
            } => {
                let pairs = value.as_object().ok_or("expected object")?;
                for (key, schema) in required {
                    let val = value
                        .get(key)
                        .ok_or_else(|| format!("missing required key '{key}'"))?;
                    schema
                        .validate(val)
                        .map_err(|e| format!("key '{key}': {e}"))?;
                }
                if !allow_extra {
                    let allowed: std::collections::HashSet<&str> =
                        required.iter().map(|(k, _)| k.as_str()).collect();
                    for (key, _) in pairs {
                        if !allowed.contains(key.as_str()) {
                            return Err(format!("unexpected key '{key}'"));
                        }
                    }
                }
                Ok(())
            }
        }
    }
}

// ─── JSON Diff ───────────────────────────────────────────────────────────────

/// A single difference between two JSON values.
#[derive(Debug, Clone, PartialEq)]
pub struct JsonDiff {
    /// Dot-separated path to the differing value.
    pub path: String,
    /// The value in the left document (None if absent).
    pub left: Option<JsonValue>,
    /// The value in the right document (None if absent).
    pub right: Option<JsonValue>,
}

/// Compute the differences between two JSON values.
pub fn json_diff(left: &JsonValue, right: &JsonValue) -> Vec<JsonDiff> {
    let mut diffs = Vec::new();
    json_diff_inner(left, right, String::new(), &mut diffs);
    diffs
}

fn json_diff_inner(left: &JsonValue, right: &JsonValue, path: String, diffs: &mut Vec<JsonDiff>) {
    if left == right {
        return;
    }
    match (left, right) {
        (JsonValue::Object(lp), JsonValue::Object(rp)) => {
            let mut seen = std::collections::HashSet::new();
            for (k, lv) in lp {
                seen.insert(k.clone());
                let child_path = if path.is_empty() {
                    k.clone()
                } else {
                    format!("{}.{}", path, k)
                };
                if let Some(rv) = rp.iter().find(|(rk, _)| rk == k).map(|(_, v)| v) {
                    json_diff_inner(lv, rv, child_path, diffs);
                } else {
                    diffs.push(JsonDiff {
                        path: child_path,
                        left: Some(lv.clone()),
                        right: None,
                    });
                }
            }
            for (k, rv) in rp {
                if !seen.contains(k) {
                    let child_path = if path.is_empty() {
                        k.clone()
                    } else {
                        format!("{}.{}", path, k)
                    };
                    diffs.push(JsonDiff {
                        path: child_path,
                        left: None,
                        right: Some(rv.clone()),
                    });
                }
            }
        }
        (JsonValue::Array(la), JsonValue::Array(ra)) => {
            let max_len = la.len().max(ra.len());
            for i in 0..max_len {
                let child_path = if path.is_empty() {
                    format!("[{i}]")
                } else {
                    format!("{path}[{i}]")
                };
                match (la.get(i), ra.get(i)) {
                    (Some(lv), Some(rv)) => json_diff_inner(lv, rv, child_path, diffs),
                    (Some(lv), None) => diffs.push(JsonDiff {
                        path: child_path,
                        left: Some(lv.clone()),
                        right: None,
                    }),
                    (None, Some(rv)) => diffs.push(JsonDiff {
                        path: child_path,
                        left: None,
                        right: Some(rv.clone()),
                    }),
                    (None, None) => {}
                }
            }
        }
        _ => {
            diffs.push(JsonDiff {
                path: if path.is_empty() {
                    "$".to_string()
                } else {
                    path
                },
                left: Some(left.clone()),
                right: Some(right.clone()),
            });
        }
    }
}

// ─── JSON Streaming Parser ───────────────────────────────────────────────────

/// A simple streaming JSON parser that processes a JSON array element by element.
///
/// Feeds chunks of text and emits complete JSON values as they are parsed.
pub struct JsonStreamParser {
    buffer: String,
    results: Vec<JsonValue>,
}

impl JsonStreamParser {
    /// Create a new streaming parser.
    pub fn new() -> Self {
        Self {
            buffer: String::new(),
            results: Vec::new(),
        }
    }

    /// Feed a chunk of text to the parser.
    /// Returns the number of complete values parsed from this chunk.
    pub fn feed(&mut self, chunk: &str) -> usize {
        self.buffer.push_str(chunk);
        let mut count = 0;
        loop {
            let trimmed = self.buffer.trim_start();
            if trimmed.is_empty() {
                break;
            }
            // Skip leading commas, brackets for array streaming
            let first = trimmed.as_bytes()[0];
            if first == b'[' || first == b',' {
                self.buffer = trimmed[1..].to_string();
                continue;
            }
            if first == b']' {
                self.buffer = trimmed[1..].to_string();
                continue;
            }
            match parse_value(trimmed) {
                Ok((val, rest)) => {
                    self.results.push(val);
                    self.buffer = rest.to_string();
                    count += 1;
                }
                Err(_) => break, // Incomplete data, wait for more
            }
        }
        count
    }

    /// Take all parsed values, leaving the internal store empty.
    pub fn take_results(&mut self) -> Vec<JsonValue> {
        std::mem::take(&mut self.results)
    }

    /// Number of complete values parsed so far.
    pub fn result_count(&self) -> usize {
        self.results.len()
    }
}
impl Default for JsonStreamParser {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Records helper ──────────────────────────────────────────────────────────

/// Convert a slice of records (each a Vec of key-value pairs) to a JSON array.
pub fn records_to_json(records: &[Vec<(String, JsonValue)>]) -> JsonValue {
    let arr: Vec<JsonValue> = records
        .iter()
        .map(|rec| JsonValue::Object(rec.clone()))
        .collect();
    JsonValue::Array(arr)
}

/// Extract records from a JSON array of objects.
pub fn json_to_records(value: &JsonValue) -> Option<Vec<Vec<(String, JsonValue)>>> {
    let arr = value.as_array()?;
    let mut records = Vec::new();
    for item in arr {
        let pairs = item.as_object()?;
        records.push(pairs.clone());
    }
    Some(records)
}

// ─── Parser helpers ──────────────────────────────────────────────────────────

fn skip_whitespace(s: &str) -> &str {
    s.trim_start_matches(|c: char| c.is_ascii_whitespace())
}

fn parse_value(s: &str) -> Result<(JsonValue, &str), String> {
    let s = skip_whitespace(s);
    if s.is_empty() {
        return Err("Unexpected end of input".to_string());
    }
    match s.as_bytes()[0] {
        b'n' => parse_null(s),
        b't' | b'f' => parse_bool(s),
        b'"' => parse_string(s),
        b'[' => parse_array(s),
        b'{' => parse_object(s),
        b'-' | b'0'..=b'9' => parse_number(s),
        c => Err(format!("Unexpected character: '{}'", c as char)),
    }
}

fn parse_null(s: &str) -> Result<(JsonValue, &str), String> {
    s.strip_prefix("null")
        .map(|rest| (JsonValue::Null, rest))
        .ok_or_else(|| format!("Expected 'null', got: {:?}", &s[..s.len().min(6)]))
}

fn parse_bool(s: &str) -> Result<(JsonValue, &str), String> {
    if let Some(rest) = s.strip_prefix("true") {
        return Ok((JsonValue::Bool(true), rest));
    }
    if let Some(rest) = s.strip_prefix("false") {
        return Ok((JsonValue::Bool(false), rest));
    }
    Err(format!("Expected boolean, got: {:?}", &s[..s.len().min(6)]))
}

fn parse_number(s: &str) -> Result<(JsonValue, &str), String> {
    let end = s
        .find(|c: char| !matches!(c, '-' | '+' | '.' | 'e' | 'E' | '0'..='9'))
        .unwrap_or(s.len());
    let num_str = &s[..end];
    let n: f64 = num_str
        .parse()
        .map_err(|e| format!("Number parse error '{}': {}", num_str, e))?;
    Ok((JsonValue::Number(n), &s[end..]))
}

fn parse_string(s: &str) -> Result<(JsonValue, &str), String> {
    // s must start with '"'
    let s = &s[1..]; // skip opening quote
    let mut result = String::new();
    let mut chars = s.char_indices();
    loop {
        match chars.next() {
            None => return Err("Unterminated string".to_string()),
            Some((_, '"')) => {
                // Find byte position after the closing quote
                let pos = s.len() - chars.as_str().len();
                return Ok((JsonValue::Str(result), &s[pos..]));
            }
            Some((_, '\\')) => match chars.next() {
                Some((_, '"')) => result.push('"'),
                Some((_, '\\')) => result.push('\\'),
                Some((_, '/')) => result.push('/'),
                Some((_, 'n')) => result.push('\n'),
                Some((_, 'r')) => result.push('\r'),
                Some((_, 't')) => result.push('\t'),
                Some((_, 'b')) => result.push('\x08'),
                Some((_, 'f')) => result.push('\x0C'),
                Some((_, 'u')) => {
                    // Basic \uXXXX support
                    let hex: String = (0..4)
                        .map(|_| chars.next().map(|(_, c)| c).unwrap_or('0'))
                        .collect();
                    let code = u32::from_str_radix(&hex, 16)
                        .map_err(|_| format!("Invalid unicode escape \\u{}", hex))?;
                    let ch = char::from_u32(code)
                        .ok_or_else(|| format!("Invalid unicode codepoint U+{:04X}", code))?;
                    result.push(ch);
                }
                Some((_, c)) => return Err(format!("Invalid escape \\{}", c)),
                None => return Err("Unterminated escape".to_string()),
            },
            Some((_, c)) => result.push(c),
        }
    }
}

fn parse_array(s: &str) -> Result<(JsonValue, &str), String> {
    // s starts with '['
    let mut s = skip_whitespace(&s[1..]);
    let mut items = Vec::new();
    if s.starts_with(']') {
        return Ok((JsonValue::Array(items), &s[1..]));
    }
    loop {
        let (val, rest) = parse_value(s)?;
        items.push(val);
        s = skip_whitespace(rest);
        if s.starts_with(']') {
            return Ok((JsonValue::Array(items), &s[1..]));
        }
        if s.starts_with(',') {
            s = skip_whitespace(&s[1..]);
        } else {
            return Err(format!(
                "Expected ',' or ']' in array, got: {:?}",
                &s[..s.len().min(5)]
            ));
        }
    }
}

fn parse_object(s: &str) -> Result<(JsonValue, &str), String> {
    // s starts with '{'
    let mut s = skip_whitespace(&s[1..]);
    let mut pairs = Vec::new();
    if s.starts_with('}') {
        return Ok((JsonValue::Object(pairs), &s[1..]));
    }
    loop {
        s = skip_whitespace(s);
        // Parse key
        if !s.starts_with('"') {
            return Err(format!(
                "Expected string key in object, got: {:?}",
                &s[..s.len().min(5)]
            ));
        }
        let (key_val, rest) = parse_string(s)?;
        let key = match key_val {
            JsonValue::Str(k) => k,
            _ => unreachable!(),
        };
        s = skip_whitespace(rest);
        if !s.starts_with(':') {
            return Err(format!(
                "Expected ':' after key, got: {:?}",
                &s[..s.len().min(5)]
            ));
        }
        s = skip_whitespace(&s[1..]);
        let (val, rest) = parse_value(s)?;
        pairs.push((key, val));
        s = skip_whitespace(rest);
        if s.starts_with('}') {
            return Ok((JsonValue::Object(pairs), &s[1..]));
        }
        if s.starts_with(',') {
            s = skip_whitespace(&s[1..]);
        } else {
            return Err(format!(
                "Expected ',' or '}}' in object, got: {:?}",
                &s[..s.len().min(5)]
            ));
        }
    }
}

// ─── Physics helpers ─────────────────────────────────────────────────────────

/// Serialize a 3D vector `[x, y, z]` to a JSON array.
pub fn vec3_to_json(v: [f64; 3]) -> JsonValue {
    JsonValue::Array(v.iter().copied().map(JsonValue::Number).collect())
}

/// Deserialize a JSON array of 3 numbers into `[f64; 3]`.
pub fn json_to_vec3(j: &JsonValue) -> Option<[f64; 3]> {
    let arr = j.as_array()?;
    if arr.len() < 3 {
        return None;
    }
    Some([arr[0].as_f64()?, arr[1].as_f64()?, arr[2].as_f64()?])
}

/// Serialize a matrix (as nested array) to JSON.
pub fn matrix_to_json(rows: &[Vec<f64>]) -> JsonValue {
    JsonValue::Array(
        rows.iter()
            .map(|row| JsonValue::Array(row.iter().copied().map(JsonValue::Number).collect()))
            .collect(),
    )
}

/// Build a JSON object `{"particles": [...]}` where each element has
/// `"position"` and `"velocity"` arrays.
pub fn particles_to_json(positions: &[[f64; 3]], velocities: &[[f64; 3]]) -> JsonValue {
    let particles: Vec<JsonValue> = positions
        .iter()
        .zip(velocities.iter())
        .map(|(pos, vel)| {
            JsonValue::Object(vec![
                ("position".to_string(), vec3_to_json(*pos)),
                ("velocity".to_string(), vec3_to_json(*vel)),
            ])
        })
        .collect();
    JsonValue::Object(vec![("particles".to_string(), JsonValue::Array(particles))])
}

// ─── JSON Streaming Writer ────────────────────────────────────────────────────

/// A streaming JSON writer that builds a JSON document incrementally.
///
/// Supports arrays and objects.  Call `begin_array`/`end_array` or
/// `begin_object`/`end_object` to open/close containers, then
/// `write_value` to add values.
pub struct JsonStreamingWriter {
    buffer: String,
    /// Stack of (is_array, needs_comma) for each open container.
    stack: Vec<(bool, bool)>,
}

impl JsonStreamingWriter {
    /// Create a new streaming writer.
    pub fn new() -> Self {
        Self {
            buffer: String::new(),
            stack: Vec::new(),
        }
    }

    /// Begin a JSON array `[`.
    pub fn begin_array(&mut self) {
        self.maybe_comma();
        self.buffer.push('[');
        self.stack.push((true, false));
    }

    /// End a JSON array `]`.
    pub fn end_array(&mut self) {
        self.buffer.push(']');
        self.stack.pop();
        if let Some(top) = self.stack.last_mut() {
            top.1 = true; // next sibling needs comma
        }
    }

    /// Begin a JSON object `{`.
    pub fn begin_object(&mut self) {
        self.maybe_comma();
        self.buffer.push('{');
        self.stack.push((false, false));
    }

    /// End a JSON object `}`.
    pub fn end_object(&mut self) {
        self.buffer.push('}');
        self.stack.pop();
        if let Some(top) = self.stack.last_mut() {
            top.1 = true;
        }
    }

    /// Write an object key (must be inside an object).
    pub fn write_key(&mut self, key: &str) {
        self.maybe_comma();
        let escaped = key.replace('\\', "\\\\").replace('"', "\\\"");
        self.buffer.push_str(&format!("\"{}\":", escaped));
        // After writing a key, the next value should NOT be preceded by a comma
        if let Some(top) = self.stack.last_mut() {
            top.1 = false;
        }
    }

    /// Write a JSON value.
    pub fn write_value(&mut self, value: &JsonValue) {
        self.maybe_comma();
        self.buffer.push_str(&value.to_json_string());
        if let Some(top) = self.stack.last_mut() {
            top.1 = true;
        }
    }

    /// Finish writing and return the complete JSON string.
    pub fn finish(self) -> String {
        self.buffer
    }

    fn maybe_comma(&mut self) {
        if let Some(&(_, needs_comma)) = self.stack.last()
            && needs_comma
        {
            self.buffer.push(',');
        }
    }
}

impl Default for JsonStreamingWriter {
    fn default() -> Self {
        Self::new()
    }
}

// ─── JSON Schema Generation ───────────────────────────────────────────────────

impl JsonSchema {
    /// Infer a `JsonSchema` from a sample value.
    ///
    /// Arrays use the schema of their first element; objects use the schemas
    /// of their values.
    pub fn infer_from(value: &JsonValue) -> Self {
        match value {
            JsonValue::Null => JsonSchema::Null,
            JsonValue::Bool(_) => JsonSchema::Bool,
            JsonValue::Number(_) => JsonSchema::Number,
            JsonValue::Str(_) => JsonSchema::StringType,
            JsonValue::Array(arr) => {
                let inner = arr.first().map(Self::infer_from).unwrap_or(JsonSchema::Any);
                JsonSchema::ArrayOf(Box::new(inner))
            }
            JsonValue::Object(pairs) => {
                let required = pairs
                    .iter()
                    .map(|(k, v)| (k.clone(), Self::infer_from(v)))
                    .collect();
                JsonSchema::ObjectWith {
                    required,
                    allow_extra: true,
                }
            }
        }
    }

    /// Serialize the schema to a human-readable string.
    pub fn to_schema_string(&self) -> String {
        match self {
            JsonSchema::Any => "any".to_string(),
            JsonSchema::Null => "null".to_string(),
            JsonSchema::Bool => "boolean".to_string(),
            JsonSchema::Number => "number".to_string(),
            JsonSchema::StringType => "string".to_string(),
            JsonSchema::ArrayOf(inner) => format!("array<{}>", inner.to_schema_string()),
            JsonSchema::ObjectWith {
                required,
                allow_extra,
            } => {
                let fields: Vec<String> = required
                    .iter()
                    .map(|(k, s)| format!("{}: {}", k, s.to_schema_string()))
                    .collect();
                let extra = if *allow_extra { ", ..." } else { "" };
                format!("object{{{}{}}}", fields.join(", "), extra)
            }
        }
    }

    /// Generate a JSON Schema draft-07 compatible object.
    pub fn to_json_schema_object(&self) -> JsonValue {
        match self {
            JsonSchema::Any => JsonValue::Object(vec![]),
            JsonSchema::Null => {
                JsonValue::Object(vec![("type".into(), JsonValue::Str("null".into()))])
            }
            JsonSchema::Bool => {
                JsonValue::Object(vec![("type".into(), JsonValue::Str("boolean".into()))])
            }
            JsonSchema::Number => {
                JsonValue::Object(vec![("type".into(), JsonValue::Str("number".into()))])
            }
            JsonSchema::StringType => {
                JsonValue::Object(vec![("type".into(), JsonValue::Str("string".into()))])
            }
            JsonSchema::ArrayOf(inner) => JsonValue::Object(vec![
                ("type".into(), JsonValue::Str("array".into())),
                ("items".into(), inner.to_json_schema_object()),
            ]),
            JsonSchema::ObjectWith {
                required,
                allow_extra,
            } => {
                let props: Vec<(String, JsonValue)> = required
                    .iter()
                    .map(|(k, s)| (k.clone(), s.to_json_schema_object()))
                    .collect();
                let req_names: Vec<JsonValue> = required
                    .iter()
                    .map(|(k, _)| JsonValue::Str(k.clone()))
                    .collect();
                JsonValue::Object(vec![
                    ("type".into(), JsonValue::Str("object".into())),
                    ("properties".into(), JsonValue::Object(props)),
                    ("required".into(), JsonValue::Array(req_names)),
                    ("additionalProperties".into(), JsonValue::Bool(*allow_extra)),
                ])
            }
        }
    }
}

// ─── JSON Patch / Diff Operations ────────────────────────────────────────────

/// A single JSON patch operation (RFC 6902 inspired).
#[derive(Debug, Clone)]
pub enum JsonPatch {
    /// Add a key-value pair to an object (or replace if key exists).
    Add {
        /// Dot-separated path to the field.
        path: String,
        /// Value to insert.
        value: JsonValue,
    },
    /// Remove a key from an object.
    Remove {
        /// Dot-separated path to the field to remove.
        path: String,
    },
    /// Replace the value at the given path.
    Replace {
        /// Dot-separated path to the field.
        path: String,
        /// New value.
        value: JsonValue,
    },
    /// Copy a value from one path to another.
    Copy {
        /// Source path.
        from: String,
        /// Destination path.
        to: String,
    },
}

impl JsonPatch {
    /// Apply this patch operation to a JSON value in-place.
    pub fn apply(&self, doc: &mut JsonValue) -> Result<(), String> {
        match self {
            JsonPatch::Add { path, value } => patch_set_key(doc, path, value.clone()),
            JsonPatch::Remove { path } => patch_remove_key(doc, path),
            JsonPatch::Replace { path, value } => patch_set_key(doc, path, value.clone()),
            JsonPatch::Copy { from, to } => {
                let val = patch_get_key(doc, from)
                    .ok_or_else(|| format!("path '{}' not found for copy", from))?
                    .clone();
                patch_set_key(doc, to, val)
            }
        }
    }

    /// Serialize this patch to a JSON object.
    pub fn to_json(&self) -> JsonValue {
        match self {
            JsonPatch::Add { path, value } => JsonValue::Object(vec![
                ("op".into(), JsonValue::Str("add".into())),
                ("path".into(), JsonValue::Str(path.clone())),
                ("value".into(), value.clone()),
            ]),
            JsonPatch::Remove { path } => JsonValue::Object(vec![
                ("op".into(), JsonValue::Str("remove".into())),
                ("path".into(), JsonValue::Str(path.clone())),
            ]),
            JsonPatch::Replace { path, value } => JsonValue::Object(vec![
                ("op".into(), JsonValue::Str("replace".into())),
                ("path".into(), JsonValue::Str(path.clone())),
                ("value".into(), value.clone()),
            ]),
            JsonPatch::Copy { from, to } => JsonValue::Object(vec![
                ("op".into(), JsonValue::Str("copy".into())),
                ("from".into(), JsonValue::Str(from.clone())),
                ("path".into(), JsonValue::Str(to.clone())),
            ]),
        }
    }
}

fn patch_set_key(doc: &mut JsonValue, key: &str, value: JsonValue) -> Result<(), String> {
    if let JsonValue::Object(pairs) = doc {
        if let Some(pos) = pairs.iter().position(|(k, _)| k == key) {
            pairs[pos].1 = value;
        } else {
            pairs.push((key.to_string(), value));
        }
        Ok(())
    } else {
        Err(format!("cannot set key '{}' on non-object", key))
    }
}

fn patch_remove_key(doc: &mut JsonValue, key: &str) -> Result<(), String> {
    if let JsonValue::Object(pairs) = doc {
        let before = pairs.len();
        pairs.retain(|(k, _)| k != key);
        if pairs.len() < before {
            Ok(())
        } else {
            Err(format!("key '{}' not found for remove", key))
        }
    } else {
        Err(format!("cannot remove key '{}' from non-object", key))
    }
}

fn patch_get_key<'a>(doc: &'a JsonValue, key: &str) -> Option<&'a JsonValue> {
    doc.get(key)
}

/// Apply a sequence of patch operations to a document.
pub fn apply_patch_sequence(doc: &mut JsonValue, patches: &[JsonPatch]) -> Result<(), String> {
    for patch in patches {
        patch.apply(doc)?;
    }
    Ok(())
}

/// Generate a sequence of `JsonPatch` operations that transforms `left` into `right`.
pub fn generate_patch(left: &JsonValue, right: &JsonValue) -> Vec<JsonPatch> {
    let diffs = json_diff(left, right);
    diffs
        .into_iter()
        .map(|diff| match (diff.left, diff.right) {
            (None, Some(v)) => JsonPatch::Add {
                path: diff.path,
                value: v,
            },
            (Some(_), None) => JsonPatch::Remove { path: diff.path },
            (Some(_), Some(v)) => JsonPatch::Replace {
                path: diff.path,
                value: v,
            },
            (None, None) => JsonPatch::Remove { path: diff.path },
        })
        .collect()
}

// ─── JSONPath Query ───────────────────────────────────────────────────────────

/// Evaluate a basic JSONPath expression against a JSON value.
///
/// Supported syntax:
/// - `$` — root document
/// - `$.a.b` — nested field access
/// - `$[n]` — array index
/// - `$.*` or `$[*]` — all children
/// - Combinations: `$.a[0].b`
pub fn jsonpath_query<'a>(root: &'a JsonValue, path: &str) -> Vec<&'a JsonValue> {
    let path = path.trim();
    if path == "$" {
        return vec![root];
    }
    if !path.starts_with('$') {
        return vec![];
    }

    // Parse path segments
    let rest = &path[1..]; // drop '$'
    let segments = parse_jsonpath_segments(rest);
    let mut current: Vec<&'a JsonValue> = vec![root];

    for seg in segments {
        let mut next = Vec::new();
        for node in current {
            match seg.as_str() {
                "*" => match node {
                    JsonValue::Object(pairs) => {
                        for (_, v) in pairs {
                            next.push(v);
                        }
                    }
                    JsonValue::Array(arr) => {
                        for v in arr {
                            next.push(v);
                        }
                    }
                    _ => {}
                },
                s if s.starts_with('[') && s.ends_with(']') => {
                    let inner = &s[1..s.len() - 1];
                    if inner == "*" {
                        if let JsonValue::Array(arr) = node {
                            for v in arr {
                                next.push(v);
                            }
                        }
                    } else if let Ok(idx) = inner.parse::<usize>()
                        && let JsonValue::Array(arr) = node
                        && let Some(v) = arr.get(idx)
                    {
                        next.push(v);
                    }
                }
                key => {
                    if let Some(v) = node.get(key) {
                        next.push(v);
                    }
                }
            }
        }
        current = next;
    }
    current
}

fn parse_jsonpath_segments(s: &str) -> Vec<String> {
    let mut segments = Vec::new();
    let mut remaining = s;

    while !remaining.is_empty() {
        if remaining.starts_with('.') {
            remaining = &remaining[1..];
            // Read until next '.' or '['
            let end = remaining.find(['.', '[']).unwrap_or(remaining.len());
            if end > 0 {
                segments.push(remaining[..end].to_string());
                remaining = &remaining[end..];
            }
        } else if remaining.starts_with('[') {
            let end = remaining.find(']').unwrap_or(remaining.len() - 1);
            segments.push(remaining[..=end].to_string());
            remaining = &remaining[(end + 1)..];
        } else {
            // Try reading as bare word
            let end = remaining.find(['.', '[']).unwrap_or(remaining.len());
            segments.push(remaining[..end].to_string());
            remaining = &remaining[end..];
        }
    }
    segments
}

// ─── JSON Compression ─────────────────────────────────────────────────────────

/// Compress a JSON string by removing insignificant whitespace.
///
/// Preserves whitespace inside string values.
pub fn json_compress(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    let mut in_string = false;
    let mut chars = input.chars().peekable();

    while let Some(c) = chars.next() {
        if in_string {
            out.push(c);
            if c == '\\' {
                // Escaped character: consume the next char verbatim
                if let Some(next) = chars.next() {
                    out.push(next);
                }
            } else if c == '"' {
                in_string = false;
            }
        } else {
            match c {
                '"' => {
                    in_string = true;
                    out.push(c);
                }
                ' ' | '\t' | '\n' | '\r' => {
                    // Skip whitespace outside strings
                }
                _ => {
                    out.push(c);
                }
            }
        }
    }
    out
}

/// Statistics about JSON compression.
#[derive(Debug, Clone)]
pub struct JsonCompressionStats {
    /// Original string length.
    pub original_len: usize,
    /// Compressed string length.
    pub compressed_len: usize,
}

impl JsonCompressionStats {
    /// Compute compression statistics for a JSON string.
    pub fn compute(input: &str) -> Self {
        let compressed = json_compress(input);
        Self {
            original_len: input.len(),
            compressed_len: compressed.len(),
        }
    }

    /// Compression ratio (compressed / original). Lower is better.
    pub fn compression_ratio(&self) -> f64 {
        if self.original_len == 0 {
            return 1.0;
        }
        self.compressed_len as f64 / self.original_len as f64
    }

    /// Bytes saved by compression.
    pub fn bytes_saved(&self) -> usize {
        self.original_len.saturating_sub(self.compressed_len)
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_null_roundtrip() {
        let j = JsonValue::Null;
        let s = j.to_json_string();
        assert_eq!(s, "null");
        let parsed = JsonValue::from_str(&s).unwrap();
        assert_eq!(parsed, JsonValue::Null);
    }

    #[test]
    fn test_bool_true() {
        let j = JsonValue::Bool(true);
        let s = j.to_json_string();
        assert_eq!(s, "true");
        assert_eq!(JsonValue::from_str(&s).unwrap(), JsonValue::Bool(true));
    }

    #[test]
    fn test_bool_false() {
        let j = JsonValue::Bool(false);
        assert_eq!(j.to_json_string(), "false");
        assert_eq!(
            JsonValue::from_str("false").unwrap(),
            JsonValue::Bool(false)
        );
    }

    #[test]
    fn test_number_integer() {
        let j = JsonValue::Number(42.0);
        let s = j.to_json_string();
        assert_eq!(s, "42");
        let parsed = JsonValue::from_str(&s).unwrap();
        assert_eq!(parsed.as_f64(), Some(42.0));
    }

    #[test]
    fn test_number_float() {
        let j = JsonValue::Number(3.125);
        let s = j.to_json_string();
        let parsed = JsonValue::from_str(&s).unwrap();
        let v = parsed.as_f64().unwrap();
        assert!((v - 3.125).abs() < 1e-10);
    }

    #[test]
    fn test_string_roundtrip() {
        let j = JsonValue::Str("hello world".to_string());
        let s = j.to_json_string();
        assert_eq!(s, "\"hello world\"");
        let parsed = JsonValue::from_str(&s).unwrap();
        assert_eq!(parsed.as_str(), Some("hello world"));
    }

    #[test]
    fn test_string_escape() {
        let j = JsonValue::Str("line1\nline2".to_string());
        let s = j.to_json_string();
        assert!(s.contains("\\n"));
        let parsed = JsonValue::from_str(&s).unwrap();
        assert_eq!(parsed.as_str(), Some("line1\nline2"));
    }

    #[test]
    fn test_array_roundtrip() {
        let j = JsonValue::Array(vec![
            JsonValue::Number(1.0),
            JsonValue::Number(2.0),
            JsonValue::Number(3.0),
        ]);
        let s = j.to_json_string();
        let parsed = JsonValue::from_str(&s).unwrap();
        let arr = parsed.as_array().unwrap();
        assert_eq!(arr.len(), 3);
        assert_eq!(arr[0].as_f64(), Some(1.0));
    }

    #[test]
    fn test_empty_array() {
        let j = JsonValue::Array(vec![]);
        assert_eq!(j.to_json_string(), "[]");
        let parsed = JsonValue::from_str("[]").unwrap();
        assert_eq!(parsed.as_array().unwrap().len(), 0);
    }

    #[test]
    fn test_object_get() {
        let j = JsonValue::Object(vec![
            ("x".to_string(), JsonValue::Number(1.0)),
            ("y".to_string(), JsonValue::Number(2.0)),
        ]);
        assert_eq!(j.get("x").unwrap().as_f64(), Some(1.0));
        assert!(j.get("z").is_none());
    }

    #[test]
    fn test_object_roundtrip() {
        let j = JsonValue::Object(vec![
            ("name".to_string(), JsonValue::Str("particle".to_string())),
            ("mass".to_string(), JsonValue::Number(1.5)),
        ]);
        let s = j.to_json_string();
        let parsed = JsonValue::from_str(&s).unwrap();
        assert_eq!(parsed.get("name").unwrap().as_str(), Some("particle"));
        assert_eq!(parsed.get("mass").unwrap().as_f64(), Some(1.5));
    }

    #[test]
    fn test_nested_object() {
        let s = r#"{"pos":{"x":1,"y":2}}"#;
        let parsed = JsonValue::from_str(s).unwrap();
        let pos = parsed.get("pos").unwrap();
        assert_eq!(pos.get("x").unwrap().as_f64(), Some(1.0));
    }

    #[test]
    fn test_get_index() {
        let j = JsonValue::Array(vec![JsonValue::Number(10.0), JsonValue::Number(20.0)]);
        assert_eq!(j.get_index(0).unwrap().as_f64(), Some(10.0));
        assert!(j.get_index(5).is_none());
    }

    #[test]
    fn test_vec3_roundtrip() {
        let v = [1.0, 2.0, 3.0f64];
        let j = vec3_to_json(v);
        let v2 = json_to_vec3(&j).unwrap();
        assert_eq!(v, v2);
    }

    #[test]
    fn test_json_to_vec3_too_short() {
        let j = JsonValue::Array(vec![JsonValue::Number(1.0)]);
        assert!(json_to_vec3(&j).is_none());
    }

    #[test]
    fn test_particles_to_json() {
        let pos = [[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]];
        let vel = [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]];
        let j = particles_to_json(&pos, &vel);
        let arr = j.get("particles").unwrap().as_array().unwrap();
        assert_eq!(arr.len(), 2);
        let p0 = &arr[0];
        let p0_pos = json_to_vec3(p0.get("position").unwrap()).unwrap();
        assert_eq!(p0_pos, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_parse_negative_number() {
        let j = JsonValue::from_str("-3.125").unwrap();
        assert!((j.as_f64().unwrap() + 3.125).abs() < 1e-10);
    }

    #[test]
    fn test_whitespace_tolerance() {
        let s = "  {  \"k\"  :  42  }  ";
        let j = JsonValue::from_str(s).unwrap();
        assert_eq!(j.get("k").unwrap().as_f64(), Some(42.0));
    }

    #[test]
    fn test_empty_object() {
        let j = JsonValue::from_str("{}").unwrap();
        assert!(j.get("anything").is_none());
    }

    #[test]
    fn test_invalid_input_error() {
        assert!(JsonValue::from_str("not_json").is_err());
    }

    // ── New tests for expanded features ─────────────────────────────────────

    #[test]
    fn test_get_path() {
        let s = r#"{"a":{"b":{"c":42}}}"#;
        let j = JsonValue::from_str(s).unwrap();
        assert_eq!(j.get_path("a.b.c").unwrap().as_f64(), Some(42.0));
        assert!(j.get_path("a.b.d").is_none());
        assert!(j.get_path("x").is_none());
    }

    #[test]
    fn test_json_len_and_is_empty() {
        let arr = JsonValue::Array(vec![JsonValue::Number(1.0), JsonValue::Number(2.0)]);
        assert_eq!(arr.len(), 2);
        assert!(!arr.is_empty());

        let empty_arr = JsonValue::Array(vec![]);
        assert_eq!(empty_arr.len(), 0);
        assert!(empty_arr.is_empty());

        assert!(JsonValue::Null.is_empty());
        assert!(!JsonValue::Number(1.0).is_empty());
    }

    #[test]
    fn test_json_merge() {
        let base = JsonValue::from_str(r#"{"a":1,"b":2}"#).unwrap();
        let overlay = JsonValue::from_str(r#"{"b":3,"c":4}"#).unwrap();
        let merged = base.merge(&overlay);
        assert_eq!(merged.get("a").unwrap().as_f64(), Some(1.0));
        assert_eq!(merged.get("b").unwrap().as_f64(), Some(3.0));
        assert_eq!(merged.get("c").unwrap().as_f64(), Some(4.0));
    }

    #[test]
    fn test_json_merge_nested() {
        let base = JsonValue::from_str(r#"{"a":{"x":1,"y":2}}"#).unwrap();
        let overlay = JsonValue::from_str(r#"{"a":{"y":3,"z":4}}"#).unwrap();
        let merged = base.merge(&overlay);
        let a = merged.get("a").unwrap();
        assert_eq!(a.get("x").unwrap().as_f64(), Some(1.0));
        assert_eq!(a.get("y").unwrap().as_f64(), Some(3.0));
        assert_eq!(a.get("z").unwrap().as_f64(), Some(4.0));
    }

    #[test]
    fn test_json_keys_values() {
        let j = JsonValue::from_str(r#"{"a":1,"b":2}"#).unwrap();
        let keys = j.keys();
        assert_eq!(keys.len(), 2);
        assert!(keys.contains(&"a"));
        assert!(keys.contains(&"b"));
        let vals = j.values();
        assert_eq!(vals.len(), 2);
    }

    #[test]
    fn test_schema_validation_number() {
        let schema = JsonSchema::Number;
        assert!(schema.validate(&JsonValue::Number(42.0)).is_ok());
        assert!(schema.validate(&JsonValue::Str("nope".into())).is_err());
    }

    #[test]
    fn test_schema_validation_object() {
        let schema = JsonSchema::ObjectWith {
            required: vec![
                ("name".into(), JsonSchema::StringType),
                ("mass".into(), JsonSchema::Number),
            ],
            allow_extra: true,
        };
        let valid = JsonValue::from_str(r#"{"name":"ball","mass":1.5,"extra":true}"#).unwrap();
        assert!(schema.validate(&valid).is_ok());

        let missing_mass = JsonValue::from_str(r#"{"name":"ball"}"#).unwrap();
        assert!(schema.validate(&missing_mass).is_err());
    }

    #[test]
    fn test_schema_no_extra_keys() {
        let schema = JsonSchema::ObjectWith {
            required: vec![("x".into(), JsonSchema::Number)],
            allow_extra: false,
        };
        let valid = JsonValue::from_str(r#"{"x":1}"#).unwrap();
        assert!(schema.validate(&valid).is_ok());

        let extra = JsonValue::from_str(r#"{"x":1,"y":2}"#).unwrap();
        assert!(schema.validate(&extra).is_err());
    }

    #[test]
    fn test_schema_array_of() {
        let schema = JsonSchema::ArrayOf(Box::new(JsonSchema::Number));
        let valid = JsonValue::from_str("[1,2,3]").unwrap();
        assert!(schema.validate(&valid).is_ok());

        let invalid = JsonValue::from_str(r#"[1,"two",3]"#).unwrap();
        assert!(schema.validate(&invalid).is_err());
    }

    #[test]
    fn test_json_diff_same() {
        let a = JsonValue::from_str(r#"{"x":1}"#).unwrap();
        let b = JsonValue::from_str(r#"{"x":1}"#).unwrap();
        let diffs = json_diff(&a, &b);
        assert!(diffs.is_empty());
    }

    #[test]
    fn test_json_diff_different_value() {
        let a = JsonValue::from_str(r#"{"x":1}"#).unwrap();
        let b = JsonValue::from_str(r#"{"x":2}"#).unwrap();
        let diffs = json_diff(&a, &b);
        assert_eq!(diffs.len(), 1);
        assert_eq!(diffs[0].path, "x");
    }

    #[test]
    fn test_json_diff_added_key() {
        let a = JsonValue::from_str(r#"{"x":1}"#).unwrap();
        let b = JsonValue::from_str(r#"{"x":1,"y":2}"#).unwrap();
        let diffs = json_diff(&a, &b);
        assert_eq!(diffs.len(), 1);
        assert_eq!(diffs[0].path, "y");
        assert!(diffs[0].left.is_none());
        assert!(diffs[0].right.is_some());
    }

    #[test]
    fn test_json_diff_removed_key() {
        let a = JsonValue::from_str(r#"{"x":1,"y":2}"#).unwrap();
        let b = JsonValue::from_str(r#"{"x":1}"#).unwrap();
        let diffs = json_diff(&a, &b);
        assert_eq!(diffs.len(), 1);
        assert_eq!(diffs[0].path, "y");
        assert!(diffs[0].left.is_some());
        assert!(diffs[0].right.is_none());
    }

    #[test]
    fn test_json_diff_array() {
        let a = JsonValue::from_str("[1,2,3]").unwrap();
        let b = JsonValue::from_str("[1,2,4]").unwrap();
        let diffs = json_diff(&a, &b);
        assert_eq!(diffs.len(), 1);
        assert!(diffs[0].path.contains("[2]"));
    }

    #[test]
    fn test_streaming_parser() {
        let mut parser = JsonStreamParser::new();
        let count = parser.feed(r#"[{"a":1},{"b":2}]"#);
        assert_eq!(count, 2);
        let results = parser.take_results();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].get("a").unwrap().as_f64(), Some(1.0));
        assert_eq!(results[1].get("b").unwrap().as_f64(), Some(2.0));
    }

    #[test]
    fn test_streaming_parser_incremental() {
        let mut parser = JsonStreamParser::new();
        parser.feed("[42,");
        assert_eq!(parser.result_count(), 1);
        parser.feed("43]");
        assert_eq!(parser.result_count(), 2);
        let results = parser.take_results();
        assert_eq!(results[0].as_f64(), Some(42.0));
        assert_eq!(results[1].as_f64(), Some(43.0));
    }

    #[test]
    fn test_records_to_json() {
        let records = vec![
            vec![
                ("id".into(), JsonValue::Number(1.0)),
                ("name".into(), JsonValue::Str("a".into())),
            ],
            vec![
                ("id".into(), JsonValue::Number(2.0)),
                ("name".into(), JsonValue::Str("b".into())),
            ],
        ];
        let j = records_to_json(&records);
        let arr = j.as_array().unwrap();
        assert_eq!(arr.len(), 2);
        assert_eq!(arr[0].get("id").unwrap().as_f64(), Some(1.0));
    }

    #[test]
    fn test_json_to_records() {
        let s = r#"[{"id":1},{"id":2}]"#;
        let j = JsonValue::from_str(s).unwrap();
        let records = json_to_records(&j).unwrap();
        assert_eq!(records.len(), 2);
    }

    #[test]
    fn test_pretty_print() {
        let j = JsonValue::from_str(r#"{"a":1,"b":[2,3]}"#).unwrap();
        let pretty = j.to_json_pretty(2);
        assert!(pretty.contains('\n'));
        assert!(pretty.contains("\"a\""));
        // Should be parseable back
        let reparsed = JsonValue::from_str(&pretty).unwrap();
        assert_eq!(reparsed.get("a").unwrap().as_f64(), Some(1.0));
    }

    #[test]
    fn test_as_bool() {
        assert_eq!(JsonValue::Bool(true).as_bool(), Some(true));
        assert_eq!(JsonValue::Number(1.0).as_bool(), None);
    }

    #[test]
    fn test_as_object() {
        let j = JsonValue::from_str(r#"{"a":1}"#).unwrap();
        assert!(j.as_object().is_some());
        assert!(JsonValue::Number(1.0).as_object().is_none());
    }

    #[test]
    fn test_matrix_to_json() {
        let m = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let j = matrix_to_json(&m);
        let arr = j.as_array().unwrap();
        assert_eq!(arr.len(), 2);
        assert_eq!(arr[0].get_index(0).unwrap().as_f64(), Some(1.0));
        assert_eq!(arr[1].get_index(1).unwrap().as_f64(), Some(4.0));
    }

    // ── JSON streaming writer tests ───────────────────────────────────────

    #[test]
    fn test_streaming_writer_empty_array() {
        let mut w = JsonStreamingWriter::new();
        w.begin_array();
        w.end_array();
        let s = w.finish();
        assert_eq!(s, "[]");
    }

    #[test]
    fn test_streaming_writer_single_value() {
        let mut w = JsonStreamingWriter::new();
        w.begin_array();
        w.write_value(&JsonValue::Number(42.0));
        w.end_array();
        let s = w.finish();
        assert_eq!(
            JsonValue::from_str(&s).unwrap(),
            JsonValue::Array(vec![JsonValue::Number(42.0)])
        );
    }

    #[test]
    fn test_streaming_writer_multiple_values() {
        let mut w = JsonStreamingWriter::new();
        w.begin_array();
        w.write_value(&JsonValue::Number(1.0));
        w.write_value(&JsonValue::Number(2.0));
        w.write_value(&JsonValue::Str("three".to_string()));
        w.end_array();
        let s = w.finish();
        let parsed = JsonValue::from_str(&s).unwrap();
        let arr = parsed.as_array().unwrap();
        assert_eq!(arr.len(), 3);
    }

    #[test]
    fn test_streaming_writer_nested_object() {
        let mut w = JsonStreamingWriter::new();
        w.begin_object();
        w.write_key("x");
        w.write_value(&JsonValue::Number(1.0));
        w.write_key("y");
        w.write_value(&JsonValue::Number(2.0));
        w.end_object();
        let s = w.finish();
        let parsed = JsonValue::from_str(&s).unwrap();
        assert_eq!(parsed.get("x").unwrap().as_f64(), Some(1.0));
    }

    // ── JSON schema generation tests ──────────────────────────────────────

    #[test]
    fn test_schema_from_number() {
        let v = JsonValue::Number(42.0);
        let schema = JsonSchema::infer_from(&v);
        assert!(schema.validate(&JsonValue::Number(1.0)).is_ok());
        assert!(schema.validate(&JsonValue::Str("nope".into())).is_err());
    }

    #[test]
    fn test_schema_from_object() {
        let j = JsonValue::from_str(r#"{"name":"alice","age":30}"#).unwrap();
        let schema = JsonSchema::infer_from(&j);
        let valid = JsonValue::from_str(r#"{"name":"bob","age":25}"#).unwrap();
        assert!(schema.validate(&valid).is_ok());
    }

    #[test]
    fn test_schema_from_array() {
        let j = JsonValue::from_str("[1,2,3]").unwrap();
        let schema = JsonSchema::infer_from(&j);
        let valid = JsonValue::from_str("[4,5,6]").unwrap();
        assert!(schema.validate(&valid).is_ok());
        let invalid = JsonValue::from_str(r#"["a","b"]"#).unwrap();
        assert!(schema.validate(&invalid).is_err());
    }

    #[test]
    fn test_schema_to_string() {
        let schema = JsonSchema::Number;
        let s = schema.to_schema_string();
        assert!(s.contains("number"));
    }

    // ── JSON patch / diff operations tests ────────────────────────────────

    #[test]
    fn test_json_patch_add() {
        let mut doc = JsonValue::from_str(r#"{"a":1}"#).unwrap();
        let patch = JsonPatch::Add {
            path: "b".to_string(),
            value: JsonValue::Number(2.0),
        };
        patch.apply(&mut doc).unwrap();
        assert_eq!(doc.get("b").unwrap().as_f64(), Some(2.0));
    }

    #[test]
    fn test_json_patch_remove() {
        let mut doc = JsonValue::from_str(r#"{"a":1,"b":2}"#).unwrap();
        let patch = JsonPatch::Remove {
            path: "b".to_string(),
        };
        patch.apply(&mut doc).unwrap();
        assert!(doc.get("b").is_none());
        assert_eq!(doc.get("a").unwrap().as_f64(), Some(1.0));
    }

    #[test]
    fn test_json_patch_replace() {
        let mut doc = JsonValue::from_str(r#"{"a":1}"#).unwrap();
        let patch = JsonPatch::Replace {
            path: "a".to_string(),
            value: JsonValue::Number(99.0),
        };
        patch.apply(&mut doc).unwrap();
        assert_eq!(doc.get("a").unwrap().as_f64(), Some(99.0));
    }

    #[test]
    fn test_json_patch_sequence() {
        let mut doc = JsonValue::from_str(r#"{"x":1}"#).unwrap();
        let patches = vec![
            JsonPatch::Add {
                path: "y".to_string(),
                value: JsonValue::Number(2.0),
            },
            JsonPatch::Replace {
                path: "x".to_string(),
                value: JsonValue::Str("hello".into()),
            },
        ];
        apply_patch_sequence(&mut doc, &patches).unwrap();
        assert_eq!(doc.get("y").unwrap().as_f64(), Some(2.0));
        assert_eq!(doc.get("x").unwrap().as_str(), Some("hello"));
    }

    // ── JSONPath basics tests ──────────────────────────────────────────────

    #[test]
    fn test_jsonpath_root() {
        let j = JsonValue::from_str(r#"{"a":1}"#).unwrap();
        let results = jsonpath_query(&j, "$");
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].get("a").unwrap().as_f64(), Some(1.0));
    }

    #[test]
    fn test_jsonpath_child() {
        let j = JsonValue::from_str(r#"{"a":{"b":42}}"#).unwrap();
        let results = jsonpath_query(&j, "$.a.b");
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].as_f64(), Some(42.0));
    }

    #[test]
    fn test_jsonpath_array_index() {
        let j = JsonValue::from_str("[10,20,30]").unwrap();
        let results = jsonpath_query(&j, "$[1]");
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].as_f64(), Some(20.0));
    }

    #[test]
    fn test_jsonpath_wildcard() {
        let j = JsonValue::from_str(r#"{"a":1,"b":2,"c":3}"#).unwrap();
        let results = jsonpath_query(&j, "$.*");
        assert_eq!(results.len(), 3);
    }

    #[test]
    fn test_jsonpath_array_wildcard() {
        let j = JsonValue::from_str("[1,2,3]").unwrap();
        let results = jsonpath_query(&j, "$[*]");
        assert_eq!(results.len(), 3);
    }

    // ── JSON compression tests ─────────────────────────────────────────────

    #[test]
    fn test_json_compress_removes_whitespace() {
        let pretty = r#"{
  "a": 1,
  "b": [
    2,
    3
  ]
}"#;
        let compressed = json_compress(pretty);
        assert!(
            !compressed.contains('\n'),
            "compressed should not contain newlines"
        );
        assert!(
            !compressed.contains("  "),
            "compressed should not contain double spaces"
        );
        let parsed = JsonValue::from_str(&compressed).unwrap();
        assert_eq!(parsed.get("a").unwrap().as_f64(), Some(1.0));
    }

    #[test]
    fn test_json_compress_roundtrip() {
        let j = JsonValue::from_str(r#"{"x":1,"y":[2,3],"z":true}"#).unwrap();
        let pretty = j.to_json_pretty(4);
        let compressed = json_compress(&pretty);
        let reparsed = JsonValue::from_str(&compressed).unwrap();
        assert_eq!(reparsed, j);
    }

    #[test]
    fn test_json_compress_stats() {
        let pretty = r#"{
  "alpha": 1,
  "beta":  2
}"#;
        let stats = JsonCompressionStats::compute(pretty);
        assert!(
            stats.original_len > stats.compressed_len,
            "compression should reduce size: {} > {}",
            stats.original_len,
            stats.compressed_len
        );
        assert!(stats.compression_ratio() < 1.0);
    }

    #[test]
    fn test_json_compress_preserves_strings_with_spaces() {
        let j = JsonValue::from_str(r#"{"msg":"hello world"}"#).unwrap();
        let pretty = j.to_json_pretty(2);
        let compressed = json_compress(&pretty);
        let reparsed = JsonValue::from_str(&compressed).unwrap();
        assert_eq!(reparsed.get("msg").unwrap().as_str(), Some("hello world"));
    }
}
