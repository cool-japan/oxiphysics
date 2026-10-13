// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Line-oriented ASCII restart encoding.
//!
//! ```text
//! # OxiPhysics restart file
//! VERSION <text>
//! TIMESTAMP <u64>
//! STEP <u64>
//! TIME <f64>
//! CRATE <text>
//! DESC <text>
//! N_PARTICLES <n>
//! BEGIN_POSITIONS / END_POSITIONS      n lines of `x y z`
//! BEGIN_VELOCITIES / END_VELOCITIES    0 or n lines of `x y z`
//! BEGIN_FORCES / END_FORCES            0 or n lines of `x y z`
//! BEGIN_MASSES / END_MASSES            0 or n lines, one value each
//! BEGIN_TYPES / END_TYPES              0 or n lines, one `u32` each
//! BEGIN_BOX / END_BOX                  3 lines of `x y z`
//! BEGIN_EXTRA_SCALAR <name> / END_EXTRA_SCALAR    repeated, one value per line
//! BEGIN_EXTRA_VECTOR <name> / END_EXTRA_VECTOR    repeated, `x y z` per line
//! ```
//!
//! Floating-point values are written with the shortest decimal form that
//! parses back to the identical `f64` (Rust's `{:?}`), so the format is
//! lossless. The free-text values (`VERSION`, `CRATE`, `DESC` and the extra
//! array names) escape `\`, line feed and carriage return as `\\`, `\n` and
//! `\r`; everything after the single separator space is kept verbatim.
//!
//! The reader is strict: malformed numbers, wrong column counts, unknown
//! lines, unterminated or duplicated sections and arrays whose length
//! disagrees with `N_PARTICLES` are reported as errors that carry the line
//! number.

use super::{RestartData, RestartMetadata};

// ── Encoder ───────────────────────────────────────────────────────────────────

/// Escape `\`, LF and CR so a free-text value stays on one line.
fn escape_text(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            other => out.push(other),
        }
    }
    out
}

/// Inverse of [`escape_text`]. Unrecognised escape sequences are kept
/// verbatim so older files containing a bare backslash still read back.
fn unescape_text(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(ch) = chars.next() {
        if ch != '\\' {
            out.push(ch);
            continue;
        }
        match chars.next() {
            Some('\\') => out.push('\\'),
            Some('n') => out.push('\n'),
            Some('r') => out.push('\r'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

fn push_vec3_lines(s: &mut String, tag: &str, vectors: &[[f64; 3]]) {
    s.push_str(&format!("BEGIN_{tag}\n"));
    for v in vectors {
        s.push_str(&format!("{:?} {:?} {:?}\n", v[0], v[1], v[2]));
    }
    s.push_str(&format!("END_{tag}\n"));
}

/// Encode `data` as ASCII text. Fails if [`RestartData::validate`] fails.
pub(super) fn encode(data: &RestartData) -> Result<String, String> {
    data.validate()?;
    let mut s = String::new();
    s.push_str("# OxiPhysics restart file\n");
    s.push_str(&format!(
        "VERSION {}\n",
        escape_text(&data.metadata.version)
    ));
    s.push_str(&format!("TIMESTAMP {}\n", data.metadata.timestamp));
    s.push_str(&format!("STEP {}\n", data.metadata.step));
    s.push_str(&format!("TIME {:?}\n", data.metadata.time));
    s.push_str(&format!(
        "CRATE {}\n",
        escape_text(&data.metadata.crate_name)
    ));
    s.push_str(&format!(
        "DESC {}\n",
        escape_text(&data.metadata.description)
    ));
    s.push_str(&format!("N_PARTICLES {}\n", data.n_particles()));
    push_vec3_lines(&mut s, "POSITIONS", &data.positions);
    push_vec3_lines(&mut s, "VELOCITIES", &data.velocities);
    push_vec3_lines(&mut s, "FORCES", &data.forces);
    s.push_str("BEGIN_MASSES\n");
    for m in &data.masses {
        s.push_str(&format!("{m:?}\n"));
    }
    s.push_str("END_MASSES\n");
    s.push_str("BEGIN_TYPES\n");
    for t in &data.types {
        s.push_str(&format!("{t}\n"));
    }
    s.push_str("END_TYPES\n");
    push_vec3_lines(&mut s, "BOX", &data.box_matrix);
    for (name, values) in &data.extra_scalars {
        s.push_str(&format!("BEGIN_EXTRA_SCALAR {}\n", escape_text(name)));
        for v in values {
            s.push_str(&format!("{v:?}\n"));
        }
        s.push_str("END_EXTRA_SCALAR\n");
    }
    for (name, vectors) in &data.extra_vectors {
        s.push_str(&format!("BEGIN_EXTRA_VECTOR {}\n", escape_text(name)));
        for v in vectors {
            s.push_str(&format!("{:?} {:?} {:?}\n", v[0], v[1], v[2]));
        }
        s.push_str("END_EXTRA_VECTOR\n");
    }
    Ok(s)
}

// ── Decoder ───────────────────────────────────────────────────────────────────

/// The fixed sections of the format.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Core {
    Positions,
    Velocities,
    Forces,
    Masses,
    Types,
    Box,
}

impl Core {
    /// Position in [`Parsed::seen`].
    fn index(self) -> usize {
        match self {
            Core::Positions => 0,
            Core::Velocities => 1,
            Core::Forces => 2,
            Core::Masses => 3,
            Core::Types => 4,
            Core::Box => 5,
        }
    }
}

/// `(section, begin marker, end marker)` for the fixed sections.
const CORE_SECTIONS: [(Core, &str, &str); 6] = [
    (Core::Positions, "BEGIN_POSITIONS", "END_POSITIONS"),
    (Core::Velocities, "BEGIN_VELOCITIES", "END_VELOCITIES"),
    (Core::Forces, "BEGIN_FORCES", "END_FORCES"),
    (Core::Masses, "BEGIN_MASSES", "END_MASSES"),
    (Core::Types, "BEGIN_TYPES", "END_TYPES"),
    (Core::Box, "BEGIN_BOX", "END_BOX"),
];

const BEGIN_EXTRA_SCALAR: &str = "BEGIN_EXTRA_SCALAR";
const END_EXTRA_SCALAR: &str = "END_EXTRA_SCALAR";
const BEGIN_EXTRA_VECTOR: &str = "BEGIN_EXTRA_VECTOR";
const END_EXTRA_VECTOR: &str = "END_EXTRA_VECTOR";

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Kind {
    Core(Core),
    ExtraScalar,
    ExtraVector,
}

/// A section that has been opened and not yet closed.
struct Section {
    kind: Kind,
    /// Extra-array name (empty for the fixed sections).
    name: String,
    /// Marker that closes this section.
    end_marker: &'static str,
    /// Line number of the `BEGIN_*` marker.
    begin_line: usize,
    /// Values of `x y z` sections (positions, velocities, forces, box, extra vectors).
    vectors: Vec<[f64; 3]>,
    /// Values of scalar sections (masses, extra scalars).
    scalars: Vec<f64>,
    /// Values of the `types` section.
    ints: Vec<u32>,
}

/// Everything collected from closed sections.
#[derive(Default)]
struct Parsed {
    positions: Vec<[f64; 3]>,
    velocities: Vec<[f64; 3]>,
    forces: Vec<[f64; 3]>,
    masses: Vec<f64>,
    types: Vec<u32>,
    box_matrix: Option<[[f64; 3]; 3]>,
    seen: [bool; 6],
    extra_scalars: Vec<(String, Vec<f64>)>,
    extra_vectors: Vec<(String, Vec<[f64; 3]>)>,
}

fn line_error(line_no: usize, message: impl std::fmt::Display) -> String {
    format!("restart ASCII line {line_no}: {message}")
}

fn parse_f64(token: &str, line_no: usize, what: &str) -> Result<f64, String> {
    token
        .parse::<f64>()
        .map_err(|e| line_error(line_no, format!("invalid {what} `{token}`: {e}")))
}

fn parse_u64(token: &str, line_no: usize, what: &str) -> Result<u64, String> {
    token
        .parse::<u64>()
        .map_err(|e| line_error(line_no, format!("invalid {what} `{token}`: {e}")))
}

fn parse_vec3(line: &str, line_no: usize, what: &str) -> Result<[f64; 3], String> {
    let mut tokens = line.split_whitespace();
    match (tokens.next(), tokens.next(), tokens.next(), tokens.next()) {
        (Some(x), Some(y), Some(z), None) => Ok([
            parse_f64(x, line_no, what)?,
            parse_f64(y, line_no, what)?,
            parse_f64(z, line_no, what)?,
        ]),
        _ => Err(line_error(
            line_no,
            format!("expected exactly 3 values for {what}, found `{line}`"),
        )),
    }
}

fn parse_single<'a>(line: &'a str, line_no: usize, what: &str) -> Result<&'a str, String> {
    let mut tokens = line.split_whitespace();
    match (tokens.next(), tokens.next()) {
        (Some(token), None) => Ok(token),
        _ => Err(line_error(
            line_no,
            format!("expected exactly 1 value for {what}, found `{line}`"),
        )),
    }
}

/// `line` is `marker` optionally followed by one whitespace character and a
/// name; return the name (verbatim, untrimmed).
fn strip_marker<'a>(line: &'a str, marker: &str) -> Option<&'a str> {
    let rest = line.strip_prefix(marker)?;
    if rest.is_empty() {
        return Some("");
    }
    rest.strip_prefix(char::is_whitespace)
}

impl Section {
    fn open(kind: Kind, name: String, end_marker: &'static str, begin_line: usize) -> Self {
        Self {
            kind,
            name,
            end_marker,
            begin_line,
            vectors: Vec::new(),
            scalars: Vec::new(),
            ints: Vec::new(),
        }
    }

    /// Recognise a `BEGIN_*` marker; `raw` is the untrimmed line (extra array
    /// names keep their whitespace), `line` the trimmed one.
    fn try_begin(raw: &str, line: &str, line_no: usize) -> Option<Self> {
        for &(core, begin, end) in &CORE_SECTIONS {
            if line == begin {
                return Some(Self::open(Kind::Core(core), String::new(), end, line_no));
            }
        }
        let raw = raw.trim_start();
        if let Some(name) = strip_marker(raw, BEGIN_EXTRA_SCALAR) {
            return Some(Self::open(
                Kind::ExtraScalar,
                unescape_text(name),
                END_EXTRA_SCALAR,
                line_no,
            ));
        }
        if let Some(name) = strip_marker(raw, BEGIN_EXTRA_VECTOR) {
            return Some(Self::open(
                Kind::ExtraVector,
                unescape_text(name),
                END_EXTRA_VECTOR,
                line_no,
            ));
        }
        None
    }

    /// Parse one data line of this section.
    fn push_line(&mut self, line: &str, line_no: usize) -> Result<(), String> {
        match self.kind {
            Kind::Core(Core::Positions) => {
                self.vectors.push(parse_vec3(line, line_no, "position")?)
            }
            Kind::Core(Core::Velocities) => {
                self.vectors.push(parse_vec3(line, line_no, "velocity")?)
            }
            Kind::Core(Core::Forces) => self.vectors.push(parse_vec3(line, line_no, "force")?),
            Kind::Core(Core::Box) => self.vectors.push(parse_vec3(line, line_no, "box row")?),
            Kind::ExtraVector => {
                self.vectors
                    .push(parse_vec3(line, line_no, "extra vector entry")?)
            }
            Kind::Core(Core::Masses) => {
                let token = parse_single(line, line_no, "mass")?;
                self.scalars.push(parse_f64(token, line_no, "mass")?);
            }
            Kind::ExtraScalar => {
                let token = parse_single(line, line_no, "extra scalar entry")?;
                self.scalars
                    .push(parse_f64(token, line_no, "extra scalar entry")?);
            }
            Kind::Core(Core::Types) => {
                let token = parse_single(line, line_no, "type")?;
                let value = token
                    .parse::<u32>()
                    .map_err(|e| line_error(line_no, format!("invalid type `{token}`: {e}")))?;
                self.ints.push(value);
            }
        }
        Ok(())
    }

    /// Move the collected values into `parsed`.
    fn finish(self, parsed: &mut Parsed, end_line: usize) -> Result<(), String> {
        match self.kind {
            Kind::Core(core) => {
                let slot = &mut parsed.seen[core.index()];
                if *slot {
                    return Err(line_error(
                        self.begin_line,
                        format!("duplicate section (second `{}` block)", self.end_marker),
                    ));
                }
                *slot = true;
                match core {
                    Core::Positions => parsed.positions = self.vectors,
                    Core::Velocities => parsed.velocities = self.vectors,
                    Core::Forces => parsed.forces = self.vectors,
                    Core::Masses => parsed.masses = self.scalars,
                    Core::Types => parsed.types = self.ints,
                    Core::Box => {
                        let rows: [[f64; 3]; 3] =
                            <[[f64; 3]; 3]>::try_from(self.vectors.as_slice()).map_err(|_| {
                                line_error(
                                    end_line,
                                    format!(
                                        "BOX section must hold exactly 3 rows, found {}",
                                        self.vectors.len()
                                    ),
                                )
                            })?;
                        parsed.box_matrix = Some(rows);
                    }
                }
            }
            Kind::ExtraScalar => parsed.extra_scalars.push((self.name, self.scalars)),
            Kind::ExtraVector => parsed.extra_vectors.push((self.name, self.vectors)),
        }
        Ok(())
    }
}

fn is_marker(line: &str) -> bool {
    line.starts_with("BEGIN_") || line.starts_with("END_")
}

/// Parse one `KEY value` header line.
fn parse_header(
    raw: &str,
    line_no: usize,
    meta: &mut RestartMetadata,
    declared_n: &mut Option<usize>,
) -> Result<(), String> {
    let line = raw.trim_start();
    let (key, value) = line
        .split_once(char::is_whitespace)
        .unwrap_or((line.trim_end(), ""));
    match key {
        "VERSION" => meta.version = unescape_text(value),
        "CRATE" => meta.crate_name = unescape_text(value),
        "DESC" => meta.description = unescape_text(value),
        "TIMESTAMP" => meta.timestamp = parse_u64(value.trim(), line_no, "TIMESTAMP")?,
        "STEP" => meta.step = parse_u64(value.trim(), line_no, "STEP")?,
        "TIME" => meta.time = parse_f64(value.trim(), line_no, "TIME")?,
        "N_PARTICLES" => {
            let n = parse_u64(value.trim(), line_no, "N_PARTICLES")?;
            let n = usize::try_from(n).map_err(|_| {
                line_error(line_no, format!("N_PARTICLES {n} does not fit in memory"))
            })?;
            *declared_n = Some(n);
        }
        _ => {
            return Err(line_error(
                line_no,
                format!("unrecognised line `{}`", raw.trim()),
            ));
        }
    }
    Ok(())
}

/// Decode ASCII restart text (a leading UTF-8 byte-order mark is ignored).
pub(super) fn decode(text: &str) -> Result<RestartData, String> {
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut metadata = RestartMetadata::new("", 0, 0, 0.0, "", "");
    let mut declared_n: Option<usize> = None;
    let mut parsed = Parsed::default();
    let mut open: Option<Section> = None;

    for (index, raw) in text.lines().enumerate() {
        let line_no = index + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        match open.take() {
            Some(mut section) => {
                if line == section.end_marker {
                    section.finish(&mut parsed, line_no)?;
                } else if is_marker(line) {
                    return Err(line_error(
                        line_no,
                        format!(
                            "unexpected `{line}` inside the section opened on line {} (expected `{}`)",
                            section.begin_line, section.end_marker
                        ),
                    ));
                } else {
                    section.push_line(line, line_no)?;
                    open = Some(section);
                }
            }
            None => match Section::try_begin(raw, line, line_no) {
                Some(section) => open = Some(section),
                None => parse_header(raw, line_no, &mut metadata, &mut declared_n)?,
            },
        }
    }
    if let Some(section) = open {
        return Err(line_error(
            section.begin_line,
            format!("section is never closed (missing `{}`)", section.end_marker),
        ));
    }

    let n = declared_n.ok_or_else(|| "restart ASCII: missing N_PARTICLES header".to_string())?;
    if parsed.positions.len() != n {
        return Err(format!(
            "restart ASCII: N_PARTICLES is {n} but the POSITIONS section holds {} entries",
            parsed.positions.len()
        ));
    }
    let optional_lengths = [
        ("VELOCITIES", parsed.velocities.len()),
        ("FORCES", parsed.forces.len()),
        ("MASSES", parsed.masses.len()),
        ("TYPES", parsed.types.len()),
    ];
    for (name, len) in optional_lengths {
        if len != 0 && len != n {
            return Err(format!(
                "restart ASCII: N_PARTICLES is {n} but the {name} section holds {len} entries (expected 0 or {n})"
            ));
        }
    }

    Ok(RestartData {
        metadata,
        positions: parsed.positions,
        velocities: parsed.velocities,
        forces: parsed.forces,
        masses: parsed.masses,
        types: parsed.types,
        box_matrix: parsed.box_matrix.unwrap_or([[0.0; 3]; 3]),
        extra_scalars: parsed.extra_scalars,
        extra_vectors: parsed.extra_vectors,
    })
}
