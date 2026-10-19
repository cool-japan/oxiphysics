// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Internals of the fixed-column IGES reader ([`super::IgesReader`]).
//!
//! An IGES file is a sequence of 80-column records.  Columns 1-72 hold the
//! data, column 73 the section letter (`S` start, `G` global, `D` directory
//! entry, `P` parameter data, `T` terminate) and columns 74-80 the sequence
//! number of the record within its section.
//!
//! - The **global section** is free-format text.  Its first two parameters
//!   (`1H,` and `1H;` by default) define the parameter delimiter and the
//!   record delimiter that every free-format record of the file uses.
//! - The **directory section** has two records per entity.  The first holds
//!   the entity type (columns 1-8), the sequence number of the first parameter
//!   data record of the entity (9-16) and the status number (65-72); the
//!   second holds the number of parameter data records (25-32) and the label
//!   (57-64).
//! - The **parameter data section** holds, per entity, free-format data in
//!   columns 1-64: the entity type number, the parameters, and the record
//!   delimiter.  Columns 65-72 repeat the sequence number of the directory
//!   entry.  The records of an entity are found through the directory entry,
//!   not by their order in the section.
//! - Free-format data may contain Hollerith constants `nHtext`, whose text can
//!   itself contain the delimiters, and Fortran-style exponents (`1.5D+2`).

use super::types::{IgesEntity, IgesEntityType, LengthUnit};

/// Zero-based column of the section letter.
const SECTION_COLUMN: usize = 72;
/// Width of the free-format data area of a global-section record.
const GLOBAL_DATA_WIDTH: usize = 72;
/// Width of the free-format data area of a parameter-data record.
const PARAMETER_DATA_WIDTH: usize = 64;
/// Index (zero-based) of the units flag in the global parameters (parameter 14).
const UNITS_FLAG_INDEX: usize = 13;
/// Index (zero-based) of the units name in the global parameters (parameter 15).
const UNITS_NAME_INDEX: usize = 14;
/// Units flag assumed when the global section leaves parameter 14 blank (IGES
/// default: inches).
const DEFAULT_UNITS_FLAG: u32 = 1;

/// What the reader extracts from an IGES file.
#[derive(Debug)]
pub(super) struct ParsedIges {
    /// The entities with their parameter data.
    pub(super) entities: Vec<IgesEntity>,
    /// Global-section parameters, parameter 1 first, exactly as written.
    pub(super) global_params: Vec<String>,
    /// Length unit of the model, or `None` when the file has no global section.
    pub(super) units: Option<LengthUnit>,
}

/// The records of the sections that carry data, in file order.
#[derive(Debug, Default)]
struct Sections<'a> {
    global: Vec<&'a str>,
    directory: Vec<&'a str>,
    parameter: Vec<&'a str>,
}

/// File up the records by the section letter in column 73.  Start and
/// terminate records carry no data and short records cannot be classified;
/// both are skipped.
fn split_sections(content: &str) -> Sections<'_> {
    let mut sections = Sections::default();
    for line in content.lines() {
        match line.as_bytes().get(SECTION_COLUMN) {
            Some(b'G') => sections.global.push(line),
            Some(b'D') => sections.directory.push(line),
            Some(b'P') => sections.parameter.push(line),
            _ => {}
        }
    }
    sections
}

/// Columns `from..to` (zero-based, end exclusive) of `line`; empty when the
/// line is shorter or the range does not fall on character boundaries.
fn columns(line: &str, from: usize, to: usize) -> &str {
    line.get(from..to.min(line.len())).unwrap_or("")
}

/// The data area of a record, padded with blanks to `width` columns, so that a
/// Hollerith constant continued on the next record keeps its length.
fn padded_data(line: &str, width: usize) -> String {
    let data = columns(line, 0, width);
    format!("{data:<width$}")
}

/// Text of a Hollerith constant `nHtext`, `None` when `token` is not one.
fn hollerith_text(token: &str) -> Option<&str> {
    let digits = token.bytes().take_while(u8::is_ascii_digit).count();
    if digits == 0 {
        return None;
    }
    token.get(digits..)?.strip_prefix('H')
}

/// Split free-format data into its parameters.
///
/// Parameters are separated by `param_delim`; the data ends at `record_delim`.
/// A Hollerith constant (`nH` followed by exactly `n` characters) is kept
/// whole, so delimiters inside the text do not split it.  Every parameter is
/// returned as written, trimmed; an omitted parameter is the empty string.
///
/// Returns the parameters and whether the record delimiter was found.
fn tokenize(
    data: &[char],
    param_delim: char,
    record_delim: char,
) -> Result<(Vec<String>, bool), String> {
    let mut tokens = Vec::new();
    let mut pos = 0usize;
    loop {
        while data.get(pos) == Some(&' ') {
            pos += 1;
        }
        let start = pos;
        let digits = data[start..]
            .iter()
            .take_while(|c| c.is_ascii_digit())
            .count();
        let after_digits = start + digits;
        if digits > 0 && data.get(after_digits) == Some(&'H') {
            let length: usize = data[start..after_digits]
                .iter()
                .collect::<String>()
                .parse()
                .map_err(|_| "Hollerith constant with an unreadable length".to_string())?;
            let end = (after_digits + 1)
                .checked_add(length)
                .filter(|end| *end <= data.len())
                .ok_or_else(|| "Hollerith constant runs past the end of the data".to_string())?;
            tokens.push(data[start..end].iter().collect::<String>());
            pos = end;
            while data.get(pos) == Some(&' ') {
                pos += 1;
            }
            match data.get(pos) {
                Some(&c) if c == param_delim => {
                    pos += 1;
                    continue;
                }
                Some(&c) if c == record_delim => return Ok((tokens, true)),
                Some(&c) => {
                    return Err(format!(
                        "unexpected `{c}` after a Hollerith constant, expected a delimiter"
                    ));
                }
                None => return Ok((tokens, false)),
            }
        }
        while pos < data.len() && data[pos] != param_delim && data[pos] != record_delim {
            pos += 1;
        }
        let token: String = data[start..pos].iter().collect();
        tokens.push(token.trim().to_string());
        match data.get(pos) {
            Some(&c) if c == record_delim => return Ok((tokens, true)),
            Some(_) => pos += 1,
            None => return Ok((tokens, false)),
        }
    }
}

/// Delimiters and parameters of the global section.
struct GlobalSection {
    param_delim: char,
    record_delim: char,
    /// Parameter 1 (the parameter delimiter) first, as written.
    params: Vec<String>,
}

/// Parse the global section: the delimiter definitions, then the parameters.
fn parse_global(records: &[&str]) -> Result<GlobalSection, String> {
    let text: String = records
        .iter()
        .map(|line| padded_data(line, GLOBAL_DATA_WIDTH))
        .collect();
    let chars: Vec<char> = text.chars().collect();

    // Parameter 1: `1H<c>`, or omitted (default `,`).
    let (param_delim, mut pos) = match chars.get(..3) {
        Some(['1', 'H', c]) => (*c, 3usize),
        _ => (',', 0usize),
    };
    if chars.get(pos) != Some(&param_delim) {
        return Err(format!(
            "malformed global section: expected the parameter delimiter `{param_delim}` after parameter 1"
        ));
    }
    pos += 1;

    // Parameter 2: `1H<c>`, or omitted (default `;`).
    let record_delim = if chars.get(pos) == Some(&param_delim) {
        ';'
    } else {
        match chars.get(pos..pos + 3) {
            Some(['1', 'H', c]) => {
                pos += 3;
                *c
            }
            _ => {
                return Err(
                    "malformed global section: parameter 2 is not a record delimiter".to_string(),
                );
            }
        }
    };
    let first_param = format!("1H{param_delim}");
    let second_param = format!("1H{record_delim}");
    // The delimiter that follows parameter 2 (or its omission).
    match chars.get(pos) {
        Some(&c) if c == param_delim => pos += 1,
        Some(&c) if c == record_delim => {
            return Ok(GlobalSection {
                param_delim,
                record_delim,
                params: vec![first_param, second_param],
            });
        }
        _ => {
            return Err(
                "malformed global section: parameter 2 is not followed by a delimiter".to_string(),
            );
        }
    }

    let (tokens, _) = tokenize(&chars[pos..], param_delim, record_delim)
        .map_err(|e| format!("malformed global section: {e}"))?;
    let mut params = Vec::with_capacity(tokens.len() + 2);
    params.push(first_param);
    params.push(second_param);
    params.extend(tokens);
    Ok(GlobalSection {
        param_delim,
        record_delim,
        params,
    })
}

/// Length unit named by an IGES units name (global parameter 15).
fn unit_from_name(name: &str) -> Option<LengthUnit> {
    match name.trim().to_ascii_uppercase().as_str() {
        "IN" | "INCH" => Some(LengthUnit::Inch),
        "MM" => Some(LengthUnit::Millimeter),
        "FT" | "FOOT" | "FEET" => Some(LengthUnit::Foot),
        "MI" | "MILE" => Some(LengthUnit::Mile),
        "M" | "METER" | "METRE" => Some(LengthUnit::Meter),
        "KM" => Some(LengthUnit::Kilometer),
        "MIL" => Some(LengthUnit::Mil),
        "UM" => Some(LengthUnit::Micrometer),
        "CM" => Some(LengthUnit::Centimeter),
        "UIN" => Some(LengthUnit::Microinch),
        _ => None,
    }
}

/// Length unit selected by the units flag (global parameter 14).
fn units_from_global(params: &[String]) -> Result<LengthUnit, String> {
    let flag_text = params.get(UNITS_FLAG_INDEX).map_or("", |s| s.trim());
    let flag = if flag_text.is_empty() {
        DEFAULT_UNITS_FLAG
    } else {
        flag_text
            .parse::<u32>()
            .map_err(|_| format!("invalid units flag `{flag_text}` in the global section"))?
    };
    let unit = match flag {
        1 => Some(LengthUnit::Inch),
        2 => Some(LengthUnit::Millimeter),
        3 => {
            let name = params.get(UNITS_NAME_INDEX).map_or("", |s| s.as_str());
            unit_from_name(hollerith_text(name).unwrap_or(name))
        }
        4 => Some(LengthUnit::Foot),
        5 => Some(LengthUnit::Mile),
        6 => Some(LengthUnit::Meter),
        7 => Some(LengthUnit::Kilometer),
        8 => Some(LengthUnit::Mil),
        9 => Some(LengthUnit::Micrometer),
        10 => Some(LengthUnit::Centimeter),
        11 => Some(LengthUnit::Microinch),
        _ => None,
    };
    unit.ok_or_else(|| format!("unsupported units flag {flag} in the global section"))
}

/// The fields of one two-record directory entry that the reader uses.
struct DirectoryEntry {
    /// Sequence number of the first record of the entry (odd, from 1).
    sequence_number: usize,
    entity_type: u32,
    /// Sequence number of the first parameter data record.
    parameter_pointer: usize,
    /// Number of parameter data records; `None` when the directory leaves the
    /// field blank, in which case the records that point back at this entry
    /// are used.
    parameter_line_count: Option<usize>,
    status: u32,
    label: String,
}

/// Parse an unsigned integer field of a directory record; blank is `None`.
fn directory_number<T: std::str::FromStr>(
    field: &str,
    what: &str,
    sequence_number: usize,
) -> Result<Option<T>, String> {
    let field = field.trim();
    if field.is_empty() {
        return Ok(None);
    }
    field.parse::<T>().map(Some).map_err(|_| {
        format!("directory entry {sequence_number}: {what} `{field}` is not a valid number")
    })
}

fn parse_directory(records: &[&str]) -> Result<Vec<DirectoryEntry>, String> {
    if !records.len().is_multiple_of(2) {
        return Err(format!(
            "the directory section holds {} records, but every entry has two",
            records.len()
        ));
    }
    let mut entries = Vec::with_capacity(records.len() / 2);
    for (index, pair) in records.chunks_exact(2).enumerate() {
        let [first, second] = pair else {
            continue;
        };
        let sequence_number = 2 * index + 1;
        let entity_type =
            directory_number::<u32>(columns(first, 0, 8), "entity type", sequence_number)?
                .ok_or_else(|| {
                    format!("directory entry {sequence_number}: the entity type is blank")
                })?;
        let parameter_pointer = directory_number::<usize>(
            columns(first, 8, 16),
            "parameter data pointer",
            sequence_number,
        )?
        .filter(|pointer| *pointer >= 1)
        .ok_or_else(|| {
            format!("directory entry {sequence_number}: entity {entity_type} has no parameter data")
        })?;
        let status = directory_number::<u32>(columns(first, 64, 72), "status", sequence_number)?
            .unwrap_or(0);
        let parameter_line_count = directory_number::<usize>(
            columns(second, 24, 32),
            "parameter line count",
            sequence_number,
        )?
        .filter(|count| *count >= 1);
        entries.push(DirectoryEntry {
            sequence_number,
            entity_type,
            parameter_pointer,
            parameter_line_count,
            status,
            label: columns(second, 56, 64).trim().to_string(),
        });
    }
    Ok(entries)
}

/// Numeric value of a free-format parameter.
///
/// Accepts integers and reals with `E` or Fortran `D` exponents (`1.5D+2`);
/// anything else (`inf`, `nan`, stray letters) is rejected.
fn parse_real(token: &str) -> Option<f64> {
    let token = token.trim();
    let plausible = !token.is_empty()
        && token.bytes().all(|b| {
            b.is_ascii_digit() || matches!(b, b'+' | b'-' | b'.' | b'E' | b'e' | b'D' | b'd')
        });
    if !plausible {
        return None;
    }
    let normalized: String = token
        .chars()
        .map(|c| if matches!(c, 'D' | 'd') { 'E' } else { c })
        .collect();
    normalized.parse::<f64>().ok().filter(|v| v.is_finite())
}

/// Build the entity of one directory entry from its parameter data records.
fn build_entity(
    entry: &DirectoryEntry,
    parameter_records: &[&str],
    param_delim: char,
    record_delim: char,
) -> Result<IgesEntity, String> {
    let seq = entry.sequence_number;
    let first = entry.parameter_pointer - 1;
    let line_count = entry.parameter_line_count.unwrap_or_else(|| {
        // No count in the directory: take the records that name this entry.
        parameter_records
            .iter()
            .skip(first)
            .take_while(|line| {
                columns(line, PARAMETER_DATA_WIDTH, SECTION_COLUMN)
                    .trim()
                    .parse::<usize>()
                    .is_ok_and(|pointer| pointer == seq)
            })
            .count()
    });
    if line_count == 0 {
        return Err(format!(
            "directory entry {seq}: entity {} has no parameter line count and no parameter data record points back at it",
            entry.entity_type
        ));
    }
    let last = first
        .checked_add(line_count)
        .filter(|last| *last <= parameter_records.len())
        .ok_or_else(|| {
            format!(
                "directory entry {seq}: its {line_count} parameter data records starting at record {} run past the {} records of the section",
                entry.parameter_pointer,
                parameter_records.len()
            )
        })?;
    let text: String = parameter_records[first..last]
        .iter()
        .map(|line| padded_data(line, PARAMETER_DATA_WIDTH))
        .collect();
    let chars: Vec<char> = text.chars().collect();
    let (tokens, terminated) = tokenize(&chars, param_delim, record_delim)
        .map_err(|e| format!("directory entry {seq}: parameter data: {e}"))?;
    if !terminated {
        return Err(format!(
            "directory entry {seq}: the parameter data lacks the record delimiter `{record_delim}`"
        ));
    }

    let (type_token, parameter_tokens) = tokens
        .split_first()
        .ok_or_else(|| format!("directory entry {seq}: the parameter data is empty"))?;
    let written_type = type_token.trim().parse::<u32>().map_err(|_| {
        format!("directory entry {seq}: the parameter data starts with `{type_token}`, not an entity type number")
    })?;
    if written_type != entry.entity_type {
        return Err(format!(
            "directory entry {seq}: the directory says entity type {} but the parameter data starts with {written_type}",
            entry.entity_type
        ));
    }

    let mut parameters = Vec::with_capacity(parameter_tokens.len());
    for (position, token) in parameter_tokens.iter().enumerate() {
        if token.is_empty() || hollerith_text(token).is_some() {
            // Omitted parameter or text constant: no numeric value.
            parameters.push(f64::NAN);
            continue;
        }
        let value = parse_real(token).ok_or_else(|| {
            format!(
                "directory entry {seq}: parameter {} `{token}` of entity {} is not a number",
                position + 1,
                entry.entity_type
            )
        })?;
        parameters.push(value);
    }

    let mut entity = IgesEntity::new(seq, IgesEntityType::from_type_number(entry.entity_type));
    entity.parameter_pointer = entry.parameter_pointer;
    entity.status = entry.status;
    entity.label = entry.label.clone();
    entity.parameters = parameters;
    entity.raw_parameters = parameter_tokens.to_vec();
    Ok(entity)
}

/// Parse the text of an IGES file.
///
/// # Errors
///
/// Returns a description of the first defect: a malformed global section, an
/// odd number of directory records, a directory field that is not a number, a
/// parameter data span outside the section, parameter data that starts with a
/// different entity type than its directory entry or lacks the record
/// delimiter, or a parameter that is neither a number, a Hollerith constant
/// nor omitted.
pub(super) fn parse(content: &str) -> Result<ParsedIges, String> {
    let sections = split_sections(content);

    let (param_delim, record_delim, global_params, units) = if sections.global.is_empty() {
        (',', ';', Vec::new(), None)
    } else {
        let global = parse_global(&sections.global)?;
        let units = units_from_global(&global.params)?;
        (
            global.param_delim,
            global.record_delim,
            global.params,
            Some(units),
        )
    };

    let directory = parse_directory(&sections.directory)?;
    let entities = directory
        .iter()
        .map(|entry| build_entity(entry, &sections.parameter, param_delim, record_delim))
        .collect::<Result<Vec<_>, _>>()?;

    Ok(ParsedIges {
        entities,
        global_params,
        units,
    })
}

#[cfg(test)]
mod tests {
    use super::super::types::IgesReader;
    use super::*;

    /// One 80-column record: 72 data columns, the section letter, and the
    /// sequence number.
    fn record(data: &str, section: char, seq: usize) -> String {
        format!("{data:<72}{section}{seq:>7}")
    }

    /// A parameter data record: 64 data columns, the directory pointer, `P`,
    /// and the sequence number.
    fn parameter_record(data: &str, de_pointer: usize, seq: usize) -> String {
        format!("{data:<64}{de_pointer:>8}P{seq:>7}")
    }

    fn directory_records(
        entity_type: u32,
        pd_pointer: usize,
        line_count: usize,
        label: &str,
        status: &str,
        seq: usize,
    ) -> [String; 2] {
        let first = format!(
            "{entity_type:>8}{pd_pointer:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{status:>8}",
            0, 0, 0, 0, 0, 0
        );
        let second = format!(
            "{entity_type:>8}{:>8}{:>8}{line_count:>8}{:>8}{:>8}{:>8}{label:>8}{:>8}",
            0, 0, 0, 0, 0, 0
        );
        [record(&first, 'D', seq), record(&second, 'D', seq + 1)]
    }

    fn global_section(units_flag: &str, units_name: &str) -> Vec<String> {
        let text = format!(
            "1H,,1H;,4HTEST,8Htest.igs,10HOxiPhysics,5H1.0.0,32,38,6,308,15,4HTEST,1.0,{units_flag},{units_name},1,0.01,15H20260101.000000,0.001,100.0,6Hauthor,3Horg,11,0;"
        );
        text.as_bytes()
            .chunks(72)
            .enumerate()
            .map(|(i, chunk)| record(&String::from_utf8_lossy(chunk), 'G', i + 1))
            .collect()
    }

    /// `(entity type, parameter text)` pairs, each parameter text in as many
    /// 64-column records as it needs, assembled into a complete IGES file.
    fn iges_file(units_flag: &str, units_name: &str, entities: &[(u32, &str)]) -> String {
        let mut lines = vec![record("generated for a test", 'S', 1)];
        lines.extend(global_section(units_flag, units_name));
        let mut directory = Vec::new();
        let mut parameters = Vec::new();
        for (index, (entity_type, text)) in entities.iter().enumerate() {
            let pd_pointer = parameters.len() + 1;
            let chunks: Vec<String> = text
                .as_bytes()
                .chunks(64)
                .map(|c| String::from_utf8_lossy(c).into_owned())
                .collect();
            directory.extend(directory_records(
                *entity_type,
                pd_pointer,
                chunks.len(),
                "",
                "0",
                2 * index + 1,
            ));
            for chunk in &chunks {
                parameters.push(parameter_record(chunk, 2 * index + 1, parameters.len() + 1));
            }
        }
        lines.extend(directory);
        lines.extend(parameters);
        lines.push(record("S      1G      3D      2P      1", 'T', 1));
        lines.join("\n")
    }

    // ── tokenizer ──────────────────────────────────────────────────────────

    fn tokens(text: &str) -> (Vec<String>, bool) {
        let chars: Vec<char> = text.chars().collect();
        tokenize(&chars, ',', ';').unwrap()
    }

    #[test]
    fn test_tokenize_plain_parameters() {
        assert_eq!(
            tokens("110,0.,0.,0.,1.,1.,1.;   "),
            (
                vec!["110", "0.", "0.", "0.", "1.", "1.", "1."]
                    .into_iter()
                    .map(String::from)
                    .collect::<Vec<_>>(),
                true
            )
        );
    }

    #[test]
    fn test_tokenize_keeps_omitted_parameters_and_hollerith_constants_whole() {
        let (parts, terminated) = tokens("406,,3HA,B,4H;,;x,2;");
        assert!(terminated);
        assert_eq!(parts, vec!["406", "", "3HA,B", "4H;,;x", "2"]);
        // A Hollerith constant may be followed by blanks before the delimiter.
        let (parts, _) = tokens("212,2HAB   ,7;");
        assert_eq!(parts, vec!["212", "2HAB", "7"]);
    }

    #[test]
    fn test_tokenize_reports_a_missing_record_delimiter_and_bad_constants() {
        assert!(!tokens("110,1.,2.").1);
        let chars: Vec<char> = "212,9HAB;".chars().collect();
        assert!(tokenize(&chars, ',', ';').is_err());
        let chars: Vec<char> = "212,2HABC,1;".chars().collect();
        assert!(tokenize(&chars, ',', ';').is_err());
    }

    #[test]
    fn test_tokenize_honours_other_delimiters() {
        let chars: Vec<char> = "110/1.5/2.5$".chars().collect();
        let (parts, terminated) = tokenize(&chars, '/', '$').unwrap();
        assert!(terminated);
        assert_eq!(parts, vec!["110", "1.5", "2.5"]);
    }

    // ── numbers ────────────────────────────────────────────────────────────

    #[test]
    fn test_parse_real_accepts_the_iges_number_forms() {
        assert_eq!(parse_real("0."), Some(0.0));
        assert_eq!(parse_real(".5"), Some(0.5));
        assert_eq!(parse_real("1.D0"), Some(1.0));
        assert_eq!(parse_real("-1.5D+2"), Some(-150.0));
        assert_eq!(parse_real("2.5d-1"), Some(0.25));
        assert_eq!(parse_real("1E3"), Some(1000.0));
        assert_eq!(parse_real("+12"), Some(12.0));
        assert_eq!(parse_real("  3.25 "), Some(3.25));
    }

    #[test]
    fn test_parse_real_rejects_everything_else() {
        for text in [
            "", "abc", "inf", "NaN", "infinity", "1.2.3", "D", "+", "1,5", "0x10",
        ] {
            assert_eq!(parse_real(text), None, "{text:?}");
        }
    }

    // ── global section ─────────────────────────────────────────────────────

    #[test]
    fn test_global_section_delimiters_and_parameters() {
        let lines = global_section("2", "2HMM");
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let global = parse_global(&refs).unwrap();
        assert_eq!((global.param_delim, global.record_delim), (',', ';'));
        assert_eq!(global.params[0], "1H,");
        assert_eq!(global.params[1], "1H;");
        assert_eq!(global.params[2], "4HTEST");
        assert_eq!(global.params[13], "2");
        assert_eq!(global.params[14], "2HMM");
        assert_eq!(global.params[17], "15H20260101.000000");
    }

    #[test]
    fn test_global_section_with_omitted_delimiter_parameters_uses_the_defaults() {
        let lines = [record(",,4HTEST,8Htest.igs;", 'G', 1)];
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let global = parse_global(&refs).unwrap();
        assert_eq!((global.param_delim, global.record_delim), (',', ';'));
        assert_eq!(global.params[2], "4HTEST");
    }

    #[test]
    fn test_global_section_with_custom_delimiters() {
        let lines = [record("1H//1H$/4HTEST/2/3.5$", 'G', 1)];
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let global = parse_global(&refs).unwrap();
        assert_eq!((global.param_delim, global.record_delim), ('/', '$'));
        assert_eq!(&global.params[2..], ["4HTEST", "2", "3.5"]);
    }

    #[test]
    fn test_units_flags() {
        let cases = [
            ("1", "2HIN", LengthUnit::Inch),
            ("2", "2HMM", LengthUnit::Millimeter),
            ("4", "2HFT", LengthUnit::Foot),
            ("5", "2HMI", LengthUnit::Mile),
            ("6", "1HM", LengthUnit::Meter),
            ("7", "2HKM", LengthUnit::Kilometer),
            ("8", "3HMIL", LengthUnit::Mil),
            ("9", "2HUM", LengthUnit::Micrometer),
            ("10", "2HCM", LengthUnit::Centimeter),
            ("11", "3HUIN", LengthUnit::Microinch),
        ];
        for (flag, name, expected) in cases {
            let parsed = parse(&iges_file(flag, name, &[])).unwrap();
            assert_eq!(parsed.units, Some(expected), "flag {flag}");
        }
    }

    #[test]
    fn test_units_flag_three_reads_the_units_name() {
        let parsed = parse(&iges_file("3", "2HCM", &[])).unwrap();
        assert_eq!(parsed.units, Some(LengthUnit::Centimeter));
        assert!(parse(&iges_file("3", "4HFURL", &[])).is_err());
    }

    #[test]
    fn test_blank_units_flag_defaults_to_inches_and_bad_flags_are_rejected() {
        let parsed = parse(&iges_file("", "", &[])).unwrap();
        assert_eq!(parsed.units, Some(LengthUnit::Inch));
        let message = parse(&iges_file("99", "", &[])).unwrap_err();
        assert!(message.contains("units flag 99"), "{message}");
        let message = parse(&iges_file("two", "", &[])).unwrap_err();
        assert!(message.contains("two"), "{message}");
    }

    #[test]
    fn test_file_without_a_global_section_has_no_units() {
        let parsed = parse("").unwrap();
        assert!(parsed.entities.is_empty());
        assert_eq!(parsed.units, None);
    }

    // ── entities ───────────────────────────────────────────────────────────

    #[test]
    fn test_line_entity_parameters_exclude_the_entity_type() {
        let parsed = parse(&iges_file("2", "2HMM", &[(110, "110,0.,0.,0.,1.,1.,1.;")])).unwrap();
        assert_eq!(parsed.entities.len(), 1);
        let entity = &parsed.entities[0];
        assert_eq!(entity.entity_type, IgesEntityType::Line);
        assert_eq!(entity.parameters, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
        assert_eq!(
            entity.raw_parameters,
            vec!["0.", "0.", "0.", "1.", "1.", "1."]
        );
        assert_eq!(entity.sequence_number, 1);
        assert_eq!(entity.parameter_pointer, 1);
    }

    #[test]
    fn test_reader_yields_the_line_of_a_110_entity() {
        let mut reader = IgesReader::new();
        reader
            .parse(&iges_file("2", "2HMM", &[(110, "110,0.,0.,0.,1.,1.,1.;")]))
            .unwrap();
        assert_eq!(reader.entity_count(), 1);
        assert_eq!(reader.units, LengthUnit::Millimeter);
        assert_eq!(
            reader.extract_lines(),
            vec![([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])]
        );
        assert!(reader.extract_points().is_empty());
        assert_eq!(reader.find_entities(IgesEntityType::Line).len(), 1);
        assert_eq!(reader.global_params[0], "1H,");
    }

    #[test]
    fn test_reader_extracts_points_lines_and_arcs_from_a_mixed_file() {
        let mut reader = IgesReader::new();
        reader
            .parse(&iges_file(
                "6",
                "1HM",
                &[
                    (116, "116,1.5D0,-2.,3.;"),
                    (100, "100,2.5,0.,0.,3.,4.,5.,0.;"),
                    (110, "110,1.,2.,3.,4.,5.,6.;"),
                    // An arc without its end point is not an arc.
                    (100, "100,2.5,0.,0.,3.,4.,5.;"),
                    // A point with an omitted coordinate is skipped.
                    (116, "116,1.,,3.;"),
                ],
            ))
            .unwrap();
        assert_eq!(reader.entity_count(), 5);
        assert_eq!(reader.units, LengthUnit::Meter);
        assert_eq!(reader.extract_points(), vec![[1.5, -2.0, 3.0]]);
        assert_eq!(
            reader.extract_lines(),
            vec![([1.0, 2.0, 3.0], [4.0, 5.0, 6.0])]
        );
        let arcs = reader.extract_circular_arcs();
        assert_eq!(arcs.len(), 1);
        assert_eq!(arcs[0].z_displacement, 2.5);
        assert_eq!(arcs[0].center, [0.0, 0.0]);
        assert_eq!(arcs[0].start, [3.0, 4.0]);
        assert_eq!(arcs[0].end, [5.0, 0.0]);
        assert!((arcs[0].radius() - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_reader_state_follows_the_last_parse() {
        let mut reader = IgesReader::new();
        reader
            .parse(&iges_file("1", "2HIN", &[(110, "110,0.,0.,0.,1.,1.,1.;")]))
            .unwrap();
        assert_eq!(reader.units, LengthUnit::Inch);
        reader.parse(&iges_file("10", "2HCM", &[])).unwrap();
        assert_eq!(reader.units, LengthUnit::Centimeter);
        assert_eq!(reader.entity_count(), 0);
        // A failing parse leaves no stale entities behind.
        assert!(
            reader
                .parse(&iges_file("6", "1HM", &[(110, "116,0.,0.,0.;")]))
                .is_err()
        );
        assert_eq!(reader.entity_count(), 0);
    }

    #[test]
    fn test_d_exponents_are_read() {
        let parsed = parse(&iges_file("6", "1HM", &[(116, "116,1.5D+2,-2.D0,.5D1,0;")])).unwrap();
        assert_eq!(parsed.entities[0].parameters, vec![150.0, -2.0, 5.0, 0.0]);
    }

    #[test]
    fn test_parameter_data_is_found_through_the_directory_pointer() {
        // The two entities' parameter records are in the opposite order of
        // their directory entries.
        let mut lines = vec![record("pointer test", 'S', 1)];
        lines.extend(global_section("6", "1HM"));
        lines.extend(directory_records(116, 2, 1, "POINT", "1", 1));
        lines.extend(directory_records(110, 1, 1, "LINE", "0", 3));
        lines.push(parameter_record("110,0.,0.,0.,1.,2.,3.;", 3, 1));
        lines.push(parameter_record("116,7.,8.,9.;", 1, 2));
        let parsed = parse(&lines.join("\n")).unwrap();
        assert_eq!(parsed.entities.len(), 2);
        assert_eq!(parsed.entities[0].entity_type, IgesEntityType::Point);
        assert_eq!(parsed.entities[0].parameters, vec![7.0, 8.0, 9.0]);
        assert_eq!(parsed.entities[0].label, "POINT");
        assert_eq!(parsed.entities[0].status, 1);
        assert_eq!(parsed.entities[1].entity_type, IgesEntityType::Line);
        assert_eq!(
            parsed.entities[1].parameters,
            vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0]
        );
        assert_eq!(parsed.entities[1].sequence_number, 3);
    }

    #[test]
    fn test_parameters_spanning_several_records_are_joined() {
        let text = format!(
            "126,{};",
            (1..=40)
                .map(|i| format!("{i}.25"))
                .collect::<Vec<_>>()
                .join(",")
        );
        assert!(text.len() > 128, "the data must need several records");
        let parsed = parse(&iges_file("6", "1HM", &[(126, &text)])).unwrap();
        let expected: Vec<f64> = (1..=40).map(|i| f64::from(i) + 0.25).collect();
        assert_eq!(parsed.entities[0].parameters, expected);
    }

    #[test]
    fn test_hollerith_and_omitted_parameters_have_no_numeric_value() {
        let parsed = parse(&iges_file("6", "1HM", &[(406, "406,2,6HA,B;C ,,3.5;")])).unwrap();
        let entity = &parsed.entities[0];
        assert_eq!(entity.raw_parameters, vec!["2", "6HA,B;C ", "", "3.5"]);
        assert_eq!(entity.parameters.len(), 4);
        assert_eq!(entity.parameters[0], 2.0);
        assert!(entity.parameters[1].is_nan());
        assert!(entity.parameters[2].is_nan());
        assert_eq!(entity.parameters[3], 3.5);
    }

    #[test]
    fn test_unknown_entity_types_keep_their_parameters() {
        let parsed = parse(&iges_file("6", "1HM", &[(314, "314,10.,20.,30.;")])).unwrap();
        assert_eq!(parsed.entities[0].entity_type, IgesEntityType::Unknown(314));
        assert_eq!(parsed.entities[0].parameters, vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn test_non_numeric_parameter_is_an_error() {
        let message =
            parse(&iges_file("6", "1HM", &[(110, "110,0.,zero,0.,1.,1.,1.;")])).unwrap_err();
        assert!(message.contains("zero"), "{message}");
        assert!(message.contains("parameter 2"), "{message}");
    }

    #[test]
    fn test_type_mismatch_between_directory_and_parameter_data_is_an_error() {
        let message = parse(&iges_file("6", "1HM", &[(110, "116,0.,0.,0.;")])).unwrap_err();
        assert!(
            message.contains("110") && message.contains("116"),
            "{message}"
        );
    }

    #[test]
    fn test_missing_record_delimiter_is_an_error() {
        let message = parse(&iges_file("6", "1HM", &[(110, "110,0.,0.,0.,1.,1.,1.")])).unwrap_err();
        assert!(message.contains("record delimiter"), "{message}");
    }

    #[test]
    fn test_blank_line_count_falls_back_to_the_records_that_point_back() {
        // Entity 1 has a count of 1 in the directory, entity 2 leaves it blank
        // and spans two records, each naming directory entry 3.
        let mut lines = vec![record("x", 'S', 1)];
        lines.extend(global_section("6", "1HM"));
        lines.extend(directory_records(116, 1, 1, "", "0", 1));
        let mut blank = directory_records(126, 2, 1, "", "0", 3);
        blank[1] = record(&format!("{:>8}{:>8}{:>8}{:>8}", 126, 0, 0, ""), 'D', 4);
        lines.extend(blank);
        lines.push(parameter_record("116,1.,2.,3.;", 1, 1));
        lines.push(parameter_record("126,1.5,2.5,", 3, 2));
        lines.push(parameter_record("3.5,4.5;", 3, 3));
        let parsed = parse(&lines.join("\n")).unwrap();
        assert_eq!(parsed.entities[0].parameters, vec![1.0, 2.0, 3.0]);
        assert_eq!(parsed.entities[1].parameters, vec![1.5, 2.5, 3.5, 4.5]);
    }

    #[test]
    fn test_blank_line_count_without_back_pointers_is_an_error() {
        let mut lines = vec![record("x", 'S', 1)];
        lines.extend(global_section("6", "1HM"));
        let mut entry = directory_records(110, 1, 1, "", "0", 1);
        entry[1] = record(&format!("{:>8}{:>8}{:>8}{:>8}", 110, 0, 0, ""), 'D', 2);
        lines.extend(entry);
        lines.push(parameter_record("110,0.,0.,0.,1.,1.,1.;", 0, 1));
        let message = parse(&lines.join("\n")).unwrap_err();
        assert!(message.contains("no parameter line count"), "{message}");
    }

    #[test]
    fn test_directory_defects_are_errors() {
        // A parameter span beyond the section.
        let mut lines = vec![record("x", 'S', 1)];
        lines.extend(global_section("6", "1HM"));
        lines.extend(directory_records(110, 5, 1, "", "0", 1));
        lines.push(parameter_record("110,0.,0.,0.,1.,1.,1.;", 1, 1));
        let message = parse(&lines.join("\n")).unwrap_err();
        assert!(message.contains("run past"), "{message}");

        // An odd number of directory records.
        let mut lines = vec![record("x", 'S', 1)];
        lines.extend(global_section("6", "1HM"));
        lines.push(directory_records(110, 1, 1, "", "0", 1)[0].clone());
        assert!(parse(&lines.join("\n")).unwrap_err().contains("two"));

        // A non-numeric entity type.
        let mut lines = vec![record("x", 'S', 1)];
        lines.extend(global_section("6", "1HM"));
        let mut entry = directory_records(110, 1, 1, "", "0", 1);
        entry[0] = record(&format!("{:>8}{:>8}", "abc", 1), 'D', 1);
        lines.extend(entry);
        assert!(parse(&lines.join("\n")).unwrap_err().contains("abc"));
    }

    #[test]
    fn test_parameter_data_honours_the_delimiters_of_the_global_section() {
        // `/` separates the parameters and `$` ends the records.
        let global_text = "1H//1H$/4HTEST/8Htest.igs/10HOxiPhysics/5H1.0.0/32/38/6/308/15/4HTEST/1.0/6/1HM/1/0.01/15H20260101.000000/0.001/100.0/6Hauthor/3Horg/11/0$";
        let mut lines = vec![record("custom delimiters", 'S', 1)];
        for (index, chunk) in global_text.as_bytes().chunks(72).enumerate() {
            lines.push(record(&String::from_utf8_lossy(chunk), 'G', index + 1));
        }
        lines.extend(directory_records(110, 1, 1, "", "0", 1));
        lines.push(parameter_record("110/0./0./0./1./2./3.$", 1, 1));
        let parsed = parse(&lines.join("\n")).unwrap();
        assert_eq!(parsed.units, Some(LengthUnit::Meter));
        assert_eq!(
            parsed.entities[0].parameters,
            vec![0.0, 0.0, 0.0, 1.0, 2.0, 3.0]
        );
        // With the default delimiters the same data is not a record at all.
        let mut wrong = lines.clone();
        wrong.pop();
        wrong.push(parameter_record("110,0.,0.,0.,1.,2.,3.;", 1, 1));
        assert!(parse(&wrong.join("\n")).is_err());
    }

    #[test]
    fn test_windows_line_endings_are_accepted() {
        let unix = iges_file("2", "2HMM", &[(110, "110,0.,0.,0.,1.,1.,1.;")]);
        let windows = unix.replace('\n', "\r\n");
        assert_ne!(unix, windows);
        let parsed = parse(&windows).unwrap();
        assert_eq!(parsed.units, Some(LengthUnit::Millimeter));
        assert_eq!(parsed.entities.len(), 1);
        assert_eq!(
            parsed.entities[0].parameters,
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        );
    }

    #[test]
    fn test_short_and_unclassifiable_records_are_ignored() {
        let mut text = iges_file("6", "1HM", &[(110, "110,0.,0.,0.,1.,1.,1.;")]);
        text.push_str("\nshort\n\u{3042}\u{3042}\u{3042}\n");
        // Column 73 falls on a multi-byte character, so no section letter.
        text.push_str(&format!("{}\u{3042}D      9\n", "\u{e9}".repeat(70)));
        let parsed = parse(&text).unwrap();
        assert_eq!(parsed.entities.len(), 1);
        assert_eq!(parsed.entities[0].parameters.len(), 6);
    }

    #[test]
    fn test_multi_byte_character_inside_a_directory_field_is_an_error_not_a_panic() {
        let mut lines = vec![record("x", 'S', 1)];
        lines.extend(global_section("6", "1HM"));
        let mut entry = directory_records(110, 1, 1, "", "0", 1);
        // Seven blanks and an `\u{e9}` straddling the end of the type field;
        // the record is padded to exactly 72 bytes so that `D` stays in column 73.
        let data = format!("       \u{e9}{:>8}", 1);
        entry[0] = format!("{data}{}D{:>7}", " ".repeat(72 - data.len()), 1);
        lines.extend(entry);
        lines.push(parameter_record("110,0.,0.,0.,1.,1.,1.;", 1, 1));
        let message = parse(&lines.join("\n")).unwrap_err();
        assert!(message.contains("directory entry 1"), "{message}");
    }
}
