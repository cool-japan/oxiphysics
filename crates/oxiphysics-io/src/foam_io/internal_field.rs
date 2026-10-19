// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! The `internalField` entry of OpenFOAM `vol*Field` files.
//!
//! A field file holds one `internalField` entry, the values in the cells, and
//! a `boundaryField` dictionary whose patches carry `value` entries of their
//! own, many of them vectors.  The internal field therefore has to be located
//! by its keyword and read up to the end of *its own* list: everything between
//! the first `(` and the last `)` of the file also contains the boundary
//! values.
//!
//! The spellings that OpenFOAM writes are understood:
//!
//! ```text
//! internalField   uniform (1 0 0);
//! internalField   nonuniform List<vector> 3((1 0 0) (2 0 0) (3 0 0));
//! internalField   nonuniform List<scalar>
//! 4
//! (
//! 1
//! 2
//! 3
//! 4
//! )
//! ;
//! internalField   nonuniform List<scalar> 100{0.5};
//! ```
//!
//! Comments (`// ...` and `/* ... */`) are ignored, and only an entry at the
//! top level of the file counts, so a `$internalField` reference inside a
//! patch is not mistaken for the entry itself.

/// Largest number of values a repeated `N{value}` entry may expand to.
///
/// The count of a repeated list is taken from the file and no data backs it
/// up, so it is bounded instead of being trusted.
const MAX_REPEATED_ENTRIES: usize = 1 << 26;

/// An `internalField` entry as it is written in the file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum InternalField {
    /// `uniform <value>;`: the text of the single value.
    Uniform(String),
    /// `nonuniform [List<type>] [N] ( ... )`: the text between the
    /// parentheses and the element count the file declares, if it does.
    Listed {
        declared: Option<usize>,
        body: String,
    },
    /// `nonuniform [List<type>] N{ value }`: one value repeated `count` times.
    Repeated { count: usize, value: String },
}

/// Locate the top-level `internalField` entry of `content`.
///
/// Returns `Ok(None)` when the file has none, in which case the caller is
/// looking at a bare list.
///
/// # Errors
///
/// Fails when the entry is not `uniform` or `nonuniform`, when its value is
/// not terminated, or when a parenthesis or brace of its list is not closed.
pub(super) fn internal_field_block(content: &str) -> Result<Option<InternalField>, String> {
    let clean = strip_comments(content);
    match find_entry(&clean, "internalField") {
        Some(after_keyword) => parse_entry(clean.get(after_keyword..).unwrap_or("")).map(Some),
        None => Ok(None),
    }
}

/// Read a bare `[N] ( ... )` list, as it appears in a `points` file or in
/// text that holds nothing but the list.
///
/// A `FoamFile` header in front of the list is skipped.  The list ends at the
/// parenthesis that matches its opening one, so text behind it does not become
/// part of it.  Returns the element count written in front of the list, if
/// there is one, and the text between the parentheses.
///
/// # Errors
///
/// Fails when there is no list or when its parenthesis is not closed.
pub(super) fn bare_list(content: &str) -> Result<(Option<usize>, String), String> {
    let clean = strip_comments(content);
    let search_from = match clean.find("FoamFile") {
        Some(pos) => {
            let brace = clean
                .get(pos..)
                .and_then(|rest| rest.find('{'))
                .map(|offset| pos + offset)
                .ok_or_else(|| "No opening brace for FoamFile".to_string())?;
            let close = matching_close(clean.get(brace..).unwrap_or(""), '{', '}')
                .ok_or_else(|| "No closing brace for FoamFile".to_string())?;
            brace + close + 1
        }
        None => 0,
    };
    let rest = clean.get(search_from..).unwrap_or("");
    let open = rest
        .find('(')
        .ok_or_else(|| "Could not find list block '(' … ')'".to_string())?;
    let close = matching_close(rest.get(open..).unwrap_or(""), '(', ')')
        .ok_or_else(|| "Unmatched '(' in list block".to_string())?;
    let declared = rest
        .get(..open)
        .and_then(|before| before.split_whitespace().next_back())
        .and_then(|token| token.parse::<usize>().ok());
    let body = rest.get(open + 1..open + close).unwrap_or("").to_string();
    Ok((declared, body))
}

impl InternalField {
    /// The scalar values of the entry.
    ///
    /// A `uniform` entry yields its single value.
    ///
    /// # Errors
    ///
    /// Fails on a token that is not a number and when the element count the
    /// file declares differs from the number of values it holds.
    pub(super) fn into_scalars(self) -> Result<Vec<f64>, String> {
        match self {
            Self::Uniform(text) => Ok(vec![parse_scalar(&text, "uniform scalar")?]),
            Self::Listed { declared, body } => parse_scalar_list(declared, &body),
            Self::Repeated { count, value } => {
                repeated(parse_scalar(&value, "repeated scalar")?, count)
            }
        }
    }

    /// The vector values of the entry.
    ///
    /// A `uniform` entry yields its single value.
    ///
    /// # Errors
    ///
    /// Fails on an element that is not an `(x y z)` triple and when the
    /// element count the file declares differs from the number of vectors it
    /// holds.
    pub(super) fn into_vectors(self) -> Result<Vec<[f64; 3]>, String> {
        match self {
            Self::Uniform(text) => Ok(vec![parse_vector(&text)?]),
            Self::Listed { declared, body } => parse_vector_list(declared, &body),
            Self::Repeated { count, value } => repeated(parse_vector(&value)?, count),
        }
    }
}

/// Parse the whitespace separated numbers of a scalar list.
pub(super) fn parse_scalar_list(declared: Option<usize>, body: &str) -> Result<Vec<f64>, String> {
    let values = body
        .split_whitespace()
        .map(|token| parse_scalar(token, "scalar"))
        .collect::<Result<Vec<f64>, String>>()?;
    check_declared(declared, values.len(), "values")?;
    Ok(values)
}

/// Parse a list of `(x y z)` triples.
pub(super) fn parse_vector_list(
    declared: Option<usize>,
    body: &str,
) -> Result<Vec<[f64; 3]>, String> {
    let mut vectors = Vec::new();
    let mut rest = body.trim_start();
    while !rest.is_empty() {
        if !rest.starts_with('(') {
            let stray: String = rest.chars().take_while(|c| !c.is_whitespace()).collect();
            return Err(format!(
                "Expected '(x y z)' in a vector list, found {stray:?}"
            ));
        }
        let close = rest
            .find(')')
            .ok_or_else(|| "Unmatched '(' in a vector list".to_string())?;
        vectors.push(parse_vector(rest.get(..=close).unwrap_or(""))?);
        rest = rest.get(close + 1..).unwrap_or("").trim_start();
    }
    check_declared(declared, vectors.len(), "vectors")?;
    Ok(vectors)
}

fn check_declared(declared: Option<usize>, found: usize, what: &str) -> Result<(), String> {
    match declared {
        Some(count) if count != found => Err(format!(
            "The list declares {count} {what} but holds {found}"
        )),
        _ => Ok(()),
    }
}

/// `count` copies of `value`, as `count{value}` stands for.
fn repeated<T: Clone>(value: T, count: usize) -> Result<Vec<T>, String> {
    if count > MAX_REPEATED_ENTRIES {
        return Err(format!(
            "A repeated list of {count} values is larger than the limit of {MAX_REPEATED_ENTRIES}"
        ));
    }
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| format!("Cannot allocate a repeated list of {count} values"))?;
    values.resize(count, value);
    Ok(values)
}

fn parse_scalar(token: &str, what: &str) -> Result<f64, String> {
    token
        .trim()
        .parse::<f64>()
        .map_err(|_| format!("Cannot parse {what}: {token:?}"))
}

/// Parse one `(x y z)` triple.
fn parse_vector(text: &str) -> Result<[f64; 3], String> {
    let trimmed = text.trim();
    let inner = trimmed
        .strip_prefix('(')
        .and_then(|rest| rest.strip_suffix(')'))
        .ok_or_else(|| format!("Cannot parse vector: {trimmed:?}"))?;
    let values = inner
        .split_whitespace()
        .map(|token| parse_scalar(token, "vector component"))
        .collect::<Result<Vec<f64>, String>>()?;
    match values.as_slice() {
        [x, y, z] => Ok([*x, *y, *z]),
        _ => Err(format!(
            "A vector needs 3 components, got {}: {trimmed:?}",
            values.len()
        )),
    }
}

/// Interpret the text that follows the `internalField` keyword.
fn parse_entry(after_keyword: &str) -> Result<InternalField, String> {
    let (mode, rest) = split_word(after_keyword.trim_start());
    match mode {
        "nonuniform" => parse_nonuniform(rest),
        "uniform" => {
            let end = rest.find(';').ok_or_else(|| {
                "internalField: the uniform value is not closed by ';'".to_string()
            })?;
            let value = rest.get(..end).unwrap_or("").trim();
            if value.is_empty() {
                Err("internalField: uniform has no value".to_string())
            } else {
                Ok(InternalField::Uniform(value.to_string()))
            }
        }
        other => Err(format!(
            "internalField: expected 'uniform' or 'nonuniform', found {other:?}"
        )),
    }
}

/// Interpret `[List<type>] [N] ( ... )` or `[List<type>] N{ value }`.
fn parse_nonuniform(after_mode: &str) -> Result<InternalField, String> {
    let mut rest = after_mode.trim_start();
    // The type name, `List<vector>`, comes before the count.
    let (word, after_word) = split_word(rest);
    if word.starts_with(|c: char| c.is_ascii_alphabetic()) {
        rest = after_word.trim_start();
    }
    let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
    let declared = match rest.get(..digits).filter(|text| !text.is_empty()) {
        Some(text) => Some(
            text.parse::<usize>()
                .map_err(|_| format!("internalField: the count {text} is too large"))?,
        ),
        None => None,
    };
    let rest = rest.get(digits..).unwrap_or("").trim_start();
    match rest.chars().next() {
        Some('(') => {
            let close = matching_close(rest, '(', ')')
                .ok_or_else(|| "internalField: unmatched '(' in the list".to_string())?;
            Ok(InternalField::Listed {
                declared,
                body: rest.get(1..close).unwrap_or("").to_string(),
            })
        }
        Some('{') => {
            let count = declared.ok_or_else(|| {
                "internalField: a '{value}' list needs an element count".to_string()
            })?;
            let close = matching_close(rest, '{', '}')
                .ok_or_else(|| "internalField: unmatched '{' in the list".to_string())?;
            Ok(InternalField::Repeated {
                count,
                value: rest.get(1..close).unwrap_or("").trim().to_string(),
            })
        }
        _ => Err("internalField: expected the list '(' after 'nonuniform'".to_string()),
    }
}

/// Split off the leading word, which ends at white space, a parenthesis, a
/// brace or a semicolon.
fn split_word(text: &str) -> (&str, &str) {
    let end = text
        .find(|c: char| c.is_whitespace() || matches!(c, '(' | '{' | ';'))
        .unwrap_or(text.len());
    text.split_at(end)
}

/// Byte offset of the delimiter that closes the one `text` starts with.
fn matching_close(text: &str, open: char, close: char) -> Option<usize> {
    let mut depth = 0usize;
    let mut quoted = false;
    let mut escaped = false;
    for (offset, c) in text.char_indices() {
        if quoted {
            match c {
                _ if escaped => escaped = false,
                '\\' => escaped = true,
                '"' => quoted = false,
                _ => {}
            }
        } else if c == '"' {
            quoted = true;
        } else if c == open {
            depth += 1;
        } else if c == close {
            depth = depth.checked_sub(1)?;
            if depth == 0 {
                return Some(offset);
            }
        }
    }
    None
}

/// Byte offset just behind the top-level `keyword` entry, if the file has one.
///
/// Only a keyword at brace depth zero that stands at the start of a statement
/// counts; `$internalField` and `internalField` inside a patch dictionary do
/// not.
fn find_entry(text: &str, keyword: &str) -> Option<usize> {
    let bytes = text.as_bytes();
    let mut depth = 0usize;
    let mut quoted = false;
    let mut escaped = false;
    for (index, &byte) in bytes.iter().enumerate() {
        if quoted {
            match byte {
                _ if escaped => escaped = false,
                b'\\' => escaped = true,
                b'"' => quoted = false,
                _ => {}
            }
            continue;
        }
        match byte {
            b'"' => quoted = true,
            b'{' => depth += 1,
            b'}' => depth = depth.saturating_sub(1),
            _ if depth == 0 && keyword_at(bytes, index, keyword) => {
                return Some(index + keyword.len());
            }
            _ => {}
        }
    }
    None
}

fn keyword_at(bytes: &[u8], index: usize, keyword: &str) -> bool {
    let starts_statement = match index.checked_sub(1).and_then(|before| bytes.get(before)) {
        None => true,
        Some(previous) => previous.is_ascii_whitespace() || matches!(previous, b';' | b'}'),
    };
    starts_statement
        && bytes.get(index..index + keyword.len()) == Some(keyword.as_bytes())
        && bytes
            .get(index + keyword.len())
            .is_some_and(u8::is_ascii_whitespace)
}

/// `text` without its `//` and `/* */` comments.
///
/// A block comment becomes one blank so that the words around it stay apart;
/// the newline that ends a line comment is kept.  String literals are copied
/// unchanged.
fn strip_comments(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' => {
                out.push(c);
                let mut escaped = false;
                for inner in chars.by_ref() {
                    out.push(inner);
                    match inner {
                        _ if escaped => escaped = false,
                        '\\' => escaped = true,
                        '"' => break,
                        _ => {}
                    }
                }
            }
            '/' if chars.peek() == Some(&'/') => {
                for skipped in chars.by_ref() {
                    if skipped == '\n' {
                        out.push('\n');
                        break;
                    }
                }
            }
            '/' if chars.peek() == Some(&'*') => {
                chars.next();
                let mut previous = '\0';
                for skipped in chars.by_ref() {
                    if previous == '*' && skipped == '/' {
                        break;
                    }
                    previous = skipped;
                }
                out.push(' ');
            }
            _ => out.push(c),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::foam_io::{parse_foam_scalar_field_internal, parse_foam_vector_field};

    /// A `U` file as `simpleFoam` writes it: four cells, an inlet with a
    /// uniform vector value, a walls patch with a non-uniform one and an
    /// outlet whose `value` list is a different length from the internal one.
    const VELOCITY_FILE: &str = "\
/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org (see docs)
    \\\\  /    A nd           | Version:  10
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volVectorField;
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   nonuniform List<vector>
4
(
(0.5 0 0)
(1 -0.25 0)
(1.5e-1 0 2)
(-3 4 5)
)
;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform (1 0 0);
    }
    outlet
    {
        type            zeroGradient;
    }
    walls
    {
        type            fixedValue;
        value           nonuniform List<vector> 2((0 0 0) (7 8 9));
    }
}


// ************************************************************************* //
";

    #[test]
    fn velocity_file_returns_only_the_cell_values() {
        let cells = parse_foam_vector_field(VELOCITY_FILE).unwrap();
        assert_eq!(
            cells,
            vec![
                [0.5, 0.0, 0.0],
                [1.0, -0.25, 0.0],
                [0.15, 0.0, 2.0],
                [-3.0, 4.0, 5.0],
            ]
        );
    }

    #[test]
    fn uniform_vector_internal_field_is_a_single_value() {
        let file = "\
FoamFile { class volVectorField; object U; }
dimensions [0 1 -1 0 0 0 0];
internalField   uniform (0 0 0);
boundaryField
{
    inlet { type fixedValue; value uniform (1 0 0); }
    outlet { type zeroGradient; }
    walls { type fixedValue; value nonuniform List<vector> 2((0 0 0) (7 8 9)); }
}
";
        let values = parse_foam_vector_field(file).unwrap();
        assert_eq!(values, vec![[0.0, 0.0, 0.0]]);
    }

    #[test]
    fn scalar_file_ignores_boundary_value_lists() {
        let file = "\
FoamFile { class volScalarField; object p; }
dimensions [0 2 -2 0 0 0 0];
internalField   nonuniform List<scalar> 3(1 2.5 -3);
boundaryField
{
    inlet { type fixedValue; value nonuniform List<scalar> 2(5 6); }
    outlet { type zeroGradient; }
}
";
        assert_eq!(
            parse_foam_scalar_field_internal(file).unwrap(),
            vec![1.0, 2.5, -3.0]
        );
    }

    #[test]
    fn scalar_file_in_the_multi_line_layout() {
        let file = "\
internalField   nonuniform List<scalar>
4
(
10
20
30
40
)
;

boundaryField
{
    outlet { type fixedValue; value nonuniform List<scalar> 1(99); }
}
";
        assert_eq!(
            parse_foam_scalar_field_internal(file).unwrap(),
            vec![10.0, 20.0, 30.0, 40.0]
        );
    }

    #[test]
    fn uniform_scalar_internal_field_is_a_single_value() {
        let file = "internalField   uniform 101325;\nboundaryField { a { type fixedValue; value nonuniform List<scalar> 2(5 6); } }\n";
        assert_eq!(
            parse_foam_scalar_field_internal(file).unwrap(),
            vec![101_325.0]
        );
        let exponent = "internalField uniform 1.5e-3;\n";
        assert_eq!(
            parse_foam_scalar_field_internal(exponent).unwrap(),
            vec![1.5e-3]
        );
    }

    #[test]
    fn empty_lists_are_valid() {
        let scalars = "internalField nonuniform List<scalar> 0();\n";
        assert!(
            parse_foam_scalar_field_internal(scalars)
                .unwrap()
                .is_empty()
        );
        let vectors = "internalField nonuniform List<vector>\n0\n(\n)\n;\n";
        assert!(parse_foam_vector_field(vectors).unwrap().is_empty());
    }

    #[test]
    fn declared_count_must_match_the_values() {
        let too_few = "internalField nonuniform List<vector> 3((1 0 0) (2 0 0));\n";
        let message = parse_foam_vector_field(too_few).unwrap_err();
        assert!(
            message.contains("declares 3") && message.contains("holds 2"),
            "{message}"
        );
        let too_many = "internalField nonuniform List<scalar> 2(1 2 3);\n";
        let message = parse_foam_scalar_field_internal(too_many).unwrap_err();
        assert!(
            message.contains("declares 2") && message.contains("holds 3"),
            "{message}"
        );
    }

    #[test]
    fn malformed_values_are_reported_with_the_token() {
        let scalar = "internalField nonuniform List<scalar> 2(1 abc);\n";
        let message = parse_foam_scalar_field_internal(scalar).unwrap_err();
        assert!(message.contains("\"abc\""), "{message}");
        let short = "internalField nonuniform List<vector> 1((1 2));\n";
        let message = parse_foam_vector_field(short).unwrap_err();
        assert!(message.contains("3 components"), "{message}");
        let stray = "internalField nonuniform List<vector> 1((1 2 3) 4);\n";
        let message = parse_foam_vector_field(stray).unwrap_err();
        assert!(message.contains("Expected '(x y z)'"), "{message}");
        let bad_uniform = "internalField uniform (1 two 3);\n";
        assert!(parse_foam_vector_field(bad_uniform).is_err());
    }

    #[test]
    fn structural_errors_are_reported() {
        for (file, fragment) in [
            ("internalField uniform 5\n", "not closed by ';'"),
            ("internalField uniform ;\n", "no value"),
            (
                "internalField #calc \"1\";\n",
                "expected 'uniform' or 'nonuniform'",
            ),
            (
                "internalField nonuniform List<scalar> 3(1 2 3\n",
                "unmatched '('",
            ),
            (
                "internalField nonuniform List<scalar> 3 1 2 3;\n",
                "expected the list '('",
            ),
            (
                "internalField nonuniform List<scalar> {1};\n",
                "needs an element count",
            ),
        ] {
            let message = parse_foam_scalar_field_internal(file).unwrap_err();
            assert!(message.contains(fragment), "{file:?}: {message}");
        }
    }

    #[test]
    fn repeated_lists_expand_to_the_declared_count() {
        let scalars = "internalField nonuniform List<scalar> 4{2.5};\n";
        assert_eq!(
            parse_foam_scalar_field_internal(scalars).unwrap(),
            vec![2.5; 4]
        );
        let vectors = "internalField nonuniform List<vector> 3{(1 0 -1)};\n";
        assert_eq!(
            parse_foam_vector_field(vectors).unwrap(),
            vec![[1.0, 0.0, -1.0]; 3]
        );
        let huge = "internalField nonuniform List<scalar> 99999999999{0};\n";
        let message = parse_foam_scalar_field_internal(huge).unwrap_err();
        assert!(message.contains("limit"), "{message}");
    }

    #[test]
    fn comments_do_not_hide_or_fake_the_entry() {
        let file = "\
/* internalField uniform (9 9 9); ) ( */
// internalField uniform (8 8 8);
FoamFile { class volVectorField; }
internalField /* the cells */ nonuniform List<vector> 1( // one cell
    (1 2 3) // its velocity )
);
boundaryField { a { type fixedValue; value uniform (4 5 6); } }
";
        assert_eq!(
            parse_foam_vector_field(file).unwrap(),
            vec![[1.0, 2.0, 3.0]]
        );
    }

    #[test]
    fn references_to_internal_field_inside_patches_are_not_the_entry() {
        let file = "\
FoamFile { class volVectorField; object U; }
internalField   uniform (0 0 1);
boundaryField
{
    inlet { type fixedValue; value $internalField; }
    internalField { type calculated; }
}
";
        assert_eq!(
            parse_foam_vector_field(file).unwrap(),
            vec![[0.0, 0.0, 1.0]]
        );
        // Without a top-level entry the file is a bare list, not a field.
        let only_reference = "boundaryField { a { value $internalField; } }\n";
        assert!(internal_field_block(only_reference).unwrap().is_none());
    }

    #[test]
    fn bare_lists_are_still_read() {
        let points = "1\n(\n    (1.0 2.0 3.0)\n)";
        assert_eq!(
            parse_foam_vector_field(points).unwrap(),
            vec![[1.0, 2.0, 3.0]]
        );
        let scalars = "3\n(\n1\n2\n3\n)\n";
        assert_eq!(
            parse_foam_scalar_field_internal(scalars).unwrap(),
            vec![1.0, 2.0, 3.0]
        );
    }

    #[test]
    fn a_bare_list_ends_at_its_own_closing_parenthesis() {
        // Text after the list, even with parentheses, is not part of it.
        let file =
            "FoamFile { class vectorField; }\n2\n(\n(1 0 0)\n(0 1 0)\n)\n// end of file (really)\n";
        let vectors = parse_foam_vector_field(file).unwrap();
        assert_eq!(vectors, vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let wrong_count = "3\n(\n(1 0 0)\n(0 1 0)\n)\n";
        let message = parse_foam_vector_field(wrong_count).unwrap_err();
        assert!(message.contains("declares 3"), "{message}");
        assert!(parse_foam_vector_field("no list here").is_err());
        assert!(parse_foam_scalar_field_internal("2\n(\n1 2\n").is_err());
    }

    #[test]
    fn internal_field_block_describes_the_entry() {
        assert_eq!(
            internal_field_block("internalField uniform 4;\n").unwrap(),
            Some(InternalField::Uniform("4".to_string()))
        );
        assert_eq!(
            internal_field_block("internalField nonuniform List<scalar> 2(1 2);\n").unwrap(),
            Some(InternalField::Listed {
                declared: Some(2),
                body: "1 2".to_string()
            })
        );
        assert_eq!(
            internal_field_block("internalField nonuniform 2{7};\n").unwrap(),
            Some(InternalField::Repeated {
                count: 2,
                value: "7".to_string()
            })
        );
        assert_eq!(internal_field_block("").unwrap(), None);
    }

    #[test]
    fn written_fields_read_back() {
        let values = vec![[1.0, 0.0, 0.0], [0.5, 0.5, 0.0], [-2.0, 1.0e-3, 4.0]];
        let patches: Vec<(&str, crate::foam_io::FoamBcType)> = vec![
            (
                "inlet",
                crate::foam_io::FoamBcType::FixedVectorValue([9.0, 8.0, 7.0]),
            ),
            ("outlet", crate::foam_io::FoamBcType::ZeroGradient),
        ];
        let text = crate::foam_io::write_vector_field("U", "[0 1 -1 0 0 0 0]", &values, &patches);
        assert_eq!(parse_foam_vector_field(&text).unwrap(), values);
        let scalar_text = crate::foam_io::FoamScalarField {
            name: "p".to_string(),
            internal_field: vec![3.0, 4.0, 5.5],
            boundary_values: vec![vec![1.0, 2.0], Vec::new()],
        }
        .to_foam_string();
        assert_eq!(
            parse_foam_scalar_field_internal(&scalar_text).unwrap(),
            vec![3.0, 4.0, 5.5]
        );
    }
}
