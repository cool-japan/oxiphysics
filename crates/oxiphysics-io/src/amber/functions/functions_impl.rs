//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::super::types::{
    AmberAngle, AmberBond, AmberDihedral, AmberEnergyFrame, AmberTopology, FfAngleParam,
    FfBondParam, LjParameter, MdcrdFrame,
};
use crate::fortran_fields::parse_real_record;

/// Extract the raw text body of a `%FLAG `name` section.
///
/// Returns the lines after `%FORMAT(...)` up to the next `%FLAG` or EOF.
pub fn parse_flag_section(content: &str, flag_name: &str) -> Option<String> {
    let target = format!("%FLAG {}", flag_name.to_uppercase());
    let lines: Vec<&str> = content.lines().collect();
    let flag_pos = lines.iter().position(|l| l.trim() == target.trim())?;
    let data_start = flag_pos + 2;
    let mut body = String::new();
    for line in lines.iter().skip(data_start) {
        let trimmed = line.trim();
        if trimmed.starts_with("%FLAG") || trimmed.starts_with("%VERSION") {
            break;
        }
        body.push_str(line);
        body.push('\n');
    }
    Some(body)
}
/// Parse Fortran-formatted real numbers (fixed-width or space-separated).
pub fn parse_fortran_reals(s: &str) -> Vec<f64> {
    s.split_whitespace()
        .filter_map(|tok| tok.replace('D', "E").replace('d', "e").parse::<f64>().ok())
        .collect()
}
/// Parse Fortran-formatted integers (space-separated).
pub fn parse_fortran_ints(s: &str) -> Vec<i64> {
    s.split_whitespace()
        .filter_map(|tok| tok.parse::<i64>().ok())
        .collect()
}
/// Parse Fortran-formatted string tokens (fixed-width 4-char columns or whitespace separated).
pub fn parse_fortran_strings(s: &str) -> Vec<String> {
    let mut result = Vec::new();
    for line in s.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let chars: Vec<char> = line.chars().collect();
        if chars.len() >= 4 {
            let mut i = 0;
            while i + 4 <= chars.len() {
                let tok: String = chars[i..i + 4].iter().collect();
                let trimmed = tok.trim().to_string();
                if !trimmed.is_empty() {
                    result.push(trimmed);
                }
                i += 4;
            }
        } else {
            for tok in line.split_whitespace() {
                result.push(tok.to_string());
            }
        }
    }
    result
}
/// Zero-based atom index designated by a prmtop atom pointer: `|pointer| / 3`.
///
/// The sign of a pointer carries no atom information (AMBER only uses it as a
/// flag in dihedral records), so it is ignored.  The result saturates at
/// `usize::MAX`, which no topology contains, instead of wrapping around to an
/// index that could designate a real atom.
fn atom_index_of_pointer(pointer: i64) -> usize {
    usize::try_from(pointer.unsigned_abs() / 3).unwrap_or(usize::MAX)
}
/// Zero-based parameter type designated by a one-based prmtop type index.
///
/// A non-positive index selects the first type; an index too large for
/// `usize` saturates, so that the table lookup of the caller misses.
fn type_slot_of_index(type_index: i64) -> usize {
    usize::try_from(type_index.saturating_sub(1).max(0)).unwrap_or(usize::MAX)
}
/// Build `AmberBond` entries from a flat integer array
/// (`3*i`, `3*j`, `type_idx` triplets).
///
/// The first two integers of each triplet are atom pointers (three times the
/// zero-based atom index, sign ignored), the third is the one-based index into
/// the bond parameter tables.  A trailing incomplete triplet is ignored.
/// [`AmberTopology::from_prmtop_str`] validates the pointers and the type
/// indices before calling this function; for unchecked input a non-positive
/// type index selects the first parameter type and an index beyond the tables
/// yields `k = 0` and `r0 = 0`.
pub fn build_bonds(ints: &[i64], force_k: &[f64], equil_r: &[f64]) -> Vec<AmberBond> {
    let mut bonds = Vec::new();
    let mut idx = 0;
    while idx + 2 < ints.len() {
        let i = atom_index_of_pointer(ints[idx]);
        let j = atom_index_of_pointer(ints[idx + 1]);
        let type_idx = type_slot_of_index(ints[idx + 2]);
        let k = force_k.get(type_idx).copied().unwrap_or(0.0);
        let r0 = equil_r.get(type_idx).copied().unwrap_or(0.0);
        bonds.push(AmberBond { i, j, k, r0 });
        idx += 3;
    }
    bonds
}
/// Build `AmberAngle` entries from a flat integer array
/// (`3*i`, `3*j`, `3*k`, `type_idx` quadruples).
///
/// The first three integers of each quadruple are atom pointers (three times
/// the zero-based atom index, sign ignored), the fourth is the one-based index
/// into the angle parameter tables.  A trailing incomplete quadruple is
/// ignored.  [`AmberTopology::from_prmtop_str`] validates the pointers and the
/// type indices before calling this function; for unchecked input a
/// non-positive type index selects the first parameter type and an index
/// beyond the tables yields `k = 0` and `theta0 = 0`.
pub fn build_angles(ints: &[i64], force_k: &[f64], equil_t: &[f64]) -> Vec<AmberAngle> {
    let mut angles = Vec::new();
    let mut idx = 0;
    while idx + 3 < ints.len() {
        let i = atom_index_of_pointer(ints[idx]);
        let j = atom_index_of_pointer(ints[idx + 1]);
        let atom_k = atom_index_of_pointer(ints[idx + 2]);
        let k_idx = type_slot_of_index(ints[idx + 3]);
        let k = force_k.get(k_idx).copied().unwrap_or(0.0);
        let theta0 = equil_t.get(k_idx).copied().unwrap_or(0.0);
        angles.push(AmberAngle {
            i,
            j,
            atom_k,
            k_idx,
            k,
            theta0,
        });
        idx += 4;
    }
    angles
}
/// Decode a prmtop atom pointer (three times the zero-based atom index).
///
/// Returns `None` for a pointer that cannot designate an atom: a negative
/// value or one that is not a multiple of three.
pub(crate) fn decode_atom_pointer(pointer: i64) -> Option<usize> {
    let pointer = u64::try_from(pointer).ok()?;
    if !pointer.is_multiple_of(3) {
        return None;
    }
    usize::try_from(pointer / 3).ok()
}
/// Read a non-negative counter from the `POINTERS` array (zero when the array
/// is too short to contain it).
pub(crate) fn pointer_count(pointers: &[i64], index: usize, name: &str) -> Result<usize, String> {
    match pointers.get(index) {
        None => Ok(0),
        Some(&value) => {
            usize::try_from(value).map_err(|_| format!("POINTERS: {name} is negative ({value})"))
        }
    }
}
/// Check the records of a flat `BONDS_*` or `ANGLES_*` prmtop section.
///
/// Every record holds `atoms` atom pointers followed by a one-based parameter
/// type index.  Returns a description of the first defect: a trailing
/// incomplete record, an atom pointer that does not designate one of the
/// `n_atoms` atoms, or a type index outside `1..=n_types`.
pub(crate) fn validate_valence_records(
    flag: &str,
    ints: &[i64],
    atoms: usize,
    n_atoms: usize,
    n_types: usize,
) -> Result<(), String> {
    let stride = atoms + 1;
    if !ints.len().is_multiple_of(stride) {
        return Err(format!(
            "{flag}: {} integers do not form whole records of {stride} values",
            ints.len()
        ));
    }
    for (record_idx, record) in ints.chunks_exact(stride).enumerate() {
        for &pointer in &record[..atoms] {
            let designates_atom = decode_atom_pointer(pointer).is_some_and(|atom| atom < n_atoms);
            if !designates_atom {
                return Err(format!(
                    "{flag}: record {record_idx} has atom pointer {pointer}, which does not designate one of the {n_atoms} atoms declared in POINTERS"
                ));
            }
        }
        let type_index = record[atoms];
        let known_type = usize::try_from(type_index).is_ok_and(|t| (1..=n_types).contains(&t));
        if !known_type {
            return Err(format!(
                "{flag}: record {record_idx} has parameter type index {type_index}, expected 1..={n_types}"
            ));
        }
    }
    Ok(())
}
/// Parse the Lennard-Jones pair parameters of a raw prmtop string.
///
/// `LENNARD_JONES_ACOEF` and `LENNARD_JONES_BCOEF` hold one entry per
/// unordered pair of atom types (`NTYPES * (NTYPES + 1) / 2` entries).  The
/// `NONBONDED_PARM_INDEX` section (`NTYPES` x `NTYPES` entries, row-major,
/// with `NTYPES` taken from `POINTERS`) maps a pair of types to the one-based
/// position of its entry.  A non-positive entry designates a 10-12 hydrogen
/// bond term that has no 12-6 parameters and is skipped.
///
/// One [`LjParameter`] is returned for every pair `type_a <= type_b` that has
/// a 12-6 term, with zero-based type indices.  The pair assignment is never
/// guessed: if `POINTERS` or `NONBONDED_PARM_INDEX` is missing, or an index
/// points outside the coefficient tables, the result is empty.
pub fn extract_lj_parameters(prmtop: &str) -> Vec<LjParameter> {
    let pointers = extract_flag_integers(prmtop, "POINTERS");
    let Some(n_types) = pointers
        .get(1)
        .and_then(|&n| usize::try_from(n).ok())
        .filter(|&n| n > 0)
    else {
        return Vec::new();
    };
    let nb_index = extract_flag_integers(prmtop, "NONBONDED_PARM_INDEX");
    let a_vals = extract_flag_floats(prmtop, "LENNARD_JONES_ACOEF");
    let b_vals = extract_flag_floats(prmtop, "LENNARD_JONES_BCOEF");
    let mut out = Vec::new();
    for type_a in 0..n_types {
        for type_b in type_a..n_types {
            let Some(&entry) = n_types
                .checked_mul(type_a)
                .and_then(|row| row.checked_add(type_b))
                .and_then(|slot| nb_index.get(slot))
            else {
                return Vec::new();
            };
            if entry <= 0 {
                continue;
            }
            let position = usize::try_from(entry - 1).ok();
            let (Some(&a_coeff), Some(&b_coeff)) = (
                position.and_then(|p| a_vals.get(p)),
                position.and_then(|p| b_vals.get(p)),
            ) else {
                return Vec::new();
            };
            out.push(LjParameter {
                type_a,
                type_b,
                a_coeff,
                b_coeff,
            });
        }
    }
    out
}
/// Parse dihedral entries from the `DIHEDRALS_INC_HYDROGEN` and
/// `DIHEDRALS_WITHOUT_HYDROGEN` sections.
///
/// Every entry is five integers: three times the zero-based indices of the
/// four atoms, then the one-based index of the torsion parameter set in the
/// `DIHEDRAL_FORCE_CONSTANT`, `DIHEDRAL_PERIODICITY` and `DIHEDRAL_PHASE`
/// tables.  The sign of the third atom pointer marks a term whose 1-4
/// non-bonded interaction is skipped ([`AmberDihedral::skip_14`]); the sign of
/// the fourth marks an improper torsion ([`AmberDihedral::is_improper`]).
///
/// AMBER stores `PK` with `E = PK * (1 + cos(PN * phi - PHASE))`, whereas the
/// returned [`AmberDihedral::v_n`] follows the `V_n / 2 * (1 + cos(..))`
/// convention of `dihedral_energy`, so `v_n` is `2 * PK`.
///
/// Entries that cannot be decoded are skipped: an atom pointer that is not a
/// multiple of three, an atom beyond the `NATOM` declared in `POINTERS` (when
/// present), a parameter index outside the three tables, or a parameter set
/// that is not finite.  Terms are returned in file order, hydrogen terms
/// first.
pub fn extract_dihedrals(prmtop: &str) -> Vec<AmberDihedral> {
    let force_constants = extract_flag_floats(prmtop, "DIHEDRAL_FORCE_CONSTANT");
    let periodicities = extract_flag_floats(prmtop, "DIHEDRAL_PERIODICITY");
    let phases = extract_flag_floats(prmtop, "DIHEDRAL_PHASE");
    let n_atoms = extract_flag_integers(prmtop, "POINTERS")
        .first()
        .and_then(|&n| usize::try_from(n).ok());
    // The sign of a dihedral atom pointer carries a flag; the magnitude
    // designates the atom.
    let atom = |pointer: i32| -> Option<usize> {
        let atom = decode_atom_pointer(i64::from(pointer).abs())?;
        n_atoms.is_none_or(|n| atom < n).then_some(atom)
    };
    let mut out = Vec::new();
    for flag in ["DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN"] {
        for chunk in extract_flag_integers(prmtop, flag).chunks_exact(5) {
            let (Some(i), Some(j), Some(k), Some(l)) = (
                atom(chunk[0]),
                atom(chunk[1]),
                atom(chunk[2]),
                atom(chunk[3]),
            ) else {
                continue;
            };
            let Some(term) = usize::try_from(chunk[4])
                .ok()
                .and_then(|t| t.checked_sub(1))
            else {
                continue;
            };
            let (Some(&pk), Some(&pn), Some(&phase)) = (
                force_constants.get(term),
                periodicities.get(term),
                phases.get(term),
            ) else {
                continue;
            };
            let periodicity = pn.abs().round();
            let v_n = 2.0 * pk;
            if !(v_n.is_finite() && phase.is_finite() && periodicity <= f64::from(i32::MAX)) {
                continue;
            }
            out.push(AmberDihedral {
                i,
                j,
                k,
                l,
                v_n,
                gamma: phase,
                n: periodicity as i32,
                is_improper: chunk[3] < 0,
                skip_14: chunk[2] < 0,
            });
        }
    }
    out
}
/// Parse a sequence of energy frames from an AMBER mdout-style string.
///
/// Only lines containing `NSTEP =` and `Etot` are considered; the rest are
/// skipped.  Returns one frame per `NSTEP` block.
pub fn parse_energy_frames(mdout: &str) -> Vec<AmberEnergyFrame> {
    let mut frames = Vec::new();
    let mut current_step: Option<u64> = None;
    let mut current_time: f64 = 0.0;
    let mut current_temp: f64 = 0.0;
    let mut e_tot = 0.0f64;
    let mut e_kin = 0.0f64;
    let mut e_pot = 0.0f64;
    let mut have_energy = false;
    for line in mdout.lines() {
        let t = line.trim();
        if t.contains("NSTEP =") {
            if let Some(step) = current_step.take() {
                if have_energy {
                    frames.push(AmberEnergyFrame {
                        step,
                        time_ps: current_time,
                        temp_k: current_temp,
                        e_tot,
                        e_kin,
                        e_pot,
                    });
                }
                have_energy = false;
            }
            current_step = parse_field(t, "NSTEP =").and_then(|s| s.parse().ok());
            current_time = parse_field(t, "TIME(PS) =")
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0);
            current_temp = parse_field(t, "TEMP(K) =")
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0);
        } else if t.starts_with("Etot") {
            e_tot = parse_field(t, "Etot   =")
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0);
            e_kin = parse_field(t, "EKtot   =")
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0);
            e_pot = parse_field(t, "EPtot   =")
                .and_then(|s| s.parse().ok())
                .unwrap_or(0.0);
            have_energy = true;
        }
    }
    if let Some(step) = current_step
        && have_energy
    {
        frames.push(AmberEnergyFrame {
            step,
            time_ps: current_time,
            temp_k: current_temp,
            e_tot,
            e_kin,
            e_pot,
        });
    }
    frames
}
/// Write a minimal AMBER restart (.inpcrd) string from positions and optional
/// velocities.
///
/// Each coordinate triplet is formatted to 12.7f in groups of 2 per line
/// (6 values per line), which matches the standard AMBER inpcrd layout.
pub fn write_inpcrd(
    title: &str,
    positions: &[[f64; 3]],
    velocities: Option<&[[f64; 3]]>,
) -> String {
    let mut out = String::new();
    out.push_str(title);
    out.push('\n');
    out.push_str(&format!("{:6}\n", positions.len()));
    let flat_pos: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    for chunk in flat_pos.chunks(6) {
        let line: String = chunk.iter().map(|v| format!("{:12.7}", v)).collect();
        out.push_str(&line);
        out.push('\n');
    }
    if let Some(vels) = velocities {
        let flat_vel: Vec<f64> = vels.iter().flat_map(|v| v.iter().copied()).collect();
        for chunk in flat_vel.chunks(6) {
            let line: String = chunk.iter().map(|v| format!("{:12.7}", v)).collect();
            out.push_str(&line);
            out.push('\n');
        }
    }
    out
}
/// Minimal built-in AMBER99SB bond parameter table.
///
/// Only the most common backbone bonds are included for testing/demo purposes.
pub fn amber99sb_bonds() -> Vec<FfBondParam> {
    vec![
        FfBondParam {
            type_a: "N".into(),
            type_b: "H".into(),
            k: 434.0,
            r0: 1.010,
        },
        FfBondParam {
            type_a: "N".into(),
            type_b: "CT".into(),
            k: 337.0,
            r0: 1.449,
        },
        FfBondParam {
            type_a: "CT".into(),
            type_b: "CT".into(),
            k: 310.0,
            r0: 1.526,
        },
        FfBondParam {
            type_a: "CT".into(),
            type_b: "HC".into(),
            k: 340.0,
            r0: 1.090,
        },
        FfBondParam {
            type_a: "CT".into(),
            type_b: "C".into(),
            k: 317.0,
            r0: 1.522,
        },
        FfBondParam {
            type_a: "C".into(),
            type_b: "O".into(),
            k: 570.0,
            r0: 1.229,
        },
        FfBondParam {
            type_a: "C".into(),
            type_b: "N".into(),
            k: 490.0,
            r0: 1.335,
        },
        FfBondParam {
            type_a: "OH".into(),
            type_b: "HO".into(),
            k: 553.0,
            r0: 0.960,
        },
    ]
}
/// Minimal built-in AMBER99SB angle parameter table.
pub fn amber99sb_angles() -> Vec<FfAngleParam> {
    vec![
        FfAngleParam {
            type_i: "H".into(),
            type_j: "N".into(),
            type_k: "H".into(),
            k: 35.0,
            theta0_deg: 120.0,
        },
        FfAngleParam {
            type_i: "H".into(),
            type_j: "N".into(),
            type_k: "CT".into(),
            k: 50.0,
            theta0_deg: 118.0,
        },
        FfAngleParam {
            type_i: "N".into(),
            type_j: "CT".into(),
            type_k: "C".into(),
            k: 63.0,
            theta0_deg: 111.2,
        },
        FfAngleParam {
            type_i: "N".into(),
            type_j: "CT".into(),
            type_k: "CT".into(),
            k: 80.0,
            theta0_deg: 111.2,
        },
        FfAngleParam {
            type_i: "CT".into(),
            type_j: "C".into(),
            type_k: "O".into(),
            k: 80.0,
            theta0_deg: 120.4,
        },
        FfAngleParam {
            type_i: "CT".into(),
            type_j: "C".into(),
            type_k: "N".into(),
            k: 70.0,
            theta0_deg: 116.6,
        },
        FfAngleParam {
            type_i: "CT".into(),
            type_j: "CT".into(),
            type_k: "HC".into(),
            k: 50.0,
            theta0_deg: 109.5,
        },
        FfAngleParam {
            type_i: "HC".into(),
            type_j: "CT".into(),
            type_k: "HC".into(),
            k: 35.0,
            theta0_deg: 109.5,
        },
    ]
}
/// Compute summary statistics over a slice of energy frames.
///
/// Returns `(mean_etot, min_etot, max_etot)`.
pub fn energy_summary(frames: &[AmberEnergyFrame]) -> (f64, f64, f64) {
    if frames.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let vals: Vec<f64> = frames.iter().map(|f| f.e_tot).collect();
    let mean = vals.iter().sum::<f64>() / vals.len() as f64;
    let min = vals.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (mean, min, max)
}
pub(super) fn parse_field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let pos = line.find(key)?;
    let rest = line[pos + key.len()..].trim_start();
    let end = rest.find([' ', '\t', '\n']).unwrap_or(rest.len());
    Some(&rest[..end])
}
pub(super) fn extract_flag_floats(prmtop: &str, flag: &str) -> Vec<f64> {
    let marker = format!("%FLAG {flag}");
    let mut in_section = false;
    let mut out = Vec::new();
    for line in prmtop.lines() {
        let t = line.trim();
        if t == marker {
            in_section = true;
            continue;
        }
        if in_section {
            if t.starts_with("%FORMAT") {
                continue;
            }
            if t.starts_with('%') {
                break;
            }
            for tok in t.split_whitespace() {
                if let Ok(v) = tok.parse::<f64>() {
                    out.push(v);
                }
            }
        }
    }
    out
}
pub(super) fn extract_flag_integers(prmtop: &str, flag: &str) -> Vec<i32> {
    let marker = format!("%FLAG {flag}");
    let mut in_section = false;
    let mut out = Vec::new();
    for line in prmtop.lines() {
        let t = line.trim();
        if t == marker {
            in_section = true;
            continue;
        }
        if in_section {
            if t.starts_with("%FORMAT") {
                continue;
            }
            if t.starts_with('%') {
                break;
            }
            for tok in t.split_whitespace() {
                if let Ok(v) = tok.parse::<i32>() {
                    out.push(v);
                }
            }
        }
    }
    out
}
/// Strip AMBER-style inline comment (everything from `;` onward).
pub fn strip_amber_comment(s: &str) -> &str {
    if let Some(pos) = s.find(';') {
        s[..pos].trim()
    } else {
        s.trim()
    }
}
/// Collect all floating-point numbers from a string, ignoring non-numeric tokens.
pub fn collect_numbers(s: &str) -> Vec<f64> {
    s.split_whitespace()
        .filter_map(|tok| tok.parse::<f64>().ok())
        .collect()
}
/// Parse a residue selection like "1", "1-10", "1,3,5", "1-5,8,10-12".
pub fn parse_residue_selection(s: &str) -> Result<Vec<usize>, String> {
    let mut result = Vec::new();
    for token in s.split(',') {
        let t = token.trim();
        if t.is_empty() {
            continue;
        }
        if let Some(dash_pos) = t.find('-') {
            let start: usize = t[..dash_pos]
                .trim()
                .parse()
                .map_err(|_| format!("Invalid range start: '{}'", &t[..dash_pos]))?;
            let end: usize = t[dash_pos + 1..]
                .trim()
                .parse()
                .map_err(|_| format!("Invalid range end: '{}'", &t[dash_pos + 1..]))?;
            if end < start {
                return Err(format!("Range end {} < start {}", end, start));
            }
            result.extend(start..=end);
        } else {
            let n: usize = t
                .parse()
                .map_err(|_| format!("Invalid residue number: '{}'", t))?;
            result.push(n);
        }
    }
    Ok(result)
}
/// Number of coordinate values per line of an MDCRD record (`10F8.3`).
const MDCRD_VALUES_PER_LINE: usize = 10;
/// Column width of an MDCRD field (`F8.3`).
const MDCRD_COLUMN_WIDTH: usize = 8;
/// Parse the periodic-box line of an MDCRD frame (`a b c`, extra values are
/// ignored).
fn parse_mdcrd_box(line: &str) -> Option<[f64; 3]> {
    match parse_real_record(line, MDCRD_COLUMN_WIDTH).ok()?.as_slice() {
        [a, b, c, ..] => Some([*a, *b, *c]),
        _ => None,
    }
}
/// Parse AMBER MDCRD trajectory data from a string.
///
/// MDCRD format: a title line, then for every frame the coordinates in
/// `10F8.3` records (ten values per line, fewer on the last line of a frame),
/// followed by one `3F8.3` box line when `has_box` is set.
///
/// Fields are cut into 8-character columns when they touch each other
/// (`1.000-100.000`) and split on blanks otherwise.  Only complete frames are
/// returned: parsing stops at the first frame whose coordinate lines are
/// truncated, do not hold exactly `3 * n_atoms` numbers, or contain something
/// that is not a number.  A final frame without its box line, or with a box
/// line of fewer than three numbers, is returned with `box_dims == None`.
/// `n_atoms == 0` yields no frames.
pub fn parse_mdcrd(s: &str, n_atoms: usize, has_box: bool) -> Vec<MdcrdFrame> {
    let Some(n_coords) = n_atoms.checked_mul(3).filter(|&n| n > 0) else {
        return Vec::new();
    };
    let coords_lines = n_coords.div_ceil(MDCRD_VALUES_PER_LINE);
    let frame_lines = coords_lines + usize::from(has_box);
    let lines: Vec<&str> = s.lines().skip(1).collect();
    let mut frames = Vec::new();
    let mut first = 0usize;
    while first + coords_lines <= lines.len() {
        let box_line = first + coords_lines;
        let mut coords: Vec<f64> = Vec::with_capacity(n_coords);
        for line in &lines[first..box_line] {
            match parse_real_record(line, MDCRD_COLUMN_WIDTH) {
                Ok(values) => coords.extend(values),
                Err(_) => return frames,
            }
        }
        if coords.len() != n_coords {
            return frames;
        }
        let box_dims = if has_box {
            lines.get(box_line).and_then(|line| parse_mdcrd_box(line))
        } else {
            None
        };
        frames.push(MdcrdFrame { coords, box_dims });
        first += frame_lines;
    }
    frames
}
/// Serialises a minimal AMBER prmtop-format string from an `AmberTopology`.
///
/// The output contains `%VERSION`, `%FLAG TITLE`, `%FLAG POINTERS`,
/// `%FLAG ATOM_NAME`, `%FLAG CHARGE`, and `%FLAG MASS` sections so that it
/// can be round-tripped through `AmberTopology::from_prmtop_str`.
pub fn write_prmtop(topo: &AmberTopology) -> String {
    let mut out = String::new();
    out.push_str("%VERSION  VERSION_STAMP = V0001.000  DATE = 01/01/26  00:00:00\n");
    out.push_str("%FLAG TITLE\n");
    out.push_str("%FORMAT(20a4)\n");
    out.push_str(&format!("{}\n", topo.title));
    out.push_str("%FLAG POINTERS\n");
    out.push_str("%FORMAT(10I8)\n");
    let n_atoms = topo.atoms.len() as i64;
    let n_types = 1i64;
    let n_bonds = topo.bonds.len() as i64;
    let n_bonds_nh = n_bonds;
    let n_angles = topo.angles.len() as i64;
    let pointers: Vec<i64> = vec![
        n_atoms, n_types, n_bonds, n_bonds_nh, n_angles, n_angles, 0, 0, 0, 0,
    ];
    for (i, p) in pointers.iter().enumerate() {
        out.push_str(&format!("{:8}", p));
        if (i + 1) % 10 == 0 {
            out.push('\n');
        }
    }
    if !pointers.len().is_multiple_of(10) {
        out.push('\n');
    }
    out.push_str("%FLAG ATOM_NAME\n");
    out.push_str("%FORMAT(20a4)\n");
    for (i, atom) in topo.atoms.iter().enumerate() {
        out.push_str(&format!("{:<4}", &atom.name[..atom.name.len().min(4)]));
        if (i + 1) % 20 == 0 {
            out.push('\n');
        }
    }
    if !topo.atoms.is_empty() && !topo.atoms.len().is_multiple_of(20) {
        out.push('\n');
    }
    pub(super) const AMBER_CHARGE_SCALE: f64 = 18.2223;
    out.push_str("%FLAG CHARGE\n");
    out.push_str("%FORMAT(5E16.8)\n");
    for (i, atom) in topo.atoms.iter().enumerate() {
        out.push_str(&format!("{:16.8E}", atom.charge * AMBER_CHARGE_SCALE));
        if (i + 1) % 5 == 0 {
            out.push('\n');
        }
    }
    if !topo.atoms.is_empty() && !topo.atoms.len().is_multiple_of(5) {
        out.push('\n');
    }
    out.push_str("%FLAG MASS\n");
    out.push_str("%FORMAT(5E16.8)\n");
    for (i, atom) in topo.atoms.iter().enumerate() {
        out.push_str(&format!("{:16.8E}", atom.mass));
        if (i + 1) % 5 == 0 {
            out.push('\n');
        }
    }
    if !topo.atoms.is_empty() && !topo.atoms.len().is_multiple_of(5) {
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amber::types::*;
    use std::str::FromStr;
    fn minimal_prmtop() -> String {
        let n_charge = -0.4157 * 18.2223;
        let h_charge = 0.2719 * 18.2223;
        format!(
            "%VERSION  VERSION_STAMP = V0001.000\n\
             %FLAG TITLE\n\
             %FORMAT(20a4)\n\
             test molecule\n\
             %FLAG POINTERS\n\
             %FORMAT(10I8)\n\
               2    0    1    0    1    0    0    0    0    0\n\
               0    0    0    0    0    0    0    0    0    0\n\
               0    0    0    0    0    0    0    0    0    0\n\
               0\n\
             %FLAG ATOM_NAME\n\
             %FORMAT(20a4)\n\
             N   H   \n\
             %FLAG CHARGE\n\
             %FORMAT(5E16.8)\n\
               {:.8E}  {:.8E}\n\
             %FLAG MASS\n\
             %FORMAT(5E16.8)\n\
                1.40100000E+01   1.00800000E+00\n\
             %FLAG AMBER_ATOM_TYPE\n\
             %FORMAT(20a4)\n\
             N   H   \n\
             %FLAG RESIDUE_LABEL\n\
             %FORMAT(20a4)\n\
             ALA \n\
             %FLAG RESIDUE_POINTER\n\
             %FORMAT(10I8)\n\
               1\n\
             %FLAG BONDS_INC_HYDROGEN\n\
             %FORMAT(10I8)\n\
               0    3    1\n\
             %FLAG BONDS_WITHOUT_HYDROGEN\n\
             %FORMAT(10I8)\n\
             %FLAG BOND_FORCE_CONSTANT\n\
             %FORMAT(5E16.8)\n\
                4.34000000E+02\n\
             %FLAG BOND_EQUIL_VALUE\n\
             %FORMAT(5E16.8)\n\
                1.01000000E+00\n\
             %FLAG ANGLES_INC_HYDROGEN\n\
             %FORMAT(10I8)\n\
             %FLAG ANGLES_WITHOUT_HYDROGEN\n\
             %FORMAT(10I8)\n\
             %FLAG ANGLE_FORCE_CONSTANT\n\
             %FORMAT(5E16.8)\n\
             %FLAG ANGLE_EQUIL_VALUE\n\
             %FORMAT(5E16.8)\n",
            n_charge, h_charge
        )
    }
    fn parse_minimal() -> AmberTopology {
        AmberTopology::from_prmtop_str(&minimal_prmtop()).unwrap()
    }
    #[test]
    fn test_parse_flag_section_found() {
        let s = "%FLAG TITLE\n%FORMAT(20a4)\nhello\n%FLAG OTHER\n%FORMAT\ndata\n";
        let result = parse_flag_section(s, "TITLE");
        assert!(result.is_some());
        assert!(result.unwrap().contains("hello"));
    }
    #[test]
    fn test_parse_flag_section_not_found() {
        let s = "%FLAG TITLE\n%FORMAT(20a4)\nhello\n";
        let result = parse_flag_section(s, "MISSING");
        assert!(result.is_none());
    }
    #[test]
    fn test_parse_fortran_reals_basic() {
        let s = "   1.23000000E+00  -4.56000000E+00\n";
        let vals = parse_fortran_reals(s);
        assert_eq!(vals.len(), 2);
        assert!((vals[0] - 1.23).abs() < 1e-6);
        assert!((vals[1] + 4.56).abs() < 1e-6);
    }
    #[test]
    fn test_parse_fortran_reals_d_exponent() {
        let s = "1.0D+00  2.5D-01\n";
        let vals = parse_fortran_reals(s);
        assert_eq!(vals.len(), 2);
        assert!((vals[1] - 0.25).abs() < 1e-9);
    }
    #[test]
    fn test_parse_fortran_ints_basic() {
        let s = "   2    0    1    0\n";
        let vals = parse_fortran_ints(s);
        assert_eq!(vals, vec![2, 0, 1, 0]);
    }
    #[test]
    fn test_parse_title() {
        let top = parse_minimal();
        assert!(top.title.contains("test molecule"), "title='{}'", top.title);
    }
    #[test]
    fn test_n_atoms() {
        let top = parse_minimal();
        assert_eq!(top.n_atoms, 2);
    }
    #[test]
    fn test_atom_names() {
        let top = parse_minimal();
        assert_eq!(top.atoms[0].name, "N");
        assert_eq!(top.atoms[1].name, "H");
    }
    #[test]
    fn test_atom_masses() {
        let top = parse_minimal();
        assert!(
            (top.atoms[0].mass - 14.01).abs() < 0.01,
            "N mass={}",
            top.atoms[0].mass
        );
        assert!(
            (top.atoms[1].mass - 1.008).abs() < 0.01,
            "H mass={}",
            top.atoms[1].mass
        );
    }
    #[test]
    fn test_atom_charges_scaled() {
        let top = parse_minimal();
        assert!(
            (top.atoms[0].charge - (-0.4157)).abs() < 0.001,
            "N charge={:.4}",
            top.atoms[0].charge
        );
    }
    #[test]
    fn test_atom_types() {
        let top = parse_minimal();
        assert_eq!(top.atoms[0].atom_type, "N");
        assert_eq!(top.atoms[1].atom_type, "H");
    }
    #[test]
    fn test_residue_names() {
        let top = parse_minimal();
        let res = top.residue_names();
        assert!(res.contains(&"ALA".to_string()), "residues={res:?}");
    }
    #[test]
    fn test_bonds_parsed() {
        let top = parse_minimal();
        assert_eq!(top.bonds.len(), 1, "bonds={}", top.bonds.len());
        assert_eq!(top.bonds[0].i, 0);
        assert_eq!(top.bonds[0].j, 1);
    }
    #[test]
    fn test_bond_force_constant() {
        let top = parse_minimal();
        assert!((top.bonds[0].k - 434.0).abs() < 1.0, "k={}", top.bonds[0].k);
    }
    #[test]
    fn test_bond_equil_value() {
        let top = parse_minimal();
        assert!(
            (top.bonds[0].r0 - 1.01).abs() < 0.01,
            "r0={}",
            top.bonds[0].r0
        );
    }
    #[test]
    fn test_total_charge() {
        let top = parse_minimal();
        let q = top.total_charge();
        assert!((q - (-0.1438)).abs() < 0.01, "total_charge={q:.4}");
    }
    #[test]
    fn test_total_mass() {
        let top = parse_minimal();
        let m = top.total_mass();
        assert!((m - 15.018).abs() < 0.1, "total_mass={m:.3}");
    }
    #[test]
    fn test_write_summary_contains_atoms() {
        let top = parse_minimal();
        let summary = top.write_summary();
        assert!(summary.contains("Atoms"), "summary={summary}");
        assert!(summary.contains('2'), "summary={summary}");
    }
    #[test]
    fn test_from_prmtop_empty_string() {
        let top = AmberTopology::from_prmtop_str("").unwrap();
        assert_eq!(top.n_atoms, 0);
        assert!(top.atoms.is_empty());
    }
    fn minimal_inpcrd() -> String {
        "test coords\n    2\n   1.0000000   2.0000000   3.0000000   4.0000000   5.0000000   6.0000000\n"
            .to_string()
    }
    #[test]
    fn test_amber_coords_parse_basic() {
        let coords = AmberCoordinates::from_str(&minimal_inpcrd()).unwrap();
        assert_eq!(coords.n_atoms, 2);
        assert_eq!(coords.coords.len(), 6);
    }
    #[test]
    fn test_amber_coords_position() {
        let coords = AmberCoordinates::from_str(&minimal_inpcrd()).unwrap();
        let p0 = coords.position(0);
        assert!((p0[0] - 1.0).abs() < 1e-6);
        assert!((p0[1] - 2.0).abs() < 1e-6);
        assert!((p0[2] - 3.0).abs() < 1e-6);
    }
    #[test]
    fn test_amber_coords_position_atom1() {
        let coords = AmberCoordinates::from_str(&minimal_inpcrd()).unwrap();
        let p1 = coords.position(1);
        assert!((p1[0] - 4.0).abs() < 1e-6);
    }
    #[test]
    fn test_amber_coords_no_velocities() {
        let coords = AmberCoordinates::from_str(&minimal_inpcrd()).unwrap();
        assert!(coords.velocities.is_none());
    }
    #[test]
    fn test_amber_coords_restart_roundtrip() {
        let coords = AmberCoordinates::from_str(&minimal_inpcrd()).unwrap();
        let restart = coords.to_restart_string();
        assert!(restart.contains("test coords"));
        assert!(restart.contains("2"));
    }
    #[test]
    fn test_amber_coords_empty() {
        let result = AmberCoordinates::from_str("");
        assert!(result.is_err());
    }
    #[test]
    fn test_amber_coords_with_velocities() {
        let s = "test\n    1\n   1.0000000   2.0000000   3.0000000   0.1000000   0.2000000   0.3000000\n";
        let coords = AmberCoordinates::from_str(s).unwrap();
        assert_eq!(coords.n_atoms, 1);
        assert!(coords.velocities.is_some());
        let vels = coords.velocities.unwrap();
        assert!((vels[0] - 0.1).abs() < 1e-6);
    }
    #[test]
    fn test_extract_lj_empty() {
        let params = extract_lj_parameters("");
        assert!(params.is_empty());
    }
    #[test]
    fn test_extract_dihedrals_empty() {
        let dihedrals = extract_dihedrals("");
        assert!(dihedrals.is_empty());
    }
    #[test]
    fn test_atom_type_names() {
        let top = parse_minimal();
        let types = top.atom_type_names();
        assert!(types.contains(&"N".to_string()));
        assert!(types.contains(&"H".to_string()));
    }
    #[test]
    fn test_n_atom_types() {
        let top = parse_minimal();
        assert_eq!(top.n_atom_types(), 2);
    }
    #[test]
    fn test_atoms_in_residue() {
        let top = parse_minimal();
        let ala_atoms = top.atoms_in_residue("ALA");
        assert_eq!(ala_atoms.len(), 2);
    }
    #[test]
    fn test_atoms_in_residue_missing() {
        let top = parse_minimal();
        let gly_atoms = top.atoms_in_residue("GLY");
        assert!(gly_atoms.is_empty());
    }
    #[test]
    fn test_extract_lj_basic() {
        // Two atom types: NTYPES * (NTYPES + 1) / 2 = 3 pair entries.
        let prmtop = lj_prmtop(
            2,
            "       1       2       2       3",
            "  1.0  2.0  3.0",
            "  0.5  1.0  1.5",
        );
        let params = extract_lj_parameters(&prmtop);
        assert_eq!(params.len(), 3);
        let pairs: Vec<(usize, usize)> = params.iter().map(|p| (p.type_a, p.type_b)).collect();
        assert_eq!(pairs, vec![(0, 0), (0, 1), (1, 1)]);
        assert!((params[0].a_coeff - 1.0).abs() < 1e-10);
        assert!((params[0].b_coeff - 0.5).abs() < 1e-10);
        // The (0, 1) pair is entry 2 of the coefficient tables.
        assert!((params[1].a_coeff - 2.0).abs() < 1e-10);
        assert!((params[1].b_coeff - 1.0).abs() < 1e-10);
        assert!((params[2].a_coeff - 3.0).abs() < 1e-10);
        assert!((params[2].b_coeff - 1.5).abs() < 1e-10);
    }
    #[test]
    fn test_lj_r_min() {
        let lj = LjParameter {
            type_a: 0,
            type_b: 0,
            a_coeff: 2.0,
            b_coeff: 2.0,
        };
        let r_min = lj.r_min();
        assert!((r_min - 2.0f64.powf(1.0 / 6.0)).abs() < 1e-10);
    }
    #[test]
    fn test_lj_epsilon() {
        let lj = LjParameter {
            type_a: 0,
            type_b: 0,
            a_coeff: 4.0,
            b_coeff: 4.0,
        };
        let eps = lj.epsilon();
        assert!((eps - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_lj_zero_b_coeff_r_min() {
        let lj = LjParameter {
            type_a: 0,
            type_b: 0,
            a_coeff: 1.0,
            b_coeff: 0.0,
        };
        assert_eq!(lj.r_min(), 0.0);
    }
    #[test]
    fn test_extract_dihedrals_basic() {
        let prmtop = dihedral_prmtop(4, "       0       3       6       9       1", "");
        let dihedrals = extract_dihedrals(&prmtop);
        assert_eq!(dihedrals.len(), 1);
        let d = &dihedrals[0];
        assert_eq!((d.i, d.j, d.k, d.l), (0, 1, 2, 3));
        assert_eq!(d.n, 3);
        assert!((d.v_n - 2.0 * 0.155555556).abs() < 1e-12, "v_n={}", d.v_n);
        assert_eq!(d.gamma, 0.0);
        assert!(!d.is_improper);
        assert!(!d.skip_14);
    }
    #[test]
    fn test_extract_dihedrals_improper_flag() {
        // A negative FOURTH atom pointer marks an improper torsion.
        let prmtop = dihedral_prmtop(4, "", "       0       3       6      -9       1");
        let dihedrals = extract_dihedrals(&prmtop);
        assert_eq!(dihedrals.len(), 1);
        let d = &dihedrals[0];
        assert!(d.is_improper);
        assert!(!d.skip_14);
        assert_eq!((d.i, d.j, d.k, d.l), (0, 1, 2, 3));
        assert_eq!(d.n, 3);
    }
    fn sample_mdout() -> &'static str {
        " NSTEP =       10  TIME(PS) =  0.010000  TEMP(K) =  298.5  PRESS =    0.0\n\
         Etot   =    -100.00  EKtot   =     50.00  EPtot   =    -150.00\n"
    }
    #[test]
    fn test_parse_energy_frames_single() {
        let frames = parse_energy_frames(sample_mdout());
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].step, 10);
        assert!((frames[0].time_ps - 0.01).abs() < 1e-6);
        assert!((frames[0].temp_k - 298.5).abs() < 1e-3);
        assert!((frames[0].e_tot - (-100.0)).abs() < 1e-6);
        assert!((frames[0].e_kin - 50.0).abs() < 1e-6);
        assert!((frames[0].e_pot - (-150.0)).abs() < 1e-6);
    }
    #[test]
    fn test_parse_energy_frames_empty() {
        let frames = parse_energy_frames("");
        assert!(frames.is_empty());
    }
    #[test]
    fn test_parse_energy_frames_two_frames() {
        let mdout = " NSTEP =        1  TIME(PS) =  0.001000  TEMP(K) =  300.0  PRESS =   0.0\n\
                     Etot   =   -200.00  EKtot   =   100.00  EPtot   =   -300.00\n\
                     \n\
                     NSTEP =        2  TIME(PS) =  0.002000  TEMP(K) =  301.0  PRESS =   0.0\n\
                     Etot   =   -190.00  EKtot   =   110.00  EPtot   =   -300.00\n";
        let frames = parse_energy_frames(mdout);
        assert_eq!(frames.len(), 2);
        assert_eq!(frames[1].step, 2);
    }
    #[test]
    fn test_energy_summary_basic() {
        let frames = vec![
            AmberEnergyFrame {
                step: 1,
                time_ps: 0.0,
                temp_k: 300.0,
                e_tot: -100.0,
                e_kin: 50.0,
                e_pot: -150.0,
            },
            AmberEnergyFrame {
                step: 2,
                time_ps: 0.1,
                temp_k: 300.0,
                e_tot: -200.0,
                e_kin: 50.0,
                e_pot: -250.0,
            },
        ];
        let (mean, min, max) = energy_summary(&frames);
        assert!((mean - (-150.0)).abs() < 1e-10);
        assert!((min - (-200.0)).abs() < 1e-10);
        assert!((max - (-100.0)).abs() < 1e-10);
    }
    #[test]
    fn test_energy_summary_empty() {
        let (mean, min, max) = energy_summary(&[]);
        assert_eq!(mean, 0.0);
        assert_eq!(min, 0.0);
        assert_eq!(max, 0.0);
    }
    #[test]
    fn test_write_inpcrd_basic() {
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let s = write_inpcrd("test restart", &positions, None);
        assert!(s.starts_with("test restart"));
        assert!(s.contains("2"));
        assert!(s.contains("1.0000000"));
    }
    #[test]
    fn test_write_inpcrd_with_velocities() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let velocities = vec![[0.1, 0.2, 0.3]];
        let s = write_inpcrd("vel test", &positions, Some(&velocities));
        assert!(s.contains("0.1000000"));
    }
    #[test]
    fn test_amber99sb_bonds_not_empty() {
        let bonds = amber99sb_bonds();
        assert!(!bonds.is_empty());
        let ct_ct = bonds.iter().find(|b| b.type_a == "CT" && b.type_b == "CT");
        assert!(ct_ct.is_some());
    }
    #[test]
    fn test_amber99sb_angles_not_empty() {
        let angles = amber99sb_angles();
        assert!(!angles.is_empty());
        for a in &angles {
            assert!(
                a.k > 0.0,
                "k should be positive for {}-{}-{}",
                a.type_i,
                a.type_j,
                a.type_k
            );
            assert!(a.theta0_deg > 0.0);
        }
    }
    // ── helpers for the prmtop section tests below ────────────────────────────

    /// Build a prmtop text from `(flag, format, body)` sections.
    fn prmtop_text(sections: &[(&str, &str, &str)]) -> String {
        let mut text = String::from("%VERSION  VERSION_STAMP = V0001.000\n");
        for (flag, format, body) in sections {
            text.push_str(&format!("%FLAG {flag}\n%FORMAT({format})\n"));
            if !body.is_empty() {
                text.push_str(body);
                text.push('\n');
            }
        }
        text
    }
    /// Topology with `natom` atoms and two angle parameter types.
    fn angle_prmtop(natom: usize, inc_h: &str, without_h: &str) -> String {
        prmtop_text(&[
            ("POINTERS", "10I8", &format!("{natom:8}")),
            ("ANGLES_INC_HYDROGEN", "10I8", inc_h),
            ("ANGLES_WITHOUT_HYDROGEN", "10I8", without_h),
            (
                "ANGLE_FORCE_CONSTANT",
                "5E16.8",
                "  5.00000000E+01  8.00000000E+01",
            ),
            (
                "ANGLE_EQUIL_VALUE",
                "5E16.8",
                "  1.90000000E+00  2.10000000E+00",
            ),
        ])
    }
    /// Topology with `natom` atoms and one bond parameter type.
    fn bond_prmtop(natom: usize, inc_h: &str, without_h: &str) -> String {
        prmtop_text(&[
            ("POINTERS", "10I8", &format!("{natom:8}")),
            ("BONDS_INC_HYDROGEN", "10I8", inc_h),
            ("BONDS_WITHOUT_HYDROGEN", "10I8", without_h),
            ("BOND_FORCE_CONSTANT", "5E16.8", "  3.40000000E+02"),
            ("BOND_EQUIL_VALUE", "5E16.8", "  1.09000000E+00"),
        ])
    }
    fn lj_prmtop(n_types: usize, nb_index: &str, a_coeff: &str, b_coeff: &str) -> String {
        prmtop_text(&[
            ("POINTERS", "10I8", &format!("{:8}{:8}", 4, n_types)),
            ("NONBONDED_PARM_INDEX", "10I8", nb_index),
            ("LENNARD_JONES_ACOEF", "5E16.8", a_coeff),
            ("LENNARD_JONES_BCOEF", "5E16.8", b_coeff),
        ])
    }
    /// Topology with `natom` atoms and two dihedral parameter sets:
    /// `PK = 0.155555556, PN = 3, PHASE = 0` and `PK = 1.1, PN = 2, PHASE = pi`.
    fn dihedral_prmtop(natom: usize, inc_h: &str, without_h: &str) -> String {
        prmtop_text(&[
            ("POINTERS", "10I8", &format!("{natom:8}")),
            ("DIHEDRALS_INC_HYDROGEN", "10I8", inc_h),
            ("DIHEDRALS_WITHOUT_HYDROGEN", "10I8", without_h),
            (
                "DIHEDRAL_FORCE_CONSTANT",
                "5E16.8",
                "  1.55555556E-01  1.10000000E+00",
            ),
            (
                "DIHEDRAL_PERIODICITY",
                "5E16.8",
                "  3.00000000E+00  2.00000000E+00",
            ),
            (
                "DIHEDRAL_PHASE",
                "5E16.8",
                "  0.00000000E+00  3.14159274E+00",
            ),
        ])
    }
    fn topology_error(text: &str) -> String {
        match AmberTopology::from_prmtop_str(text) {
            Ok(_) => panic!("the topology must be rejected:\n{text}"),
            Err(message) => message,
        }
    }

    #[test]
    fn test_write_prmtop_output_passes_the_strict_reader() {
        let text = write_prmtop(&parse_minimal());
        let again = AmberTopology::from_prmtop_str(&text).unwrap();
        assert_eq!(again.n_atoms, 2);
        assert_eq!(again.atoms[0].name, "N");
        assert_eq!(again.atoms[1].name, "H");
        assert!((again.atoms[0].mass - 14.01).abs() < 1e-6);
    }

    // ── angles (third atom) ───────────────────────────────────────────────────

    #[test]
    fn test_angle_third_atom_is_parsed() {
        let text = prmtop_text(&[
            ("POINTERS", "10I8", "       3"),
            (
                "ANGLES_WITHOUT_HYDROGEN",
                "10I8",
                "       0       3       6       1",
            ),
            ("ANGLE_FORCE_CONSTANT", "5E16.8", "  6.30000000E+01"),
            ("ANGLE_EQUIL_VALUE", "5E16.8", "  1.94000000E+00"),
        ]);
        let top = AmberTopology::from_prmtop_str(&text).unwrap();
        assert_eq!(top.angles.len(), 1);
        let angle = &top.angles[0];
        assert_eq!((angle.i, angle.j, angle.atom_k), (0, 1, 2));
        assert_eq!(angle.k_idx, 0);
        assert!((angle.k - 63.0).abs() < 1e-12, "k={}", angle.k);
        assert!(
            (angle.theta0 - 1.94).abs() < 1e-12,
            "theta0={}",
            angle.theta0
        );
    }
    #[test]
    fn test_angle_sections_combine_in_file_order_with_their_own_types() {
        let text = angle_prmtop(
            4,
            "       3       0       9       2",
            "       0       3       6       1       6       9       3       1",
        );
        let top = AmberTopology::from_prmtop_str(&text).unwrap();
        let atoms: Vec<(usize, usize, usize, usize)> = top
            .angles
            .iter()
            .map(|a| (a.i, a.j, a.atom_k, a.k_idx))
            .collect();
        assert_eq!(atoms, vec![(1, 0, 3, 1), (0, 1, 2, 0), (2, 3, 1, 0)]);
        assert!((top.angles[0].k - 80.0).abs() < 1e-9);
        assert!((top.angles[0].theta0 - 2.1).abs() < 1e-9);
        assert!((top.angles[1].k - 50.0).abs() < 1e-9);
        assert!((top.angles[2].theta0 - 1.9).abs() < 1e-9);
    }
    #[test]
    fn test_build_angles_reads_the_third_pointer() {
        let angles = build_angles(&[0, 3, 6, 1, 9, 6, 3, 2], &[10.0, 20.0], &[1.0, 2.0]);
        assert_eq!(angles.len(), 2);
        assert_eq!((angles[0].i, angles[0].j, angles[0].atom_k), (0, 1, 2));
        assert_eq!((angles[1].i, angles[1].j, angles[1].atom_k), (3, 2, 1));
        assert_eq!(angles[1].k_idx, 1);
        assert!((angles[1].k - 20.0).abs() < 1e-12);
        assert!((angles[1].theta0 - 2.0).abs() < 1e-12);
        // A trailing incomplete quadruple is ignored.
        assert_eq!(build_angles(&[0, 3, 6, 1, 0, 3], &[1.0], &[1.0]).len(), 1);
    }
    #[test]
    fn test_build_angles_pointer_sign_is_ignored_and_extremes_do_not_panic() {
        // The sign of an angle pointer carries no meaning: |p| / 3 designates
        // the atom, so a negative pointer cannot wrap around to a huge index.
        let angles = build_angles(&[-3, 3, -6, 1], &[10.0], &[1.0]);
        assert_eq!(angles.len(), 1);
        assert_eq!((angles[0].i, angles[0].j, angles[0].atom_k), (1, 1, 2));
        // The extreme 64-bit values must not overflow (`i64::MIN - 1` and
        // `i64::MIN.abs()` both would).
        let extreme = build_angles(&[i64::MIN, i64::MAX, i64::MIN, i64::MIN], &[10.0], &[1.0]);
        assert_eq!(extreme.len(), 1);
        assert!(extreme[0].i >= 1 << 30 && extreme[0].j >= 1 << 30 && extreme[0].atom_k >= 1 << 30);
        // A non-positive type index selects the first parameter type ...
        assert_eq!(extreme[0].k_idx, 0);
        assert!((extreme[0].k - 10.0).abs() < 1e-12);
        // ... and an index beyond the tables yields zero parameters.
        let beyond = build_angles(&[0, 3, 6, i64::MAX], &[10.0], &[1.0]);
        assert_eq!((beyond[0].k, beyond[0].theta0), (0.0, 0.0));
    }
    #[test]
    fn test_build_bonds_pointer_and_type_arithmetic_never_overflows() {
        let bonds = build_bonds(&[0, 3, 1, 6, 9, 2], &[300.0, 400.0], &[1.0, 2.0]);
        assert_eq!(bonds.len(), 2);
        assert_eq!((bonds[0].i, bonds[0].j), (0, 1));
        assert_eq!((bonds[1].i, bonds[1].j), (2, 3));
        assert!((bonds[1].k - 400.0).abs() < 1e-12);
        assert!((bonds[1].r0 - 2.0).abs() < 1e-12);
        let extreme = build_bonds(&[i64::MIN, -3, i64::MIN], &[300.0], &[1.0]);
        assert_eq!(extreme.len(), 1);
        assert_eq!(extreme[0].j, 1);
        assert!((extreme[0].k - 300.0).abs() < 1e-12);
        // A trailing incomplete triplet is ignored.
        assert_eq!(build_bonds(&[0, 3, 1, 6], &[1.0], &[1.0]).len(), 1);
    }
    #[test]
    fn test_angle_atoms_beyond_natom_are_rejected() {
        // Each of the three atom pointers is checked in turn.
        for ints in [
            "      12       3       6       1",
            "       0      12       6       1",
            "       0       3      12       1",
        ] {
            let message = topology_error(&angle_prmtop(3, "", ints));
            assert!(
                message.contains("ANGLES_WITHOUT_HYDROGEN") && message.contains("pointer 12"),
                "{message}"
            );
        }
        let message = topology_error(&angle_prmtop(3, "       0       3       9       1", ""));
        assert!(message.contains("ANGLES_INC_HYDROGEN"), "{message}");
        // The same record is fine once the topology has enough atoms.
        let top = AmberTopology::from_prmtop_str(&angle_prmtop(
            4,
            "       0       3       9       1",
            "",
        ))
        .unwrap();
        assert_eq!(top.angles[0].atom_k, 3);
    }
    #[test]
    fn test_angle_malformed_records_are_rejected() {
        // negative pointer, pointer that is not a multiple of three
        for ints in [
            "       0       3      -3       1",
            "       0       3       4       1",
        ] {
            let message = topology_error(&angle_prmtop(3, "", ints));
            assert!(message.contains("atom pointer"), "{message}");
        }
        // type index 0 and past the end of the tables
        for ints in [
            "       0       3       6       0",
            "       0       3       6       3",
        ] {
            let message = topology_error(&angle_prmtop(3, "", ints));
            assert!(message.contains("parameter type index"), "{message}");
        }
        // incomplete trailing record
        let message = topology_error(&angle_prmtop(3, "", "       0       3       6"));
        assert!(message.contains("whole records"), "{message}");
    }
    #[test]
    fn test_bond_atoms_and_types_are_validated() {
        let top = AmberTopology::from_prmtop_str(&bond_prmtop(2, "       0       3       1", ""))
            .unwrap();
        assert_eq!((top.bonds[0].i, top.bonds[0].j), (0, 1));
        let message = topology_error(&bond_prmtop(2, "", "       0       6       1"));
        assert!(message.contains("BONDS_WITHOUT_HYDROGEN"), "{message}");
        let message = topology_error(&bond_prmtop(2, "       0       3       2", ""));
        assert!(message.contains("parameter type index"), "{message}");
        let message = topology_error(&bond_prmtop(2, "       0       3", ""));
        assert!(message.contains("whole records"), "{message}");
    }
    #[test]
    fn test_prmtop_counters_are_sanity_checked() {
        let negative = prmtop_text(&[("POINTERS", "10I8", "      -1")]);
        assert!(topology_error(&negative).contains("NATOM"));
        // An atom count that no input of this size could describe must not
        // allocate a gigantic atom table.
        let huge = prmtop_text(&[("POINTERS", "10I8", "99999999999999")]);
        assert!(topology_error(&huge).contains("more than"));
    }
    #[test]
    fn test_validate_topology_checks_the_third_angle_atom() {
        let atom = |name: &str| AmberAtom {
            name: name.to_string(),
            residue_name: "R".to_string(),
            charge: 0.0,
            mass: 12.0,
            atom_type: "C".to_string(),
        };
        let topo = AmberTopology {
            title: "t".to_string(),
            atoms: vec![atom("A"), atom("B")],
            bonds: Vec::new(),
            angles: vec![AmberAngle {
                i: 0,
                j: 1,
                atom_k: 5,
                k_idx: 0,
                k: 50.0,
                theta0: 1.9,
            }],
            n_atoms: 2,
            n_bonds: 0,
        };
        let issues = crate::amber::validate_topology(&topo);
        assert_eq!(issues.len(), 1, "issues: {issues:?}");
        assert!(issues[0].contains("Angle 0") && issues[0].contains("k=5"));
    }

    // ── Lennard-Jones pair table ──────────────────────────────────────────────

    #[test]
    fn test_extract_lj_three_types_follow_the_nonbonded_index() {
        // AMBER numbering of the 6 unordered pairs of 3 types.
        let prmtop = lj_prmtop(
            3,
            "       1       2       4       2       3       5       4       5       6",
            "  1.0  2.0  3.0  4.0  5.0  6.0",
            "  0.1  0.2  0.3  0.4  0.5  0.6",
        );
        let params = extract_lj_parameters(&prmtop);
        let pairs: Vec<(usize, usize)> = params.iter().map(|p| (p.type_a, p.type_b)).collect();
        assert_eq!(pairs, vec![(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]);
        let a: Vec<f64> = params.iter().map(|p| p.a_coeff).collect();
        let b: Vec<f64> = params.iter().map(|p| p.b_coeff).collect();
        assert_eq!(a, vec![1.0, 2.0, 4.0, 3.0, 5.0, 6.0]);
        assert_eq!(b, vec![0.1, 0.2, 0.4, 0.3, 0.5, 0.6]);
    }
    #[test]
    fn test_extract_lj_skips_hydrogen_bond_pairs() {
        let prmtop = lj_prmtop(
            2,
            "       1      -1      -1       2",
            "  1.0  2.0",
            "  0.5  1.5",
        );
        let params = extract_lj_parameters(&prmtop);
        let pairs: Vec<(usize, usize)> = params.iter().map(|p| (p.type_a, p.type_b)).collect();
        assert_eq!(pairs, vec![(0, 0), (1, 1)]);
        assert!((params[1].a_coeff - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_extract_lj_never_guesses_the_pair_assignment() {
        let tables = ("  1.0  2.0  3.0", "  0.5  1.0  1.5");
        // The tables alone (no POINTERS / NONBONDED_PARM_INDEX) used to be
        // mislabelled as the diagonal (i, i) types.
        let tables_only = prmtop_text(&[
            ("LENNARD_JONES_ACOEF", "5E16.8", tables.0),
            ("LENNARD_JONES_BCOEF", "5E16.8", tables.1),
        ]);
        assert!(extract_lj_parameters(&tables_only).is_empty());
        // POINTERS but no index.
        let no_index = prmtop_text(&[
            ("POINTERS", "10I8", "       4       2"),
            ("LENNARD_JONES_ACOEF", "5E16.8", tables.0),
            ("LENNARD_JONES_BCOEF", "5E16.8", tables.1),
        ]);
        assert!(extract_lj_parameters(&no_index).is_empty());
        // Index shorter than NTYPES * NTYPES.
        let short = lj_prmtop(2, "       1       2       2", tables.0, tables.1);
        assert!(extract_lj_parameters(&short).is_empty());
        // Index pointing past the coefficient tables.
        let beyond = lj_prmtop(2, "       1       2       2       4", tables.0, tables.1);
        assert!(extract_lj_parameters(&beyond).is_empty());
        // Zero types.
        let zero = lj_prmtop(0, "", tables.0, tables.1);
        assert!(extract_lj_parameters(&zero).is_empty());
    }

    // ── dihedrals ─────────────────────────────────────────────────────────────

    #[test]
    fn test_extract_dihedrals_skip_14_flag() {
        // A negative THIRD atom pointer means "skip the 1-4 interaction".
        let prmtop = dihedral_prmtop(4, "       0       3      -6       9       1", "");
        let d = &extract_dihedrals(&prmtop)[0];
        assert!(d.skip_14);
        assert!(!d.is_improper);
        assert_eq!((d.i, d.j, d.k, d.l), (0, 1, 2, 3));
    }
    #[test]
    fn test_extract_dihedrals_both_flags_and_parameter_lookup() {
        let prmtop = dihedral_prmtop(
            4,
            "       0       3      -6      -9       2",
            "       0       3       6       9       2       3       0       6       9       1",
        );
        let dihedrals = extract_dihedrals(&prmtop);
        assert_eq!(dihedrals.len(), 3);
        // hydrogen section first, file order inside a section
        let d = &dihedrals[0];
        assert!(d.skip_14 && d.is_improper);
        assert_eq!(d.n, 2);
        assert!((d.v_n - 2.2).abs() < 1e-9, "v_n={}", d.v_n);
        assert!((d.gamma - 3.14159274).abs() < 1e-9, "gamma={}", d.gamma);
        assert!(!dihedrals[1].skip_14 && !dihedrals[1].is_improper);
        assert_eq!(dihedrals[1].n, 2);
        assert_eq!((dihedrals[2].i, dihedrals[2].j), (1, 0));
        assert_eq!(dihedrals[2].n, 3);
    }
    #[test]
    fn test_extract_dihedrals_reproduce_the_amber_energy() {
        // AMBER: E = PK * (1 + cos(PN * phi - PHASE)).
        let prmtop = dihedral_prmtop(4, "       0       3       6       9       2", "");
        let d = &extract_dihedrals(&prmtop)[0];
        let phi = std::f64::consts::FRAC_PI_2;
        let amber = 1.1 * (1.0 + (2.0 * phi - 3.14159274_f64).cos());
        let ours = crate::amber::dihedral_energy(d.v_n, f64::from(d.n), phi, d.gamma);
        assert!((ours - amber).abs() < 1e-9, "ours={ours}, amber={amber}");
        assert!((ours - 2.2).abs() < 1e-6);
    }
    #[test]
    fn test_extract_dihedrals_skips_undecodable_entries() {
        let prmtop = dihedral_prmtop(
            4,
            // parameter set 3 does not exist, set 0 is not one-based
            "       0       3       6       9       3       0       3       6       9       0",
            // atom 4 does not exist; pointer 4 is not a multiple of three;
            // the last record is valid
            "       0       3       6      12       1       0       3       6       4       1       0       3       6       9       1",
        );
        let dihedrals = extract_dihedrals(&prmtop);
        assert_eq!(dihedrals.len(), 1);
        assert_eq!((dihedrals[0].k, dihedrals[0].l), (2, 3));
        // Without POINTERS the atom count is unknown, so only the parameter
        // index and the pointer arithmetic are checked.
        let no_pointers = "%FLAG DIHEDRALS_INC_HYDROGEN\n%FORMAT(10I8)\n0 3 6 12 1\n%FLAG DIHEDRAL_FORCE_CONSTANT\n%FORMAT(5E16.8)\n1.0\n%FLAG DIHEDRAL_PERIODICITY\n%FORMAT(5E16.8)\n2.0\n%FLAG DIHEDRAL_PHASE\n%FORMAT(5E16.8)\n0.0\n";
        assert_eq!(extract_dihedrals(no_pointers).len(), 1);
    }
    #[test]
    fn test_extract_dihedrals_negative_periodicity_is_its_magnitude() {
        let prmtop = prmtop_text(&[
            ("DIHEDRALS_INC_HYDROGEN", "10I8", "0 3 6 9 1"),
            ("DIHEDRAL_FORCE_CONSTANT", "5E16.8", "0.5"),
            ("DIHEDRAL_PERIODICITY", "5E16.8", "-2.0"),
            ("DIHEDRAL_PHASE", "5E16.8", "0.0"),
        ]);
        let dihedrals = extract_dihedrals(&prmtop);
        assert_eq!(dihedrals[0].n, 2);
        assert!((dihedrals[0].v_n - 1.0).abs() < 1e-12);
    }

    // ── MDCRD ─────────────────────────────────────────────────────────────────

    /// Coordinates of one frame as `10F8.3` lines (plus a `3F8.3` box line).
    fn mdcrd_frame(coords: &[f64], box_dims: Option<[f64; 3]>) -> String {
        let mut text = String::new();
        for line in coords.chunks(10) {
            text.extend(line.iter().map(|v| format!("{v:8.3}")));
            text.push('\n');
        }
        if let Some(dims) = box_dims {
            text.extend(dims.iter().map(|v| format!("{v:8.3}")));
            text.push('\n');
        }
        text
    }
    #[test]
    fn test_mdcrd_boxless_final_frame_does_not_panic() {
        let frames = parse_mdcrd("t\n1 2 3 4 5 6\n", 2, true);
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].coords, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(frames[0].box_dims, None);
    }
    #[test]
    fn test_mdcrd_frames_with_boxes() {
        // 5 atoms = 15 coordinates = two lines per frame, plus the box line.
        let mut text = String::from("title\n");
        let frame =
            |f: usize| -> Vec<f64> { (0..15).map(|i| (f * 100 + i) as f64 * 0.5).collect() };
        for f in 0..3 {
            text.push_str(&mdcrd_frame(&frame(f), Some([30.0 + f as f64, 40.0, 50.0])));
        }
        let frames = parse_mdcrd(&text, 5, true);
        assert_eq!(frames.len(), 3);
        for (f, parsed) in frames.iter().enumerate() {
            assert_eq!(parsed.coords, frame(f));
            assert_eq!(parsed.box_dims, Some([30.0 + f as f64, 40.0, 50.0]));
            assert_eq!(parsed.n_atoms(), 5);
        }
        // Dropping the last box line keeps the frame and loses only the box.
        let truncated = text
            .trim_end()
            .rsplit_once('\n')
            .map(|(head, _)| head.to_string());
        let frames = parse_mdcrd(&truncated.unwrap(), 5, true);
        assert_eq!(frames.len(), 3);
        assert_eq!(frames[2].box_dims, None);
        assert_eq!(frames[1].box_dims, Some([31.0, 40.0, 50.0]));
    }
    #[test]
    fn test_mdcrd_two_thousand_frames_parse_completely() {
        // 4 atoms = 12 coordinates = lines of 10 + 2 values.  Four-digit values
        // fill their F8.3 column, so neighbouring fields touch.
        let coords =
            |f: usize| -> Vec<f64> { (0..12).map(|i| f as f64 + i as f64 * 0.5).collect() };
        let mut text = String::from("2000 frames\n");
        for f in 0..2000 {
            text.push_str(&mdcrd_frame(
                &coords(f),
                Some([50.0 + f as f64, 60.0, 70.0]),
            ));
        }
        let frames = parse_mdcrd(&text, 4, true);
        assert_eq!(frames.len(), 2000);
        for f in [0, 1, 999, 1000, 1998, 1999] {
            assert_eq!(frames[f].coords, coords(f), "frame {f}");
            assert_eq!(
                frames[f].box_dims,
                Some([50.0 + f as f64, 60.0, 70.0]),
                "frame {f}"
            );
        }
        assert!(
            text.lines()
                .any(|l| l.split_whitespace().any(|t| t.len() > 8))
        );
    }
    #[test]
    fn test_mdcrd_adjacent_negative_fields_parse() {
        let coords = [1.0, -100.0, 3.0, -100.0, -200.5, 6.0];
        let text = format!("t\n{}", mdcrd_frame(&coords, None));
        assert!(text.contains("1.000-100.000"), "{text}");
        let frames = parse_mdcrd(&text, 2, false);
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].coords, coords.to_vec());
        assert_eq!(frames[0].position(1), Some([-100.0, -200.5, 6.0]));
    }
    #[test]
    fn test_mdcrd_touching_box_fields() {
        let text = format!(
            "t\n{}",
            mdcrd_frame(&[1.0, 2.0, 3.0], Some([1000.0, 1000.5, 1001.0]))
        );
        assert!(text.contains("1000.0001000.5001001.000"), "{text}");
        let frames = parse_mdcrd(&text, 1, true);
        assert_eq!(frames[0].box_dims, Some([1000.0, 1000.5, 1001.0]));
    }
    #[test]
    fn test_mdcrd_short_box_line_keeps_the_frame() {
        let frames = parse_mdcrd("t\n   1.000   2.000   3.000\n  50.000  60.000\n", 1, true);
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].box_dims, None);
    }
    #[test]
    fn test_mdcrd_stops_at_an_incomplete_frame() {
        // 5 atoms = two lines per frame, no box.
        let coords: Vec<f64> = (0..15).map(f64::from).collect();
        let mut text = String::from("t\n");
        for _ in 0..3 {
            text.push_str(&mdcrd_frame(&coords, None));
        }
        // a fourth frame with only its first line
        text.push_str(&mdcrd_frame(&coords[..10], None));
        assert_eq!(parse_mdcrd(&text, 5, false).len(), 3);
        // a fourth frame whose second line is short
        let mut short = String::from("t\n");
        short.push_str(&mdcrd_frame(&coords, None));
        short.push_str(&mdcrd_frame(&coords[..12], None));
        assert_eq!(parse_mdcrd(&short, 5, false).len(), 1);
        // a second frame that holds an extra value
        let mut long = String::from("t\n");
        long.push_str(&mdcrd_frame(&coords, None));
        long.push_str(&mdcrd_frame(&[coords.as_slice(), &[99.0]].concat(), None));
        assert_eq!(parse_mdcrd(&long, 5, false).len(), 1);
    }
    #[test]
    fn test_mdcrd_stops_at_garbage() {
        let good = mdcrd_frame(&[1.0, 2.0, 3.0], None);
        let text = format!("t\n{good}   1.000   x.000   3.000\n{good}");
        let frames = parse_mdcrd(&text, 1, false);
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].coords, vec![1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_mdcrd_degenerate_atom_counts() {
        let text = format!("t\n{}", mdcrd_frame(&[1.0, 2.0, 3.0], None));
        assert!(parse_mdcrd(&text, 0, false).is_empty());
        assert!(parse_mdcrd(&text, 0, true).is_empty());
        assert!(parse_mdcrd(&text, usize::MAX, false).is_empty());
        assert!(parse_mdcrd(&text, usize::MAX / 3, true).is_empty());
        assert!(parse_mdcrd("", 1, false).is_empty());
    }
    #[test]
    fn test_mdcrd_accepts_crlf_line_endings() {
        let frames = parse_mdcrd(
            "t\r\n   1.000   2.000   3.000\r\n  50.000  60.000  70.000\r\n",
            1,
            true,
        );
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].box_dims, Some([50.0, 60.0, 70.0]));
    }

    // ── mutation sweeps ───────────────────────────────────────────────────────

    #[test]
    fn test_mutated_prmtop_text_never_yields_out_of_range_terms() {
        use crate::fortran_fields::test_support::{SplitMix64, mutate_text};
        let seeds = [
            minimal_prmtop(),
            angle_prmtop(
                4,
                "       3       0       9       2",
                "       0       3       6       1       6       9       3       1",
            ),
            bond_prmtop(3, "       0       3       1", "       3       6       1"),
            dihedral_prmtop(
                4,
                "       0       3      -6      -9       2",
                "       0       3       6       9       2       3       0       6       9       1",
            ),
            lj_prmtop(
                3,
                "       1       2       4       2       3       5       4       5       6",
                "  1.0  2.0  3.0  4.0  5.0  6.0",
                "  0.1  0.2  0.3  0.4  0.5  0.6",
            ),
        ];
        let mut rng = SplitMix64::new(0xA2B3_0001);
        let mut accepted = 0usize;
        let mut rejected = 0usize;
        for seed in &seeds {
            for round in 0..1200 {
                let text = mutate_text(seed, &mut rng, 1 + round % 4);
                match AmberTopology::from_prmtop_str(&text) {
                    Ok(top) => {
                        accepted += 1;
                        assert_eq!(top.atoms.len(), top.n_atoms);
                        for bond in &top.bonds {
                            assert!(bond.i < top.n_atoms && bond.j < top.n_atoms, "{text}");
                        }
                        for angle in &top.angles {
                            assert!(
                                angle.i < top.n_atoms
                                    && angle.j < top.n_atoms
                                    && angle.atom_k < top.n_atoms,
                                "{text}"
                            );
                        }
                    }
                    Err(_) => rejected += 1,
                }
                let pointers = extract_flag_integers(&text, "POINTERS");
                let n_atoms = pointers.first().and_then(|&n| usize::try_from(n).ok());
                for d in extract_dihedrals(&text) {
                    if let Some(n) = n_atoms {
                        assert!(d.i < n && d.j < n && d.k < n && d.l < n, "{text}");
                    }
                    assert!(
                        d.v_n.is_finite() && d.gamma.is_finite() && d.n >= 0,
                        "{text}"
                    );
                }
                let n_types = pointers
                    .get(1)
                    .and_then(|&n| usize::try_from(n).ok())
                    .unwrap_or(0);
                for pair in extract_lj_parameters(&text) {
                    assert!(
                        pair.type_a <= pair.type_b && pair.type_b < n_types,
                        "{text}"
                    );
                }
            }
        }
        assert!(accepted > 0, "no damaged topology was accepted");
        assert!(rejected > 0, "no damaged topology was rejected");
    }
    #[test]
    fn test_mutated_mdcrd_text_keeps_whole_frames() {
        use crate::fortran_fields::test_support::{SplitMix64, mutate_text};
        let mut seed = String::from("title\n");
        for frame in 0..3 {
            let coords: Vec<f64> = (0..12)
                .map(|i| f64::from(frame) * 100.0 - 150.0 + f64::from(i) * 0.5)
                .collect();
            seed.push_str(&mdcrd_frame(&coords, Some([40.0, 41.0, 42.0])));
        }
        let mut rng = SplitMix64::new(0x3DC2_0001);
        let mut frames_seen = 0usize;
        for n_atoms in [1usize, 2, 4, 5] {
            for has_box in [false, true] {
                for round in 0..300 {
                    let text = mutate_text(&seed, &mut rng, 1 + round % 4);
                    for frame in parse_mdcrd(&text, n_atoms, has_box) {
                        frames_seen += 1;
                        assert_eq!(frame.coords.len(), 3 * n_atoms, "{text}");
                        assert!(has_box || frame.box_dims.is_none(), "{text}");
                    }
                }
            }
        }
        assert!(frames_seen > 0, "no damaged trajectory yielded a frame");
    }
}
