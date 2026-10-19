// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! U.S. NRC Regulatory Guide 1.99, Revision 2 (May 1988), "Radiation
//! Embrittlement of Reactor Vessel Materials", Regulatory Positions 1.1 and 1.2.
//!
//! * **Position 1.1** — the adjustment of the reference temperature
//!   `ΔRT_NDT = CF · f^(0.28 − 0.10·log₁₀ f)` (Equation 2).  The chemistry
//!   factor `CF` is read from Table 1 (welds) or Table 2 (base metal: plates and
//!   forgings) as a function of the best-estimate copper and nickel contents,
//!   and `f` is the fast fluence in units of 10¹⁹ n/cm² (E > 1 MeV).  Linear
//!   interpolation in the tables is permitted by the guide.
//! * **Position 1.2** — the decrease of the Charpy upper-shelf energy, read from
//!   Figure 2 as a function of copper content and fluence (again with linear
//!   interpolation).
//!
//! Source: U.S. NRC, Regulatory Guide 1.99, Revision 2, May 1988 (NRC ADAMS
//! accession ML003740284).  Tables 1 and 2 are transcribed from the guide (the
//! two tables were checked cell by cell against an independent copy of the
//! text, which also appears in 10 CFR 50.61).  Figure 2 is a chart,
//! not a table: it was digitized from the printed figure.  Every curve of the
//! figure is a straight line on its log-log axes, so each is stored by its
//! value at 10¹⁸ n/cm² and its slope; the digitized lines reproduce the printed
//! curves to within about ±0.9 % of their value (±0.2 % rms).

/// Lower end of the fluence range of the guide (Figure 1), in n/cm²
/// (E > 1 MeV).  At or below this fluence no embrittlement is predicted.
pub const RG199_MIN_FLUENCE_N_CM2: f64 = 1.0e17;

/// Product form of a reactor pressure vessel beltline material.
///
/// It selects the chemistry-factor table of Regulatory Position 1.1 and the
/// copper labelling of the upper-shelf curves of Position 1.2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RpvProductForm {
    /// Plates and forgings (Table 2).
    BaseMetal,
    /// Weld metal (Table 1).
    Weld,
}

/// Number of copper rows (0.00 … 0.40 wt-% in steps of 0.01) of the tables.
const CU_ROWS: usize = 41;
/// Number of nickel columns (0 … 1.2 wt-% in steps of 0.2) of the tables.
const NI_COLS: usize = 7;
/// Copper spacing of the table rows (wt-%).
const CU_STEP: f64 = 0.01;
/// Nickel spacing of the table columns (wt-%).
const NI_STEP: f64 = 0.20;
/// Largest copper content of the tables (wt-%).
const CU_MAX: f64 = 0.40;
/// Largest nickel content of the tables (wt-%).
const NI_MAX: f64 = 1.20;

type ChemistryFactorTable = [[u16; NI_COLS]; CU_ROWS];

/// Table 1 of RG 1.99 Rev. 2: chemistry factor for welds, °F.
///
/// Rows: copper 0.00 … 0.40 wt-% in steps of 0.01; columns: nickel 0, 0.2, … 1.2 wt-%.
static WELD_CF_F: ChemistryFactorTable = [
    [20, 20, 20, 20, 20, 20, 20],        // Cu 0.00
    [20, 20, 20, 20, 20, 20, 20],        // Cu 0.01
    [21, 26, 27, 27, 27, 27, 27],        // Cu 0.02
    [22, 35, 41, 41, 41, 41, 41],        // Cu 0.03
    [24, 43, 54, 54, 54, 54, 54],        // Cu 0.04
    [26, 49, 67, 68, 68, 68, 68],        // Cu 0.05
    [29, 52, 77, 82, 82, 82, 82],        // Cu 0.06
    [32, 55, 85, 95, 95, 95, 95],        // Cu 0.07
    [36, 58, 90, 106, 108, 108, 108],    // Cu 0.08
    [40, 61, 94, 115, 122, 122, 122],    // Cu 0.09
    [44, 65, 97, 122, 133, 135, 135],    // Cu 0.10
    [49, 68, 101, 130, 144, 148, 148],   // Cu 0.11
    [52, 72, 103, 135, 153, 161, 161],   // Cu 0.12
    [58, 76, 106, 139, 162, 172, 176],   // Cu 0.13
    [61, 79, 109, 142, 168, 182, 188],   // Cu 0.14
    [66, 84, 112, 146, 175, 191, 200],   // Cu 0.15
    [70, 88, 115, 149, 178, 199, 211],   // Cu 0.16
    [75, 92, 119, 151, 184, 207, 221],   // Cu 0.17
    [79, 95, 122, 154, 187, 214, 230],   // Cu 0.18
    [83, 100, 126, 157, 191, 220, 238],  // Cu 0.19
    [88, 104, 129, 160, 194, 223, 245],  // Cu 0.20
    [92, 108, 133, 164, 197, 229, 252],  // Cu 0.21
    [97, 112, 137, 167, 200, 232, 257],  // Cu 0.22
    [101, 117, 140, 169, 203, 236, 263], // Cu 0.23
    [105, 121, 144, 173, 206, 239, 268], // Cu 0.24
    [110, 126, 148, 176, 209, 243, 272], // Cu 0.25
    [113, 130, 151, 180, 212, 246, 276], // Cu 0.26
    [119, 134, 155, 184, 216, 249, 280], // Cu 0.27
    [122, 138, 160, 187, 218, 251, 284], // Cu 0.28
    [128, 142, 164, 191, 222, 254, 287], // Cu 0.29
    [131, 146, 167, 194, 225, 257, 290], // Cu 0.30
    [136, 151, 172, 198, 228, 260, 293], // Cu 0.31
    [140, 155, 175, 202, 231, 263, 296], // Cu 0.32
    [144, 160, 180, 205, 234, 266, 299], // Cu 0.33
    [149, 164, 184, 209, 238, 269, 302], // Cu 0.34
    [153, 168, 187, 212, 241, 272, 305], // Cu 0.35
    [158, 172, 191, 216, 245, 275, 308], // Cu 0.36
    [162, 177, 196, 220, 248, 278, 311], // Cu 0.37
    [166, 182, 200, 223, 250, 281, 314], // Cu 0.38
    [171, 185, 203, 227, 254, 285, 317], // Cu 0.39
    [175, 189, 207, 231, 257, 288, 320], // Cu 0.40
];

/// Table 2 of RG 1.99 Rev. 2: chemistry factor for base metal (plates and
/// forgings), °F.
///
/// Rows: copper 0.00 … 0.40 wt-% in steps of 0.01; columns: nickel 0, 0.2, … 1.2 wt-%.
static BASE_CF_F: ChemistryFactorTable = [
    [20, 20, 20, 20, 20, 20, 20],        // Cu 0.00
    [20, 20, 20, 20, 20, 20, 20],        // Cu 0.01
    [20, 20, 20, 20, 20, 20, 20],        // Cu 0.02
    [20, 20, 20, 20, 20, 20, 20],        // Cu 0.03
    [22, 26, 26, 26, 26, 26, 26],        // Cu 0.04
    [25, 31, 31, 31, 31, 31, 31],        // Cu 0.05
    [28, 37, 37, 37, 37, 37, 37],        // Cu 0.06
    [31, 43, 44, 44, 44, 44, 44],        // Cu 0.07
    [34, 48, 51, 51, 51, 51, 51],        // Cu 0.08
    [37, 53, 58, 58, 58, 58, 58],        // Cu 0.09
    [41, 58, 65, 65, 67, 67, 67],        // Cu 0.10
    [45, 62, 72, 74, 77, 77, 77],        // Cu 0.11
    [49, 67, 79, 83, 86, 86, 86],        // Cu 0.12
    [53, 71, 85, 91, 96, 96, 96],        // Cu 0.13
    [57, 75, 91, 100, 105, 106, 106],    // Cu 0.14
    [61, 80, 99, 110, 115, 117, 117],    // Cu 0.15
    [65, 84, 104, 118, 123, 125, 125],   // Cu 0.16
    [69, 88, 110, 127, 132, 135, 135],   // Cu 0.17
    [73, 92, 115, 134, 141, 144, 144],   // Cu 0.18
    [78, 97, 120, 142, 150, 154, 154],   // Cu 0.19
    [82, 102, 125, 149, 159, 164, 165],  // Cu 0.20
    [86, 107, 129, 155, 167, 172, 174],  // Cu 0.21
    [91, 112, 134, 161, 176, 181, 184],  // Cu 0.22
    [95, 117, 138, 167, 184, 190, 194],  // Cu 0.23
    [100, 121, 143, 172, 191, 199, 204], // Cu 0.24
    [104, 126, 148, 176, 199, 208, 214], // Cu 0.25
    [109, 130, 151, 180, 205, 216, 221], // Cu 0.26
    [114, 134, 155, 184, 211, 225, 230], // Cu 0.27
    [119, 138, 160, 187, 216, 233, 239], // Cu 0.28
    [124, 142, 164, 191, 221, 241, 248], // Cu 0.29
    [129, 146, 167, 194, 225, 249, 257], // Cu 0.30
    [134, 151, 172, 198, 228, 255, 266], // Cu 0.31
    [139, 155, 175, 202, 231, 260, 274], // Cu 0.32
    [144, 160, 180, 205, 234, 264, 282], // Cu 0.33
    [149, 164, 184, 209, 238, 268, 290], // Cu 0.34
    [153, 168, 187, 212, 241, 272, 298], // Cu 0.35
    [158, 173, 191, 216, 245, 275, 303], // Cu 0.36
    [162, 177, 196, 220, 248, 278, 308], // Cu 0.37
    [166, 182, 200, 223, 250, 281, 313], // Cu 0.38
    [171, 185, 203, 227, 254, 285, 317], // Cu 0.39
    [175, 189, 207, 231, 257, 288, 320], // Cu 0.40
];

/// Bilinear interpolation of a chemistry-factor table (°F).
///
/// Copper and nickel are clamped to the range of the tables (0 … 0.40 and
/// 0 … 1.20 wt-%); a NaN input propagates to the result.
fn interpolate_chemistry_factor(
    table: &ChemistryFactorTable,
    copper_wt: f64,
    nickel_wt: f64,
) -> f64 {
    if copper_wt.is_nan() || nickel_wt.is_nan() {
        return f64::NAN;
    }
    let row_pos = copper_wt.clamp(0.0, CU_MAX) / CU_STEP;
    let col_pos = nickel_wt.clamp(0.0, NI_MAX) / NI_STEP;
    // Lower cell corner; the clamped upper edge belongs to the last cell.
    let row = (row_pos.floor() as usize).min(CU_ROWS - 2);
    let col = (col_pos.floor() as usize).min(NI_COLS - 2);
    let row_frac = (row_pos - row as f64).clamp(0.0, 1.0);
    let col_frac = (col_pos - col as f64).clamp(0.0, 1.0);
    let at = |r: usize, c: usize| {
        table
            .get(r)
            .and_then(|line| line.get(c))
            .map_or(f64::NAN, |&value| f64::from(value))
    };
    let lower = at(row, col) * (1.0 - col_frac) + at(row, col + 1) * col_frac;
    let upper = at(row + 1, col) * (1.0 - col_frac) + at(row + 1, col + 1) * col_frac;
    lower * (1.0 - row_frac) + upper * row_frac
}

/// Chemistry factor CF in °F from RG 1.99 Rev. 2 Table 1 (welds) or Table 2
/// (base metal), bilinearly interpolated in copper and nickel content.
///
/// `copper_wt` and `nickel_wt` are the best-estimate contents in weight-percent.
/// Values beyond the tabulated range (copper above 0.40 %, nickel above 1.20 %)
/// are clamped to the table edge; the guide requires such cases to be justified
/// by data.  The factor never drops below 20 °F, the floor of both tables.
pub fn chemistry_factor_fahrenheit(form: RpvProductForm, copper_wt: f64, nickel_wt: f64) -> f64 {
    let table = match form {
        RpvProductForm::Weld => &WELD_CF_F,
        RpvProductForm::BaseMetal => &BASE_CF_F,
    };
    interpolate_chemistry_factor(table, copper_wt, nickel_wt)
}

/// Chemistry factor CF in °C: [`chemistry_factor_fahrenheit`] converted as a
/// temperature *difference* (multiplied by 5/9, no offset).
pub fn chemistry_factor_for(form: RpvProductForm, copper_wt: f64, nickel_wt: f64) -> f64 {
    chemistry_factor_fahrenheit(form, copper_wt, nickel_wt) * 5.0 / 9.0
}

/// The fluence factor `f^(0.28 − 0.10·log₁₀ f)` of Equation 2 / Figure 1, with
/// the fluence `f` taken in units of 10¹⁹ n/cm² (E > 1 MeV).
///
/// Returns 0 for fluences at or below [`RG199_MIN_FLUENCE_N_CM2`], the lower end
/// of the guide's range.  Above 10²⁰ n/cm² the guide is not validated (the
/// factor even decreases beyond about 2.5·10²⁰ n/cm²), so results there are an
/// extrapolation.
pub fn rg199_fluence_factor(fluence_n_cm2: f64) -> f64 {
    if fluence_n_cm2 <= RG199_MIN_FLUENCE_N_CM2 {
        return 0.0;
    }
    let f = fluence_n_cm2 / 1.0e19;
    f.powf(0.28 - 0.10 * f.log10())
}

/// One straight line of Figure 2 on its log-log axes:
/// `decrease [%] = pct_at_1e18 · (fluence / 10¹⁸ n/cm²)^slope`.
#[derive(Debug, Clone, Copy)]
struct UseLine {
    /// Decrease in upper-shelf energy at 10¹⁸ n/cm² (percent).
    pct_at_1e18: f64,
    /// Slope on the log-log axes.
    slope: f64,
}

impl UseLine {
    /// Decrease (percent) at `fluence_n_cm2`.
    fn percent_at(&self, fluence_n_cm2: f64) -> f64 {
        self.pct_at_1e18 * (fluence_n_cm2 / 1.0e18).powf(self.slope)
    }
}

/// Copper contents (wt-%) of the plotted curves, in the *base-metal* labelling
/// of Figure 2, ascending: 0.10 … 0.35.  The weld labels of the figure are the
/// base-metal labels minus [`USE_WELD_COPPER_OFFSET`].
const USE_BASE_COPPER: [f64; 6] = [0.10, 0.15, 0.20, 0.25, 0.30, 0.35];

/// Copper spacing (wt-%) between consecutive curves of Figure 2.
const USE_COPPER_STEP: f64 = 0.05;

/// Copper label difference between the base-metal and weld scales of Figure 2:
/// the curve labelled 0.35 % (base metal) is the curve labelled 0.30 % (weld).
const USE_WELD_COPPER_OFFSET: f64 = 0.05;

/// Digitized curves of Figure 2, in the order of [`USE_BASE_COPPER`].
///
/// The two highest-copper curves end on the "UPPER LIMIT" line ([`USE_UPPER_LIMIT`])
/// where they meet it (near 7·10¹⁸ and 3·10¹⁹ n/cm²) and follow it from there on.
const USE_CURVES: [UseLine; 6] = [
    UseLine {
        pct_at_1e18: 10.98,
        slope: 0.2355,
    },
    UseLine {
        pct_at_1e18: 14.02,
        slope: 0.2327,
    },
    UseLine {
        pct_at_1e18: 16.94,
        slope: 0.2327,
    },
    UseLine {
        pct_at_1e18: 20.09,
        slope: 0.2292,
    },
    UseLine {
        pct_at_1e18: 23.13,
        slope: 0.2274,
    },
    UseLine {
        pct_at_1e18: 26.00,
        slope: 0.2260,
    },
];

/// The "UPPER LIMIT" line of Figure 2: no curve rises above it.
const USE_UPPER_LIMIT: UseLine = UseLine {
    pct_at_1e18: 30.18,
    slope: 0.1479,
};

/// Predicted decrease of the Charpy upper-shelf energy, in percent of the
/// unirradiated value (RG 1.99 Rev. 2, Regulatory Position 1.2, Figure 2).
///
/// * The figure holds six curves labelled with the copper content of the
///   material (0.10 … 0.35 wt-% for base metal, 0.05 … 0.30 wt-% for welds); the
///   result is interpolated linearly in copper between the two enclosing
///   curves, as the guide permits.  Copper beyond the plotted curves is clamped
///   to the nearest curve (for low copper this is the conservative choice).
/// * Every curve is capped by the "UPPER LIMIT" line.
/// * Figure 2 plots fluences from 10¹⁸ to 6·10¹⁹ n/cm².  Above 6·10¹⁹ and between
///   [`RG199_MIN_FLUENCE_N_CM2`] and 10¹⁸ n/cm² the straight lines of the figure
///   are continued (an extrapolation); at or below
///   [`RG199_MIN_FLUENCE_N_CM2`] the result is 0.
///
/// The result lies in 0 … 100.
pub fn upper_shelf_energy_decrease_percent(
    form: RpvProductForm,
    copper_wt: f64,
    fluence_n_cm2: f64,
) -> f64 {
    if fluence_n_cm2 <= RG199_MIN_FLUENCE_N_CM2 {
        return 0.0;
    }
    let base_copper = match form {
        RpvProductForm::BaseMetal => copper_wt,
        RpvProductForm::Weld => copper_wt + USE_WELD_COPPER_OFFSET,
    };
    let lowest = USE_BASE_COPPER[0];
    let highest = USE_BASE_COPPER[USE_BASE_COPPER.len() - 1];
    let position = (base_copper.clamp(lowest, highest) - lowest) / USE_COPPER_STEP;
    // Index of the lower enclosing curve (the highest curve closes the last cell).
    let lower = (position.floor() as usize).min(USE_CURVES.len() - 2);
    let fraction = (position - lower as f64).clamp(0.0, 1.0);
    let upper_limit = USE_UPPER_LIMIT.percent_at(fluence_n_cm2);
    let at_curve = |index: usize| {
        USE_CURVES.get(index).map_or(f64::NAN, |line| {
            line.percent_at(fluence_n_cm2).min(upper_limit)
        })
    };
    let percent = at_curve(lower) * (1.0 - fraction) + at_curve(lower + 1) * fraction;
    percent.clamp(0.0, 100.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nuclear_materials_advanced::{RpvCharpy, RpvSteelGrade, chemistry_factor};

    const WELD: RpvProductForm = RpvProductForm::Weld;
    const BASE: RpvProductForm = RpvProductForm::BaseMetal;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    // ── Tables 1 and 2 ──────────────────────────────────────────────────────

    #[test]
    fn tables_have_the_20_f_floor_and_grow_with_copper_and_nickel() {
        for (name, table) in [("weld", &WELD_CF_F), ("base", &BASE_CF_F)] {
            for (i, row) in table.iter().enumerate() {
                assert!(row.iter().all(|&v| v >= 20), "{name} row {i} below 20 °F");
                assert!(
                    row.windows(2).all(|w| w[1] >= w[0]),
                    "{name} row {i} must not decrease with nickel"
                );
            }
            assert!(
                table.windows(2).all(|rows| rows[1]
                    .iter()
                    .zip(rows[0].iter())
                    .all(|(next, previous)| next >= previous)),
                "{name} table must not decrease with copper in any nickel column"
            );
        }
        // Both tables end in the same row (0.40 % Cu).
        assert_eq!(WELD_CF_F[CU_ROWS - 1], BASE_CF_F[CU_ROWS - 1]);
    }

    #[test]
    fn spot_values_match_the_published_weld_table() {
        // (Cu, Ni, CF in °F) read from Table 1 of RG 1.99 Rev. 2.
        let published = [
            (0.00, 0.0, 20.0),
            (0.02, 0.2, 26.0),
            (0.05, 0.4, 67.0),
            (0.10, 0.6, 122.0),
            (0.14, 0.8, 168.0),
            (0.20, 0.8, 194.0),
            (0.25, 1.0, 243.0),
            (0.30, 1.2, 290.0),
            (0.35, 0.6, 212.0),
            (0.40, 1.2, 320.0),
        ];
        for (cu, ni, cf) in published {
            let got = chemistry_factor_fahrenheit(WELD, cu, ni);
            assert!(
                close(got, cf, 1e-9),
                "weld CF({cu}, {ni}) = {got}, table {cf}"
            );
        }
    }

    #[test]
    fn spot_values_match_the_published_base_metal_table() {
        // (Cu, Ni, CF in °F) read from Table 2 of RG 1.99 Rev. 2.
        let published = [
            (0.00, 0.0, 20.0),
            (0.04, 0.2, 26.0),
            (0.05, 0.0, 25.0),
            (0.10, 0.4, 65.0),
            (0.18, 0.8, 141.0),
            (0.20, 0.8, 159.0),
            (0.30, 0.6, 194.0),
            (0.35, 1.2, 298.0),
            (0.40, 1.2, 320.0),
        ];
        for (cu, ni, cf) in published {
            let got = chemistry_factor_fahrenheit(BASE, cu, ni);
            assert!(
                close(got, cf, 1e-9),
                "base CF({cu}, {ni}) = {got}, table {cf}"
            );
        }
    }

    #[test]
    fn old_invented_fit_was_far_below_the_table_for_a_weld() {
        // Cu 0.20 % / Ni 0.80 % weld: Table 1 gives 194 °F = 107.8 °C; the former
        // fit `(Cu - 0.10) * 250 + Ni * 20` returned 41 °C.
        let cf_c = chemistry_factor_for(WELD, 0.20, 0.80);
        assert!(close(cf_c, 194.0 * 5.0 / 9.0, 1e-9), "CF = {cf_c} °C");
        assert!(cf_c > 100.0);
    }

    #[test]
    fn chemistry_factor_interpolates_bilinearly_between_table_entries() {
        // Midway in copper between rows 0.20 and 0.21 at Ni = 0.8.
        let mid_cu = chemistry_factor_fahrenheit(WELD, 0.205, 0.80);
        assert!(close(mid_cu, (194.0 + 197.0) / 2.0, 1e-9), "{mid_cu}");
        // Midway in nickel between 0 and 0.2 at Cu = 0.20.
        let mid_ni = chemistry_factor_fahrenheit(WELD, 0.20, 0.10);
        assert!(close(mid_ni, (88.0 + 104.0) / 2.0, 1e-9), "{mid_ni}");
        // Centre of the cell (Cu 0.20..0.21, Ni 0..0.2): mean of the 4 corners.
        let centre = chemistry_factor_fahrenheit(WELD, 0.205, 0.10);
        assert!(
            close(centre, (88.0 + 104.0 + 92.0 + 108.0) / 4.0, 1e-9),
            "{centre}"
        );
        // Cu 0.2025 (a quarter of the way from row 0.20 to row 0.21) and Ni 0.15
        // (three quarters of the way from column 0 to column 0.2), base metal:
        // row 0.20 -> (82, 102), row 0.21 -> (86, 107).
        let quarter = chemistry_factor_fahrenheit(BASE, 0.2025, 0.15);
        let lower_row = 82.0 + 0.75 * (102.0 - 82.0);
        let upper_row = 86.0 + 0.75 * (107.0 - 86.0);
        let by_hand = lower_row + 0.25 * (upper_row - lower_row);
        assert!(close(by_hand, 98.1875, 1e-12));
        assert!(close(quarter, by_hand, 1e-9), "{quarter} vs {by_hand}");
        // Table values themselves are reproduced exactly on the grid.
        assert!(close(
            chemistry_factor_fahrenheit(WELD, 0.29, 0.6),
            191.0,
            1e-9
        ));
    }

    #[test]
    fn chemistry_factor_is_continuous_across_the_grid_lines() {
        let eps = 1e-9;
        for &cu in &[0.03, 0.10, 0.20, 0.29, 0.39] {
            for &ni in &[0.2, 0.6, 1.0] {
                for form in [WELD, BASE] {
                    let below = chemistry_factor_fahrenheit(form, cu - eps, ni - eps);
                    let above = chemistry_factor_fahrenheit(form, cu + eps, ni + eps);
                    assert!(close(below, above, 1e-4), "{form:?} Cu {cu} Ni {ni}");
                }
            }
        }
    }

    #[test]
    fn chemistry_factor_clamps_outside_the_table_and_propagates_nan() {
        assert_eq!(
            chemistry_factor_fahrenheit(WELD, 0.9, 3.0),
            chemistry_factor_fahrenheit(WELD, 0.40, 1.20)
        );
        assert_eq!(
            chemistry_factor_fahrenheit(BASE, -0.1, -1.0),
            chemistry_factor_fahrenheit(BASE, 0.0, 0.0)
        );
        assert!(close(
            chemistry_factor_fahrenheit(BASE, f64::INFINITY, 0.0),
            175.0,
            1e-9
        ));
        assert!(chemistry_factor_fahrenheit(WELD, f64::NAN, 0.5).is_nan());
        assert!(chemistry_factor_fahrenheit(WELD, 0.2, f64::NAN).is_nan());
    }

    #[test]
    fn chemistry_factor_in_celsius_is_the_difference_conversion() {
        for (cu, ni) in [(0.0, 0.0), (0.2, 0.8), (0.35, 1.0)] {
            for form in [WELD, BASE] {
                let f = chemistry_factor_fahrenheit(form, cu, ni);
                assert!(close(
                    chemistry_factor_for(form, cu, ni),
                    f * 5.0 / 9.0,
                    1e-12
                ));
            }
        }
        // The 20 °F floor is 11.1 °C: there is never a zero chemistry factor.
        assert!(close(chemistry_factor(0.0, 0.0), 20.0 * 5.0 / 9.0, 1e-12));
        assert!(close(chemistry_factor(0.05, 0.0), 25.0 * 5.0 / 9.0, 1e-12));
    }

    #[test]
    fn weld_and_base_metal_tables_differ_where_the_guide_says_so() {
        // Weld metal is more sensitive to copper at moderate copper content.
        assert!(chemistry_factor_for(WELD, 0.20, 0.8) > chemistry_factor_for(BASE, 0.20, 0.8));
        assert!(chemistry_factor_for(WELD, 0.05, 0.4) > chemistry_factor_for(BASE, 0.05, 0.4));
    }

    // ── Fluence factor (Equation 2 / Figure 1) ──────────────────────────────

    #[test]
    fn fluence_factor_matches_figure_1() {
        // Figure 1: 0.42 at 1e18, 1.0 at 1e19, 1.5 at 1e20 n/cm².
        assert!(close(rg199_fluence_factor(1.0e19), 1.0, 1e-12));
        assert!(close(
            rg199_fluence_factor(1.0e18),
            0.1_f64.powf(0.38),
            1e-12
        ));
        assert!(close(rg199_fluence_factor(1.0e18), 0.417, 0.005));
        assert!(close(
            rg199_fluence_factor(1.0e20),
            10.0_f64.powf(0.18),
            1e-12
        ));
        assert!(close(rg199_fluence_factor(1.0e20), 1.514, 0.005));
        // f = 0.01 (1e17 n/cm²): 0.01^0.48 = 0.1096, the lowest point of Figure 1.
        assert!(close(rg199_fluence_factor(1.000_000_1e17), 0.1096, 1e-3));
    }

    #[test]
    fn fluence_factor_is_zero_at_or_below_the_range_and_monotone_inside_it() {
        assert_eq!(rg199_fluence_factor(0.0), 0.0);
        assert_eq!(rg199_fluence_factor(-5.0), 0.0);
        assert_eq!(rg199_fluence_factor(RG199_MIN_FLUENCE_N_CM2), 0.0);
        let mut previous = 0.0;
        let mut fluence = 1.1e17;
        while fluence <= 1.0e20 {
            let factor = rg199_fluence_factor(fluence);
            assert!(factor > previous, "fluence factor must rise at {fluence}");
            previous = factor;
            fluence *= 1.07;
        }
    }

    // ── Figure 2 (upper-shelf energy decrease) ──────────────────────────────

    #[test]
    fn use_curves_start_at_the_values_printed_at_1e18() {
        // Values read at the left edge of the plotted curves (10^18 n/cm²).
        let figure = [
            (0.10, 11.0),
            (0.15, 14.0),
            (0.20, 17.0),
            (0.25, 20.0),
            (0.30, 23.0),
            (0.35, 26.0),
        ];
        for (copper, percent) in figure {
            let got = upper_shelf_energy_decrease_percent(BASE, copper, 1.0e18 * 1.000_000_1);
            assert!(
                close(got, percent, 0.5),
                "base Cu {copper}: {got} vs {percent}"
            );
        }
        assert!(close(USE_UPPER_LIMIT.percent_at(1.0e18), 30.0, 0.5));
    }

    #[test]
    fn use_curves_end_at_the_values_printed_at_6e19() {
        // Right edge of the figure (6e19 n/cm²), read off the printed curves.
        let figure = [
            (0.10, 28.8),
            (0.15, 36.4),
            (0.20, 44.0),
            (0.25, 51.3),
            (0.30, 55.3), // merged with the upper limit
            (0.35, 55.3), // merged with the upper limit
        ];
        for (copper, percent) in figure {
            let got = upper_shelf_energy_decrease_percent(BASE, copper, 6.0e19);
            assert!(
                close(got, percent, 1.0),
                "base Cu {copper}: {got} vs {percent}"
            );
        }
        assert!(close(USE_UPPER_LIMIT.percent_at(6.0e19), 55.3, 1.0));
    }

    #[test]
    fn use_curves_at_the_1e19_grid_line() {
        let figure = [(0.10, 18.9), (0.15, 24.0), (0.20, 28.9), (0.25, 34.1)];
        for (copper, percent) in figure {
            let got = upper_shelf_energy_decrease_percent(BASE, copper, 1.0e19);
            assert!(
                close(got, percent, 0.5),
                "base Cu {copper}: {got} vs {percent}"
            );
        }
    }

    #[test]
    fn weld_labels_are_the_base_metal_labels_minus_five_hundredths() {
        for &fluence in &[2.0e18, 1.0e19, 4.0e19] {
            for &weld_cu in &[0.05, 0.10, 0.17, 0.25, 0.30] {
                let weld = upper_shelf_energy_decrease_percent(WELD, weld_cu, fluence);
                let base = upper_shelf_energy_decrease_percent(BASE, weld_cu + 0.05, fluence);
                assert!(close(weld, base, 1e-12), "Cu {weld_cu}, f {fluence}");
            }
        }
    }

    #[test]
    fn use_decrease_grows_with_copper_and_fluence_and_stays_below_the_upper_limit() {
        let mut fluence = 2.0e17;
        while fluence < 1.0e20 {
            let limit = USE_UPPER_LIMIT.percent_at(fluence);
            let mut previous = 0.0;
            let mut copper = 0.0;
            while copper <= 0.45 {
                let percent = upper_shelf_energy_decrease_percent(BASE, copper, fluence);
                assert!(percent >= previous - 1e-12, "monotone in Cu at {fluence}");
                assert!(
                    percent <= limit + 1e-12,
                    "{percent} above upper limit {limit}"
                );
                previous = percent;
                copper += 0.0125;
            }
            fluence *= 1.3;
        }
        for &copper in &[0.0, 0.12, 0.2, 0.33, 0.5] {
            let mut previous = 0.0;
            let mut fluence = 1.5e17;
            while fluence < 1.0e20 {
                let percent = upper_shelf_energy_decrease_percent(BASE, copper, fluence);
                assert!(
                    percent >= previous - 1e-12,
                    "monotone in fluence, Cu {copper}"
                );
                previous = percent;
                fluence *= 1.1;
            }
        }
    }

    #[test]
    fn use_decrease_interpolates_linearly_in_copper_between_curves() {
        let fluence = 1.0e19;
        let low = upper_shelf_energy_decrease_percent(BASE, 0.15, fluence);
        let high = upper_shelf_energy_decrease_percent(BASE, 0.20, fluence);
        let mid = upper_shelf_energy_decrease_percent(BASE, 0.175, fluence);
        assert!(close(mid, 0.5 * (low + high), 1e-12));
        let quarter = upper_shelf_energy_decrease_percent(BASE, 0.1625, fluence);
        assert!(close(quarter, 0.75 * low + 0.25 * high, 1e-12));
    }

    #[test]
    fn upper_limit_caps_the_two_highest_copper_curves() {
        // Curve "0.35 % base / 0.30 % weld" meets the upper limit near 7e18 n/cm².
        let after = upper_shelf_energy_decrease_percent(BASE, 0.35, 2.0e19);
        assert!(close(after, USE_UPPER_LIMIT.percent_at(2.0e19), 1e-12));
        let before = upper_shelf_energy_decrease_percent(BASE, 0.35, 3.0e18);
        assert!(before < USE_UPPER_LIMIT.percent_at(3.0e18) - 0.5);
        // Copper beyond the plotted curves is clamped to the highest curve.
        let beyond = upper_shelf_energy_decrease_percent(BASE, 0.60, 2.0e19);
        assert!(close(beyond, after, 1e-12));
    }

    #[test]
    fn use_decrease_for_low_copper_is_clamped_to_the_lowest_curve() {
        let lowest = upper_shelf_energy_decrease_percent(BASE, 0.10, 5.0e18);
        assert!(close(
            upper_shelf_energy_decrease_percent(BASE, 0.0, 5.0e18),
            lowest,
            1e-12
        ));
        assert!(close(
            upper_shelf_energy_decrease_percent(BASE, -1.0, 5.0e18),
            lowest,
            1e-12
        ));
        let weld_lowest = upper_shelf_energy_decrease_percent(WELD, 0.05, 5.0e18);
        assert!(close(
            upper_shelf_energy_decrease_percent(WELD, 0.0, 5.0e18),
            weld_lowest,
            1e-12
        ));
    }

    #[test]
    fn use_decrease_has_no_jump_at_1e18_and_vanishes_at_the_range_limit() {
        for &copper in &[0.1, 0.2, 0.3] {
            let below = upper_shelf_energy_decrease_percent(BASE, copper, 1.0e18 * (1.0 - 1e-9));
            let above = upper_shelf_energy_decrease_percent(BASE, copper, 1.0e18 * (1.0 + 1e-9));
            assert!(close(below, above, 1e-6), "Cu {copper}: {below} vs {above}");
            assert!(below > 0.0);
            assert_eq!(
                upper_shelf_energy_decrease_percent(BASE, copper, RG199_MIN_FLUENCE_N_CM2),
                0.0
            );
            assert_eq!(upper_shelf_energy_decrease_percent(BASE, copper, 0.0), 0.0);
        }
    }

    #[test]
    fn use_decrease_never_exceeds_100_percent() {
        let huge = upper_shelf_energy_decrease_percent(BASE, 0.35, 1.0e40);
        assert!(close(huge, 100.0, 1e-12));
    }

    // ── RpvCharpy ───────────────────────────────────────────────────────────

    fn charpy(fluence: f64, copper: f64, nickel: f64) -> RpvCharpy {
        RpvCharpy::new(
            RpvSteelGrade::A508Class3,
            120.0,
            -20.0,
            fluence,
            copper,
            nickel,
            0.01,
        )
    }

    #[test]
    fn rtndt_shift_is_the_chemistry_factor_times_the_fluence_factor() {
        // At 1e19 n/cm² the fluence factor is 1: the shift equals CF.
        let weld = charpy(1.0e19, 0.20, 0.80).with_product_form(WELD);
        assert!(close(weld.rtndt_shift(), 194.0 * 5.0 / 9.0, 1e-9));
        let base = charpy(1.0e19, 0.20, 0.80);
        assert_eq!(base.product_form, BASE);
        assert!(close(base.rtndt_shift(), 159.0 * 5.0 / 9.0, 1e-9));
        // At 1e18 n/cm² the factor is 0.1^0.38.
        let low = charpy(1.0e18, 0.20, 0.80);
        assert!(close(
            low.rtndt_shift(),
            159.0 * 5.0 / 9.0 * 0.1_f64.powf(0.38),
            1e-9
        ));
    }

    #[test]
    fn rtndt_shift_has_no_cutoff_at_1e18() {
        // Figure 1 reaches down to 1e17 n/cm²; the old 1e18 cutoff removed 40 % of CF.
        let just_below = charpy(0.999e18, 0.20, 0.80).rtndt_shift();
        let just_above = charpy(1.001e18, 0.20, 0.80).rtndt_shift();
        assert!(
            just_below > 0.0 && close(just_below, just_above, 0.2),
            "{just_below} vs {just_above}"
        );
        assert_eq!(charpy(1.0e17, 0.20, 0.80).rtndt_shift(), 0.0);
        assert!(charpy(2.0e17, 0.20, 0.80).rtndt_shift() > 0.0);
    }

    #[test]
    fn rtndt_shift_of_a_low_copper_steel_is_not_zero() {
        // 20 °F floor: even copper-free steel embrittles.
        let shift = charpy(1.0e19, 0.0, 0.0).rtndt_shift();
        assert!(close(shift, 20.0 * 5.0 / 9.0, 1e-9), "{shift}");
    }

    #[test]
    fn use_drop_applies_the_figure_2_percentage_to_the_unirradiated_energy() {
        // Base metal, Cu 0.20 %, 1e19 n/cm²: about 29 % of 120 J.
        let record = charpy(1.0e19, 0.20, 0.7);
        let expected = 120.0 * upper_shelf_energy_decrease_percent(BASE, 0.20, 1.0e19) / 100.0;
        assert!(close(record.use_drop(), expected, 1e-12));
        assert!(
            close(record.use_drop(), 120.0 * 0.289, 0.6),
            "{}",
            record.use_drop()
        );
        assert!(close(
            record.use_irradiated(),
            120.0 - record.use_drop(),
            1e-12
        ));
        // The weld curve with the same label (Cu 0.15 %) gives the same decrease.
        let weld = charpy(1.0e19, 0.15, 0.7).with_product_form(WELD);
        assert!(close(weld.use_drop(), record.use_drop(), 1e-12));
    }

    #[test]
    fn use_drop_is_independent_of_nickel_and_phosphorus() {
        let a = RpvCharpy::new(
            RpvSteelGrade::A533GradeB,
            100.0,
            0.0,
            5.0e18,
            0.2,
            0.2,
            0.001,
        );
        let b = RpvCharpy::new(
            RpvSteelGrade::A533GradeB,
            100.0,
            0.0,
            5.0e18,
            0.2,
            1.0,
            0.03,
        );
        assert_eq!(a.use_drop(), b.use_drop());
    }
}
