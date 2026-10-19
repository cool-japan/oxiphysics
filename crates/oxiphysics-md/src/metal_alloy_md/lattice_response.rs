// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lattice response of a cubic EAM crystal: static pressure, elastic constants
//! and quasi-harmonic thermal expansion.
//!
//! Everything is derived from the energy of a periodic [`AlloySystem`]
//! supercell of the element's own crystal structure (FCC or BCC, see
//! [`EamParams::structure`]) that is deformed homogeneously.  A shear
//! deforms the periodic box together with the atoms (a tilted, monoclinic
//! [`SimBox`]), so the periodic images follow the strain instead of leaving
//! a mismatch across the boundary.
//!
//! - [`lattice_pressure_eam`] — static pressure at the tabulated lattice
//!   constant, the derivative of the energy with respect to the volume.
//! - [`elastic_constants_eam`] — `C11`, `C12`, `C44` from central second
//!   differences of the strain energy.
//! - [`quasi_harmonic_lattice_constant`] and [`thermal_expansion_coefficient`]
//!   — lattice constant and linear expansion from the minimum of the Helmholtz
//!   free energy `F(a, T) = E_static(a) / N + F_vib(a, T)` of an Einstein
//!   crystal whose frequency follows the lattice constant.

use super::functions::{EV_PER_ANGSTROM_AMU_TO_ANGSTROM_PER_PS2, generate_crystal_positions};
use super::types::{AlloySystem, CrystalStructure, EamParams, MixingRule, SimBox};
use crate::error::MdError;

#[cfg(test)]
mod tests;

/// Conversion from eV/Å³ to GPa (CODATA 2018: 1 eV/Å³ = 160.2176634 GPa).
const EV_PER_ANGSTROM3_TO_GPA: f64 = 160.217_663_4;
/// Boltzmann constant in eV/K (CODATA 2018).
const KB_EV_PER_K: f64 = 8.617_333_262e-5;
/// Reduced Planck constant in eV·ps (CODATA 2018).
const HBAR_EV_PS: f64 = 6.582_119_569e-4;
/// Strain amplitude `epsilon` of the central differences of the elastic constants.
const ELASTIC_STRAIN: f64 = 1.0e-3;
/// Largest relative disagreement between the elastic constants obtained with
/// the strain amplitudes `epsilon` and `2 epsilon` that is accepted as smooth.
const ELASTIC_SMOOTHNESS_TOLERANCE: f64 = 1.0e-2;
/// Floor (eV/Å³) of the scale against which the smoothness check is made, so
/// that constants which vanish by symmetry do not trigger it.
const ELASTIC_SMOOTHNESS_SCALE_FLOOR: f64 = 1.0e-3;
/// Relative lattice-constant step of the five-point central difference of the
/// pressure.
const PRESSURE_STEP: f64 = 5.0e-4;
/// Displacement (Å) of the single atom of the Einstein force constant.
const EINSTEIN_DISPLACEMENT: f64 = 1.0e-3;
/// Largest relative deviation from `lattice_a` that the lattice-constant search
/// of the quasi-harmonic model explores.
const MAX_LATTICE_SHIFT: f64 = 0.03;
/// Relative spacing of the lattice constants scanned before the refinement.
const SCAN_STEP: f64 = 0.0025;
/// Relative width at which the golden-section refinement stops.
const REFINEMENT_TOLERANCE: f64 = 1.0e-9;
/// Upper bound of the iterations of the golden-section refinement.
const MAX_REFINEMENT_ITERATIONS: usize = 200;
/// Largest supercell (conventional cells per edge) that is accepted; the
/// neighbour search is `O(N²)` in the number of atoms.
const MAX_SUPERCELL_CELLS: usize = 12;

/// Upper-triangular deformation gradient
/// `F = [[xx, xy, 0], [0, yy, 0], [0, 0, zz]]` of the periodic supercell.
///
/// The box vectors become `a = (xx L, 0, 0)`, `b = (xy L, yy L, 0)` and
/// `c = (0, 0, zz L)`, and every atom moves to `F r`.  A pure `xy` is the simple
/// shear `x -> x + xy y`.
#[derive(Debug, Clone, Copy)]
struct Deformation {
    xx: f64,
    yy: f64,
    zz: f64,
    xy: f64,
}

impl Deformation {
    /// No deformation.
    const IDENTITY: Self = Self {
        xx: 1.0,
        yy: 1.0,
        zz: 1.0,
        xy: 0.0,
    };

    /// Uniaxial strain `epsilon` along x.
    fn uniaxial(epsilon: f64) -> Self {
        Self {
            xx: 1.0 + epsilon,
            ..Self::IDENTITY
        }
    }

    /// Equibiaxial strain `epsilon` along x and y.
    fn biaxial(epsilon: f64) -> Self {
        Self {
            xx: 1.0 + epsilon,
            yy: 1.0 + epsilon,
            ..Self::IDENTITY
        }
    }

    /// Simple shear `x -> x + gamma y` (engineering shear strain `gamma`).
    fn simple_shear(gamma: f64) -> Self {
        Self {
            xy: gamma,
            ..Self::IDENTITY
        }
    }
}

/// Number of atoms of the conventional cubic cell, or an error for a structure
/// without cubic symmetry.
fn atoms_per_cell(structure: CrystalStructure) -> Result<usize, MdError> {
    match structure {
        CrystalStructure::FCC => Ok(4),
        CrystalStructure::BCC => Ok(2),
        CrystalStructure::HCP => Err(MdError::General(
            "the cubic lattice response is defined for FCC and BCC crystals only, \
             not for HCP"
                .to_string(),
        )),
    }
}

/// Check that the parameters describe a usable cubic crystal and return the
/// number of atoms of its conventional cell.
fn validate_params(params: &EamParams) -> Result<usize, MdError> {
    let atoms = atoms_per_cell(params.structure)?;
    for (name, value) in [
        ("lattice_a", params.lattice_a),
        ("cutoff", params.cutoff),
        ("re", params.re),
    ] {
        if !(value.is_finite() && value > 0.0) {
            return Err(MdError::General(format!(
                "EAM parameter `{name}` must be positive and finite, got {value}"
            )));
        }
    }
    Ok(atoms)
}

/// Smallest number of conventional cells per supercell edge for which the
/// minimum-image convention of the periodic supercell is exact.
///
/// The brute-force neighbour search of [`AlloySystem`] only sees pairs closer
/// than half of the box edge, so `n a` has to exceed twice the cutoff.  The
/// lattice-constant search of the quasi-harmonic model shrinks the lattice by
/// up to 3 %, which is included.
///
/// # Errors
/// [`MdError::General`] for an HCP crystal, non-positive or non-finite
/// `lattice_a`, `cutoff` or `re`, or a cutoff so long that the supercell would
/// exceed 12 conventional cells per edge.
pub fn minimum_supercell_cells(params: &EamParams) -> Result<usize, MdError> {
    validate_params(params)?;
    let shrunk_edge = params.lattice_a * (1.0 - MAX_LATTICE_SHIFT);
    let cells = (2.0 * params.cutoff / shrunk_edge).floor() + 1.0;
    if cells > MAX_SUPERCELL_CELLS as f64 {
        return Err(MdError::General(format!(
            "cutoff {} Å needs a supercell of {cells} cells per edge at lattice constant \
             {} Å, more than the supported {MAX_SUPERCELL_CELLS}",
            params.cutoff, params.lattice_a
        )));
    }
    Ok(cells as usize)
}

/// Check that `n_cells` is a supercell edge for which the minimum-image
/// convention holds over the whole lattice-constant search range.
fn check_supercell(params: &EamParams, n_cells: usize) -> Result<(), MdError> {
    let minimum = minimum_supercell_cells(params)?;
    if n_cells < minimum {
        return Err(MdError::General(format!(
            "n_cells = {n_cells} is too small: the periodic supercell edge n_cells * a must \
             exceed twice the cutoff ({} Å), which needs at least {minimum} cells for \
             a = {} Å",
            params.cutoff, params.lattice_a
        )));
    }
    if n_cells > MAX_SUPERCELL_CELLS {
        return Err(MdError::General(format!(
            "n_cells = {n_cells} exceeds the supported maximum of {MAX_SUPERCELL_CELLS}"
        )));
    }
    Ok(())
}

/// Periodic supercell of `n_cells`³ conventional cells with lattice constant
/// `lattice_a`, deformed by `deformation`.
fn crystal_system(
    params: &EamParams,
    lattice_a: f64,
    n_cells: usize,
    deformation: Deformation,
) -> AlloySystem {
    let edge = n_cells as f64 * lattice_a;
    let sim_box = SimBox::new(
        edge * deformation.xx,
        edge * deformation.yy,
        edge * deformation.zz,
    )
    .with_xy_tilt(edge * deformation.xy);
    let mut system = AlloySystem::new(vec![params.clone()], sim_box, MixingRule::Geometric, 0.001);
    for position in generate_crystal_positions(lattice_a, n_cells, params.structure) {
        let deformed = [
            deformation.xx * position[0] + deformation.xy * position[1],
            deformation.yy * position[1],
            deformation.zz * position[2],
        ];
        system.add_atom(system.sim_box.wrap(deformed), 0);
    }
    system
}

/// Potential energy per atom (eV) of a non-empty system.
fn energy_per_atom(system: &AlloySystem) -> f64 {
    system.potential_energy() / system.num_atoms().max(1) as f64
}

/// Static pressure (GPa) of the perfect crystal at its tabulated lattice
/// constant `lattice_a`: `P = -dE/dV` of the energy per atom with respect to the
/// atomic volume, by a five-point central difference (error `O(h^4)`) in the
/// lattice constant.  Positive values mean the crystal is compressed (it would
/// expand if released).
///
/// A pressure that is not close to zero means `lattice_a` is not the
/// equilibrium lattice constant of the potential; the elastic constants of
/// [`elastic_constants_eam`] then belong to a pre-stressed state and the
/// quasi-harmonic expansion of [`thermal_expansion_coefficient`] is not defined
/// (see [`EamParams::calibrated_to_zero_pressure`]).
///
/// # Errors
/// See [`minimum_supercell_cells`].
pub fn lattice_pressure_eam(params: &EamParams) -> Result<f64, MdError> {
    let atoms = validate_params(params)?;
    let n_cells = minimum_supercell_cells(params)?;
    let a = params.lattice_a;
    let step = PRESSURE_STEP * a;
    let energy_at = |lattice_a: f64| {
        energy_per_atom(&crystal_system(
            params,
            lattice_a,
            n_cells,
            Deformation::IDENTITY,
        ))
    };
    let de_da = (8.0 * (energy_at(a + step) - energy_at(a - step))
        - (energy_at(a + 2.0 * step) - energy_at(a - 2.0 * step)))
        / (12.0 * step);
    let dv_da = 3.0 * a * a / atoms as f64;
    Ok(-de_da / dv_da * EV_PER_ANGSTROM3_TO_GPA)
}

impl EamParams {
    /// A copy of these parameters whose pair-repulsion amplitude `a_coeff` is
    /// chosen so that [`EamParams::lattice_a`] is the zero-pressure lattice
    /// constant of the potential ([`lattice_pressure_eam`] returns zero).
    ///
    /// The static energy is exactly linear in `a_coeff`, hence so is the pressure,
    /// and two pressure evaluations give the root.  The tabulated presets
    /// ([`EamParams::copper`], ...) are not self-consistent: at their `lattice_a`
    /// they are under a tensile pressure of 29 GPa (aluminium) to 124 GPa
    /// (nickel), their strain-energy curvatures `C11` and `(C11 + 2 C12) / 3` are
    /// negative, and they have no quasi-harmonic equilibrium.  The calibrated copy
    /// is a stable crystal at `lattice_a`, for which [`elastic_constants_eam`]
    /// returns the physical elastic constants and
    /// [`thermal_expansion_coefficient`] is defined.
    ///
    /// # Errors
    /// See [`minimum_supercell_cells`]; also [`MdError::General`] if the pressure
    /// does not depend on `a_coeff` or the zero-pressure amplitude is not
    /// positive.
    pub fn calibrated_to_zero_pressure(&self) -> Result<Self, MdError> {
        let pressure_for = |a_coeff: f64| {
            let mut trial = self.clone();
            trial.a_coeff = a_coeff;
            lattice_pressure_eam(&trial)
        };
        let reference = pressure_for(self.a_coeff)?;
        let slope = pressure_for(self.a_coeff + 1.0)? - reference;
        if !(slope.is_finite() && slope > 0.0) {
            return Err(MdError::General(format!(
                "the static pressure does not increase with the repulsion amplitude \
                 (slope {slope} GPa per eV)"
            )));
        }
        let a_coeff = self.a_coeff - reference / slope;
        if !(a_coeff.is_finite() && a_coeff > 0.0) {
            return Err(MdError::General(format!(
                "no positive repulsion amplitude gives zero pressure at lattice_a = {} Å \
                 (solution {a_coeff} eV)",
                self.lattice_a
            )));
        }
        let mut calibrated = self.clone();
        calibrated.a_coeff = a_coeff;
        Ok(calibrated)
    }
}

/// Elastic constants `(C11, C12, C44)` in GPa of a cubic EAM crystal.
///
/// The crystal has the element's own structure ([`EamParams::structure`]: FCC
/// or BCC) and lattice constant [`EamParams::lattice_a`]; the periodic
/// supercell is the smallest one that satisfies the minimum-image convention
/// ([`minimum_supercell_cells`]).  With `V0` the volume of the supercell and
/// `epsilon = 1e-3`, the constants are central second differences of the
/// strain energy (`E+` and `E-` are the energies for `+epsilon` and
/// `-epsilon`, `E0` is the undeformed energy):
///
/// - uniaxial strain `e_xx = epsilon`: `E+ + E- - 2 E0 = V0 C11 epsilon²`,
/// - equibiaxial strain `e_xx = e_yy = epsilon`:
///   `E+ + E- - 2 E0 = 2 V0 (C11 + C12) epsilon²`, so `C12` is that value
///   *minus* `C11`,
/// - simple shear `x -> x + gamma y`, `gamma = 2 epsilon` (engineering shear
///   strain): `E+ + E- - 2 E0 = V0 C44 gamma²`.  The periodic box is sheared
///   together with the atoms (a tilted [`SimBox`]).
///
/// The result is `C_ij = (1/V0) d²E / d e_i d e_j` at the lattice constant
/// `lattice_a`.  These coincide with the physical elastic constants only when
/// `lattice_a` is the zero-pressure lattice constant of the potential; check
/// [`lattice_pressure_eam`].  Under a residual pressure `P` the bulk modulus
/// of the strain energy obeys `(C11 + 2 C12) / 3 = V d²E/dV² - 2 P / 3`.
///
/// The strain amplitudes `epsilon` and `2 epsilon` are both evaluated; if they
/// disagree by more than 1 % of the largest constant the energy is not smooth
/// under strain (a neighbour shell sits on the cutoff radius), and an error is
/// returned instead of a meaningless number.
///
/// # Errors
/// [`MdError::General`] for an HCP crystal (the cubic constants are not
/// defined), invalid parameters (see [`minimum_supercell_cells`]), or an energy
/// that is not smooth under strain.
pub fn elastic_constants_eam(params: &EamParams) -> Result<(f64, f64, f64), MdError> {
    validate_params(params)?;
    let n_cells = minimum_supercell_cells(params)?;
    let fine = elastic_constants_at_strain(params, n_cells, ELASTIC_STRAIN);
    let coarse = elastic_constants_at_strain(params, n_cells, 2.0 * ELASTIC_STRAIN);
    let scale = fine
        .0
        .abs()
        .max(fine.1.abs())
        .max(fine.2.abs())
        .max(ELASTIC_SMOOTHNESS_SCALE_FLOOR);
    let worst = (fine.0 - coarse.0)
        .abs()
        .max((fine.1 - coarse.1).abs())
        .max((fine.2 - coarse.2).abs());
    let smooth = worst.is_finite() && worst <= ELASTIC_SMOOTHNESS_TOLERANCE * scale;
    if !smooth {
        return Err(MdError::General(format!(
            "the EAM energy is not smooth under a strain of {ELASTIC_STRAIN}: the elastic \
             constants change by {:.3} GPa between strain amplitudes {ELASTIC_STRAIN} and \
             {}; a neighbour shell probably sits on the cutoff radius {} Å",
            worst * EV_PER_ANGSTROM3_TO_GPA,
            2.0 * ELASTIC_STRAIN,
            params.cutoff
        )));
    }
    Ok((
        fine.0 * EV_PER_ANGSTROM3_TO_GPA,
        fine.1 * EV_PER_ANGSTROM3_TO_GPA,
        fine.2 * EV_PER_ANGSTROM3_TO_GPA,
    ))
}

/// `(C11, C12, C44)` in eV/Å³ from the central differences with strain
/// amplitude `epsilon` (see [`elastic_constants_eam`]).
fn elastic_constants_at_strain(
    params: &EamParams,
    n_cells: usize,
    epsilon: f64,
) -> (f64, f64, f64) {
    let a = params.lattice_a;
    let energy = |deformation: Deformation| {
        crystal_system(params, a, n_cells, deformation).potential_energy()
    };
    let volume = (n_cells as f64 * a).powi(3);
    let reference = energy(Deformation::IDENTITY);
    let second_difference =
        |plus: Deformation, minus: Deformation| energy(plus) + energy(minus) - 2.0 * reference;
    let c11 = second_difference(
        Deformation::uniaxial(epsilon),
        Deformation::uniaxial(-epsilon),
    ) / (volume * epsilon * epsilon);
    let c11_plus_c12 = second_difference(
        Deformation::biaxial(epsilon),
        Deformation::biaxial(-epsilon),
    ) / (2.0 * volume * epsilon * epsilon);
    let gamma = 2.0 * epsilon;
    let c44 = second_difference(
        Deformation::simple_shear(gamma),
        Deformation::simple_shear(-gamma),
    ) / (volume * gamma * gamma);
    (c11, c11_plus_c12 - c11, c44)
}

/// Einstein force constant `k = d²E / du²` (eV/Å²) of one atom of the perfect
/// crystal of lattice constant `lattice_a`, displaced by `u` along x while all
/// other atoms stay on their sites.
///
/// It is the central difference of the analytic force on the displaced atom,
/// `k = -(F(+u) - F(-u)) / (2 u)`, and contains the pair terms as well as the
/// embedding terms of the displaced atom and of all its neighbours.  By cubic
/// site symmetry it is the same along y and z.  Differencing the force (a first
/// derivative) instead of the energy keeps the round-off at about 1e-13
/// relative; an energy second difference at this displacement is good to only
/// about 1e-9, which would make the minimum of the free energy noisy at the
/// 1e-5 Å level.
fn einstein_force_constant(
    params: &EamParams,
    lattice_a: f64,
    n_cells: usize,
) -> Result<f64, MdError> {
    let mut system = crystal_system(params, lattice_a, n_cells, Deformation::IDENTITY);
    let origin_x = system
        .atoms
        .first()
        .map(|atom| atom.position[0])
        .ok_or_else(|| MdError::General("the supercell contains no atoms".to_string()))?;
    let mut force_at = |shift: f64| -> Result<f64, MdError> {
        if let Some(atom) = system.atoms.first_mut() {
            atom.position[0] = origin_x + shift;
        }
        system.compute_forces();
        system
            .atoms
            .first()
            .map(|atom| atom.force[0])
            .ok_or_else(|| MdError::General("the supercell contains no atoms".to_string()))
    };
    let plus = force_at(EINSTEIN_DISPLACEMENT)?;
    let minus = force_at(-EINSTEIN_DISPLACEMENT)?;
    Ok(-(plus - minus) / (2.0 * EINSTEIN_DISPLACEMENT))
}

/// Einstein angular frequency (rad/ps) of an atom of mass `mass` (amu) with
/// force constant `force_constant` (eV/Å²), or `None` if the force constant is
/// not positive (a dynamically unstable lattice).
///
/// `omega² = k / m` with the conversion `1 eV / (Å amu) = 9648.533 Å/ps²`.
fn einstein_frequency(force_constant: f64, mass: f64) -> Option<f64> {
    (force_constant.is_finite() && force_constant > 0.0)
        .then(|| (force_constant * EV_PER_ANGSTROM_AMU_TO_ANGSTROM_PER_PS2 / mass).sqrt())
}

/// Helmholtz free energy per atom (eV) of three independent quantum harmonic
/// oscillators of angular frequency `omega` (rad/ps) at `temperature` (K):
///
/// ```text
/// F_vib = 3 [ hbar omega / 2 + k_B T ln(1 - exp(-hbar omega / k_B T)) ]
/// ```
///
/// including the zero-point energy; for `hbar omega << k_B T` it reduces to the
/// classical `3 k_B T ln(hbar omega / k_B T)`.
fn vibrational_free_energy(omega: f64, temperature: f64) -> f64 {
    let kt = KB_EV_PER_K * temperature;
    let x = HBAR_EV_PS * omega / kt;
    3.0 * kt * (0.5 * x + (-(-x).exp_m1()).ln())
}

/// Free energy per atom (eV) at lattice constant `lattice_a` and `temperature`,
/// `E_static(a) / N + F_vib(omega_E(a), T)`; `+inf` where the Einstein force
/// constant is not positive.
fn free_energy_per_atom(
    params: &EamParams,
    lattice_a: f64,
    n_cells: usize,
    temperature: f64,
) -> f64 {
    let static_energy = energy_per_atom(&crystal_system(
        params,
        lattice_a,
        n_cells,
        Deformation::IDENTITY,
    ));
    let omega = einstein_force_constant(params, lattice_a, n_cells)
        .ok()
        .and_then(|k| einstein_frequency(k, params.mass));
    match omega {
        Some(omega) => static_energy + vibrational_free_energy(omega, temperature),
        None => f64::INFINITY,
    }
}

/// Sentence fragment that reports the static pressure at `lattice_a` for error
/// messages.
fn static_pressure_note(params: &EamParams) -> String {
    match lattice_pressure_eam(params) {
        Ok(pressure) => format!("the static pressure at lattice_a is {pressure:.1} GPa"),
        Err(_) => "the static pressure could not be evaluated".to_string(),
    }
}

/// Grid index of the local minimum of `values` closest to `centre` that has a
/// neighbour on both sides, or `None`.
fn nearest_interior_minimum(values: &[f64], centre: usize) -> Option<usize> {
    let is_minimum = |index: usize| -> bool {
        if index == 0 || index + 1 >= values.len() {
            return false;
        }
        let (below, here, above) = (values[index - 1], values[index], values[index + 1]);
        here.is_finite() && here <= below && here <= above && (here < below || here < above)
    };
    (0..values.len()).find_map(|distance| {
        [centre.checked_sub(distance), Some(centre + distance)]
            .into_iter()
            .flatten()
            .find(|&index| is_minimum(index))
    })
}

/// Golden-section search for the minimum of `f` in `[lo, hi]`.
///
/// Only comparisons of function values are used, so `+inf` values (unstable
/// regions) are handled.  Stops when the interval is narrower than `tolerance`.
fn golden_section_minimum(f: &impl Fn(f64) -> f64, lo: f64, hi: f64, tolerance: f64) -> f64 {
    const INVERSE_GOLDEN_RATIO: f64 = 0.618_033_988_749_894_9;
    let (mut lo, mut hi) = (lo, hi);
    let mut c = hi - INVERSE_GOLDEN_RATIO * (hi - lo);
    let mut d = lo + INVERSE_GOLDEN_RATIO * (hi - lo);
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..MAX_REFINEMENT_ITERATIONS {
        if hi - lo <= tolerance {
            break;
        }
        if fc < fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - INVERSE_GOLDEN_RATIO * (hi - lo);
            fc = f(c);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + INVERSE_GOLDEN_RATIO * (hi - lo);
            fd = f(d);
        }
    }
    0.5 * (lo + hi)
}

/// Lattice constant (Å) that minimises the quasi-harmonic Helmholtz free energy
/// at `temperature` (K).
///
/// The crystal is a periodic supercell of `n_cells`³ conventional cells of the
/// element's structure ([`EamParams::structure`]).  Its free energy per atom is
///
/// ```text
/// F(a, T) = E_static(a) / N + 3 [ hbar omega_E(a) / 2 + k_B T ln(1 - exp(-hbar omega_E(a) / k_B T)) ]
/// ```
///
/// (the classical limit is `3 k_B T ln(hbar omega_E / k_B T)`), where the
/// Einstein frequency follows from the force constant of one atom displaced
/// against its fixed neighbours, `omega_E² = k(a) / m`.  By the identity
/// `<omega²> = k / m` this is the root-mean-square phonon frequency of the
/// crystal.  Anharmonicity beyond the volume dependence of `omega_E` is not
/// included.
///
/// The search is local: the free energy is scanned in steps of 0.25 % of
/// `lattice_a` within +-3 % and the minimum closest to `lattice_a` is refined
/// by a golden-section search.  `lattice_a` therefore has to be near the
/// potential's equilibrium lattice constant (see [`lattice_pressure_eam`]).
///
/// # Errors
/// [`MdError::General`] if the temperature is not positive and finite, the
/// mass is not positive, the structure is HCP, `n_cells` is below
/// [`minimum_supercell_cells`] (or above 12), the Einstein force constant is
/// not positive at `lattice_a`, or the free energy has no minimum within +-3 %
/// of `lattice_a`.
pub fn quasi_harmonic_lattice_constant(
    params: &EamParams,
    temperature: f64,
    n_cells: usize,
) -> Result<f64, MdError> {
    validate_params(params)?;
    if !(params.mass.is_finite() && params.mass > 0.0) {
        return Err(MdError::General(format!(
            "atomic mass must be positive and finite, got {}",
            params.mass
        )));
    }
    if !(temperature.is_finite() && temperature > 0.0) {
        return Err(MdError::General(format!(
            "temperature must be positive and finite, got {temperature} K"
        )));
    }
    check_supercell(params, n_cells)?;
    let a0 = params.lattice_a;
    let steps = (MAX_LATTICE_SHIFT / SCAN_STEP).round() as usize;
    let grid: Vec<f64> = (0..=2 * steps)
        .map(|i| a0 * (1.0 + (i as f64 - steps as f64) * SCAN_STEP))
        .collect();
    let free_energy = |a: f64| free_energy_per_atom(params, a, n_cells, temperature);
    let values: Vec<f64> = grid.iter().map(|&a| free_energy(a)).collect();
    if values.get(steps).is_none_or(|value| !value.is_finite()) {
        return Err(MdError::General(format!(
            "the Einstein force constant is not positive at lattice_a = {a0} Å, the lattice is \
             dynamically unstable ({})",
            static_pressure_note(params)
        )));
    }
    let Some(index) = nearest_interior_minimum(&values, steps) else {
        return Err(MdError::General(format!(
            "the quasi-harmonic free energy at {temperature} K has no minimum within {:.0} % of \
             lattice_a = {a0} Å ({}); the potential is not in equilibrium there",
            100.0 * MAX_LATTICE_SHIFT,
            static_pressure_note(params)
        )));
    };
    Ok(golden_section_minimum(
        &free_energy,
        grid[index - 1],
        grid[index + 1],
        REFINEMENT_TOLERANCE * a0,
    ))
}

/// Linear thermal expansion coefficient `alpha` (1/K) of a cubic EAM crystal
/// between `t_low` and `t_high` (K) in the quasi-harmonic Einstein model.
///
/// The equilibrium lattice constants `a(T_low)` and `a(T_high)` come from
/// [`quasi_harmonic_lattice_constant`] (minimum of the Helmholtz free energy,
/// including the zero-point energy) and
///
/// ```text
/// alpha = (a(T_high) - a(T_low)) / ( (a(T_high) + a(T_low)) / 2 * (T_high - T_low) )
/// ```
///
/// is their secant slope over the interval.  `n_cells` is the number of
/// conventional cells per supercell edge, at least
/// [`minimum_supercell_cells`].
///
/// Each lattice constant is located to about `1e-8` to `1e-7` Å, so the secant
/// is reliable to better than 1 % for intervals of 1 K and wider; a narrower
/// interval measures the numerical resolution rather than the expansion.
///
/// # Errors
/// [`MdError::General`] unless `0 < t_low < t_high` (both finite), and for every
/// error of [`quasi_harmonic_lattice_constant`]; in particular a potential that
/// has no equilibrium lattice near `lattice_a` yields an error instead of a
/// made-up value.
pub fn thermal_expansion_coefficient(
    params: &EamParams,
    t_low: f64,
    t_high: f64,
    n_cells: usize,
) -> Result<f64, MdError> {
    if !(t_low.is_finite() && t_high.is_finite() && t_low > 0.0 && t_high > t_low) {
        return Err(MdError::General(format!(
            "thermal expansion needs 0 < t_low < t_high, got t_low = {t_low} K, \
             t_high = {t_high} K"
        )));
    }
    let a_low = quasi_harmonic_lattice_constant(params, t_low, n_cells)?;
    let a_high = quasi_harmonic_lattice_constant(params, t_high, n_cells)?;
    Ok((a_high - a_low) / (0.5 * (a_low + a_high) * (t_high - t_low)))
}
