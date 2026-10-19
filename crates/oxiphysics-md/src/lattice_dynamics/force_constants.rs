// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Harmonic force constants and the dynamical matrix `D(q)` of a crystal.
//!
//! In the harmonic approximation the energy of a crystal with `N` atoms per
//! primitive cell is `E = ½ Σ u_α(0 s) Φ^{ss'}_{αβ}(R) u_β(R s')`, a sum over
//! the atoms `s` of the home cell and the atoms `s'` of the cell displaced by
//! the lattice translation `R`.  The real-space *force-constant blocks*
//! `Φ^{ss'}(R)` ([`ForceConstantBlock`]) are the second derivatives of that
//! energy.  Translation invariance of the crystal relates them by the
//! *acoustic sum rule*
//!
//! ```text
//! Σ_{s', R} Φ^{ss'}_{αβ}(R) = 0        for every s, α, β,
//! ```
//!
//! and the symmetry of second derivatives by `Φ^{ss'}(R) = Φ^{s's}(−R)ᵀ`.
//!
//! The Bloch transform with the lattice phase `e^{i q·R}` and the mass
//! weighting gives the Hermitian `3N × 3N` dynamical matrix
//!
//! ```text
//! D^{ss'}_{αβ}(q) = (M_s M_s')^{-1/2} Σ_R Φ^{ss'}_{αβ}(R) e^{i q·R},
//! ```
//!
//! whose eigenvalues are the squared phonon frequencies `ω²(q)`.  Because of
//! the acoustic sum rule three of the `3N` branches vanish at `q = 0`.
//!
//! The phase uses lattice translations only (not the positions of the basis
//! atoms inside the cell), so `D(q)` is periodic in reciprocal space and the
//! frequencies are the same as in the position-phase convention.

use crate::error::MdError;
use oxiphysics_core::numerical_methods::symmetric_eigen_n;

use super::functions::{norm3, scale3};
use super::types::LatticeVectors;

/// Eigenvalues `ω²` of `D(q)` whose magnitude does not exceed this fraction of the
/// size of the contributions that add up to the matrix entries (before they
/// cancel) are indistinguishable from rounding noise and are reported as
/// exactly zero.  The acoustic branches at Γ, where the on-site block cancels
/// the neighbour blocks, are the main case.
const EIGENVALUE_NOISE_FLOOR: f64 = 1.0e-12;

/// Divisions per reciprocal-lattice vector of the default
/// [`DynamicalMatrix::max_frequency`] grid.  The value is divisible by 2, 3 and
/// 4, so the grid contains the usual zone-boundary points.
const MAX_FREQUENCY_GRID_DIVISIONS: usize = 12;

/// Relative tolerance of the consistency checks on force-constant blocks.
const BLOCK_MATCH_TOLERANCE: f64 = 1.0e-8;

/// Acoustic branches count as vanishing at Γ below this fraction of the largest
/// frequency of the crystal.
const ACOUSTIC_FREQUENCY_TOLERANCE: f64 = 1.0e-6;

/// Lattice translations `n₁a₁ + n₂a₂ + n₃a₃` with `|nᵢ|` up to this range are
/// searched for neighbours by [`DynamicalMatrix::new_cubic`].
const NEIGHBOUR_SEARCH_RANGE: i32 = 3;

/// Neighbour distances that agree within this relative tolerance form a shell.
const SHELL_TOLERANCE: f64 = 1.0e-9;

/// One real-space force-constant block `Φ^{ss'}(R)`.
///
/// The block couples atom `s` of the home cell with atom `s'` of the cell
/// displaced by the lattice translation `R`: `Φ_{αβ} = ∂²E / ∂u_α(0 s) ∂u_β(R s')`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ForceConstantBlock {
    /// Lattice translation `R` from the home cell to the cell of the second
    /// atom, in Cartesian components and the length unit of the lattice vectors.
    pub translation: [f64; 3],
    /// Basis index `s` of the first atom, in the home cell.
    pub first: usize,
    /// Basis index `s'` of the second atom, in the cell displaced by `R`.
    pub second: usize,
    /// Second derivatives `Φ_{αβ}` (energy / length²).
    pub phi: [[f64; 3]; 3],
}

impl ForceConstantBlock {
    /// Whether the block is finite in every component.
    fn is_finite(&self) -> bool {
        self.translation.iter().all(|v| v.is_finite())
            && self.phi.iter().flatten().all(|v| v.is_finite())
    }

    /// Whether `other` is the transposed partner `Φ^{s's}(−R)ᵀ` of this block.
    fn is_partner(&self, other: &Self, phi_scale: f64) -> bool {
        let translation_scale = 1.0 + norm3(self.translation);
        let opposite = (0..3).all(|axis| {
            (self.translation[axis] + other.translation[axis]).abs()
                <= BLOCK_MATCH_TOLERANCE * translation_scale
        });
        opposite
            && other.first == self.second
            && other.second == self.first
            && (0..3).all(|a| {
                (0..3).all(|b| {
                    (self.phi[a][b] - other.phi[b][a]).abs() <= BLOCK_MATCH_TOLERANCE * phi_scale
                })
            })
    }
}

/// Mass-weighted dynamical matrix of a crystal from its interatomic force
/// constants.
///
/// The crystal is described by the primitive [`LatticeVectors`], the masses of
/// the `N` basis atoms and the real-space blocks `Φ^{ss'}(R)`.  The Bloch matrix
/// `D(q)` ([`DynamicalMatrix::matrix_at`]) is Hermitian and `3N × 3N`; its
/// eigenvalues give the phonon frequencies ([`DynamicalMatrix::at_q`]).  See
/// the module documentation for the formulas.
#[derive(Debug, Clone)]
pub struct DynamicalMatrix {
    /// Number of atoms per primitive cell.
    pub n_atoms: usize,
    /// Atomic masses, length `n_atoms` (a consistent mass unit with the
    /// force constants).
    pub masses: Vec<f64>,
    /// Real-space force-constant blocks, complete in the sense of
    /// [`DynamicalMatrix::new`]: every block has its transposed partner.
    pub force_constants: Vec<ForceConstantBlock>,
    /// Primitive lattice vectors.
    pub lattice: LatticeVectors,
}

/// Hermitian Bloch matrix `D(q) = Re + i Im` of a crystal at one wave vector.
///
/// Rows and columns are ordered `(atom, α)`, so both parts are `3N × 3N`.
#[derive(Debug, Clone)]
pub struct BlochMatrix {
    /// Real part `Re D`, a symmetric matrix.
    pub real: Vec<Vec<f64>>,
    /// Imaginary part `Im D`, an antisymmetric matrix.
    pub imag: Vec<Vec<f64>>,
    /// Largest sum, over a row, of the absolute values of the terms that add up
    /// to the entries.  It bounds the spectral radius from above and sets the
    /// scale of the rounding errors of the entries and of the eigenvalues.
    magnitude: f64,
}

impl DynamicalMatrix {
    /// Build a dynamical matrix from basis masses and real-space force
    /// constants.
    ///
    /// # Errors
    /// [`MdError::General`] if there are no masses or one is not finite and
    /// positive, a block refers to a basis atom that does not exist or holds a
    /// non-finite number, or a block lacks its transposed partner
    /// `Φ^{s's}(−R)ᵀ` (the symmetry of second derivatives; a block set without
    /// it has no real spectrum).  The acoustic sum rule is *not* required, see
    /// [`DynamicalMatrix::acoustic_sum_rule_violation`] and
    /// [`DynamicalMatrix::enforce_acoustic_sum_rule`].
    pub fn new(
        masses: Vec<f64>,
        force_constants: Vec<ForceConstantBlock>,
        lattice: LatticeVectors,
    ) -> Result<Self, MdError> {
        if masses.is_empty() {
            return Err(MdError::General(
                "DynamicalMatrix needs at least one basis atom".into(),
            ));
        }
        if let Some(atom) = masses.iter().position(|&m| !(m.is_finite() && m > 0.0)) {
            return Err(MdError::General(format!(
                "DynamicalMatrix: mass of basis atom {atom} must be finite and positive, got {}",
                masses[atom]
            )));
        }
        let n_atoms = masses.len();
        for (index, block) in force_constants.iter().enumerate() {
            if block.first >= n_atoms || block.second >= n_atoms {
                return Err(MdError::General(format!(
                    "DynamicalMatrix: block {index} couples basis atoms {} and {} but there are only {n_atoms}",
                    block.first, block.second
                )));
            }
            if !block.is_finite() {
                return Err(MdError::General(format!(
                    "DynamicalMatrix: block {index} contains a non-finite number"
                )));
            }
        }
        let phi_scale = force_constants
            .iter()
            .flat_map(|block| block.phi.iter().flatten())
            .fold(0.0_f64, |scale, value| scale.max(value.abs()));
        for (index, block) in force_constants.iter().enumerate() {
            if !force_constants
                .iter()
                .any(|other| block.is_partner(other, phi_scale))
            {
                return Err(MdError::General(format!(
                    "DynamicalMatrix: block {index} (basis atoms {} -> {}, translation {:?}) has no transposed partner at the opposite translation",
                    block.first, block.second, block.translation
                )));
            }
        }
        Ok(Self {
            n_atoms,
            masses,
            force_constants,
            lattice,
        })
    }

    /// Construct the dynamical matrix of a mono-atomic crystal whose atoms are
    /// joined to their nearest neighbours by central springs.
    ///
    /// Every neighbour at the lattice translation `R` is coupled by a spring of
    /// stiffness `k_spring` along the bond direction `r̂ = R / |R|`, i.e. by
    /// the block `Φ(R) = −k r̂ r̂ᵀ`; the on-site block `Φ(0) = k Σ_R r̂ r̂ᵀ`
    /// satisfies the acoustic sum rule.  For a simple cubic lattice this gives
    /// `ω²(q) = (2k/m) Σ_α (1 − cos q_α a)` per polarisation `α`.
    ///
    /// Neighbours are taken shell by shell, nearest first, until at least
    /// `coordination` of them are connected; a shell is never split, which
    /// keeps the point-group symmetry of the lattice.  The standard
    /// coordinations select exactly the nearest shell: 6 (simple cubic),
    /// 8 (body-centred cubic) and 12 (face-centred cubic).  `coordination = 0`
    /// gives uncoupled atoms with vanishing frequencies.  Only translations
    /// `n₁a₁ + n₂a₂ + n₃a₃` with `|nᵢ| ≤ 3` are searched, so the primitive
    /// vectors should be reasonably reduced.
    ///
    /// # Arguments
    /// * `k_spring`    – spring constant (energy / length²); a negative value
    ///   describes an unstable lattice and yields negative frequencies
    /// * `mass`        – atomic mass
    /// * `coordination`– number of nearest neighbours (e.g. 6 for simple cubic)
    /// * `lattice`     – primitive lattice vectors
    pub fn new_cubic(
        k_spring: f64,
        mass: f64,
        coordination: usize,
        lattice: LatticeVectors,
    ) -> Self {
        let mut blocks = Vec::new();
        let mut on_site = [[0.0_f64; 3]; 3];
        for translation in neighbour_translations(&lattice, coordination) {
            let direction = scale3(translation, 1.0 / norm3(translation));
            let mut phi = [[0.0_f64; 3]; 3];
            for (a, row) in phi.iter_mut().enumerate() {
                for (b, entry) in row.iter_mut().enumerate() {
                    let projector = direction[a] * direction[b];
                    *entry = -k_spring * projector;
                    on_site[a][b] += k_spring * projector;
                }
            }
            blocks.push(ForceConstantBlock {
                translation,
                first: 0,
                second: 0,
                phi,
            });
        }
        blocks.push(ForceConstantBlock {
            translation: [0.0; 3],
            first: 0,
            second: 0,
            phi: on_site,
        });
        Self {
            n_atoms: 1,
            masses: vec![mass],
            force_constants: blocks,
            lattice,
        }
    }

    /// Length scale of the lattice, the longest primitive vector.
    fn length_scale(&self) -> f64 {
        [self.lattice.a1, self.lattice.a2, self.lattice.a3]
            .iter()
            .map(|&vector| norm3(vector))
            .fold(0.0_f64, f64::max)
    }

    /// Whether `block` couples an atom with itself in the home cell.
    fn is_on_site(&self, block: &ForceConstantBlock) -> bool {
        block.first == block.second
            && norm3(block.translation) <= BLOCK_MATCH_TOLERANCE * self.length_scale()
    }

    /// Largest absolute entry of the acoustic sum `Σ_{s', R} Φ^{ss'}(R)` over
    /// all basis atoms `s`.
    ///
    /// It vanishes for a translation-invariant crystal; otherwise the acoustic
    /// branches do not reach zero at Γ.
    pub fn acoustic_sum_rule_violation(&self) -> f64 {
        let mut sums = vec![[[0.0_f64; 3]; 3]; self.masses.len()];
        for block in &self.force_constants {
            if let Some(sum) = sums.get_mut(block.first) {
                for (a, row) in sum.iter_mut().enumerate() {
                    for (b, entry) in row.iter_mut().enumerate() {
                        *entry += block.phi[a][b];
                    }
                }
            }
        }
        sums.iter()
            .flatten()
            .flatten()
            .fold(0.0_f64, |worst, value| worst.max(value.abs()))
    }

    /// Replace the on-site blocks `Φ^{ss}(0)` by `−Σ_{(s', R) ≠ (s, 0)} Φ^{ss'}(R)`
    /// so that the acoustic sum rule holds.
    ///
    /// The on-site block is kept symmetric (and the dynamical matrix
    /// Hermitian) by using the symmetric part of the sum, which equals the sum
    /// itself for any translation-invariant set of force constants.
    pub fn enforce_acoustic_sum_rule(&mut self) {
        let n_atoms = self.masses.len();
        let mut sums = vec![[[0.0_f64; 3]; 3]; n_atoms];
        for block in &self.force_constants {
            if self.is_on_site(block) {
                continue;
            }
            if let Some(sum) = sums.get_mut(block.first) {
                for (a, row) in sum.iter_mut().enumerate() {
                    for (b, entry) in row.iter_mut().enumerate() {
                        *entry += block.phi[a][b];
                    }
                }
            }
        }
        let kept: Vec<ForceConstantBlock> = self
            .force_constants
            .iter()
            .filter(|block| !self.is_on_site(block))
            .copied()
            .collect();
        self.force_constants = kept;
        for (atom, sum) in sums.iter().enumerate() {
            let mut phi = [[0.0_f64; 3]; 3];
            for (a, row) in phi.iter_mut().enumerate() {
                for (b, entry) in row.iter_mut().enumerate() {
                    *entry = -0.5 * (sum[a][b] + sum[b][a]);
                }
            }
            self.force_constants.push(ForceConstantBlock {
                translation: [0.0; 3],
                first: atom,
                second: atom,
                phi,
            });
        }
    }

    /// The Hermitian dynamical matrix `D(q)` at the Cartesian wave vector `q`.
    ///
    /// Only the Hermitian part `(D + D†)/2` is returned, so a block list that
    /// was assembled by hand without [`DynamicalMatrix::new`]'s checks still
    /// has a real spectrum.  Returns `None` if a mass is not finite and
    /// positive.  Blocks that refer to a missing basis atom are ignored.
    pub fn matrix_at(&self, q: [f64; 3]) -> Option<BlochMatrix> {
        let n_atoms = self.masses.len();
        if self.masses.iter().any(|&m| !(m.is_finite() && m > 0.0)) {
            return None;
        }
        let dim = 3 * n_atoms;
        let mut real = vec![vec![0.0_f64; dim]; dim];
        let mut imag = vec![vec![0.0_f64; dim]; dim];
        let mut row_magnitude = vec![0.0_f64; dim];
        for block in &self.force_constants {
            if block.first >= n_atoms || block.second >= n_atoms {
                continue;
            }
            let weight = 1.0 / (self.masses[block.first] * self.masses[block.second]).sqrt();
            let phase = q[0] * block.translation[0]
                + q[1] * block.translation[1]
                + q[2] * block.translation[2];
            let (sin_phase, cos_phase) = phase.sin_cos();
            for a in 0..3 {
                for b in 0..3 {
                    let value = weight * block.phi[a][b];
                    real[3 * block.first + a][3 * block.second + b] += value * cos_phase;
                    imag[3 * block.first + a][3 * block.second + b] += value * sin_phase;
                    row_magnitude[3 * block.first + a] += value.abs();
                }
            }
        }
        for row in 0..dim {
            for col in (row + 1)..dim {
                let symmetric = 0.5 * (real[row][col] + real[col][row]);
                real[row][col] = symmetric;
                real[col][row] = symmetric;
                let antisymmetric = 0.5 * (imag[row][col] - imag[col][row]);
                imag[row][col] = antisymmetric;
                imag[col][row] = -antisymmetric;
            }
            imag[row][row] = 0.0;
        }
        let magnitude = row_magnitude.into_iter().fold(0.0_f64, f64::max);
        Some(BlochMatrix {
            real,
            imag,
            magnitude,
        })
    }

    /// Squared phonon frequencies `ω²`, the `3N` eigenvalues of `D(q)` in
    /// ascending order.
    ///
    /// The Hermitian matrix `A + iB` is diagonalised through the real symmetric
    /// `6N × 6N` matrix `[[A, −B], [B, A]]`, whose spectrum holds every
    /// eigenvalue twice.  Eigenvalues below `10⁻¹²` of the size of the terms
    /// that add up to the matrix entries cannot be told from rounding noise and
    /// are returned as exactly zero, which makes the acoustic branches at Γ
    /// exactly zero.  A negative eigenvalue is a genuine instability of the
    /// lattice.  All entries are NaN if a mass is not finite and positive.
    pub fn squared_frequencies_at(&self, q: [f64; 3]) -> Vec<f64> {
        let n_modes = 3 * self.masses.len();
        let Some(matrix) = self.matrix_at(q) else {
            return vec![f64::NAN; n_modes];
        };
        let mut embedded = vec![vec![0.0_f64; 2 * n_modes]; 2 * n_modes];
        for row in 0..n_modes {
            for col in 0..n_modes {
                embedded[row][col] = matrix.real[row][col];
                embedded[row][col + n_modes] = -matrix.imag[row][col];
                embedded[row + n_modes][col] = matrix.imag[row][col];
                embedded[row + n_modes][col + n_modes] = matrix.real[row][col];
            }
        }
        let (doubled, _) = symmetric_eigen_n(&embedded);
        let floor = EIGENVALUE_NOISE_FLOOR * matrix.magnitude;
        (0..n_modes)
            .map(|mode| {
                let eigenvalue = 0.5 * (doubled[2 * mode] + doubled[2 * mode + 1]);
                if eigenvalue.abs() <= floor {
                    0.0
                } else {
                    eigenvalue
                }
            })
            .collect()
    }

    /// Phonon frequencies `ω = sign(ω²) √|ω²|` at the Cartesian wave vector
    /// `q`, `3N` values in ascending order.
    ///
    /// Unstable modes (`ω² < 0`) are reported as negative frequencies.  The
    /// wave vector is in radians per length unit of the lattice vectors: the
    /// zone boundary of a simple cubic lattice with lattice constant `a` lies
    /// at `q = π/a`.  See [`DynamicalMatrix::squared_frequencies_at`] for the
    /// numerical details.
    ///
    /// # Arguments
    /// * `q` – wave vector (rad / length)
    pub fn at_q(&self, q: [f64; 3]) -> Vec<f64> {
        self.squared_frequencies_at(q)
            .into_iter()
            .map(|eigenvalue| {
                if eigenvalue == 0.0 {
                    0.0
                } else {
                    eigenvalue.signum() * eigenvalue.abs().sqrt()
                }
            })
            .collect()
    }

    /// Phonon frequencies at the zone centre Γ (`q = 0`).
    ///
    /// Three of them, the acoustic branches, vanish when the force constants
    /// obey the acoustic sum rule; the remaining `3N − 3` are the optical
    /// frequencies of the crystal.
    pub fn phonon_frequencies_gamma(&self) -> Vec<f64> {
        self.at_q([0.0; 3])
    }

    /// Whether the three lowest branches vanish at Γ, i.e. stay below
    /// `10⁻⁶` of the largest frequency of the crystal.
    pub fn acoustic_at_gamma(&self) -> bool {
        let tolerance = ACOUSTIC_FREQUENCY_TOLERANCE * self.max_frequency();
        self.phonon_frequencies_gamma()
            .iter()
            .take(3)
            .all(|&omega| omega.abs() <= tolerance)
    }

    /// Largest phonon frequency found on a uniform `n × n × n` grid of wave
    /// vectors `q = (i b₁ + j b₂ + k b₃) / n` (`0 ≤ i, j, k < n`) covering one
    /// reciprocal cell, or `0` if no frequency is positive.
    ///
    /// The grid contains the zone-boundary points whose reduced coordinates
    /// are multiples of `1/n`, so the maximum is exact for models whose
    /// extremum lies there (nearest-neighbour central forces on the cubic
    /// lattices) and a lower bound otherwise.  A degenerate lattice is probed
    /// at Γ only.
    pub fn max_frequency_on_grid(&self, divisions: usize) -> f64 {
        let n = divisions.max(1);
        let volume = self.lattice.volume();
        if !(volume.is_finite() && volume.abs() > f64::MIN_POSITIVE) {
            return top_frequency(&self.phonon_frequencies_gamma());
        }
        let (b1, b2, b3) = self.lattice.reciprocal_vectors();
        let mut largest = 0.0_f64;
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let (fi, fj, fk) = (
                        i as f64 / n as f64,
                        j as f64 / n as f64,
                        k as f64 / n as f64,
                    );
                    let q = [
                        fi * b1[0] + fj * b2[0] + fk * b3[0],
                        fi * b1[1] + fj * b2[1] + fk * b3[1],
                        fi * b1[2] + fj * b2[2] + fk * b3[2],
                    ];
                    largest = largest.max(top_frequency(&self.at_q(q)));
                }
            }
        }
        largest
    }

    /// Largest phonon frequency of the crystal, sampled on a 12 × 12 × 12 grid
    /// of the reciprocal cell (see [`DynamicalMatrix::max_frequency_on_grid`]).
    pub fn max_frequency(&self) -> f64 {
        self.max_frequency_on_grid(MAX_FREQUENCY_GRID_DIVISIONS)
    }
}

/// Highest entry of an ascending frequency list (`0` for an empty list).
fn top_frequency(frequencies: &[f64]) -> f64 {
    frequencies.last().copied().unwrap_or(0.0)
}

/// Lattice translations of the nearest-neighbour shells that together hold at
/// least `coordination` sites.
fn neighbour_translations(lattice: &LatticeVectors, coordination: usize) -> Vec<[f64; 3]> {
    if coordination == 0 {
        return Vec::new();
    }
    let length_scale = [lattice.a1, lattice.a2, lattice.a3]
        .iter()
        .map(|&vector| norm3(vector))
        .fold(0.0_f64, f64::max);
    let mut candidates: Vec<[f64; 3]> = Vec::new();
    for n1 in -NEIGHBOUR_SEARCH_RANGE..=NEIGHBOUR_SEARCH_RANGE {
        for n2 in -NEIGHBOUR_SEARCH_RANGE..=NEIGHBOUR_SEARCH_RANGE {
            for n3 in -NEIGHBOUR_SEARCH_RANGE..=NEIGHBOUR_SEARCH_RANGE {
                let translation = [0, 1, 2].map(|axis| {
                    f64::from(n1) * lattice.a1[axis]
                        + f64::from(n2) * lattice.a2[axis]
                        + f64::from(n3) * lattice.a3[axis]
                });
                if norm3(translation) > BLOCK_MATCH_TOLERANCE * length_scale {
                    candidates.push(translation);
                }
            }
        }
    }
    // Nearest first.  Translations of one shell differ in length only by
    // rounding noise, so their relative order is irrelevant: whole shells are
    // taken below.
    candidates.sort_by(|a, b| norm3(*a).total_cmp(&norm3(*b)));
    let mut selected: Vec<[f64; 3]> = Vec::new();
    let mut start = 0;
    while start < candidates.len() && selected.len() < coordination {
        let shell_length = norm3(candidates[start]);
        let mut end = start;
        while end < candidates.len()
            && (norm3(candidates[end]) - shell_length).abs() <= SHELL_TOLERANCE * shell_length
        {
            end += 1;
        }
        selected.extend_from_slice(&candidates[start..end]);
        start = end;
    }
    selected
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice_dynamics::{DiatomicChain, PhononDispersion};
    use std::f64::consts::PI;

    const K: f64 = 1.7;
    const M: f64 = 2.3;

    fn assert_close(actual: &[f64], expected: &[f64], tol: f64, context: &str) {
        assert_eq!(actual.len(), expected.len(), "{context}: mode count");
        for (mode, (a, e)) in actual.iter().zip(expected).enumerate() {
            assert!(
                (a - e).abs() <= tol * (1.0 + e.abs()),
                "{context}: mode {mode}: {actual:?} vs {expected:?}"
            );
        }
    }

    fn sorted(mut values: Vec<f64>) -> Vec<f64> {
        values.sort_by(f64::total_cmp);
        values
    }

    fn dot(q: [f64; 3], r: [f64; 3]) -> f64 {
        q[0] * r[0] + q[1] * r[1] + q[2] * r[2]
    }

    /// Nearest-neighbour translations of the simple cubic lattice.
    fn sc_neighbours(a: f64) -> Vec<[f64; 3]> {
        let mut out = Vec::new();
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let mut r = [0.0; 3];
                r[axis] = sign * a;
                out.push(r);
            }
        }
        out
    }

    /// Nearest-neighbour translations of the face-centred cubic lattice: the 12
    /// vectors `(a/2)(±1, ±1, 0)` and their permutations.
    fn fcc_neighbours(a: f64) -> Vec<[f64; 3]> {
        let mut out = Vec::new();
        for zero_axis in 0..3 {
            let (u, v) = ((zero_axis + 1) % 3, (zero_axis + 2) % 3);
            for s1 in [-0.5, 0.5] {
                for s2 in [-0.5, 0.5] {
                    let mut r = [0.0; 3];
                    r[u] = s1 * a;
                    r[v] = s2 * a;
                    out.push(r);
                }
            }
        }
        out
    }

    /// Nearest-neighbour translations of the body-centred cubic lattice.
    fn bcc_neighbours(a: f64) -> Vec<[f64; 3]> {
        let mut out = Vec::new();
        for s1 in [-0.5, 0.5] {
            for s2 in [-0.5, 0.5] {
                for s3 in [-0.5, 0.5] {
                    out.push([s1 * a, s2 * a, s3 * a]);
                }
            }
        }
        out
    }

    #[test]
    fn simple_cubic_acoustic_branches_vanish_at_gamma() {
        let dm = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(1.3));
        assert_eq!(dm.phonon_frequencies_gamma(), vec![0.0; 3]);
        assert!(dm.acoustic_at_gamma());
        assert!(dm.acoustic_sum_rule_violation() < 1.0e-12);
        // Six neighbour blocks and the on-site block 2k * identity.
        assert_eq!(dm.force_constants.len(), 7);
        let on_site = dm
            .force_constants
            .iter()
            .find(|block| dm.is_on_site(block))
            .expect("on-site block");
        for a in 0..3 {
            for b in 0..3 {
                let expected = if a == b { 2.0 * K } else { 0.0 };
                assert!(
                    (on_site.phi[a][b] - expected).abs() < 1.0e-12,
                    "{on_site:?}"
                );
            }
        }
    }

    #[test]
    fn simple_cubic_dispersion_is_the_analytic_nearest_neighbour_formula() {
        let a = 1.3;
        let dm = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(a));
        let points = [
            [0.7 / a, -1.1 / a, 2.3 / a],
            [PI / a, 0.0, 0.0],
            [PI / a, PI / a, 0.0],
            [PI / a, PI / a, PI / a],
            [0.3 / a, 0.0, 0.0],
        ];
        for q in points {
            // The matrix is diagonal: omega_alpha^2 = (2k/m)(1 - cos q_alpha a).
            let expected = sorted(
                (0..3)
                    .map(|i| (2.0 * K / M * (1.0 - (q[i] * a).cos())).sqrt())
                    .collect(),
            );
            assert_close(&dm.at_q(q), &expected, 1.0e-9, &format!("q = {q:?}"));
        }
        // X point: two transverse modes at rest, the longitudinal one at 2 sqrt(k/m).
        let x = dm.at_q([PI / a, 0.0, 0.0]);
        assert_close(&x, &[0.0, 0.0, 2.0 * (K / M).sqrt()], 1.0e-9, "X");
    }

    #[test]
    fn long_wavelength_limit_is_the_sound_velocity() {
        let a = 1.3;
        let dm = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(a));
        let sound_velocity = a * (K / M).sqrt();
        let q = 1.0e-3 / a;
        let omega = dm.at_q([q, 0.0, 0.0]);
        assert_eq!(omega[0], 0.0, "transverse branch");
        assert!(
            (omega[2] / q - sound_velocity).abs() < 1.0e-5 * sound_velocity,
            "omega / q = {} vs v = {sound_velocity}",
            omega[2] / q
        );
        // Linear dispersion: doubling q doubles omega.
        let doubled = dm.at_q([2.0 * q, 0.0, 0.0]);
        assert!((doubled[2] / omega[2] - 2.0).abs() < 1.0e-6);
    }

    #[test]
    fn fcc_zone_boundary_frequencies_match_the_analytic_values() {
        let a = 1.9;
        let dm = DynamicalMatrix::new_cubic(K, M, 12, LatticeVectors::fcc(a));
        assert_eq!(dm.force_constants.len(), 13);
        assert_eq!(dm.phonon_frequencies_gamma(), vec![0.0; 3]);
        assert!(dm.acoustic_sum_rule_violation() < 1.0e-12);
        // X = (2 pi / a)(1, 0, 0): omega_T^2 = 4k/m (twice), omega_L^2 = 8k/m.
        let x = dm.at_q([2.0 * PI / a, 0.0, 0.0]);
        let (t, l) = ((4.0 * K / M).sqrt(), (8.0 * K / M).sqrt());
        assert_close(&x, &[t, t, l], 1.0e-9, "X");
        // L = (pi / a)(1, 1, 1): omega_T^2 = 2k/m (twice), omega_L^2 = 8k/m.
        let point_l = dm.at_q([PI / a, PI / a, PI / a]);
        let t = (2.0 * K / M).sqrt();
        assert_close(&point_l, &[t, t, l], 1.0e-9, "L");
    }

    #[test]
    fn bcc_models_the_known_soft_branch_and_h_point() {
        // Nearest-neighbour central forces alone leave a flat zero branch along
        // Gamma-N for the bcc lattice; at H all three modes are degenerate.
        let a = 1.4;
        let dm = DynamicalMatrix::new_cubic(K, M, 8, LatticeVectors::bcc(a));
        assert_eq!(dm.force_constants.len(), 9);
        let h = dm.at_q([2.0 * PI / a, 0.0, 0.0]);
        let expected = (16.0 * K / (3.0 * M)).sqrt();
        assert_close(&h, &[expected; 3], 1.0e-9, "H");
        for qx in [0.3, 0.9, 1.7] {
            let along_n = dm.at_q([qx / a, qx / a, 0.0]);
            assert_eq!(along_n[0], 0.0, "qx = {qx}: {along_n:?}");
            assert!(along_n[1] > 0.0 && along_n[2] > 0.0);
        }
    }

    #[test]
    fn trace_of_the_dynamical_matrix_is_the_neighbour_sum_on_every_lattice() {
        // tr D(q) = sum_alpha omega_alpha^2 = (k/m) sum_R (1 - cos q.R), an
        // identity that involves no eigenvalue solver.
        let a = 1.1;
        let cases = [
            (LatticeVectors::cubic(a), 6, sc_neighbours(a)),
            (LatticeVectors::fcc(a), 12, fcc_neighbours(a)),
            (LatticeVectors::bcc(a), 8, bcc_neighbours(a)),
        ];
        for (lattice, coordination, neighbours) in cases {
            let dm = DynamicalMatrix::new_cubic(K, M, coordination, lattice);
            assert_eq!(dm.force_constants.len(), coordination + 1);
            for q in [[0.4, -0.9, 1.3], [2.0, 0.5, -0.7], [-1.6, 2.2, 0.1]] {
                let squared = dm.squared_frequencies_at(q);
                let trace: f64 = squared.iter().sum();
                let expected: f64 = neighbours
                    .iter()
                    .map(|&r| K / M * (1.0 - dot(q, r).cos()))
                    .sum();
                assert!(
                    (trace - expected).abs() < 1.0e-9 * (1.0 + expected),
                    "z = {coordination}, q = {q:?}: trace {trace} vs {expected}"
                );
                assert!(squared.iter().all(|&w| w >= 0.0), "{squared:?}");
            }
        }
    }

    #[test]
    fn frequencies_are_even_in_q_and_periodic_in_the_reciprocal_lattice() {
        let dm = DynamicalMatrix::new_cubic(K, M, 12, LatticeVectors::fcc(1.9));
        let (b1, b2, b3) = dm.lattice.reciprocal_vectors();
        let q = [0.37, -1.21, 0.88];
        let reference = dm.at_q(q);
        assert_close(
            &dm.at_q([-q[0], -q[1], -q[2]]),
            &reference,
            1.0e-9,
            "inversion",
        );
        let shifts = [
            b1,
            b2,
            b3,
            [b1[0] + b2[0], b1[1] + b2[1], b1[2] + b2[2]],
            [b3[0] - b1[0], b3[1] - b1[1], b3[2] - b1[2]],
        ];
        for g in shifts {
            let shifted = [q[0] + g[0], q[1] + g[1], q[2] + g[2]];
            assert_close(&dm.at_q(shifted), &reference, 1.0e-9, &format!("G = {g:?}"));
        }
    }

    #[test]
    fn frequencies_scale_with_the_square_root_of_stiffness_over_mass() {
        let lattice = LatticeVectors::cubic(1.0);
        let q = [0.9, 1.4, -0.6];
        let base = DynamicalMatrix::new_cubic(K, M, 6, lattice.clone()).at_q(q);
        let stiffer = DynamicalMatrix::new_cubic(4.0 * K, M, 6, lattice.clone()).at_q(q);
        let heavier = DynamicalMatrix::new_cubic(K, 4.0 * M, 6, lattice).at_q(q);
        let doubled: Vec<f64> = base.iter().map(|w| 2.0 * w).collect();
        let halved: Vec<f64> = base.iter().map(|w| 0.5 * w).collect();
        assert_close(&stiffer, &doubled, 1.0e-9, "4k");
        assert_close(&heavier, &halved, 1.0e-9, "4m");
    }

    #[test]
    fn coordination_selects_whole_neighbour_shells() {
        let neighbours = |lattice: LatticeVectors, coordination: usize| {
            DynamicalMatrix::new_cubic(K, M, coordination, lattice)
                .force_constants
                .len()
                - 1
        };
        // Simple cubic shells hold 6, 12, 8 and 6 sites.
        let cubic = LatticeVectors::cubic(1.0);
        assert_eq!(neighbours(cubic.clone(), 1), 6);
        assert_eq!(neighbours(cubic.clone(), 6), 6);
        assert_eq!(neighbours(cubic.clone(), 7), 18);
        assert_eq!(neighbours(cubic.clone(), 18), 18);
        assert_eq!(neighbours(cubic.clone(), 19), 26);
        // fcc: 12, 6, 24; bcc: 8, 6, 12.
        assert_eq!(neighbours(LatticeVectors::fcc(1.0), 12), 12);
        assert_eq!(neighbours(LatticeVectors::fcc(1.0), 13), 18);
        assert_eq!(neighbours(LatticeVectors::bcc(1.0), 8), 8);
        assert_eq!(neighbours(LatticeVectors::bcc(1.0), 9), 14);
        // No neighbours: free atoms.
        let free = DynamicalMatrix::new_cubic(K, M, 0, cubic);
        assert_eq!(free.force_constants.len(), 1);
        assert_eq!(free.at_q([0.7, 0.2, 0.1]), vec![0.0; 3]);
        assert_eq!(free.max_frequency(), 0.0);
    }

    #[test]
    fn a_second_neighbour_shell_stiffens_the_zone_boundary() {
        // Coordination 18 adds the 12 next-nearest sites at a sqrt(2) to the 6
        // nearest ones.  Check the trace against the explicit lattice sum and
        // the longitudinal X-point frequency against the nearest-neighbour model.
        let a = 1.0;
        let dm = DynamicalMatrix::new_cubic(K, M, 18, LatticeVectors::cubic(a));
        assert_eq!(dm.force_constants.len(), 19);
        assert!(dm.acoustic_sum_rule_violation() < 1.0e-12);
        let q = [PI / a, 0.0, 0.0];
        let first_shell: Vec<[f64; 3]> = sc_neighbours(a);
        let mut second_shell = Vec::new();
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            for s1 in [-1.0, 1.0] {
                for s2 in [-1.0, 1.0] {
                    let mut r = [0.0; 3];
                    r[i] = s1 * a;
                    r[j] = s2 * a;
                    second_shell.push(r);
                }
            }
        }
        let trace: f64 = dm.squared_frequencies_at(q).iter().sum();
        let expected: f64 = first_shell
            .iter()
            .chain(second_shell.iter())
            .map(|&r| K / M * (1.0 - dot(q, r).cos()))
            .sum();
        assert!(
            (trace - expected).abs() < 1.0e-9 * (1.0 + expected),
            "{trace} vs {expected}"
        );
        let nearest_only = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(a));
        assert!(dm.at_q(q)[2] > nearest_only.at_q(q)[2]);
    }

    #[test]
    fn negative_spring_constant_reports_unstable_modes_as_negative_frequencies() {
        let a = 1.0;
        let dm = DynamicalMatrix::new_cubic(-K, M, 6, LatticeVectors::cubic(a));
        let at_x = dm.at_q([PI / a, 0.0, 0.0]);
        assert_close(
            &at_x,
            &[-2.0 * (K / M).sqrt(), 0.0, 0.0],
            1.0e-9,
            "unstable X",
        );
        assert_eq!(dm.max_frequency(), 0.0);
    }

    #[test]
    fn max_frequency_samples_the_brillouin_zone_not_just_gamma() {
        let a = 1.3;
        let simple_cubic = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(a));
        // All three branches meet at 2 sqrt(k/m) at the corner (pi/a)(1, 1, 1).
        let expected = 2.0 * (K / M).sqrt();
        assert!((simple_cubic.max_frequency() - expected).abs() < 1.0e-9);
        assert_eq!(
            simple_cubic.max_frequency_on_grid(1),
            0.0,
            "grid of Gamma only"
        );
        // fcc: the grid includes X, where omega_L = sqrt(8 k / m).
        let fcc = DynamicalMatrix::new_cubic(K, M, 12, LatticeVectors::fcc(a));
        assert!(fcc.max_frequency() >= (8.0 * K / M).sqrt() - 1.0e-9);
        // The old Gamma-only estimate was zero for an acoustic crystal.
        assert_eq!(fcc.phonon_frequencies_gamma().last().copied(), Some(0.0));
    }

    /// Diatomic chain along x: masses `m1`, `m2` alternate, neighbours are joined
    /// by springs `k`, and the transverse cell edges are large and uncoupled.
    fn diatomic_chain_dynamical_matrix(m1: f64, m2: f64, k: f64, a: f64) -> DynamicalMatrix {
        let spring = [[-k, 0.0, 0.0], [0.0; 3], [0.0; 3]];
        let block = |first: usize, second: usize, x: f64| ForceConstantBlock {
            translation: [x, 0.0, 0.0],
            first,
            second,
            phi: spring,
        };
        // Bond inside the cell (0 - 1) and across the cell boundary (1 - 0'),
        // each listed with its transposed partner.
        let blocks = vec![
            block(0, 1, 0.0),
            block(1, 0, 0.0),
            block(1, 0, a),
            block(0, 1, -a),
        ];
        let lattice = LatticeVectors::new([a, 0.0, 0.0], [0.0, 3.0 * a, 0.0], [0.0, 0.0, 3.0 * a]);
        let mut dm = DynamicalMatrix::new(vec![m1, m2], blocks, lattice).expect("valid chain");
        dm.enforce_acoustic_sum_rule();
        dm
    }

    #[test]
    fn diatomic_chain_matches_the_analytic_acoustic_and_optical_branches() {
        let (m1, m2, k, a) = (1.0, 3.0, 2.5, 1.4);
        let dm = diatomic_chain_dynamical_matrix(m1, m2, k, a);
        assert_eq!(dm.n_atoms, 2);
        assert!(dm.acoustic_sum_rule_violation() < 1.0e-12);
        let mut chain = DiatomicChain::new(m1, m2, k);
        chain.a = a;
        for fraction in [0.0, 0.1, 0.25, 0.5, 0.8, 1.0] {
            let q = fraction * PI / a;
            let omega = dm.at_q([q, 0.0, 0.0]);
            assert_eq!(omega.len(), 6);
            // y and z motion is not coupled: four zero modes.
            assert_close(
                &omega[..4],
                &[0.0; 4],
                1.0e-9,
                &format!("transverse, q = {q}"),
            );
            assert_close(
                &omega[4..],
                &[chain.acoustic_frequency(q), chain.optical_frequency(q)],
                1.0e-9,
                &format!("longitudinal, q = {q}"),
            );
        }
        // Optical gap at Gamma: omega^2 = k (1/m1 + 1/m2) * 2.
        let gamma = dm.phonon_frequencies_gamma();
        assert_close(
            &gamma[4..],
            &[0.0, (2.0 * k * (1.0 / m1 + 1.0 / m2)).sqrt()],
            1.0e-9,
            "Gamma",
        );
        assert!(dm.acoustic_at_gamma());
    }

    #[test]
    fn dynamical_matrix_of_the_chain_is_hermitian_with_the_expected_phases() {
        let (m1, m2, k, a) = (1.0, 3.0, 2.5, 1.4);
        let dm = diatomic_chain_dynamical_matrix(m1, m2, k, a);
        let q = 0.6 / a;
        let BlochMatrix { real, imag, .. } = dm.matrix_at([q, 0.0, 0.0]).expect("positive masses");
        let weight = 1.0 / (m1 * m2).sqrt();
        // D^{01}_{xx} = -k (1 + exp(-i q a)) / sqrt(m1 m2)
        assert!((real[0][3] + k * weight * (1.0 + (q * a).cos())).abs() < 1.0e-12);
        assert!((imag[0][3] - k * weight * (q * a).sin()).abs() < 1.0e-12);
        for row in 0..6 {
            for col in 0..6 {
                assert!(
                    (real[row][col] - real[col][row]).abs() < 1.0e-14,
                    "Re symmetric"
                );
                assert!(
                    (imag[row][col] + imag[col][row]).abs() < 1.0e-14,
                    "Im antisymmetric"
                );
            }
        }
        // On-site terms 2k / m.
        assert!((real[0][0] - 2.0 * k / m1).abs() < 1.0e-12);
        assert!((real[3][3] - 2.0 * k / m2).abs() < 1.0e-12);
    }

    #[test]
    fn dispersion_path_of_a_simple_cubic_crystal_reaches_the_analytic_points() {
        let dm = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(1.0));
        let dispersion = PhononDispersion::compute_sc(dm, 10);
        assert!(dispersion.acoustic_modes_at_gamma(1.0e-9));
        let top = 2.0 * (K / M).sqrt();
        // Gamma -> X -> M -> Gamma with 10 points per segment.
        assert_close(&dispersion.frequencies[10], &[0.0, 0.0, top], 1.0e-9, "X");
        assert_close(&dispersion.frequencies[20], &[0.0, top, top], 1.0e-9, "M");
        assert_close(
            &dispersion.frequencies[30],
            &[0.0; 3],
            1.0e-9,
            "closing Gamma",
        );
        assert!((dispersion.max_frequency() - top).abs() < 1.0e-9);
    }

    #[test]
    fn acoustic_check_ignores_optical_branches_of_a_diatomic_crystal() {
        let dm = diatomic_chain_dynamical_matrix(1.0, 3.0, 2.5, 1.4);
        let dispersion = PhononDispersion::compute_sc(dm, 4);
        let at_gamma = &dispersion.frequencies[0];
        assert!(at_gamma[5] > 1.0, "optical branch at Gamma: {at_gamma:?}");
        assert!(dispersion.acoustic_modes_at_gamma(1.0e-9));
    }

    #[test]
    fn violated_acoustic_sum_rule_is_detected_and_repaired() {
        let mut dm = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(1.0));
        let on_site = (0..dm.force_constants.len())
            .find(|&i| dm.is_on_site(&dm.force_constants[i]))
            .expect("on-site block");
        dm.force_constants[on_site].phi[0][0] *= 1.1;
        assert!(dm.acoustic_sum_rule_violation() > 0.1);
        assert!(!dm.acoustic_at_gamma());
        let gamma = dm.phonon_frequencies_gamma();
        assert!(
            (gamma[2] - (0.2 * K / M).sqrt()).abs() < 1.0e-9,
            "{gamma:?}"
        );

        dm.enforce_acoustic_sum_rule();
        assert!(dm.acoustic_sum_rule_violation() < 1.0e-12);
        assert!(dm.acoustic_at_gamma());
        assert_eq!(dm.force_constants.len(), 7, "one on-site block per atom");
    }

    fn block(first: usize, second: usize, x: f64, value: f64) -> ForceConstantBlock {
        ForceConstantBlock {
            translation: [x, 0.0, 0.0],
            first,
            second,
            phi: [[value, 0.0, 0.0], [0.0; 3], [0.0; 3]],
        }
    }

    #[test]
    fn constructor_rejects_inconsistent_input() {
        let lattice = LatticeVectors::cubic(1.0);
        let good = vec![
            block(0, 0, 0.0, 1.0),
            block(0, 0, 1.0, -0.5),
            block(0, 0, -1.0, -0.5),
        ];
        assert!(DynamicalMatrix::new(vec![1.0], good.clone(), lattice.clone()).is_ok());

        let new = |masses: Vec<f64>, blocks: Vec<ForceConstantBlock>| {
            DynamicalMatrix::new(masses, blocks, lattice.clone())
        };
        // masses
        assert!(new(vec![], good.clone()).is_err());
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(new(vec![bad], good.clone()).is_err(), "mass {bad}");
        }
        // basis index out of range
        let mut out_of_range = good.clone();
        out_of_range.push(block(0, 3, 0.0, 0.0));
        assert!(new(vec![1.0], out_of_range).is_err());
        // non-finite numbers
        let mut not_finite = good.clone();
        not_finite[0].phi[1][1] = f64::NAN;
        assert!(new(vec![1.0], not_finite).is_err());
        let mut far = good.clone();
        far[1].translation[2] = f64::INFINITY;
        assert!(new(vec![1.0], far).is_err());
        // a coupling without its transposed partner has no real spectrum
        assert!(new(vec![1.0], vec![block(0, 0, 1.0, -0.5)]).is_err());
        // a partner at the right place but with the wrong sign
        let mut wrong = good.clone();
        wrong[2].phi[0][0] = 0.5;
        assert!(new(vec![1.0], wrong).is_err());
        // an asymmetric on-site block is its own partner and must be symmetric
        let mut skew = good;
        skew[0].phi[0][1] = 0.3;
        assert!(new(vec![1.0], skew).is_err());
    }

    #[test]
    fn invalid_masses_in_hand_built_structs_give_nan_instead_of_panicking() {
        let mut dm = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(1.0));
        dm.masses[0] = 0.0;
        assert!(dm.matrix_at([0.1, 0.2, 0.3]).is_none());
        let omega = dm.at_q([0.1, 0.2, 0.3]);
        assert_eq!(omega.len(), 3);
        assert!(omega.iter().all(|w| w.is_nan()));
        assert!(!dm.acoustic_at_gamma());
        // Blocks that name a missing basis atom are skipped.
        let mut stray = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(1.0));
        stray.force_constants.push(block(0, 5, 0.0, 9.0));
        let reference = DynamicalMatrix::new_cubic(K, M, 6, LatticeVectors::cubic(1.0));
        assert_close(
            &stray.at_q([0.4, 0.5, 0.6]),
            &reference.at_q([0.4, 0.5, 0.6]),
            1.0e-12,
            "stray block",
        );
    }
}
