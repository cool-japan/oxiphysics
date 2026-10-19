// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Post-Hartree-Fock electron correlation on a closed-shell canonical RHF
//! reference: MP2 (with spin-component split), SCS-MP2, singlet CIS with
//! oscillator strengths, and the closed-shell CCS energy.
//!
//! All quantities are in atomic units (Hartree, Bohr).  Spatial molecular
//! orbitals are indexed `0..n_mo`, the first `n_occ` doubly occupied, ordered by
//! ascending orbital energy.  Two-electron integrals use chemists' notation
//! (pq|rs) = ∫∫ φ_p(1) φ_q(1) r₁₂⁻¹ φ_r(2) φ_s(2).

use super::error::QuantumChemistryError;
use super::integrals::dipole_integral;
use super::types::BasisFunction;
use super::types_2::{HartreeFockSolver, Mp2Result};
use super::types_3::{CisState, HartreeFockResult, ScfStatus};
use oxiphysics_core::numerical_methods::symmetric_eigen_n;

type QcResult<T> = Result<T, QuantumChemistryError>;

/// Electron correlation methods: MP2, CIS, and coupled cluster singles (CCS).
#[derive(Debug, Clone)]
pub struct ElectronCorrelation {
    /// Number of doubly occupied spatial orbitals.
    pub n_occ: usize,
    /// Number of virtual spatial orbitals.
    pub n_virt: usize,
    /// Canonical orbital energies ε_p \[Hartree\], length `n_occ + n_virt`.
    pub orbital_energies: Vec<f64>,
    /// Full MO-basis two-electron integral tensor (pq|rs), row-major with index
    /// `((p·n + q)·n + r)·n + s`, `n = n_occ + n_virt`.
    pub mo_eri: Vec<f64>,
    /// HF reference energy \[Hartree\].
    pub e_hf: f64,
    /// MO-basis dipole integrals ⟨p| r_k |q⟩ for k = x, y, z, each row-major
    /// `n × n` \[Bohr\].  Required by [`cis_states`](Self::cis_states).
    pub mo_dipole: Option<[Vec<f64>; 3]>,
}

impl ElectronCorrelation {
    /// Create an electron-correlation solver from canonical orbital energies and
    /// the full MO two-electron tensor (see [`mo_eri`](Self::mo_eri) for the layout).
    ///
    /// # Errors
    /// [`QuantumChemistryError::DimensionMismatch`] if `orbital_energies` does not
    /// have `n_occ + n_virt` entries or `mo_eri` does not have `(n_occ + n_virt)⁴`.
    pub fn new(
        n_occ: usize,
        n_virt: usize,
        orbital_energies: Vec<f64>,
        mo_eri: Vec<f64>,
        e_hf: f64,
    ) -> QcResult<Self> {
        let ec = Self {
            n_occ,
            n_virt,
            orbital_energies,
            mo_eri,
            e_hf,
            mo_dipole: None,
        };
        ec.validate()?;
        Ok(ec)
    }

    /// Attach MO-basis dipole integrals ⟨p| r_k |q⟩ (each row-major `n × n`).
    ///
    /// # Errors
    /// [`QuantumChemistryError::DimensionMismatch`] if a component has the wrong length.
    pub fn with_mo_dipoles(mut self, dipoles: [Vec<f64>; 3]) -> QcResult<Self> {
        let n = self.n_mo();
        for component in &dipoles {
            if component.len() != n * n {
                return Err(QuantumChemistryError::DimensionMismatch {
                    what: "MO dipole integrals",
                    expected: n * n,
                    found: component.len(),
                });
            }
        }
        self.mo_dipole = Some(dipoles);
        Ok(self)
    }

    /// Build the correlation solver from a converged closed-shell RHF calculation.
    ///
    /// The AO integrals (μν|λσ) of `solver` are transformed to the MO basis with
    /// four quarter transformations (O(n⁵)) using `result.mo_coefficients`.  When
    /// the solver carries its basis functions (it was built with
    /// [`HartreeFockSolver::from_atoms`] or [`HartreeFockSolver::from_basis`]), the
    /// AO dipole integrals are transformed as well (enabling CIS oscillator
    /// strengths); otherwise `mo_dipole` is left `None` and can be supplied with
    /// [`with_mo_dipoles`](Self::with_mo_dipoles).
    ///
    /// # Errors
    /// * [`QuantumChemistryError::ScfNotConverged`] if `result.status` is not converged.
    /// * [`QuantumChemistryError::OpenShellReference`] for an odd electron count.
    /// * [`QuantumChemistryError::DimensionMismatch`] if the coefficient matrix,
    ///   orbital energies or AO integrals do not match `solver.n_basis`.
    pub fn from_hartree_fock(
        solver: &HartreeFockSolver,
        result: &HartreeFockResult,
    ) -> QcResult<Self> {
        if result.status != ScfStatus::Converged {
            return Err(QuantumChemistryError::ScfNotConverged);
        }
        if !solver.n_electrons.is_multiple_of(2) {
            return Err(QuantumChemistryError::OpenShellReference(
                solver.n_electrons,
            ));
        }
        let n = solver.n_basis;
        let n_occ = solver.n_electrons / 2;
        check_len("MO coefficient rows", n, result.mo_coefficients.len())?;
        for row in &result.mo_coefficients {
            check_len("MO coefficient columns", n, row.len())?;
        }
        check_len("orbital energies", n, result.orbital_energies.len())?;
        if n_occ > n {
            return Err(QuantumChemistryError::DimensionMismatch {
                what: "occupied orbitals vs basis functions",
                expected: n,
                found: n_occ,
            });
        }
        let n_pair = n * (n + 1) / 2;
        check_len(
            "packed AO ERIs",
            n_pair * (n_pair + 1) / 2,
            solver.eri.len(),
        )?;

        // Full AO tensor (μν|λσ).
        let mut ao = vec![0.0_f64; n * n * n * n];
        for mu in 0..n {
            for nu in 0..n {
                for la in 0..n {
                    for si in 0..n {
                        ao[((mu * n + nu) * n + la) * n + si] =
                            solver.eri[solver.eri_index(mu, nu, la, si)];
                    }
                }
            }
        }
        let mo_eri = ao_to_mo_eri(&ao, &result.mo_coefficients, n);

        let mut ec = Self::new(
            n_occ,
            n - n_occ,
            result.orbital_energies.clone(),
            mo_eri,
            result.energy,
        )?;
        if solver.basis_functions.len() == n
            && let Ok(ao_dip) = ao_dipole_integrals(&solver.basis_functions, [0.0; 3])
        {
            let c = &result.mo_coefficients;
            let dipoles = [
                ao_to_mo_one_electron(&ao_dip[0], c, n),
                ao_to_mo_one_electron(&ao_dip[1], c, n),
                ao_to_mo_one_electron(&ao_dip[2], c, n),
            ];
            ec = ec.with_mo_dipoles(dipoles)?;
        }
        Ok(ec)
    }

    /// Number of spatial molecular orbitals.
    pub fn n_mo(&self) -> usize {
        self.n_occ + self.n_virt
    }

    /// MO integral (pq|rs).  Indices must be below [`n_mo`](Self::n_mo); call
    /// only after validation (all public methods validate first).
    fn eri(&self, p: usize, q: usize, r: usize, s: usize) -> f64 {
        let n = self.n_mo();
        self.mo_eri[((p * n + q) * n + r) * n + s]
    }

    /// Check that all stored arrays have consistent sizes.
    fn validate(&self) -> QcResult<()> {
        let n = self.n_mo();
        check_len("orbital energies", n, self.orbital_energies.len())?;
        check_len(
            "MO two-electron integrals",
            n * n * n * n,
            self.mo_eri.len(),
        )?;
        if let Some(d) = &self.mo_dipole {
            for component in d {
                check_len("MO dipole integrals", n * n, component.len())?;
            }
        }
        Ok(())
    }

    /// Second-order Møller–Plesset correlation energy (closed shell).
    ///
    /// E_MP2 = −Σ_{ijab} (ia|jb) [2(ia|jb) − (ib|ja)] / (ε_a + ε_b − ε_i − ε_j)
    ///
    /// split into opposite-spin E_OS = −Σ (ia|jb)² / Δ and same-spin
    /// E_SS = −Σ (ia|jb)[(ia|jb) − (ib|ja)] / Δ (E_MP2 = E_OS + E_SS).
    ///
    /// # Errors
    /// Inconsistent dimensions, or a non-positive denominator (non-aufbau reference).
    pub fn mp2_energy(&self) -> QcResult<Mp2Result> {
        self.validate()?;
        let o = self.n_occ;
        let n = self.n_mo();
        let eps = &self.orbital_energies;
        let mut e_os = 0.0;
        let mut e_ss = 0.0;
        for i in 0..o {
            for j in 0..o {
                for a in o..n {
                    for b in o..n {
                        let denom = eps[a] + eps[b] - eps[i] - eps[j];
                        if denom <= 0.0 {
                            return Err(QuantumChemistryError::NonPositiveDenominator(denom));
                        }
                        let iajb = self.eri(i, a, j, b);
                        let ibja = self.eri(i, b, j, a);
                        e_os -= iajb * iajb / denom;
                        e_ss -= iajb * (iajb - ibja) / denom;
                    }
                }
            }
        }
        let e_mp2 = e_os + e_ss;
        Ok(Mp2Result {
            e_hf: self.e_hf,
            e_mp2,
            e_total: self.e_hf + e_mp2,
            e_ss,
            e_os,
        })
    }

    /// Spin-component-scaled MP2 total energy E_HF + c_os E_OS + c_ss E_SS
    /// (Grimme 2003 uses c_os = 6/5, c_ss = 1/3).
    ///
    /// # Errors
    /// Same as [`mp2_energy`](Self::mp2_energy).
    pub fn scs_mp2_energy(&self, c_os: f64, c_ss: f64) -> QcResult<f64> {
        let mp2 = self.mp2_energy()?;
        Ok(self.e_hf + c_os * mp2.e_os + c_ss * mp2.e_ss)
    }

    /// Lowest `n_states` singlet CIS excited states.
    ///
    /// Diagonalises A_{ia,jb} = δ_ij δ_ab (ε_a − ε_i) + 2(ia|jb) − (ij|ab) for
    /// spin-adapted singlet amplitudes X_ia (Σ X² = 1).  The transition dipole is
    /// ⟨0|μ|S⟩ = √2 Σ_ia X_ia ⟨i|r|a⟩ and the oscillator strength
    /// f = ⅔ ΔE |⟨0|μ|S⟩|² (length gauge).
    ///
    /// # Errors
    /// Inconsistent dimensions, or [`QuantumChemistryError::MissingData`] when no
    /// MO dipole integrals are available.
    pub fn cis_states(&self, n_states: usize) -> QcResult<Vec<CisState>> {
        self.validate()?;
        let dipoles = self
            .mo_dipole
            .as_ref()
            .ok_or(QuantumChemistryError::MissingData("MO dipole integrals"))?;
        let o = self.n_occ;
        let v = self.n_virt;
        let n = self.n_mo();
        let dim = o * v;
        if dim == 0 || n_states == 0 {
            return Ok(Vec::new());
        }
        let eps = &self.orbital_energies;
        let mut a_mat = vec![vec![0.0_f64; dim]; dim];
        for i in 0..o {
            for a in 0..v {
                let ia = i * v + a;
                for j in 0..o {
                    for b in 0..v {
                        let jb = j * v + b;
                        let mut val =
                            2.0 * self.eri(i, o + a, j, o + b) - self.eri(i, j, o + a, o + b);
                        if ia == jb {
                            val += eps[o + a] - eps[i];
                        }
                        a_mat[ia][jb] = val;
                    }
                }
            }
        }
        let (evals, evecs) = symmetric_eigen_n(&a_mat);
        let mut states = Vec::with_capacity(n_states.min(dim));
        for (energy, x) in evals.iter().zip(evecs.iter()).take(n_states) {
            let mut transition = [0.0_f64; 3];
            for (k, comp) in dipoles.iter().enumerate() {
                let mut t = 0.0;
                for i in 0..o {
                    for a in 0..v {
                        t += x[i * v + a] * comp[i * n + o + a];
                    }
                }
                transition[k] = std::f64::consts::SQRT_2 * t;
            }
            let t2 = transition.iter().map(|t| t * t).sum::<f64>();
            let (dominant, _) = x.iter().enumerate().fold((0, 0.0_f64), |best, (idx, c)| {
                if c.abs() > best.1 {
                    (idx, c.abs())
                } else {
                    best
                }
            });
            states.push(CisState {
                excitation_energy: *energy,
                oscillator_strength: 2.0 / 3.0 * energy * t2,
                dominant_occ: dominant / v,
                dominant_virt: dominant % v,
                amplitudes: x.clone(),
            });
        }
        Ok(states)
    }

    /// T1 diagnostic of Lee & Taylor, √(Σ_ia t_ia² / n_occ) for closed-shell
    /// spatial amplitudes (equivalently ‖t‖ / √N_el with spin-orbital amplitudes).
    pub fn ccs_t1_diagnostic(&self, t1_amplitudes: &[f64]) -> f64 {
        if t1_amplitudes.is_empty() {
            return 0.0;
        }
        let norm_sq: f64 = t1_amplitudes.iter().map(|t| t * t).sum();
        let n_e = (self.n_occ as f64).max(1.0);
        (norm_sq / n_e).sqrt()
    }

    /// Closed-shell coupled-cluster-singles correlation energy for amplitudes
    /// `t1[i·n_virt + a]` on a canonical HF reference (f_ia = 0):
    ///
    /// E_CCS = Σ_{ia,jb} [2(ia|jb) − (ib|ja)] t_ia t_jb
    ///
    /// # Errors
    /// [`QuantumChemistryError::DimensionMismatch`] if `t1_amplitudes` does not
    /// have `n_occ · n_virt` entries.
    pub fn ccs_energy(&self, t1_amplitudes: &[f64]) -> QcResult<f64> {
        self.validate()?;
        let o = self.n_occ;
        let v = self.n_virt;
        check_len("T1 amplitudes", o * v, t1_amplitudes.len())?;
        let mut e = 0.0;
        for i in 0..o {
            for a in 0..v {
                let t_ia = t1_amplitudes[i * v + a];
                for j in 0..o {
                    for b in 0..v {
                        let t_jb = t1_amplitudes[j * v + b];
                        let k = 2.0 * self.eri(i, o + a, j, o + b) - self.eri(i, o + b, j, o + a);
                        e += k * t_ia * t_jb;
                    }
                }
            }
        }
        Ok(e)
    }
}

/// Return a dimension-mismatch error unless `found == expected`.
fn check_len(what: &'static str, expected: usize, found: usize) -> QcResult<()> {
    if expected == found {
        Ok(())
    } else {
        Err(QuantumChemistryError::DimensionMismatch {
            what,
            expected,
            found,
        })
    }
}

/// Transform a full AO tensor (μν|λσ) to the MO basis (pq|rs) with four
/// quarter transformations; `c[p][μ]` = C_μp.
fn ao_to_mo_eri(ao: &[f64], c: &[Vec<f64>], n: usize) -> Vec<f64> {
    // Transform the index at `position` (0 = outermost) of a row-major n⁴ tensor.
    let transform = |input: &[f64], position: u32| -> Vec<f64> {
        let mut out = vec![0.0_f64; n * n * n * n];
        let stride = n.pow(3 - position);
        for (idx, slot) in out.iter_mut().enumerate() {
            let p = (idx / stride) % n;
            let base = idx - p * stride;
            let mut acc = 0.0;
            for (mu, &c_mu) in c[p].iter().enumerate().take(n) {
                acc += c_mu * input[base + mu * stride];
            }
            *slot = acc;
        }
        out
    };
    let t1 = transform(ao, 0);
    let t2 = transform(&t1, 1);
    let t3 = transform(&t2, 2);
    transform(&t3, 3)
}

/// Transform a one-electron AO matrix to the MO basis: M_pq = Σ C_μp A_μν C_νq,
/// returned row-major `n × n`.
fn ao_to_mo_one_electron(ao: &[Vec<f64>], c: &[Vec<f64>], n: usize) -> Vec<f64> {
    let mut out = vec![0.0_f64; n * n];
    for p in 0..n {
        for q in 0..n {
            let mut acc = 0.0;
            for mu in 0..n {
                for nu in 0..n {
                    acc += c[p][mu] * ao[mu][nu] * c[q][nu];
                }
            }
            out[p * n + q] = acc;
        }
    }
    out
}

/// AO dipole integrals ⟨μ| (r − origin)_k |ν⟩ \[Bohr\] for contracted Cartesian
/// Gaussians of any angular momentum (McMurchie–Davidson, see
/// [`dipole_integral`]).
///
/// # Errors
/// Never fails for a well-formed basis; the `Result` is kept for API stability.
pub fn ao_dipole_integrals(
    basis: &[BasisFunction],
    origin: [f64; 3],
) -> QcResult<[Vec<Vec<f64>>; 3]> {
    let n = basis.len();
    let mut d = [
        vec![vec![0.0_f64; n]; n],
        vec![vec![0.0_f64; n]; n],
        vec![vec![0.0_f64; n]; n],
    ];
    for i in 0..n {
        for j in i..n {
            let acc = dipole_integral(&basis[i], &basis[j], origin);
            for k in 0..3 {
                d[k][i][j] = acc[k];
                d[k][j][i] = acc[k];
            }
        }
    }
    Ok(d)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// H₂ at R = 1.4 bohr with ζ = 1.24 (Szabo & Ostlund, ch. 3 and 6).
    fn h2_reference() -> (HartreeFockSolver, HartreeFockResult) {
        let r_ang = 1.4 * super::super::functions::BOHR_TO_ANGSTROM;
        let solver = HartreeFockSolver::from_atoms(&[1, 1], &[[0.0; 3], [0.0, 0.0, r_ang]], 0)
            .expect("H2 solver");
        let result = solver.run_scf();
        assert_eq!(result.status, ScfStatus::Converged);
        (solver, result)
    }

    #[test]
    fn h2_hf_reference_matches_szabo_ostlund() {
        let (solver, result) = h2_reference();
        assert!((solver.overlap[0][1] - 0.6593).abs() < 1e-4);
        assert!(
            (result.energy + 1.1167).abs() < 5e-4,
            "E_HF {}",
            result.energy
        );
        assert!((result.orbital_energies[0] + 0.578).abs() < 2e-3);
        assert!((result.orbital_energies[1] - 0.670).abs() < 2e-3);
    }

    #[test]
    fn h2_mo_integrals_match_szabo_ostlund() {
        let (solver, result) = h2_reference();
        let ec = ElectronCorrelation::from_hartree_fock(&solver, &result).expect("MP2 setup");
        assert_eq!((ec.n_occ, ec.n_virt), (1, 1));
        assert!((ec.eri(0, 0, 0, 0) - 0.6746).abs() < 1e-3);
        assert!((ec.eri(0, 0, 1, 1) - 0.6636).abs() < 1e-3);
        assert!((ec.eri(0, 1, 0, 1) - 0.1813).abs() < 1e-3);
        assert!((ec.eri(1, 1, 1, 1) - 0.6975).abs() < 1e-3);
        // 8-fold permutational symmetry survives the transformation.
        assert!((ec.eri(0, 1, 0, 1) - ec.eri(1, 0, 1, 0)).abs() < 1e-12);
        assert!((ec.eri(0, 0, 1, 1) - ec.eri(1, 1, 0, 0)).abs() < 1e-12);
    }

    #[test]
    fn h2_mp2_energy_is_negative_and_matches_reference() {
        let (solver, result) = h2_reference();
        let ec = ElectronCorrelation::from_hartree_fock(&solver, &result).expect("MP2 setup");
        let mp2 = ec.mp2_energy().expect("MP2");
        assert!(mp2.e_mp2 < 0.0, "E_MP2 must be negative: {}", mp2.e_mp2);
        assert!((mp2.e_mp2 + 0.0132).abs() < 3e-4, "E_MP2 {}", mp2.e_mp2);
        // Two-orbital closed form: K² / (2(ε₁ − ε₂)).
        let k = ec.eri(0, 1, 0, 1);
        let eps = &ec.orbital_energies;
        let closed = k * k / (2.0 * (eps[0] - eps[1]));
        assert!((mp2.e_mp2 - closed).abs() < 1e-12);
        // One occupied orbital: no same-spin pairs.
        assert!(mp2.e_ss.abs() < 1e-14);
        assert!((mp2.e_os - mp2.e_mp2).abs() < 1e-14);
        assert!((mp2.e_total - mp2.e_hf - mp2.e_mp2).abs() < 1e-14);
        // SCS-MP2 with unit scaling reproduces MP2.
        let scs = ec.scs_mp2_energy(1.0, 1.0).expect("SCS");
        assert!((scs - mp2.e_total).abs() < 1e-14);
    }

    #[test]
    fn h2_cis_singlet_and_oscillator_strength() {
        let (solver, result) = h2_reference();
        let ec = ElectronCorrelation::from_hartree_fock(&solver, &result).expect("CIS setup");
        let states = ec.cis_states(3).expect("CIS");
        assert_eq!(states.len(), 1);
        let s = &states[0];
        let eps = &ec.orbital_energies;
        let expected = eps[1] - eps[0] + 2.0 * ec.eri(0, 1, 0, 1) - ec.eri(0, 0, 1, 1);
        assert!((s.excitation_energy - expected).abs() < 1e-12);
        assert!((s.excitation_energy - 0.947).abs() < 5e-3);
        // Minimal-basis transition dipole along the bond:
        // ⟨σg|z|σu⟩ = (z_B − z_A) / (2 √(1 − S²)).
        let s12 = solver.overlap[0][1];
        let mu = 1.4 / (2.0 * (1.0 - s12 * s12).sqrt());
        let f_expected = 4.0 / 3.0 * s.excitation_energy * mu * mu;
        assert!(
            (s.oscillator_strength - f_expected).abs() < 1e-6,
            "f {} vs {f_expected}",
            s.oscillator_strength
        );
        assert_eq!((s.dominant_occ, s.dominant_virt), (0, 0));
    }

    #[test]
    fn ccs_energy_uses_exchange_corrected_integrals() {
        let (solver, result) = h2_reference();
        let ec = ElectronCorrelation::from_hartree_fock(&solver, &result).expect("setup");
        assert_eq!(ec.ccs_energy(&[0.0]).expect("CCS"), 0.0);
        let t = 0.1;
        let e = ec.ccs_energy(&[t]).expect("CCS");
        // [2(gu|gu) − (gu|gu)] t² = K t²
        assert!((e - ec.eri(0, 1, 0, 1) * t * t).abs() < 1e-14);
        assert!(matches!(
            ec.ccs_energy(&[0.1, 0.2]),
            Err(QuantumChemistryError::DimensionMismatch { .. })
        ));
    }

    #[test]
    fn missing_or_inconsistent_inputs_are_errors() {
        assert!(matches!(
            ElectronCorrelation::new(2, 2, vec![-0.6, -0.5, 0.3, 0.4], vec![0.1; 4], -1.0),
            Err(QuantumChemistryError::DimensionMismatch { .. })
        ));
        assert!(matches!(
            ElectronCorrelation::new(1, 1, vec![-0.5], vec![0.0; 16], -1.0),
            Err(QuantumChemistryError::DimensionMismatch { .. })
        ));
        let ec = ElectronCorrelation::new(1, 1, vec![-0.5, 0.5], vec![0.1; 16], -1.0)
            .expect("valid sizes");
        assert_eq!(
            ec.cis_states(1).err(),
            Some(QuantumChemistryError::MissingData("MO dipole integrals"))
        );
    }

    #[test]
    fn non_aufbau_reference_is_rejected() {
        let ec = ElectronCorrelation::new(1, 1, vec![0.5, -0.5], vec![0.1; 16], -1.0)
            .expect("valid sizes");
        assert!(matches!(
            ec.mp2_energy(),
            Err(QuantumChemistryError::NonPositiveDenominator(_))
        ));
    }

    #[test]
    fn unconverged_reference_is_rejected() {
        let r_ang = 1.4 * super::super::functions::BOHR_TO_ANGSTROM;
        let mut solver = HartreeFockSolver::from_atoms(&[1, 1], &[[0.0; 3], [0.0, 0.0, r_ang]], 0)
            .expect("H2 solver");
        solver.max_iter = 1;
        let result = solver.run_scf();
        assert_eq!(result.status, ScfStatus::NotConverged);
        assert_eq!(
            ElectronCorrelation::from_hartree_fock(&solver, &result).err(),
            Some(QuantumChemistryError::ScfNotConverged)
        );
    }

    #[test]
    fn dipole_integrals_of_s_functions() {
        let a = BasisFunction::sto3g_s([0.0, 0.0, 0.0], 1.24);
        let b = BasisFunction::sto3g_s([0.0, 0.0, 1.4], 1.24);
        let d = ao_dipole_integrals(&[a.clone(), b.clone()], [0.0; 3]).expect("s dipoles");
        // ⟨A|z|A⟩ = 0, ⟨B|z|B⟩ = 1.4, ⟨A|z|B⟩ = S_AB · 0.7 (equal exponents).
        assert!(d[2][0][0].abs() < 1e-12);
        assert!((d[2][1][1] - 1.4).abs() < 1e-6);
        assert!((d[2][0][1] - a.overlap_with(&b) * 0.7).abs() < 1e-10);
        assert!(d[0][0][1].abs() < 1e-14);
        // A p_z function on the same centre as `a`: only ⟨s|z|p_z⟩ survives,
        // and the p_z density is symmetric about the origin.
        let p = BasisFunction::sto3g_p([0.0; 3], 1.72, super::super::types::CartesianAxis::Z);
        let d = ao_dipole_integrals(&[a, p], [0.0; 3]).expect("s/p dipoles");
        assert!(d[2][0][1].abs() > 0.1, "⟨s|z|p_z⟩ = {}", d[2][0][1]);
        assert!(d[0][0][1].abs() < 1e-14 && d[1][0][1].abs() < 1e-14);
        assert!(d[2][1][1].abs() < 1e-14);
        assert!((d[2][0][1] - d[2][1][0]).abs() < 1e-15);
    }
}
