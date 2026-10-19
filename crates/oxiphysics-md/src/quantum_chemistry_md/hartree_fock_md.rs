// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hartree-Fock Born-Oppenheimer molecular dynamics.
//!
//! The electronic energy at every geometry is the converged restricted
//! Hartree-Fock/STO-3G energy of [`HartreeFockSolver`].  Nuclear forces are the
//! analytic nuclear-repulsion forces plus central finite differences of the
//! electronic energy (analytic HF gradients are not available in the solver).
//! All quantities are in Hartree atomic units: lengths in Bohr, energies in
//! Hartree, masses in electron masses and time in ħ/E_h ≈ 24.19 as.

use super::dist;
use super::error::QuantumChemistryMdError;
use crate::quantum_chemistry::{HartreeFockSolver, ScfStatus};

/// Bohr → Ångström conversion (inverse of the factor used by
/// [`HartreeFockSolver::from_atoms`]).
const BOHR_PER_ANGSTROM: f64 = 1.889_725_988_6;

/// Conversion factor from unified atomic mass units to electron masses.
pub const AMU_TO_ELECTRON_MASS: f64 = 1_822.888_486;

/// Forces cached at the end of a time step, tagged with the configuration
/// they belong to so that a later change of the public fields invalidates them.
#[derive(Debug, Clone)]
struct ForceCache {
    positions: Vec<[f64; 3]>,
    nuclear_charges: Vec<f64>,
    total_charge: i32,
    forces: Vec<[f64; 3]>,
}

/// Hartree-Fock Born-Oppenheimer molecular dynamics driver.
///
/// At each geometry a restricted Hartree-Fock/STO-3G SCF calculation is
/// converged; the nuclei move on the resulting Born-Oppenheimer surface with
/// the velocity-Verlet integrator.  Supported elements are those with an
/// STO-3G basis in [`HartreeFockSolver::from_atoms`] (H, C, N, O, F, S), and
/// the electron count must be even (closed shell).
#[derive(Debug, Clone)]
pub struct HartreeFockMd {
    /// Nuclear charges Z_I for each atom (must be integral atomic numbers).
    pub nuclear_charges: Vec<f64>,
    /// Cartesian positions of nuclei \[Bohr\].
    pub positions: Vec<[f64; 3]>,
    /// Velocities of nuclei \[Bohr / a.u. time\].
    pub velocities: Vec<[f64; 3]>,
    /// Nuclear masses \[electron masses\] (multiply amu by [`AMU_TO_ELECTRON_MASS`]).
    pub masses: Vec<f64>,
    /// Net molecular charge (0 for a neutral molecule).
    pub total_charge: i32,
    /// SCF energy convergence threshold \[Hartree\]; it must be far below
    /// `fd_step` × the force accuracy required.
    pub scf_energy_threshold: f64,
    /// Maximum number of SCF iterations per energy evaluation.
    pub scf_max_iter: usize,
    /// Central finite-difference step for the electronic forces \[Bohr\].
    pub fd_step: f64,
    cache: Option<ForceCache>,
}

impl HartreeFockMd {
    /// Create a new Hartree-Fock MD driver for a neutral molecule.
    ///
    /// # Arguments
    /// * `nuclear_charges` – Atomic numbers as floats.
    /// * `positions` – Initial nuclear positions in Bohr.
    /// * `masses` – Nuclear masses in electron masses.
    ///
    /// # Errors
    /// [`QuantumChemistryMdError::DimensionMismatch`] if the three inputs differ
    /// in length, [`QuantumChemistryMdError::UnsupportedNuclearCharge`] for a
    /// charge that is not a supported atomic number, and
    /// [`QuantumChemistryMdError::InvalidParameter`] for a non-positive or
    /// non-finite mass or a non-finite coordinate.
    pub fn new(
        nuclear_charges: Vec<f64>,
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
    ) -> Result<Self, QuantumChemistryMdError> {
        let n = positions.len();
        let md = Self {
            nuclear_charges,
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
            total_charge: 0,
            scf_energy_threshold: 1e-11,
            scf_max_iter: 200,
            fd_step: 1e-4,
            cache: None,
        };
        md.validate()?;
        Ok(md)
    }

    /// Set the net molecular charge (builder style).
    pub fn with_total_charge(mut self, total_charge: i32) -> Self {
        self.total_charge = total_charge;
        self.cache = None;
        self
    }

    /// Converged Hartree-Fock electronic energy in Hartree (nuclear repulsion
    /// excluded).  Zero when the system has no electrons.
    ///
    /// # Errors
    /// Invalid configuration, unsupported element, open shell or SCF
    /// non-convergence.
    pub fn scf_energy(&self) -> Result<f64, QuantumChemistryMdError> {
        self.validate()?;
        let atomic_numbers = self.atomic_numbers()?;
        let n_electrons = self.electron_count(&atomic_numbers)?;
        self.electronic_energy_at(&atomic_numbers, n_electrons, &self.positions)
    }

    /// Compute the nuclear repulsion energy Σ_{I<J} Z_I Z_J / R_IJ in Hartree.
    pub fn nuclear_repulsion(&self) -> f64 {
        let n = self.positions.len().min(self.nuclear_charges.len());
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.positions[i], &self.positions[j]).max(1e-6);
                e += self.nuclear_charges[i] * self.nuclear_charges[j] / r;
            }
        }
        e
    }

    /// Nuclear forces F_I = −∂E_BO/∂R_I in Hartree/Bohr.
    ///
    /// The nuclear-repulsion part is analytic; the electronic part is the
    /// central finite difference −(E_el(R + h) − E_el(R − h)) / 2h with
    /// h = [`fd_step`](Self::fd_step), costing 6N SCF calculations.
    ///
    /// # Errors
    /// As for [`scf_energy`](Self::scf_energy), at any displaced geometry.
    pub fn forces(&self) -> Result<Vec<[f64; 3]>, QuantumChemistryMdError> {
        self.validate()?;
        let atomic_numbers = self.atomic_numbers()?;
        let n_electrons = self.electron_count(&atomic_numbers)?;
        self.forces_at(&atomic_numbers, n_electrons, &self.positions)
    }

    /// Energy gradient ∂E_BO/∂R_I in Hartree/Bohr (the negative of [`forces`](Self::forces)).
    ///
    /// # Errors
    /// As for [`forces`](Self::forces).
    pub fn gradient(&self) -> Result<Vec<[f64; 3]>, QuantumChemistryMdError> {
        Ok(self
            .forces()?
            .into_iter()
            .map(|f| [-f[0], -f[1], -f[2]])
            .collect())
    }

    /// Advance the system by one velocity-Verlet time step.
    ///
    /// v(t+dt/2) = v(t) + dt F(t) / 2m, R(t+dt) = R(t) + dt v(t+dt/2),
    /// v(t+dt) = v(t+dt/2) + dt F(t+dt) / 2m.
    ///
    /// The forces at the end of the step are cached and reused at the start of
    /// the next one (one force evaluation per step) as long as the positions,
    /// charges and total charge are unchanged.  On error the state is left
    /// untouched.
    ///
    /// # Arguments
    /// * `dt` – Time step in atomic units (1 a.u. ≈ 24.2 as).
    ///
    /// # Errors
    /// [`QuantumChemistryMdError::InvalidParameter`] for a non-positive or
    /// non-finite `dt`, otherwise as for [`forces`](Self::forces).
    pub fn step(&mut self, dt: f64) -> Result<(), QuantumChemistryMdError> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "dt",
                value: dt,
            });
        }
        self.validate()?;
        let atomic_numbers = self.atomic_numbers()?;
        let n_electrons = self.electron_count(&atomic_numbers)?;
        let f_start = match self.cached_forces() {
            Some(f) => f,
            None => self.forces_at(&atomic_numbers, n_electrons, &self.positions)?,
        };
        let mut velocities = self.velocities.clone();
        let mut positions = self.positions.clone();
        for ((vel, pos), (f, &m)) in velocities
            .iter_mut()
            .zip(positions.iter_mut())
            .zip(f_start.iter().zip(self.masses.iter()))
        {
            for k in 0..3 {
                vel[k] += 0.5 * dt * f[k] / m;
                pos[k] += dt * vel[k];
            }
        }
        let f_end = self.forces_at(&atomic_numbers, n_electrons, &positions)?;
        for ((vel, f), &m) in velocities
            .iter_mut()
            .zip(f_end.iter())
            .zip(self.masses.iter())
        {
            for k in 0..3 {
                vel[k] += 0.5 * dt * f[k] / m;
            }
        }
        self.positions = positions;
        self.velocities = velocities;
        self.cache = Some(ForceCache {
            positions: self.positions.clone(),
            nuclear_charges: self.nuclear_charges.clone(),
            total_charge: self.total_charge,
            forces: f_end,
        });
        Ok(())
    }

    /// Nuclear kinetic energy Σ_I ½ m_I |v_I|² in Hartree.
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities
            .iter()
            .zip(self.masses.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }

    /// Total Born-Oppenheimer energy (SCF electronic + nuclear repulsion) in Hartree.
    ///
    /// # Errors
    /// As for [`scf_energy`](Self::scf_energy).
    pub fn total_energy(&self) -> Result<f64, QuantumChemistryMdError> {
        Ok(self.scf_energy()? + self.nuclear_repulsion())
    }

    // -- private helpers --

    fn validate(&self) -> Result<(), QuantumChemistryMdError> {
        let n = self.positions.len();
        for (what, found) in [
            ("nuclear_charges", self.nuclear_charges.len()),
            ("masses", self.masses.len()),
            ("velocities", self.velocities.len()),
        ] {
            if found != n {
                return Err(QuantumChemistryMdError::DimensionMismatch {
                    what,
                    expected: n,
                    found,
                });
            }
        }
        if let Some(&m) = self.masses.iter().find(|m| !(m.is_finite() && **m > 0.0)) {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "mass",
                value: m,
            });
        }
        if let Some(x) = self
            .positions
            .iter()
            .flat_map(|p| p.iter())
            .find(|x| !x.is_finite())
        {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "position",
                value: *x,
            });
        }
        if !(self.fd_step.is_finite() && self.fd_step > 0.0) {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "fd_step",
                value: self.fd_step,
            });
        }
        self.atomic_numbers().map(|_| ())
    }

    /// Integral atomic numbers with an STO-3G basis.
    fn atomic_numbers(&self) -> Result<Vec<u32>, QuantumChemistryMdError> {
        self.nuclear_charges
            .iter()
            .enumerate()
            .map(|(index, &charge)| {
                let rounded = charge.round();
                let supported = charge.is_finite()
                    && (charge - rounded).abs() < 1e-9
                    && (1.0..=118.0).contains(&rounded)
                    && HartreeFockSolver::from_atoms(&[rounded as u32], &[[0.0; 3]], 0).is_some();
                if supported {
                    Ok(rounded as u32)
                } else {
                    Err(QuantumChemistryMdError::UnsupportedNuclearCharge { index, charge })
                }
            })
            .collect()
    }

    /// Number of electrons (neutral-atom count minus the total charge),
    /// validated as non-negative and even.
    fn electron_count(&self, atomic_numbers: &[u32]) -> Result<usize, QuantumChemistryMdError> {
        let mut neutral: i64 = 0;
        for (index, &z) in atomic_numbers.iter().enumerate() {
            let solver = HartreeFockSolver::from_atoms(&[z], &[[0.0; 3]], 0).ok_or(
                QuantumChemistryMdError::UnsupportedNuclearCharge {
                    index,
                    charge: z as f64,
                },
            )?;
            neutral += solver.n_electrons as i64;
        }
        let electrons = neutral - i64::from(self.total_charge);
        if electrons < 0 {
            return Err(QuantumChemistryMdError::InvalidElectronCount {
                electrons,
                capacity: 0,
            });
        }
        if electrons % 2 != 0 {
            return Err(QuantumChemistryMdError::OpenShell(electrons));
        }
        usize::try_from(electrons).map_err(|_| QuantumChemistryMdError::InvalidElectronCount {
            electrons,
            capacity: 0,
        })
    }

    /// Converged RHF electronic energy at `positions` (Bohr).
    fn electronic_energy_at(
        &self,
        atomic_numbers: &[u32],
        n_electrons: usize,
        positions: &[[f64; 3]],
    ) -> Result<f64, QuantumChemistryMdError> {
        if n_electrons == 0 {
            return Ok(0.0);
        }
        let positions_ang: Vec<[f64; 3]> = positions
            .iter()
            .map(|p| {
                [
                    p[0] / BOHR_PER_ANGSTROM,
                    p[1] / BOHR_PER_ANGSTROM,
                    p[2] / BOHR_PER_ANGSTROM,
                ]
            })
            .collect();
        let mut solver = HartreeFockSolver::from_atoms(atomic_numbers, &positions_ang, 0).ok_or(
            QuantumChemistryMdError::DimensionMismatch {
                what: "positions",
                expected: atomic_numbers.len(),
                found: positions.len(),
            },
        )?;
        let capacity = 2 * solver.n_basis;
        if n_electrons > capacity {
            return Err(QuantumChemistryMdError::InvalidElectronCount {
                electrons: n_electrons as i64,
                capacity,
            });
        }
        solver.n_electrons = n_electrons;
        solver.energy_threshold = self.scf_energy_threshold;
        solver.max_iter = self.scf_max_iter.max(2);
        let result = solver.run_scf();
        if result.status != ScfStatus::Converged || !result.energy.is_finite() {
            return Err(QuantumChemistryMdError::ScfNotConverged);
        }
        Ok(result.energy - result.nuclear_repulsion)
    }

    /// Analytic nuclear-repulsion forces plus finite-difference electronic forces.
    fn forces_at(
        &self,
        atomic_numbers: &[u32],
        n_electrons: usize,
        positions: &[[f64; 3]],
    ) -> Result<Vec<[f64; 3]>, QuantumChemistryMdError> {
        let n = positions.len();
        let mut forces = vec![[0.0f64; 3]; n];
        // Nuclear repulsion: F_I = Σ_J Z_I Z_J (R_I − R_J)/R_IJ³.
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&positions[i], &positions[j]);
                if r < 1e-12 {
                    continue;
                }
                let coeff = self.nuclear_charges[i] * self.nuclear_charges[j] / (r * r * r);
                for k in 0..3 {
                    let fk = coeff * (positions[i][k] - positions[j][k]);
                    forces[i][k] += fk;
                    forces[j][k] -= fk;
                }
            }
        }
        if n_electrons == 0 {
            return Ok(forces);
        }
        let h = self.fd_step;
        let mut displaced = positions.to_vec();
        for i in 0..n {
            for k in 0..3 {
                let original = displaced[i][k];
                displaced[i][k] = original + h;
                let e_plus = self.electronic_energy_at(atomic_numbers, n_electrons, &displaced)?;
                displaced[i][k] = original - h;
                let e_minus = self.electronic_energy_at(atomic_numbers, n_electrons, &displaced)?;
                displaced[i][k] = original;
                forces[i][k] -= (e_plus - e_minus) / (2.0 * h);
            }
        }
        Ok(forces)
    }

    /// Forces cached by the previous step, if they belong to the current state.
    fn cached_forces(&self) -> Option<Vec<[f64; 3]>> {
        self.cache.as_ref().and_then(|c| {
            (c.positions == self.positions
                && c.nuclear_charges == self.nuclear_charges
                && c.total_charge == self.total_charge)
                .then(|| c.forces.clone())
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const PROTON_MASS: f64 = 1_836.152_673;

    fn h2(r: f64) -> HartreeFockMd {
        HartreeFockMd::new(
            vec![1.0, 1.0],
            vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0]],
            vec![PROTON_MASS, PROTON_MASS],
        )
        .expect("valid H2")
    }

    #[test]
    fn test_hf_nuclear_repulsion_single_atom() {
        let hf = HartreeFockMd::new(vec![8.0], vec![[0.0, 0.0, 0.0]], vec![16.0]).expect("O atom");
        assert_eq!(hf.nuclear_repulsion(), 0.0);
    }

    #[test]
    fn test_hf_nuclear_repulsion_two_atoms() {
        let e = h2(1.0).nuclear_repulsion();
        assert!((e - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hf_rejects_unsupported_or_fractional_charges() {
        let err = HartreeFockMd::new(vec![1.5], vec![[0.0; 3]], vec![1.0]);
        assert!(matches!(
            err,
            Err(QuantumChemistryMdError::UnsupportedNuclearCharge { index: 0, .. })
        ));
        let err = HartreeFockMd::new(
            vec![1.0, 26.0],
            vec![[0.0; 3], [2.0, 0.0, 0.0]],
            vec![1.0, 1.0],
        );
        assert!(matches!(
            err,
            Err(QuantumChemistryMdError::UnsupportedNuclearCharge { index: 1, .. })
        ));
        let err = HartreeFockMd::new(vec![1.0], vec![[0.0; 3]], vec![1.0, 2.0]);
        assert!(matches!(
            err,
            Err(QuantumChemistryMdError::DimensionMismatch { .. })
        ));
    }

    #[test]
    fn test_hf_open_shell_rejected() {
        // H2+ has one electron.
        let hf = h2(2.0).with_total_charge(1);
        assert_eq!(hf.scf_energy(), Err(QuantumChemistryMdError::OpenShell(1)));
    }

    #[test]
    fn test_hf_scf_energy_depends_on_geometry() {
        // RHF/STO-3G (ζ = 1.24) for H2 at the textbook geometry R = 1.4 bohr:
        // E_total = −1.1167 Hartree (Szabo & Ostlund, Table 3.6).
        let e_eq = h2(1.4).total_energy().expect("scf");
        assert!((e_eq - (-1.1167)).abs() < 2e-3, "E(1.4) = {e_eq}");
        let e_el_eq = h2(1.4).scf_energy().expect("scf");
        let e_el_far = h2(3.0).scf_energy().expect("scf");
        assert!(e_el_eq < 0.0);
        // Electronic energy alone is more negative at short distance.
        assert!(e_el_eq < e_el_far);
        // The total energy is bound relative to stretched H2.
        assert!(e_eq < h2(3.0).total_energy().expect("scf"));
    }

    #[test]
    fn test_hf_total_energy_is_scf_plus_repulsion() {
        let hf = h2(1.4);
        let e_total = hf.total_energy().expect("scf");
        let e_scf = hf.scf_energy().expect("scf");
        assert!((e_total - (e_scf + hf.nuclear_repulsion())).abs() < 1e-12);
    }

    #[test]
    fn test_hf_forces_are_negative_energy_gradient() {
        let hf = h2(1.8);
        let f = hf.forces().expect("forces");
        let h = 1e-3;
        let e_plus = h2(1.8 + h).total_energy().expect("scf");
        let e_minus = h2(1.8 - h).total_energy().expect("scf");
        let de_dr = (e_plus - e_minus) / (2.0 * h);
        // Atom 1 sits at +x, so F_1x = −dE/dR.
        assert!((f[1][0] + de_dr).abs() < 1e-5, "{} vs {}", f[1][0], -de_dr);
        // Newton's third law and no transverse components.
        assert!((f[0][0] + f[1][0]).abs() < 1e-6);
        assert!(f[0][1].abs() < 1e-6 && f[1][2].abs() < 1e-6);
        let g = hf.gradient().expect("gradient");
        assert!((g[1][0] + f[1][0]).abs() < 1e-15);
    }

    #[test]
    fn test_hf_bond_restoring_forces() {
        // Stretched H2 is pulled together, compressed H2 is pushed apart.
        let stretched = h2(2.2).forces().expect("forces");
        assert!(stretched[1][0] < 0.0 && stretched[0][0] > 0.0);
        let compressed = h2(1.0).forces().expect("forces");
        assert!(compressed[1][0] > 0.0 && compressed[0][0] < 0.0);
    }

    #[test]
    fn test_hf_bare_nuclei_repel_and_conserve_energy() {
        // H2^{2+}: two bare protons released at rest must fly apart.
        let mut hf = h2(1.4).with_total_charge(2);
        assert_eq!(hf.scf_energy().expect("no electrons"), 0.0);
        let e0 = hf.kinetic_energy() + hf.nuclear_repulsion();
        let r0 = dist(&hf.positions[0], &hf.positions[1]);
        for _ in 0..400 {
            hf.step(0.05).expect("step");
        }
        let r1 = dist(&hf.positions[0], &hf.positions[1]);
        assert!(r1 > r0 + 1e-3, "r0 {r0} r1 {r1}");
        assert!(hf.velocities[0][0] < 0.0 && hf.velocities[1][0] > 0.0);
        let e1 = hf.kinetic_energy() + hf.nuclear_repulsion();
        assert!((e1 - e0).abs() < 1e-6, "drift {}", e1 - e0);
    }

    #[test]
    fn test_hf_md_conserves_energy_with_electrons() {
        let mut hf = h2(1.6);
        let e0 = hf.kinetic_energy() + hf.total_energy().expect("scf");
        for _ in 0..5 {
            hf.step(5.0).expect("step");
        }
        let e1 = hf.kinetic_energy() + hf.total_energy().expect("scf");
        // The bond contracts towards the RHF minimum and energy is conserved.
        assert!(dist(&hf.positions[0], &hf.positions[1]) < 1.6);
        assert!((e1 - e0).abs() < 1e-5, "drift {}", e1 - e0);
    }

    #[test]
    fn test_hf_step_rejects_bad_dt_and_keeps_state() {
        let mut hf = h2(1.4);
        let before = hf.positions.clone();
        assert!(hf.step(-1.0).is_err());
        assert!(hf.step(f64::NAN).is_err());
        assert_eq!(hf.positions, before);
    }

    #[test]
    fn test_hf_force_cache_invalidated_by_position_change() {
        let mut hf = h2(1.4).with_total_charge(2);
        hf.step(0.1).expect("step");
        assert!(hf.cached_forces().is_some());
        hf.positions[1][0] += 0.5;
        assert!(hf.cached_forces().is_none());
    }

    #[test]
    fn test_hf_single_atom_forces_zero() {
        // Neutral He is not in the STO-3G table; use a closed-shell H^- anion.
        let hf = HartreeFockMd::new(vec![1.0], vec![[0.0; 3]], vec![PROTON_MASS])
            .expect("H")
            .with_total_charge(-1);
        let f = hf.forces().expect("forces");
        for c in f[0] {
            assert!(c.abs() < 1e-8);
        }
    }

    #[test]
    fn test_hf_forces_vanish_at_large_separation_without_electrons() {
        let hf = h2(1000.0).with_total_charge(2);
        let f = hf.forces().expect("forces");
        assert!(f[0][0].abs() < 1e-5);
    }
}
