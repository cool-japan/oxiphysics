// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Real-space grid Kohn-Sham DFT Born-Oppenheimer molecular dynamics.
//!
//! The Kohn-Sham equations of a closed-shell, spin-unpolarised molecule are
//! solved on a [`UniformGrid`] with homogeneous Dirichlet boundaries:
//!
//! ```text
//! [−½∇² + V_ext + V_H[n] + v_xc[n] − a₀ K̂] φ_i = ε_i φ_i,   n = 2 Σ_{i∈occ} |φ_i|²
//! ```
//!
//! * −½∇²: fourth-order finite differences.
//! * V_ext: every nucleus is a Gaussian charge of width σ
//!   ([`DensityFunctionalMd::nuclear_smearing`]), V_I(r) = −Z_I erf(|r−R_I|/σ)/|r−R_I|,
//!   which removes the Coulomb singularity from the grid and becomes the bare
//!   −Z_I/|r−R_I| for |r−R_I| ≫ σ.  This is an all-electron treatment with a
//!   regularised nucleus, not a norm-conserving pseudopotential.
//! * V_H: free-space Poisson equation solved by zero-padded FFT convolution
//!   (oxifft).
//! * v_xc: the functional selected by [`XcFunctional`] (see the `xc` module);
//!   for B3LYP the exact-exchange operator K̂ is applied with the occupied
//!   orbitals of the previous iteration.
//! * Diagonalisation: LOBPCG with a kinetic-energy preconditioner.
//! * Self-consistency: Pulay (DIIS) density mixing on top of the linear
//!   step n ← (1−α) n_in + α n_out (plain linear mixing when the history is 0).
//!
//! Ionic forces are Hellmann-Feynman forces, exact derivatives of the discrete
//! total energy at self-consistency (the grid does not move with the atoms, so
//! there are no Pulay terms).  Atomic units throughout.

use super::error::QuantumChemistryMdError;
use super::ks_eigensolver::{EigenPairs, lobpcg};
use super::real_space_grid::{FreeSpacePoisson, KineticPreconditioner, UniformGrid, erf};
use super::xc;
use super::{XcFunctional, dist};

/// 1 e·Bohr in Debye.
const E_BOHR_TO_DEBYE: f64 = 2.541_746_473;
/// 2/√π.
const TWO_OVER_SQRT_PI: f64 = 1.128_379_167_095_512_6;

/// Components of the Kohn-Sham total energy \[Hartree\].
#[derive(Debug, Clone, Default, PartialEq)]
pub struct KsEnergies {
    /// Non-interacting kinetic energy T_s = 2 Σ_i ⟨φ_i|−½∇²|φ_i⟩.
    pub kinetic: f64,
    /// Electron-nucleus energy ∫ V_ext n dV.
    pub external: f64,
    /// Hartree energy ½ ∫ V_H n dV.
    pub hartree: f64,
    /// Semi-local exchange energy of the functional.
    pub exchange: f64,
    /// Exact-exchange contribution a₀ E_x^HF (non-zero only for hybrids).
    pub exact_exchange: f64,
    /// Correlation energy of the functional.
    pub correlation: f64,
    /// Ion-ion repulsion Σ_{I<J} Z_I Z_J / R_IJ.
    pub ion_ion: f64,
    /// Sum of all terms.
    pub total: f64,
}

/// Result of one Kohn-Sham iteration.
#[derive(Debug, Clone, PartialEq)]
pub struct KsStepReport {
    /// Energy of the output orbitals and density.
    pub energies: KsEnergies,
    /// ∫ |n_out − n_in| dV (electrons).
    pub density_residual: f64,
    /// Kohn-Sham eigenvalues of all computed bands \[Hartree\].
    pub eigenvalues: Vec<f64>,
    /// LOBPCG iterations used for the diagonalisation.
    pub eigensolver_iterations: usize,
    /// Largest eigensolver residual ‖Hφ − εφ‖ of the occupied bands.
    pub eigensolver_residual: f64,
}

/// Result of a converged self-consistent-field run.
#[derive(Debug, Clone, PartialEq)]
pub struct KsScfReport {
    /// Energy at self-consistency.
    pub energies: KsEnergies,
    /// Number of Kohn-Sham iterations performed.
    pub steps: usize,
    /// Final ∫ |n_out − n_in| dV.
    pub density_residual: f64,
    /// Kohn-Sham eigenvalues \[Hartree\].
    pub eigenvalues: Vec<f64>,
}

/// Grid-dependent precomputed operators.
#[derive(Debug, Clone)]
struct GridOperators {
    poisson: FreeSpacePoisson,
    preconditioner: KineticPreconditioner,
}

impl GridOperators {
    fn new(grid: &UniformGrid) -> Self {
        Self {
            poisson: FreeSpacePoisson::new(grid),
            preconditioner: KineticPreconditioner::new(grid),
        }
    }

    fn matches(&self, grid: &UniformGrid) -> bool {
        self.poisson.grid() == grid && self.preconditioner.grid() == grid
    }
}

/// Forces cached by [`DensityFunctionalMd::step`] together with the positions
/// they belong to.
#[derive(Debug, Clone)]
struct ForceCache {
    positions: Vec<[f64; 3]>,
    forces: Vec<[f64; 3]>,
}

/// Grid Kohn-Sham DFT Born-Oppenheimer MD driver (see the module documentation).
#[derive(Debug, Clone)]
pub struct DensityFunctionalMd {
    /// Exchange-correlation functional choice.
    pub functional: XcFunctional,
    /// Real-space grid.
    pub grid: UniformGrid,
    /// Electron density on the grid \[e/Bohr³\] (storage order of [`UniformGrid`]).
    pub density: Vec<f64>,
    /// Nuclear charges Z_I.
    pub nuclear_charges: Vec<f64>,
    /// Ionic positions \[Bohr\].
    pub ionic_positions: Vec<[f64; 3]>,
    /// Ionic velocities \[Bohr/a.u.\].
    pub ionic_velocities: Vec<[f64; 3]>,
    /// Ionic masses \[electron masses\] (needed only by [`step`](Self::step)).
    pub ionic_masses: Vec<f64>,
    /// Net charge of the system; the electron count is Σ Z_I − `total_charge`.
    pub total_charge: i32,
    /// Width σ of the Gaussian nuclear charges \[Bohr\].
    pub nuclear_smearing: f64,
    /// Density-mixing parameter α ∈ (0, 1].
    pub mixing: f64,
    /// Number of previous iterations used by Pulay (DIIS) density mixing;
    /// 0 or 1 selects plain linear mixing n ← (1−α) n_in + α n_out.
    pub pulay_history: usize,
    /// Unoccupied bands computed in addition to the occupied ones.
    pub extra_bands: usize,
    /// Residual tolerance of the eigensolver (unit-norm grid vectors).
    pub eigensolver_tolerance: f64,
    /// Maximum LOBPCG iterations per Kohn-Sham step.
    pub eigensolver_max_iter: usize,
    /// SCF convergence threshold on ∫|n_out − n_in| dV and on the energy change.
    pub scf_tolerance: f64,
    /// Maximum Kohn-Sham iterations of [`run_scf`](Self::run_scf).
    pub scf_max_steps: usize,
    orbitals: Vec<Vec<f64>>,
    eigenvalues: Vec<f64>,
    output_density: Vec<f64>,
    mixing_history: Vec<(Vec<f64>, Vec<f64>)>,
    operators: Option<GridOperators>,
    force_cache: Option<ForceCache>,
}

impl DensityFunctionalMd {
    /// Construct a DFT driver without nuclei.
    ///
    /// # Arguments
    /// * `functional` – XC functional to use.
    /// * `grid` – Real-space grid.
    /// * `density` – Initial electron density on `grid` (may be all zero).
    ///
    /// The nuclear smearing defaults to 1.5 grid spacings, the mixing to 0.3.
    ///
    /// # Errors
    /// [`QuantumChemistryMdError::DimensionMismatch`] if `density` does not
    /// match the grid, [`QuantumChemistryMdError::InvalidParameter`] for a
    /// negative or non-finite density value.
    pub fn new(
        functional: XcFunctional,
        grid: UniformGrid,
        density: Vec<f64>,
    ) -> Result<Self, QuantumChemistryMdError> {
        let smearing = 1.5 * grid.spacing;
        let md = Self {
            functional,
            grid,
            density,
            nuclear_charges: Vec::new(),
            ionic_positions: Vec::new(),
            ionic_velocities: Vec::new(),
            ionic_masses: Vec::new(),
            total_charge: 0,
            nuclear_smearing: smearing,
            mixing: 0.3,
            pulay_history: 6,
            extra_bands: 2,
            eigensolver_tolerance: 1e-7,
            eigensolver_max_iter: 400,
            scf_tolerance: 1e-7,
            scf_max_steps: 300,
            orbitals: Vec::new(),
            eigenvalues: Vec::new(),
            output_density: Vec::new(),
            mixing_history: Vec::new(),
            operators: None,
            force_cache: None,
        };
        md.validate()?;
        Ok(md)
    }

    /// Add nuclei (builder style).
    ///
    /// # Errors
    /// [`QuantumChemistryMdError::DimensionMismatch`] if the three slices
    /// differ in length, [`QuantumChemistryMdError::InvalidParameter`] for a
    /// non-positive charge or mass or a non-finite position.
    pub fn with_ions(
        mut self,
        nuclear_charges: Vec<f64>,
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
    ) -> Result<Self, QuantumChemistryMdError> {
        self.ionic_velocities = vec![[0.0; 3]; positions.len()];
        self.nuclear_charges = nuclear_charges;
        self.ionic_positions = positions;
        self.ionic_masses = masses;
        self.force_cache = None;
        self.validate()?;
        if self.ionic_masses.len() != self.ionic_positions.len() {
            return Err(QuantumChemistryMdError::DimensionMismatch {
                what: "ionic_masses",
                expected: self.ionic_positions.len(),
                found: self.ionic_masses.len(),
            });
        }
        Ok(self)
    }

    /// Set the net charge of the system (builder style).
    pub fn with_total_charge(mut self, total_charge: i32) -> Self {
        self.total_charge = total_charge;
        self
    }

    /// Number of electrons Σ Z_I − `total_charge`.
    ///
    /// # Errors
    /// [`QuantumChemistryMdError::InvalidElectronCount`] if it is negative or
    /// not an integer, [`QuantumChemistryMdError::OpenShell`] if it is odd.
    pub fn electron_count(&self) -> Result<usize, QuantumChemistryMdError> {
        let n: f64 = self.nuclear_charges.iter().sum::<f64>() - f64::from(self.total_charge);
        let rounded = n.round();
        if !n.is_finite() || (n - rounded).abs() > 1e-9 || rounded < 0.0 {
            return Err(QuantumChemistryMdError::InvalidElectronCount {
                electrons: rounded as i64,
                capacity: self.grid.n_points(),
            });
        }
        let electrons = rounded as i64;
        if electrons % 2 != 0 {
            return Err(QuantumChemistryMdError::OpenShell(electrons));
        }
        usize::try_from(electrons).map_err(|_| QuantumChemistryMdError::InvalidElectronCount {
            electrons,
            capacity: self.grid.n_points(),
        })
    }

    /// Occupied and extra Kohn-Sham orbitals of the last step, normalised as
    /// ∫ φ² dV = 1.
    pub fn orbitals(&self) -> &[Vec<f64>] {
        &self.orbitals
    }

    /// Kohn-Sham eigenvalues of the last step \[Hartree\].
    pub fn eigenvalues(&self) -> &[f64] {
        &self.eigenvalues
    }

    /// Total number of electrons ∫ n dV in [`density`](Self::density).
    pub fn integrated_density(&self) -> f64 {
        self.density.iter().sum::<f64>() * self.grid.dv()
    }

    /// Semi-local exchange energy of [`density`](Self::density) in Hartree.
    ///
    /// LDA: Slater; PBE: PBE exchange; B3LYP: 0.8 Slater + 0.72 ΔB88.  The
    /// exact-exchange part of a hybrid depends on the orbitals and is given by
    /// [`exact_exchange_energy`](Self::exact_exchange_energy).
    pub fn exchange_energy(&self) -> f64 {
        self.xc_on_grid(&self.density).exchange
    }

    /// Correlation energy of [`density`](Self::density) in Hartree
    /// (LDA: VWN5; PBE: PBE; B3LYP: 0.19 VWN-RPA + 0.81 LYP).
    pub fn correlation_energy(&self) -> f64 {
        self.xc_on_grid(&self.density).correlation
    }

    /// Hartree-Fock exchange energy E_x^HF = −Σ_{ij} (ij|ji) of the occupied
    /// orbitals of the last Kohn-Sham step (unscaled by the hybrid fraction).
    ///
    /// # Errors
    /// Invalid electron count.
    pub fn exact_exchange_energy(&self) -> Result<f64, QuantumChemistryMdError> {
        let n_occ = self.electron_count()? / 2;
        if self.orbitals.len() < n_occ {
            return Ok(0.0);
        }
        let operators = self.operators_for_grid();
        Ok(self.hf_exchange_energy(&operators.poisson, &self.orbitals[..n_occ]))
    }

    /// Hartree energy ½ ∫ V_H[n] n dV of [`density`](Self::density).
    pub fn hartree_energy(&self) -> f64 {
        let operators = self.operators_for_grid();
        let v_h = operators.poisson.potential(&self.density);
        0.5 * dot(&v_h, &self.density) * self.grid.dv()
    }

    /// External (electron-nucleus) potential on the grid \[Hartree\].
    pub fn external_potential(&self) -> Vec<f64> {
        let sigma = self.nuclear_smearing;
        self.grid
            .points()
            .iter()
            .map(|g| {
                self.nuclear_charges
                    .iter()
                    .zip(self.ionic_positions.iter())
                    .map(|(&z, r)| -z * smeared_coulomb(dist(g, r), sigma))
                    .sum()
            })
            .collect()
    }

    /// Ion-ion repulsion energy Σ_{I<J} Z_I Z_J / R_IJ \[Hartree\].
    pub fn ion_ion_energy(&self) -> f64 {
        let n = self.ionic_positions.len().min(self.nuclear_charges.len());
        let mut e = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.ionic_positions[i], &self.ionic_positions[j]);
                if r > 0.0 {
                    e += self.nuclear_charges[i] * self.nuclear_charges[j] / r;
                }
            }
        }
        e
    }

    /// Perform one Kohn-Sham self-consistency iteration.
    ///
    /// Builds V_eff = V_ext + V_H[n_in] + v_xc[n_in] (plus −a₀K̂ for hybrids)
    /// from the current density, solves for the lowest
    /// N/2 + [`extra_bands`](Self::extra_bands) orbitals, forms
    /// n_out = 2 Σ_occ |φ_i|² and mixes n ← (1−α) n_in + α n_out, accelerated by
    /// Pulay mixing over the last [`pulay_history`](Self::pulay_history) steps.
    ///
    /// # Errors
    /// Invalid configuration or electron count, too many bands for the grid,
    /// or eigensolver non-convergence (the state is then unchanged).
    pub fn kohn_sham_step(&mut self) -> Result<KsStepReport, QuantumChemistryMdError> {
        self.validate()?;
        let n_electrons = self.electron_count()?;
        let operators = self.take_operators();
        let result = self.kohn_sham_step_with(&operators, n_electrons);
        self.operators = Some(operators);
        result
    }

    /// Iterate [`kohn_sham_step`](Self::kohn_sham_step) to self-consistency.
    ///
    /// Converged when ∫|n_out − n_in| dV and the change of the total energy
    /// are both below [`scf_tolerance`](Self::scf_tolerance).  On convergence
    /// [`density`](Self::density) is set to the output density of the final
    /// orbitals, so that [`forces`](Self::forces) are the Hellmann-Feynman
    /// forces of the converged state.
    ///
    /// # Errors
    /// [`QuantumChemistryMdError::ScfNotConverged`] after
    /// [`scf_max_steps`](Self::scf_max_steps) iterations, otherwise as for
    /// [`kohn_sham_step`](Self::kohn_sham_step).
    pub fn run_scf(&mut self) -> Result<KsScfReport, QuantumChemistryMdError> {
        self.mixing_history.clear();
        let mut previous_energy = f64::INFINITY;
        for step in 1..=self.scf_max_steps {
            let report = self.kohn_sham_step()?;
            let de = (report.energies.total - previous_energy).abs();
            previous_energy = report.energies.total;
            if report.density_residual < self.scf_tolerance && de < self.scf_tolerance {
                self.density = self.output_density.clone();
                return Ok(KsScfReport {
                    energies: report.energies,
                    steps: step,
                    density_residual: report.density_residual,
                    eigenvalues: report.eigenvalues,
                });
            }
        }
        Err(QuantumChemistryMdError::ScfNotConverged)
    }

    /// Hellmann-Feynman forces on the ions in Hartree/Bohr.
    ///
    /// F_I = −∫ n(r) ∂V_I(r)/∂R_I dV + Σ_{J≠I} Z_I Z_J (R_I − R_J)/R_IJ³
    ///     = Z_I Σ_g n_g dV φ'(u)/u (R_I − r_g) + ion-ion,
    ///
    /// with u = |r_g − R_I| and φ(u) = erf(u/σ)/u (φ' < 0, so the electron
    /// density attracts the nucleus).  For σ → 0 the electronic term is the
    /// bare −Z_I Σ_g n_g dV (R_I − r_g)/u³.  The forces use
    /// [`density`](Self::density) as it is; call [`run_scf`](Self::run_scf)
    /// first for Born-Oppenheimer forces.
    ///
    /// # Errors
    /// Inconsistent ionic data or density size.
    pub fn forces(&self) -> Result<Vec<[f64; 3]>, QuantumChemistryMdError> {
        self.validate()?;
        let dv = self.grid.dv();
        let sigma = self.nuclear_smearing;
        let points = self.grid.points();
        let mut forces = vec![[0.0f64; 3]; self.ionic_positions.len()];
        for ((f, r_i), &z_i) in forces
            .iter_mut()
            .zip(self.ionic_positions.iter())
            .zip(self.nuclear_charges.iter())
        {
            for (g, &n_g) in points.iter().zip(self.density.iter()) {
                if n_g == 0.0 {
                    continue;
                }
                let u = dist(g, r_i);
                let coeff = z_i * n_g * dv * smeared_coulomb_derivative_over_u(u, sigma);
                for k in 0..3 {
                    f[k] += coeff * (r_i[k] - g[k]);
                }
            }
        }
        let n = self.ionic_positions.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist(&self.ionic_positions[i], &self.ionic_positions[j]);
                if r < 1e-12 {
                    continue;
                }
                let coeff = self.nuclear_charges[i] * self.nuclear_charges[j] / (r * r * r);
                for k in 0..3 {
                    let fk = coeff * (self.ionic_positions[i][k] - self.ionic_positions[j][k]);
                    forces[i][k] += fk;
                    forces[j][k] -= fk;
                }
            }
        }
        Ok(forces)
    }

    /// Born-Oppenheimer velocity-Verlet step of length `dt` (atomic units).
    ///
    /// The electronic structure is converged with [`run_scf`](Self::run_scf)
    /// at the new positions; the forces at the end of a step are cached and
    /// reused at the start of the next one.  Returns the SCF report at the new
    /// positions.  On error the ionic state, density and orbitals are restored.
    ///
    /// # Errors
    /// [`QuantumChemistryMdError::InvalidParameter`] for a bad `dt` or mass,
    /// [`QuantumChemistryMdError::DimensionMismatch`] for missing masses,
    /// otherwise as for [`run_scf`](Self::run_scf).
    pub fn step(&mut self, dt: f64) -> Result<KsScfReport, QuantumChemistryMdError> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "dt",
                value: dt,
            });
        }
        if self.ionic_masses.len() != self.ionic_positions.len() {
            return Err(QuantumChemistryMdError::DimensionMismatch {
                what: "ionic_masses",
                expected: self.ionic_positions.len(),
                found: self.ionic_masses.len(),
            });
        }
        if let Some(&m) = self
            .ionic_masses
            .iter()
            .find(|m| !(m.is_finite() && **m > 0.0))
        {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "ionic mass",
                value: m,
            });
        }
        let saved = (
            self.ionic_positions.clone(),
            self.ionic_velocities.clone(),
            self.density.clone(),
            self.orbitals.clone(),
            self.eigenvalues.clone(),
        );
        let result = self.verlet_step(dt);
        if result.is_err() {
            let (positions, velocities, density, orbitals, eigenvalues) = saved;
            self.ionic_positions = positions;
            self.ionic_velocities = velocities;
            self.density = density;
            self.orbitals = orbitals;
            self.eigenvalues = eigenvalues;
        }
        result
    }

    /// Ionic kinetic energy Σ ½ m v² \[Hartree\].
    pub fn ionic_kinetic_energy(&self) -> f64 {
        self.ionic_velocities
            .iter()
            .zip(self.ionic_masses.iter())
            .map(|(v, &m)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }

    /// Electronic dipole moment −∫ r n(r) dV in Debye (electron charge −1).
    pub fn dipole_moment(&self) -> [f64; 3] {
        let dv = self.grid.dv();
        let mut d = [0.0f64; 3];
        for (g, &n_g) in self.grid.points().iter().zip(self.density.iter()) {
            for k in 0..3 {
                d[k] -= n_g * g[k] * dv;
            }
        }
        d.map(|x| x * E_BOHR_TO_DEBYE)
    }

    /// Total dipole moment Σ_I Z_I R_I − ∫ r n(r) dV in Debye.
    pub fn total_dipole_moment(&self) -> [f64; 3] {
        let mut d = self.dipole_moment();
        for (&z, r) in self.nuclear_charges.iter().zip(self.ionic_positions.iter()) {
            for k in 0..3 {
                d[k] += z * r[k] * E_BOHR_TO_DEBYE;
            }
        }
        d
    }

    // -- private helpers --

    fn validate(&self) -> Result<(), QuantumChemistryMdError> {
        let n_points = self.grid.n_points();
        if self.density.len() != n_points {
            return Err(QuantumChemistryMdError::DimensionMismatch {
                what: "density",
                expected: n_points,
                found: self.density.len(),
            });
        }
        if let Some(&x) = self.density.iter().find(|x| !(x.is_finite() && **x >= 0.0)) {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "density",
                value: x,
            });
        }
        let n_ions = self.ionic_positions.len();
        for (what, found) in [
            ("nuclear_charges", self.nuclear_charges.len()),
            ("ionic_velocities", self.ionic_velocities.len()),
        ] {
            if found != n_ions {
                return Err(QuantumChemistryMdError::DimensionMismatch {
                    what,
                    expected: n_ions,
                    found,
                });
            }
        }
        if let Some(&z) = self
            .nuclear_charges
            .iter()
            .find(|z| !(z.is_finite() && **z > 0.0))
        {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "nuclear charge",
                value: z,
            });
        }
        if let Some(&x) = self
            .ionic_positions
            .iter()
            .flat_map(|p| p.iter())
            .find(|x| !x.is_finite())
        {
            return Err(QuantumChemistryMdError::InvalidParameter {
                name: "ionic position",
                value: x,
            });
        }
        for (name, value, ok) in [
            (
                "nuclear_smearing",
                self.nuclear_smearing,
                self.nuclear_smearing.is_finite() && self.nuclear_smearing > 0.0,
            ),
            (
                "mixing",
                self.mixing,
                self.mixing.is_finite() && self.mixing > 0.0 && self.mixing <= 1.0,
            ),
            (
                "eigensolver_tolerance",
                self.eigensolver_tolerance,
                self.eigensolver_tolerance.is_finite() && self.eigensolver_tolerance > 0.0,
            ),
        ] {
            if !ok {
                return Err(QuantumChemistryMdError::InvalidParameter { name, value });
            }
        }
        Ok(())
    }

    /// Grid operators, rebuilt if the grid changed.
    fn take_operators(&mut self) -> GridOperators {
        match self.operators.take() {
            Some(ops) if ops.matches(&self.grid) => ops,
            _ => GridOperators::new(&self.grid),
        }
    }

    /// Cached operators when valid, otherwise freshly built ones.
    fn operators_for_grid(&self) -> std::borrow::Cow<'_, GridOperators> {
        match &self.operators {
            Some(ops) if ops.matches(&self.grid) => std::borrow::Cow::Borrowed(ops),
            _ => std::borrow::Cow::Owned(GridOperators::new(&self.grid)),
        }
    }

    fn kohn_sham_step_with(
        &mut self,
        operators: &GridOperators,
        n_electrons: usize,
    ) -> Result<KsStepReport, QuantumChemistryMdError> {
        let n_points = self.grid.n_points();
        let dv = self.grid.dv();
        let n_occ = n_electrons / 2;
        let n_bands = if n_occ == 0 {
            0
        } else {
            n_occ + self.extra_bands
        };
        if n_bands > n_points {
            return Err(QuantumChemistryMdError::InvalidElectronCount {
                electrons: n_electrons as i64,
                capacity: n_points,
            });
        }
        let n_in = self.density.clone();
        let v_ext = self.external_potential();
        let a0 = xc::exact_exchange_fraction(&self.functional);

        let (orbitals, eigenvalues, iterations, residual) = if n_occ == 0 {
            (Vec::new(), Vec::new(), 0, 0.0)
        } else {
            let v_h = operators.poisson.potential(&n_in);
            let xc_in = self.xc_on_grid(&n_in);
            let v_loc: Vec<f64> = v_ext
                .iter()
                .zip(v_h.iter())
                .zip(xc_in.potential.iter())
                .map(|((a, b), c)| a + b + c)
                .collect();
            let initial = self.initial_block(n_bands);
            let reference: Vec<Vec<f64>> = if a0 > 0.0 {
                if self.orbitals.len() >= n_occ && self.orbitals.iter().all(|o| o.len() == n_points)
                {
                    self.orbitals[..n_occ].to_vec()
                } else {
                    // Bootstrap the exchange operator from the semi-local problem.
                    let pairs =
                        self.diagonalise(operators, &v_loc, 0.0, &[], initial.clone(), n_occ)?;
                    physical_orbitals(&pairs.vectors, dv)[..n_occ].to_vec()
                }
            } else {
                Vec::new()
            };
            let pairs = self.diagonalise(operators, &v_loc, a0, &reference, initial, n_occ)?;
            (
                physical_orbitals(&pairs.vectors, dv),
                pairs.values,
                pairs.iterations,
                pairs.max_residual,
            )
        };

        let mut n_out = vec![0.0f64; n_points];
        for phi in orbitals.iter().take(n_occ) {
            for (n, p) in n_out.iter_mut().zip(phi.iter()) {
                *n += 2.0 * p * p;
            }
        }
        let energies = self.energy_terms(
            operators,
            &orbitals[..n_occ.min(orbitals.len())],
            &n_out,
            &v_ext,
        );
        let density_residual: f64 = n_out
            .iter()
            .zip(n_in.iter())
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            * dv;
        self.density = self.mix_density(n_in, &n_out);
        self.output_density = n_out;
        self.orbitals = orbitals;
        self.eigenvalues = eigenvalues.clone();
        Ok(KsStepReport {
            energies,
            density_residual,
            eigenvalues,
            eigensolver_iterations: iterations,
            eigensolver_residual: residual,
        })
    }

    /// Solve for the lowest bands of −½∇² + v_loc − a₀ K̂[reference].
    ///
    /// Returns the Euclidean-normalised Ritz pairs of the whole block.
    fn diagonalise(
        &self,
        operators: &GridOperators,
        v_loc: &[f64],
        a0: f64,
        reference: &[Vec<f64>],
        initial: Vec<Vec<f64>>,
        n_wanted: usize,
    ) -> Result<EigenPairs, QuantumChemistryMdError> {
        let grid = &self.grid;
        let apply = |block: &[Vec<f64>]| -> Vec<Vec<f64>> {
            block
                .iter()
                .map(|v| {
                    let lap = grid.laplacian(v);
                    let mut hv: Vec<f64> = lap
                        .iter()
                        .zip(v.iter())
                        .zip(v_loc.iter())
                        .map(|((l, x), u)| -0.5 * l + u * x)
                        .collect();
                    if a0 > 0.0 {
                        for phi in reference {
                            let pair: Vec<f64> =
                                phi.iter().zip(v.iter()).map(|(a, b)| a * b).collect();
                            let v_pair = operators.poisson.potential(&pair);
                            for ((h, p), w) in hv.iter_mut().zip(phi.iter()).zip(v_pair.iter()) {
                                *h -= a0 * p * w;
                            }
                        }
                    }
                    hv
                })
                .collect()
        };
        let precondition =
            |r: &[f64], lambda: f64| operators.preconditioner.apply(r, (-lambda).max(0.0) + 0.5);
        lobpcg(
            apply,
            precondition,
            initial,
            n_wanted,
            self.eigensolver_tolerance,
            self.eigensolver_max_iter,
        )
    }

    /// Starting block: the previous orbitals when compatible, otherwise
    /// Gaussian-type functions on the nuclei with a deterministic perturbation
    /// that breaks every spatial symmetry.
    fn initial_block(&self, n_bands: usize) -> Vec<Vec<f64>> {
        let n_points = self.grid.n_points();
        let sqrt_dv = self.grid.dv().sqrt();
        if self.orbitals.len() == n_bands && self.orbitals.iter().all(|o| o.len() == n_points) {
            return self
                .orbitals
                .iter()
                .map(|o| o.iter().map(|x| x * sqrt_dv).collect())
                .collect();
        }
        let centres: Vec<[f64; 3]> = if self.ionic_positions.is_empty() {
            let p0 = self.grid.point(0);
            let p1 = self.grid.point(n_points - 1);
            vec![[
                0.5 * (p0[0] + p1[0]),
                0.5 * (p0[1] + p1[1]),
                0.5 * (p0[2] + p1[2]),
            ]]
        } else {
            self.ionic_positions.clone()
        };
        let points = self.grid.points();
        let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
        (0..n_bands)
            .map(|band| {
                points
                    .iter()
                    .map(|g| {
                        let mut value = 0.0;
                        for (atom, c) in centres.iter().enumerate() {
                            let d = [g[0] - c[0], g[1] - c[1], g[2] - c[2]];
                            let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                            let angular = match band % 4 {
                                0 => 1.0,
                                1 => d[0],
                                2 => d[1],
                                _ => d[2],
                            };
                            // Alternate the relative phase of the centres for
                            // every second shell of trial functions.
                            let sign = if (band / 4) % 2 == 1 && atom % 2 == 1 {
                                -1.0
                            } else {
                                1.0
                            };
                            value += sign * angular * (-0.5 * r2).exp();
                        }
                        state = state
                            .wrapping_mul(6_364_136_223_846_793_005)
                            .wrapping_add(1_442_695_040_888_963_407);
                        let noise = ((state >> 11) as f64 / (1u64 << 53) as f64) - 0.5;
                        value + 1e-3 * noise
                    })
                    .collect()
            })
            .collect()
    }

    /// Kohn-Sham total energy of `occupied` orbitals with density `n`.
    fn energy_terms(
        &self,
        operators: &GridOperators,
        occupied: &[Vec<f64>],
        n: &[f64],
        v_ext: &[f64],
    ) -> KsEnergies {
        let dv = self.grid.dv();
        let kinetic: f64 = occupied
            .iter()
            .map(|phi| {
                let lap = self.grid.laplacian(phi);
                2.0 * (-0.5) * dot(phi, &lap) * dv
            })
            .sum();
        let external = dot(v_ext, n) * dv;
        let v_h = operators.poisson.potential(n);
        let hartree = 0.5 * dot(&v_h, n) * dv;
        let xc_out = self.xc_on_grid(n);
        let a0 = xc::exact_exchange_fraction(&self.functional);
        let exact_exchange = if a0 > 0.0 {
            a0 * self.hf_exchange_energy(&operators.poisson, occupied)
        } else {
            0.0
        };
        let ion_ion = self.ion_ion_energy();
        let total = kinetic
            + external
            + hartree
            + xc_out.exchange
            + exact_exchange
            + xc_out.correlation
            + ion_ion;
        KsEnergies {
            kinetic,
            external,
            hartree,
            exchange: xc_out.exchange,
            exact_exchange,
            correlation: xc_out.correlation,
            ion_ion,
            total,
        }
    }

    /// E_x^HF = −Σ_{i,j} ∫ φ_iφ_j V[φ_iφ_j] dV for closed-shell orbitals.
    fn hf_exchange_energy(&self, poisson: &FreeSpacePoisson, occupied: &[Vec<f64>]) -> f64 {
        let dv = self.grid.dv();
        let mut e = 0.0;
        for (i, phi_i) in occupied.iter().enumerate() {
            for (j, phi_j) in occupied.iter().enumerate().skip(i) {
                let pair: Vec<f64> = phi_i.iter().zip(phi_j.iter()).map(|(a, b)| a * b).collect();
                let v_pair = poisson.potential(&pair);
                let integral = dot(&pair, &v_pair) * dv;
                e -= if i == j { integral } else { 2.0 * integral };
            }
        }
        e
    }

    /// Semi-local XC energies and potential of density `n` on the grid.
    pub(super) fn xc_on_grid(&self, n: &[f64]) -> XcOnGrid {
        let dv = self.grid.dv();
        let gga = xc::is_gga(&self.functional);
        let gradient: Vec<Vec<f64>> = if gga {
            (0..3).map(|axis| self.grid.derivative(n, axis)).collect()
        } else {
            Vec::new()
        };
        let mut exchange = 0.0;
        let mut correlation = 0.0;
        let mut potential = vec![0.0f64; n.len()];
        let mut d_sigma = vec![0.0f64; if gga { n.len() } else { 0 }];
        for (idx, &rho) in n.iter().enumerate() {
            let sigma = if gga {
                gradient.iter().map(|g| g[idx] * g[idx]).sum()
            } else {
                0.0
            };
            let e = xc::energy_density(&self.functional, rho, sigma);
            exchange += e.exchange * dv;
            correlation += e.correlation * dv;
            let d = xc::derivatives(&self.functional, rho, sigma);
            potential[idx] = d.d_rho;
            if gga {
                d_sigma[idx] = d.d_sigma;
            }
        }
        if gga {
            // v = ∂f/∂ρ − 2 Σ_k D_k(∂f/∂σ · D_k ρ): the exact derivative of the
            // discrete energy because the central difference D_k is antisymmetric.
            for (axis, g) in gradient.iter().enumerate() {
                let flux: Vec<f64> = d_sigma.iter().zip(g.iter()).map(|(a, b)| a * b).collect();
                let div = self.grid.derivative(&flux, axis);
                for (v, dvv) in potential.iter_mut().zip(div.iter()) {
                    *v -= 2.0 * dvv;
                }
            }
        }
        XcOnGrid {
            exchange,
            correlation,
            potential,
        }
    }

    /// Next input density from the current input and output densities.
    ///
    /// Pulay/DIIS: with residuals F_i = n_out,i − n_in,i of the stored steps,
    /// minimise ‖Σ c_i F_i‖ subject to Σ c_i = 1 and return
    /// Σ c_i (n_in,i + α F_i); falls back to linear mixing when the history is
    /// too short or the DIIS system is singular.  Negative values are clipped.
    fn mix_density(&mut self, n_in: Vec<f64>, n_out: &[f64]) -> Vec<f64> {
        let alpha = self.mixing;
        let residual: Vec<f64> = n_out.iter().zip(n_in.iter()).map(|(o, i)| o - i).collect();
        let linear: Vec<f64> = n_in
            .iter()
            .zip(residual.iter())
            .map(|(i, f)| (i + alpha * f).max(0.0))
            .collect();
        if self.pulay_history < 2 {
            self.mixing_history.clear();
            return linear;
        }
        if self
            .mixing_history
            .first()
            .is_some_and(|(n, _)| n.len() != n_in.len())
        {
            self.mixing_history.clear();
        }
        self.mixing_history.push((n_in, residual));
        while self.mixing_history.len() > self.pulay_history {
            self.mixing_history.remove(0);
        }
        let k = self.mixing_history.len();
        if k < 2 {
            return linear;
        }
        let mut system = vec![vec![0.0f64; k + 1]; k + 1];
        for i in 0..k {
            for j in i..k {
                let b = dot(&self.mixing_history[i].1, &self.mixing_history[j].1);
                system[i][j] = b;
                system[j][i] = b;
            }
            system[i][k] = 1.0;
            system[k][i] = 1.0;
        }
        let mut rhs = vec![0.0f64; k + 1];
        rhs[k] = 1.0;
        match solve_dense(system, rhs) {
            Some(c) => {
                let mut next = vec![0.0f64; n_out.len()];
                for (ci, (n_i, f_i)) in c.iter().zip(self.mixing_history.iter()) {
                    for ((x, a), f) in next.iter_mut().zip(n_i.iter()).zip(f_i.iter()) {
                        *x += ci * (a + alpha * f);
                    }
                }
                next.iter().map(|x| x.max(0.0)).collect()
            }
            None => {
                self.mixing_history.clear();
                linear
            }
        }
    }

    fn verlet_step(&mut self, dt: f64) -> Result<KsScfReport, QuantumChemistryMdError> {
        let f_start = match &self.force_cache {
            Some(cache) if cache.positions == self.ionic_positions => cache.forces.clone(),
            _ => {
                self.run_scf()?;
                self.forces()?
            }
        };
        let masses = self.ionic_masses.clone();
        for ((pos, vel), (f, &m)) in self
            .ionic_positions
            .iter_mut()
            .zip(self.ionic_velocities.iter_mut())
            .zip(f_start.iter().zip(masses.iter()))
        {
            for k in 0..3 {
                vel[k] += 0.5 * dt * f[k] / m;
                pos[k] += dt * vel[k];
            }
        }
        self.force_cache = None;
        let report = self.run_scf()?;
        let f_end = self.forces()?;
        for ((vel, f), &m) in self
            .ionic_velocities
            .iter_mut()
            .zip(f_end.iter())
            .zip(masses.iter())
        {
            for k in 0..3 {
                vel[k] += 0.5 * dt * f[k] / m;
            }
        }
        self.force_cache = Some(ForceCache {
            positions: self.ionic_positions.clone(),
            forces: f_end,
        });
        Ok(report)
    }
}

/// Semi-local XC quantities on the grid.
pub(super) struct XcOnGrid {
    /// Semi-local exchange energy \[Hartree\].
    pub(super) exchange: f64,
    /// Correlation energy \[Hartree\].
    pub(super) correlation: f64,
    /// v_xc on the grid \[Hartree\].
    pub(super) potential: Vec<f64>,
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

/// Solve the dense system `a x = b` by Gaussian elimination with partial
/// pivoting; `None` if it is numerically singular.
fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    let scale = a
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0f64, |m, x| m.max(x.abs()));
    if !(scale.is_finite() && scale > 0.0) {
        return None;
    }
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| {
            a[i][col]
                .abs()
                .partial_cmp(&a[j][col].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        if a[pivot][col].abs() <= 1e-14 * scale {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in (col + 1)..n {
            let factor = a[row][col] / a[col][col];
            if factor != 0.0 {
                for c in col..n {
                    a[row][c] -= factor * a[col][c];
                }
                b[row] -= factor * b[col];
            }
        }
    }
    let mut x = vec![0.0f64; n];
    for row in (0..n).rev() {
        let mut acc = b[row];
        for c in (row + 1)..n {
            acc -= a[row][c] * x[c];
        }
        x[row] = acc / a[row][row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// Convert Euclidean-normalised grid vectors to ∫φ² dV = 1 normalisation.
fn physical_orbitals(vectors: &[Vec<f64>], dv: f64) -> Vec<Vec<f64>> {
    let scale = 1.0 / dv.sqrt();
    vectors
        .iter()
        .map(|v| v.iter().map(|x| x * scale).collect())
        .collect()
}

/// φ(u) = erf(u/σ)/u, the potential of a unit Gaussian charge of width σ.
pub(super) fn smeared_coulomb(u: f64, sigma: f64) -> f64 {
    let y = u / sigma;
    if y < 1e-2 {
        let y2 = y * y;
        TWO_OVER_SQRT_PI / sigma * (1.0 - y2 / 3.0 + y2 * y2 / 10.0 - y2 * y2 * y2 / 42.0)
    } else {
        erf(y) / u
    }
}

/// φ'(u)/u for φ(u) = erf(u/σ)/u (finite at u = 0).
pub(super) fn smeared_coulomb_derivative_over_u(u: f64, sigma: f64) -> f64 {
    let y = u / sigma;
    let s3 = sigma * sigma * sigma;
    if y < 1e-2 {
        let y2 = y * y;
        TWO_OVER_SQRT_PI / s3 * (-2.0 / 3.0 + 0.4 * y2 - y2 * y2 / 7.0 + y2 * y2 * y2 / 27.0)
    } else {
        (TWO_OVER_SQRT_PI / sigma * (-y * y).exp() * u - erf(y)) / (u * u * u)
    }
}
