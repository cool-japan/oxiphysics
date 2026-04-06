//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{HBAR, KB};

/// One-dimensional diatomic harmonic chain.
///
/// Alternating masses m1 and m2 connected by springs of constant `k_spring`
/// with unit-cell spacing `a`.
#[derive(Debug, Clone)]
pub struct DiatomicChain {
    /// Mass of atom type 1 (must be ≤ m2 for conventional notation).
    pub m1: f64,
    /// Mass of atom type 2 (m2 ≥ m1).
    pub m2: f64,
    /// Spring constant.
    pub k_spring: f64,
    /// Lattice spacing (full unit-cell length).
    pub a: f64,
}
impl DiatomicChain {
    /// Construct a new diatomic chain.  No constraint on ordering of masses.
    pub fn new(m1: f64, m2: f64, k_spring: f64) -> Self {
        Self {
            m1,
            m2,
            k_spring,
            a: 1.0,
        }
    }
    /// Optical branch frequency at wave vector q.
    pub fn optical_frequency(&self, q: f64) -> f64 {
        let m1 = self.m1;
        let m2 = self.m2;
        let k = self.k_spring;
        let s = (q * self.a / 2.0).sin();
        let sum_inv = 1.0 / m1 + 1.0 / m2;
        let disc = sum_inv * sum_inv - 4.0 * s * s / (m1 * m2);
        let disc = disc.max(0.0);
        (k * (sum_inv + disc.sqrt())).sqrt()
    }
    /// Acoustic branch frequency at wave vector q.
    pub fn acoustic_frequency(&self, q: f64) -> f64 {
        let m1 = self.m1;
        let m2 = self.m2;
        let k = self.k_spring;
        let s = (q * self.a / 2.0).sin();
        let sum_inv = 1.0 / m1 + 1.0 / m2;
        let disc = sum_inv * sum_inv - 4.0 * s * s / (m1 * m2);
        let disc = disc.max(0.0);
        (k * (sum_inv - disc.sqrt())).max(0.0).sqrt()
    }
    /// Bandgap: difference between optical and acoustic branches at zone boundary.
    ///
    /// ω_op(π/a) - ω_ac(π/a)
    pub fn bandgap(&self) -> f64 {
        let q_bz = PI / self.a;
        let op = self.optical_frequency(q_bz);
        let ac = self.acoustic_frequency(q_bz);
        (op - ac).abs()
    }
}
/// Phonon dispersion along a high-symmetry path in the Brillouin zone.
///
/// Supports simple cubic and FCC Brillouin zone paths:
/// - Simple cubic: Γ(0,0,0) → X(π/a,0,0) → M(π/a,π/a,0) → Γ(0,0,0)
/// - FCC:          Γ(0,0,0) → X(2π/a,0,0) → W(2π/a,π/a,0) → Γ(0,0,0)
///
/// Each segment is sampled at `n_points` uniformly spaced k-points.
#[derive(Debug, Clone)]
pub struct PhononDispersion {
    /// Dynamical matrix used to compute frequencies.
    pub dynmat: DynamicalMatrix,
    /// Number of k-points per segment.
    pub n_points: usize,
    /// Computed k-points (each is \[kx, ky, kz\]).
    pub kpoints: Vec<[f64; 3]>,
    /// Phonon frequencies at each k-point (one Vec`f64` per k-point).
    pub frequencies: Vec<Vec<f64>>,
}
impl PhononDispersion {
    /// Compute the phonon dispersion along the simple-cubic Γ→X→M→Γ path.
    ///
    /// # Arguments
    /// * `dynmat`   – dynamical matrix
    /// * `n_points` – number of k-points per segment
    pub fn compute_sc(dynmat: DynamicalMatrix, n_points: usize) -> Self {
        let a = norm3(dynmat.lattice.a1).max(f64::EPSILON);
        let bz = PI / a;
        let gamma = [0.0_f64, 0.0, 0.0];
        let x = [bz, 0.0, 0.0];
        let m = [bz, bz, 0.0];
        let mut kpoints = Vec::new();
        let mut frequencies = Vec::new();
        for i in 0..n_points {
            let t = i as f64 / n_points as f64;
            let k = interpolate3(gamma, x, t);
            let freqs = dynmat.at_q(k);
            kpoints.push(k);
            frequencies.push(freqs);
        }
        for i in 0..n_points {
            let t = i as f64 / n_points as f64;
            let k = interpolate3(x, m, t);
            let freqs = dynmat.at_q(k);
            kpoints.push(k);
            frequencies.push(freqs);
        }
        for i in 0..=n_points {
            let t = i as f64 / n_points as f64;
            let k = interpolate3(m, gamma, t);
            let freqs = dynmat.at_q(k);
            kpoints.push(k);
            frequencies.push(freqs);
        }
        Self {
            dynmat,
            n_points,
            kpoints,
            frequencies,
        }
    }
    /// Maximum phonon frequency over the entire dispersion (rad/s or consistent unit).
    pub fn max_frequency(&self) -> f64 {
        self.frequencies
            .iter()
            .flat_map(|fs| fs.iter().cloned())
            .fold(0.0_f64, f64::max)
    }
    /// Minimum nonzero phonon frequency (useful for identifying acoustic branches).
    pub fn min_nonzero_frequency(&self) -> f64 {
        self.frequencies
            .iter()
            .flat_map(|fs| fs.iter().cloned())
            .filter(|&w| w > 1e-10)
            .fold(f64::INFINITY, f64::min)
    }
    /// All frequencies flattened into a single Vec.
    pub fn all_frequencies(&self) -> Vec<f64> {
        self.frequencies
            .iter()
            .flat_map(|v| v.iter().cloned())
            .collect()
    }
    /// Number of k-points in the path.
    pub fn n_kpoints(&self) -> usize {
        self.kpoints.len()
    }
    /// Check that acoustic modes vanish at Γ (q = 0).
    ///
    /// Returns true if frequencies at the first k-point are all < tolerance.
    pub fn acoustic_modes_at_gamma(&self, tol: f64) -> bool {
        if let Some(fs) = self.frequencies.first() {
            fs.iter().all(|&w| w < tol)
        } else {
            true
        }
    }
}
/// One-dimensional monatomic harmonic chain unit cell.
///
/// Atoms of equal mass `mass` are connected by springs of constant `spring_const`
/// with equilibrium spacing `equilibrium_spacing`.
#[derive(Debug, Clone)]
pub struct UnitCell1D {
    /// Atomic mass (amu or kg).
    pub mass: f64,
    /// Spring constant (N/m or eV/Å²).
    pub spring_const: f64,
    /// Equilibrium lattice spacing (Å or m).
    pub equilibrium_spacing: f64,
}
impl UnitCell1D {
    /// Construct a new 1D unit cell.
    pub fn new(mass: f64, spring_const: f64, equilibrium_spacing: f64) -> Self {
        Self {
            mass,
            spring_const,
            equilibrium_spacing,
        }
    }
    /// Acoustic dispersion relation for the 1D chain.
    ///
    /// ω(q) = 2 sqrt(k/m) |sin(q a / 2)|
    pub fn acoustic_dispersion(&self, q: f64) -> f64 {
        let omega_max = 2.0 * (self.spring_const / self.mass).sqrt();
        omega_max * (q * self.equilibrium_spacing / 2.0).sin().abs()
    }
    /// Group velocity ∂ω/∂q for the 1D acoustic branch.
    ///
    /// v_g = a sqrt(k/m) |cos(q a / 2)|
    pub fn group_velocity(&self, q: f64) -> f64 {
        let v0 = self.equilibrium_spacing * (self.spring_const / self.mass).sqrt();
        v0 * (q * self.equilibrium_spacing / 2.0).cos().abs()
    }
    /// Maximum frequency at the Brillouin-zone boundary (q = π/a).
    ///
    /// ω_BZ = 2 sqrt(k/m)
    pub fn zone_boundary_frequency(&self) -> f64 {
        2.0 * (self.spring_const / self.mass).sqrt()
    }
}
/// Mode Grüneisen parameters and macroscopic Grüneisen parameter.
///
/// The mode Grüneisen parameter characterises how the phonon frequency changes
/// with volume:
///
/// γ_i = −(V / ω_i) ∂ω_i / ∂V ≈ −(V / ω_i) (ω_i(V+δV) − ω_i(V−δV)) / (2δV)
///
/// The macroscopic Grüneisen parameter is the mode-average weighted by
/// the heat capacity:
///
/// γ = Σ_i γ_i C_i / Σ_i C_i
#[derive(Debug, Clone)]
pub struct GruneisenParameter {
    /// Phonon frequencies at volume V (rad/s).
    pub freqs_v: Vec<f64>,
    /// Phonon frequencies at volume V + δV.
    pub freqs_vp: Vec<f64>,
    /// Phonon frequencies at volume V − δV.
    pub freqs_vm: Vec<f64>,
    /// Reference volume V.
    pub volume: f64,
    /// Volume perturbation δV.
    pub delta_v: f64,
    /// Temperature for heat-capacity weighting (K).
    pub temperature: f64,
}
impl GruneisenParameter {
    /// Construct from frequency sets at three volumes.
    ///
    /// # Arguments
    /// * `freqs_v`     – frequencies at volume V
    /// * `freqs_vp`    – frequencies at volume V + δV
    /// * `freqs_vm`    – frequencies at volume V − δV
    /// * `volume`      – reference volume V
    /// * `delta_v`     – volume step δV
    /// * `temperature` – temperature (K)
    pub fn new(
        freqs_v: Vec<f64>,
        freqs_vp: Vec<f64>,
        freqs_vm: Vec<f64>,
        volume: f64,
        delta_v: f64,
        temperature: f64,
    ) -> Self {
        Self {
            freqs_v,
            freqs_vp,
            freqs_vm,
            volume,
            delta_v,
            temperature,
        }
    }
    /// Compute mode Grüneisen parameter γ_i for each mode.
    ///
    /// γ_i = −(V / ω_i) · (ω_{i,V+} − ω_{i,V−}) / (2 δV)
    pub fn mode_gruneisen(&self) -> Vec<f64> {
        let n = self
            .freqs_v
            .len()
            .min(self.freqs_vp.len())
            .min(self.freqs_vm.len());
        (0..n)
            .map(|i| {
                let wi = self.freqs_v[i];
                if wi.abs() < 1e-10 || self.delta_v.abs() < 1e-30 {
                    return 0.0;
                }
                let dw_dv = (self.freqs_vp[i] - self.freqs_vm[i]) / (2.0 * self.delta_v);
                -(self.volume / wi) * dw_dv
            })
            .collect()
    }
    /// Macroscopic Grüneisen parameter γ = Σ γ_i C_i / Σ C_i.
    ///
    /// Uses the Einstein single-mode heat capacity for weighting.
    pub fn macroscopic_gruneisen(&self) -> f64 {
        let mode_g = self.mode_gruneisen();
        let n = mode_g.len();
        if n == 0 {
            return 0.0;
        }
        let mut sum_gci = 0.0_f64;
        let mut sum_ci = 0.0_f64;
        for (i, &gamma_i) in mode_g.iter().enumerate() {
            let wi = self.freqs_v[i];
            if wi <= 0.0 || self.temperature <= 0.0 {
                continue;
            }
            let x = HBAR * wi / (KB * self.temperature);
            if x > 700.0 {
                continue;
            }
            let ex = x.exp();
            let denom = (ex - 1.0) * (ex - 1.0);
            if denom < 1e-300 {
                continue;
            }
            let ci = KB * x * x * ex / denom;
            sum_gci += gamma_i * ci;
            sum_ci += ci;
        }
        if sum_ci.abs() < 1e-300 {
            0.0
        } else {
            sum_gci / sum_ci
        }
    }
    /// Thermal pressure contribution dP/dT|_V = γ C_v / V.
    ///
    /// # Arguments
    /// * `cv` – total constant-volume heat capacity (J/K)
    pub fn thermal_pressure_coefficient(&self, cv: f64) -> f64 {
        if self.volume <= 0.0 {
            return 0.0;
        }
        self.macroscopic_gruneisen() * cv / self.volume
    }
}
/// Full 3N × 3N mass-weighted dynamical matrix from interatomic force constants.
///
/// The dynamical matrix D(q) is constructed from the force constant matrix Φ:
///
/// D_{αβ}^{ss'}(q) = (1/√(M_s M_s')) Σ_l Φ_{αβ}^{ss'}(0l) exp(i q · R_l)
///
/// where the sum is over lattice vectors R_l, s/s' are sublattice indices,
/// and α/β are Cartesian components.
///
/// For a mono-atomic lattice with isotropic spring constants this reduces to
/// a diagonal matrix whose diagonal elements are the phonon frequencies squared.
#[derive(Debug, Clone)]
pub struct DynamicalMatrix {
    /// Number of atoms per unit cell.
    pub n_atoms: usize,
    /// Atomic masses (kg or amu), length n_atoms.
    pub masses: Vec<f64>,
    /// Force constant matrix (3N × 3N, row-major flat array).
    pub force_constants: Vec<f64>,
    /// Lattice vectors.
    pub lattice: LatticeVectors,
}
impl DynamicalMatrix {
    /// Construct the dynamical matrix for a mono-atomic cubic crystal.
    ///
    /// Uses the nearest-neighbour isotropic spring constant `k_spring` and the
    /// given mass to build a simple diagonal Γ-point dynamical matrix.
    ///
    /// # Arguments
    /// * `k_spring`    – harmonic spring constant (energy / length²)
    /// * `mass`        – atomic mass
    /// * `coordination`– number of nearest neighbours (e.g. 6 for simple cubic)
    /// * `lattice`     – primitive lattice vectors
    pub fn new_cubic(
        k_spring: f64,
        mass: f64,
        coordination: usize,
        lattice: LatticeVectors,
    ) -> Self {
        let d = k_spring * coordination as f64 / mass;
        let mut fc = vec![0.0_f64; 9];
        fc[0] = d;
        fc[4] = d;
        fc[8] = d;
        Self {
            n_atoms: 1,
            masses: vec![mass],
            force_constants: fc,
            lattice,
        }
    }
    /// Compute eigenvalues (ω²) of the 3×3 Γ-point dynamical matrix via
    /// the diagonal approximation.
    ///
    /// Returns three phonon frequencies ω_i = √(max(0, D_ii)).
    pub fn phonon_frequencies_gamma(&self) -> Vec<f64> {
        if self.force_constants.len() < 9 {
            return vec![0.0; 3];
        }
        let fc = &self.force_constants;
        vec![
            fc[0].max(0.0).sqrt(),
            fc[4].max(0.0).sqrt(),
            fc[8].max(0.0).sqrt(),
        ]
    }
    /// Evaluate the (diagonal) dynamical matrix at wave vector q.
    ///
    /// For the isotropic cubic model the q-dependence enters through a
    /// structure factor:
    ///
    /// D_diag(q) = D_0 · (1 − cos(q · a)) averaged over nearest-neighbour bonds.
    ///
    /// This simplified formula handles the Γ-point (q → 0) limit correctly.
    ///
    /// # Arguments
    /// * `q` – wave vector as 3-vector (rad / lattice_parameter units)
    pub fn at_q(&self, q: [f64; 3]) -> Vec<f64> {
        if self.force_constants.len() < 9 {
            return vec![0.0; 3];
        }
        let a = norm3(self.lattice.a1).max(f64::EPSILON);
        let sx = (q[0] * a / 2.0).sin();
        let sy = (q[1] * a / 2.0).sin();
        let sz = (q[2] * a / 2.0).sin();
        let factor = (4.0 / 3.0) * (sx * sx + sy * sy + sz * sz);
        let d0 = self.force_constants[0];
        let freq_sq = (d0 * factor).max(0.0);
        vec![freq_sq.sqrt(), freq_sq.sqrt(), freq_sq.sqrt()]
    }
    /// Acoustic mode at Γ: ω(Γ) = 0 for acoustic branches.
    ///
    /// Verifies that the acoustic branch vanishes at q = 0.
    pub fn acoustic_at_gamma(&self) -> bool {
        let omegas = self.at_q([0.0, 0.0, 0.0]);
        omegas.iter().all(|&w| w.abs() < 1e-8)
    }
    /// Estimate the maximum phonon frequency (zone boundary).
    pub fn max_frequency(&self) -> f64 {
        self.phonon_frequencies_gamma()
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max)
    }
}
/// Phonon thermodynamic properties from a set of normal-mode frequencies.
///
/// Computes in the harmonic approximation:
/// - Zero-point energy E_ZPE = Σ_i ℏω_i / 2
/// - Internal energy U(T) = Σ_i ℏω_i \[n_BE(ω_i, T) + 1/2\]
/// - Heat capacity C_v(T) = ∂U/∂T
/// - Entropy S(T) = Σ_i k_B \[(n+1) ln(n+1) − n ln(n)\], n = n_BE(ω_i, T)
/// - Helmholtz free energy F(T) = U(T) − T S(T)
#[derive(Debug, Clone)]
pub struct ThermodynamicProperties {
    /// Phonon frequencies ω_i (rad/s).
    pub frequencies: Vec<f64>,
    /// Temperature (K).
    pub temperature: f64,
}
impl ThermodynamicProperties {
    /// Construct from a set of phonon frequencies and temperature.
    ///
    /// Acoustic modes near zero (ω < `tol`) are excluded automatically.
    ///
    /// # Arguments
    /// * `frequencies`  – phonon frequencies (rad/s)
    /// * `temperature`  – temperature (K)
    /// * `tol`          – frequency threshold below which modes are excluded (rad/s)
    pub fn new(frequencies: Vec<f64>, temperature: f64, tol: f64) -> Self {
        let filtered: Vec<f64> = frequencies.into_iter().filter(|&w| w > tol).collect();
        Self {
            frequencies: filtered,
            temperature,
        }
    }
    /// Zero-point energy E_ZPE = Σ_i ℏω_i / 2 (J per unit cell).
    pub fn zero_point_energy(&self) -> f64 {
        self.frequencies.iter().map(|&w| 0.5 * HBAR * w).sum()
    }
    /// Bose-Einstein occupation number for mode i.
    fn n_be(&self, omega: f64) -> f64 {
        if self.temperature <= 0.0 || omega <= 0.0 {
            return 0.0;
        }
        let x = HBAR * omega / (KB * self.temperature);
        if x > 700.0 {
            return 0.0;
        }
        1.0 / (x.exp() - 1.0)
    }
    /// Internal energy U(T) = Σ_i ℏω_i \[n(ω_i) + 1/2\] (J per unit cell).
    pub fn internal_energy(&self) -> f64 {
        self.frequencies
            .iter()
            .map(|&w| HBAR * w * (self.n_be(w) + 0.5))
            .sum()
    }
    /// Constant-volume heat capacity C_v(T) = ∂U/∂T (J/K per unit cell).
    ///
    /// C_v = Σ_i k_B x_i² exp(x_i) / (exp(x_i) − 1)²,  x_i = ℏω_i / k_B T
    pub fn heat_capacity(&self) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        self.frequencies
            .iter()
            .map(|&w| {
                let x = HBAR * w / (KB * self.temperature);
                if x > 700.0 || x <= 0.0 {
                    return 0.0;
                }
                let ex = x.exp();
                let denom = (ex - 1.0) * (ex - 1.0);
                if denom < 1e-300 {
                    return 0.0;
                }
                KB * x * x * ex / denom
            })
            .sum()
    }
    /// Entropy S(T) = Σ_i k_B \[(n+1) ln(n+1) − n ln(n)\] (J/K per unit cell).
    pub fn entropy(&self) -> f64 {
        if self.temperature <= 0.0 {
            return 0.0;
        }
        self.frequencies
            .iter()
            .map(|&w| {
                let n = self.n_be(w);
                if n < 1e-300 {
                    return 0.0;
                }
                KB * ((n + 1.0) * (n + 1.0).ln() - n * n.ln())
            })
            .sum()
    }
    /// Helmholtz free energy F(T) = U(T) − T S(T) (J per unit cell).
    pub fn helmholtz_free_energy(&self) -> f64 {
        self.internal_energy() - self.temperature * self.entropy()
    }
    /// Debye temperature estimate: θ_D = ℏ ω_max / k_B.
    pub fn debye_temperature(&self) -> f64 {
        let omega_max = self.frequencies.iter().cloned().fold(0.0_f64, f64::max);
        HBAR * omega_max / KB
    }
    /// Heat capacity in units of R (Dulong-Petit limit = 3).
    pub fn heat_capacity_per_r(&self) -> f64 {
        self.heat_capacity() / KB
    }
}
/// Quasi-harmonic approximation (QHA) and thermal expansion coefficient.
///
/// The QHA accounts for anharmonicity by computing the Helmholtz free energy
/// F(V, T) = E_static(V) + F_phonon(V, T) at multiple volumes and minimising
/// to find the equilibrium volume V₀(T).
///
/// The thermal expansion coefficient is:
/// α = −(1/V)(∂V/∂T)_P = γ C_v / (V K_T)
///
/// where K_T is the isothermal bulk modulus.
#[derive(Debug, Clone)]
pub struct AnharmonicCorrections {
    /// Volumes at which calculations were performed.
    pub volumes: Vec<f64>,
    /// Static (electronic) energies at each volume (J).
    pub static_energies: Vec<f64>,
    /// Phonon Helmholtz free energies at each volume and the stored temperature (J).
    pub phonon_free_energies: Vec<f64>,
    /// Temperature (K).
    pub temperature: f64,
    /// Macroscopic Grüneisen parameter (dimensionless).
    pub gruneisen: f64,
    /// Isothermal bulk modulus K_T (Pa).
    pub bulk_modulus: f64,
}
impl AnharmonicCorrections {
    /// Construct an AnharmonicCorrections object.
    ///
    /// # Arguments
    /// * `volumes`             – volumes (m³ per unit cell)
    /// * `static_energies`     – DFT total energies (J)
    /// * `phonon_free_energies`– phonon Helmholtz free energies (J)
    /// * `temperature`         – temperature (K)
    /// * `gruneisen`           – macroscopic Grüneisen parameter
    /// * `bulk_modulus`        – isothermal bulk modulus (Pa)
    pub fn new(
        volumes: Vec<f64>,
        static_energies: Vec<f64>,
        phonon_free_energies: Vec<f64>,
        temperature: f64,
        gruneisen: f64,
        bulk_modulus: f64,
    ) -> Self {
        Self {
            volumes,
            static_energies,
            phonon_free_energies,
            temperature,
            gruneisen,
            bulk_modulus,
        }
    }
    /// Total Helmholtz free energy F(V, T) = E_static + F_phonon at each volume (J).
    pub fn total_free_energies(&self) -> Vec<f64> {
        self.static_energies
            .iter()
            .zip(self.phonon_free_energies.iter())
            .map(|(e, f)| e + f)
            .collect()
    }
    /// Equilibrium volume V₀(T) from minimum of total F(V).
    ///
    /// Finds the volume corresponding to the minimum free energy by
    /// simple parabolic interpolation.
    pub fn equilibrium_volume(&self) -> f64 {
        let fv = self.total_free_energies();
        if fv.is_empty() {
            return 0.0;
        }
        let (min_idx, _) = fv
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap_or((0, &0.0));
        if min_idx > 0 && min_idx < fv.len() - 1 {
            let va = self.volumes[min_idx - 1];
            let vb = self.volumes[min_idx];
            let vc = self.volumes[min_idx + 1];
            let fa = fv[min_idx - 1];
            let fb = fv[min_idx];
            let fc = fv[min_idx + 1];
            let denom = 2.0 * (fa - 2.0 * fb + fc);
            if denom.abs() > 1e-300 {
                let _v0 = vb - (vc - va) * (fa - fb) / denom;
                return _v0.max(va).min(vc);
            }
        }
        self.volumes.get(min_idx).cloned().unwrap_or(0.0)
    }
    /// Linear thermal expansion coefficient α (K⁻¹).
    ///
    /// α = γ C_v / (3 V K_T)  (1/3 factor for linear vs volumetric)
    ///
    /// # Arguments
    /// * `cv` – heat capacity C_v at the reference temperature (J/K)
    pub fn thermal_expansion_coefficient(&self, cv: f64) -> f64 {
        let v0 = self.equilibrium_volume();
        if v0 <= 0.0 || self.bulk_modulus <= 0.0 {
            return 0.0;
        }
        self.gruneisen * cv / (3.0 * v0 * self.bulk_modulus)
    }
    /// Volumetric thermal expansion coefficient β = 3α (K⁻¹).
    ///
    /// # Arguments
    /// * `cv` – heat capacity C_v (J/K)
    pub fn volumetric_expansion(&self, cv: f64) -> f64 {
        3.0 * self.thermal_expansion_coefficient(cv)
    }
    /// Pressure correction ΔP(T) = −∂F/∂V at the equilibrium volume (Pa).
    ///
    /// Uses a finite-difference approximation over the stored volumes.
    pub fn pressure_correction(&self) -> f64 {
        let fv = self.total_free_energies();
        let n = fv.len();
        if n < 2 {
            return 0.0;
        }
        let mid = n / 2;
        if mid == 0 {
            return 0.0;
        }
        let df = fv[mid] - fv[mid - 1];
        let dv = self.volumes[mid] - self.volumes[mid - 1];
        if dv.abs() < 1e-30 {
            return 0.0;
        }
        -df / dv
    }
}
/// Primitive lattice vectors of a 3D crystal.
#[derive(Clone, Debug, PartialEq)]
pub struct LatticeVectors {
    /// First primitive vector a1.
    pub a1: [f64; 3],
    /// Second primitive vector a2.
    pub a2: [f64; 3],
    /// Third primitive vector a3.
    pub a3: [f64; 3],
}
impl LatticeVectors {
    /// Construct `LatticeVectors` from three primitive vectors.
    pub fn new(a1: [f64; 3], a2: [f64; 3], a3: [f64; 3]) -> Self {
        Self { a1, a2, a3 }
    }
    /// Simple cubic lattice with lattice parameter `a`.
    pub fn cubic(a: f64) -> Self {
        Self {
            a1: [a, 0.0, 0.0],
            a2: [0.0, a, 0.0],
            a3: [0.0, 0.0, a],
        }
    }
    /// Face-centred cubic lattice with conventional lattice parameter `a`.
    ///
    /// Primitive vectors for FCC:
    /// - a1 = a/2 (0, 1, 1)
    /// - a2 = a/2 (1, 0, 1)
    /// - a3 = a/2 (1, 1, 0)
    pub fn fcc(a: f64) -> Self {
        let h = a * 0.5;
        Self {
            a1: [0.0, h, h],
            a2: [h, 0.0, h],
            a3: [h, h, 0.0],
        }
    }
    /// Body-centred cubic lattice with conventional lattice parameter `a`.
    ///
    /// Primitive vectors for BCC:
    /// - a1 = a/2 (−1, 1, 1)
    /// - a2 = a/2 (1, −1, 1)
    /// - a3 = a/2 (1, 1, −1)
    pub fn bcc(a: f64) -> Self {
        let h = a * 0.5;
        Self {
            a1: [-h, h, h],
            a2: [h, -h, h],
            a3: [h, h, -h],
        }
    }
    /// Volume of the unit cell: `V = a1 · (a2 × a3)`.
    pub fn volume(&self) -> f64 {
        dot3(self.a1, cross3(self.a2, self.a3))
    }
    /// Compute reciprocal lattice vectors `b1`, `b2`, `b3`.
    ///
    /// `b1 = 2π (a2 × a3) / (a1 · a2 × a3)`
    pub fn reciprocal_vectors(&self) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let vol = self.volume();
        let inv = 2.0 * PI / vol;
        let b1 = scale3(cross3(self.a2, self.a3), inv);
        let b2 = scale3(cross3(self.a3, self.a1), inv);
        let b3 = scale3(cross3(self.a1, self.a2), inv);
        (b1, b2, b3)
    }
}
/// Phonon density of states computed by direct histogram sampling.
///
/// All phonon frequencies from the dispersion are binned into a histogram.
/// The result is normalised so that ∫ g(ω) dω = 1.
#[derive(Debug, Clone)]
pub struct PhononDOS {
    /// Histogram bin boundaries (n_bins + 1 values).
    pub bin_edges: Vec<f64>,
    /// Normalised DOS g(ω) in each bin.
    pub dos: Vec<f64>,
    /// Number of bins.
    pub n_bins: usize,
}
impl PhononDOS {
    /// Compute the phonon DOS from a flat list of frequencies.
    ///
    /// # Arguments
    /// * `frequencies` – slice of phonon frequencies
    /// * `n_bins`      – number of histogram bins
    pub fn from_frequencies(frequencies: &[f64], n_bins: usize) -> Self {
        if frequencies.is_empty() || n_bins == 0 {
            return Self {
                bin_edges: vec![0.0; n_bins + 1],
                dos: vec![0.0; n_bins],
                n_bins,
            };
        }
        let omega_max = frequencies.iter().cloned().fold(0.0_f64, f64::max);
        let omega_min = 0.0_f64;
        let dw = (omega_max - omega_min) / n_bins as f64;
        if dw <= 0.0 {
            return Self {
                bin_edges: vec![0.0; n_bins + 1],
                dos: vec![0.0; n_bins],
                n_bins,
            };
        }
        let mut hist = vec![0usize; n_bins];
        for &w in frequencies {
            if w < omega_min || w > omega_max {
                continue;
            }
            let bin = ((w - omega_min) / dw).floor() as usize;
            let bin = bin.min(n_bins - 1);
            hist[bin] += 1;
        }
        let total = frequencies.len() as f64;
        let dos: Vec<f64> = hist.iter().map(|&c| c as f64 / (total * dw)).collect();
        let bin_edges: Vec<f64> = (0..=n_bins).map(|i| omega_min + i as f64 * dw).collect();
        Self {
            bin_edges,
            dos,
            n_bins,
        }
    }
    /// Evaluate g(ω) at a given frequency by linear interpolation between bin centres.
    pub fn evaluate(&self, omega: f64) -> f64 {
        if self.n_bins == 0 || self.bin_edges.len() < 2 {
            return 0.0;
        }
        let dw = self.bin_edges[1] - self.bin_edges[0];
        if dw <= 0.0 {
            return 0.0;
        }
        let idx = ((omega - self.bin_edges[0]) / dw).floor() as isize;
        if idx < 0 || idx as usize >= self.n_bins {
            return 0.0;
        }
        self.dos[idx as usize]
    }
    /// Normalisation check: ∫ g(ω) dω ≈ 1.
    pub fn norm(&self) -> f64 {
        if self.bin_edges.len() < 2 {
            return 0.0;
        }
        let dw = self.bin_edges[1] - self.bin_edges[0];
        self.dos.iter().sum::<f64>() * dw
    }
    /// Maximum DOS value (peak of the spectrum).
    pub fn peak(&self) -> f64 {
        self.dos.iter().cloned().fold(0.0_f64, f64::max)
    }
    /// Frequency of the DOS peak (bin centre of the tallest bin).
    pub fn peak_frequency(&self) -> f64 {
        let idx = self
            .dos
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);
        if self.bin_edges.len() > idx + 1 {
            0.5 * (self.bin_edges[idx] + self.bin_edges[idx + 1])
        } else {
            0.0
        }
    }
}
/// Driven harmonic oscillator with damping.
///
/// Equation of motion: m ẍ + c ẋ + k x = F₀ cos(ωt)
#[derive(Debug, Clone)]
pub struct ForcedOscillator {
    /// Mass (kg or amu).
    pub mass: f64,
    /// Spring constant (N/m or eV/Å²).
    pub spring: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
}
impl ForcedOscillator {
    /// Construct a new forced oscillator.
    pub fn new(mass: f64, spring: f64, damping: f64) -> Self {
        Self {
            mass,
            spring,
            damping,
        }
    }
    /// Natural (undamped) angular frequency ω₀ = sqrt(k/m).
    pub fn natural_frequency(&self) -> f64 {
        if self.mass <= 0.0 {
            return 0.0;
        }
        (self.spring / self.mass).sqrt()
    }
    /// Damped angular frequency ω_d = sqrt(ω₀² − (c/2m)²).
    pub fn damped_frequency(&self) -> f64 {
        let w0 = self.natural_frequency();
        let gamma = self.damping / (2.0 * self.mass.max(f64::EPSILON));
        let d = w0 * w0 - gamma * gamma;
        d.max(0.0).sqrt()
    }
    /// Steady-state amplitude for harmonic forcing F₀ cos(ωt).
    ///
    /// X = F₀ / sqrt((k − mω²)² + (cω)²)
    pub fn steady_state_amplitude(&self, f0: f64, omega: f64) -> f64 {
        let k_eff = self.spring - self.mass * omega * omega;
        let c_eff = self.damping * omega;
        let denom = (k_eff * k_eff + c_eff * c_eff).sqrt();
        if denom < f64::EPSILON {
            return f64::INFINITY;
        }
        f0 / denom
    }
    /// Quality factor Q = m ω₀ / c = sqrt(mk) / c.
    pub fn quality_factor(&self) -> f64 {
        let denom = self.damping;
        if denom < f64::EPSILON {
            return f64::INFINITY;
        }
        (self.mass * self.spring).sqrt() / denom
    }
}
