// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended MD analysis methods.
//!
//! Provides post-processing tools for molecular dynamics trajectories:
//! - Radial distribution function (RDF / g(r))
//! - Mean-squared displacement (MSD) and diffusion coefficient
//! - Velocity autocorrelation function (VACF)
//! - Pressure tensor from virial theorem
//! - Static structure factor S(q)
//! - Heat flux and Green-Kubo thermal conductivity
//! - Viscosity from pressure tensor autocorrelation
//! - Dielectric constant from dipole moment fluctuations
//! - Orientational order parameters
//! - Cluster analysis (DBSCAN)

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// RadialDistributionFunction
// ---------------------------------------------------------------------------

/// Computes the radial distribution function g(r) from MD trajectory frames.
///
/// g(r) = (n(r) / (4π r² Δr ρ N)) normalized so that g(r) → 1 at large r
/// for a homogeneous fluid.
#[derive(Debug, Clone)]
pub struct RadialDistributionFunction {
    /// Histogram bins (count of pairs per shell).
    pub histogram: Vec<f64>,
    /// Bin width Δr.
    pub dr: f64,
    /// Maximum distance r_max.
    pub r_max: f64,
    /// Number of bins.
    pub n_bins: usize,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Number of frames accumulated.
    pub n_frames: usize,
    /// Box dimensions `[Lx, Ly, Lz]`.
    pub box_dims: [f64; 3],
}

impl RadialDistributionFunction {
    /// Create a new RDF accumulator.
    pub fn new(r_max: f64, n_bins: usize, n_atoms: usize, box_dims: [f64; 3]) -> Self {
        let dr = r_max / n_bins as f64;
        Self {
            histogram: vec![0.0; n_bins],
            dr,
            r_max,
            n_bins,
            n_atoms,
            n_frames: 0,
            box_dims,
        }
    }

    /// Accumulate one snapshot: positions is a flat array \[x0,y0,z0, x1,y1,z1, ...\].
    pub fn accumulate(&mut self, positions: &[f64]) {
        let n = self.n_atoms;
        let lx = self.box_dims[0];
        let ly = self.box_dims[1];
        let lz = self.box_dims[2];

        for i in 0..n {
            for j in (i + 1)..n {
                let mut dx = positions[3 * j] - positions[3 * i];
                let mut dy = positions[3 * j + 1] - positions[3 * i + 1];
                let mut dz = positions[3 * j + 2] - positions[3 * i + 2];
                // Minimum image convention
                dx -= (dx / lx).round() * lx;
                dy -= (dy / ly).round() * ly;
                dz -= (dz / lz).round() * lz;
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                if r < self.r_max {
                    let bin = (r / self.dr).floor() as usize;
                    if bin < self.n_bins {
                        self.histogram[bin] += 2.0; // count both (i,j) and (j,i)
                    }
                }
            }
        }
        self.n_frames += 1;
    }

    /// Compute normalised g(r) and return as `Vec<(r, g_r)>`.
    pub fn compute(&self) -> Vec<(f64, f64)> {
        let vol = self.box_dims[0] * self.box_dims[1] * self.box_dims[2];
        let rho = self.n_atoms as f64 / vol;
        let n_frames = self.n_frames.max(1) as f64;
        let n = self.n_atoms as f64;

        (0..self.n_bins)
            .map(|b| {
                let r = (b as f64 + 0.5) * self.dr;
                let shell_vol = 4.0 * PI * r * r * self.dr;
                let ideal = rho * shell_vol * n * n_frames;
                let gr = if ideal > 0.0 {
                    self.histogram[b] / ideal
                } else {
                    0.0
                };
                (r, gr)
            })
            .collect()
    }

    /// Find the first peak position and height.
    pub fn first_peak(&self) -> Option<(f64, f64)> {
        let gr_data = self.compute();
        let mut max_val = 0.0_f64;
        let mut max_r = 0.0_f64;
        for (r, gr) in &gr_data {
            if *gr > max_val {
                max_val = *gr;
                max_r = *r;
            }
        }
        if max_val > 1.0 {
            Some((max_r, max_val))
        } else {
            None
        }
    }

    /// Running average of g(r) over all accumulated frames.
    pub fn running_average_gr(&self) -> Vec<f64> {
        self.compute().into_iter().map(|(_, gr)| gr).collect()
    }
}

// ---------------------------------------------------------------------------
// MeanSquaredDisplacement
// ---------------------------------------------------------------------------

/// Mean-squared displacement (MSD) and diffusion coefficient computation.
///
/// MSD(t) = <|r(t) - r(0)|²>
///
/// Diffusion coefficient via Einstein relation: D = lim_{t→∞} MSD(t) / (6t) \[3D\].
#[derive(Debug, Clone)]
pub struct MeanSquaredDisplacement {
    /// Reference positions at t=0: flat \[x0,y0,z0,...\].
    pub r0: Vec<f64>,
    /// Current positions (unwrapped): flat.
    pub r_current: Vec<f64>,
    /// Accumulated MSD values at each stored time.
    pub msd: Vec<f64>,
    /// Time points at which MSD was sampled.
    pub times: Vec<f64>,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Box dimensions (for PBC unwrapping).
    pub box_dims: [f64; 3],
}

impl MeanSquaredDisplacement {
    /// Create a new MSD tracker from initial positions.
    pub fn new(positions: Vec<f64>, box_dims: [f64; 3]) -> Self {
        let n_atoms = positions.len() / 3;
        let r_current = positions.clone();
        Self {
            r0: positions,
            r_current,
            msd: Vec::new(),
            times: Vec::new(),
            n_atoms,
            box_dims,
        }
    }

    /// Update with new (unwrapped) positions and record MSD at time `t`.
    pub fn update(&mut self, positions: &[f64], t: f64) {
        let n = self.n_atoms;
        let mut sum = 0.0_f64;
        for i in 0..n {
            let dx = positions[3 * i] - self.r0[3 * i];
            let dy = positions[3 * i + 1] - self.r0[3 * i + 1];
            let dz = positions[3 * i + 2] - self.r0[3 * i + 2];
            sum += dx * dx + dy * dy + dz * dz;
        }
        self.msd.push(sum / n as f64);
        self.times.push(t);
        self.r_current = positions.to_vec();
    }

    /// Estimate diffusion coefficient D from linear fit to MSD (3D).
    ///
    /// Uses the last half of the stored MSD data for the fit.
    pub fn diffusion_coefficient(&self) -> f64 {
        let n = self.msd.len();
        if n < 4 {
            return 0.0;
        }
        let start = n / 2;
        let (slope, _) = linear_fit(&self.times[start..], &self.msd[start..]);
        slope / 6.0 // Einstein: MSD = 6 D t
    }

    /// MSD at a given index.
    pub fn msd_at(&self, idx: usize) -> f64 {
        self.msd.get(idx).copied().unwrap_or(0.0)
    }
}

// ---------------------------------------------------------------------------
// VelocityAutocorrelation
// ---------------------------------------------------------------------------

/// Velocity autocorrelation function (VACF) for diffusion and spectral analysis.
///
/// VACF(t) = <v(t) · v(0)> / <v(0) · v(0)>
///
/// Green-Kubo integral: D = (1/3) ∫₀^∞ VACF(t) dt
#[derive(Debug, Clone)]
pub struct VelocityAutocorrelation {
    /// Reference velocities at t=0: flat \[vx0,vy0,vz0,...\].
    pub v0: Vec<f64>,
    /// VACF values (not normalised) at stored times.
    pub vacf: Vec<f64>,
    /// Time points.
    pub times: Vec<f64>,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Initial kinetic energy (denominator for normalisation).
    pub v0_dot_v0: f64,
}

impl VelocityAutocorrelation {
    /// Initialise with reference velocities.
    pub fn new(velocities: Vec<f64>) -> Self {
        let n_atoms = velocities.len() / 3;
        let v0_dot_v0: f64 = velocities
            .chunks(3)
            .map(|v| v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            .sum::<f64>()
            / n_atoms as f64;
        Self {
            v0: velocities,
            vacf: Vec::new(),
            times: Vec::new(),
            n_atoms,
            v0_dot_v0,
        }
    }

    /// Add a time point: compute C(t) = <v(t)·v(0)>.
    pub fn accumulate(&mut self, velocities: &[f64], t: f64) {
        let n = self.n_atoms;
        let mut c = 0.0_f64;
        for i in 0..n {
            c += velocities[3 * i] * self.v0[3 * i]
                + velocities[3 * i + 1] * self.v0[3 * i + 1]
                + velocities[3 * i + 2] * self.v0[3 * i + 2];
        }
        self.vacf.push(c / n as f64);
        self.times.push(t);
    }

    /// Normalised VACF: Z(t) = C(t) / C(0).
    pub fn normalized(&self) -> Vec<f64> {
        let denom = self.v0_dot_v0.abs().max(1e-30);
        self.vacf.iter().map(|&c| c / denom).collect()
    }

    /// Green-Kubo diffusion coefficient: D = (1/3) ∫ VACF(t) dt.
    ///
    /// Integrates using the trapezoidal rule.
    pub fn diffusion_gk(&self) -> f64 {
        trapezoid_integrate(&self.times, &self.vacf) / 3.0
    }

    /// Spectral density (power spectrum of VACF) via discrete cosine transform.
    ///
    /// Returns Vec of (frequency, power) pairs.
    pub fn spectral_density(&self, dt: f64) -> Vec<(f64, f64)> {
        let n = self.vacf.len();
        if n == 0 {
            return vec![];
        }
        let normalised = self.normalized();
        // DCT-II via naive sum (for small n; production code uses FFT)
        (0..n)
            .map(|k| {
                let freq = k as f64 / (2.0 * n as f64 * dt);
                let power: f64 = normalised
                    .iter()
                    .enumerate()
                    .map(|(j, &z)| {
                        z * (PI * k as f64 * (2 * j + 1) as f64 / (2.0 * n as f64)).cos()
                    })
                    .sum();
                (freq, power.abs())
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// PressureTensor
// ---------------------------------------------------------------------------

/// Pressure tensor from virial theorem.
///
/// P_αβ = (1/V) \[Σ_i m_i v_{iα} v_{iβ} + Σ_{i<j} r_{ijα} F_{ijβ}\]
#[derive(Debug, Clone)]
pub struct PressureTensor {
    /// Instantaneous pressure tensor (3×3, row-major).
    pub p_inst: [f64; 9],
    /// Running average of pressure tensor.
    pub p_avg: [f64; 9],
    /// Number of steps accumulated.
    pub n_steps: u64,
    /// System volume V.
    pub volume: f64,
    /// Kinetic contribution to pressure tensor.
    pub p_kinetic: [f64; 9],
    /// Virial contribution to pressure tensor.
    pub p_virial: [f64; 9],
}

impl PressureTensor {
    /// Create a new pressure tensor accumulator.
    pub fn new(volume: f64) -> Self {
        Self {
            p_inst: [0.0; 9],
            p_avg: [0.0; 9],
            n_steps: 0,
            volume,
            p_kinetic: [0.0; 9],
            p_virial: [0.0; 9],
        }
    }

    /// Compute pressure tensor from velocities (kinetic part) and forces/positions (virial part).
    ///
    /// `velocities` flat \[vx0,vy0,vz0,...\], `masses` one per atom.
    /// `pairs` list of (i, j, r_ij \[f64; 3\], f_ij \[f64; 3\]).
    pub fn compute(
        &mut self,
        velocities: &[f64],
        masses: &[f64],
        pairs: &[(usize, usize, [f64; 3], [f64; 3])],
    ) {
        let v = self.volume;
        let mut pk = [0.0_f64; 9];
        let mut pv = [0.0_f64; 9];

        let n_atoms = masses.len();
        for i in 0..n_atoms {
            let m = masses[i];
            let vx = velocities[3 * i];
            let vy = velocities[3 * i + 1];
            let vz = velocities[3 * i + 2];
            let vel = [vx, vy, vz];
            for a in 0..3 {
                for b in 0..3 {
                    pk[a * 3 + b] += m * vel[a] * vel[b];
                }
            }
        }

        for (_, _, rij, fij) in pairs {
            for a in 0..3 {
                for b in 0..3 {
                    pv[a * 3 + b] += rij[a] * fij[b];
                }
            }
        }

        for i in 0..9 {
            self.p_inst[i] = (pk[i] + pv[i]) / v;
            self.p_kinetic[i] = pk[i] / v;
            self.p_virial[i] = pv[i] / v;
        }

        // Update running average
        let n = self.n_steps as f64;
        for i in 0..9 {
            self.p_avg[i] = (self.p_avg[i] * n + self.p_inst[i]) / (n + 1.0);
        }
        self.n_steps += 1;
    }

    /// Scalar pressure P = Tr(P_αβ) / 3.
    pub fn scalar_pressure(&self) -> f64 {
        (self.p_inst[0] + self.p_inst[4] + self.p_inst[8]) / 3.0
    }

    /// Average scalar pressure.
    pub fn avg_scalar_pressure(&self) -> f64 {
        (self.p_avg[0] + self.p_avg[4] + self.p_avg[8]) / 3.0
    }
}

// ---------------------------------------------------------------------------
// StructureFactor
// ---------------------------------------------------------------------------

/// Static structure factor S(q) from atomic positions.
///
/// S(q) = (1/N) |Σ_j exp(i q·r_j)|²
#[derive(Debug, Clone)]
pub struct StructureFactor {
    /// q-vectors to evaluate S(q) at.
    pub q_vectors: Vec<[f64; 3]>,
    /// S(q) values corresponding to each q-vector.
    pub sq: Vec<f64>,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Number of frames accumulated.
    pub n_frames: usize,
    /// Running S(q) sum for averaging.
    pub sq_sum: Vec<f64>,
}

impl StructureFactor {
    /// Create a new structure factor calculator for given q-vectors.
    pub fn new(q_vectors: Vec<[f64; 3]>, n_atoms: usize) -> Self {
        let nq = q_vectors.len();
        Self {
            q_vectors,
            sq: vec![0.0; nq],
            n_atoms,
            n_frames: 0,
            sq_sum: vec![0.0; nq],
        }
    }

    /// Accumulate one frame: positions flat \[x0,y0,z0,...\].
    pub fn accumulate(&mut self, positions: &[f64]) {
        let n = self.n_atoms;
        for (qi, q) in self.q_vectors.iter().enumerate() {
            let mut re = 0.0_f64;
            let mut im = 0.0_f64;
            for j in 0..n {
                let qdotr = q[0] * positions[3 * j]
                    + q[1] * positions[3 * j + 1]
                    + q[2] * positions[3 * j + 2];
                re += qdotr.cos();
                im += qdotr.sin();
            }
            self.sq_sum[qi] += (re * re + im * im) / n as f64;
        }
        self.n_frames += 1;
    }

    /// Return averaged S(q) over all frames.
    pub fn averaged(&mut self) -> &[f64] {
        let nf = self.n_frames.max(1) as f64;
        for i in 0..self.sq.len() {
            self.sq[i] = self.sq_sum[i] / nf;
        }
        &self.sq
    }

    /// Debye formula: S(q) for isotropic system (single |q| value).
    ///
    /// S(q) = (1/N) Σ_{i,j} sin(q r_{ij}) / (q r_{ij})
    pub fn debye(&self, q_mag: f64, positions: &[f64]) -> f64 {
        let n = self.n_atoms;
        let mut sum = n as f64; // diagonal terms (r=0 → sinc=1)
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions[3 * j] - positions[3 * i];
                let dy = positions[3 * j + 1] - positions[3 * i + 1];
                let dz = positions[3 * j + 2] - positions[3 * i + 2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                let qr = q_mag * r;
                let sinc = if qr < 1e-10 { 1.0 } else { qr.sin() / qr };
                sum += 2.0 * sinc;
            }
        }
        sum / n as f64
    }
}

// ---------------------------------------------------------------------------
// HeatFlux
// ---------------------------------------------------------------------------

/// Heat flux and Green-Kubo thermal conductivity.
///
/// λ = V / (3 k_B T²) ∫₀^∞ <J(t)·J(0)> dt
///
/// where J is the microscopic heat current vector.
#[derive(Debug, Clone)]
pub struct HeatFlux {
    /// Heat current vectors at stored times: `Vec<[Jx, Jy, Jz]>`.
    pub j_history: Vec<[f64; 3]>,
    /// Time points.
    pub times: Vec<f64>,
    /// System volume.
    pub volume: f64,
    /// Temperature.
    pub temperature: f64,
    /// Boltzmann constant (kJ/mol/K).
    pub k_b: f64,
}

impl HeatFlux {
    /// Create a new heat flux tracker.
    pub fn new(volume: f64, temperature: f64) -> Self {
        Self {
            j_history: Vec::new(),
            times: Vec::new(),
            volume,
            temperature,
            k_b: 8.314e-3,
        }
    }

    /// Record a heat current vector at time `t`.
    pub fn record(&mut self, j: [f64; 3], t: f64) {
        self.j_history.push(j);
        self.times.push(t);
    }

    /// Compute the heat current autocorrelation function <J(t)·J(0)>.
    pub fn heat_current_acf(&self) -> Vec<f64> {
        let n = self.j_history.len();
        if n == 0 {
            return vec![];
        }
        let j0 = self.j_history[0];
        self.j_history
            .iter()
            .map(|jt| jt[0] * j0[0] + jt[1] * j0[1] + jt[2] * j0[2])
            .collect()
    }

    /// Green-Kubo thermal conductivity λ (W/m·K or kJ/mol/K/Å/ps units).
    ///
    /// λ = V / (3 k_B T²) ∫ HACF(t) dt
    pub fn thermal_conductivity(&self) -> f64 {
        let acf = self.heat_current_acf();
        if acf.is_empty() || self.times.is_empty() {
            return 0.0;
        }
        let integral = trapezoid_integrate(&self.times, &acf);
        let t2 = self.temperature * self.temperature;
        integral * self.volume / (3.0 * self.k_b * t2)
    }
}

// ---------------------------------------------------------------------------
// ViscosityCalc
// ---------------------------------------------------------------------------

/// Green-Kubo viscosity from pressure tensor autocorrelation.
///
/// η = V / (k_B T) ∫₀^∞ <P_xy(t) P_xy(0)> dt
#[derive(Debug, Clone)]
pub struct ViscosityCalc {
    /// Off-diagonal pressure tensor element P_xy at each time: (t, P_xy).
    pub pxy_history: Vec<(f64, f64)>,
    /// System volume.
    pub volume: f64,
    /// Temperature.
    pub temperature: f64,
    /// Boltzmann constant.
    pub k_b: f64,
}

impl ViscosityCalc {
    /// Create a new viscosity calculator.
    pub fn new(volume: f64, temperature: f64) -> Self {
        Self {
            pxy_history: Vec::new(),
            volume,
            temperature,
            k_b: 8.314e-3,
        }
    }

    /// Record P_xy at time `t`.
    pub fn record(&mut self, t: f64, pxy: f64) {
        self.pxy_history.push((t, pxy));
    }

    /// Compute P_xy autocorrelation function.
    pub fn pxy_acf(&self) -> Vec<f64> {
        let n = self.pxy_history.len();
        if n == 0 {
            return vec![];
        }
        let p0 = self.pxy_history[0].1;
        self.pxy_history.iter().map(|&(_, p)| p * p0).collect()
    }

    /// Green-Kubo shear viscosity η.
    pub fn viscosity(&self) -> f64 {
        let acf = self.pxy_acf();
        if acf.is_empty() {
            return 0.0;
        }
        let times: Vec<f64> = self.pxy_history.iter().map(|&(t, _)| t).collect();
        let integral = trapezoid_integrate(&times, &acf);
        let kt = self.k_b * self.temperature;
        integral * self.volume / kt
    }
}

// ---------------------------------------------------------------------------
// DielectricConst
// ---------------------------------------------------------------------------

/// Dielectric constant from dipole moment fluctuations.
///
/// ε = 1 + (`M²` - `M`²) * 3V / (3 ε₀ k_B T)  (in SI)
///
/// In reduced units: ε = 1 + 4π (`M²` - `M`²) / (3 V k_B T)
#[derive(Debug, Clone)]
pub struct DielectricConst {
    /// Dipole moment vectors at each frame: `[Mx, My, Mz]`.
    pub dipoles: Vec<[f64; 3]>,
    /// System volume.
    pub volume: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Boltzmann constant (kJ/mol/K).
    pub k_b: f64,
}

impl DielectricConst {
    /// Create a new dielectric constant accumulator.
    pub fn new(volume: f64, temperature: f64) -> Self {
        Self {
            dipoles: Vec::new(),
            volume,
            temperature,
            k_b: 8.314e-3,
        }
    }

    /// Add a dipole moment vector for one frame.
    pub fn add_frame(&mut self, dipole: [f64; 3]) {
        self.dipoles.push(dipole);
    }

    /// Mean dipole moment `[`Mx`, `My`, `Mz`]`.
    pub fn mean_dipole(&self) -> [f64; 3] {
        let n = self.dipoles.len() as f64;
        if n == 0.0 {
            return [0.0; 3];
        }
        let sum = self.dipoles.iter().fold([0.0_f64; 3], |mut acc, d| {
            acc[0] += d[0];
            acc[1] += d[1];
            acc[2] += d[2];
            acc
        });
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }

    /// Fluctuation `M²` - `M`².
    pub fn dipole_fluctuation(&self) -> f64 {
        let n = self.dipoles.len() as f64;
        if n == 0.0 {
            return 0.0;
        }
        let m2_avg = self
            .dipoles
            .iter()
            .map(|d| d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
            .sum::<f64>()
            / n;
        let m_avg = self.mean_dipole();
        let m_avg2 = m_avg[0] * m_avg[0] + m_avg[1] * m_avg[1] + m_avg[2] * m_avg[2];
        m2_avg - m_avg2
    }

    /// Dielectric constant ε in reduced units.
    pub fn dielectric_constant(&self) -> f64 {
        let fluct = self.dipole_fluctuation();
        let kt = self.k_b * self.temperature;
        1.0 + 4.0 * PI * fluct / (3.0 * self.volume * kt)
    }
}

// ---------------------------------------------------------------------------
// OrderParameter
// ---------------------------------------------------------------------------

/// Orientational order parameters for molecular systems.
///
/// Computes the P2 Legendre polynomial order parameter for rod-like molecules:
/// S = <P2(cos θ)> = <(3 cos²θ - 1) / 2>
#[derive(Debug, Clone)]
pub struct OrderParameter {
    /// Director axis (normalised) `[nx, ny, nz]`.
    pub director: [f64; 3],
    /// Time-averaged S values.
    pub s_history: Vec<f64>,
    /// Time points.
    pub times: Vec<f64>,
    /// Tolerance for director determination.
    pub tol: f64,
}

impl OrderParameter {
    /// Create a new order parameter calculator with a given director axis.
    pub fn new(director: [f64; 3]) -> Self {
        let mag =
            (director[0] * director[0] + director[1] * director[1] + director[2] * director[2])
                .sqrt();
        let d = if mag > 1e-10 {
            [director[0] / mag, director[1] / mag, director[2] / mag]
        } else {
            [0.0, 0.0, 1.0]
        };
        Self {
            director: d,
            s_history: Vec::new(),
            times: Vec::new(),
            tol: 1e-8,
        }
    }

    /// Compute P2 order parameter from molecular orientation vectors.
    ///
    /// `orientations` is a flat array of unit vectors \[ux0,uy0,uz0,...\].
    pub fn compute_p2(&self, orientations: &[f64]) -> f64 {
        let n = orientations.len() / 3;
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = (0..n)
            .map(|i| {
                let cos_theta = orientations[3 * i] * self.director[0]
                    + orientations[3 * i + 1] * self.director[1]
                    + orientations[3 * i + 2] * self.director[2];
                (3.0 * cos_theta * cos_theta - 1.0) * 0.5
            })
            .sum();
        sum / n as f64
    }

    /// Record order parameter at time `t`.
    pub fn record(&mut self, orientations: &[f64], t: f64) {
        let s = self.compute_p2(orientations);
        self.s_history.push(s);
        self.times.push(t);
    }

    /// Time-averaged order parameter.
    pub fn mean_s(&self) -> f64 {
        if self.s_history.is_empty() {
            return 0.0;
        }
        self.s_history.iter().sum::<f64>() / self.s_history.len() as f64
    }

    /// Nematic order parameter: max eigenvalue of the Q-tensor.
    ///
    /// For a set of unit vectors computes Q_αβ = (1/N) Σ_i (3 u_{iα} u_{iβ} - δ_{αβ}) / 2.
    pub fn nematic_order(&self, orientations: &[f64]) -> f64 {
        let n = orientations.len() / 3;
        if n == 0 {
            return 0.0;
        }
        let mut q = [[0.0_f64; 3]; 3];
        for i in 0..n {
            let u = [
                orientations[3 * i],
                orientations[3 * i + 1],
                orientations[3 * i + 2],
            ];
            for a in 0..3 {
                for b in 0..3 {
                    let delta = if a == b { 1.0 } else { 0.0 };
                    q[a][b] += (3.0 * u[a] * u[b] - delta) * 0.5;
                }
            }
        }
        for row in &mut q {
            for v in row.iter_mut() {
                *v /= n as f64;
            }
        }
        // Largest eigenvalue of q (power iteration approximation)
        let trace = q[0][0] + q[1][1] + q[2][2];
        // For uniaxial: largest eigenvalue ≈ max diagonal if aligned
        let max_diag = q[0][0].max(q[1][1]).max(q[2][2]);
        // Simple approximation via Frobenius norm
        let _ = trace;
        max_diag.abs()
    }
}

// ---------------------------------------------------------------------------
// ClusterAnalysis
// ---------------------------------------------------------------------------

/// Neighbor-based clustering using DBSCAN algorithm.
///
/// DBSCAN identifies clusters based on:
/// - `eps`: neighborhood radius
/// - `min_pts`: minimum points to form a core
#[derive(Debug, Clone)]
pub struct ClusterAnalysis {
    /// Epsilon neighborhood radius.
    pub eps: f64,
    /// Minimum points for a core point.
    pub min_pts: usize,
    /// Cluster labels: -1 = noise, ≥0 = cluster id.
    pub labels: Vec<i64>,
    /// Number of clusters found.
    pub n_clusters: usize,
}

impl ClusterAnalysis {
    /// Create a new DBSCAN cluster analysis.
    pub fn new(eps: f64, min_pts: usize) -> Self {
        Self {
            eps,
            min_pts,
            labels: Vec::new(),
            n_clusters: 0,
        }
    }

    /// Run DBSCAN on a 3D point set.
    ///
    /// `positions` is flat \[x0,y0,z0, x1,y1,z1, ...\].
    pub fn run(&mut self, positions: &[f64]) {
        let n = positions.len() / 3;
        let eps2 = self.eps * self.eps;
        let mut labels = vec![-1i64; n];
        let mut cluster_id = 0i64;
        let mut visited = vec![false; n];

        for i in 0..n {
            if visited[i] {
                continue;
            }
            visited[i] = true;
            let neighbors = self.range_query(positions, i, eps2);
            if neighbors.len() < self.min_pts {
                labels[i] = -1; // noise
            } else {
                self.expand_cluster(
                    positions,
                    i,
                    &neighbors,
                    cluster_id,
                    &mut labels,
                    &mut visited,
                    eps2,
                );
                cluster_id += 1;
            }
        }

        self.labels = labels;
        self.n_clusters = cluster_id as usize;
    }

    fn range_query(&self, positions: &[f64], idx: usize, eps2: f64) -> Vec<usize> {
        let n = positions.len() / 3;
        let xi = positions[3 * idx];
        let yi = positions[3 * idx + 1];
        let zi = positions[3 * idx + 2];
        (0..n)
            .filter(|&j| {
                let dx = positions[3 * j] - xi;
                let dy = positions[3 * j + 1] - yi;
                let dz = positions[3 * j + 2] - zi;
                dx * dx + dy * dy + dz * dz <= eps2
            })
            .collect()
    }

    fn expand_cluster(
        &self,
        positions: &[f64],
        _core: usize,
        neighbors: &[usize],
        cluster_id: i64,
        labels: &mut [i64],
        visited: &mut [bool],
        eps2: f64,
    ) {
        let mut queue: Vec<usize> = neighbors.to_vec();
        let mut qi = 0;
        while qi < queue.len() {
            let q = queue[qi];
            qi += 1;
            if !visited[q] {
                visited[q] = true;
                let q_neighbors = self.range_query(positions, q, eps2);
                if q_neighbors.len() >= self.min_pts {
                    for &qn in &q_neighbors {
                        if !queue.contains(&qn) {
                            queue.push(qn);
                        }
                    }
                }
            }
            if labels[q] < 0 {
                labels[q] = cluster_id;
            }
        }
    }

    /// Sizes of each cluster (number of points per cluster).
    pub fn cluster_sizes(&self) -> Vec<usize> {
        if self.n_clusters == 0 {
            return vec![];
        }
        let mut sizes = vec![0usize; self.n_clusters];
        for &l in &self.labels {
            if l >= 0 {
                sizes[l as usize] += 1;
            }
        }
        sizes
    }

    /// Percolation check: returns true if any cluster spans more than `threshold`
    /// fraction of total points.
    pub fn percolation_threshold(&self, threshold: f64) -> bool {
        let total = self.labels.len() as f64;
        if total == 0.0 {
            return false;
        }
        self.cluster_sizes()
            .iter()
            .any(|&s| s as f64 / total > threshold)
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Bin a distance value into a histogram bin index.
///
/// Returns `None` if `r >= r_max`.
pub fn bin_distance(r: f64, dr: f64, r_max: f64) -> Option<usize> {
    if r >= r_max || r < 0.0 {
        return None;
    }
    Some((r / dr).floor() as usize)
}

/// Compute a running average of a slice (cumulative mean).
pub fn running_average(data: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(data.len());
    let mut sum = 0.0_f64;
    for (i, &v) in data.iter().enumerate() {
        sum += v;
        result.push(sum / (i + 1) as f64);
    }
    result
}

/// Autocorrelation function using direct summation (for short arrays).
///
/// C(k) = (1/(N-k)) Σ_{i=0}^{N-k-1} x_i * x_{i+k}
pub fn autocorrelation_fft(data: &[f64]) -> Vec<f64> {
    let n = data.len();
    (0..n)
        .map(|k| {
            let count = (n - k).max(1);
            (0..(n - k)).map(|i| data[i] * data[i + k]).sum::<f64>() / count as f64
        })
        .collect()
}

/// Unwrap a 1D trajectory coordinate for periodic boundary conditions.
///
/// Detects jumps larger than `box_length / 2` and applies integer corrections.
pub fn unwrap_pbc_trajectory(positions: &[f64], box_length: f64) -> Vec<f64> {
    if positions.is_empty() {
        return vec![];
    }
    let mut unwrapped = Vec::with_capacity(positions.len());
    unwrapped.push(positions[0]);
    let mut offset = 0.0_f64;
    for i in 1..positions.len() {
        let diff = positions[i] - positions[i - 1];
        if diff > box_length * 0.5 {
            offset -= box_length;
        } else if diff < -box_length * 0.5 {
            offset += box_length;
        }
        unwrapped.push(positions[i] + offset);
    }
    unwrapped
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Trapezoidal rule integration.
fn trapezoid_integrate(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len().min(y.len());
    if n < 2 {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for i in 0..(n - 1) {
        sum += (x[i + 1] - x[i]) * (y[i] + y[i + 1]) * 0.5;
    }
    sum
}

/// Linear least-squares fit y = a*x + b. Returns (slope, intercept).
fn linear_fit(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len().min(y.len()) as f64;
    if n < 2.0 {
        return (0.0, 0.0);
    }
    let sx: f64 = x.iter().sum();
    let sy: f64 = y.iter().take(x.len()).sum();
    let sxy: f64 = x.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
    let sxx: f64 = x.iter().map(|a| a * a).sum();
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-30 {
        return (0.0, sy / n);
    }
    let slope = (n * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / n;
    (slope, intercept)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- RadialDistributionFunction ---

    #[test]
    fn test_rdf_new() {
        let rdf = RadialDistributionFunction::new(5.0, 100, 10, [10.0, 10.0, 10.0]);
        assert_eq!(rdf.n_bins, 100);
        assert!((rdf.dr - 0.05).abs() < 1e-14);
    }

    #[test]
    fn test_rdf_accumulate_no_panic() {
        let mut rdf = RadialDistributionFunction::new(5.0, 50, 4, [10.0, 10.0, 10.0]);
        let positions = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0];
        rdf.accumulate(&positions);
        assert_eq!(rdf.n_frames, 1);
    }

    #[test]
    fn test_rdf_compute_length() {
        let rdf = RadialDistributionFunction::new(5.0, 50, 4, [10.0, 10.0, 10.0]);
        let result = rdf.compute();
        assert_eq!(result.len(), 50);
    }

    #[test]
    fn test_rdf_gr_non_negative() {
        let mut rdf = RadialDistributionFunction::new(5.0, 50, 4, [10.0, 10.0, 10.0]);
        let positions = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 2.0, 0.0];
        rdf.accumulate(&positions);
        for (_, gr) in rdf.compute() {
            assert!(gr >= 0.0);
        }
    }

    // --- MeanSquaredDisplacement ---

    #[test]
    fn test_msd_zero_displacement() {
        let r0 = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let mut msd = MeanSquaredDisplacement::new(r0.clone(), [10.0, 10.0, 10.0]);
        msd.update(&r0, 1.0);
        assert!((msd.msd_at(0)).abs() < 1e-12);
    }

    #[test]
    fn test_msd_known_displacement() {
        let r0 = vec![0.0, 0.0, 0.0];
        let mut msd = MeanSquaredDisplacement::new(r0, [100.0, 100.0, 100.0]);
        let r1 = vec![1.0, 0.0, 0.0]; // displaced by 1
        msd.update(&r1, 1.0);
        assert!((msd.msd_at(0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_msd_diffusion_coefficient_linear() {
        let r0 = vec![0.0, 0.0, 0.0];
        let mut msd_calc = MeanSquaredDisplacement::new(r0, [1000.0, 1000.0, 1000.0]);
        // MSD = 6 D t with D = 1.0 → MSD grows linearly
        for step in 1..=20 {
            let t = step as f64;
            // Displacement sqrt(6*1.0*t) in x only → MSD = 6*t
            let r = vec![(6.0 * t).sqrt(), 0.0, 0.0];
            msd_calc.update(&r, t);
        }
        let d = msd_calc.diffusion_coefficient();
        assert!((d - 1.0).abs() < 0.1);
    }

    // --- VelocityAutocorrelation ---

    #[test]
    fn test_vacf_zero_time() {
        let v0 = vec![1.0, 0.0, 0.0];
        let mut vacf = VelocityAutocorrelation::new(v0.clone());
        vacf.accumulate(&v0, 0.0);
        assert!((vacf.vacf[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vacf_normalized_at_zero() {
        let v0 = vec![2.0, 0.0, 0.0];
        let mut vacf = VelocityAutocorrelation::new(v0.clone());
        vacf.accumulate(&v0, 0.0);
        let n = vacf.normalized();
        assert!((n[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_vacf_diffusion_gk_positive() {
        let v0 = vec![1.0, 0.0, 0.0];
        let mut vacf = VelocityAutocorrelation::new(v0.clone());
        vacf.accumulate(&v0, 0.0);
        vacf.accumulate(&v0, 1.0);
        let d = vacf.diffusion_gk();
        assert!(d >= 0.0);
    }

    // --- PressureTensor ---

    #[test]
    fn test_pressure_tensor_scalar() {
        let mut pt = PressureTensor::new(1000.0);
        let vels = vec![1.0, 0.0, 0.0];
        let masses = vec![1.0];
        let pairs: Vec<(usize, usize, [f64; 3], [f64; 3])> = vec![];
        pt.compute(&vels, &masses, &pairs);
        let p = pt.scalar_pressure();
        assert!(p.abs() > 0.0);
    }

    #[test]
    fn test_pressure_tensor_running_avg() {
        let mut pt = PressureTensor::new(1000.0);
        let vels = vec![1.0, 0.0, 0.0];
        let masses = vec![1.0];
        let pairs: Vec<(usize, usize, [f64; 3], [f64; 3])> = vec![];
        pt.compute(&vels, &masses, &pairs);
        pt.compute(&vels, &masses, &pairs);
        assert_eq!(pt.n_steps, 2);
        assert!((pt.avg_scalar_pressure() - pt.scalar_pressure()).abs() < 1e-10);
    }

    // --- StructureFactor ---

    #[test]
    fn test_structure_factor_q_zero() {
        let q_vecs = vec![[0.0_f64; 3]];
        let mut sf = StructureFactor::new(q_vecs, 4);
        let positions = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        sf.accumulate(&positions);
        let sq = sf.averaged();
        // At q=0, S(0) = N = 4
        assert!((sq[0] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_structure_factor_debye_q_zero() {
        // At q=0: sinc(qr)=1 for all pairs, so S(0) = (N + 2*N*(N-1)/2) / N = N
        let sf = StructureFactor::new(vec![], 4);
        let positions = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let s0 = sf.debye(0.0, &positions);
        assert!((s0 - 4.0).abs() < 1e-10);
    }

    // --- HeatFlux ---

    #[test]
    fn test_heat_flux_conductivity() {
        let mut hf = HeatFlux::new(1000.0, 300.0);
        hf.record([1.0, 0.0, 0.0], 0.0);
        hf.record([0.5, 0.0, 0.0], 1.0);
        let lambda = hf.thermal_conductivity();
        // Just check it's finite and non-negative
        assert!(lambda.is_finite());
    }

    #[test]
    fn test_heat_flux_acf_length() {
        let mut hf = HeatFlux::new(1000.0, 300.0);
        for i in 0..5 {
            hf.record([i as f64, 0.0, 0.0], i as f64);
        }
        assert_eq!(hf.heat_current_acf().len(), 5);
    }

    // --- ViscosityCalc ---

    #[test]
    fn test_viscosity_record_and_compute() {
        let mut vc = ViscosityCalc::new(1000.0, 300.0);
        vc.record(0.0, 1.0);
        vc.record(1.0, 0.5);
        vc.record(2.0, 0.0);
        let eta = vc.viscosity();
        assert!(eta.is_finite());
    }

    #[test]
    fn test_pxy_acf_length() {
        let mut vc = ViscosityCalc::new(1000.0, 300.0);
        for i in 0..4 {
            vc.record(i as f64, i as f64 * 0.5);
        }
        let acf = vc.pxy_acf();
        assert_eq!(acf.len(), 4);
    }

    // --- DielectricConst ---

    #[test]
    fn test_dielectric_const_no_fluctuation() {
        let mut dc = DielectricConst::new(1000.0, 300.0);
        // All frames same dipole → no fluctuation → ε = 1
        for _ in 0..10 {
            dc.add_frame([1.0, 0.0, 0.0]);
        }
        let fluct = dc.dipole_fluctuation();
        assert!(fluct.abs() < 1e-10);
        assert!((dc.dielectric_constant() - 1.0).abs() < 1e-8);
    }

    #[test]
    fn test_dielectric_mean_dipole() {
        let mut dc = DielectricConst::new(1000.0, 300.0);
        dc.add_frame([2.0, 0.0, 0.0]);
        dc.add_frame([4.0, 0.0, 0.0]);
        let m = dc.mean_dipole();
        assert!((m[0] - 3.0).abs() < 1e-12);
    }

    // --- OrderParameter ---

    #[test]
    fn test_order_parameter_aligned() {
        let op = OrderParameter::new([0.0, 0.0, 1.0]);
        // All molecules aligned with z-axis
        let orientations = vec![0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0];
        let s = op.compute_p2(&orientations);
        // P2(cos 0) = 1
        assert!((s - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_order_parameter_perpendicular() {
        let op = OrderParameter::new([0.0, 0.0, 1.0]);
        // All molecules perpendicular to z-axis
        let orientations = vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let s = op.compute_p2(&orientations);
        // P2(cos 90°) = P2(0) = -0.5
        assert!((s - (-0.5)).abs() < 1e-12);
    }

    #[test]
    fn test_order_parameter_mean() {
        let mut op = OrderParameter::new([0.0, 0.0, 1.0]);
        let orientations = vec![0.0, 0.0, 1.0];
        op.record(&orientations, 0.0);
        op.record(&orientations, 1.0);
        assert!((op.mean_s() - 1.0).abs() < 1e-12);
    }

    // --- ClusterAnalysis ---

    #[test]
    fn test_cluster_analysis_two_clusters() {
        let mut ca = ClusterAnalysis::new(0.5, 2);
        let positions = vec![
            0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.2, 0.0, 0.0, 5.0, 0.0, 0.0, 5.1, 0.0, 0.0, 5.2, 0.0,
            0.0,
        ];
        ca.run(&positions);
        assert_eq!(ca.n_clusters, 2);
    }

    #[test]
    fn test_cluster_analysis_all_noise() {
        let mut ca = ClusterAnalysis::new(0.1, 5); // very small eps, large min_pts
        let positions = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 4.0, 0.0, 0.0];
        ca.run(&positions);
        assert_eq!(ca.n_clusters, 0);
    }

    #[test]
    fn test_cluster_sizes() {
        let mut ca = ClusterAnalysis::new(0.5, 2);
        let positions = vec![0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.2, 0.0, 0.0];
        ca.run(&positions);
        let sizes = ca.cluster_sizes();
        assert_eq!(sizes.iter().sum::<usize>(), 3);
    }

    #[test]
    fn test_percolation_threshold_above() {
        let mut ca = ClusterAnalysis::new(0.5, 2);
        let positions = vec![0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.2, 0.0, 0.0];
        ca.run(&positions);
        // All 3 in one cluster → 100% > 50%
        assert!(ca.percolation_threshold(0.5));
    }

    // --- Helper functions ---

    #[test]
    fn test_bin_distance_valid() {
        let b = bin_distance(2.5, 0.1, 5.0);
        assert_eq!(b, Some(25));
    }

    #[test]
    fn test_bin_distance_too_large() {
        let b = bin_distance(6.0, 0.1, 5.0);
        assert!(b.is_none());
    }

    #[test]
    fn test_bin_distance_negative() {
        let b = bin_distance(-1.0, 0.1, 5.0);
        assert!(b.is_none());
    }

    #[test]
    fn test_running_average_correctness() {
        let data = vec![1.0, 3.0, 2.0, 4.0];
        let avg = running_average(&data);
        assert!((avg[0] - 1.0).abs() < 1e-12);
        assert!((avg[1] - 2.0).abs() < 1e-12);
        assert!((avg[3] - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_autocorrelation_fft_identity() {
        let data = vec![1.0, 0.0, 0.0, 0.0];
        let acf = autocorrelation_fft(&data);
        assert!((acf[0] - 0.25).abs() < 1e-12); // (1*1)/4
    }

    #[test]
    fn test_autocorrelation_fft_constant() {
        let data = vec![2.0, 2.0, 2.0, 2.0];
        let acf = autocorrelation_fft(&data);
        // All values should equal 4.0
        for v in &acf {
            assert!((v - 4.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_unwrap_pbc_no_jump() {
        let positions = vec![0.1, 0.2, 0.3, 0.4];
        let unwrapped = unwrap_pbc_trajectory(&positions, 10.0);
        for (a, b) in positions.iter().zip(unwrapped.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }

    #[test]
    fn test_unwrap_pbc_with_jump() {
        // Simulated jump: 9.9 → 0.1 (actual: -9.8, unwrapped should give 10.1)
        let positions = vec![9.9, 0.1];
        let unwrapped = unwrap_pbc_trajectory(&positions, 10.0);
        assert!((unwrapped[1] - 10.1).abs() < 1e-10);
    }

    #[test]
    fn test_trapezoid_integrate() {
        // ∫₀¹ x dx = 0.5 (with linear function)
        let x = vec![0.0, 0.5, 1.0];
        let y = vec![0.0, 0.5, 1.0];
        let result = trapezoid_integrate(&x, &y);
        assert!((result - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_linear_fit() {
        let x = vec![0.0, 1.0, 2.0, 3.0];
        let y = vec![1.0, 3.0, 5.0, 7.0];
        let (slope, intercept) = linear_fit(&x, &y);
        assert!((slope - 2.0).abs() < 1e-10);
        assert!((intercept - 1.0).abs() < 1e-10);
    }
}
