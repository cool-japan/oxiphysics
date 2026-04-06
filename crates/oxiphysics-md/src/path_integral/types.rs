//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{HBAR, KB};

/// Propagator for the free ring polymer (spring part only) in normal-mode space.
///
/// Each normal mode k (k > 0) oscillates at frequency 2*omega_P*sin(k*pi/P).
/// The centroid (k=0) undergoes free translation.
pub struct BeadPropagator;
impl BeadPropagator {
    /// Normal mode frequencies for a ring polymer of P beads at temperature T.
    ///
    /// omega_k = 2 * omega_P * sin(k * pi / P)   for k = 0, 1, ..., P-1
    pub fn normal_mode_frequencies(n_beads: usize, _mass: f64, temperature: f64) -> Vec<f64> {
        let p = n_beads as f64;
        let omega_p = p * KB * temperature / HBAR;
        (0..n_beads)
            .map(|k| {
                let kf = k as f64;
                2.0 * omega_p * (kf * std::f64::consts::PI / p).sin()
            })
            .collect()
    }
    /// Propagate a single normal mode by dt using exact (harmonic) evolution.
    ///
    /// For mode k > 0 with frequency omega_k:
    ///   q(t+dt) = q(t)*cos(omega_k*dt) + (v(t)/omega_k)*sin(omega_k*dt)
    ///   v(t+dt) = -q(t)*omega_k*sin(omega_k*dt) + v(t)*cos(omega_k*dt)
    ///
    /// For centroid (k=0): free particle motion.
    pub fn propagate_mode(q: &mut [f64; 3], v: &mut [f64; 3], omega_k: f64, dt: f64) {
        if omega_k.abs() < 1e-30 {
            for d in 0..3 {
                q[d] += v[d] * dt;
            }
        } else {
            let cos_wdt = (omega_k * dt).cos();
            let sin_wdt = (omega_k * dt).sin();
            for d in 0..3 {
                let q_old = q[d];
                let v_old = v[d];
                q[d] = q_old * cos_wdt + (v_old / omega_k) * sin_wdt;
                v[d] = -q_old * omega_k * sin_wdt + v_old * cos_wdt;
            }
        }
    }
}
/// Multi-atom PIMD system: a collection of `PimdAtom` ring polymers.
#[derive(Debug, Clone)]
pub struct PimdState {
    /// Atoms in the system (each a ring polymer).
    pub atoms: Vec<PimdAtom>,
    /// System temperature (K).
    pub temperature: f64,
    /// Number of Trotter beads per atom.
    pub n_beads: usize,
}
impl PimdState {
    /// Quantum kinetic energy via virial estimator.
    ///
    /// KE_vir = (3*N/2)*k_B*T + (1/(2*P)) * sum_{a,i} (r_{a,i} - r_bar_a) . F_{a,i}
    pub fn quantum_kinetic_energy(&self) -> f64 {
        let n = self.atoms.len();
        let p = self.n_beads as f64;
        let classical_part = 1.5 * n as f64 * KB * self.temperature;
        let mut virial = 0.0f64;
        for atom in &self.atoms {
            let c = atom.centroid();
            for i in 0..self.n_beads {
                for d in 0..3 {
                    virial += (atom.replicas[i][d] - c[d]) * atom.forces[i][d];
                }
            }
        }
        classical_part + virial / (2.0 * p)
    }
    /// Centroid positions of all atoms.
    pub fn centroid_positions(&self) -> Vec<[f64; 3]> {
        self.atoms.iter().map(|a| a.centroid()).collect()
    }
}
/// Thermostatted Ring Polymer MD integrator.
///
/// Combines the exact free ring-polymer propagator (in normal-mode space)
/// with PILE-L friction for the stochastic/thermostat part and a velocity
/// Verlet step for the physical forces.
#[derive(Debug, Clone)]
pub struct TrpmdIntegrator {
    /// Number of ring-polymer beads.
    pub n_beads: usize,
    /// Simulation temperature (K).
    pub temperature: f64,
    /// Normal-mode propagator.
    pub propagator: NormalModePropagator,
    /// PILE-L thermostat.
    pub thermostat: PileLThermostat,
}
impl TrpmdIntegrator {
    /// Create a TRPMD integrator.
    pub fn new(n_beads: usize, temperature: f64, gamma0: f64) -> Self {
        TrpmdIntegrator {
            n_beads,
            temperature,
            propagator: NormalModePropagator::new(n_beads, temperature),
            thermostat: PileLThermostat::new(n_beads, temperature, gamma0),
        }
    }
    /// One TRPMD step for a single ring polymer.
    ///
    /// Algorithm (BAOAB-like):
    /// 1. Half thermostat (O) in normal-mode space
    /// 2. Half force kick (B) in bead space
    /// 3. Free propagation (A) in normal-mode space
    /// 4. Half force kick (B) in bead space
    /// 5. Half thermostat (O) in normal-mode space
    ///
    /// `forces` — external forces on each bead (length n_beads).
    pub fn step(
        &self,
        pos: &mut Vec<[f64; 3]>,
        vel: &mut Vec<[f64; 3]>,
        forces: &[[f64; 3]],
        mass: f64,
        dt: f64,
    ) {
        let c = &self.propagator.c_matrix;
        let mut mode_vel = beads_to_normal_modes(vel, c);
        self.thermostat.apply(&mut mode_vel, mass, dt * 0.5);
        *vel = normal_modes_to_beads(&mode_vel, c);
        for i in 0..self.n_beads {
            for d in 0..3 {
                vel[i][d] += 0.5 * dt * forces[i][d] / mass;
            }
        }
        let (new_pos, new_vel) = self.propagator.propagate(pos, vel, mass, dt);
        *pos = new_pos;
        *vel = new_vel;
        for i in 0..self.n_beads {
            for d in 0..3 {
                vel[i][d] += 0.5 * dt * forces[i][d] / mass;
            }
        }
        let mut mode_vel2 = beads_to_normal_modes(vel, c);
        self.thermostat.apply(&mut mode_vel2, mass, dt * 0.5);
        *vel = normal_modes_to_beads(&mode_vel2, c);
    }
}
/// A single atom represented as a ring polymer for PIMD.
///
/// Holds `n_beads` replica positions, fictitious masses (all equal to `mass`),
/// and the current (physical) forces on each replica.
#[derive(Debug, Clone)]
pub struct PimdAtom {
    /// Replica bead positions.
    pub replicas: Vec<[f64; 3]>,
    /// Per-replica masses (all equal to `mass`).
    pub masses: Vec<f64>,
    /// Per-replica physical forces (updated externally).
    pub forces: Vec<[f64; 3]>,
    /// Number of Trotter beads.
    pub n_beads: usize,
}
impl PimdAtom {
    /// Create a new `PimdAtom` with all replicas at the origin.
    pub fn new(n_beads: usize, mass: f64) -> Self {
        Self {
            replicas: vec![[0.0; 3]; n_beads],
            masses: vec![mass; n_beads],
            forces: vec![[0.0; 3]; n_beads],
            n_beads,
        }
    }
    /// Centroid position: average over all replicas.
    pub fn centroid(&self) -> [f64; 3] {
        let mut c = [0.0f64; 3];
        for r in &self.replicas {
            for d in 0..3 {
                c[d] += r[d];
            }
        }
        let inv_p = 1.0 / self.n_beads as f64;
        [c[0] * inv_p, c[1] * inv_p, c[2] * inv_p]
    }
    /// Harmonic spring forces between adjacent replicas at temperature `T`.
    ///
    /// The spring constant is k = m * omega_P^2 where omega_P = P*k_B*T/hbar.
    /// F_i = m * omega_P^2 * (r_{i-1} + r_{i+1} - 2*r_i)
    pub fn spring_forces(&self, t: f64) -> Vec<[f64; 3]> {
        let p = self.n_beads;
        let mass = self.masses[0];
        let omega_p = p as f64 * KB * t / HBAR;
        let k = mass * omega_p * omega_p;
        let mut fs = vec![[0.0f64; 3]; p];
        for i in 0..p {
            let prev = (i + p - 1) % p;
            let next = (i + 1) % p;
            for d in 0..3 {
                fs[i][d] = k
                    * (self.replicas[prev][d] + self.replicas[next][d] - 2.0 * self.replicas[i][d]);
            }
        }
        fs
    }
}
/// PILE-L (Path Integral Langevin Equation - Least-coupling) thermostat.
///
/// Each normal mode k is coupled to a Langevin bath with friction γ_k.
/// The centroid mode uses the user-supplied γ₀, while higher modes use
/// γ_k = 2 ω_k (critical damping).
///
/// Reference: Ceriotti, Parrinello, Markland & Manolopoulos (2010).
#[derive(Debug, Clone)]
pub struct PileLThermostat {
    /// Number of beads.
    pub n_beads: usize,
    /// Centroid friction γ₀ (rad s^-1 or ps^-1 in MD units).
    pub gamma0: f64,
    /// Per-mode friction γ_k (length n_beads).
    pub gamma: Vec<f64>,
    /// Temperature (K).
    pub temperature: f64,
    /// Normal-mode frequencies.
    pub omega: Vec<f64>,
}
impl PileLThermostat {
    /// Create a PILE-L thermostat.
    ///
    /// * `n_beads` — number of ring-polymer beads P
    /// * `temperature` — simulation temperature (K)
    /// * `gamma0` — centroid friction (ps^-1 in typical MD units)
    pub fn new(n_beads: usize, temperature: f64, gamma0: f64) -> Self {
        let omega = normal_mode_frequencies(n_beads, temperature);
        let mut gamma = vec![0.0f64; n_beads];
        gamma[0] = gamma0;
        for k in 1..n_beads {
            gamma[k] = 2.0 * omega[k];
        }
        PileLThermostat {
            n_beads,
            gamma0,
            gamma,
            temperature,
            omega,
        }
    }
    /// Apply the O-U step to a single normal-mode velocity component.
    ///
    /// v → v exp(-γ_k dt/2) + σ_k ξ, where σ_k = sqrt(k_B T / m) * sqrt(1 - exp(-γ_k dt)).
    fn ou_step(gamma_k: f64, temperature: f64, mass: f64, v: f64, dt: f64) -> f64 {
        let decay = (-gamma_k * dt).exp();
        let sigma = ((1.0 - decay * decay) * KB * temperature / mass).sqrt();
        let mut rng = rand::rng();
        use rand::RngExt;
        let xi: f64 = rng.random_range(-1.0_f64..1.0_f64);
        decay * v + sigma * xi
    }
    /// Apply PILE-L thermostat to normal-mode velocities of a single particle.
    ///
    /// Performs the stochastic velocity update for each mode k.
    pub fn apply(&self, mode_vel: &mut [[f64; 3]], mass: f64, dt: f64) {
        for k in 0..self.n_beads {
            for d in 0..3 {
                mode_vel[k][d] =
                    Self::ou_step(self.gamma[k], self.temperature, mass, mode_vel[k][d], dt);
            }
        }
    }
}
/// Centroid Molecular Dynamics (CMD) propagator.
///
/// In CMD the physical force acts only on the centroid degree of freedom,
/// while the internal modes are thermostatted to sample the quantum
/// distribution.  The centroid trajectory approximates real-time quantum
/// dynamics for computing time-correlation functions.
///
/// Reference: Cao & Voth, J. Chem. Phys. 99, 10070 (1993).
pub struct CentroidMD;
impl CentroidMD {
    /// Propagate the centroid position and velocity by dt under the centroid force.
    ///
    /// The centroid force is the mean of the per-bead external forces:
    /// F_c = (1/P) * sum_i F_ext_i
    pub fn step(ring: &mut RingPolymer, per_bead_forces: &[[f64; 3]], dt: f64) {
        let p = ring.n_beads as f64;
        let mut f_c = [0.0_f64; 3];
        for i in 0..ring.n_beads {
            for d in 0..3 {
                f_c[d] += per_bead_forces[i][d];
            }
        }
        for d in 0..3 {
            f_c[d] /= p;
        }
        let inv_m = 1.0 / ring.mass;
        let mut cv = ring.centroid_velocity();
        for d in 0..3 {
            cv[d] += f_c[d] * inv_m * 0.5 * dt;
        }
        let mut cpos = ring.centroid();
        for d in 0..3 {
            cpos[d] += cv[d] * dt;
        }
        for d in 0..3 {
            cv[d] += f_c[d] * inv_m * 0.5 * dt;
        }
        let old_c = ring.centroid();
        for i in 0..ring.n_beads {
            for d in 0..3 {
                ring.positions[i][d] += cpos[d] - old_c[d];
            }
            ring.velocities[i].copy_from_slice(&cv);
        }
    }
    /// Centroid mean-square displacement from initial centroid position.
    pub fn centroid_msd(ring: &RingPolymer, initial_centroid: [f64; 3]) -> f64 {
        let c = ring.centroid();
        (0..3).map(|d| (c[d] - initial_centroid[d]).powi(2)).sum()
    }
}
/// Ring Polymer Molecular Dynamics (RPMD) propagator.
///
/// RPMD uses the ring polymer Hamiltonian to compute approximate
/// quantum dynamical properties. The ring polymer evolves under both
/// the inter-bead springs and the physical potential.
///
/// Reference: Craig & Manolopoulos, J. Chem. Phys. 121, 3368 (2004).
pub struct RpmdStep;
impl RpmdStep {
    /// Single RPMD step: identical to velocity-Verlet but with per-bead
    /// physical forces.
    ///
    /// In RPMD the internal (spring) modes are NOT thermostatted -- only
    /// the physical forces and springs drive the dynamics.
    pub fn step(ring: &mut RingPolymer, physical_forces: &[[f64; 3]], temperature: f64, dt: f64) {
        PimdStep::velocity_verlet_per_bead(ring, physical_forces, temperature, dt);
    }
    /// Multiple RPMD steps with constant forces (useful for free-ring propagation).
    pub fn multi_step(
        ring: &mut RingPolymer,
        physical_forces: &[[f64; 3]],
        temperature: f64,
        dt: f64,
        n_steps: usize,
    ) {
        for _ in 0..n_steps {
            Self::step(ring, physical_forces, temperature, dt);
        }
    }
    /// Compute the RPMD Hamiltonian (for energy conservation checks).
    ///
    /// H = KE_classical + V_spring + (1/P) * sum V_phys(r_i)
    pub fn hamiltonian(ring: &RingPolymer, temperature: f64, bead_potentials: &[f64]) -> f64 {
        let ke = ring.classical_kinetic_energy();
        let v_spring = ring.spring_potential_energy(temperature);
        let v_phys = ring.potential_energy_estimator(bead_potentials);
        ke + v_spring + v_phys
    }
}
/// Driver for a single PIMD velocity-Verlet step.
pub struct PimdStep;
impl PimdStep {
    /// Perform one velocity-Verlet step for the ring polymer.
    pub fn velocity_verlet(
        ring: &mut RingPolymer,
        external_force: [f64; 3],
        temperature: f64,
        dt: f64,
    ) {
        let n = ring.n_beads;
        let inv_m = 1.0 / ring.bead_mass;
        let f_spring = ring.spring_forces(temperature);
        for i in 0..n {
            for d in 0..3 {
                let f_total = f_spring[i][d] + external_force[d];
                ring.velocities[i][d] += f_total * inv_m * 0.5 * dt;
            }
        }
        for i in 0..n {
            for d in 0..3 {
                ring.positions[i][d] += ring.velocities[i][d] * dt;
            }
        }
        let f_spring_new = ring.spring_forces(temperature);
        for i in 0..n {
            for d in 0..3 {
                let f_total = f_spring_new[i][d] + external_force[d];
                ring.velocities[i][d] += f_total * inv_m * 0.5 * dt;
            }
        }
    }
    /// Perform one velocity-Verlet step with per-bead external forces.
    pub fn velocity_verlet_per_bead(
        ring: &mut RingPolymer,
        external_forces: &[[f64; 3]],
        temperature: f64,
        dt: f64,
    ) {
        assert_eq!(external_forces.len(), ring.n_beads);
        let n = ring.n_beads;
        let inv_m = 1.0 / ring.bead_mass;
        let f_spring = ring.spring_forces(temperature);
        for i in 0..n {
            for d in 0..3 {
                let f_total = f_spring[i][d] + external_forces[i][d];
                ring.velocities[i][d] += f_total * inv_m * 0.5 * dt;
            }
        }
        for i in 0..n {
            for d in 0..3 {
                ring.positions[i][d] += ring.velocities[i][d] * dt;
            }
        }
        let f_spring_new = ring.spring_forces(temperature);
        for i in 0..n {
            for d in 0..3 {
                let f_total = f_spring_new[i][d] + external_forces[i][d];
                ring.velocities[i][d] += f_total * inv_m * 0.5 * dt;
            }
        }
    }
}
/// Path Integral Langevin Equation (PILE) thermostat.
///
/// Applies Langevin friction and noise to each bead independently.
/// The centroid mode can use a different friction coefficient than
/// the internal modes for optimal sampling.
///
/// Reference: Ceriotti, Parrinello, Markland & Manolopoulos,
/// J. Chem. Phys. 133, 124104 (2010).
pub struct PileThermostat {
    /// Friction coefficient for centroid mode (1/s).
    pub gamma_centroid: f64,
    /// Friction coefficient for internal modes (1/s).
    /// If None, uses the optimal value omega_P.
    pub gamma_internal: Option<f64>,
    /// Temperature (K).
    pub temperature: f64,
}
impl PileThermostat {
    /// Create a PILE thermostat.
    pub fn new(gamma_centroid: f64, temperature: f64) -> Self {
        Self {
            gamma_centroid,
            gamma_internal: None,
            temperature,
        }
    }
    /// Create a PILE thermostat with explicit internal friction.
    pub fn with_internal_friction(
        gamma_centroid: f64,
        gamma_internal: f64,
        temperature: f64,
    ) -> Self {
        Self {
            gamma_centroid,
            gamma_internal: Some(gamma_internal),
            temperature,
        }
    }
    /// Optimal internal-mode friction coefficient: omega_P = P * k_B * T / hbar.
    pub fn optimal_gamma_internal(&self, n_beads: usize) -> f64 {
        let p = n_beads as f64;
        p * KB * self.temperature / HBAR
    }
    /// Apply the Langevin thermostat to the centroid velocity.
    ///
    /// Uses the velocity-rescale form:
    /// v_new = c1 * v_old + c2 * R
    /// where c1 = exp(-gamma*dt), c2 = sqrt((1-c1^2)*k_B*T/m), R ~ N(0,1).
    ///
    /// `noise` should be a sample from the standard normal distribution N(0,1).
    pub fn apply_centroid(&self, centroid_vel: &mut [f64; 3], mass: f64, dt: f64, noise: [f64; 3]) {
        let c1 = (-self.gamma_centroid * dt).exp();
        let sigma = (KB * self.temperature / mass).sqrt();
        let c2 = ((1.0 - c1 * c1).max(0.0)).sqrt() * sigma;
        for d in 0..3 {
            centroid_vel[d] = c1 * centroid_vel[d] + c2 * noise[d];
        }
    }
    /// Apply the thermostat to an internal (non-centroid) bead velocity.
    pub fn apply_internal(
        &self,
        vel: &mut [f64; 3],
        mass: f64,
        n_beads: usize,
        dt: f64,
        noise: [f64; 3],
    ) {
        let gamma = self
            .gamma_internal
            .unwrap_or_else(|| self.optimal_gamma_internal(n_beads));
        let c1 = (-gamma * dt).exp();
        let sigma = (KB * self.temperature / mass).sqrt();
        let c2 = ((1.0 - c1 * c1).max(0.0)).sqrt() * sigma;
        for d in 0..3 {
            vel[d] = c1 * vel[d] + c2 * noise[d];
        }
    }
}
/// Classical ring polymer representing a single quantum particle.
///
/// A particle of mass `mass` is replicated into `n_beads` beads arranged
/// cyclically.  The inter-bead harmonic springs reproduce quantum
/// fluctuations at temperature `T`.
#[derive(Debug, Clone)]
pub struct RingPolymer {
    /// Number of Trotter beads P.
    pub n_beads: usize,
    /// Physical particle mass (kg).
    pub mass: f64,
    /// Bead positions (m), indexed `[bead][x/y/z]`.
    pub positions: Vec<[f64; 3]>,
    /// Bead velocities (m s^-1).
    pub velocities: Vec<[f64; 3]>,
    /// Effective bead mass for centroid-mode dynamics.
    pub bead_mass: f64,
}
impl RingPolymer {
    /// Create a new ring polymer with all beads placed at `center`.
    pub fn new(n_beads: usize, mass: f64, center: [f64; 3]) -> Self {
        assert!(n_beads >= 1, "ring polymer needs at least 1 bead");
        let positions = vec![center; n_beads];
        let velocities = vec![[0.0_f64; 3]; n_beads];
        Self {
            n_beads,
            mass,
            positions,
            velocities,
            bead_mass: mass,
        }
    }
    /// Harmonic spring forces between adjacent beads.
    ///
    /// omega_P = P * k_B * T / hbar
    /// F_i = -m * omega_P^2 * (2*r_i - r_{i-1} - r_{i+1})
    pub fn spring_forces(&self, temperature: f64) -> Vec<[f64; 3]> {
        let p = self.n_beads as f64;
        let omega_p = p * KB * temperature / HBAR;
        let k = self.mass * omega_p * omega_p;
        let mut forces = vec![[0.0_f64; 3]; self.n_beads];
        for i in 0..self.n_beads {
            let prev = (i + self.n_beads - 1) % self.n_beads;
            let next = (i + 1) % self.n_beads;
            for d in 0..3 {
                forces[i][d] = -k
                    * (2.0 * self.positions[i][d]
                        - self.positions[prev][d]
                        - self.positions[next][d]);
            }
        }
        forces
    }
    /// Spring potential energy: sum of (1/2)*m*omega_P^2 * |r_{i+1} - r_i|^2.
    pub fn spring_potential_energy(&self, temperature: f64) -> f64 {
        let p = self.n_beads as f64;
        let omega_p = p * KB * temperature / HBAR;
        let k = self.mass * omega_p * omega_p;
        let sum: f64 = (0..self.n_beads)
            .map(|i| {
                let j = (i + 1) % self.n_beads;
                (0..3)
                    .map(|d| (self.positions[j][d] - self.positions[i][d]).powi(2))
                    .sum::<f64>()
            })
            .sum();
        0.5 * k * sum
    }
    /// Centroid position -- mean over all beads.
    pub fn centroid(&self) -> [f64; 3] {
        let mut c = [0.0_f64; 3];
        for pos in &self.positions {
            for d in 0..3 {
                c[d] += pos[d];
            }
        }
        let inv_p = 1.0 / (self.n_beads as f64);
        [c[0] * inv_p, c[1] * inv_p, c[2] * inv_p]
    }
    /// Centroid velocity -- mean velocity over all beads.
    pub fn centroid_velocity(&self) -> [f64; 3] {
        let mut cv = [0.0_f64; 3];
        for vel in &self.velocities {
            for d in 0..3 {
                cv[d] += vel[d];
            }
        }
        let inv_p = 1.0 / (self.n_beads as f64);
        [cv[0] * inv_p, cv[1] * inv_p, cv[2] * inv_p]
    }
    /// Radius of gyration of the ring polymer (quantum spread).
    ///
    /// R_g = sqrt( (1/P) sum_i |r_i - r_bar|^2 )
    pub fn radius_of_gyration(&self) -> f64 {
        let c = self.centroid();
        let msd: f64 = self
            .positions
            .iter()
            .map(|r| (0..3).map(|d| (r[d] - c[d]).powi(2)).sum::<f64>())
            .sum::<f64>();
        (msd / (self.n_beads as f64)).sqrt()
    }
    /// Convert Cartesian bead coordinates to staging coordinates.
    ///
    /// u_0 = centroid = (1/P) sum r_k
    /// u_k = r_k - r_{k-1}  for k = 1 ... P-1
    pub fn to_staging(&self) -> Vec<[f64; 3]> {
        let mut staging = vec![[0.0_f64; 3]; self.n_beads];
        staging[0] = self.centroid();
        for k in 1..self.n_beads {
            for d in 0..3 {
                staging[k][d] = self.positions[k][d] - self.positions[k - 1][d];
            }
        }
        staging
    }
    /// Convert staging coordinates back to Cartesian bead positions.
    pub fn from_staging(staging: &[[f64; 3]], n_beads: usize) -> Vec<[f64; 3]> {
        assert_eq!(staging.len(), n_beads, "staging length mismatch");
        if n_beads == 0 {
            return vec![];
        }
        if n_beads == 1 {
            return vec![staging[0]];
        }
        let p = n_beads as f64;
        let mut partial_sums = vec![[0.0_f64; 3]; n_beads];
        for k in 1..n_beads {
            for d in 0..3 {
                partial_sums[k][d] = partial_sums[k - 1][d] + staging[k][d];
            }
        }
        let mut r0 = [0.0_f64; 3];
        for d in 0..3 {
            let sum_sk: f64 = (1..n_beads).map(|k| partial_sums[k][d]).sum();
            r0[d] = staging[0][d] - sum_sk / p;
        }
        let mut positions = vec![[0.0_f64; 3]; n_beads];
        positions[0] = r0;
        for k in 1..n_beads {
            for d in 0..3 {
                positions[k][d] = r0[d] + partial_sums[k][d];
            }
        }
        positions
    }
    /// Convert bead coordinates to normal mode coordinates via matrix transform.
    ///
    /// Uses an explicit P x P orthogonal transform matrix built from
    /// the DFT basis. The transform stores the forward DFT coefficients
    /// (complex, stored as pairs) and the inverse recovers positions exactly.
    ///
    /// For P beads, modes\[k\]\[d\] stores the k-th mode amplitude for dimension d.
    /// The transform matrix T\[k\]\[n\] is computed on the fly.
    pub fn to_normal_modes(&self) -> Vec<[f64; 3]> {
        let p = self.n_beads;
        let t_matrix = build_normal_mode_matrix(p);
        let mut modes = vec![[0.0_f64; 3]; p];
        for k in 0..p {
            for d in 0..3 {
                let mut sum = 0.0;
                for n in 0..p {
                    sum += t_matrix[k * p + n] * self.positions[n][d];
                }
                modes[k][d] = sum;
            }
        }
        modes
    }
    /// Convert normal mode coordinates back to Cartesian bead positions.
    pub fn from_normal_modes(modes: &[[f64; 3]], n_beads: usize) -> Vec<[f64; 3]> {
        let t_matrix = build_normal_mode_matrix(n_beads);
        let mut positions = vec![[0.0_f64; 3]; n_beads];
        for n in 0..n_beads {
            for d in 0..3 {
                let mut sum = 0.0;
                for k in 0..n_beads {
                    sum += t_matrix[k * n_beads + n] * modes[k][d];
                }
                positions[n][d] = sum;
            }
        }
        positions
    }
    /// Primitive estimator for the quantum kinetic energy.
    ///
    /// KE_prim = 3*P*k_B*T/2 - (m*omega_P^2/2) * sum_i |r_{i+1} - r_i|^2
    pub fn primitive_kinetic_energy(&self, temperature: f64) -> f64 {
        let p = self.n_beads as f64;
        let omega_p = p * KB * temperature / HBAR;
        let classical_ke = 1.5 * p * KB * temperature;
        let spring_sum: f64 = (0..self.n_beads)
            .map(|i| {
                let j = (i + 1) % self.n_beads;
                (0..3)
                    .map(|d| (self.positions[j][d] - self.positions[i][d]).powi(2))
                    .sum::<f64>()
            })
            .sum();
        classical_ke - 0.5 * self.mass * omega_p * omega_p * spring_sum
    }
    /// Virial kinetic energy estimator.
    ///
    /// KE_vir = (3/2)*k_B*T + (1/(2P)) * sum_i (r_i - r_bar) . F_ext_i
    ///
    /// where F_ext_i is the external (physical) force on bead i.
    /// The virial estimator has lower variance than the primitive estimator.
    pub fn virial_kinetic_energy(&self, temperature: f64, external_forces: &[[f64; 3]]) -> f64 {
        let c = self.centroid();
        let p = self.n_beads as f64;
        let mut virial_sum = 0.0_f64;
        for i in 0..self.n_beads {
            for d in 0..3 {
                virial_sum += (self.positions[i][d] - c[d]) * external_forces[i][d];
            }
        }
        1.5 * KB * temperature + virial_sum / (2.0 * p)
    }
    /// Thermodynamic estimator for the potential energy.
    ///
    /// V_est = (1/P) * sum_i V(r_i)
    ///
    /// where V(r_i) is the physical potential evaluated at each bead position.
    pub fn potential_energy_estimator(&self, bead_potentials: &[f64]) -> f64 {
        assert_eq!(bead_potentials.len(), self.n_beads);
        bead_potentials.iter().sum::<f64>() / (self.n_beads as f64)
    }
    /// Classical kinetic energy of the beads (for thermostat monitoring).
    pub fn classical_kinetic_energy(&self) -> f64 {
        let mut ke = 0.0_f64;
        for vel in &self.velocities {
            for d in 0..3 {
                ke += 0.5 * self.bead_mass * vel[d] * vel[d];
            }
        }
        ke
    }
    /// Instantaneous temperature from classical kinetic energy.
    ///
    /// T = 2*KE / (3*P*k_B)
    pub fn instantaneous_temperature(&self) -> f64 {
        let ke = self.classical_kinetic_energy();
        2.0 * ke / (3.0 * self.n_beads as f64 * KB)
    }
}
/// Nosé-Hoover-style PIMD thermostat for internal (non-centroid) modes.
///
/// Each internal mode k has its own thermostat with coupling frequency omega_k.
/// The centroid mode uses a separate, softer thermostat.
pub struct PimdThermostat {
    /// Friction coefficients for each mode.
    pub gamma: Vec<f64>,
    /// Temperature (K).
    pub temperature: f64,
    /// Mass of the physical particle (kg).
    pub mass: f64,
}
impl PimdThermostat {
    /// Create a PIMD thermostat with optimal mode-dependent friction.
    ///
    /// gamma_k = 2 * omega_k  (critically damped)
    pub fn new_optimal(n_beads: usize, mass: f64, temperature: f64) -> Self {
        let freqs = StagingTransform::staging_frequencies(n_beads, temperature);
        let gamma: Vec<f64> = freqs.iter().map(|&w| 2.0 * w).collect();
        Self {
            gamma,
            temperature,
            mass,
        }
    }
    /// Apply Langevin thermostat to a single normal mode velocity.
    ///
    /// v_new = c1 * v + c2 * noise
    /// where c1 = exp(-gamma_k * dt), c2 = sqrt((1-c1^2) * k_B*T/m)
    pub fn apply_mode(&self, vel: &mut [f64; 3], mode_index: usize, dt: f64, noise: [f64; 3]) {
        if mode_index >= self.gamma.len() {
            return;
        }
        let g = self.gamma[mode_index];
        let c1 = (-g * dt).exp();
        let sigma = (KB * self.temperature / self.mass).sqrt();
        let c2 = ((1.0 - c1 * c1).max(0.0)).sqrt() * sigma;
        for d in 0..3 {
            vel[d] = c1 * vel[d] + c2 * noise[d];
        }
    }
}
/// Free-particle ring-polymer propagator in normal-mode space.
///
/// For each non-centroid mode k with frequency ω_k, the exact harmonic
/// propagator over time step dt is applied to (q_k, p_k/m):
///
///   q_k(t+dt) = q_k cos(ω_k dt) + p_k/(m ω_k) sin(ω_k dt)
///   p_k(t+dt) = -m ω_k q_k sin(ω_k dt) + p_k cos(ω_k dt)
///
/// The centroid mode (k=0) is propagated with free-particle motion.
#[derive(Debug, Clone)]
pub struct NormalModePropagator {
    /// Number of beads.
    pub n_beads: usize,
    /// Normal-mode frequencies (rad s^-1).
    pub omega: Vec<f64>,
    /// Transformation matrix C (n_beads × n_beads).
    pub c_matrix: Vec<Vec<f64>>,
}
impl NormalModePropagator {
    /// Build a normal-mode propagator for temperature `t` and `n_beads` beads.
    pub fn new(n_beads: usize, temperature: f64) -> Self {
        NormalModePropagator {
            n_beads,
            omega: normal_mode_frequencies(n_beads, temperature),
            c_matrix: normal_mode_matrix(n_beads),
        }
    }
    /// Apply the free ring-polymer propagator to a single particle.
    ///
    /// `pos` and `vel` are in bead representation (length P each).
    /// `mass` is the physical particle mass.
    /// Returns updated (pos, vel) in bead space.
    pub fn propagate(
        &self,
        pos: &[[f64; 3]],
        vel: &[[f64; 3]],
        mass: f64,
        dt: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let p = self.n_beads;
        let mut q = beads_to_normal_modes(pos, &self.c_matrix);
        let mut v = beads_to_normal_modes(vel, &self.c_matrix);
        for k in 0..p {
            let omega_k = self.omega[k];
            if omega_k < 1e-30 {
                for d in 0..3 {
                    q[k][d] += v[k][d] * dt;
                }
            } else {
                let cos_w = (omega_k * dt).cos();
                let sin_w = (omega_k * dt).sin();
                let inv_mw = 1.0 / (mass * omega_k);
                for d in 0..3 {
                    let q_new = q[k][d] * cos_w + v[k][d] * inv_mw * sin_w;
                    let v_new = -mass * omega_k * q[k][d] * sin_w + v[k][d] * cos_w;
                    q[k][d] = q_new;
                    v[k][d] = v_new;
                }
            }
        }
        let new_pos = normal_modes_to_beads(&q, &self.c_matrix);
        let new_vel = normal_modes_to_beads(&v, &self.c_matrix);
        (new_pos, new_vel)
    }
}
/// Ring Polymer Contraction: evaluate expensive potential on a contracted ring.
///
/// For a ring with P beads, the contracted ring of P' < P beads is obtained
/// by taking the P' lowest normal modes.  This allows the expensive part of
/// the potential (e.g. ab initio) to be evaluated on fewer replicas.
#[derive(Debug, Clone)]
pub struct RingPolymerContraction {
    /// Full bead number P.
    pub p_full: usize,
    /// Contracted bead number P'.
    pub p_contracted: usize,
}
impl RingPolymerContraction {
    /// Create a ring-polymer contraction scheme.
    pub fn new(p_full: usize, p_contracted: usize) -> Self {
        assert!(p_contracted <= p_full, "p_contracted must be <= p_full");
        RingPolymerContraction {
            p_full,
            p_contracted,
        }
    }
    /// Contract bead positions from P to P' by truncating normal modes.
    ///
    /// Returns P' bead positions representing the contracted ring.
    pub fn contract(&self, beads: &[[f64; 3]], temperature: f64) -> Vec<[f64; 3]> {
        let c_full = normal_mode_matrix(self.p_full);
        let modes_full = beads_to_normal_modes(beads, &c_full);
        let p = self.p_contracted;
        let c_cont = normal_mode_matrix(p);
        let omega_full = normal_mode_frequencies(self.p_full, temperature);
        let omega_cont = normal_mode_frequencies(p, temperature);
        let mut modes_cont = vec![[0.0f64; 3]; p];
        let centroid_scale = (p as f64 / self.p_full as f64).sqrt();
        for d in 0..3 {
            modes_cont[0][d] = modes_full[0][d] * centroid_scale;
        }
        for k in 1..p.min(self.p_full) {
            let scale = if omega_full[k] > 1e-30 {
                omega_cont[k] / omega_full[k]
            } else {
                1.0
            };
            for d in 0..3 {
                modes_cont[k][d] = modes_full[k][d] * scale;
            }
        }
        normal_modes_to_beads(&modes_cont, &c_cont)
    }
}
/// Quantum tunneling corrections using the Feynman-Hibbs effective potential.
///
/// The Feynman-Hibbs effective potential replaces the bare potential V(x)
/// by an effective potential that includes leading quantum corrections:
///
/// V_FH(x) = V(x) + (hbar^2 / (24 * m * k_B * T)) * d^2V/dx^2
///
/// This is the Wigner-Kirkwood expansion to order hbar^2.
///
/// Reference: Feynman & Hibbs, "Quantum Mechanics and Path Integrals" (1965).
pub struct FeynmanHibbsCorrection {
    /// Particle mass (kg).
    pub mass: f64,
    /// Temperature (K).
    pub temperature: f64,
}
impl FeynmanHibbsCorrection {
    /// Create a new Feynman-Hibbs correction calculator.
    pub fn new(mass: f64, temperature: f64) -> Self {
        Self { mass, temperature }
    }
    /// Feynman-Hibbs prefactor: hbar^2 / (24 * m * k_B * T).
    ///
    /// This multiplies the second derivative of the potential.
    pub fn prefactor(&self) -> f64 {
        let hbar = HBAR;
        hbar * hbar / (24.0 * self.mass * KB * self.temperature)
    }
    /// Feynman-Hibbs effective potential given V(x) and V''(x).
    ///
    /// V_eff = V + prefactor * V''
    pub fn effective_potential(&self, v: f64, v_second_deriv: f64) -> f64 {
        v + self.prefactor() * v_second_deriv
    }
    /// Quantum de Broglie thermal wavelength.
    ///
    /// lambda_th = hbar * sqrt(2*pi / (m * k_B * T))
    pub fn thermal_de_broglie_wavelength(&self) -> f64 {
        let hbar = HBAR;
        hbar * (2.0 * std::f64::consts::PI / (self.mass * KB * self.temperature)).sqrt()
    }
    /// Estimate quantum tunneling probability through a barrier of height E_b
    /// and width a using the semiclassical (WKB) approximation.
    ///
    /// T_WKB ≈ exp(-2*a/hbar * sqrt(2*m*E_b))
    pub fn wkb_tunneling_probability(&self, barrier_height: f64, barrier_width: f64) -> f64 {
        let hbar = HBAR;
        let exponent = -2.0 * barrier_width / hbar * (2.0 * self.mass * barrier_height).sqrt();
        exponent.exp().min(1.0)
    }
}
/// Proper staging transformation with fictitious bead masses.
///
/// The staging transformation decouples the inter-bead springs by introducing
/// fictitious masses m_k = (k/(k+1)) * m for k = 1..P-1, making each staging
/// mode evolve independently.
///
/// Reference: Tuckerman et al., J. Chem. Phys. 99, 2796 (1993).
pub struct StagingTransform;
impl StagingTransform {
    /// Compute the fictitious staging masses for a ring polymer of P beads.
    ///
    /// m_0 = m  (centroid mass, the full physical mass)
    /// m_k = (k+1)/k * m  for k = 1..P-1
    pub fn fictitious_masses(n_beads: usize, mass: f64) -> Vec<f64> {
        let mut masses = vec![mass; n_beads];
        for k in 1..n_beads {
            masses[k] = mass * (k as f64 + 1.0) / (k as f64);
        }
        masses
    }
    /// Convert Cartesian bead coordinates to staging coordinates.
    ///
    /// u_0 = x_0  (first bead)
    /// u_k = x_k - ((k-1)*x_{k+1} + x_0) / k   for k = 1..P-1
    ///
    /// This is the form that diagonalises the spring matrix exactly.
    pub fn to_staging_proper(positions: &[[f64; 3]], n_beads: usize) -> Vec<[f64; 3]> {
        if n_beads == 0 {
            return vec![];
        }
        if n_beads == 1 {
            return positions.to_vec();
        }
        let mut u = vec![[0.0_f64; 3]; n_beads];
        u[0] = positions[0];
        let p = n_beads;
        for k in 1..p {
            for d in 0..3 {
                let frac = k as f64 / (k as f64 + 1.0);
                let x_next = if k + 1 < p {
                    positions[k + 1][d]
                } else {
                    positions[0][d]
                };
                u[k][d] = positions[k][d] - frac * x_next - (1.0 - frac) * positions[0][d];
            }
        }
        u
    }
    /// Spring frequencies in staging coordinates.
    ///
    /// omega_k = sqrt(k*(k+1)) * omega_P  for k = 1..P-1
    /// omega_0 = 0  (centroid)
    pub fn staging_frequencies(n_beads: usize, temperature: f64) -> Vec<f64> {
        let p = n_beads as f64;
        let omega_p = p * KB * temperature / HBAR;
        let mut freqs = vec![0.0_f64; n_beads];
        for k in 1..n_beads {
            let kf = k as f64;
            freqs[k] = (kf * (kf + 1.0)).sqrt() * omega_p;
        }
        freqs
    }
}
