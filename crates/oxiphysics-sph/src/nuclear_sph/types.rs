//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{NUCLEON_MASS, RHO_0};

/// Nuclear density profile model.
///
/// Parametrized as a Woods-Saxon distribution:
/// `ρ(r) = ρ₀ / (1 + exp((r - R) / a))`
#[derive(Debug, Clone)]
pub struct WoodsSaxonProfile {
    /// Central density (nuclear saturation density ρ₀) in fm⁻³.
    pub rho_central: f64,
    /// Nuclear radius R = r₀ A^{1/3} in fm.
    pub radius_fm: f64,
    /// Diffuseness parameter a in fm.
    pub diffuseness_fm: f64,
    /// Mass number A.
    pub mass_number: u32,
    /// Atomic number Z.
    pub atomic_number: u32,
}
impl WoodsSaxonProfile {
    /// Construct a Woods-Saxon profile for a nucleus with mass number A and charge Z.
    pub fn new(mass_number: u32, atomic_number: u32) -> Self {
        let a = mass_number as f64;
        let r0 = 1.2;
        Self {
            rho_central: RHO_0,
            radius_fm: r0 * a.powf(1.0 / 3.0),
            diffuseness_fm: 0.524,
            mass_number,
            atomic_number,
        }
    }
    /// Evaluate nuclear density at radial distance r (fm).
    pub fn density(&self, r_fm: f64) -> f64 {
        self.rho_central / (1.0 + ((r_fm - self.radius_fm) / self.diffuseness_fm).exp())
    }
    /// Proton density ρₚ(r) = Z/A · ρ(r).
    pub fn proton_density(&self, r_fm: f64) -> f64 {
        self.density(r_fm) * (self.atomic_number as f64) / (self.mass_number as f64)
    }
    /// Neutron density ρₙ(r) = N/A · ρ(r).
    pub fn neutron_density(&self, r_fm: f64) -> f64 {
        let n = (self.mass_number - self.atomic_number) as f64;
        self.density(r_fm) * n / (self.mass_number as f64)
    }
    /// Root-mean-square radius (fm).
    pub fn rms_radius(&self) -> f64 {
        let n_steps = 1000usize;
        let r_max = self.radius_fm * 5.0;
        let dr = r_max / n_steps as f64;
        let mut num = 0.0_f64;
        let mut den = 0.0_f64;
        for i in 0..n_steps {
            let r = (i as f64 + 0.5) * dr;
            let rho = self.density(r);
            let dv = 4.0 * PI * r * r * dr;
            num += r * r * rho * dv;
            den += rho * dv;
        }
        if den < 1e-30 { 0.0 } else { (num / den).sqrt() }
    }
    /// Skin thickness (neutron skin) t = R_n - R_p in fm.
    pub fn neutron_skin_thickness(&self) -> f64 {
        let z = self.atomic_number as f64;
        let a = self.mass_number as f64;
        let t_approx = 0.9 * (a - 2.0 * z) / a;
        t_approx.max(0.0)
    }
}
/// Parametrization set for the Skyrme nuclear EOS.
///
/// Controls compressibility (stiff vs soft) and symmetry energy.
#[derive(Debug, Clone)]
pub struct SkyrmeParams {
    /// Coefficient A in MeV (typically negative, ~-356 MeV for K=200).
    pub a_mev: f64,
    /// Coefficient B in MeV (positive, determines saturation).
    pub b_mev: f64,
    /// Power law exponent γ (sigma).
    pub gamma: f64,
    /// Nuclear compressibility K₀ in MeV.
    pub k0_mev: f64,
    /// Symmetry energy coefficient S₀ in MeV at saturation.
    pub s0_mev: f64,
    /// Slope parameter L of symmetry energy in MeV.
    pub l_mev: f64,
}
impl SkyrmeParams {
    /// Soft EOS (K₀ ≈ 200 MeV).
    ///
    /// Parameters A and B satisfy the saturation conditions:
    /// P(ρ₀) = 0  and  u(ρ₀) = E/A ≈ −16 MeV.
    pub fn soft() -> Self {
        Self {
            a_mev: -224.0,
            b_mev: 208.0,
            gamma: 7.0 / 6.0,
            k0_mev: 200.0,
            s0_mev: 31.6,
            l_mev: 60.0,
        }
    }
    /// Stiff EOS (K₀ ≈ 380 MeV).
    pub fn stiff() -> Self {
        Self {
            a_mev: -124.0,
            b_mev: 71.0,
            gamma: 2.0,
            k0_mev: 380.0,
            s0_mev: 31.6,
            l_mev: 80.0,
        }
    }
    /// Standard soft EOS with K₀ = 240 MeV (NL3-like).
    pub fn standard() -> Self {
        Self {
            a_mev: -209.2,
            b_mev: 156.4,
            gamma: 4.0 / 3.0,
            k0_mev: 240.0,
            s0_mev: 31.6,
            l_mev: 58.0,
        }
    }
    /// Potential energy per particle from the Skyrme EOS (isospin-symmetric).
    ///
    /// `u(ρ) = A/2 · (ρ/ρ₀) + B/(γ+1) · (ρ/ρ₀)^γ`  \[MeV\]
    pub fn potential_energy_per_particle(&self, rho: f64) -> f64 {
        let x = rho / RHO_0;
        self.a_mev * 0.5 * x + self.b_mev / (self.gamma + 1.0) * x.powf(self.gamma)
    }
    /// Pressure from the Skyrme EOS (isospin-symmetric).
    ///
    /// `P(ρ) = ρ² d(u/ρ)/dρ`  \[MeV/fm³\]
    pub fn pressure(&self, rho: f64) -> f64 {
        let x = rho / RHO_0;
        let du_dx = self.a_mev * 0.5
            + self.b_mev * self.gamma / (self.gamma + 1.0) * x.powf(self.gamma - 1.0);
        rho * rho * du_dx / RHO_0
    }
    /// Symmetry energy E_sym(ρ) = S₀ + L/3 · (ρ-ρ₀)/ρ₀  \[MeV\] (linear approx).
    pub fn symmetry_energy(&self, rho: f64) -> f64 {
        self.s0_mev + self.l_mev / 3.0 * (rho - RHO_0) / RHO_0
    }
    /// Isospin-asymmetry correction to energy per particle.
    ///
    /// `δu = E_sym(ρ) · δ²`  where δ = (ρₙ - ρₚ) / ρ.
    pub fn asymmetry_energy(&self, rho: f64, rho_n: f64, rho_p: f64) -> f64 {
        if rho < 1e-20 {
            return 0.0;
        }
        let delta = (rho_n - rho_p) / rho;
        self.symmetry_energy(rho) * delta * delta
    }
}
/// Minimum spanning tree cluster finder for SPH particles.
///
/// Groups nuclear particles that are closer than `r_cut` (in fm).
#[derive(Debug, Clone)]
pub struct NuclearClusterFinder {
    /// Clustering radius cutoff in fm.
    pub r_cut: f64,
    /// Minimum density threshold ρ_min for a particle to join a cluster.
    pub rho_min: f64,
}
impl NuclearClusterFinder {
    /// Create a cluster finder with given cutoff radius and minimum density.
    pub fn new(r_cut: f64, rho_min: f64) -> Self {
        Self { r_cut, rho_min }
    }
    /// Assign cluster labels to all particles using a simple union-find approach.
    ///
    /// Returns a vector of cluster IDs (same length as `particles`).
    pub fn find_clusters(&self, particles: &[NuclearParticle]) -> Vec<usize> {
        let n = particles.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut Vec<usize>, mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
            let ra = find(parent, a);
            let rb = find(parent, b);
            if ra != rb {
                parent[ra] = rb;
            }
        }
        for i in 0..n {
            if particles[i].rho < self.rho_min {
                continue;
            }
            for j in (i + 1)..n {
                if particles[j].rho < self.rho_min {
                    continue;
                }
                let r = len3(sub3(particles[i].pos, particles[j].pos));
                if r < self.r_cut {
                    union(&mut parent, i, j);
                }
            }
        }
        let mut label_map = std::collections::HashMap::new();
        let mut next_label = 0usize;
        let mut labels = vec![0usize; n];
        for i in 0..n {
            let root = find(&mut parent, i);
            let label = *label_map.entry(root).or_insert_with(|| {
                let l = next_label;
                next_label += 1;
                l
            });
            labels[i] = label;
        }
        labels
    }
    /// Compute mass and charge of each cluster.
    ///
    /// Returns a vector of `(A, Z, centroid_pos)` tuples.
    pub fn cluster_properties(&self, particles: &[NuclearParticle]) -> Vec<(u32, u32, [f64; 3])> {
        let labels = self.find_clusters(particles);
        let n_clusters = labels.iter().copied().max().unwrap_or(0) + 1;
        let mut mass = vec![0u32; n_clusters];
        let mut charge = vec![0u32; n_clusters];
        let mut centroid = vec![[0.0_f64; 3]; n_clusters];
        let mut count = vec![0u32; n_clusters];
        for (i, p) in particles.iter().enumerate() {
            let cl = labels[i];
            mass[cl] += 1;
            if p.nucleon_type == NucleonType::Proton {
                charge[cl] += 1;
            }
            centroid[cl] = add3(centroid[cl], p.pos);
            count[cl] += 1;
        }
        (0..n_clusters)
            .map(|cl| {
                let c = count[cl].max(1) as f64;
                (mass[cl], charge[cl], scale3(centroid[cl], 1.0 / c))
            })
            .collect()
    }
}
/// Relativistic SPH state variables for a single particle.
#[derive(Debug, Clone)]
pub struct RelativisticSphParticle {
    /// Position x^μ (spatial part) in fm.
    pub pos: [f64; 3],
    /// Coordinate velocity dx/dt (in units of c).
    pub vel: [f64; 3],
    /// Lorentz factor γ.
    pub lorentz_factor: f64,
    /// Rest-frame baryon number density ρ₀ in fm⁻³.
    pub rho0: f64,
    /// Specific internal energy ε in MeV.
    pub epsilon: f64,
    /// Pressure P in MeV/fm³.
    pub pressure: f64,
    /// SPH smoothing length h in fm.
    pub h: f64,
    /// Baryon number weight.
    pub baryon_number: f64,
}
impl RelativisticSphParticle {
    /// Create a relativistic SPH particle at rest.
    pub fn at_rest(pos: [f64; 3], rho0: f64, h: f64, baryon_number: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            lorentz_factor: 1.0,
            rho0,
            epsilon: 0.0,
            pressure: 0.0,
            h,
            baryon_number,
        }
    }
    /// Update the Lorentz factor from the current velocity.
    pub fn update_lorentz_factor(&mut self) {
        let v2 = dot3(self.vel, self.vel);
        self.lorentz_factor = (1.0 - v2.min(1.0 - 1e-10)).powf(-0.5);
    }
    /// Total energy density in the lab frame e = γ² (ε + P) ρ₀ - P \[MeV/fm³\].
    pub fn energy_density_lab(&self) -> f64 {
        let gamma2 = self.lorentz_factor * self.lorentz_factor;
        gamma2 * (self.epsilon + self.pressure / self.rho0.max(1e-30)) * self.rho0 - self.pressure
    }
    /// Four-velocity u^i = γ v^i (spatial components only).
    pub fn four_velocity_spatial(&self) -> [f64; 3] {
        scale3(self.vel, self.lorentz_factor)
    }
}
/// A nuclear SPH particle (pseudo-particle / test particle).
#[derive(Debug, Clone)]
pub struct NuclearParticle {
    /// Position in fm.
    pub pos: [f64; 3],
    /// Momentum per unit mass (velocity × γ for relativistic, or just velocity).
    pub vel: [f64; 3],
    /// Force accumulator in MeV/fm.
    pub force: [f64; 3],
    /// SPH smoothing length h in fm.
    pub h: f64,
    /// Local baryon number density ρ in fm⁻³ (computed by SPH sum).
    pub rho: f64,
    /// Internal energy per particle u in MeV.
    pub u: f64,
    /// Baryon number carried by this pseudo-particle (usually 1/N_test).
    pub baryon_number: f64,
    /// Nucleon type: proton or neutron.
    pub nucleon_type: NucleonType,
    /// Isospin projection τ₃: +1/2 for proton, −1/2 for neutron.
    pub isospin: f64,
    /// Whether this particle has been absorbed into a fragment.
    pub in_fragment: bool,
}
impl NuclearParticle {
    /// Create a new nuclear particle.
    pub fn new(
        pos: [f64; 3],
        vel: [f64; 3],
        h: f64,
        baryon_number: f64,
        nucleon_type: NucleonType,
    ) -> Self {
        let isospin = match nucleon_type {
            NucleonType::Proton => 0.5,
            NucleonType::Neutron => -0.5,
        };
        Self {
            pos,
            vel,
            force: [0.0; 3],
            h,
            rho: RHO_0,
            u: 0.0,
            baryon_number,
            nucleon_type,
            isospin,
            in_fragment: false,
        }
    }
    /// Lorentz factor γ for the particle's velocity (c = 1).
    pub fn lorentz_factor(&self) -> f64 {
        let v2 = dot3(self.vel, self.vel);
        (1.0 - v2).max(1e-10).powf(-0.5)
    }
    /// Relativistic momentum p = γ m v per unit nucleon mass.
    pub fn relativistic_momentum(&self) -> [f64; 3] {
        let gamma = self.lorentz_factor();
        scale3(self.vel, gamma * NUCLEON_MASS)
    }
    /// Kinetic energy T = (γ - 1) m c²  \[MeV\].
    pub fn kinetic_energy(&self) -> f64 {
        let gamma = self.lorentz_factor();
        (gamma - 1.0) * NUCLEON_MASS * self.baryon_number
    }
}
/// High-level nuclear SPH simulation driver.
#[derive(Debug)]
pub struct NuclearSphSimulation {
    /// List of SPH particles.
    pub particles: Vec<NuclearParticle>,
    /// Skyrme EOS parameters.
    pub eos: SkyrmeParams,
    /// Pauli exclusion force parameters.
    pub pauli: PauliParams,
    /// Global SPH smoothing length h in fm.
    pub h: f64,
    /// Current simulation time in fm/c.
    pub time: f64,
    /// Current time step Δt in fm/c.
    pub dt: f64,
    /// Step counter.
    pub step: u64,
    /// Artificial viscosity α coefficient.
    pub alpha_visc: f64,
    /// Artificial viscosity β coefficient.
    pub beta_visc: f64,
}
impl NuclearSphSimulation {
    /// Create a new nuclear SPH simulation.
    pub fn new(eos: SkyrmeParams, h: f64, dt: f64) -> Self {
        Self {
            particles: Vec::new(),
            eos,
            pauli: PauliParams::default(),
            h,
            time: 0.0,
            dt,
            step: 0,
            alpha_visc: 1.0,
            beta_visc: 2.0,
        }
    }
    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: NuclearParticle) {
        self.particles.push(p);
    }
    /// Recompute baryon densities for all particles via SPH sum.
    pub fn update_densities(&mut self) {
        let positions: Vec<[f64; 3]> = self.particles.iter().map(|p| p.pos).collect();
        let baryons: Vec<f64> = self.particles.iter().map(|p| p.baryon_number).collect();
        let h = self.h;
        for i in 0..self.particles.len() {
            let mut rho = 0.0;
            for j in 0..self.particles.len() {
                let r = len3(sub3(positions[i], positions[j]));
                rho += baryons[j] * gaussian_kernel(r, h);
            }
            self.particles[i].rho = rho;
        }
    }
    /// Recompute forces on all particles.
    pub fn update_forces(&mut self) {
        let n = self.particles.len();
        let snapshot: Vec<NuclearParticle> = self.particles.clone();
        let h = self.h;
        for i in 0..n {
            let neighbours: Vec<NuclearParticle> = snapshot
                .iter()
                .enumerate()
                .filter(|(j, _)| *j != i)
                .map(|(_, p)| p.clone())
                .collect();
            let f_skyrme = skyrme_force(&snapshot[i], &neighbours, &self.eos, h);
            let f_visc = artificial_viscosity_force(
                &snapshot[i],
                &neighbours,
                h,
                self.alpha_visc,
                self.beta_visc,
            );
            let mut f_pauli = [0.0_f64; 3];
            for nb in &neighbours {
                let same = nb.nucleon_type == snapshot[i].nucleon_type;
                let r_ij = sub3(snapshot[i].pos, nb.pos);
                let p_ij = sub3(
                    scale3(snapshot[i].vel, NUCLEON_MASS),
                    scale3(nb.vel, NUCLEON_MASS),
                );
                let fp = pauli_force(&self.pauli, r_ij, p_ij, same);
                f_pauli = add3(f_pauli, fp);
            }
            self.particles[i].force = add3(add3(f_skyrme, f_visc), f_pauli);
        }
    }
    /// Advance the simulation by one time step (velocity Verlet).
    pub fn step_forward(&mut self) {
        for p in &mut self.particles {
            verlet_kick(p, self.dt);
        }
        for p in &mut self.particles {
            verlet_drift(p, self.dt);
        }
        self.update_densities();
        self.update_forces();
        for p in &mut self.particles {
            verlet_kick(p, self.dt);
        }
        self.time += self.dt;
        self.step += 1;
    }
    /// Total kinetic energy of all particles in MeV.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }
    /// Total baryon number (should be conserved).
    pub fn total_baryon_number(&self) -> f64 {
        self.particles.iter().map(|p| p.baryon_number).sum()
    }
    /// Number of protons and neutrons.
    pub fn proton_neutron_count(&self) -> (usize, usize) {
        let p = self
            .particles
            .iter()
            .filter(|par| par.nucleon_type == NucleonType::Proton)
            .count();
        let n = self
            .particles
            .iter()
            .filter(|par| par.nucleon_type == NucleonType::Neutron)
            .count();
        (p, n)
    }
}
/// Parameters controlling the quantum Pauli exclusion force.
///
/// The force prevents identical nucleons from occupying the same
/// phase-space cell — a key ingredient in QMD and pBUU models.
#[derive(Debug, Clone)]
pub struct PauliParams {
    /// Phase-space cell size in (fm · MeV/c)³.
    pub phase_space_cell: f64,
    /// Pauli potential strength in MeV.
    pub strength_mev: f64,
    /// Width parameter q₀ in fm⁻¹ for position-space Gaussian.
    pub q0_fm_inv: f64,
    /// Width parameter p₀ in MeV/c for momentum-space Gaussian.
    pub p0_mev_c: f64,
}
/// Parameters for the QMD nuclear potential.
///
/// V_QMD = t₁ δ(r_ij) + t₂ δ²(r_ij)  (in Skyrme form for SPH)
#[derive(Debug, Clone)]
pub struct QmdParams {
    /// Gaussian wavepacket width L in fm².
    pub wavepacket_width_fm2: f64,
    /// Two-body interaction t₁ in MeV·fm³.
    pub t1_mev_fm3: f64,
    /// Three-body interaction t₃ in MeV·fm^{3γ}.
    pub t3_mev: f64,
    /// Power γ for three-body density dependence.
    pub gamma: f64,
    /// Symmetry energy coefficient cs in MeV.
    pub cs_mev: f64,
    /// Coulomb interaction coefficient e²/(4πε₀) in MeV·fm.
    pub coulomb_coeff: f64,
    /// Momentum-dependent interaction coefficient in MeV.
    pub momentum_dep_coeff: f64,
}
/// Identity of a nuclear SPH pseudo-particle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NucleonType {
    /// Proton (charge +1).
    Proton,
    /// Neutron (charge 0).
    Neutron,
}
/// Represents a fission fragment with mass, charge, and kinetic energy.
#[derive(Debug, Clone)]
pub struct FissionFragment {
    /// Fragment mass number A.
    pub mass_number: u32,
    /// Fragment charge Z.
    pub charge: u32,
    /// Centre of mass position in fm.
    pub pos: [f64; 3],
    /// Centre of mass velocity (units of c).
    pub vel: [f64; 3],
    /// Excitation energy E* in MeV.
    pub excitation_energy_mev: f64,
    /// Deformation parameter β₂.
    pub deformation: f64,
}
impl FissionFragment {
    /// Create a fission fragment.
    pub fn new(mass_number: u32, charge: u32, pos: [f64; 3], vel: [f64; 3]) -> Self {
        Self {
            mass_number,
            charge,
            pos,
            vel,
            excitation_energy_mev: 0.0,
            deformation: 0.0,
        }
    }
    /// Total kinetic energy TKE in MeV (non-relativistic).
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = dot3(self.vel, self.vel);
        0.5 * (self.mass_number as f64) * NUCLEON_MASS * v2
    }
    /// Coulomb repulsion energy between two fragments at their separation.
    pub fn coulomb_energy_with(&self, other: &FissionFragment) -> f64 {
        let r_vec = sub3(self.pos, other.pos);
        let r = len3(r_vec);
        if r < 1e-10 {
            return 0.0;
        }
        1.44 * (self.charge as f64) * (other.charge as f64) / r
    }
    /// Acceleration due to Coulomb force from another fragment \[c²/fm\].
    pub fn coulomb_acceleration(&self, other: &FissionFragment) -> [f64; 3] {
        let r_vec = sub3(self.pos, other.pos);
        let r = len3(r_vec);
        if r < 1e-10 {
            return [0.0; 3];
        }
        let f_mag = 1.44 * (self.charge as f64) * (other.charge as f64) / (r * r * r);
        let m = (self.mass_number as f64) * NUCLEON_MASS;
        scale3(r_vec, f_mag / m)
    }
}
/// Scission point geometry for symmetric and asymmetric fission.
#[derive(Debug, Clone)]
pub struct ScissionConfig {
    /// Heavy fragment mass number A_H.
    pub a_heavy: u32,
    /// Light fragment mass number A_L.
    pub a_light: u32,
    /// Nuclear charge of heavy fragment Z_H.
    pub z_heavy: u32,
    /// Nuclear charge of light fragment Z_L.
    pub z_light: u32,
    /// Scission separation distance d in fm.
    pub separation_fm: f64,
    /// Total kinetic energy release TKE in MeV.
    pub tke_mev: f64,
}
impl ScissionConfig {
    /// Create symmetric fission configuration.
    pub fn symmetric(a_compound: u32, z_compound: u32) -> Self {
        let a_h = a_compound / 2;
        let a_l = a_compound - a_h;
        let z_h = z_compound / 2;
        let z_l = z_compound - z_h;
        let tke = 0.1071 * (z_compound as f64).powi(2) / ((a_compound as f64).powf(1.0 / 3.0));
        Self {
            a_heavy: a_h,
            a_light: a_l,
            z_heavy: z_h,
            z_light: z_l,
            separation_fm: 1.2 * (a_h as f64).powf(1.0 / 3.0)
                + 1.2 * (a_l as f64).powf(1.0 / 3.0)
                + 2.0,
            tke_mev: tke,
        }
    }
    /// Generate the two fission fragments from this scission configuration.
    pub fn generate_fragments(&self) -> (FissionFragment, FissionFragment) {
        let half_sep = self.separation_fm * 0.5;
        let heavy =
            FissionFragment::new(self.a_heavy, self.z_heavy, [0.0, 0.0, half_sep], [0.0; 3]);
        let light =
            FissionFragment::new(self.a_light, self.z_light, [0.0, 0.0, -half_sep], [0.0; 3]);
        (heavy, light)
    }
}
/// Initial condition builder for a heavy-ion collision event.
#[derive(Debug, Clone)]
pub struct HeavyIonCollision {
    /// Projectile nucleus profile.
    pub projectile: WoodsSaxonProfile,
    /// Target nucleus profile.
    pub target: WoodsSaxonProfile,
    /// Impact parameter b in fm.
    pub impact_parameter_fm: f64,
    /// Beam energy per nucleon in MeV (lab frame).
    pub beam_energy_mev: f64,
    /// Number of SPH test particles per nucleon.
    pub test_particles: u32,
}
impl HeavyIonCollision {
    /// Construct a collision system.
    pub fn new(
        projectile: WoodsSaxonProfile,
        target: WoodsSaxonProfile,
        impact_parameter_fm: f64,
        beam_energy_mev: f64,
        test_particles: u32,
    ) -> Self {
        Self {
            projectile,
            target,
            impact_parameter_fm,
            beam_energy_mev,
            test_particles,
        }
    }
    /// Beam velocity β = v/c from lab-frame kinetic energy.
    pub fn beam_velocity(&self) -> f64 {
        let t = self.beam_energy_mev;
        let m = NUCLEON_MASS;
        let gamma = 1.0 + t / m;
        let beta2 = 1.0 - 1.0 / (gamma * gamma);
        beta2.max(0.0).sqrt()
    }
    /// Lorentz γ factor for the beam.
    pub fn beam_lorentz_factor(&self) -> f64 {
        1.0 + self.beam_energy_mev / NUCLEON_MASS
    }
    /// Longitudinal Lorentz contraction factor.
    pub fn lorentz_contraction(&self) -> f64 {
        self.beam_lorentz_factor()
    }
    /// Estimate total nucleon-nucleon cross-section in mb (empirical).
    pub fn nn_cross_section_mb(&self) -> f64 {
        let e_mev = self.beam_energy_mev;
        if e_mev < 300.0 {
            17.0 + 32.0 * (1.0 - e_mev / 300.0).powi(2)
        } else {
            40.0 + 30.0 / (1.0 + (e_mev / 1000.0).powi(2))
        }
    }
    /// Number of participating nucleons N_part (Glauber model, approximation).
    pub fn n_part_glauber(&self) -> f64 {
        let b = self.impact_parameter_fm;
        let r_p = self.projectile.radius_fm;
        let r_t = self.target.radius_fm;
        let sigma_nn = self.nn_cross_section_mb() * 0.1 / (PI);
        let r_sum = r_p + r_t;
        if b >= r_sum {
            return 0.0;
        }
        let a_p = self.projectile.mass_number as f64;
        let a_t = self.target.mass_number as f64;
        let overlap = (r_sum - b).max(0.0) / r_sum;
        let t_overlap = overlap * overlap * PI * sigma_nn;
        a_p + a_t - (a_p * (-a_t * t_overlap).exp() + a_t * (-a_p * t_overlap).exp())
    }
}
