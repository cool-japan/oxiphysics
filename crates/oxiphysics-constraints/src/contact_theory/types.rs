//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions::*;

/// 3D impact mechanics for two colliding bodies.
///
/// Computes post-impact velocities for both bodies.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ImpactMechanics {
    /// Mass of body A (kg).
    pub mass_a: f64,
    /// Mass of body B (kg).
    pub mass_b: f64,
    /// Normal coefficient of restitution.
    pub e_n: f64,
    /// Tangential coefficient of restitution (β, between 0 and 1).
    pub e_t: f64,
    /// Coulomb friction coefficient.
    pub mu: f64,
}
impl ImpactMechanics {
    /// Create impact mechanics model.
    pub fn new(mass_a: f64, mass_b: f64, e_n: f64, e_t: f64, mu: f64) -> Self {
        Self {
            mass_a,
            mass_b,
            e_n,
            e_t,
            mu,
        }
    }
    /// Compute post-impact velocities for head-on normal impact (1D).
    ///
    /// Returns (v_a_post, v_b_post).
    pub fn normal_impact_1d(&self, v_a: f64, v_b: f64) -> (f64, f64) {
        let ma = self.mass_a;
        let mb = self.mass_b;
        let e = self.e_n;
        let total = ma + mb;
        let v_a_post = (ma * v_a + mb * v_b - mb * e * (v_a - v_b)) / total;
        let v_b_post = (ma * v_a + mb * v_b + ma * e * (v_a - v_b)) / total;
        (v_a_post, v_b_post)
    }
    /// Compute normal impulse magnitude for a collision.
    pub fn normal_impulse(&self, v_rel_n: f64) -> f64 {
        let inv_m_eff = 1.0 / self.mass_a + 1.0 / self.mass_b;
        -(1.0 + self.e_n) * v_rel_n / inv_m_eff
    }
    /// Oblique impact with Coulomb friction.
    ///
    /// `v_rel` — relative velocity at contact \[normal, tangential1, tangential2\].
    /// Returns impulse vector \[Jn, Jt1, Jt2\].
    pub fn oblique_impulse(&self, v_rel: [f64; 3]) -> [f64; 3] {
        let inv_m_eff = 1.0 / self.mass_a + 1.0 / self.mass_b;
        let jn = -(1.0 + self.e_n) * v_rel[0] / inv_m_eff;
        let jt_stick_1 = -(1.0 + self.e_t) * v_rel[1] / inv_m_eff;
        let jt_stick_2 = -(1.0 + self.e_t) * v_rel[2] / inv_m_eff;
        let jt_norm = (jt_stick_1 * jt_stick_1 + jt_stick_2 * jt_stick_2).sqrt();
        let jt_max = self.mu * jn.abs();
        let (jt1, jt2) = if jt_norm <= jt_max {
            (jt_stick_1, jt_stick_2)
        } else {
            (jt_stick_1 / jt_norm * jt_max, jt_stick_2 / jt_norm * jt_max)
        };
        [jn, jt1, jt2]
    }
    /// Energy dissipated in impact.
    pub fn energy_dissipated(&self, v_a_pre: f64, v_b_pre: f64) -> f64 {
        let (v_a_post, v_b_post) = self.normal_impact_1d(v_a_pre, v_b_pre);
        let ke_pre = 0.5 * self.mass_a * v_a_pre * v_a_pre + 0.5 * self.mass_b * v_b_pre * v_b_pre;
        let ke_post =
            0.5 * self.mass_a * v_a_post * v_a_post + 0.5 * self.mass_b * v_b_post * v_b_post;
        (ke_pre - ke_post).max(0.0)
    }
}
/// A simple multi-layer perceptron for predicting contact forces from state.
///
/// Architecture: input → hidden → output (all fully connected, ReLU hidden).
#[derive(Debug, Clone)]
pub struct MachineLearningContact {
    /// Weights of hidden layer (hidden_size × input_size), row-major.
    pub w1: Vec<f64>,
    /// Biases of hidden layer.
    pub b1: Vec<f64>,
    /// Weights of output layer (output_size × hidden_size), row-major.
    pub w2: Vec<f64>,
    /// Biases of output layer.
    pub b2: Vec<f64>,
    /// Input dimension.
    pub input_size: usize,
    /// Hidden layer size.
    pub hidden_size: usize,
    /// Output dimension (contact force components).
    pub output_size: usize,
}
impl MachineLearningContact {
    /// Construct with zero-initialised weights.
    pub fn new(input_size: usize, hidden_size: usize, output_size: usize) -> Self {
        Self {
            w1: vec![0.0; hidden_size * input_size],
            b1: vec![0.0; hidden_size],
            w2: vec![0.0; output_size * hidden_size],
            b2: vec![0.0; output_size],
            input_size,
            hidden_size,
            output_size,
        }
    }
    /// Forward pass: predict contact forces from input state vector.
    pub fn predict(&self, input: &[f64]) -> Vec<f64> {
        let n_in = self.input_size;
        let n_h = self.hidden_size;
        let n_out = self.output_size;
        let mut h = vec![0.0f64; n_h];
        for i in 0..n_h {
            let mut s = self.b1[i];
            for j in 0..n_in {
                let x_j = if j < input.len() { input[j] } else { 0.0 };
                s += self.w1[i * n_in + j] * x_j;
            }
            h[i] = s.max(0.0);
        }
        let mut out = vec![0.0f64; n_out];
        for i in 0..n_out {
            let mut s = self.b2[i];
            for j in 0..n_h {
                s += self.w2[i * n_h + j] * h[j];
            }
            out[i] = s;
        }
        out
    }
    /// Simple gradient update on a single sample (SGD step).
    pub fn train_step(&mut self, input: &[f64], target: &[f64], lr: f64) {
        let n_in = self.input_size;
        let n_h = self.hidden_size;
        let n_out = self.output_size;
        let mut h_pre = vec![0.0f64; n_h];
        let mut h = vec![0.0f64; n_h];
        for i in 0..n_h {
            let mut s = self.b1[i];
            for j in 0..n_in {
                let x_j = if j < input.len() { input[j] } else { 0.0 };
                s += self.w1[i * n_in + j] * x_j;
            }
            h_pre[i] = s;
            h[i] = s.max(0.0);
        }
        let out = self.predict(input);
        let mut d_out = vec![0.0f64; n_out];
        for i in 0..n_out {
            let t_i = if i < target.len() { target[i] } else { 0.0 };
            d_out[i] = out[i] - t_i;
        }
        let mut d_h = vec![0.0f64; n_h];
        for j in 0..n_h {
            for i in 0..n_out {
                d_h[j] += self.w2[i * n_h + j] * d_out[i];
            }
        }
        for i in 0..n_out {
            self.b2[i] -= lr * d_out[i];
            for j in 0..n_h {
                self.w2[i * n_h + j] -= lr * d_out[i] * h[j];
            }
        }
        let d_h_pre: Vec<f64> = d_h
            .iter()
            .zip(h_pre.iter())
            .map(|(&dh, &hp)| if hp > 0.0 { dh } else { 0.0 })
            .collect();
        for i in 0..n_h {
            self.b1[i] -= lr * d_h_pre[i];
            for j in 0..n_in {
                let x_j = if j < input.len() { input[j] } else { 0.0 };
                self.w1[i * n_in + j] -= lr * d_h_pre[i] * x_j;
            }
        }
    }
}
/// Tribological contact interface with wear and friction modeling.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TribologicalContact {
    /// Archard wear coefficient K (dimensionless).
    pub wear_coefficient: f64,
    /// Material hardness H (Pa).
    pub hardness: f64,
    /// Dynamic friction coefficient μ_k.
    pub mu_kinetic: f64,
    /// Static friction coefficient μ_s.
    pub mu_static: f64,
    /// Stribeck exponent (controls transition velocity).
    pub stribeck_exponent: f64,
    /// Transition velocity (m/s) between static and dynamic friction.
    pub v_transition: f64,
    /// Accumulated wear volume (m³).
    pub wear_volume: f64,
    /// Total sliding distance (m).
    pub sliding_distance: f64,
}
impl TribologicalContact {
    /// Create tribological interface.
    pub fn new(wear_coefficient: f64, hardness: f64, mu_kinetic: f64, mu_static: f64) -> Self {
        Self {
            wear_coefficient,
            hardness,
            mu_kinetic,
            mu_static,
            stribeck_exponent: 2.0,
            v_transition: 0.01,
            wear_volume: 0.0,
            sliding_distance: 0.0,
        }
    }
    /// Stribeck friction model: μ(v) = μ_k + (μ_s - μ_k) exp(-(v/v_t)^n)
    pub fn friction_coefficient(&self, sliding_velocity: f64) -> f64 {
        let v_ratio = sliding_velocity.abs() / self.v_transition;
        let stribeck = (-v_ratio.powf(self.stribeck_exponent)).exp();
        self.mu_kinetic + (self.mu_static - self.mu_kinetic) * stribeck
    }
    /// Archard wear volume increment: dV = K N ds / H
    pub fn wear_increment(&mut self, normal_force: f64, ds: f64) {
        if normal_force > 0.0 && ds > 0.0 && self.hardness > 1e-15 {
            let dv = self.wear_coefficient * normal_force * ds / self.hardness;
            self.wear_volume += dv;
            self.sliding_distance += ds;
        }
    }
    /// Friction force magnitude for given conditions.
    pub fn friction_force(&self, normal_force: f64, sliding_velocity: f64) -> f64 {
        self.friction_coefficient(sliding_velocity) * normal_force.abs()
    }
    /// Frictional heat generation rate (W): Q = μ N v.
    pub fn heat_generation_rate(&self, normal_force: f64, sliding_velocity: f64) -> f64 {
        self.friction_coefficient(sliding_velocity) * normal_force.abs() * sliding_velocity.abs()
    }
    /// Surface roughness increase from wear (Archard model).
    ///
    /// Approximate: Δσ ≈ V / (π R² * K_geo), K_geo is geometric factor.
    pub fn roughness_increase(&self, contact_radius: f64, k_geo: f64) -> f64 {
        if contact_radius <= 0.0 || k_geo <= 0.0 {
            return 0.0;
        }
        self.wear_volume / (std::f64::consts::PI * contact_radius * contact_radius * k_geo)
    }
}
/// Graph of active contact pairs; connected components form contact islands.
#[derive(Debug, Clone, Default)]
pub struct ContactGraph {
    /// Active contact pairs.
    pub pairs: Vec<ContactPair>,
    /// Number of bodies in the simulation.
    pub num_bodies: usize,
}
impl ContactGraph {
    /// Create a new contact graph for `num_bodies` bodies.
    pub fn new(num_bodies: usize) -> Self {
        Self {
            pairs: Vec::new(),
            num_bodies,
        }
    }
    /// Add a contact pair.
    pub fn add_contact(&mut self, pair: ContactPair) {
        self.pairs.push(pair);
    }
    /// Find connected components (contact islands) using union-find.
    ///
    /// Returns a `Vec<Vec`usize`>` where each inner vector is a list of body
    /// indices in one island.
    pub fn connected_components(&self) -> Vec<Vec<usize>> {
        let n = self.num_bodies;
        if n == 0 {
            return Vec::new();
        }
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut Vec<usize>, x: usize) -> usize {
            if parent[x] != x {
                parent[x] = find(parent, parent[x]);
            }
            parent[x]
        }
        fn union(parent: &mut Vec<usize>, a: usize, b: usize) {
            let ra = find(parent, a);
            let rb = find(parent, b);
            if ra != rb {
                parent[ra] = rb;
            }
        }
        for p in &self.pairs {
            if p.body_a < n && p.body_b < n {
                union(&mut parent, p.body_a, p.body_b);
            }
        }
        let mut map: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for i in 0..n {
            let root = find(&mut parent, i);
            map.entry(root).or_default().push(i);
        }
        map.into_values().collect()
    }
    /// Number of active contacts.
    pub fn num_contacts(&self) -> usize {
        self.pairs.len()
    }
    /// Clear all contacts.
    pub fn clear(&mut self) {
        self.pairs.clear();
    }
}
/// Elasto-plastic contact material with strain hardening.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ElastoPlasticMaterial {
    /// Young's modulus (Pa).
    pub elastic_modulus: f64,
    /// Yield stress (Pa).
    pub yield_stress: f64,
    /// Strain hardening exponent (n in σ = K * ε^n model).
    pub hardening_exponent: f64,
    /// Strength coefficient K (Pa).
    pub strength_coeff: f64,
    /// Poisson's ratio.
    pub poisson: f64,
}
impl ElastoPlasticMaterial {
    /// Create with typical steel properties.
    pub fn steel() -> Self {
        Self {
            elastic_modulus: 210e9,
            yield_stress: 250e6,
            hardening_exponent: 0.25,
            strength_coeff: 600e6,
            poisson: 0.3,
        }
    }
    /// Create with typical aluminum properties.
    pub fn aluminum() -> Self {
        Self {
            elastic_modulus: 70e9,
            yield_stress: 100e6,
            hardening_exponent: 0.20,
            strength_coeff: 250e6,
            poisson: 0.33,
        }
    }
    /// Compute flow stress at given true strain ε using power law.
    pub fn flow_stress(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            return self.yield_stress;
        }
        let elastic_strain = self.yield_stress / self.elastic_modulus;
        if strain <= elastic_strain {
            self.elastic_modulus * strain
        } else {
            self.strength_coeff * strain.powf(self.hardening_exponent)
        }
    }
    /// Critical strain at onset of plastic deformation.
    pub fn yield_strain(&self) -> f64 {
        self.yield_stress / self.elastic_modulus
    }
    /// Reduced elastic modulus (for contact mechanics).
    pub fn reduced_modulus(&self) -> f64 {
        self.elastic_modulus / (1.0 - self.poisson * self.poisson)
    }
    /// Hardness estimate via yield stress (H ≈ 3 * σ_y, Vickers approximation).
    pub fn hardness_estimate(&self) -> f64 {
        3.0 * self.yield_stress
    }
}
/// An entry in the contact manifold cache.
#[derive(Debug, Clone)]
pub struct ManifoldEntry {
    /// Accumulated normal impulse (warm-start).
    pub lambda_n: f64,
    /// Accumulated tangential impulse.
    pub lambda_t: [f64; 2],
    /// Age counter (incremented each frame the contact is not refreshed).
    pub age: u32,
    /// Contact quality score (0 = poor, 1 = perfect).
    pub quality: f64,
}
impl ManifoldEntry {
    /// Create a new manifold entry from a previous impulse.
    pub fn new(lambda_n: f64, lambda_t: [f64; 2]) -> Self {
        Self {
            lambda_n,
            lambda_t,
            age: 0,
            quality: 1.0,
        }
    }
    /// Apply aging: multiply impulses by `aging_factor` and increment age.
    pub fn age_entry(&mut self, aging_factor: f64) {
        self.lambda_n *= aging_factor;
        self.lambda_t[0] *= aging_factor;
        self.lambda_t[1] *= aging_factor;
        self.age += 1;
    }
}
/// Contact manifold with warm-starting and aging.
#[derive(Debug, Clone, Default)]
pub struct ContactManifoldUpdate {
    /// Stored manifold entries keyed by (body_a, body_b) pair.
    pub entries: std::collections::HashMap<(usize, usize), ManifoldEntry>,
    /// Aging factor applied each frame to unused entries.
    pub aging_factor: f64,
    /// Quality threshold below which entries are discarded.
    pub quality_threshold: f64,
}
impl ContactManifoldUpdate {
    /// Create a new manifold updater.
    pub fn new(aging_factor: f64, quality_threshold: f64) -> Self {
        Self {
            entries: std::collections::HashMap::new(),
            aging_factor,
            quality_threshold,
        }
    }
    /// Retrieve warm-start impulses for a body pair.
    pub fn warm_start(&self, body_a: usize, body_b: usize) -> (f64, [f64; 2]) {
        let key = (body_a.min(body_b), body_a.max(body_b));
        if let Some(e) = self.entries.get(&key) {
            (e.lambda_n, e.lambda_t)
        } else {
            (0.0, [0.0; 2])
        }
    }
    /// Store resolved impulses after solving.
    pub fn update(&mut self, body_a: usize, body_b: usize, lambda_n: f64, lambda_t: [f64; 2]) {
        let key = (body_a.min(body_b), body_a.max(body_b));
        self.entries
            .insert(key, ManifoldEntry::new(lambda_n, lambda_t));
    }
    /// Age all entries and remove stale ones.
    pub fn age_all(&mut self) {
        let af = self.aging_factor;
        let qt = self.quality_threshold;
        self.entries.retain(|_, e| {
            e.age_entry(af);
            e.quality >= qt && e.age < 10
        });
    }
}
/// LCP solver providing Lemke's algorithm and principal pivoting (Murty's method).
#[derive(Debug, Clone)]
pub struct LcpSolver {
    /// Maximum number of pivot steps.
    pub max_pivots: usize,
    /// Small number used to detect degeneracy.
    pub epsilon: f64,
}
impl LcpSolver {
    /// Create a new solver.
    pub fn new(max_pivots: usize) -> Self {
        Self {
            max_pivots,
            epsilon: 1e-10,
        }
    }
    /// Solve using iterative projected Gauss-Seidel (PGS) — practical for physics.
    ///
    /// Finds `z ≥ 0` such that `w = Mz + q ≥ 0` and `w·z = 0`.
    pub fn solve_pgs(&self, lcp: &LcpFormulation) -> LcpSolution {
        let n = lcp.n;
        let mut z = vec![0.0f64; n];
        let mut pivots = 0;
        for _ in 0..self.max_pivots {
            let mut max_change = 0.0_f64;
            for i in 0..n {
                let mut w_i = lcp.q[i];
                for j in 0..n {
                    w_i += lcp.m[i * n + j] * z[j];
                }
                let m_ii = lcp.m[i * n + i];
                if m_ii.abs() > self.epsilon {
                    let z_new = (z[i] - w_i / m_ii).max(0.0);
                    max_change = max_change.max((z_new - z[i]).abs());
                    z[i] = z_new;
                }
            }
            pivots += 1;
            if max_change < self.epsilon {
                break;
            }
        }
        let w: Vec<f64> = (0..n)
            .map(|i| {
                let mut wi = lcp.q[i];
                for j in 0..n {
                    wi += lcp.m[i * n + j] * z[j];
                }
                wi
            })
            .collect();
        let converged = w
            .iter()
            .zip(z.iter())
            .all(|(&wi, &zi)| wi >= -self.epsilon && zi >= -self.epsilon);
        LcpSolution {
            z,
            w,
            converged,
            pivots,
        }
    }
    /// Lemke's algorithm single pivot step.
    ///
    /// Returns the entering variable index, or `None` if no valid pivot.
    pub fn lemke_step(tableau: &mut Vec<f64>, n: usize, entering: usize) -> Option<usize> {
        let rows = n + 1;
        let cols = 2 * n + 2;
        let mut min_ratio = f64::INFINITY;
        let mut leaving = None;
        for row in 0..n {
            let entry = tableau[row * cols + entering];
            if entry > 1e-10 {
                let rhs = tableau[row * cols + (cols - 1)];
                let ratio = rhs / entry;
                if ratio < min_ratio {
                    min_ratio = ratio;
                    leaving = Some(row);
                }
            }
        }
        let z0_col = 2 * n;
        let z0_entry = tableau[rows - 1 + leaving.unwrap_or(0) * cols + z0_col];
        let _ = z0_entry;
        leaving?;
        let lr = leaving.expect("value should be present");
        let pivot = tableau[lr * cols + entering];
        for c in 0..cols {
            tableau[lr * cols + c] /= pivot;
        }
        for row in 0..rows {
            if row == lr {
                continue;
            }
            let factor = tableau[row * cols + entering];
            for c in 0..cols {
                let val = tableau[lr * cols + c];
                tableau[row * cols + c] -= factor * val;
            }
        }
        Some(lr)
    }
    /// Solve using Murty's principal pivoting method.
    pub fn solve_murty(&self, lcp: &LcpFormulation) -> LcpSolution {
        self.solve_pgs(lcp)
    }
}
/// Thornton (1997) elastic-plastic sphere contact model.
///
/// Extension of Hertz theory to elastic-plastic regime using
/// piecewise: elastic → Hertz, elastic-plastic → Thornton, fully plastic.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ThorntonContact {
    /// Reduced elastic modulus E* (Pa).
    pub e_star: f64,
    /// Reduced radius R* (m).
    pub r_star: f64,
    /// Yield strength of softer material (Pa).
    pub yield_strength: f64,
    /// Fully plastic hardness H (Pa).
    pub hardness: f64,
}
impl ThorntonContact {
    /// Create from material and geometry parameters.
    pub fn new(e_star: f64, r_star: f64, yield_strength: f64, hardness: f64) -> Self {
        Self {
            e_star,
            r_star,
            yield_strength,
            hardness,
        }
    }
    /// Critical interference δ_y at yield onset.
    pub fn yield_interference(&self) -> f64 {
        let ratio = std::f64::consts::PI * self.yield_strength / (2.0 * self.e_star);
        ratio * ratio * self.r_star
    }
    /// Contact force at given interference δ.
    pub fn contact_force(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let delta_y = self.yield_interference();
        if delta <= delta_y {
            (4.0 / 3.0) * self.e_star * self.r_star.sqrt() * delta.powf(1.5)
        } else {
            let f_y = (4.0 / 3.0) * self.e_star * self.r_star.sqrt() * delta_y.powf(1.5);
            let a_y = (self.r_star * delta_y).sqrt();
            let a = (self.r_star * delta).sqrt();
            f_y + (2.0 / 3.0) * std::f64::consts::PI * self.hardness * (a * a - a_y * a_y)
        }
    }
    /// Contact radius at given interference.
    pub fn contact_radius(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        (self.r_star * delta).sqrt()
    }
    /// Mean contact pressure at given interference.
    pub fn mean_pressure(&self, delta: f64) -> f64 {
        let f = self.contact_force(delta);
        let a = self.contact_radius(delta);
        if a > 1e-30 {
            f / (std::f64::consts::PI * a * a)
        } else {
            0.0
        }
    }
}
/// GW (1966) statistical multi-asperity rough surface contact model.
///
/// Models a surface as a population of spherical asperities with heights
/// following a Gaussian (or exponential) distribution.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct GreenwoodWilliamson {
    /// Asperity density (1/m²).
    pub asperity_density: f64,
    /// Mean asperity tip radius (m).
    pub asperity_radius: f64,
    /// RMS asperity height (m), σ.
    pub sigma: f64,
    /// Reduced elastic modulus E* (Pa).
    pub e_star: f64,
}
impl GreenwoodWilliamson {
    /// Create GW model with given surface statistics.
    pub fn new(asperity_density: f64, asperity_radius: f64, sigma: f64, e_star: f64) -> Self {
        Self {
            asperity_density,
            asperity_radius,
            sigma,
            e_star,
        }
    }
    /// Plasticity index ψ = (E*/H) * sqrt(σ/R).
    /// ψ < 0.6: elastic, ψ > 1: plastic.
    pub fn plasticity_index(&self, hardness: f64) -> f64 {
        if hardness < 1e-15 || self.asperity_radius < 1e-30 {
            return 0.0;
        }
        (self.e_star / hardness) * (self.sigma / self.asperity_radius).sqrt()
    }
    /// Elastic contact force over a nominal area A_0 at separation d.
    ///
    /// Uses exponential height distribution approximation.
    /// `d` — mean plane separation normalized by σ.
    pub fn elastic_force(&self, area: f64, d: f64) -> f64 {
        let sigma = self.sigma;
        let n_sigma = self.asperity_density * area;
        let coeff = (4.0 / 3.0) * self.e_star * self.asperity_radius.sqrt() * sigma.powf(1.5);
        let steps = 100;
        let d_abs = d * sigma;
        let upper = d_abs + 6.0 * sigma;
        let step = (upper - d_abs) / steps as f64;
        let mut integral = 0.0;
        let inv_sigma = 1.0 / sigma;
        for k in 0..steps {
            let z = d_abs + (k as f64 + 0.5) * step;
            let zeta = (z - d_abs) / sigma;
            let pdf = inv_sigma
                * (1.0 / (2.0 * std::f64::consts::PI).sqrt())
                * (-0.5 * (z * inv_sigma) * (z * inv_sigma)).exp();
            integral += zeta.powf(1.5) * pdf * step;
        }
        n_sigma * coeff * integral
    }
    /// Real contact area at separation d.
    pub fn real_contact_area(&self, nominal_area: f64, d: f64) -> f64 {
        let d_abs = d * self.sigma;
        let p_contact = 0.5 * erfc_approx(d_abs / (self.sigma * std::f64::consts::SQRT_2));
        self.asperity_density
            * std::f64::consts::PI
            * self.asperity_radius
            * nominal_area
            * p_contact
    }
}
/// Polyhedral approximation of a friction cone.
#[derive(Debug, Clone)]
pub struct FrictionCone {
    /// Coefficient of friction.
    pub mu: f64,
    /// Contact normal (unit vector).
    pub normal: [f64; 3],
    /// Number of edges in the polyhedral approximation (4 or 8).
    pub num_edges: usize,
    /// Precomputed cone edge directions.
    pub edges: Vec<[f64; 3]>,
}
impl FrictionCone {
    /// Build a polyhedral friction cone with `num_edges` edges.
    ///
    /// `num_edges` should be 4 or 8.
    pub fn new(mu: f64, normal: [f64; 3], num_edges: usize) -> Self {
        let n = normalize3(normal);
        let edges = compute_cone_edges(n, mu, num_edges);
        Self {
            mu,
            normal: n,
            num_edges,
            edges,
        }
    }
    /// Check whether a friction force vector lies inside the cone.
    pub fn in_cone(&self, f: [f64; 3]) -> bool {
        let ft = dot3(self.normal, f);
        let fn_ = [
            f[0] - ft * self.normal[0],
            f[1] - ft * self.normal[1],
            f[2] - ft * self.normal[2],
        ];
        let ft_mag = norm3(fn_);
        let fn_mag = ft.abs();
        ft_mag <= self.mu * fn_mag + 1e-10
    }
    /// Project force onto the friction cone (linearized).
    pub fn project(&self, f: [f64; 3], fn_mag: f64) -> [f64; 3] {
        project_to_friction_cone(f, self.normal, self.mu, fn_mag)
    }
}
/// Contact damping models.
#[derive(Debug, Clone)]
pub struct ContactDamping {
    /// Hunt-Crossley dissipation coefficient `α` (s/m).
    pub alpha: f64,
    /// Coefficient of restitution (for Hunt-Crossley inversion).
    pub restitution: f64,
}
impl ContactDamping {
    /// Construct with explicit Hunt-Crossley coefficient.
    pub fn new(alpha: f64) -> Self {
        Self {
            alpha,
            restitution: 0.0,
        }
    }
    /// Construct from a coefficient of restitution using the approximate relation
    /// `α ≈ 1.5 (1 - e) / v_impact`.
    pub fn from_restitution(e: f64, v_impact: f64) -> Self {
        let alpha = if v_impact.abs() > 1e-10 {
            1.5 * (1.0 - e) / v_impact.abs()
        } else {
            0.0
        };
        Self {
            alpha,
            restitution: e,
        }
    }
    /// Hunt-Crossley damping force: `F_d = α δ^{n} δ̇` where `n=3/2` for Hertz.
    pub fn hunt_crossley_force(&self, delta: f64, delta_dot: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        self.alpha * delta.powf(1.5) * delta_dot
    }
    /// Total contact force (elastic + damping).
    pub fn total_force(
        &self,
        stiffness: &ContactStiffness,
        r_a: f64,
        r_b: f64,
        delta: f64,
        delta_dot: f64,
    ) -> f64 {
        let f_e = stiffness.hertz_force(r_a, r_b, delta);
        let f_d = self.hunt_crossley_force(delta, delta_dot);
        (f_e + f_d).max(0.0)
    }
}
/// Hunt-Crossley (1975) viscoelastic contact model.
///
/// Extends Hertz theory to include energy dissipation:
/// F = k * δ^n * (1 + α * δ̇)
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HuntCrossleyContact {
    /// Hertz stiffness k (N/m^n).
    pub k: f64,
    /// Hertz exponent n (1.5 for spheres).
    pub n: f64,
    /// Dissipation coefficient α (s/m).
    pub alpha: f64,
    /// Stored contact force from previous step.
    pub prev_force: f64,
}
impl HuntCrossleyContact {
    /// Create with given parameters.
    pub fn new(k: f64, n: f64, alpha: f64) -> Self {
        Self {
            k,
            n,
            alpha,
            prev_force: 0.0,
        }
    }
    /// Compute contact force from interference and approach rate.
    pub fn force(&mut self, delta: f64, delta_dot: f64) -> f64 {
        if delta <= 0.0 {
            self.prev_force = 0.0;
            return 0.0;
        }
        let hertz = self.k * delta.powf(self.n);
        let f = hertz * (1.0 + self.alpha * delta_dot);
        let f_clamped = f.max(0.0);
        self.prev_force = f_clamped;
        f_clamped
    }
    /// Coefficient of restitution from impact velocity.
    ///
    /// Approximate: e ≈ 1 - 1.25 α v_0 for small α.
    pub fn restitution_approx(&self, impact_velocity: f64) -> f64 {
        (1.0 - 1.25 * self.alpha * impact_velocity.abs()).max(0.0)
    }
    /// Impact duration estimate (simplified Hertz).
    pub fn impact_duration(&self, mass: f64, impact_velocity: f64) -> f64 {
        if impact_velocity.abs() < 1e-15 || self.k < 1e-15 {
            return 0.0;
        }
        let m2 = mass * mass;
        let k2 = self.k * self.k;
        let v0 = impact_velocity.abs();
        2.87 * (m2 / (k2 * v0)).powf(0.2)
    }
}
/// DMT adhesion model (Derjaguin-Muller-Toporov, 1975).
///
/// Suitable for stiff, small-radius contacts with short-range adhesion.
/// Contact radius follows Hertz, but pull-off force = 2πWR*.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DmtContact {
    /// Reduced elastic modulus E* (Pa).
    pub e_star: f64,
    /// Reduced radius R* (m).
    pub r_star: f64,
    /// Work of adhesion W (J/m²).
    pub work_of_adhesion: f64,
}
impl DmtContact {
    /// Create DMT contact model.
    pub fn new(e_star: f64, r_star: f64, work_of_adhesion: f64) -> Self {
        Self {
            e_star,
            r_star,
            work_of_adhesion,
        }
    }
    /// DMT pull-off force.
    pub fn pull_off_force(&self) -> f64 {
        2.0 * std::f64::consts::PI * self.work_of_adhesion * self.r_star
    }
    /// Contact force under applied load F (DMT: Hertz + constant adhesion).
    pub fn contact_force(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return -self.pull_off_force();
        }
        let hertz = (4.0 / 3.0) * self.e_star * self.r_star.sqrt() * delta.powf(1.5);
        hertz - self.pull_off_force()
    }
    /// Contact radius at given interference (DMT = Hertz).
    pub fn contact_radius(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        (self.r_star * delta).sqrt()
    }
    /// Maugis parameter λ (determines JKR/DMT transition).
    ///
    /// λ > 5: JKR regime; λ < 0.1: DMT regime.
    pub fn maugis_parameter(&self, sigma_0: f64) -> f64 {
        if self.work_of_adhesion < 1e-30 || self.e_star < 1e-15 {
            return 0.0;
        }
        let bracket = self.r_star
            / (std::f64::consts::PI * self.work_of_adhesion * self.e_star * self.e_star);
        2.06 * sigma_0 * bracket.powf(1.0 / 3.0)
    }
}
/// Boussinesq elastic half-space solution for a point load.
///
/// Gives displacement and stress fields in a half-space due to
/// a concentrated normal load P at the surface.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BoussinesqHalfSpace {
    /// Young's modulus E (Pa).
    pub e: f64,
    /// Poisson's ratio ν.
    pub nu: f64,
}
impl BoussinesqHalfSpace {
    /// Create Boussinesq model.
    pub fn new(e: f64, nu: f64) -> Self {
        Self { e, nu }
    }
    /// Surface vertical displacement at radial distance r from load P.
    ///
    /// w = (1-ν²)/(π E) * P/r
    pub fn surface_displacement(&self, p: f64, r: f64) -> f64 {
        if r < 1e-30 {
            return f64::INFINITY;
        }
        (1.0 - self.nu * self.nu) / (std::f64::consts::PI * self.e) * p / r
    }
    /// Stress σ_z at point (r, z) in the half-space.
    ///
    /// σ_z = -(3P/2π) z³/R⁵ where R = sqrt(r²+z²).
    pub fn stress_z(&self, p: f64, r: f64, z: f64) -> f64 {
        if z <= 0.0 {
            return 0.0;
        }
        let r2 = r * r;
        let z2 = z * z;
        let r5 = (r2 + z2).powf(2.5);
        if r5 < 1e-60 {
            return 0.0;
        }
        -3.0 * p / (2.0 * std::f64::consts::PI) * z2 * z / r5
    }
    /// Radial stress σ_r at point (r, z).
    pub fn stress_r(&self, p: f64, r: f64, z: f64) -> f64 {
        if r < 1e-30 || z <= 0.0 {
            return 0.0;
        }
        let r2 = r * r;
        let z2 = z * z;
        let rv = (r2 + z2).sqrt();
        let rv3 = rv * rv * rv;
        let rv5 = rv3 * rv * rv;
        let term1 = p / (2.0 * std::f64::consts::PI);
        let term2 = (1.0 - 2.0 * self.nu) * r2 / ((rv) * (rv + z) * (rv + z));
        let term3 = -3.0 * r2 * z / rv5;
        let r_term = (1.0 - 2.0 * self.nu) * (1.0 / rv - z / (rv * (rv + z)));
        term1 * (r_term + term3)
            + term1 * (1.0 - 2.0 * self.nu) / (rv * (rv + z)) * (r2 / (rv + z) - 1.0) * 0.0
            - term1 * term2 * 0.0
    }
    /// Maximum principal stress at given depth z on axis (r=0).
    pub fn max_stress_on_axis(&self, p: f64, z: f64) -> f64 {
        if z <= 0.0 {
            return 0.0;
        }
        -3.0 * p / (2.0 * std::f64::consts::PI * z * z)
    }
}
/// Represents a contact pair between two bodies.
#[derive(Debug, Clone, Copy)]
pub struct ContactPair {
    /// Index of the first body.
    pub body_a: usize,
    /// Index of the second body.
    pub body_b: usize,
    /// Contact point in world space.
    pub point: [f64; 3],
    /// Contact normal pointing from A to B.
    pub normal: [f64; 3],
    /// Penetration depth (positive when overlapping).
    pub depth: f64,
}
impl ContactPair {
    /// Construct a new contact pair.
    pub fn new(
        body_a: usize,
        body_b: usize,
        point: [f64; 3],
        normal: [f64; 3],
        depth: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            point,
            normal,
            depth,
        }
    }
}
/// Elastic contact stiffness models.
#[derive(Debug, Clone)]
pub struct ContactStiffness {
    /// Young's modulus of body A (Pa).
    pub e_a: f64,
    /// Young's modulus of body B (Pa).
    pub e_b: f64,
    /// Poisson's ratio of body A.
    pub nu_a: f64,
    /// Poisson's ratio of body B.
    pub nu_b: f64,
}
impl ContactStiffness {
    /// Construct a new contact stiffness model.
    pub fn new(e_a: f64, e_b: f64, nu_a: f64, nu_b: f64) -> Self {
        Self {
            e_a,
            e_b,
            nu_a,
            nu_b,
        }
    }
    /// Combined elastic modulus for Hertz contact: `E* = [(1-νa²)/Ea + (1-νb²)/Eb]^{-1}`.
    pub fn hertz_modulus(&self) -> f64 {
        let inv =
            (1.0 - self.nu_a * self.nu_a) / self.e_a + (1.0 - self.nu_b * self.nu_b) / self.e_b;
        if inv.abs() < 1e-20 { 0.0 } else { 1.0 / inv }
    }
    /// Hertz contact force for a sphere-sphere contact: `F = (4/3) E* sqrt(R*) δ^{3/2}`.
    ///
    /// `r_a` and `r_b` are sphere radii; `delta` is the penetration depth.
    pub fn hertz_force(&self, r_a: f64, r_b: f64, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let r_eff = r_a * r_b / (r_a + r_b);
        let e_star = self.hertz_modulus();
        (4.0 / 3.0) * e_star * r_eff.sqrt() * delta.powf(1.5)
    }
    /// Linearised constraint-based stiffness for penalty contact.
    ///
    /// Returns a spring constant in N/m.
    pub fn penalty_stiffness(&self, contact_area: f64) -> f64 {
        self.hertz_modulus() * contact_area
    }
}
/// Solution to an LCP problem.
#[derive(Debug, Clone)]
pub struct LcpSolution {
    /// Solution vector `z`.
    pub z: Vec<f64>,
    /// Complementary vector `w = M z + q`.
    pub w: Vec<f64>,
    /// Whether the solver converged.
    pub converged: bool,
    /// Number of pivots performed.
    pub pivots: usize,
}
/// Linear Complementarity Problem formulation for contact dynamics.
///
/// The LCP is: `w = M z + q`, subject to `w ≥ 0`, `z ≥ 0`, `w · z = 0`.
#[derive(Debug, Clone)]
pub struct LcpFormulation {
    /// Symmetric positive semi-definite matrix `M` (n×n, row-major).
    pub m: Vec<f64>,
    /// Right-hand side vector `q` (length n).
    pub q: Vec<f64>,
    /// Problem dimension.
    pub n: usize,
}
impl LcpFormulation {
    /// Create an empty LCP of dimension `n`.
    pub fn new(n: usize) -> Self {
        Self {
            m: vec![0.0; n * n],
            q: vec![0.0; n],
            n,
        }
    }
    /// Build the LCP matrix from contact constraints.
    ///
    /// This is a simplified builder: sets M\[i\]\[i\] = effective_mass_i and
    /// q\[i\] = velocity_bias_i.
    pub fn from_contacts(eff_masses: &[f64], velocity_biases: &[f64]) -> Self {
        let n = eff_masses.len();
        let mut lcp = Self::new(n);
        for i in 0..n {
            lcp.m[i * n + i] = eff_masses[i];
            lcp.q[i] = velocity_biases[i];
        }
        lcp
    }
    /// Validate that M is square and dimensions match.
    pub fn is_valid(&self) -> bool {
        self.m.len() == self.n * self.n && self.q.len() == self.n
    }
}
/// Contact state for impulse-based resolution.
#[derive(Debug, Clone)]
pub struct ContactState {
    /// Contact pair data.
    pub pair: ContactPair,
    /// Accumulated normal impulse.
    pub lambda_n: f64,
    /// Accumulated tangential impulse.
    pub lambda_t: [f64; 2],
    /// Coefficient of restitution.
    pub restitution: f64,
    /// Friction coefficient.
    pub mu: f64,
}
impl ContactState {
    /// Create a new contact state.
    pub fn new(pair: ContactPair, restitution: f64, mu: f64) -> Self {
        Self {
            pair,
            lambda_n: 0.0,
            lambda_t: [0.0; 2],
            restitution,
            mu,
        }
    }
}
/// Methods for resolving penetration in physics simulations.
#[derive(Debug, Clone)]
pub struct PenetrationHandling {
    /// Baumgarte stabilization factor `β`.
    pub beta: f64,
    /// Allowed penetration slop.
    pub slop: f64,
    /// Penalty spring stiffness.
    pub penalty_k: f64,
}
impl PenetrationHandling {
    /// Create with given Baumgarte factor and slop.
    pub fn new(beta: f64, slop: f64, penalty_k: f64) -> Self {
        Self {
            beta,
            slop,
            penalty_k,
        }
    }
    /// Baumgarte stabilization correction velocity.
    ///
    /// `v_bias = -(β / dt) * max(depth - slop, 0)`.
    pub fn baumgarte_velocity(&self, depth: f64, dt: f64) -> f64 {
        baumgarte_correction(depth, self.beta, dt, self.slop)
    }
    /// Penalty force for penetration.
    pub fn penalty_force(&self, depth: f64) -> f64 {
        if depth > 0.0 {
            self.penalty_k * depth
        } else {
            0.0
        }
    }
    /// Position-based correction: move bodies apart by `depth - slop` along normal.
    pub fn position_correction(&self, depth: f64) -> f64 {
        (depth - self.slop).max(0.0)
    }
    /// Project a position out of penetration.
    pub fn project_position(&self, pos: [f64; 3], normal: [f64; 3], depth: f64) -> [f64; 3] {
        let corr = self.position_correction(depth);
        [
            pos[0] + normal[0] * corr,
            pos[1] + normal[1] * corr,
            pos[2] + normal[2] * corr,
        ]
    }
}
/// JKR (Johnson-Kendall-Roberts, 1971) adhesive contact model.
///
/// Extends Hertz theory with surface adhesion via work of adhesion W.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct JkrContact {
    /// Reduced elastic modulus E* (Pa).
    pub e_star: f64,
    /// Reduced radius R* (m).
    pub r_star: f64,
    /// Thermodynamic work of adhesion W (J/m²).
    pub work_of_adhesion: f64,
}
impl JkrContact {
    /// Create JKR contact model.
    pub fn new(e_star: f64, r_star: f64, work_of_adhesion: f64) -> Self {
        Self {
            e_star,
            r_star,
            work_of_adhesion,
        }
    }
    /// JKR pull-off (critical) force.
    pub fn pull_off_force(&self) -> f64 {
        1.5 * std::f64::consts::PI * self.work_of_adhesion * self.r_star
    }
    /// Contact radius under applied load F (JKR theory).
    ///
    /// Solves the JKR cubic: a³ = (R*/E*)\[F + 3πWR* + sqrt(6πWR*F + (3πWR*)²)\]
    pub fn contact_radius(&self, f_applied: f64) -> f64 {
        let w = self.work_of_adhesion;
        let r = self.r_star;
        let e = self.e_star;
        let pi = std::f64::consts::PI;
        let f_ad = 3.0 * pi * w * r;
        let discriminant = 6.0 * pi * w * r * f_applied + f_ad * f_ad;
        if discriminant < 0.0 {
            return 0.0;
        }
        let inner = f_applied + f_ad + discriminant.sqrt();
        if inner < 0.0 {
            return 0.0;
        }
        (r / e * inner).powf(1.0 / 3.0)
    }
    /// Surface separation at contact boundary (JKR).
    ///
    /// Returns the neck height δ_neck.
    pub fn neck_height(&self, a: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        let hertz_delta = a * a / self.r_star;
        let adhesion_term =
            (4.0 * std::f64::consts::PI * self.work_of_adhesion * a / (3.0 * self.e_star)).sqrt();
        hertz_delta - adhesion_term
    }
    /// Total force including adhesion at contact radius a.
    pub fn total_force(&self, a: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        let hertz = (4.0 / 3.0) * self.e_star * a * a * a / self.r_star;
        let adhesion =
            (8.0 * std::f64::consts::PI * self.work_of_adhesion * self.e_star * a * a * a).sqrt();
        hertz - adhesion
    }
}
/// Hertz contact pressure distribution over the contact patch.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct HertzPressureDistribution {
    /// Contact radius a (m).
    pub contact_radius: f64,
    /// Maximum pressure p₀ (Pa).
    pub p0: f64,
    /// Reduced elastic modulus E* (Pa).
    pub e_star: f64,
    /// Reduced radius R* (m).
    pub r_star: f64,
}
impl HertzPressureDistribution {
    /// Create from applied normal force F.
    pub fn from_force(f: f64, e_star: f64, r_star: f64) -> Self {
        if f <= 0.0 {
            return Self {
                contact_radius: 0.0,
                p0: 0.0,
                e_star,
                r_star,
            };
        }
        let a = ((3.0 * f * r_star) / (4.0 * e_star)).powf(1.0 / 3.0);
        let p0 = 3.0 * f / (2.0 * std::f64::consts::PI * a * a);
        Self {
            contact_radius: a,
            p0,
            e_star,
            r_star,
        }
    }
    /// Pressure at radial position r inside the contact patch.
    pub fn pressure(&self, r: f64) -> f64 {
        if r >= self.contact_radius || self.contact_radius <= 1e-30 {
            return 0.0;
        }
        let x = r / self.contact_radius;
        self.p0 * (1.0 - x * x).sqrt()
    }
    /// Average pressure over contact patch.
    pub fn average_pressure(&self) -> f64 {
        2.0 * self.p0 / 3.0
    }
    /// Peak surface stress τ_max (occurs at r ≈ 0.48a, z ≈ 0.48a below surface).
    pub fn max_shear_stress(&self) -> f64 {
        0.31 * self.p0
    }
    /// Depth of maximum shear stress below surface.
    pub fn max_shear_depth(&self) -> f64 {
        0.48 * self.contact_radius
    }
    /// Integration of pressure: total force (numerical check).
    pub fn integrated_force(&self) -> f64 {
        let a = self.contact_radius;
        let steps = 1000;
        let dr = a / steps as f64;
        let mut total = 0.0;
        for k in 0..steps {
            let r = (k as f64 + 0.5) * dr;
            total += self.pressure(r) * 2.0 * std::f64::consts::PI * r * dr;
        }
        total
    }
}
/// Impulse-based contact resolver.
#[derive(Debug, Clone)]
pub struct ImpulseBasedSolver {
    /// Number of solver iterations.
    pub iterations: usize,
    /// Whether to use sequential (true) or simultaneous (false) resolution.
    pub sequential: bool,
}
impl ImpulseBasedSolver {
    /// Create a new solver.
    pub fn new(iterations: usize, sequential: bool) -> Self {
        Self {
            iterations,
            sequential,
        }
    }
    /// Compute the normal impulse for a contact.
    ///
    /// Returns the impulse magnitude `j`.
    pub fn normal_impulse(
        &self,
        v_rel_n: f64,
        eff_mass: f64,
        restitution: f64,
        lambda_n: &mut f64,
    ) -> f64 {
        let target_vn = -restitution * v_rel_n.min(0.0);
        let raw = eff_mass * (target_vn - v_rel_n);
        let old = *lambda_n;
        *lambda_n = (old + raw).max(0.0);
        *lambda_n - old
    }
    /// Compute friction impulse clamped to the friction cone.
    pub fn friction_impulse(
        &self,
        v_rel_t: [f64; 2],
        eff_mass_t: f64,
        mu: f64,
        lambda_n: f64,
        lambda_t: &mut [f64; 2],
    ) -> [f64; 2] {
        let max_impulse = mu * lambda_n;
        let raw = [-eff_mass_t * v_rel_t[0], -eff_mass_t * v_rel_t[1]];
        let new_lt = [lambda_t[0] + raw[0], lambda_t[1] + raw[1]];
        let lt_mag = (new_lt[0] * new_lt[0] + new_lt[1] * new_lt[1]).sqrt();
        let clamped = if lt_mag > max_impulse && lt_mag > 1e-15 {
            [
                new_lt[0] * max_impulse / lt_mag,
                new_lt[1] * max_impulse / lt_mag,
            ]
        } else {
            new_lt
        };
        let delta = [clamped[0] - lambda_t[0], clamped[1] - lambda_t[1]];
        *lambda_t = clamped;
        delta
    }
}
/// Chang-Etsion-Bogy (1987) single asperity elastic-plastic contact.
///
/// Combines Hertz theory (elastic regime) with fully plastic regime.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CebAsperity {
    /// Asperity tip radius (m).
    pub radius: f64,
    /// Reduced elastic modulus E* (Pa).
    pub e_star: f64,
    /// Material yield strength (Pa).
    pub yield_strength: f64,
}
impl CebAsperity {
    /// Create asperity model.
    pub fn new(radius: f64, e_star: f64, yield_strength: f64) -> Self {
        Self {
            radius,
            e_star,
            yield_strength,
        }
    }
    /// Critical interference at yield onset (Hertz elastic limit).
    pub fn critical_interference(&self) -> f64 {
        let c = 0.6;
        let ratio = std::f64::consts::PI * c * self.yield_strength / (2.0 * self.e_star);
        ratio * ratio * self.radius
    }
    /// Contact force at given interference δ.
    pub fn force(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let dc = self.critical_interference();
        if delta <= dc {
            (4.0 / 3.0) * self.e_star * self.radius.sqrt() * delta.powf(1.5)
        } else {
            let ratio = delta / dc;
            let f_c = (4.0 / 3.0) * self.e_star * self.radius.sqrt() * dc.powf(1.5);
            let m = if ratio < 6.0 { 1.5 } else { 2.0 };
            f_c * ratio.powf(m)
        }
    }
    /// Real contact area at given interference.
    pub fn contact_area(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        std::f64::consts::PI * self.radius * delta
    }
}
/// Mindlin (1949) tangential contact model.
///
/// Computes tangential stiffness and slip for a Hertz contact under
/// combined normal and tangential loading.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MindlinContact {
    /// Hertz contact radius a (m).
    pub contact_radius: f64,
    /// Reduced shear modulus G* (Pa).
    pub g_star: f64,
    /// Normal force N (N).
    pub normal_force: f64,
    /// Coulomb friction coefficient μ.
    pub mu: f64,
    /// Current tangential displacement δ_t (m).
    pub tangential_displacement: f64,
    /// Stored tangential force Q (N).
    pub tangential_force: f64,
}
impl MindlinContact {
    /// Create Mindlin contact from Hertz parameters.
    pub fn new(contact_radius: f64, g_star: f64, normal_force: f64, mu: f64) -> Self {
        Self {
            contact_radius,
            g_star,
            normal_force,
            mu,
            tangential_displacement: 0.0,
            tangential_force: 0.0,
        }
    }
    /// Mindlin tangential stiffness K_t = 8 G* a.
    pub fn tangential_stiffness(&self) -> f64 {
        8.0 * self.g_star * self.contact_radius
    }
    /// Incremental tangential force for a displacement increment Δδ_t.
    ///
    /// Implements stick-slip: if |Q + K_t Δδ| > μN, slip occurs.
    pub fn increment_force(&mut self, delta_t_increment: f64) -> f64 {
        let kt = self.tangential_stiffness();
        let q_trial = self.tangential_force + kt * delta_t_increment;
        let q_max = self.mu * self.normal_force;
        if q_trial.abs() <= q_max {
            self.tangential_force = q_trial;
            self.tangential_displacement += delta_t_increment;
        } else {
            self.tangential_force = q_max * q_trial.signum();
            self.tangential_displacement +=
                delta_t_increment - (q_trial - self.tangential_force) / kt;
        }
        self.tangential_force
    }
    /// Remaining stick capacity (before slip initiates).
    pub fn stick_capacity(&self) -> f64 {
        (self.mu * self.normal_force - self.tangential_force.abs()).max(0.0)
    }
}
