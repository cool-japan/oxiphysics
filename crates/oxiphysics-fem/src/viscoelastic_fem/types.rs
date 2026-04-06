//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
/// One term of a Prony series: modulus E_i and relaxation time tau_i.
#[derive(Debug, Clone)]
pub struct PronyTerm {
    /// Partial modulus for this Maxwell arm (Pa or normalised).
    pub modulus: f64,
    /// Relaxation time (seconds).
    pub relaxation_time: f64,
}
/// Internal variables h_i associated with each Prony arm at one integration point.
#[derive(Debug, Clone)]
pub struct InternalVariable {
    /// Current value of each h_i.
    pub values: Vec<f64>,
    /// Number of Prony terms.
    pub n_terms: usize,
}
impl InternalVariable {
    /// Create zeroed internal variables for `n_terms` Prony arms.
    pub fn new(n_terms: usize) -> Self {
        Self {
            values: vec![0.0; n_terms],
            n_terms,
        }
    }
    /// Exact-integration update over time increment `dt` for a strain increment `delta_eps`.
    pub fn update(&mut self, strain_increment: f64, material: &ViscoelasticMaterial, dt: f64) {
        for (i, term) in material.prony_terms.iter().enumerate() {
            let xi = (-dt / term.relaxation_time).exp();
            let coeff = if dt.abs() < 1e-15 {
                term.modulus
            } else {
                term.modulus * (1.0 - xi) * term.relaxation_time / dt
            };
            self.values[i] = self.values[i] * xi + coeff * strain_increment;
        }
    }
    /// Update with time-temperature superposition shift factor.
    pub fn update_with_tts(
        &mut self,
        strain_increment: f64,
        material: &ViscoelasticMaterial,
        dt: f64,
        shift_factor: f64,
    ) {
        let dt_reduced = if shift_factor.abs() > 1e-30 {
            dt / shift_factor
        } else {
            dt
        };
        self.update(strain_increment, material, dt_reduced);
    }
    /// Total internal stress contribution: sum of all h_i.
    pub fn total_stress(&self) -> f64 {
        self.values.iter().sum()
    }
    /// Reset all internal variables to zero.
    pub fn reset(&mut self) {
        for v in &mut self.values {
            *v = 0.0;
        }
    }
}
/// Generalized Maxwell model in shear modulus form.
///
/// G(t) = G_inf + sum_i G_i * exp(-t / tau_i)
#[derive(Debug, Clone)]
pub struct ShearProny {
    /// Long-term shear modulus G_inf (Pa).
    pub g_inf: f64,
    /// Prony terms: (G_i, tau_i).
    pub terms: Vec<(f64, f64)>,
}
impl ShearProny {
    /// Create from G_inf and (G_i, tau_i) pairs.
    pub fn new(g_inf: f64, terms: Vec<(f64, f64)>) -> Self {
        Self { g_inf, terms }
    }
    /// Typical rubber (3-term Prony series for shear).
    pub fn rubber() -> Self {
        Self {
            g_inf: 0.5e6,
            terms: vec![(2.0e6, 0.01), (1.0e6, 0.1), (0.5e6, 1.0)],
        }
    }
    /// Shear relaxation modulus G(t).
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        self.g_inf
            + self
                .terms
                .iter()
                .map(|(g_i, tau_i)| g_i * (-t / tau_i).exp())
                .sum::<f64>()
    }
    /// Instantaneous shear modulus G_0 = G_inf + sum G_i.
    pub fn instantaneous_modulus(&self) -> f64 {
        self.g_inf + self.terms.iter().map(|(g_i, _)| g_i).sum::<f64>()
    }
    /// Storage shear modulus G'(omega).
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        self.g_inf
            + self
                .terms
                .iter()
                .map(|(g_i, tau_i)| {
                    let wt = omega * tau_i;
                    g_i * wt * wt / (1.0 + wt * wt)
                })
                .sum::<f64>()
    }
    /// Loss shear modulus G''(omega).
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        self.terms
            .iter()
            .map(|(g_i, tau_i)| {
                let wt = omega * tau_i;
                g_i * wt / (1.0 + wt * wt)
            })
            .sum::<f64>()
    }
    /// Loss tangent tan(delta) = G'' / G'.
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        let g_prime = self.storage_modulus(omega);
        if g_prime.abs() < 1e-30 {
            return 0.0;
        }
        self.loss_modulus(omega) / g_prime
    }
    /// Peak loss tangent frequency (approximate — scan over decades).
    pub fn peak_loss_frequency(&self) -> f64 {
        if self.terms.is_empty() {
            return 0.0;
        }
        let dominant = self
            .terms
            .iter()
            .max_by(|(g1, _), (g2, _)| g1.partial_cmp(g2).unwrap_or(std::cmp::Ordering::Equal));
        if let Some((_, tau)) = dominant {
            if tau.abs() > 1e-30 { 1.0 / tau } else { 0.0 }
        } else {
            0.0
        }
    }
}
/// Generalised Maxwell (Prony series) viscoelastic material.
///
/// The relaxation modulus is
/// ```text
/// E(t) = E_inf + sum_i E_i * exp(-t / tau_i)
/// ```
#[derive(Debug, Clone)]
pub struct ViscoelasticMaterial {
    /// Long-term (equilibrium) modulus.
    pub e_inf: f64,
    /// Prony series terms.
    pub prony_terms: Vec<PronyTerm>,
    /// Poisson ratio (used for extension to 3-D; not used in 1-D scalar FEM here).
    pub nu: f64,
}
impl ViscoelasticMaterial {
    /// Typical polymer with three Prony arms.
    pub fn polymer() -> Self {
        Self {
            e_inf: 0.5e9,
            prony_terms: vec![
                PronyTerm {
                    modulus: 1.5e9,
                    relaxation_time: 1.0,
                },
                PronyTerm {
                    modulus: 0.8e9,
                    relaxation_time: 10.0,
                },
                PronyTerm {
                    modulus: 0.4e9,
                    relaxation_time: 100.0,
                },
            ],
            nu: 0.35,
        }
    }
    /// Rubber-like material with broader spectrum.
    pub fn rubber() -> Self {
        Self {
            e_inf: 1.0e6,
            prony_terms: vec![
                PronyTerm {
                    modulus: 5.0e6,
                    relaxation_time: 0.01,
                },
                PronyTerm {
                    modulus: 3.0e6,
                    relaxation_time: 0.1,
                },
                PronyTerm {
                    modulus: 1.0e6,
                    relaxation_time: 1.0,
                },
                PronyTerm {
                    modulus: 0.5e6,
                    relaxation_time: 10.0,
                },
            ],
            nu: 0.49,
        }
    }
    /// Instantaneous (glassy) modulus E_0 = E_inf + sum(E_i).
    pub fn instantaneous_modulus(&self) -> f64 {
        self.e_inf + self.prony_terms.iter().map(|t| t.modulus).sum::<f64>()
    }
    /// Relaxation modulus E(t) = E_inf + sum E_i * exp(-t/tau_i).
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        self.e_inf
            + self
                .prony_terms
                .iter()
                .map(|term| term.modulus * (-t / term.relaxation_time).exp())
                .sum::<f64>()
    }
    /// Loss modulus E''(omega) = sum E_i * omega*tau_i / (1 + (omega*tau_i)^2).
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        self.prony_terms
            .iter()
            .map(|term| {
                let wt = omega * term.relaxation_time;
                term.modulus * wt / (1.0 + wt * wt)
            })
            .sum()
    }
    /// Storage modulus E'(omega) = E_inf + sum E_i * (omega*tau_i)^2 / (1 + (omega*tau_i)^2).
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        self.e_inf
            + self
                .prony_terms
                .iter()
                .map(|term| {
                    let wt = omega * term.relaxation_time;
                    term.modulus * wt * wt / (1.0 + wt * wt)
                })
                .sum::<f64>()
    }
    /// Loss tangent tan(delta) = E''(omega) / E'(omega).
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        let e_storage = self.storage_modulus(omega);
        if e_storage.abs() < 1e-30 {
            return 0.0;
        }
        self.loss_modulus(omega) / e_storage
    }
    /// Number of Prony terms.
    pub fn num_terms(&self) -> usize {
        self.prony_terms.len()
    }
    /// Relaxation modulus normalized by instantaneous modulus.
    pub fn normalized_relaxation(&self, t: f64) -> f64 {
        let e0 = self.instantaneous_modulus();
        if e0.abs() < 1e-30 {
            return 1.0;
        }
        self.relaxation_modulus(t) / e0
    }
}
/// Kelvin-Voigt model (spring and dashpot in parallel).
///
/// Constitutive equation: sigma = E * epsilon + eta * d(epsilon)/dt
#[derive(Debug, Clone)]
pub struct KelvinVoigtModel {
    /// Spring modulus (Pa).
    pub modulus: f64,
    /// Dashpot viscosity (Pa·s).
    pub viscosity: f64,
}
impl KelvinVoigtModel {
    /// Create a new Kelvin-Voigt model.
    pub fn new(modulus: f64, viscosity: f64) -> Self {
        Self { modulus, viscosity }
    }
    /// Retardation time tau = eta / E (s).
    pub fn retardation_time(&self) -> f64 {
        if self.modulus.abs() < 1e-30 {
            return f64::INFINITY;
        }
        self.viscosity / self.modulus
    }
    /// Creep compliance J(t) = (1/E) * (1 - exp(-t/tau)).
    pub fn creep_compliance(&self, t: f64) -> f64 {
        if self.modulus.abs() < 1e-30 {
            return 0.0;
        }
        let tau = self.retardation_time();
        if tau.is_infinite() || tau.abs() < 1e-30 {
            return 1.0 / self.modulus;
        }
        (1.0 / self.modulus) * (1.0 - (-t / tau).exp())
    }
    /// Strain response to constant applied stress sigma_0.
    ///
    /// epsilon(t) = (sigma_0 / E) * (1 - exp(-t/tau))
    pub fn creep_strain(&self, sigma_0: f64, t: f64) -> f64 {
        sigma_0 * self.creep_compliance(t)
    }
    /// Relaxation: for KV model, relaxation is not well-defined (infinite relaxation time).
    /// Returns the instantaneous stress for given strain (elastic part only).
    pub fn instantaneous_stress(&self, epsilon: f64) -> f64 {
        self.modulus * epsilon
    }
    /// Storage modulus for Kelvin-Voigt model.
    ///
    /// E'(omega) = E (constant for KV model)
    pub fn storage_modulus(&self, _omega: f64) -> f64 {
        self.modulus
    }
    /// Loss modulus for Kelvin-Voigt model.
    ///
    /// E''(omega) = eta * omega
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        self.viscosity * omega
    }
    /// Loss tangent tan(delta) = E'' / E'.
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        if self.modulus.abs() < 1e-30 {
            return f64::INFINITY;
        }
        self.viscosity * omega / self.modulus
    }
}
/// Standard Linear Solid (SLS / Zener model).
///
/// Combination of a spring E1 in series with a Kelvin-Voigt unit (E2, eta).
///
/// Relaxation modulus:
/// E(t) = E1 * E2 / (E1 + E2) + E1^2 / (E1 + E2) * exp(-t / tau_r)
///
/// where tau_r = eta * (E1 + E2) / (E1 * E2)
#[derive(Debug, Clone)]
pub struct StandardLinearSolid {
    /// Spring modulus E1 (parallel branch, Pa).
    pub e1: f64,
    /// Spring modulus E2 (series spring in Maxwell arm, Pa).
    pub e2: f64,
    /// Dashpot viscosity eta (Maxwell arm, Pa·s).
    pub eta: f64,
}
impl StandardLinearSolid {
    /// Create a new Standard Linear Solid.
    pub fn new(e1: f64, e2: f64, eta: f64) -> Self {
        Self { e1, e2, eta }
    }
    /// Long-term (rubbery) modulus E_inf = E1 * E2 / (E1 + E2).
    pub fn equilibrium_modulus(&self) -> f64 {
        let denom = self.e1 + self.e2;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        self.e1 * self.e2 / denom
    }
    /// Instantaneous (glassy) modulus E_0 = E1.
    pub fn instantaneous_modulus(&self) -> f64 {
        self.e1
    }
    /// Stress relaxation time tau_r = eta * (E1 + E2) / (E1 * E2).
    pub fn relaxation_time(&self) -> f64 {
        let prod = self.e1 * self.e2;
        if prod.abs() < 1e-30 {
            return f64::INFINITY;
        }
        self.eta * (self.e1 + self.e2) / prod
    }
    /// Relaxation modulus E(t).
    ///
    /// E(t) = E_inf + (E_0 - E_inf) * exp(-t / tau_r)
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        let e_inf = self.equilibrium_modulus();
        let e0 = self.instantaneous_modulus();
        let tau_r = self.relaxation_time();
        if tau_r.is_infinite() || tau_r.abs() < 1e-30 {
            return e0;
        }
        e_inf + (e0 - e_inf) * (-t / tau_r).exp()
    }
    /// Creep retardation time tau_c = eta / E2.
    pub fn retardation_time(&self) -> f64 {
        if self.e2.abs() < 1e-30 {
            return f64::INFINITY;
        }
        self.eta / self.e2
    }
    /// Creep compliance J(t) = (1/E_inf) - (1/E_inf - 1/E_0) * exp(-t/tau_c).
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let e_inf = self.equilibrium_modulus();
        let e0 = self.instantaneous_modulus();
        if e_inf.abs() < 1e-30 || e0.abs() < 1e-30 {
            return 0.0;
        }
        let tau_c = self.retardation_time();
        if tau_c.is_infinite() || tau_c.abs() < 1e-30 {
            return 1.0 / e_inf;
        }
        1.0 / e_inf - (1.0 / e_inf - 1.0 / e0) * (-t / tau_c).exp()
    }
    /// Storage modulus E'(omega) for SLS.
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let e_inf = self.equilibrium_modulus();
        let e0 = self.instantaneous_modulus();
        let tau_r = self.relaxation_time();
        if tau_r.is_infinite() || tau_r.abs() < 1e-30 {
            return e_inf;
        }
        let wt2 = (omega * tau_r).powi(2);
        e_inf + (e0 - e_inf) * wt2 / (1.0 + wt2)
    }
    /// Loss modulus E''(omega) for SLS.
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let e0 = self.instantaneous_modulus();
        let e_inf = self.equilibrium_modulus();
        let tau_r = self.relaxation_time();
        if tau_r.is_infinite() || tau_r.abs() < 1e-30 {
            return 0.0;
        }
        let wt = omega * tau_r;
        (e0 - e_inf) * wt / (1.0 + wt * wt)
    }
    /// Loss tangent tan(delta) = E'' / E'.
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        let e_prime = self.storage_modulus(omega);
        if e_prime.abs() < 1e-30 {
            return 0.0;
        }
        self.loss_modulus(omega) / e_prime
    }
}
/// State at a single integration point: stress vector and internal variables.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ViscoelasticState {
    /// Stress components (Voigt notation).
    pub stress: Vec<f64>,
    /// Internal (hidden) variables; one `Vec`f64` per Prony arm.
    pub internal_vars: Vec<Vec<f64>>,
}
impl ViscoelasticState {
    /// Create zeroed state with `n_stress` stress components and `n_prony` Prony arms,
    /// each arm having `n_stress` components.
    #[allow(dead_code)]
    pub fn new(n_stress: usize, n_prony: usize) -> Self {
        Self {
            stress: vec![0.0; n_stress],
            internal_vars: vec![vec![0.0; n_stress]; n_prony],
        }
    }
}
/// Assembled 1-D viscoelastic beam made of `ViscoelasticElement1D` elements.
#[derive(Debug)]
pub struct ViscoelasticBeam1D {
    /// Individual elements.
    pub elements: Vec<ViscoelasticElement1D>,
    /// Number of nodes (= n_elements + 1).
    pub n_nodes: usize,
}
impl ViscoelasticBeam1D {
    /// Create a uniform beam of `n_elements` equal-length elements.
    pub fn new(n_elements: usize, length: f64, area: f64, material: ViscoelasticMaterial) -> Self {
        let elem_len = length / n_elements as f64;
        let elements = (0..n_elements)
            .map(|_| ViscoelasticElement1D::new(elem_len, area, material.clone()))
            .collect();
        Self {
            elements,
            n_nodes: n_elements + 1,
        }
    }
    /// Assemble global stiffness matrix (n_nodes × n_nodes).
    pub fn assemble_stiffness(&self, dt: f64) -> Vec<Vec<f64>> {
        let n = self.n_nodes;
        let mut k = vec![vec![0.0; n]; n];
        for (e, elem) in self.elements.iter().enumerate() {
            let ke = elem.stiffness_matrix(dt);
            for i in 0..2 {
                for j in 0..2 {
                    k[e + i][e + j] += ke[i][j];
                }
            }
        }
        k
    }
    /// Assemble internal force vector.
    pub fn assemble_internal_forces(&self) -> Vec<f64> {
        let n = self.n_nodes;
        let mut f = vec![0.0; n];
        for (e, elem) in self.elements.iter().enumerate() {
            let fe = elem.internal_force_vector();
            f[e] += fe[0];
            f[e + 1] += fe[1];
        }
        f
    }
    /// Apply Dirichlet boundary conditions by the penalty method.
    pub fn apply_dirichlet(&self, k: &mut Vec<Vec<f64>>, f: &mut Vec<f64>, dofs: &[(usize, f64)]) {
        let penalty = {
            let max_k = k
                .iter()
                .flat_map(|row| row.iter())
                .cloned()
                .fold(0.0_f64, f64::max);
            max_k * 1.0e14
        };
        for &(dof, val) in dofs {
            k[dof][dof] += penalty;
            f[dof] += penalty * val;
        }
    }
    /// Solve one time step.
    pub fn solve_step(&mut self, forces: &[f64], dt: f64) -> Vec<f64> {
        let mut k = self.assemble_stiffness(dt);
        let mut f = forces.to_vec();
        self.apply_dirichlet(&mut k, &mut f, &[(0, 0.0)]);
        let u = solve_tridiagonal_system(&k, &f);
        for (e, elem) in self.elements.iter_mut().enumerate() {
            let du = u[e + 1] - u[e];
            let d_eps = du / elem.length;
            elem.stress_update(d_eps, dt);
        }
        u
    }
    /// Solve one time step with temperature-dependent properties.
    pub fn solve_step_with_temperature(
        &mut self,
        forces: &[f64],
        dt: f64,
        wlf: &WlfParameters,
        temperature: f64,
    ) -> Vec<f64> {
        let shift = wlf.shift_factor(temperature);
        let dt_reduced = if shift.abs() > 1e-30 { dt / shift } else { dt };
        let mut k = self.assemble_stiffness(dt_reduced);
        let mut f = forces.to_vec();
        self.apply_dirichlet(&mut k, &mut f, &[(0, 0.0)]);
        let u = solve_tridiagonal_system(&k, &f);
        for (e, elem) in self.elements.iter_mut().enumerate() {
            let du = u[e + 1] - u[e];
            let d_eps = du / elem.length;
            elem.stress_update_with_tts(d_eps, dt, shift);
        }
        u
    }
    /// Total strain energy in the beam.
    pub fn total_strain_energy(&self) -> f64 {
        self.elements.iter().map(|e| e.strain_energy()).sum()
    }
    /// Maximum absolute stress in any element.
    pub fn max_stress(&self) -> f64 {
        self.elements
            .iter()
            .map(|e| e.stress().abs())
            .fold(0.0f64, f64::max)
    }
}
/// Single 1-D bar/truss element with viscoelastic constitutive response.
#[derive(Debug, Clone)]
pub struct ViscoelasticElement1D {
    /// Element length (m).
    pub length: f64,
    /// Cross-sectional area (m²).
    pub area: f64,
    /// Material definition.
    pub material: ViscoelasticMaterial,
    /// Internal variables at the single integration point.
    pub internal_vars: InternalVariable,
    /// Current total strain.
    pub(super) strain: f64,
    /// Current Cauchy stress.
    pub(super) stress: f64,
}
impl ViscoelasticElement1D {
    /// Create a new element.
    pub fn new(length: f64, area: f64, material: ViscoelasticMaterial) -> Self {
        let n = material.prony_terms.len();
        Self {
            length,
            area,
            material,
            internal_vars: InternalVariable::new(n),
            strain: 0.0,
            stress: 0.0,
        }
    }
    /// Algorithmic tangent modulus for the current time step size `dt`.
    pub fn algorithmic_modulus(&self, dt: f64) -> f64 {
        let sum: f64 = self
            .material
            .prony_terms
            .iter()
            .map(|term| {
                if dt.abs() < 1e-15 {
                    term.modulus
                } else {
                    let xi = (-dt / term.relaxation_time).exp();
                    term.modulus * (1.0 - xi) * term.relaxation_time / dt
                }
            })
            .sum();
        self.material.e_inf + sum
    }
    /// Element stiffness matrix (2×2) based on algorithmic tangent.
    pub fn stiffness_matrix(&self, dt: f64) -> [[f64; 2]; 2] {
        let k = self.algorithmic_modulus(dt) * self.area / self.length;
        [[k, -k], [-k, k]]
    }
    /// Update stress for a given strain increment and return new stress.
    pub fn stress_update(&mut self, strain_increment: f64, dt: f64) -> f64 {
        self.internal_vars
            .update(strain_increment, &self.material, dt);
        self.strain += strain_increment;
        self.stress =
            self.material.e_inf * self.strain + self.internal_vars.values.iter().sum::<f64>();
        self.stress
    }
    /// Update stress with time-temperature superposition.
    pub fn stress_update_with_tts(
        &mut self,
        strain_increment: f64,
        dt: f64,
        shift_factor: f64,
    ) -> f64 {
        self.internal_vars
            .update_with_tts(strain_increment, &self.material, dt, shift_factor);
        self.strain += strain_increment;
        self.stress =
            self.material.e_inf * self.strain + self.internal_vars.values.iter().sum::<f64>();
        self.stress
    }
    /// Current strain.
    pub fn strain(&self) -> f64 {
        self.strain
    }
    /// Current stress.
    pub fn stress(&self) -> f64 {
        self.stress
    }
    /// Elastic strain energy in the element.
    pub fn strain_energy(&self) -> f64 {
        0.5 * self.stress * self.strain * self.area * self.length
    }
    /// Internal (element) force vector \[f1, f2\].
    pub fn internal_force_vector(&self) -> [f64; 2] {
        let f = self.stress * self.area;
        [-f, f]
    }
}
/// WLF (Williams-Landel-Ferry) shift factor parameters.
///
/// log(a_T) = −C1 (T − T_ref) / (C2 + (T − T_ref))
#[allow(dead_code)]
#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct WLFShift {
    /// WLF constant C1.
    pub C1: f64,
    /// WLF constant C2.
    pub C2: f64,
    /// Reference temperature T_ref.
    pub T_ref: f64,
}
impl WLFShift {
    /// log10 shift factor: log10(a_T) = −C1 (T − T_ref) / (C2 + (T − T_ref)).
    #[allow(dead_code)]
    #[allow(non_snake_case)]
    pub fn shift_factor(&self, T: f64) -> f64 {
        let dt = T - self.T_ref;
        let denom = self.C2 + dt;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        -self.C1 * dt / denom
    }
}
/// Arrhenius time-temperature superposition parameters.
#[derive(Debug, Clone)]
pub struct ArrheniusParameters {
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
    /// Activation energy / gas constant (K).
    pub delta_h_over_r: f64,
}
impl ArrheniusParameters {
    /// Compute the shift factor a_T using Arrhenius equation.
    pub fn shift_factor(&self, temperature: f64) -> f64 {
        let exponent = self.delta_h_over_r * (1.0 / temperature - 1.0 / self.t_ref);
        exponent.exp()
    }
}
/// Prony series relaxation function: E(t) = E_inf + Σ E_i exp(-t/tau_i).
#[allow(dead_code)]
#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct PronyRelaxation {
    /// Long-term equilibrium modulus.
    pub E_inf: f64,
    /// Prony series terms.
    pub elements: Vec<PronyElement>,
}
impl PronyRelaxation {
    /// Relaxation modulus E(t) = E_inf + Σ E_i exp(-t/tau_i).
    #[allow(dead_code)]
    #[allow(non_snake_case)]
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        self.E_inf
            + self
                .elements
                .iter()
                .map(|e| e.E_i * (-t / e.tau_i).exp())
                .sum::<f64>()
    }
}
/// Maxwell viscoelastic model (spring + dashpot in series).
///
/// Stress relaxation:
/// sigma(t) = E * epsilon_0 * exp(-t * E / eta)
///          = sigma_0 * exp(-t / tau)
/// where tau = eta / E.
#[derive(Debug, Clone)]
pub struct MaxwellModel {
    /// Spring stiffness (Pa).
    pub modulus: f64,
    /// Dashpot viscosity (Pa·s).
    pub viscosity: f64,
}
impl MaxwellModel {
    /// Create a new Maxwell model.
    pub fn new(modulus: f64, viscosity: f64) -> Self {
        Self { modulus, viscosity }
    }
    /// Relaxation time tau = eta / E (s).
    pub fn relaxation_time(&self) -> f64 {
        if self.modulus.abs() < 1e-30 {
            return f64::INFINITY;
        }
        self.viscosity / self.modulus
    }
    /// Stress relaxation at time t for initial strain epsilon_0.
    ///
    /// sigma(t) = E * epsilon_0 * exp(-t / tau)
    pub fn stress_relaxation(&self, epsilon_0: f64, t: f64) -> f64 {
        let tau = self.relaxation_time();
        if tau.is_infinite() || tau.abs() < 1e-30 {
            return 0.0;
        }
        self.modulus * epsilon_0 * (-t / tau).exp()
    }
    /// Relaxation modulus E(t) = E * exp(-t / tau).
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        let tau = self.relaxation_time();
        if tau.is_infinite() || tau.abs() < 1e-30 {
            return 0.0;
        }
        self.modulus * (-t / tau).exp()
    }
    /// Creep compliance J(t) = t / eta + 1 / E.
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let j0 = if self.modulus.abs() > 1e-30 {
            1.0 / self.modulus
        } else {
            0.0
        };
        let j_creep = if self.viscosity.abs() > 1e-30 {
            t / self.viscosity
        } else {
            0.0
        };
        j0 + j_creep
    }
    /// Storage modulus E'(omega) for Maxwell model.
    ///
    /// E'(omega) = E * (omega*tau)^2 / (1 + (omega*tau)^2)
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let tau = self.relaxation_time();
        if tau.is_infinite() || tau.abs() < 1e-30 {
            return 0.0;
        }
        let wt = omega * tau;
        self.modulus * wt * wt / (1.0 + wt * wt)
    }
    /// Loss modulus E''(omega) for Maxwell model.
    ///
    /// E''(omega) = E * omega*tau / (1 + (omega*tau)^2)
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let tau = self.relaxation_time();
        if tau.is_infinite() || tau.abs() < 1e-30 {
            return 0.0;
        }
        let wt = omega * tau;
        self.modulus * wt / (1.0 + wt * wt)
    }
    /// Loss tangent tan(delta) = E'' / E'.
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        let e_prime = self.storage_modulus(omega);
        if e_prime.abs() < 1e-30 {
            return f64::INFINITY;
        }
        self.loss_modulus(omega) / e_prime
    }
}
/// A single Prony series term: modulus E_i and relaxation time tau_i.
#[allow(dead_code)]
#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct PronyElement {
    /// Partial modulus for this Maxwell arm.
    pub E_i: f64,
    /// Relaxation time tau_i (s).
    pub tau_i: f64,
}
/// A 4-node finite element with a Prony-series viscoelastic constitutive law.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PronyFemElement {
    /// Global node indices for the 4 corner nodes.
    pub nodes: [usize; 4],
    /// Prony series material.
    pub prony: PronyRelaxation,
    /// Integration-point state.
    pub state: ViscoelasticState,
    /// Time-step size (s).
    pub dt: f64,
}
impl PronyFemElement {
    /// Create a new element with zeroed state.
    #[allow(dead_code)]
    pub fn new(nodes: [usize; 4], prony: PronyRelaxation, dt: f64) -> Self {
        let n_prony = prony.elements.len();
        let state = ViscoelasticState::new(1, n_prony);
        Self {
            nodes,
            prony,
            state,
            dt,
        }
    }
    /// Update stresses for the given strain increment vector.
    ///
    /// Uses exact-integration (algorithmic tangent) for the 1-D case per
    /// component.  Returns the updated stress vector.
    #[allow(dead_code)]
    pub fn update_stresses(&mut self, strain_inc: &[f64]) -> Vec<f64> {
        let n_comp = strain_inc.len();
        if self.state.stress.len() != n_comp {
            self.state.stress = vec![0.0; n_comp];
        }
        if self.state.internal_vars.len() != self.prony.elements.len() {
            self.state.internal_vars = vec![vec![0.0; n_comp]; self.prony.elements.len()];
        }
        for k in 0..self.prony.elements.len() {
            if self.state.internal_vars[k].len() != n_comp {
                self.state.internal_vars[k] = vec![0.0; n_comp];
            }
        }
        let dt = self.dt;
        for (k, elem) in self.prony.elements.iter().enumerate() {
            let xi = (-dt / elem.tau_i).exp();
            let coeff = if dt.abs() < 1e-15 {
                elem.E_i
            } else {
                elem.E_i * (1.0 - xi) * elem.tau_i / dt
            };
            for c in 0..n_comp {
                self.state.internal_vars[k][c] =
                    self.state.internal_vars[k][c] * xi + coeff * strain_inc[c];
            }
        }
        for c in 0..n_comp {
            let h_sum: f64 = self.state.internal_vars.iter().map(|iv| iv[c]).sum();
            self.state.stress[c] += self.prony.E_inf * strain_inc[c] + h_sum
                - self.state.internal_vars.iter().map(|iv| iv[c]).sum::<f64>()
                + h_sum;
        }
        for c in 0..n_comp {
            self.state.stress[c] = self.prony.E_inf * strain_inc[c]
                + self.state.internal_vars.iter().map(|iv| iv[c]).sum::<f64>();
        }
        self.state.stress.clone()
    }
}
/// Kelvin-Voigt creep compliance:
/// J(t) = J_0 + J_inf (1 − exp(−t/tau))
#[allow(dead_code)]
#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct Creep {
    /// Instantaneous compliance J_0 = 1/E_0.
    pub J0: f64,
    /// Long-term creep compliance amplitude J_inf.
    pub J_inf: f64,
    /// Retardation time tau.
    pub tau: f64,
}
impl Creep {
    /// Creep compliance J(t) = J_0 + J_inf (1 − exp(−t/tau)).
    #[allow(dead_code)]
    #[allow(non_snake_case)]
    pub fn compliance(&self, t: f64) -> f64 {
        self.J0 + self.J_inf * (1.0 - (-t / self.tau).exp())
    }
}
/// WLF (Williams-Landel-Ferry) time-temperature superposition parameters.
#[derive(Debug, Clone)]
pub struct WlfParameters {
    /// Reference temperature T_ref (K or °C).
    pub t_ref: f64,
    /// WLF constant C1.
    pub c1: f64,
    /// WLF constant C2.
    pub c2: f64,
}
impl WlfParameters {
    /// Create WLF parameters.
    pub fn new(t_ref: f64, c1: f64, c2: f64) -> Self {
        Self { t_ref, c1, c2 }
    }
    /// Typical WLF parameters for amorphous polymers.
    pub fn typical_polymer() -> Self {
        Self {
            t_ref: 25.0,
            c1: 17.44,
            c2: 51.6,
        }
    }
    /// Compute the shift factor log10(a_T) using WLF equation.
    ///
    /// log10(a_T) = -C1 * (T - T_ref) / (C2 + (T - T_ref))
    pub fn log_shift_factor(&self, temperature: f64) -> f64 {
        let dt = temperature - self.t_ref;
        let denom = self.c2 + dt;
        if denom.abs() < 1e-10 {
            return 0.0;
        }
        -self.c1 * dt / denom
    }
    /// Compute the shift factor a_T.
    pub fn shift_factor(&self, temperature: f64) -> f64 {
        10.0_f64.powf(self.log_shift_factor(temperature))
    }
    /// Compute reduced time: t_reduced = t / a_T.
    pub fn reduced_time(&self, real_time: f64, temperature: f64) -> f64 {
        let at = self.shift_factor(temperature);
        if at.abs() < 1e-30 {
            return real_time;
        }
        real_time / at
    }
    /// Compute relaxation modulus at a given temperature and time.
    pub fn relaxation_modulus_at_temp(
        &self,
        material: &ViscoelasticMaterial,
        time: f64,
        temperature: f64,
    ) -> f64 {
        let t_reduced = self.reduced_time(time, temperature);
        material.relaxation_modulus(t_reduced)
    }
}
/// Result of a frequency sweep (complex modulus analysis).
#[derive(Debug, Clone)]
pub struct FrequencySweepPoint {
    /// Angular frequency (rad/s).
    pub omega: f64,
    /// Storage modulus E' (Pa).
    pub storage: f64,
    /// Loss modulus E'' (Pa).
    pub loss: f64,
    /// Loss tangent tan(delta).
    pub tan_delta: f64,
    /// Complex modulus magnitude |E*| (Pa).
    pub magnitude: f64,
}
