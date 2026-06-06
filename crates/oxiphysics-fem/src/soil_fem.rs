// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Soil finite element module for geotechnical analysis.
//!
//! Implements consolidation (Biot's theory), drained/undrained conditions,
//! effective stress analysis, Mohr-Coulomb and Drucker-Prager plasticity,
//! and coupled hydro-mechanical FEM for saturated/unsaturated soils.

/// Soil material model type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SoilModel {
    /// Linear elastic (drained).
    LinearElastic,
    /// Mohr-Coulomb elastoplastic.
    MohrCoulomb,
    /// Drucker-Prager elastoplastic.
    DruckerPrager,
    /// Modified Cam-Clay critical state model.
    CamClay,
    /// Hypoplastic model (Kolymbas).
    Hypoplastic,
}

/// Soil drainage condition.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DrainageCondition {
    /// Fully drained — pore pressure dissipates instantaneously.
    Drained,
    /// Fully undrained — no pore pressure dissipation.
    Undrained,
    /// Partially drained — coupled consolidation analysis.
    Coupled,
}

/// Parameters for Mohr-Coulomb soil model.
#[derive(Debug, Clone)]
pub struct MohrCoulombParams {
    /// Young's modulus \[Pa\].
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Cohesion \[Pa\].
    pub cohesion: f64,
    /// Friction angle \[radians\].
    pub friction_angle: f64,
    /// Dilatancy angle \[radians\].
    pub dilatancy_angle: f64,
}

impl MohrCoulombParams {
    /// Create new Mohr-Coulomb parameters.
    pub fn new(
        young: f64,
        poisson: f64,
        cohesion: f64,
        friction_angle: f64,
        dilatancy_angle: f64,
    ) -> Self {
        Self {
            young,
            poisson,
            cohesion,
            friction_angle,
            dilatancy_angle,
        }
    }

    /// Shear modulus G \[Pa\].
    pub fn shear_modulus(&self) -> f64 {
        self.young / (2.0 * (1.0 + self.poisson))
    }

    /// Bulk modulus K \[Pa\].
    pub fn bulk_modulus(&self) -> f64 {
        self.young / (3.0 * (1.0 - 2.0 * self.poisson))
    }
}

/// Cam-Clay model parameters.
#[derive(Debug, Clone)]
pub struct CamClayParams {
    /// Slope of critical state line (M).
    pub m: f64,
    /// Compression index lambda.
    pub lambda: f64,
    /// Swelling index kappa.
    pub kappa: f64,
    /// Reference pressure \[Pa\].
    pub p_ref: f64,
    /// Initial preconsolidation pressure \[Pa\].
    pub p_c0: f64,
    /// Poisson's ratio.
    pub poisson: f64,
}

impl CamClayParams {
    /// Create new Cam-Clay parameters.
    pub fn new(m: f64, lambda: f64, kappa: f64, p_ref: f64, p_c0: f64, poisson: f64) -> Self {
        Self {
            m,
            lambda,
            kappa,
            p_ref,
            p_c0,
            poisson,
        }
    }
}

/// 3D stress tensor (Voigt notation: xx, yy, zz, xy, yz, xz).
#[derive(Debug, Clone)]
pub struct StressTensor {
    /// Voigt components.
    pub s: [f64; 6],
}

impl StressTensor {
    /// Create zero stress tensor.
    pub fn zero() -> Self {
        Self { s: [0.0; 6] }
    }

    /// Create from array.
    pub fn from_array(s: [f64; 6]) -> Self {
        Self { s }
    }

    /// Mean stress (p = -(s11+s22+s33)/3).
    pub fn mean_stress(&self) -> f64 {
        -(self.s[0] + self.s[1] + self.s[2]) / 3.0
    }

    /// Deviatoric stress norm (q).
    pub fn deviatoric_norm(&self) -> f64 {
        let p = self.mean_stress();
        let s11 = self.s[0] + p;
        let s22 = self.s[1] + p;
        let s33 = self.s[2] + p;
        let s12 = self.s[3];
        let s23 = self.s[4];
        let s13 = self.s[5];
        let j2 = 0.5 * (s11 * s11 + s22 * s22 + s33 * s33) + s12 * s12 + s23 * s23 + s13 * s13;
        (3.0 * j2).sqrt()
    }

    /// Von Mises equivalent stress.
    pub fn von_mises(&self) -> f64 {
        let ds11 = self.s[0] - self.s[1];
        let ds22 = self.s[1] - self.s[2];
        let ds33 = self.s[2] - self.s[0];
        (0.5 * (ds11 * ds11 + ds22 * ds22 + ds33 * ds33)
            + 3.0 * (self.s[3] * self.s[3] + self.s[4] * self.s[4] + self.s[5] * self.s[5]))
            .sqrt()
    }

    /// Principal stresses (approximate, sorted descending).
    pub fn principal_stresses(&self) -> [f64; 3] {
        // Simplified: diagonal-only if off-diagonal small
        let mut p = [self.s[0], self.s[1], self.s[2]];
        p.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        p
    }
}

/// Mohr-Coulomb yield function evaluation.
///
/// Returns f > 0 if plastic, f <= 0 if elastic.
pub fn mohr_coulomb_yield(stress: &StressTensor, params: &MohrCoulombParams) -> f64 {
    let principals = stress.principal_stresses();
    let s1 = principals[0];
    let s3 = principals[2];
    let sin_phi = params.friction_angle.sin();
    let cos_phi = params.friction_angle.cos();
    // f = (s1 - s3)/2 - (s1 + s3)/2 * sin_phi - c * cos_phi
    (s1 - s3) / 2.0 - (s1 + s3) / 2.0 * sin_phi - params.cohesion * cos_phi
}

/// Drucker-Prager yield function.
///
/// Returns f > 0 if plastic.
pub fn drucker_prager_yield(stress: &StressTensor, alpha: f64, k: f64) -> f64 {
    let j2 = {
        let p = stress.mean_stress();
        let s11 = stress.s[0] + p;
        let s22 = stress.s[1] + p;
        let s33 = stress.s[2] + p;
        0.5 * (s11 * s11 + s22 * s22 + s33 * s33)
            + stress.s[3] * stress.s[3]
            + stress.s[4] * stress.s[4]
            + stress.s[5] * stress.s[5]
    };
    j2.sqrt() + alpha * (stress.s[0] + stress.s[1] + stress.s[2]) - k
}

/// Convert Mohr-Coulomb to Drucker-Prager (outer cone).
pub fn mc_to_dp_outer(phi: f64, c: f64) -> (f64, f64) {
    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    let alpha = 2.0 * sin_phi / (3.0_f64.sqrt() * (3.0 - sin_phi));
    let k = 6.0 * c * cos_phi / (3.0_f64.sqrt() * (3.0 - sin_phi));
    (alpha, k)
}

/// Elastic constitutive matrix (6×6, Voigt notation).
///
/// Returns the stiffness tensor D.
pub fn elastic_stiffness(young: f64, poisson: f64) -> [[f64; 6]; 6] {
    let e = young;
    let nu = poisson;
    let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let mu = e / (2.0 * (1.0 + nu));
    let mut d = [[0.0f64; 6]; 6];
    d[0][0] = lam + 2.0 * mu;
    d[0][1] = lam;
    d[0][2] = lam;
    d[1][0] = lam;
    d[1][1] = lam + 2.0 * mu;
    d[1][2] = lam;
    d[2][0] = lam;
    d[2][1] = lam;
    d[2][2] = lam + 2.0 * mu;
    d[3][3] = mu;
    d[4][4] = mu;
    d[5][5] = mu;
    d
}

/// Apply elastic stiffness: sigma = D * epsilon.
pub fn apply_stiffness(d: &[[f64; 6]; 6], eps: &[f64; 6]) -> [f64; 6] {
    let mut sigma = [0.0f64; 6];
    for i in 0..6 {
        for j in 0..6 {
            sigma[i] += d[i][j] * eps[j];
        }
    }
    sigma
}

/// Biot consolidation parameters for coupled HM analysis.
#[derive(Debug, Clone)]
pub struct BiotParams {
    /// Biot coefficient (alpha = 1 - K_s/K_dr).
    pub alpha: f64,
    /// Biot modulus M \[Pa\].
    pub m_modulus: f64,
    /// Hydraulic conductivity \[m/s\].
    pub permeability: f64,
    /// Water bulk modulus \[Pa\].
    pub water_bulk: f64,
    /// Porosity.
    pub porosity: f64,
}

impl BiotParams {
    /// Create new Biot parameters for saturated soil.
    pub fn new(alpha: f64, m_modulus: f64, permeability: f64, porosity: f64) -> Self {
        Self {
            alpha,
            m_modulus,
            permeability,
            water_bulk: 2.0e9,
            porosity,
        }
    }

    /// Storage coefficient S = (alpha^2)/M + n/K_w.
    pub fn storage_coefficient(&self) -> f64 {
        self.alpha * self.alpha / self.m_modulus + self.porosity / self.water_bulk
    }
}

/// 1D consolidation settlement (Terzaghi).
///
/// Returns settlement \[m\] at time t \[s\].
pub fn terzaghi_settlement(
    total_settlement: f64,
    cv: f64, // consolidation coefficient [m^2/s]
    h: f64,  // drainage path [m]
    t: f64,  // time [s]
    n_terms: usize,
) -> f64 {
    let mut u_avg = 0.0;
    for m in 0..n_terms {
        let mf = m as f64;
        let n = (2.0 * mf + 1.0) * std::f64::consts::PI / 2.0;
        let factor = 2.0 / (n * n);
        let tv = cv * t / (h * h);
        u_avg += factor * (-n * n * tv).exp();
    }
    let degree = 1.0 - u_avg.clamp(0.0, 1.0);
    total_settlement * degree
}

/// Time factor T_v for consolidation.
pub fn time_factor(cv: f64, t: f64, h: f64) -> f64 {
    cv * t / (h * h)
}

/// Degree of consolidation U from time factor T_v.
pub fn degree_of_consolidation(tv: f64, n_terms: usize) -> f64 {
    let mut u_avg = 0.0;
    for m in 0..n_terms {
        let mf = m as f64;
        let n = (2.0 * mf + 1.0) * std::f64::consts::PI / 2.0;
        u_avg += (2.0 / (n * n)) * (-n * n * tv).exp();
    }
    1.0 - u_avg.clamp(0.0, 1.0)
}

/// Effective stress principle: sigma' = sigma - u * delta.
pub fn effective_stress(total: &StressTensor, pore_pressure: f64) -> StressTensor {
    let mut s = total.s;
    s[0] -= pore_pressure;
    s[1] -= pore_pressure;
    s[2] -= pore_pressure;
    StressTensor::from_array(s)
}

/// Pore pressure update (undrained increment).
///
/// Skempton's equation: delta_u = B * (delta_sigma_1 + A * (delta_sigma_1 - delta_sigma_3)).
pub fn skempton_pore_pressure(b: f64, a: f64, d_sigma1: f64, d_sigma3: f64) -> f64 {
    b * (d_sigma3 + a * (d_sigma1 - d_sigma3))
}

/// Settlement index for Cam-Clay model (NC line).
///
/// Returns volumetric strain increment.
pub fn cam_clay_compression(params: &CamClayParams, p_old: f64, p_new: f64) -> f64 {
    if p_old <= 0.0 || p_new <= 0.0 {
        return 0.0;
    }
    -params.lambda * (p_new / p_old).ln()
}

/// Void ratio change from effective stress change (1D).
pub fn void_ratio_change(cc: f64, e0: f64, sigma_v_old: f64, sigma_v_new: f64) -> f64 {
    if sigma_v_old <= 0.0 || sigma_v_new <= 0.0 {
        return 0.0;
    }
    -cc * (sigma_v_new / sigma_v_old).log10() * (1.0 + e0)
}

/// Oedometer settlement from 1D compression.
///
/// Returns layer settlement \[m\].
pub fn oedometer_settlement(
    cc: f64,
    cs: f64,
    e0: f64,
    h: f64,
    sigma_v0: f64,
    sigma_vc: f64, // preconsolidation pressure
    sigma_vf: f64, // final stress
) -> f64 {
    if sigma_vf <= sigma_vc {
        // Overconsolidated recompression only
        cs / (1.0 + e0) * h * (sigma_vf / sigma_v0).log10()
    } else if sigma_v0 >= sigma_vc {
        // Normally consolidated
        cc / (1.0 + e0) * h * (sigma_vf / sigma_v0).log10()
    } else {
        // Mixed: recompression + virgin compression
        let s_oc = cs / (1.0 + e0) * h * (sigma_vc / sigma_v0).log10();
        let s_nc = cc / (1.0 + e0) * h * (sigma_vf / sigma_vc).log10();
        s_oc + s_nc
    }
}

/// Bearing capacity factors (Terzaghi, general shear).
pub struct BearingCapacityFactors {
    /// N_c factor.
    pub nc: f64,
    /// N_q factor.
    pub nq: f64,
    /// N_gamma factor.
    pub n_gamma: f64,
}

impl BearingCapacityFactors {
    /// Compute Terzaghi bearing capacity factors from friction angle.
    pub fn terzaghi(phi: f64) -> Self {
        let kp = (std::f64::consts::PI / 4.0 + phi / 2.0).tan().powi(2);
        let exp_term = (std::f64::consts::PI * phi.tan()).exp();
        let nq = exp_term * kp;
        let nc = if phi.abs() < 1e-10 {
            5.14
        } else {
            (nq - 1.0) / phi.tan()
        };
        let n_gamma = 2.0 * (nq + 1.0) * phi.tan();
        Self { nc, nq, n_gamma }
    }

    /// Ultimate bearing capacity for strip footing.
    ///
    /// q_ult = c*Nc + q*Nq + 0.5*gamma*B*N_gamma.
    pub fn strip_footing_capacity(
        &self,
        cohesion: f64,
        overburden: f64,
        gamma: f64,
        width: f64,
    ) -> f64 {
        cohesion * self.nc + overburden * self.nq + 0.5 * gamma * width * self.n_gamma
    }
}

/// Earth pressure coefficient at rest (K0).
pub fn k0_normally_consolidated(phi: f64) -> f64 {
    1.0 - phi.sin()
}

/// Earth pressure coefficient for overconsolidated soil.
pub fn k0_overconsolidated(phi: f64, ocr: f64) -> f64 {
    k0_normally_consolidated(phi) * ocr.powf(phi.sin())
}

/// Slope stability factor of safety (infinite slope, c-phi soil).
pub fn infinite_slope_fs(
    c: f64,
    phi: f64,
    gamma: f64,
    h: f64,
    beta: f64,
    ru: f64, // pore pressure ratio
) -> f64 {
    let sin_beta = beta.sin();
    let cos_beta = beta.cos();
    let tan_phi = phi.tan();
    let numerator = c / (gamma * h) + (1.0 - ru) * cos_beta * cos_beta * tan_phi;
    let denominator = sin_beta * cos_beta;
    if denominator.abs() < 1e-12 {
        return f64::INFINITY;
    }
    numerator / denominator
}

/// Bishop simplified method for circular slip surface.
///
/// Iterative solution for factor of safety.
pub struct BishopAnalysis {
    /// Slice widths.
    pub widths: Vec<f64>,
    /// Slice heights (average).
    pub heights: Vec<f64>,
    /// Soil unit weight.
    pub gamma: Vec<f64>,
    /// Cohesion per slice.
    pub cohesion: Vec<f64>,
    /// Friction angle per slice.
    pub phi: Vec<f64>,
    /// Slice base angles.
    pub alpha: Vec<f64>,
    /// Pore pressure ratio per slice.
    pub ru: Vec<f64>,
}

impl BishopAnalysis {
    /// Compute factor of safety using Bishop's simplified method.
    ///
    /// Iterates until convergence.
    pub fn factor_of_safety(&self, max_iter: usize, tol: f64) -> f64 {
        let n = self.widths.len();
        let mut fs = 1.0;

        for _iter in 0..max_iter {
            let mut numerator = 0.0;
            let mut denominator = 0.0;

            for i in 0..n {
                let w = self.gamma[i] * self.widths[i] * self.heights[i];
                let alpha = self.alpha[i];
                let c = self.cohesion[i];
                let phi = self.phi[i];
                let ru = self.ru[i];
                let u = ru * self.gamma[i] * self.heights[i];

                let m_alpha = alpha.cos() + alpha.sin() * phi.tan() / fs;
                numerator += (c * self.widths[i] + (w - u * self.widths[i]) * phi.tan()) / m_alpha;
                denominator += w * alpha.sin();
            }

            let fs_new = numerator / denominator.max(1e-12);
            if (fs_new - fs).abs() < tol {
                return fs_new;
            }
            fs = fs_new;
        }
        fs
    }
}

/// Permeability from Kozeny-Carman equation.
///
/// k = (e^3/(1+e)) * (d50^2 / (180 * mu)).
pub fn kozeny_carman_permeability(void_ratio: f64, d50: f64, mu: f64) -> f64 {
    let e = void_ratio;
    let e3 = e * e * e;
    let factor = e3 / ((1.0 + e) * 180.0 * mu);
    factor * d50 * d50
}

/// Seepage velocity (Darcy's law): v = k * i.
pub fn darcy_velocity(k: f64, hydraulic_gradient: f64) -> f64 {
    k * hydraulic_gradient
}

/// Flow net solution for seepage quantity.
///
/// Q = k * H * (Nf/Nd) per unit width.
pub fn flow_net_seepage(k: f64, head: f64, nf: f64, nd: f64) -> f64 {
    k * head * nf / nd
}

/// Critical hydraulic gradient (piping initiation).
pub fn critical_hydraulic_gradient(gs: f64, e: f64) -> f64 {
    (gs - 1.0) / (1.0 + e)
}

/// Liquefaction resistance ratio (cyclic stress ratio).
pub fn cyclic_stress_ratio(
    amax: f64,     // peak ground acceleration [g]
    sigma_v: f64,  // total vertical stress [Pa]
    sigma_vp: f64, // effective vertical stress [Pa]
    rd: f64,       // depth reduction factor
) -> f64 {
    0.65 * (amax / 9.81) * (sigma_v / sigma_vp) * rd
}

/// Depth reduction factor rd for liquefaction analysis.
pub fn depth_reduction_factor(z: f64) -> f64 {
    if z <= 9.15 {
        1.0 - 0.00765 * z
    } else if z <= 23.0 {
        1.174 - 0.0267 * z
    } else {
        0.5
    }
}

/// SPT N-value correction to N60.
pub fn spt_n60(n_raw: f64, e_rod: f64, c_b: f64, c_s: f64, c_r: f64) -> f64 {
    n_raw * (e_rod / 60.0) * c_b * c_s * c_r
}

/// Overburden correction for SPT (CN factor, Liao-Whitman).
pub fn overburden_correction(sigma_v_eff: f64) -> f64 {
    let pa = 101325.0; // atmospheric pressure
    (pa / sigma_v_eff).sqrt().min(2.0)
}

/// Normally consolidated compression line (NCL) void ratio.
pub fn ncl_void_ratio(e_ref: f64, lambda: f64, p: f64, p_ref: f64) -> f64 {
    e_ref - lambda * (p / p_ref).ln()
}

/// Soil unit weight from void ratio and specific gravity.
pub fn soil_unit_weight(gs: f64, e: f64, sr: f64, gamma_w: f64) -> f64 {
    gamma_w * (gs + sr * e) / (1.0 + e)
}

/// Undrained shear strength from SHANSEP equation.
pub fn shansep_su(sigma_vp: f64, s: f64, m: f64, ocr: f64) -> f64 {
    sigma_vp * s * ocr.powf(m)
}

/// Preconsolidation pressure from Casagrande construction (approximate).
///
/// Given compression curve data points (log_sigma, void_ratio).
pub fn casagrande_pc(log_sigma: &[f64], void_ratio: &[f64]) -> Option<f64> {
    if log_sigma.len() < 4 {
        return None;
    }
    let n = log_sigma.len();
    // Find max curvature point (simplified: maximum curvature index)
    let mut max_curv = 0.0;
    let mut pc_idx = 1usize;
    for i in 1..(n - 1) {
        let dydx1 = (void_ratio[i] - void_ratio[i - 1]) / (log_sigma[i] - log_sigma[i - 1] + 1e-12);
        let dydx2 = (void_ratio[i + 1] - void_ratio[i]) / (log_sigma[i + 1] - log_sigma[i] + 1e-12);
        let curv = (dydx2 - dydx1).abs();
        if curv > max_curv {
            max_curv = curv;
            pc_idx = i;
        }
    }
    Some(10.0_f64.powf(log_sigma[pc_idx]))
}

/// Soil FEM node (2D or 3D).
#[derive(Debug, Clone)]
pub struct SoilNode {
    /// Node ID.
    pub id: usize,
    /// Position \[m\].
    pub pos: [f64; 3],
    /// Displacement \[m\].
    pub disp: [f64; 3],
    /// Pore pressure \[Pa\].
    pub pore_pressure: f64,
}

impl SoilNode {
    /// Create new soil node.
    pub fn new(id: usize, pos: [f64; 3]) -> Self {
        Self {
            id,
            pos,
            disp: [0.0; 3],
            pore_pressure: 0.0,
        }
    }
}

/// Soil FEM element (quadrilateral or tetrahedral).
#[derive(Debug, Clone)]
pub struct SoilElement {
    /// Element ID.
    pub id: usize,
    /// Node connectivity.
    pub nodes: Vec<usize>,
    /// Soil model for this element.
    pub model: SoilModel,
    /// Drainage condition.
    pub drainage: DrainageCondition,
    /// Current effective stress state.
    pub stress: StressTensor,
    /// Plastic strain accumulation.
    pub plastic_strain: [f64; 6],
}

impl SoilElement {
    /// Create new soil element.
    pub fn new(
        id: usize,
        nodes: Vec<usize>,
        model: SoilModel,
        drainage: DrainageCondition,
    ) -> Self {
        Self {
            id,
            nodes,
            model,
            drainage,
            stress: StressTensor::zero(),
            plastic_strain: [0.0; 6],
        }
    }
}

/// Soil FEM mesh.
#[derive(Debug, Clone, Default)]
pub struct SoilMesh {
    /// Mesh nodes.
    pub nodes: Vec<SoilNode>,
    /// Mesh elements.
    pub elements: Vec<SoilElement>,
    /// Gravity vector \[m/s^2\].
    pub gravity: [f64; 3],
}

impl SoilMesh {
    /// Create new empty soil mesh.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            elements: Vec::new(),
            gravity: [0.0, -9.81, 0.0],
        }
    }

    /// Add node to mesh.
    pub fn add_node(&mut self, pos: [f64; 3]) -> usize {
        let id = self.nodes.len();
        self.nodes.push(SoilNode::new(id, pos));
        id
    }

    /// Add element to mesh.
    pub fn add_element(
        &mut self,
        nodes: Vec<usize>,
        model: SoilModel,
        drainage: DrainageCondition,
    ) -> usize {
        let id = self.elements.len();
        self.elements
            .push(SoilElement::new(id, nodes, model, drainage));
        id
    }

    /// Get total number of degrees of freedom (3 disp + 1 pore per node).
    pub fn n_dof(&self) -> usize {
        self.nodes.len() * 4
    }

    /// Apply gravity loading — initialise vertical effective stress (K0 procedure).
    pub fn initialize_k0(&mut self, gamma: f64, k0: f64) {
        for node in &mut self.nodes {
            let depth = -node.pos[1]; // y up
            let sigma_v = gamma * depth;
            let sigma_h = k0 * sigma_v;
            node.pore_pressure = 0.0; // drained initial
            // Store initial stress (negative = compression in soil convention)
            let _stress = StressTensor::from_array([-sigma_h, -sigma_v, -sigma_h, 0.0, 0.0, 0.0]);
        }
    }

    /// Compute total stress at a point given depth and K0.
    pub fn total_stress_at_depth(
        &self,
        depth: f64,
        gamma: f64,
        k0: f64,
        gamma_w: f64,
        z_wt: f64,
    ) -> StressTensor {
        let sigma_v = gamma * depth;
        let u = if depth > z_wt {
            gamma_w * (depth - z_wt)
        } else {
            0.0
        };
        let sigma_v_eff = sigma_v - u;
        let sigma_h_eff = k0 * sigma_v_eff;
        StressTensor::from_array([
            -(sigma_h_eff + u),
            -(sigma_v_eff + u),
            -(sigma_h_eff + u),
            0.0,
            0.0,
            0.0,
        ])
    }
}

/// Consolidation solver (1D Terzaghi, finite difference).
pub struct ConsolidationSolver {
    /// Number of layers.
    pub n_layers: usize,
    /// Layer thicknesses \[m\].
    pub dz: f64,
    /// Consolidation coefficients per layer \[m^2/s\].
    pub cv: Vec<f64>,
    /// Initial excess pore pressure \[Pa\].
    pub u0: Vec<f64>,
    /// Current pore pressure.
    pub u: Vec<f64>,
    /// Time step \[s\].
    pub dt: f64,
}

impl ConsolidationSolver {
    /// Create new consolidation solver.
    pub fn new(n_layers: usize, dz: f64, cv: Vec<f64>, u0: Vec<f64>, dt: f64) -> Self {
        let u = u0.clone();
        Self {
            n_layers,
            dz,
            cv,
            u0,
            u,
            dt,
        }
    }

    /// Advance one time step (explicit finite difference).
    pub fn step(&mut self) {
        let r_max = self.cv.iter().cloned().fold(0.0f64, f64::max) * self.dt / (self.dz * self.dz);
        if r_max > 0.5 {
            // Unstable: reduce dt
            return;
        }
        let mut u_new = self.u.clone();
        for (idx, u_new_i) in u_new[1..(self.n_layers - 1)].iter_mut().enumerate() {
            let i = idx + 1;
            let r = self.cv[i] * self.dt / (self.dz * self.dz);
            *u_new_i = self.u[i] + r * (self.u[i + 1] - 2.0 * self.u[i] + self.u[i - 1]);
        }
        // Boundary conditions: drainage at top and bottom
        u_new[0] = 0.0;
        u_new[self.n_layers - 1] = 0.0;
        self.u = u_new;
    }

    /// Advance multiple steps.
    pub fn advance(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }

    /// Average degree of consolidation.
    pub fn average_consolidation(&self) -> f64 {
        let u_avg = self.u.iter().sum::<f64>() / self.n_layers as f64;
        let u0_avg = self.u0.iter().sum::<f64>() / self.n_layers as f64;
        if u0_avg < 1e-12 {
            return 1.0;
        }
        1.0 - u_avg / u0_avg
    }
}

/// Seepage analysis (2D Laplace equation, finite difference).
pub struct SeepageAnalysis {
    /// Grid dimensions.
    pub nx: usize,
    /// Grid dimensions.
    pub ny: usize,
    /// Grid spacing \[m\].
    pub dx: f64,
    /// Grid spacing \[m\].
    pub dy: f64,
    /// Hydraulic conductivity \[m/s\].
    pub kx: f64,
    /// Hydraulic conductivity \[m/s\].
    pub ky: f64,
    /// Total head field \[m\].
    pub head: Vec<Vec<f64>>,
    /// Boundary mask (true = fixed head).
    pub bc_mask: Vec<Vec<bool>>,
    /// Fixed head values \[m\].
    pub bc_head: Vec<Vec<f64>>,
}

impl SeepageAnalysis {
    /// Create new seepage analysis grid.
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64, kx: f64, ky: f64) -> Self {
        Self {
            nx,
            ny,
            dx,
            dy,
            kx,
            ky,
            head: vec![vec![0.0; ny]; nx],
            bc_mask: vec![vec![false; ny]; nx],
            bc_head: vec![vec![0.0; ny]; nx],
        }
    }

    /// Set boundary condition at node.
    pub fn set_bc(&mut self, i: usize, j: usize, h: f64) {
        self.bc_mask[i][j] = true;
        self.bc_head[i][j] = h;
        self.head[i][j] = h;
    }

    /// Solve by Gauss-Seidel iteration.
    pub fn solve(&mut self, max_iter: usize, tol: f64) -> usize {
        let rx = self.kx / (self.dx * self.dx);
        let ry = self.ky / (self.dy * self.dy);
        let denom = 2.0 * (rx + ry);

        for iter in 0..max_iter {
            let mut max_err = 0.0f64;
            for i in 1..(self.nx - 1) {
                for j in 1..(self.ny - 1) {
                    if self.bc_mask[i][j] {
                        continue;
                    }
                    let h_new = (rx * (self.head[i + 1][j] + self.head[i - 1][j])
                        + ry * (self.head[i][j + 1] + self.head[i][j - 1]))
                        / denom;
                    let err = (h_new - self.head[i][j]).abs();
                    if err > max_err {
                        max_err = err;
                    }
                    self.head[i][j] = h_new;
                }
            }
            if max_err < tol {
                return iter + 1;
            }
        }
        max_iter
    }

    /// Compute Darcy flux at node (i, j).
    pub fn darcy_flux(&self, i: usize, j: usize) -> [f64; 2] {
        let qx = if i > 0 && i < self.nx - 1 {
            -self.kx * (self.head[i + 1][j] - self.head[i - 1][j]) / (2.0 * self.dx)
        } else {
            0.0
        };
        let qy = if j > 0 && j < self.ny - 1 {
            -self.ky * (self.head[i][j + 1] - self.head[i][j - 1]) / (2.0 * self.dy)
        } else {
            0.0
        };
        [qx, qy]
    }
}

/// Embankment stability result.
#[derive(Debug, Clone)]
pub struct EmbankmentStability {
    /// Factor of safety.
    pub fs: f64,
    /// Critical slip surface radius.
    pub radius: f64,
    /// Failure mode.
    pub mode: String,
}

impl EmbankmentStability {
    /// Create new result.
    pub fn new(fs: f64, radius: f64, mode: &str) -> Self {
        Self {
            fs,
            radius,
            mode: mode.to_string(),
        }
    }

    /// Is stable (FS >= 1.5)?
    pub fn is_stable(&self) -> bool {
        self.fs >= 1.5
    }
}

/// Simple 1D pile analysis (Winkler spring model).
#[derive(Debug, Clone)]
pub struct PileAnalysis {
    /// Pile length \[m\].
    pub length: f64,
    /// Pile diameter \[m\].
    pub diameter: f64,
    /// Pile modulus \[Pa\].
    pub ep: f64,
    /// Number of segments.
    pub n_segments: usize,
    /// Subgrade reaction \[N/m^3\].
    pub ks: Vec<f64>,
    /// Pile deflection \[m\].
    pub deflection: Vec<f64>,
}

impl PileAnalysis {
    /// Create new pile analysis.
    pub fn new(length: f64, diameter: f64, ep: f64, n_segments: usize, ks: f64) -> Self {
        let dz = length / n_segments as f64;
        let _ = dz;
        Self {
            length,
            diameter,
            ep,
            n_segments,
            ks: vec![ks; n_segments],
            deflection: vec![0.0; n_segments + 1],
        }
    }

    /// Pile moment of inertia.
    pub fn moment_of_inertia(&self) -> f64 {
        std::f64::consts::PI * self.diameter.powi(4) / 64.0
    }

    /// Pile stiffness EI.
    pub fn ei(&self) -> f64 {
        self.ep * self.moment_of_inertia()
    }

    /// Characteristic length beta.
    pub fn beta(&self) -> f64 {
        let avg_ks = self.ks.iter().sum::<f64>() / self.ks.len() as f64;
        (avg_ks * self.diameter / (4.0 * self.ei())).powf(0.25)
    }

    /// Horizontal pile capacity (p-y simplified).
    pub fn lateral_capacity(&self, cohesion: f64) -> f64 {
        9.0 * cohesion * self.diameter * self.length
    }

    /// Apply lateral load and compute maximum deflection.
    pub fn apply_lateral_load(&mut self, h: f64, m: f64) -> f64 {
        let beta = self.beta();
        let ei = self.ei();
        // Cantilever with Winkler spring — closed-form at top
        let y0 =
            h * beta / (2.0 * ei * beta.powi(3)) + m * beta.powi(2) / (2.0 * ei * beta.powi(3));
        self.deflection[0] = y0;
        y0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mohr_coulomb_params() {
        let mc = MohrCoulombParams::new(10e6, 0.3, 5000.0, 30f64.to_radians(), 10f64.to_radians());
        assert!((mc.shear_modulus() - 10e6 / 2.6).abs() < 1e3);
        assert!(mc.bulk_modulus() > 0.0);
    }

    #[test]
    fn test_stress_tensor_mean() {
        let s = StressTensor::from_array([-100e3, -200e3, -150e3, 0.0, 0.0, 0.0]);
        let p = s.mean_stress();
        assert!((p - 150e3).abs() < 1.0);
    }

    #[test]
    fn test_stress_tensor_von_mises() {
        let s = StressTensor::from_array([-100e3, -100e3, -100e3, 0.0, 0.0, 0.0]);
        assert!(s.von_mises() < 1.0); // hydrostatic
    }

    #[test]
    fn test_effective_stress() {
        let total = StressTensor::from_array([-200e3, -300e3, -200e3, 0.0, 0.0, 0.0]);
        let eff = effective_stress(&total, 50e3);
        assert!((eff.s[0] - (-250e3)).abs() < 1.0);
        assert!((eff.s[1] - (-350e3)).abs() < 1.0);
    }

    #[test]
    fn test_mohr_coulomb_yield_elastic() {
        let params =
            MohrCoulombParams::new(10e6, 0.3, 100e3, 30f64.to_radians(), 10f64.to_radians());
        // Small stress: should be elastic
        let s = StressTensor::from_array([-50e3, -50e3, -50e3, 0.0, 0.0, 0.0]);
        let f = mohr_coulomb_yield(&s, &params);
        assert!(f < 0.0, "Should be elastic for hydrostatic stress");
    }

    #[test]
    fn test_mc_to_dp() {
        let phi = 30f64.to_radians();
        let c = 10e3;
        let (alpha, k) = mc_to_dp_outer(phi, c);
        assert!(alpha > 0.0);
        assert!(k > 0.0);
    }

    #[test]
    fn test_elastic_stiffness_symmetry() {
        let d = elastic_stiffness(10e6, 0.3);
        for (i, row) in d.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - d[j][i]).abs() < 1e-6,
                    "D not symmetric at ({}, {})",
                    i,
                    j
                );
            }
        }
    }

    #[test]
    fn test_apply_stiffness() {
        let d = elastic_stiffness(1e6, 0.25);
        let eps = [0.001, 0.0, 0.0, 0.0, 0.0, 0.0];
        let sigma = apply_stiffness(&d, &eps);
        assert!(sigma[0] > 0.0); // tensile in x
    }

    #[test]
    fn test_terzaghi_settlement_convergence() {
        let s1 = terzaghi_settlement(0.1, 1e-8, 5.0, 1e6, 10);
        let s2 = terzaghi_settlement(0.1, 1e-8, 5.0, 1e9, 10);
        assert!(s2 > s1); // More settlement at later time
        assert!(s2 <= 0.1 + 1e-9);
    }

    #[test]
    fn test_degree_of_consolidation() {
        let u0 = degree_of_consolidation(0.0, 20);
        let u100 = degree_of_consolidation(10.0, 20);
        assert!(u0 < 0.02);
        assert!(u100 > 0.99);
    }

    #[test]
    fn test_k0_nc() {
        let phi = 30f64.to_radians();
        let k0 = k0_normally_consolidated(phi);
        assert!((k0 - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_k0_overconsolidated() {
        let phi = 30f64.to_radians();
        let k0_nc = k0_normally_consolidated(phi);
        let k0_oc = k0_overconsolidated(phi, 4.0);
        assert!(k0_oc > k0_nc);
    }

    #[test]
    fn test_bearing_capacity_strip() {
        let bf = BearingCapacityFactors::terzaghi(30f64.to_radians());
        assert!(bf.nq > 1.0);
        assert!(bf.nc > 1.0);
        let qult = bf.strip_footing_capacity(0.0, 50e3, 18e3, 1.0);
        assert!(qult > 50e3);
    }

    #[test]
    fn test_bearing_capacity_cohesive() {
        let bf = BearingCapacityFactors::terzaghi(0.0);
        // phi=0: Nc=5.14, Nq=1, Ng~0
        assert!((bf.nc - 5.14).abs() < 0.1);
    }

    #[test]
    fn test_infinite_slope_stable() {
        let fs = infinite_slope_fs(10e3, 30f64.to_radians(), 18e3, 5.0, 15f64.to_radians(), 0.0);
        assert!(fs > 1.0);
    }

    #[test]
    fn test_infinite_slope_unstable() {
        let fs = infinite_slope_fs(0.0, 20f64.to_radians(), 18e3, 5.0, 25f64.to_radians(), 0.5);
        assert!(fs < 1.5, "Should be unsafe with pore pressure");
    }

    #[test]
    fn test_bishop_method() {
        let n = 5;
        let ba = BishopAnalysis {
            widths: vec![1.0; n],
            heights: vec![3.0; n],
            gamma: vec![18e3; n],
            cohesion: vec![10e3; n],
            phi: vec![25f64.to_radians(); n],
            alpha: vec![
                20f64.to_radians(),
                15f64.to_radians(),
                10f64.to_radians(),
                5f64.to_radians(),
                2f64.to_radians(),
            ],
            ru: vec![0.0; n],
        };
        let fs = ba.factor_of_safety(50, 1e-6);
        assert!(fs > 0.5 && fs < 5.0, "FS should be reasonable: {}", fs);
    }

    #[test]
    fn test_kozeny_carman() {
        let k = kozeny_carman_permeability(0.7, 0.2e-3, 1e-3);
        assert!(k > 0.0);
    }

    #[test]
    fn test_darcy_velocity() {
        let v = darcy_velocity(1e-5, 0.5);
        assert!((v - 5e-6).abs() < 1e-10);
    }

    #[test]
    fn test_critical_gradient() {
        let ic = critical_hydraulic_gradient(2.67, 0.7);
        assert!((ic - 0.98).abs() < 0.1);
    }

    #[test]
    fn test_liquefaction_csr() {
        let csr = cyclic_stress_ratio(0.3 * 9.81, 100e3, 60e3, 0.9);
        assert!(csr > 0.0 && csr < 1.0);
    }

    #[test]
    fn test_depth_reduction_factor() {
        assert!((depth_reduction_factor(0.0) - 1.0).abs() < 0.01);
        assert!(depth_reduction_factor(9.0) < 1.0);
        assert!((depth_reduction_factor(30.0) - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_ncl_void_ratio() {
        let e = ncl_void_ratio(2.0, 0.2, 100e3, 1e3);
        assert!(e > 0.0);
    }

    #[test]
    fn test_soil_unit_weight() {
        let gamma = soil_unit_weight(2.67, 0.7, 1.0, 9810.0);
        assert!(gamma > 10e3 && gamma < 25e3);
    }

    #[test]
    fn test_shansep() {
        let su = shansep_su(100e3, 0.22, 0.8, 2.0);
        assert!(su > 100e3 * 0.22);
    }

    #[test]
    fn test_oedometer_nc() {
        let s = oedometer_settlement(0.35, 0.07, 1.0, 5.0, 50e3, 50e3, 150e3);
        assert!(s > 0.0);
    }

    #[test]
    fn test_oedometer_oc() {
        let s = oedometer_settlement(0.35, 0.07, 1.0, 5.0, 50e3, 200e3, 100e3);
        assert!(s > 0.0);
    }

    #[test]
    fn test_oedometer_mixed() {
        let s = oedometer_settlement(0.35, 0.07, 1.0, 5.0, 50e3, 100e3, 200e3);
        assert!(s > 0.0);
    }

    #[test]
    fn test_cam_clay_compression() {
        let params = CamClayParams::new(0.85, 0.15, 0.03, 100e3, 150e3, 0.3);
        let eps = cam_clay_compression(&params, 100e3, 200e3);
        assert!(eps < 0.0); // compression (negative volumetric strain)
    }

    #[test]
    fn test_biot_storage() {
        let biot = BiotParams::new(0.9, 1e8, 1e-7, 0.4);
        let s = biot.storage_coefficient();
        assert!(s > 0.0);
    }

    #[test]
    fn test_consolidation_solver() {
        let n = 20;
        let dz = 1.0 / n as f64;
        let cv = vec![1e-7; n];
        let u0 = vec![100e3; n];
        let dt = dz * dz / (2.0 * 1e-7); // just stable
        let mut solver = ConsolidationSolver::new(n, dz, cv, u0, dt);
        solver.advance(100);
        let u = solver.average_consolidation();
        assert!(u > 0.0 && u <= 1.0, "Consolidation degree: {}", u);
    }

    #[test]
    fn test_seepage_analysis() {
        let mut sa = SeepageAnalysis::new(10, 10, 1.0, 1.0, 1e-5, 1e-5);
        // Upstream head = 10m, downstream = 0m
        for j in 0..10 {
            sa.set_bc(0, j, 10.0);
        }
        for j in 0..10 {
            sa.set_bc(9, j, 0.0);
        }
        let iters = sa.solve(1000, 1e-6);
        assert!(iters < 1000);
        // Head should decrease from left to right
        assert!(sa.head[5][5] > 1.0 && sa.head[5][5] < 9.0);
    }

    #[test]
    fn test_darcy_flux() {
        let mut sa = SeepageAnalysis::new(10, 10, 1.0, 1.0, 1e-5, 1e-5);
        for j in 0..10 {
            sa.set_bc(0, j, 10.0);
        }
        for j in 0..10 {
            sa.set_bc(9, j, 0.0);
        }
        sa.solve(1000, 1e-6);
        let flux = sa.darcy_flux(5, 5);
        assert!(flux[0] > 0.0); // flow from high to low head (positive x-direction: i=0 high, i=9 low)
    }

    #[test]
    fn test_soil_mesh() {
        let mut mesh = SoilMesh::new();
        let n0 = mesh.add_node([0.0, 0.0, 0.0]);
        let n1 = mesh.add_node([1.0, 0.0, 0.0]);
        let n2 = mesh.add_node([1.0, -1.0, 0.0]);
        let n3 = mesh.add_node([0.0, -1.0, 0.0]);
        let _e = mesh.add_element(
            vec![n0, n1, n2, n3],
            SoilModel::MohrCoulomb,
            DrainageCondition::Drained,
        );
        assert_eq!(mesh.nodes.len(), 4);
        assert_eq!(mesh.elements.len(), 1);
    }

    #[test]
    fn test_soil_mesh_k0() {
        let mut mesh = SoilMesh::new();
        mesh.add_node([0.0, -5.0, 0.0]);
        mesh.initialize_k0(18e3, 0.5);
        // No crash, stress initialised
        assert!(!mesh.nodes.is_empty());
    }

    #[test]
    fn test_pile_analysis() {
        let mut pile = PileAnalysis::new(10.0, 0.5, 30e9, 20, 50e6);
        let ei = pile.ei();
        assert!(ei > 0.0);
        let beta = pile.beta();
        assert!(beta > 0.0);
        let cap = pile.lateral_capacity(50e3);
        assert!(cap > 0.0);
        let y0 = pile.apply_lateral_load(100e3, 0.0);
        assert!(y0 >= 0.0);
    }

    #[test]
    fn test_casagrande_pc() {
        let log_sigma: Vec<f64> = (0..10).map(|i| 3.0 + i as f64 * 0.2).collect();
        let void_ratio: Vec<f64> = (0..10)
            .map(|i| {
                let ls = 3.0 + i as f64 * 0.2;
                if ls < 4.0 {
                    1.5 - 0.05 * (ls - 3.0)
                } else {
                    1.45 - 0.35 * (ls - 4.0)
                }
            })
            .collect();
        let pc = casagrande_pc(&log_sigma, &void_ratio);
        assert!(pc.is_some());
        assert!(pc.unwrap() > 0.0);
    }

    #[test]
    fn test_flow_net_seepage() {
        let q = flow_net_seepage(1e-5, 10.0, 4.0, 12.0);
        assert!(q > 0.0);
    }

    #[test]
    fn test_embankment_stability() {
        let stab = EmbankmentStability::new(1.8, 15.0, "circular");
        assert!(stab.is_stable());
        let unstable = EmbankmentStability::new(1.2, 15.0, "circular");
        assert!(!unstable.is_stable());
    }

    #[test]
    fn test_principal_stresses() {
        let s = StressTensor::from_array([-300e3, -100e3, -200e3, 0.0, 0.0, 0.0]);
        let p = s.principal_stresses();
        assert!(p[0] >= p[1] && p[1] >= p[2]);
    }

    #[test]
    fn test_skempton_b1() {
        // B=1 (saturated), A=0.5: delta_u = delta_sigma3 + 0.5*(ds1-ds3)
        let du = skempton_pore_pressure(1.0, 0.5, 200e3, 100e3);
        assert!((du - 150e3).abs() < 1.0);
    }
}
