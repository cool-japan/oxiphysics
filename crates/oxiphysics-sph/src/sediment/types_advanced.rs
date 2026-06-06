//\! Advanced sediment transport types
//\!
//\! Contains KynchSettling, TurbidityCurrentSph, ShieldsParam, and other
//\! advanced sediment transport models.

use super::functions::*;
use super::types::SedimentTransportMode;

/// Kynch batch settling theory for concentrated suspensions.
///
/// Implements the Kynch flux theory: ∂φ/∂t + ∂F(φ)/∂z = 0
/// where F(φ) = φ v_s(φ) is the sediment flux.
pub struct KynchSettling {
    /// Settling velocity at dilute limit v_s0 (m/s).
    pub v_s0: f64,
    /// Maximum packing concentration φ_max (−).
    pub phi_max: f64,
    /// Richardson-Zaki exponent n (typically 4.65 for laminar regime).
    pub rz_exponent: f64,
}
impl KynchSettling {
    /// Construct a Kynch settling model.
    pub fn new(v_s0: f64, phi_max: f64, rz_exponent: f64) -> Self {
        KynchSettling {
            v_s0,
            phi_max,
            rz_exponent,
        }
    }
    /// Richardson-Zaki hindered settling velocity: v_s(φ) = v_s0 (1-φ/φ_max)^n.
    pub fn hindered_velocity(&self, phi: f64) -> f64 {
        if phi >= self.phi_max {
            return 0.0;
        }
        self.v_s0 * (1.0 - phi / self.phi_max).powf(self.rz_exponent)
    }
    /// Sediment flux F(φ) = φ v_s(φ).
    pub fn flux(&self, phi: f64) -> f64 {
        phi * self.hindered_velocity(phi)
    }
    /// Characteristic speed dF/dφ (wave propagation speed in Kynch theory).
    pub fn characteristic_speed(&self, phi: f64) -> f64 {
        let dphi = 1e-6;
        if phi + dphi >= self.phi_max {
            return 0.0;
        }
        (self.flux(phi + dphi) - self.flux(phi - dphi.min(phi))) / (2.0 * dphi)
    }
    /// Critical concentration at which flux is maximum.
    pub fn critical_concentration(&self) -> f64 {
        let n_pts = 100usize;
        let mut best_phi = 0.0;
        let mut best_f = 0.0;
        for k in 0..n_pts {
            let phi = self.phi_max * (k as f64 + 0.5) / n_pts as f64;
            let f = self.flux(phi);
            if f > best_f {
                best_f = f;
                best_phi = phi;
            }
        }
        best_phi
    }
}
/// Turbidity current model: density-driven underflow carrying suspended sediment.
///
/// Implements layer-averaged equations following Parker (1986) framework.
pub struct TurbidityCurrentSph {
    /// Layer thickness h (m).
    pub layer_thickness: f64,
    /// Depth-averaged suspended concentration C (−).
    pub concentration: f64,
    /// Depth-averaged current velocity u (m/s).
    pub velocity: f64,
    /// Bed slope S.
    pub bed_slope: f64,
    /// Drag coefficient C_D.
    pub drag_coeff: f64,
    /// Entrainment coefficient e_w.
    pub entrainment_coeff: f64,
    /// Sediment settling velocity w_s (m/s).
    pub settling_velocity: f64,
    /// Fluid density (ambient) ρ_a (kg/m³).
    pub rho_ambient: f64,
    /// Sediment density ρ_s (kg/m³).
    pub rho_sediment: f64,
}
impl TurbidityCurrentSph {
    /// Construct a turbidity current model.
    pub fn new(
        layer_thickness: f64,
        concentration: f64,
        velocity: f64,
        bed_slope: f64,
        drag_coeff: f64,
        entrainment_coeff: f64,
        settling_velocity: f64,
        rho_ambient: f64,
        rho_sediment: f64,
    ) -> Self {
        TurbidityCurrentSph {
            layer_thickness,
            concentration,
            velocity,
            bed_slope,
            drag_coeff,
            entrainment_coeff,
            settling_velocity,
            rho_ambient,
            rho_sediment,
        }
    }
    /// Reduced gravity g' = g (ρ_m - ρ_a) / ρ_a.
    pub fn reduced_gravity(&self) -> f64 {
        let rho_m = self.mixture_density();
        G * (rho_m - self.rho_ambient) / self.rho_ambient.max(1e-300)
    }
    /// Mixture density ρ_m = ρ_a + C (ρ_s - ρ_a).
    pub fn mixture_density(&self) -> f64 {
        self.rho_ambient + self.concentration * (self.rho_sediment - self.rho_ambient)
    }
    /// Densimetric Froude number Fr_d = u / sqrt(g' h).
    pub fn froude_number(&self) -> f64 {
        let gp = self.reduced_gravity();
        let h = self.layer_thickness;
        if gp < 1e-300 || h < 1e-300 {
            return 0.0;
        }
        self.velocity / (gp * h).sqrt()
    }
    /// Bed shear stress τ_b = ρ_m C_D u² (Pa).
    pub fn bed_shear_stress(&self) -> f64 {
        self.mixture_density() * self.drag_coeff * self.velocity * self.velocity
    }
    /// Momentum equation: du/dt = g' S - C_D u²/h - g' ∂h/∂x contributions.
    pub fn du_dt(&self) -> f64 {
        let gp = self.reduced_gravity();
        let driving = gp * self.bed_slope;
        let friction =
            self.drag_coeff * self.velocity * self.velocity / self.layer_thickness.max(1e-300);
        driving - friction
    }
    /// Concentration equation: dC/dt = (E - D) / h.
    ///
    /// E = erosion rate (volumetric), D = deposition rate = w_s * C.
    pub fn dc_dt(&self, erosion_rate: f64) -> f64 {
        let deposition = self.settling_velocity * self.concentration;
        (erosion_rate - deposition) / self.layer_thickness.max(1e-300)
    }
    /// Advance one time step using explicit Euler.
    pub fn advance(&mut self, erosion_rate: f64, dt: f64) {
        let du = self.du_dt() * dt;
        let dc = self.dc_dt(erosion_rate) * dt;
        self.velocity += du;
        self.concentration = (self.concentration + dc).max(0.0);
    }
}
/// Shields parameter for bed-load transport initiation.
pub struct ShieldsParam {
    /// Grain diameter d (m).
    pub grain_diameter: f64,
    /// Sediment density ρ_s (kg/m³).
    pub rho_s: f64,
    /// Fluid density ρ_f (kg/m³).
    pub rho_f: f64,
    /// Critical Shields parameter θ_cr (typically ~0.047 for sand).
    pub theta_cr: f64,
}
impl ShieldsParam {
    /// Create a Shields parameter object.
    pub fn new(grain_diameter: f64, rho_s: f64, rho_f: f64, theta_cr: f64) -> Self {
        ShieldsParam {
            grain_diameter,
            rho_s,
            rho_f,
            theta_cr,
        }
    }
    /// Compute Shields stress θ from bed shear stress τ_b.
    pub fn theta(&self, tau_b: f64) -> f64 {
        shields_stress(tau_b, self.grain_diameter, self.rho_s, self.rho_f)
    }
    /// Is the bed-load motion initiated for bed shear stress `tau_b`?
    pub fn is_mobile(&self, tau_b: f64) -> bool {
        self.theta(tau_b) > self.theta_cr
    }
    /// Excess Shields stress (θ − θ_cr), clamped to 0.
    pub fn excess_shields(&self, tau_b: f64) -> f64 {
        (self.theta(tau_b) - self.theta_cr).max(0.0)
    }
    /// Dimensionless grain size D_* = d * \[(s-1)*g/ν²\]^{1/3}.
    pub fn dimensionless_grain_size(&self, nu: f64) -> f64 {
        let s = self.rho_s / self.rho_f;
        let d = self.grain_diameter;
        d * ((s - 1.0) * G / (nu * nu)).powf(1.0 / 3.0)
    }
}
/// Time-evolving scour around a pier or abutment.
///
/// Tracks the scour hole geometry and applies feedback to the flow field.
pub struct TimeEvolvingScour {
    /// Current scour depth (m).
    pub depth: f64,
    /// Equilibrium scour depth y_se (m).
    pub equilibrium_depth: f64,
    /// Time scale of scour T_s (s).
    pub time_scale: f64,
    /// Current time (s).
    pub time: f64,
}
impl TimeEvolvingScour {
    /// Construct a time-evolving scour model.
    pub fn new(equilibrium_depth: f64, time_scale: f64) -> Self {
        TimeEvolvingScour {
            depth: 0.0,
            equilibrium_depth,
            time_scale,
            time: 0.0,
        }
    }
    /// Temporal scour depth y_s(t) = y_se (1 - exp(-t/T_s)).
    pub fn scour_at_time(&self, t: f64) -> f64 {
        self.equilibrium_depth * (1.0 - (-t / self.time_scale.max(1e-300)).exp())
    }
    /// Advance scour evolution.
    pub fn advance(&mut self, dt: f64) {
        self.time += dt;
        self.depth = self.scour_at_time(self.time);
    }
    /// Rate of scour dDs/dt = (y_se - y_s) / T_s.
    pub fn scour_rate(&self) -> f64 {
        (self.equilibrium_depth - self.depth) / self.time_scale.max(1e-300)
    }
}
/// One-line shoreline evolution model (Pelnard-Considère, 1956).
///
/// Solves ∂y_s/∂t = K_2 ∂²y_s/∂x² + sources where y_s is
/// the shoreline position and K_2 is the longshore diffusivity.
pub struct ShorelineEvolution {
    /// Shoreline position y_s (m) at each along-shore cell.
    pub position: Vec<f64>,
    /// Along-shore grid spacing Δx (m).
    pub dx: f64,
    /// Longshore diffusivity K_2 = K H_b^{5/2} / (8 (D_c+D_b)) (m²/s).
    pub diffusivity: f64,
    /// Closure depth D_c + D_b (m).
    pub closure_depth: f64,
    /// CERC coefficient K.
    pub cerc_k: f64,
}
impl ShorelineEvolution {
    /// Construct a shoreline evolution model.
    pub fn new(
        position: Vec<f64>,
        dx: f64,
        diffusivity: f64,
        closure_depth: f64,
        cerc_k: f64,
    ) -> Self {
        ShorelineEvolution {
            position,
            dx,
            diffusivity,
            closure_depth,
            cerc_k,
        }
    }
    /// Advance shoreline one time step using explicit diffusion.
    pub fn advance(&mut self, dt: f64) {
        let n = self.position.len();
        if n < 3 {
            return;
        }
        let mut new_pos = self.position.clone();
        for (i, np_i) in new_pos.iter_mut().enumerate().take(n - 1).skip(1) {
            let d2y = (self.position[i + 1] - 2.0 * self.position[i] + self.position[i - 1])
                / (self.dx * self.dx);
            *np_i += dt * self.diffusivity * d2y;
        }
        self.position = new_pos;
    }
    /// Maximum stable time step (CFL condition): Δt = Δx² / (2 K_2).
    pub fn max_dt(&self) -> f64 {
        if self.diffusivity < 1e-300 {
            return f64::INFINITY;
        }
        self.dx * self.dx / (2.0 * self.diffusivity)
    }
    /// Longshore drift (CERC): Q_ls = K H_b^{5/2} sin(2α_b).
    ///
    /// `h_b` breaking wave height (m), `alpha_b` wave angle at breaking (rad).
    pub fn longshore_drift(&self, h_b: f64, alpha_b: f64) -> f64 {
        let _g = G;
        let factor = 0.5644;
        self.cerc_k * factor * h_b.powf(2.5) * (2.0 * alpha_b).sin()
    }
    /// Net shoreline advance in a single time step from a point source (m).
    pub fn point_source_advance(&self, q_source: f64, dy_strip: f64, dt: f64) -> f64 {
        if self.closure_depth < 1e-300 || dy_strip < 1e-300 {
            return 0.0;
        }
        q_source * dt / (self.closure_depth * dy_strip)
    }
}
/// Log-normal grain size distribution.
///
/// Characterises a mixed-grain-size bed for fractional transport computations.
pub struct GrainSizeDistribution {
    /// Median grain diameter d₅₀ (m).
    pub d50: f64,
    /// Geometric standard deviation σ_g (−).
    pub sigma_g: f64,
}
impl GrainSizeDistribution {
    /// Construct a log-normal GSD.
    pub fn new(d50: f64, sigma_g: f64) -> Self {
        GrainSizeDistribution { d50, sigma_g }
    }
    /// Grain diameter at percentile p ∈ \[0, 1\] (log-normal).
    pub fn diameter_at_percentile(&self, p: f64) -> f64 {
        use std::f64::consts::SQRT_2;
        let p_clamp = p.clamp(0.001, 0.999);
        let t = (2.0 * p_clamp - 1.0).abs();
        let sign = if p_clamp > 0.5 { 1.0 } else { -1.0 };
        let sigma_ln = self.sigma_g.ln().abs();
        let z = sign * SQRT_2 * erfinv_approx(t);
        self.d50 * (sigma_ln * z).exp()
    }
    /// d₁₆ (16th percentile), roughly d50 / σ_g.
    pub fn d16(&self) -> f64 {
        if self.sigma_g < 1e-300 {
            return self.d50;
        }
        self.d50 / self.sigma_g
    }
    /// d₈₄ (84th percentile), roughly d50 * σ_g.
    pub fn d84(&self) -> f64 {
        self.d50 * self.sigma_g
    }
    /// Sorting coefficient ψ = (d₈₄/d₁₆)^0.5.
    pub fn sorting_coeff(&self) -> f64 {
        let d84 = self.d84();
        let d16 = self.d16();
        if d16 < 1e-300 {
            return 1.0;
        }
        (d84 / d16).sqrt()
    }
}
/// Sediment-water mixture density model.
///
/// Computes the bulk density of a sediment-water mixture given
/// volumetric concentrations of multiple grain-size fractions.
pub struct MixtureDensity {
    /// Sediment density ρ_s (kg/m³).
    pub rho_s: f64,
    /// Fluid density ρ_f (kg/m³).
    pub rho_f: f64,
}
impl MixtureDensity {
    /// Construct a mixture density model.
    pub fn new(rho_s: f64, rho_f: f64) -> Self {
        MixtureDensity { rho_s, rho_f }
    }
    /// Mixture density ρ_m = ρ_f + C (ρ_s − ρ_f) for total volumetric concentration C.
    pub fn density(&self, concentration: f64) -> f64 {
        self.rho_f + concentration * (self.rho_s - self.rho_f)
    }
    /// Bulk density for a multi-fraction suspension.
    ///
    /// `fractions` is a slice of (volumetric concentration, grain density) pairs.
    pub fn multi_fraction_density(&self, fractions: &[(f64, f64)]) -> f64 {
        let c_total: f64 = fractions.iter().map(|(c, _)| c).sum();
        let rho_sed: f64 = if c_total > 1e-300 {
            fractions.iter().map(|(c, rho)| c * rho).sum::<f64>() / c_total
        } else {
            self.rho_s
        };
        self.rho_f + c_total * (rho_sed - self.rho_f)
    }
    /// Submerged specific gravity (s − 1) = (ρ_s − ρ_f) / ρ_f.
    pub fn submerged_specific_gravity(&self) -> f64 {
        (self.rho_s - self.rho_f) / self.rho_f.max(1e-300)
    }
}
/// Grain-size sorting and armoring model for mixed-grain-size beds.
///
/// Tracks multiple grain-size fractions with fractional transport and
/// bed surface armouring (hiding/exposure effects).
pub struct ArmoringModel {
    /// Number of grain-size fractions.
    pub n_fractions: usize,
    /// Grain diameter of each fraction d_i (m).
    pub diameters: Vec<f64>,
    /// Volume fraction of each class in subsurface p_i (−).
    pub subsurface_fractions: Vec<f64>,
    /// Volume fraction of each class in surface layer F_i (−).
    pub surface_fractions: Vec<f64>,
    /// Active layer thickness L_a (m).
    pub active_layer: f64,
}
impl ArmoringModel {
    /// Construct an armoring model with uniform subsurface fractions.
    pub fn new(diameters: Vec<f64>, active_layer: f64) -> Self {
        let n = diameters.len();
        let frac = if n > 0 { 1.0 / n as f64 } else { 0.0 };
        ArmoringModel {
            n_fractions: n,
            subsurface_fractions: vec![frac; n],
            surface_fractions: vec![frac; n],
            diameters,
            active_layer,
        }
    }
    /// Geometric mean diameter d_m of the surface layer.
    pub fn surface_mean_diameter(&self) -> f64 {
        if self.surface_fractions.is_empty() {
            return 0.0;
        }
        let log_dm: f64 = self
            .surface_fractions
            .iter()
            .zip(self.diameters.iter())
            .map(|(f, d)| f * d.max(1e-300).ln())
            .sum();
        log_dm.exp()
    }
    /// Hiding-exposure correction factor (Egiazaroff, 1965).
    ///
    /// ξ_i = \[log(19) / log(19 d_i/d_m)\]² (ratio of critical Shields for fraction i to mean).
    pub fn hiding_factor(&self, frac_idx: usize) -> f64 {
        let d_m = self.surface_mean_diameter();
        if d_m < 1e-300 || frac_idx >= self.n_fractions {
            return 1.0;
        }
        let d_i = self.diameters[frac_idx];
        let log_arg = 19.0 * d_i / d_m;
        if log_arg <= 1.0 {
            return 1.0;
        }
        let log19 = 19.0_f64.ln();
        (log19 / log_arg.ln()).powi(2)
    }
    /// Fractional transport rate q_i using MPM with hiding correction.
    ///
    /// Returns transport rate (m²/s) for fraction `i` under bed shear τ_b.
    pub fn fractional_transport(&self, frac_idx: usize, tau_b: f64, theta_cr_ref: f64) -> f64 {
        if frac_idx >= self.n_fractions {
            return 0.0;
        }
        let d_i = self.diameters[frac_idx];
        let xi_i = self.hiding_factor(frac_idx);
        let theta_cr_i = theta_cr_ref * xi_i;
        let fi = self.surface_fractions[frac_idx];
        let q_i = mpm_bedload(
            shields_stress(tau_b, d_i, RHO_SEDIMENT, RHO_WATER),
            theta_cr_i,
            d_i,
            RHO_SEDIMENT,
            RHO_WATER,
        );
        fi * q_i
    }
    /// Update surface layer fractions from active layer exchange.
    ///
    /// `delta_z` is the bed change (negative = erosion). Simplistic model.
    pub fn update_surface(&mut self, delta_z: f64) {
        if delta_z < 0.0 {
            let mix = (-delta_z / self.active_layer.max(1e-300)).min(1.0);
            for i in 0..self.n_fractions {
                self.surface_fractions[i] =
                    (1.0 - mix) * self.surface_fractions[i] + mix * self.subsurface_fractions[i];
            }
        }
        let total: f64 = self.surface_fractions.iter().sum();
        if total > 1e-300 {
            for f in &mut self.surface_fractions {
                *f /= total;
            }
        }
    }
}
/// Bedform (dune) geometry and migration model.
///
/// Implements the van Rijn (1984) bedform predictor and dune migration
/// velocities based on sediment transport theory.
pub struct DuneMigration {
    /// Dune height H_d (m).
    pub height: f64,
    /// Dune length L_d (m).
    pub length: f64,
    /// Bed-load transport rate q_b (m²/s).
    pub bed_load: f64,
    /// Porosity of bed material n (−).
    pub porosity: f64,
    /// Current bed elevation profile η (m).
    pub bed_elevation: Vec<f64>,
    /// Grid spacing Δx (m).
    pub dx: f64,
}
impl DuneMigration {
    /// Construct a dune migration model.
    pub fn new(height: f64, length: f64, bed_load: f64, porosity: f64, dx: f64) -> Self {
        DuneMigration {
            height,
            length,
            bed_load,
            porosity,
            bed_elevation: Vec::new(),
            dx,
        }
    }
    /// van Rijn dune height prediction: H_d = 0.11 d (H/d)^0.3 (1-e^{-0.5T_*}) (25-T_*).
    ///
    /// T_* = (θ - θ_cr) / θ_cr, d grain diameter (m), H water depth (m).
    pub fn van_rijn_height(&self, depth: f64, d50: f64, theta: f64, theta_cr: f64) -> f64 {
        let t_star = if theta_cr > 1e-300 {
            (theta - theta_cr).max(0.0) / theta_cr
        } else {
            0.0
        };
        if t_star >= 25.0 || d50 < 1e-300 {
            return 0.0;
        }
        0.11 * d50 * (depth / d50).powf(0.3) * (1.0 - (-0.5 * t_star).exp()) * (25.0 - t_star)
    }
    /// van Rijn dune length: L_d ≈ 7.3 H.
    pub fn van_rijn_length(&self, depth: f64) -> f64 {
        7.3 * depth
    }
    /// Dune migration velocity c_d = q_b / ((1-n) H_d / 2).
    pub fn migration_velocity(&self) -> f64 {
        let denom = (1.0 - self.porosity) * self.height / 2.0;
        if denom < 1e-300 {
            return 0.0;
        }
        self.bed_load / denom
    }
    /// Advance bed profile by one time step using kinematic wave equation.
    pub fn advance_profile(&mut self, dt: f64) {
        let n = self.bed_elevation.len();
        if n < 2 {
            return;
        }
        let c = self.migration_velocity();
        let mut new_bed = self.bed_elevation.clone();
        for (i, nb_i) in new_bed.iter_mut().enumerate().take(n).skip(1) {
            *nb_i = self.bed_elevation[i]
                - dt * c * (self.bed_elevation[i] - self.bed_elevation[i - 1]) / self.dx;
        }
        self.bed_elevation = new_bed;
    }
    /// Dune form drag coefficient C_f_dune (van Rijn).
    ///
    /// C_f = 0.5 (H_d/H) (H_d/L_d).
    pub fn form_drag_coeff(&self, depth: f64) -> f64 {
        if depth < 1e-300 || self.length < 1e-300 {
            return 0.0;
        }
        0.5 * (self.height / depth) * (self.height / self.length)
    }
}
/// Multi-fraction Exner equation for bed evolution.
///
/// ∂η/∂t = -1/(1−n) Σ_i ∂q_i/∂x
///
/// where q_i is the bed-load flux of fraction i.
pub struct MultiFractionExner {
    /// Bed elevation grid η (m).
    pub eta: Vec<f64>,
    /// Grid spacing Δx (m).
    pub dx: f64,
    /// Bed porosity n (−).
    pub porosity: f64,
    /// Number of grain-size fractions.
    pub n_fractions: usize,
}
impl MultiFractionExner {
    /// Construct a multi-fraction Exner model.
    pub fn new(eta: Vec<f64>, dx: f64, porosity: f64, n_fractions: usize) -> Self {
        MultiFractionExner {
            eta,
            dx,
            porosity,
            n_fractions,
        }
    }
    /// Update bed elevation from fractional bed-load fluxes.
    ///
    /// `q_fractions[i][j]` = flux of fraction i at grid node j (m²/s).
    pub fn update(&mut self, q_fractions: &[Vec<f64>], dt: f64) {
        let n = self.eta.len();
        let scale = 1.0 / (1.0 - self.porosity);
        for j in 1..n {
            let mut total_dq = 0.0;
            for q_frac in q_fractions {
                if j < q_frac.len() {
                    total_dq += (q_frac[j] - q_frac[j - 1]) / self.dx;
                }
            }
            self.eta[j] -= dt * scale * total_dq;
        }
    }
    /// Total volume change (m²) from initial elevation.
    pub fn volume_change(&self, eta_0: &[f64]) -> f64 {
        self.eta
            .iter()
            .zip(eta_0)
            .map(|(e, e0)| (e - e0) * self.dx)
            .sum()
    }
}
/// Extended Shields parameter with slope correction.
///
/// On a sloping bed the critical Shields parameter is modified by the
/// ratio of angle of repose to local bed slope (Soulsby & Damgaard, 2005).
pub struct ShieldsSlope {
    /// Base critical Shields parameter θ_cr0 (horizontal bed).
    pub theta_cr0: f64,
    /// Angle of repose φ_r (radians).
    pub angle_of_repose: f64,
}
impl ShieldsSlope {
    /// Construct a slope-corrected Shields parameter.
    pub fn new(theta_cr0: f64, angle_of_repose_deg: f64) -> Self {
        ShieldsSlope {
            theta_cr0,
            angle_of_repose: angle_of_repose_deg.to_radians(),
        }
    }
    /// Critical Shields parameter on a bed sloping at angle β (radians).
    ///
    /// θ_cr(β) = θ_cr0 * (cos β − sin β / tan φ_r) for down-slope.
    pub fn critical_shields_downslope(&self, beta: f64) -> f64 {
        let tan_phi = self.angle_of_repose.tan().max(1e-300);
        self.theta_cr0 * (beta.cos() - beta.sin() / tan_phi).max(0.0)
    }
    /// Critical Shields for up-slope transport.
    ///
    /// θ_cr_up = θ_cr0 * (cos β + sin β / tan φ_r).
    pub fn critical_shields_upslope(&self, beta: f64) -> f64 {
        let tan_phi = self.angle_of_repose.tan().max(1e-300);
        self.theta_cr0 * (beta.cos() + beta.sin() / tan_phi)
    }
    /// Slope correction factor κ = θ_cr(β) / θ_cr0.
    pub fn slope_correction(&self, beta: f64) -> f64 {
        let tan_phi = self.angle_of_repose.tan().max(1e-300);
        (beta.cos() - beta.sin() / tan_phi).max(0.0)
    }
}
/// Turbulent sediment diffusivity model.
pub struct TurbulentDiffusion {
    /// Schmidt number Sc (ratio of momentum to mass diffusivity).
    pub schmidt_number: f64,
    /// Turbulent kinematic viscosity ν_t (m²/s).
    pub turbulent_viscosity: f64,
}
impl TurbulentDiffusion {
    /// Construct a turbulent diffusion model.
    pub fn new(schmidt_number: f64, turbulent_viscosity: f64) -> Self {
        TurbulentDiffusion {
            schmidt_number,
            turbulent_viscosity,
        }
    }
    /// Turbulent sediment diffusivity ε_s = ν_t / Sc.
    pub fn diffusivity(&self) -> f64 {
        self.turbulent_viscosity / self.schmidt_number.max(1e-10)
    }
    /// Turbulent diffusion flux −ε_s * ∇C.
    pub fn diffusion_flux(&self, grad_c: [f64; 3]) -> [f64; 3] {
        let eps = self.diffusivity();
        scale3(grad_c, -eps)
    }
    /// Peclet number Pe = u*L / ε_s for advection vs. diffusion.
    pub fn peclet_number(&self, velocity: f64, length_scale: f64) -> f64 {
        let eps = self.diffusivity();
        if eps < 1e-300 {
            return 0.0;
        }
        velocity * length_scale / eps
    }
}
/// Flocculation model for cohesive sediment (clay, silt).
///
/// Implements the Krone (1962) and Winterwerp (1998) flocculation framework
/// with aggregate growth and breakup kinetics.
pub struct FlocculationModel {
    /// Primary particle diameter d_p (m).
    pub primary_diameter: f64,
    /// Maximum floc diameter d_f_max (m).
    pub floc_diameter_max: f64,
    /// Fractal dimension n_f (typically 2.0 for cohesive sediment).
    pub fractal_dimension: f64,
    /// Aggregation rate coefficient k_a (m³/s).
    pub aggregation_rate: f64,
    /// Breakup rate coefficient k_b (s^{n_b}).
    pub breakup_rate: f64,
    /// Breakup exponent n_b.
    pub breakup_exponent: f64,
}
impl FlocculationModel {
    /// Construct a flocculation model with default parameters for fine cohesive sediment.
    pub fn new(primary_diameter: f64) -> Self {
        FlocculationModel {
            primary_diameter,
            floc_diameter_max: primary_diameter * 1000.0,
            fractal_dimension: 2.0,
            aggregation_rate: 1e-15,
            breakup_rate: 1e-4,
            breakup_exponent: 0.5,
        }
    }
    /// Number of primary particles per floc: n_p = (d_f/d_p)^{n_f}.
    pub fn particles_per_floc(&self, floc_diameter: f64) -> f64 {
        if self.primary_diameter < 1e-300 {
            return 1.0;
        }
        (floc_diameter / self.primary_diameter).powf(self.fractal_dimension)
    }
    /// Effective floc density: ρ_f = ρ_w + (ρ_s - ρ_w) (d_p/d_f)^{3-n_f}.
    pub fn floc_density(&self, floc_diameter: f64, rho_s: f64, rho_w: f64) -> f64 {
        if floc_diameter < 1e-300 {
            return rho_s;
        }
        let ratio = (self.primary_diameter / floc_diameter).powf(3.0 - self.fractal_dimension);
        rho_w + (rho_s - rho_w) * ratio
    }
    /// Floc settling velocity using Stokes law with effective floc density.
    pub fn floc_settling_velocity(
        &self,
        floc_diameter: f64,
        rho_s: f64,
        rho_w: f64,
        mu: f64,
    ) -> f64 {
        let rho_f = self.floc_density(floc_diameter, rho_s, rho_w);
        settling_velocity_stokes(floc_diameter, rho_f, rho_w, mu)
    }
    /// Aggregation rate G_a = k_a G C² where G is shear rate and C is concentration.
    pub fn aggregation_rate_fn(&self, shear_rate: f64, concentration: f64) -> f64 {
        self.aggregation_rate * shear_rate * concentration * concentration
    }
    /// Breakup rate G_b = k_b G^{n_b} C.
    pub fn breakup_rate_fn(&self, shear_rate: f64, concentration: f64) -> f64 {
        self.breakup_rate * shear_rate.powf(self.breakup_exponent) * concentration
    }
    /// Equilibrium floc size from aggregation-breakup balance.
    pub fn equilibrium_diameter(&self, shear_rate: f64) -> f64 {
        if shear_rate < 1e-300 {
            return self.floc_diameter_max;
        }
        let ratio = self.breakup_rate / (self.aggregation_rate * shear_rate).max(1e-300);
        let d_eq = self.primary_diameter * ratio.powf(1.0 / self.fractal_dimension);
        d_eq.clamp(self.primary_diameter, self.floc_diameter_max)
    }
}
/// A fluid particle augmented with sediment transport quantities.
#[derive(Clone, Debug)]
pub struct SedimentParticle {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Suspended sediment volumetric concentration (−).
    pub concentration: f64,
    /// Settling velocity (m/s, positive downward).
    pub settling_velocity: f64,
    /// Bed-load flux magnitude (m²/s).
    pub bed_load_flux: f64,
    /// Grain diameter (m).
    pub grain_diameter: f64,
    /// SPH smoothing length (m).
    pub smoothing_length: f64,
}
impl SedimentParticle {
    /// Create a new sediment particle.
    pub fn new(pos: [f64; 3], vel: [f64; 3], density: f64, grain_diameter: f64) -> Self {
        let mu = RHO_WATER * NU_WATER;
        let vs = settling_velocity_stokes(grain_diameter, RHO_SEDIMENT, RHO_WATER, mu);
        SedimentParticle {
            pos,
            vel,
            density,
            pressure: 0.0,
            concentration: 0.0,
            settling_velocity: vs,
            bed_load_flux: 0.0,
            grain_diameter,
            smoothing_length: grain_diameter * 10.0,
        }
    }
    /// Update the settling velocity using the current grain diameter.
    pub fn update_settling_velocity(&mut self) {
        let mu = self.density * NU_WATER;
        self.settling_velocity =
            settling_velocity_stokes(self.grain_diameter, RHO_SEDIMENT, self.density, mu);
    }
    /// Kinetic energy of the particle.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.density * dot3(self.vel, self.vel)
    }
}
/// CERC formula for longshore sediment transport.
///
/// Q_ls = K * H_b^{5/2} / (8 * ρ_w * g^{1/2} * (tan β))
///
/// where H_b is breaking wave height, K is empirical coefficient, tan β is
/// beach slope.
pub struct CoastalMorphology {
    /// CERC coefficient K (typically 0.39).
    pub cerc_coeff: f64,
    /// Water density (kg/m³).
    pub rho_water: f64,
    /// Beach slope tangent tan β.
    pub beach_slope: f64,
    /// Shoreline position x (m).
    pub shoreline: Vec<f64>,
    /// Along-shore grid spacing dy (m).
    pub dy: f64,
}
impl CoastalMorphology {
    /// Construct a coastal morphology model.
    pub fn new(cerc_coeff: f64, beach_slope: f64, shoreline: Vec<f64>, dy: f64) -> Self {
        CoastalMorphology {
            cerc_coeff,
            rho_water: RHO_WATER,
            beach_slope,
            shoreline,
            dy,
        }
    }
    /// Longshore sediment transport rate Q_ls (m³/s) from breaking wave height H_b.
    pub fn longshore_transport(&self, h_b: f64) -> f64 {
        if h_b < 0.0 || self.beach_slope < 1e-300 {
            return 0.0;
        }
        self.cerc_coeff * h_b.powf(2.5) / (8.0 * self.rho_water * G.sqrt() * self.beach_slope)
    }
    /// Update shoreline position from gradient of longshore transport.
    ///
    /// ∂x_s/∂t = -1/(D_c + D_b) * ∂Q/∂y
    ///
    /// `q_ls` transport at each along-shore cell, `dt` time step (s),
    /// `closure_depth` D_c + D_b (m).
    pub fn update_shoreline(&mut self, q_ls: &[f64], dt: f64, closure_depth: f64) {
        let n = self.shoreline.len();
        if n == 0 || q_ls.len() != n || closure_depth < 1e-300 {
            return;
        }
        for i in 1..n - 1 {
            let dq_dy = (q_ls[i + 1] - q_ls[i - 1]) / (2.0 * self.dy);
            self.shoreline[i] -= dt * dq_dy / closure_depth;
        }
    }
}
/// SPH-based bed-load flux calculator for unstructured particle simulations.
pub struct SphBedLoadFlux {
    /// Grain diameter (m).
    pub d50: f64,
    /// Critical Shields parameter θ_cr.
    pub theta_cr: f64,
    /// MPM coefficient (default 8.0).
    pub mpm_coeff: f64,
    /// Sediment density (kg/m³).
    pub rho_s: f64,
    /// Fluid density (kg/m³).
    pub rho_f: f64,
}
impl SphBedLoadFlux {
    /// Construct an SPH bed-load flux calculator.
    pub fn new(d50: f64, theta_cr: f64, rho_s: f64, rho_f: f64) -> Self {
        SphBedLoadFlux {
            d50,
            theta_cr,
            mpm_coeff: 8.0,
            rho_s,
            rho_f,
        }
    }
    /// Compute bed-load flux vector at particle i from local bed shear vector.
    pub fn flux_vector(&self, tau_b_vec: [f64; 3]) -> [f64; 3] {
        let tau_mag = len3(tau_b_vec);
        let theta = shields_stress(tau_mag, self.d50, self.rho_s, self.rho_f);
        if theta <= self.theta_cr {
            return [0.0; 3];
        }
        let q = mpm_bedload(theta, self.theta_cr, self.d50, self.rho_s, self.rho_f);
        if tau_mag < 1e-300 {
            return [0.0; 3];
        }
        scale3(
            [
                tau_b_vec[0] / tau_mag,
                tau_b_vec[1] / tau_mag,
                tau_b_vec[2] / tau_mag,
            ],
            q,
        )
    }
    /// Bed-load reference concentration at a near-bed SPH particle.
    ///
    /// Van Rijn reference concentration: C_a = 0.015 (d/a) T*^{1.5} / D*^{0.3}.
    pub fn reference_concentration(&self, theta: f64, nu: f64) -> f64 {
        let t_star = ((theta - self.theta_cr) / self.theta_cr.max(1e-300)).max(0.0);
        let s = self.rho_s / self.rho_f;
        let d_star = self.d50 * ((s - 1.0) * G / (nu * nu)).powf(1.0 / 3.0);
        if d_star < 1e-300 {
            return 0.0;
        }
        let a = self.d50 * 0.5;
        0.015 * (self.d50 / a) * t_star.powf(1.5) / d_star.powf(0.3)
    }
}
/// SPH particle with full suspended sediment transport state.
#[derive(Clone, Debug)]
pub struct SuspendedSedimentParticle {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Suspended sediment concentration (volumetric, −).
    pub concentration: f64,
    /// Grain diameter d₅₀ (m).
    pub d50: f64,
    /// Settling velocity w_s (m/s).
    pub ws: f64,
    /// Turbulent diffusivity ε_s (m²/s).
    pub diffusivity: f64,
    /// SPH smoothing length (m).
    pub h: f64,
    /// Is this a near-bed particle?
    pub is_near_bed: bool,
}
impl SuspendedSedimentParticle {
    /// Construct a suspended sediment particle.
    pub fn new(pos: [f64; 3], vel: [f64; 3], density: f64, d50: f64, h: f64) -> Self {
        let mu = RHO_WATER * NU_WATER;
        let ws = settling_velocity_stokes(d50, RHO_SEDIMENT, RHO_WATER, mu);
        SuspendedSedimentParticle {
            pos,
            vel,
            density,
            pressure: 0.0,
            concentration: 0.0,
            d50,
            ws,
            diffusivity: 1e-5,
            h,
            is_near_bed: false,
        }
    }
    /// Rouse number Ro = w_s / (κ u*).
    pub fn rouse_number(&self, shear_velocity: f64) -> f64 {
        let kappa = 0.41;
        if shear_velocity < 1e-300 {
            return 0.0;
        }
        self.ws / (kappa * shear_velocity)
    }
    /// Transport mode based on Rouse number.
    pub fn transport_mode(&self, shear_velocity: f64) -> SedimentTransportMode {
        let ro = self.rouse_number(shear_velocity);
        if ro > 2.5 {
            SedimentTransportMode::BedLoad
        } else if ro > 1.2 {
            SedimentTransportMode::SaltationLoad
        } else if ro > 0.8 {
            SedimentTransportMode::SuspendedLoad
        } else {
            SedimentTransportMode::WashLoad
        }
    }
}
/// Morphodynamic time-step controller.
///
/// Implements the Courant condition for bed-level change and limits the
/// morphodynamic time step to ensure stability of the Exner equation.
pub struct MorphodynamicTimestep {
    /// Maximum allowed relative bed change per step (−).
    pub max_rel_change: f64,
    /// Grid spacing Δx (m).
    pub dx: f64,
    /// Morphological acceleration factor f_morph.
    pub f_morph: f64,
}
impl MorphodynamicTimestep {
    /// Construct a morphodynamic time-step controller.
    pub fn new(max_rel_change: f64, dx: f64, f_morph: f64) -> Self {
        MorphodynamicTimestep {
            max_rel_change,
            dx,
            f_morph,
        }
    }
    /// Maximum morphodynamic time step from Courant condition.
    ///
    /// Δt_morph ≤ (1−n) max_rel_change η_max Δx / (f_morph max|q_b|)
    pub fn max_dt(&self, eta: &[f64], q_b: &[f64], porosity: f64) -> f64 {
        let eta_max = eta
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max)
            .abs()
            .max(1e-300);
        let q_max = q_b
            .iter()
            .cloned()
            .map(|q| q.abs())
            .fold(0.0_f64, f64::max)
            .max(1e-300);
        (1.0 - porosity) * self.max_rel_change * eta_max * self.dx / (self.f_morph * q_max)
    }
    /// Acoustic CFL criterion for bed-load wave propagation.
    ///
    /// The bed-load celerity c_b = dq_b/dη; estimated as q_b / η.
    pub fn bedload_cfl_dt(&self, eta: f64, q_b: f64, porosity: f64) -> f64 {
        let c_b = if eta.abs() > 1e-300 {
            q_b.abs() / (eta.abs() * (1.0 - porosity))
        } else {
            0.0
        };
        if c_b < 1e-300 {
            return f64::INFINITY;
        }
        self.dx / (self.f_morph * c_b)
    }
}
