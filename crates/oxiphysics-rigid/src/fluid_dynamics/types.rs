use super::functions::*;
// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// Advanced VIV lock-in prediction using Strouhal-Reynolds correlation.
#[derive(Debug, Clone, Copy)]
pub struct VIVPredictor {
    /// Cylinder diameter \[m\].
    pub diameter: f64,
    /// Fluid density \[kg/m³\].
    pub rho: f64,
    /// Fluid kinematic viscosity \[m²/s\].
    pub nu: f64,
    /// Structure natural frequency \[Hz\].
    pub f_natural: f64,
    /// Structural mass ratio m* = m / (ρ D² L).
    pub mass_ratio: f64,
    /// Damping ratio ζ.
    pub damping_ratio: f64,
}
impl VIVPredictor {
    /// Create a VIV predictor.
    pub fn new(
        diameter: f64,
        rho: f64,
        nu: f64,
        f_natural: f64,
        mass_ratio: f64,
        damping_ratio: f64,
    ) -> Self {
        Self {
            diameter,
            rho,
            nu,
            f_natural,
            mass_ratio,
            damping_ratio,
        }
    }
    /// Strouhal number from Reynolds number (Williamson correlation for cylinder).
    pub fn strouhal_from_re(&self, re: f64) -> f64 {
        if re < 250.0 {
            0.2
        } else if re < 2e5 {
            0.198 * (1.0 - 19.7 / re)
        } else {
            0.27
        }
    }
    /// Reynolds number at flow velocity U.
    pub fn reynolds(&self, u: f64) -> f64 {
        u * self.diameter / self.nu
    }
    /// Shedding frequency at U \[Hz\].
    pub fn shedding_frequency(&self, u: f64) -> f64 {
        let re = self.reynolds(u);
        let st = self.strouhal_from_re(re);
        st * u / self.diameter
    }
    /// Reduced velocity Ur = U / (f_n * D).
    pub fn reduced_velocity(&self, u: f64) -> f64 {
        u / (self.f_natural * self.diameter)
    }
    /// Synchronization (lock-in) bandwidth: 4.5 ≤ Ur ≤ 8 (typical for circular cylinder).
    pub fn is_lock_in(&self, u: f64) -> bool {
        let ur = self.reduced_velocity(u);
        (4.5..=8.0).contains(&ur)
    }
    /// Scruton number Sc = 2π m* ζ (mass-damping parameter).
    pub fn scruton_number(&self) -> f64 {
        2.0 * PI * self.mass_ratio * self.damping_ratio
    }
    /// Maximum amplitude ratio A/D from Griffith correlation: A/D ≈ 1.3/√Sc.
    pub fn max_amplitude_ratio(&self) -> f64 {
        let sc = self.scruton_number();
        if sc < 1e-10 {
            return 10.0;
        }
        1.3 / sc.sqrt()
    }
    /// Peak transverse force coefficient CL_peak (Gopalkrishnan model).
    pub fn peak_cl(&self) -> f64 {
        let a_d = self.max_amplitude_ratio().min(1.0);
        0.7 * (1.0 - a_d).max(0.0)
    }
}
/// Kelvin-Helmholtz instability linear analysis at a shear interface.
///
/// Two fluid layers with different velocities and densities.
#[derive(Debug, Clone, Copy)]
pub struct KelvinHelmholtz {
    /// Velocity of upper layer U₁ \[m/s\].
    pub u1: f64,
    /// Velocity of lower layer U₂ \[m/s\].
    pub u2: f64,
    /// Density of upper layer ρ₁ \[kg/m³\].
    pub rho1: f64,
    /// Density of lower layer ρ₂ \[kg/m³\].
    pub rho2: f64,
    /// Surface tension σ \[N/m\] (stabilizing).
    pub surface_tension: f64,
    /// Gravitational acceleration g \[m/s²\].
    pub gravity: f64,
}
impl KelvinHelmholtz {
    /// Create a Kelvin-Helmholtz analysis.
    pub fn new(u1: f64, u2: f64, rho1: f64, rho2: f64, surface_tension: f64, gravity: f64) -> Self {
        Self {
            u1,
            u2,
            rho1,
            rho2,
            surface_tension,
            gravity,
        }
    }
    /// Atwood number A = (ρ₂ - ρ₁) / (ρ₁ + ρ₂).
    pub fn atwood_number(&self) -> f64 {
        (self.rho2 - self.rho1) / (self.rho1 + self.rho2)
    }
    /// KH growth rate ω_i for wavenumber k \[1/m\].
    /// ω_i² = k² ρ₁ρ₂ (U₁-U₂)² / (ρ₁+ρ₂)² - (σk³/(ρ₁+ρ₂) + A*g*k).
    pub fn growth_rate(&self, k: f64) -> f64 {
        let rho_sum = self.rho1 + self.rho2;
        let du = self.u1 - self.u2;
        let kh_term = k * k * self.rho1 * self.rho2 * du * du / (rho_sum * rho_sum);
        let st_term = self.surface_tension * k * k * k / rho_sum;
        let grav_term = self.atwood_number() * self.gravity * k;
        let omega_i_sq = kh_term - st_term - grav_term;
        if omega_i_sq <= 0.0 {
            0.0
        } else {
            omega_i_sq.sqrt()
        }
    }
    /// Critical velocity difference (minimum ΔU for instability at wavenumber k).
    pub fn critical_velocity_difference(&self, k: f64) -> f64 {
        let rho_sum = self.rho1 + self.rho2;
        let denom = k * k * self.rho1 * self.rho2;
        if denom < 1e-14 {
            return f64::INFINITY;
        }
        let stabilizing = (self.surface_tension * k * k * k / rho_sum
            + self.atwood_number() * self.gravity * k)
            * rho_sum
            * rho_sum;
        (stabilizing / denom).max(0.0).sqrt()
    }
    /// Most unstable wavenumber k_max (peak growth rate).
    /// Approximate: k_max ≈ (g (ρ₂-ρ₁) / (2σ))^{1/2}.
    pub fn most_unstable_wavenumber(&self) -> f64 {
        if self.surface_tension < 1e-14 {
            return f64::INFINITY;
        }
        let num = self.gravity * (self.rho2 - self.rho1).abs();
        let den = 2.0 * self.surface_tension;
        (num / den).sqrt()
    }
    /// Phase velocity at wavenumber k: c = (ρ₁U₁ + ρ₂U₂) / (ρ₁+ρ₂).
    pub fn phase_velocity(&self, _k: f64) -> f64 {
        (self.rho1 * self.u1 + self.rho2 * self.u2) / (self.rho1 + self.rho2)
    }
}
/// Turbulent boundary layer analysis (log-law, integral methods).
#[derive(Debug, Clone, Copy)]
pub struct TurbulentBoundaryLayer {
    /// Free-stream velocity U_inf \[m/s\].
    pub u_inf: f64,
    /// Kinematic viscosity ν \[m²/s\].
    pub nu: f64,
    /// Distance from leading edge x \[m\].
    pub x: f64,
    /// Fluid density ρ \[kg/m³\].
    pub rho: f64,
}
impl TurbulentBoundaryLayer {
    /// Create a turbulent boundary layer.
    pub fn new(u_inf: f64, nu: f64, x: f64, rho: f64) -> Self {
        Self { u_inf, nu, x, rho }
    }
    /// Reynolds number Rex = U_inf * x / ν.
    pub fn reynolds_x(&self) -> f64 {
        self.u_inf * self.x / self.nu
    }
    /// Boundary layer thickness δ₉₉ (1/7 power law approximation).
    /// δ ≈ 0.37 x Rex^{-1/5}.
    pub fn thickness(&self) -> f64 {
        let rex = self.reynolds_x();
        if rex < 1.0 {
            return 0.0;
        }
        0.37 * self.x * rex.powf(-0.2)
    }
    /// Displacement thickness δ* for turbulent BL (1/7 power law: δ* = δ/8).
    pub fn displacement_thickness(&self) -> f64 {
        self.thickness() / 8.0
    }
    /// Momentum thickness θ for turbulent BL (1/7 power law: θ = 7δ/72).
    pub fn momentum_thickness(&self) -> f64 {
        7.0 * self.thickness() / 72.0
    }
    /// Shape factor H = δ* / θ.
    pub fn shape_factor(&self) -> f64 {
        let dt = self.displacement_thickness();
        let mt = self.momentum_thickness();
        if mt < 1e-14 {
            return 0.0;
        }
        dt / mt
    }
    /// Skin friction coefficient Cf (Schlichting formula for turbulent BL).
    /// Cf ≈ 0.0592 Rex^{-1/5}.
    pub fn skin_friction_coefficient(&self) -> f64 {
        let rex = self.reynolds_x();
        if rex < 1.0 {
            return 0.0;
        }
        0.0592 * rex.powf(-0.2)
    }
    /// Wall shear stress τ_w = 0.5 ρ U² Cf \[Pa\].
    pub fn wall_shear_stress(&self) -> f64 {
        0.5 * self.rho * self.u_inf * self.u_inf * self.skin_friction_coefficient()
    }
    /// Friction velocity u_τ = √(τ_w / ρ) \[m/s\].
    pub fn friction_velocity(&self) -> f64 {
        (self.wall_shear_stress() / self.rho).max(0.0).sqrt()
    }
    /// Log-law velocity at wall-normal distance y \[m\].
    /// u/u_τ = (1/κ) ln(y u_τ / ν) + B, κ = 0.41, B = 5.1.
    pub fn log_law_velocity(&self, y: f64) -> f64 {
        let u_tau = self.friction_velocity();
        if u_tau < 1e-14 {
            return 0.0;
        }
        let y_plus = y * u_tau / self.nu;
        if y_plus < 5.0 {
            y_plus * u_tau
        } else {
            let kappa = 0.41;
            let b = 5.1;
            (y_plus.ln() / kappa + b) * u_tau
        }
    }
    /// Viscous sublayer thickness y_v = 5 ν / u_τ \[m\].
    pub fn viscous_sublayer_thickness(&self) -> f64 {
        let u_tau = self.friction_velocity();
        if u_tau < 1e-14 {
            return f64::INFINITY;
        }
        5.0 * self.nu / u_tau
    }
    /// Total drag coefficient over plate length x.
    pub fn drag_coefficient(&self) -> f64 {
        let rex = self.reynolds_x();
        if rex < 1.0 {
            return 0.0;
        }
        0.074 * rex.powf(-0.2)
    }
}
/// Vortex-induced vibration (VIV) model for bluff bodies.
#[derive(Debug, Clone, Copy)]
pub struct VortexInducedVibration {
    /// Strouhal number St (typically 0.2 for circular cylinder).
    pub strouhal: f64,
    /// Diameter of cylinder \[m\].
    pub diameter: f64,
    /// Fluid density \[kg/m³\].
    pub rho: f64,
    /// Lift coefficient at lock-in CL.
    pub cl: f64,
    /// Lock-in bandwidth ratio (±Δv/v_s).
    pub lock_in_bandwidth: f64,
}
impl VortexInducedVibration {
    /// Create a VIV model.
    pub fn new(strouhal: f64, diameter: f64, rho: f64, cl: f64, lock_in_bandwidth: f64) -> Self {
        Self {
            strouhal,
            diameter,
            rho,
            cl,
            lock_in_bandwidth,
        }
    }
    /// Vortex shedding frequency at flow velocity v \[Hz\].
    pub fn shedding_frequency(&self, v: f64) -> f64 {
        self.strouhal * v / self.diameter
    }
    /// Check if lock-in occurs for a structure with natural frequency f_n \[Hz\].
    pub fn is_lock_in(&self, v: f64, f_natural: f64) -> bool {
        let f_shed = self.shedding_frequency(v);
        let ratio = f_shed / f_natural;
        (ratio - 1.0).abs() <= self.lock_in_bandwidth
    }
    /// Transverse VIV force amplitude \[N/m\] at velocity v \[m/s\].
    pub fn transverse_force_per_length(&self, v: f64) -> f64 {
        0.5 * self.rho * v * v * self.diameter * self.cl
    }
    /// Instantaneous transverse force at time t \[N/m\].
    pub fn force_per_length(&self, v: f64, t: f64) -> f64 {
        let omega = 2.0 * PI * self.shedding_frequency(v);
        let amp = self.transverse_force_per_length(v);
        amp * (omega * t).sin()
    }
    /// Reduced velocity Ur = v / (f_n * D).
    pub fn reduced_velocity(&self, v: f64, f_natural: f64) -> f64 {
        v / (f_natural * self.diameter)
    }
}
/// Vortex panel method for lifting airfoil (with Kutta condition).
#[derive(Debug, Clone)]
pub struct VortexPanelMethod {
    /// Panels defining the airfoil surface.
    pub panels: Vec<Panel>,
    /// Free-stream velocity \[m/s\].
    pub u_inf: f64,
    /// Angle of attack \[rad\].
    pub alpha: f64,
    /// Vortex strengths γ at each panel after solve.
    pub gamma: Vec<f64>,
}
impl VortexPanelMethod {
    /// Create a vortex panel method.
    pub fn new(panels: Vec<Panel>, u_inf: f64, alpha: f64) -> Self {
        let n = panels.len();
        Self {
            panels,
            u_inf,
            alpha,
            gamma: vec![0.0; n],
        }
    }
    /// Total circulation Γ = Σ γᵢ * Δsᵢ.
    pub fn total_circulation(&self) -> f64 {
        self.panels
            .iter()
            .zip(self.gamma.iter())
            .map(|(p, &g)| g * p.length())
            .sum()
    }
    /// Kutta-Joukowski lift per unit span: L = ρ * V * Γ.
    pub fn lift_per_span(&self, rho: f64) -> f64 {
        rho * self.u_inf * self.total_circulation()
    }
    /// Lift coefficient Cl = L / (0.5 ρ V² c), where c = chord.
    pub fn lift_coefficient(&self, rho: f64, chord: f64) -> f64 {
        let q = 0.5 * rho * self.u_inf * self.u_inf;
        if q < 1e-14 || chord < 1e-14 {
            return 0.0;
        }
        self.lift_per_span(rho) / (q * chord)
    }
    /// Assign uniform vortex strength (simplified, for testing).
    pub fn set_uniform_gamma(&mut self, g: f64) {
        for v in &mut self.gamma {
            *v = g;
        }
    }
}
/// A single flat panel for 2-D panel method (source/vortex/doublet).
#[derive(Debug, Clone)]
pub struct Panel {
    /// Start point of panel (x1, y1).
    pub p1: [f64; 2],
    /// End point of panel (x2, y2).
    pub p2: [f64; 2],
    /// Panel strength (source strength σ, vortex γ, or doublet μ).
    pub strength: f64,
}
impl Panel {
    /// Create a new panel.
    pub fn new(p1: [f64; 2], p2: [f64; 2]) -> Self {
        Self {
            p1,
            p2,
            strength: 0.0,
        }
    }
    /// Panel length.
    pub fn length(&self) -> f64 {
        let dx = self.p2[0] - self.p1[0];
        let dy = self.p2[1] - self.p1[1];
        (dx * dx + dy * dy).sqrt()
    }
    /// Panel center point (collocation point).
    pub fn center(&self) -> [f64; 2] {
        [
            0.5 * (self.p1[0] + self.p2[0]),
            0.5 * (self.p1[1] + self.p2[1]),
        ]
    }
    /// Outward unit normal (perpendicular to panel, pointing outward).
    pub fn normal(&self) -> [f64; 2] {
        let dx = self.p2[0] - self.p1[0];
        let dy = self.p2[1] - self.p1[1];
        let l = self.length().max(1e-15);
        [dy / l, -dx / l]
    }
    /// Panel tangent unit vector.
    pub fn tangent(&self) -> [f64; 2] {
        let dx = self.p2[0] - self.p1[0];
        let dy = self.p2[1] - self.p1[1];
        let l = self.length().max(1e-15);
        [dx / l, dy / l]
    }
    /// Source panel influence coefficient: normal velocity at point p due to this panel.
    /// Uses 2-D logarithmic source kernel.
    pub fn source_influence(&self, p: [f64; 2]) -> f64 {
        let xc = self.center()[0];
        let yc = self.center()[1];
        let dx = p[0] - xc;
        let dy = p[1] - yc;
        let r2 = dx * dx + dy * dy;
        if r2 < 1e-20 {
            return 0.5;
        }
        let n = self.normal();
        let l = self.length();
        let r = r2.sqrt();
        l / (2.0 * PI) * (n[0] * dx + n[1] * dy) / r2 * l.min(r)
    }
    /// Vortex panel influence (velocity potential induced at field point).
    pub fn vortex_influence(&self, p: [f64; 2]) -> f64 {
        let xc = self.center()[0];
        let yc = self.center()[1];
        let dx = p[0] - xc;
        let dy = p[1] - yc;
        let r2 = (dx * dx + dy * dy).max(1e-20);
        let t = self.tangent();
        let l = self.length();
        -l / (2.0 * PI) * (t[0] * dy - t[1] * dx) / r2
    }
}
/// Hydrodynamic body with added mass, radiation damping, and Froude-Krylov forces.
#[derive(Debug, Clone)]
pub struct HydrodynamicBody {
    /// Added mass tensor (6×6 for 6-DoF) — stored as diagonal for simplicity.
    pub added_mass: [f64; 6],
    /// Radiation damping diagonal \[N·s/m\].
    pub radiation_damping: [f64; 6],
    /// Water density \[kg/m³\].
    pub rho_water: f64,
    /// Body volume \[m³\].
    pub volume: f64,
    /// Froude-Krylov force coefficient (fraction of pressure force).
    pub fk_coeff: f64,
}
impl HydrodynamicBody {
    /// Create a hydrodynamic body.
    pub fn new(
        added_mass: [f64; 6],
        radiation_damping: [f64; 6],
        rho_water: f64,
        volume: f64,
        fk_coeff: f64,
    ) -> Self {
        Self {
            added_mass,
            radiation_damping,
            rho_water,
            volume,
            fk_coeff,
        }
    }
    /// Create a sphere of given radius with theoretical added mass.
    pub fn sphere(radius: f64, rho_water: f64) -> Self {
        let volume = 4.0 / 3.0 * PI * radius.powi(3);
        let m_a = 0.5 * rho_water * volume;
        let added_mass = [m_a, m_a, m_a, 0.0, 0.0, 0.0];
        let damping = [0.0; 6];
        Self::new(added_mass, damping, rho_water, volume, 0.97)
    }
    /// Buoyancy force \[N\] (vertical, upward positive).
    pub fn buoyancy(&self) -> f64 {
        const G: f64 = 9.81;
        self.rho_water * G * self.volume
    }
    /// Added mass force in direction i for acceleration a_i.
    pub fn added_mass_force(&self, accel: [f64; 6]) -> [f64; 6] {
        let mut f = [0.0; 6];
        for i in 0..6 {
            f[i] = -self.added_mass[i] * accel[i];
        }
        f
    }
    /// Radiation damping force for velocity v.
    pub fn radiation_force(&self, vel: [f64; 6]) -> [f64; 6] {
        let mut f = [0.0; 6];
        for i in 0..6 {
            f[i] = -self.radiation_damping[i] * vel[i];
        }
        f
    }
    /// Froude-Krylov force (diffraction excluded, simplified).
    /// `wave_pressure_grad` is the incident wave pressure gradient \[Pa/m\].
    pub fn froude_krylov_force(&self, wave_pressure_grad: f64) -> f64 {
        self.fk_coeff * wave_pressure_grad * self.volume
    }
}
/// Morison equation for force on a cylinder in waves.
#[derive(Debug, Clone, Copy)]
pub struct MorisonEquation {
    /// Drag coefficient Cd.
    pub cd: f64,
    /// Inertia coefficient Cm.
    pub cm: f64,
    /// Cylinder diameter \[m\].
    pub diameter: f64,
    /// Cylinder length per unit \[m\].
    pub length: f64,
    /// Fluid density \[kg/m³\].
    pub rho: f64,
}
impl MorisonEquation {
    /// Create a Morison equation model.
    pub fn new(cd: f64, cm: f64, diameter: f64, length: f64, rho: f64) -> Self {
        Self {
            cd,
            cm,
            diameter,
            length,
            rho,
        }
    }
    /// Cross-sectional area \[m²\].
    pub fn cross_section(&self) -> f64 {
        PI / 4.0 * self.diameter * self.diameter
    }
    /// Total inline force per unit length \[N/m\] from Morison equation.
    /// `u_w` = water particle velocity \[m/s\], `a_w` = acceleration \[m/s²\].
    pub fn force_per_length(&self, u_w: f64, a_w: f64) -> f64 {
        let f_drag = 0.5 * self.rho * self.cd * self.diameter * u_w * u_w.abs();
        let f_inertia = self.rho * self.cm * self.cross_section() * a_w;
        f_drag + f_inertia
    }
    /// Total force over cylinder length.
    pub fn total_force(&self, u_w: f64, a_w: f64) -> f64 {
        self.force_per_length(u_w, a_w) * self.length
    }
    /// Keulegan-Carpenter number KC = U * T / D.
    /// `u_max` = max velocity amplitude \[m/s\], `period` = wave period \[s\].
    pub fn kc_number(&self, u_max: f64, period: f64) -> f64 {
        u_max * period / self.diameter
    }
    /// Check if flow is drag-dominated (KC > 40).
    pub fn is_drag_dominated(&self, u_max: f64, period: f64) -> bool {
        self.kc_number(u_max, period) > 40.0
    }
}
/// Extended Morison analysis including drag force, inertia force separation.
#[derive(Debug, Clone, Copy)]
pub struct ExtendedMorison {
    /// Drag coefficient Cd.
    pub cd: f64,
    /// Added mass coefficient Ca = Cm - 1.
    pub ca: f64,
    /// Cylinder diameter \[m\].
    pub diameter: f64,
    /// Fluid density \[kg/m³\].
    pub rho: f64,
}
impl ExtendedMorison {
    /// Create an extended Morison model.
    pub fn new(cd: f64, ca: f64, diameter: f64, rho: f64) -> Self {
        Self {
            cd,
            ca,
            diameter,
            rho,
        }
    }
    /// Cross-sectional area \[m²\].
    pub fn area(&self) -> f64 {
        PI / 4.0 * self.diameter * self.diameter
    }
    /// Inertia coefficient Cm = 1 + Ca.
    pub fn cm(&self) -> f64 {
        1.0 + self.ca
    }
    /// Drag force per unit length \[N/m\].
    pub fn drag_force(&self, u: f64) -> f64 {
        0.5 * self.rho * self.cd * self.diameter * u * u.abs()
    }
    /// Inertia (Froude-Krylov + added mass) force per unit length \[N/m\].
    pub fn inertia_force(&self, a_fluid: f64) -> f64 {
        self.rho * self.cm() * self.area() * a_fluid
    }
    /// Added mass force per unit length \[N/m\] (structure acceleration a_s).
    pub fn added_mass_force(&self, a_structure: f64) -> f64 {
        -self.rho * self.ca * self.area() * a_structure
    }
    /// Total inline force per unit length \[N/m\].
    pub fn total_force(&self, u_fluid: f64, a_fluid: f64, a_structure: f64) -> f64 {
        let u_rel = u_fluid;
        self.drag_force(u_rel) + self.inertia_force(a_fluid) + self.added_mass_force(a_structure)
    }
    /// Maximum force per unit length for regular wave (Airy theory).
    /// `h_wave` = wave height \[m\], `T` = period \[s\], `depth` = water depth \[m\].
    pub fn max_force_regular_wave(&self, h_wave: f64, period: f64, depth: f64) -> f64 {
        let omega = 2.0 * PI / period;
        let k = omega * omega / 9.81;
        let _kd = k * depth;
        let u_max = PI * h_wave / period;
        let a_max = omega * u_max;
        let f_drag = self.drag_force(u_max);
        let f_inertia = self.inertia_force(a_max);
        f_drag.max(f_inertia)
    }
}
/// Ship hydrodynamics: resistance, wave-making, and seakeeping.
#[derive(Debug, Clone)]
pub struct ShipHydrodynamics {
    /// Ship length between perpendiculars L \[m\].
    pub length: f64,
    /// Ship beam B \[m\].
    pub beam: f64,
    /// Ship draft T \[m\].
    pub draft: f64,
    /// Block coefficient Cb (≈ 0.5–0.85).
    pub block_coeff: f64,
    /// Water density ρ \[kg/m³\].
    pub rho: f64,
    /// Kinematic viscosity ν \[m²/s\].
    pub nu: f64,
}
impl ShipHydrodynamics {
    /// Create a ship hydrodynamics model.
    pub fn new(length: f64, beam: f64, draft: f64, block_coeff: f64, rho: f64, nu: f64) -> Self {
        Self {
            length,
            beam,
            draft,
            block_coeff,
            rho,
            nu,
        }
    }
    /// Displacement volume ∇ = Cb * L * B * T \[m³\].
    pub fn displacement_volume(&self) -> f64 {
        self.block_coeff * self.length * self.beam * self.draft
    }
    /// Displacement mass Δ = ρ ∇ \[kg\].
    pub fn displacement_mass(&self) -> f64 {
        self.rho * self.displacement_volume()
    }
    /// Froude number Fr = V / √(g L).
    pub fn froude_number(&self, speed: f64) -> f64 {
        let g = 9.81;
        speed / (g * self.length).sqrt()
    }
    /// Reynolds number Re = V L / ν.
    pub fn reynolds_number(&self, speed: f64) -> f64 {
        speed * self.length / self.nu
    }
    /// ITTC 1957 friction coefficient Cf = 0.075 / (log10(Re) - 2)².
    pub fn friction_coefficient(&self, speed: f64) -> f64 {
        let re = self.reynolds_number(speed);
        if re < 1e5 {
            return 0.0;
        }
        let log_re = re.log10();
        0.075 / (log_re - 2.0).powi(2)
    }
    /// Wetted surface area S (Denny approximation).
    pub fn wetted_surface(&self) -> f64 {
        self.length * (1.7 * self.draft + self.block_coeff * self.beam)
    }
    /// Frictional resistance Rf = 0.5 ρ V² S Cf \[N\].
    pub fn frictional_resistance(&self, speed: f64) -> f64 {
        let cf = self.friction_coefficient(speed);
        let s = self.wetted_surface();
        0.5 * self.rho * speed * speed * s * cf
    }
    /// Wave-making resistance coefficient Cw (Michell integral approximation).
    /// Uses simple polynomial fit for Froude number range 0.1–0.5.
    pub fn wave_resistance_coefficient(&self, speed: f64) -> f64 {
        let fr = self.froude_number(speed);
        if !(0.05..=0.6).contains(&fr) {
            return 0.0;
        }
        let a = self.block_coeff;
        let fr_peak = 0.28 + 0.05 * a;
        let amplitude = 3e-3 * a;
        amplitude * (-(((fr - fr_peak) / 0.1).powi(2))).exp()
    }
    /// Wave-making resistance Rw \[N\].
    pub fn wave_resistance(&self, speed: f64) -> f64 {
        let cw = self.wave_resistance_coefficient(speed);
        let s = self.wetted_surface();
        0.5 * self.rho * speed * speed * s * cw
    }
    /// Total resistance Rt = Rf + Rw (simplified) \[N\].
    pub fn total_resistance(&self, speed: f64) -> f64 {
        self.frictional_resistance(speed) + self.wave_resistance(speed)
    }
    /// Effective power Pe = Rt * V \[W\].
    pub fn effective_power(&self, speed: f64) -> f64 {
        self.total_resistance(speed) * speed
    }
    /// Added mass for surge (longitudinal) motion: ma ≈ 0.07 Δ.
    pub fn added_mass_surge(&self) -> f64 {
        0.07 * self.displacement_mass()
    }
    /// Added mass for sway (transverse) motion: ma ≈ 0.8 Δ for full forms.
    pub fn added_mass_sway(&self) -> f64 {
        0.8 * self.displacement_mass()
    }
    /// Added mass for heave motion: ma ≈ ρ π/4 B² L Cm.
    pub fn added_mass_heave(&self) -> f64 {
        self.rho * PI / 4.0 * self.beam * self.beam * self.length * 1.0
    }
}
/// Cavitation index and inception analysis.
#[derive(Debug, Clone, Copy)]
pub struct CavitationAnalysis {
    /// Free-stream velocity \[m/s\].
    pub v_inf: f64,
    /// Free-stream pressure \[Pa\].
    pub p_inf: f64,
    /// Vapor pressure of liquid \[Pa\] (at operating temperature).
    pub p_vapor: f64,
    /// Liquid density \[kg/m³\].
    pub rho: f64,
}
impl CavitationAnalysis {
    /// Create a cavitation analysis model.
    pub fn new(v_inf: f64, p_inf: f64, p_vapor: f64, rho: f64) -> Self {
        Self {
            v_inf,
            p_inf,
            p_vapor,
            rho,
        }
    }
    /// Thoma (cavitation) number σ = (p_inf - p_v) / (0.5 ρ V²).
    pub fn thoma_number(&self) -> f64 {
        let dyn_pressure = 0.5 * self.rho * self.v_inf * self.v_inf;
        if dyn_pressure < 1e-14 {
            return f64::INFINITY;
        }
        (self.p_inf - self.p_vapor) / dyn_pressure
    }
    /// Local cavitation index σ_loc at a point with pressure p_loc.
    pub fn local_cavitation_index(&self, p_loc: f64) -> f64 {
        let dyn_pressure = 0.5 * self.rho * self.v_inf * self.v_inf;
        if dyn_pressure < 1e-14 {
            return f64::INFINITY;
        }
        (p_loc - self.p_vapor) / dyn_pressure
    }
    /// Inception cavitation: returns true if cavitation starts (σ < σ_incipient).
    /// σ_incipient ≈ |Cp_min| (minimum pressure coefficient).
    pub fn inception_cavitation(&self, cp_min: f64) -> bool {
        self.thoma_number() < cp_min.abs()
    }
    /// Critical velocity for cavitation onset: V_c = √(2(p_inf - p_v) / (ρ |Cp_min|)).
    pub fn critical_velocity(&self, cp_min: f64) -> f64 {
        let denom = self.rho * cp_min.abs();
        if denom < 1e-14 {
            return f64::INFINITY;
        }
        (2.0 * (self.p_inf - self.p_vapor) / denom).max(0.0).sqrt()
    }
    /// Bubble radius growth rate (Rayleigh-Plesset simplified).
    /// ṘR̈ = (p_b - p_inf) / ρ.
    pub fn bubble_growth_rate(&self, r: f64, p_bubble: f64) -> f64 {
        if r < 1e-15 {
            return 0.0;
        }
        ((p_bubble - self.p_inf) / self.rho).max(0.0).sqrt()
    }
    /// Collapse time (Rayleigh collapse): t_c = 0.915 R₀ √(ρ / (p_inf - p_v)).
    pub fn collapse_time(&self, r0: f64) -> f64 {
        let dp = self.p_inf - self.p_vapor;
        if dp < 1e-14 {
            return f64::INFINITY;
        }
        0.915 * r0 * (self.rho / dp).sqrt()
    }
}
/// Wake / slipstream model for velocity deficit behind a body.
#[derive(Debug, Clone, Copy)]
pub struct WakeEffect {
    /// Thrust/drag coefficient Ct.
    pub ct: f64,
    /// Rotor/body diameter \[m\].
    pub diameter: f64,
    /// Wake expansion factor k.
    pub expansion_factor: f64,
}
impl WakeEffect {
    /// Create a wake model.
    pub fn new(ct: f64, diameter: f64, expansion_factor: f64) -> Self {
        Self {
            ct,
            diameter,
            expansion_factor,
        }
    }
    /// Wake velocity deficit factor at downstream distance x \[m\] (Jensen model).
    /// Returns u_wake / u_inf.
    pub fn velocity_ratio(&self, x: f64, freestream_v: f64) -> f64 {
        let _ = freestream_v;
        let r = self.diameter / 2.0;
        let r_wake = r + self.expansion_factor * x;
        let area_ratio = (r / r_wake).powi(2);
        1.0 - (1.0 - (1.0 - self.ct).sqrt()) * area_ratio
    }
    /// Wake velocity at distance x \[m\].
    pub fn wake_velocity(&self, x: f64, freestream_v: f64) -> f64 {
        self.velocity_ratio(x, freestream_v) * freestream_v
    }
    /// Blockage ratio for array of turbines: fraction of cross-section blocked.
    pub fn blockage_ratio(&self, array_spacing: f64) -> f64 {
        let a_body = PI / 4.0 * self.diameter * self.diameter;
        let a_array = array_spacing * array_spacing;
        a_body / a_array
    }
}
/// Vortex Lattice Method solver for steady aerodynamics.
#[derive(Debug, Clone)]
pub struct VortexLatticeSolver {
    /// Wing panels (horseshoe vortices).
    pub panels: Vec<HorseshoeVortex>,
    /// Free-stream velocity vector \[m/s\].
    pub v_inf: [f64; 3],
    /// Air density \[kg/m³\].
    pub rho: f64,
}
impl VortexLatticeSolver {
    /// Create a VLM solver.
    pub fn new(panels: Vec<HorseshoeVortex>, v_inf: [f64; 3], rho: f64) -> Self {
        Self { panels, v_inf, rho }
    }
    /// Compute AIC (Aerodynamic Influence Coefficient) matrix entry:
    /// normal velocity at control point i due to unit-strength horseshoe j.
    pub fn aic_entry(&self, i: usize, j: usize) -> f64 {
        if i >= self.panels.len() || j >= self.panels.len() {
            return 0.0;
        }
        let mut panel_j = self.panels[j].clone();
        panel_j.gamma = 1.0;
        let cp = self.panels[i].control_point;
        let v = panel_j.total_velocity(cp);
        let n = self.panels[i].normal;
        dot3(v, n)
    }
    /// RHS vector: -V_inf · n_i for each panel.
    pub fn rhs(&self) -> Vec<f64> {
        self.panels
            .iter()
            .map(|p| -dot3(self.v_inf, p.normal))
            .collect()
    }
    /// Total lift \[N\] from Kutta-Joukowski: L = ρ V_inf × Γ * span.
    pub fn total_lift(&self) -> f64 {
        let v_mag = norm3(self.v_inf);
        self.panels
            .iter()
            .map(|p| {
                let dl = [p.b[0] - p.a[0], p.b[1] - p.a[1], p.b[2] - p.a[2]];
                let span_len = norm3(dl);
                self.rho * v_mag * p.gamma * span_len
            })
            .sum()
    }
}
/// Water entry slam force model (von Karman / Wagner theory).
#[derive(Debug, Clone, Copy)]
pub struct SlamForce {
    /// Deadrise angle β \[rad\].
    pub deadrise_angle: f64,
    /// Water density \[kg/m³\].
    pub rho_water: f64,
    /// Entry width (beam) at waterline \[m\].
    pub beam: f64,
}
impl SlamForce {
    /// Create a slam force model.
    pub fn new(deadrise_angle: f64, rho_water: f64, beam: f64) -> Self {
        Self {
            deadrise_angle,
            rho_water,
            beam,
        }
    }
    /// Von Karman added mass coefficient C_K for wedge entry.
    pub fn added_mass_coeff(&self) -> f64 {
        PI / (2.0 * self.deadrise_angle.tan() * self.deadrise_angle.tan())
    }
    /// Wagner impact pressure coefficient (peak).
    pub fn wagner_pressure_coeff(&self) -> f64 {
        PI * PI / (4.0 * (self.deadrise_angle.tan()).powi(2))
    }
    /// Maximum impact pressure \[Pa\] for entry velocity v \[m/s\].
    pub fn max_impact_pressure(&self, v: f64) -> f64 {
        0.5 * self.rho_water * v * v * self.wagner_pressure_coeff()
    }
    /// Impact force \[N\] for entry velocity v \[m/s\] and submergence depth s \[m\].
    pub fn impact_force(&self, v: f64, s: f64) -> f64 {
        let ca = self.added_mass_coeff();
        self.rho_water * ca * v * v * self.beam * s
    }
    /// Entry velocity for given momentum impulse and mass.
    pub fn entry_velocity_from_impulse(&self, impulse: f64, mass: f64) -> f64 {
        impulse / mass
    }
}
/// Source panel method for 2-D potential flow around a body.
#[derive(Debug, Clone)]
pub struct SourcePanelMethod {
    /// Panels making up the body.
    pub panels: Vec<Panel>,
    /// Free-stream velocity magnitude \[m/s\].
    pub u_inf: f64,
    /// Free-stream angle of attack \[rad\].
    pub alpha: f64,
}
impl SourcePanelMethod {
    /// Create a source panel method model.
    pub fn new(panels: Vec<Panel>, u_inf: f64, alpha: f64) -> Self {
        Self {
            panels,
            u_inf,
            alpha,
        }
    }
    /// Create a circle approximated by n panels.
    pub fn circle(radius: f64, n_panels: usize, u_inf: f64, alpha: f64) -> Self {
        let mut panels = Vec::with_capacity(n_panels);
        for i in 0..n_panels {
            let theta1 = 2.0 * PI * i as f64 / n_panels as f64;
            let theta2 = 2.0 * PI * (i + 1) as f64 / n_panels as f64;
            let p1 = [radius * theta1.cos(), radius * theta1.sin()];
            let p2 = [radius * theta2.cos(), radius * theta2.sin()];
            panels.push(Panel::new(p1, p2));
        }
        Self::new(panels, u_inf, alpha)
    }
    /// Assemble the influence coefficient matrix A (n×n).
    pub fn influence_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.panels.len();
        let mut a = vec![vec![0.0f64; n]; n];
        for (i, a_row) in a.iter_mut().enumerate() {
            let cp = self.panels[i].center();
            let ni = self.panels[i].normal();
            for (j, a_ij) in a_row.iter_mut().enumerate() {
                if i == j {
                    *a_ij = 0.5;
                } else {
                    let inf = self.panels[j].source_influence(cp);
                    *a_ij = inf * (ni[0] * 1.0 + ni[1] * 1.0).signum().abs();
                    *a_ij = self.panels[j].source_influence(cp);
                }
            }
        }
        a
    }
    /// Free-stream normal velocity at each panel.
    pub fn freestream_rhs(&self) -> Vec<f64> {
        let ca = self.alpha.cos();
        let sa = self.alpha.sin();
        self.panels
            .iter()
            .map(|p| {
                let n = p.normal();
                -(self.u_inf * ca * n[0] + self.u_inf * sa * n[1])
            })
            .collect()
    }
    /// Pressure coefficient from velocity at panel center.
    /// Cp = 1 - (v/u_inf)².
    pub fn cp_at_panel(&self, v_tangential: f64) -> f64 {
        if self.u_inf.abs() < 1e-14 {
            return 0.0;
        }
        1.0 - (v_tangential / self.u_inf).powi(2)
    }
}
/// Drag force models for different flow regimes.
#[derive(Debug, Clone, Copy)]
pub struct DragModel {
    /// Fluid dynamic viscosity \[Pa·s\].
    pub viscosity: f64,
    /// Fluid density \[kg/m³\].
    pub rho: f64,
    /// Particle/body diameter \[m\].
    pub diameter: f64,
}
impl DragModel {
    /// Create a drag model.
    pub fn new(viscosity: f64, rho: f64, diameter: f64) -> Self {
        Self {
            viscosity,
            rho,
            diameter,
        }
    }
    /// Reynolds number for velocity v \[m/s\].
    pub fn reynolds(&self, v: f64) -> f64 {
        self.rho * v * self.diameter / self.viscosity
    }
    /// Stokes drag (Re << 1): F = 3π μ d v.
    pub fn stokes_drag(&self, v: f64) -> f64 {
        3.0 * PI * self.viscosity * self.diameter * v
    }
    /// Newton (turbulent) drag: F = 0.5 * Cd * ρ * A * v².
    /// Uses Cd = 0.44 for a sphere.
    pub fn newton_drag(&self, v: f64) -> f64 {
        let cd = 0.44;
        let a = PI / 4.0 * self.diameter * self.diameter;
        0.5 * cd * self.rho * a * v * v
    }
    /// Schiller-Naumann drag coefficient for intermediate Re.
    pub fn schiller_naumann_cd(&self, v: f64) -> f64 {
        let re = self.reynolds(v).max(1e-10);
        if re < 1000.0 {
            24.0 / re * (1.0 + 0.15 * re.powf(0.687))
        } else {
            0.44
        }
    }
    /// Schiller-Naumann drag force.
    pub fn intermediate_drag(&self, v: f64) -> f64 {
        let cd = self.schiller_naumann_cd(v);
        let a = PI / 4.0 * self.diameter * self.diameter;
        0.5 * cd * self.rho * a * v * v
    }
    /// Select appropriate drag model based on Reynolds number.
    pub fn drag_force(&self, v: f64) -> f64 {
        let re = self.reynolds(v.abs());
        if re < 1.0 {
            self.stokes_drag(v.abs())
        } else {
            self.intermediate_drag(v.abs())
        }
    }
}
/// Rayleigh-Taylor instability: heavy fluid above light fluid.
#[derive(Debug, Clone, Copy)]
pub struct RayleighTaylor {
    /// Density of upper (heavy) fluid ρ_H \[kg/m³\].
    pub rho_heavy: f64,
    /// Density of lower (light) fluid ρ_L \[kg/m³\].
    pub rho_light: f64,
    /// Surface tension σ \[N/m\].
    pub surface_tension: f64,
    /// Gravitational acceleration g \[m/s²\].
    pub gravity: f64,
}
impl RayleighTaylor {
    /// Create a Rayleigh-Taylor analysis.
    pub fn new(rho_heavy: f64, rho_light: f64, surface_tension: f64, gravity: f64) -> Self {
        Self {
            rho_heavy,
            rho_light,
            surface_tension,
            gravity,
        }
    }
    /// Atwood number A = (ρ_H - ρ_L) / (ρ_H + ρ_L).
    pub fn atwood_number(&self) -> f64 {
        (self.rho_heavy - self.rho_light) / (self.rho_heavy + self.rho_light)
    }
    /// Growth rate σ(k) = √(A g k - σ k³ / (ρ_H+ρ_L)).
    pub fn growth_rate(&self, k: f64) -> f64 {
        let rho_sum = self.rho_heavy + self.rho_light;
        let unstable = self.atwood_number() * self.gravity * k;
        let stable = self.surface_tension * k * k * k / rho_sum;
        let sigma_sq = unstable - stable;
        if sigma_sq <= 0.0 {
            0.0
        } else {
            sigma_sq.sqrt()
        }
    }
    /// Critical wavenumber k_c: above k_c, interface is stable.
    pub fn critical_wavenumber(&self) -> f64 {
        let rho_sum = self.rho_heavy + self.rho_light;
        if self.surface_tension < 1e-14 {
            return f64::INFINITY;
        }
        (self.atwood_number() * self.gravity * rho_sum / self.surface_tension)
            .max(0.0)
            .sqrt()
    }
    /// Most unstable wavenumber k_max = k_c / √3.
    pub fn most_unstable_wavenumber(&self) -> f64 {
        self.critical_wavenumber() / 3.0_f64.sqrt()
    }
    /// Maximum growth rate at most unstable mode.
    pub fn max_growth_rate(&self) -> f64 {
        self.growth_rate(self.most_unstable_wavenumber())
    }
    /// Bubble rise velocity (classical RT: v_b ≈ 0.23 √(A g λ_c)).
    pub fn bubble_velocity(&self) -> f64 {
        let kc = self.critical_wavenumber();
        if kc < 1e-14 {
            return 0.0;
        }
        let lambda_c = 2.0 * PI / kc;
        0.23 * (self.atwood_number() * self.gravity * lambda_c)
            .max(0.0)
            .sqrt()
    }
}
/// Aerodynamic torque model: pitching moment and damping.
#[derive(Debug, Clone, Copy)]
pub struct AerodynamicTorque {
    /// Pitching moment coefficient slope dCm/dα.
    pub cm_alpha: f64,
    /// Aerodynamic damping coefficient Cmq.
    pub cm_q: f64,
    /// Reference area \[m²\].
    pub ref_area: f64,
    /// Reference chord \[m\].
    pub chord: f64,
    /// Air density \[kg/m³\].
    pub rho: f64,
}
impl AerodynamicTorque {
    /// Create an aerodynamic torque model.
    pub fn new(cm_alpha: f64, cm_q: f64, ref_area: f64, chord: f64, rho: f64) -> Self {
        Self {
            cm_alpha,
            cm_q,
            ref_area,
            chord,
            rho,
        }
    }
    /// Pitching moment \[N·m\] for angle of attack α \[rad\] and velocity v \[m/s\].
    pub fn pitching_moment(&self, alpha: f64, v: f64) -> f64 {
        let q = 0.5 * self.rho * v * v;
        self.cm_alpha * alpha * q * self.ref_area * self.chord
    }
    /// Aerodynamic damping torque for pitch rate q \[rad/s\] and velocity v \[m/s\].
    pub fn damping_torque(&self, pitch_rate: f64, v: f64) -> f64 {
        if v < 1e-10 {
            return 0.0;
        }
        let q = 0.5 * self.rho * v * v;
        let cm_q_dim = pitch_rate * self.chord / (2.0 * v);
        self.cm_q * cm_q_dim * q * self.ref_area * self.chord
    }
    /// Total aerodynamic moment \[N·m\].
    pub fn total_moment(&self, alpha: f64, pitch_rate: f64, v: f64) -> f64 {
        self.pitching_moment(alpha, v) + self.damping_torque(pitch_rate, v)
    }
}
/// Aerodynamic body with lift/drag coefficient look-up tables vs angle of attack.
#[derive(Debug, Clone)]
pub struct AerodynamicBody {
    /// Reference area \[m²\].
    pub ref_area: f64,
    /// Reference chord \[m\].
    pub chord: f64,
    /// Angle of attack table \[rad\].
    pub aoa_table: Vec<f64>,
    /// Lift coefficient table (Cl vs alpha).
    pub cl_table: Vec<f64>,
    /// Drag coefficient table (Cd vs alpha).
    pub cd_table: Vec<f64>,
    /// Air density \[kg/m³\].
    pub rho: f64,
}
impl AerodynamicBody {
    /// Create an aerodynamic body with look-up tables.
    pub fn new(
        ref_area: f64,
        chord: f64,
        aoa_table: Vec<f64>,
        cl_table: Vec<f64>,
        cd_table: Vec<f64>,
        rho: f64,
    ) -> Self {
        Self {
            ref_area,
            chord,
            aoa_table,
            cl_table,
            cd_table,
            rho,
        }
    }
    /// Create with default flat plate coefficients.
    pub fn flat_plate(ref_area: f64, chord: f64, rho: f64) -> Self {
        let n = 19;
        let aoas: Vec<f64> = (0..n)
            .map(|i| (-45.0 + i as f64 * 5.0).to_radians())
            .collect();
        let cls: Vec<f64> = aoas.iter().map(|&a| 2.0 * PI * a).collect();
        let cds: Vec<f64> = aoas
            .iter()
            .map(|&a| 0.01 + 2.0 * (a.sin()).powi(2))
            .collect();
        Self::new(ref_area, chord, aoas, cls, cds, rho)
    }
    /// Interpolate coefficient from table at given alpha \[rad\].
    fn interp_coeff(table_x: &[f64], table_y: &[f64], x: f64) -> f64 {
        let n = table_x.len();
        if n == 0 {
            return 0.0;
        }
        if x <= table_x[0] {
            return table_y[0];
        }
        if x >= table_x[n - 1] {
            return table_y[n - 1];
        }
        for i in 1..n {
            if x <= table_x[i] {
                let t = (x - table_x[i - 1]) / (table_x[i] - table_x[i - 1]);
                return table_y[i - 1] + t * (table_y[i] - table_y[i - 1]);
            }
        }
        table_y[n - 1]
    }
    /// Lift coefficient at angle of attack α \[rad\].
    pub fn cl(&self, alpha: f64) -> f64 {
        Self::interp_coeff(&self.aoa_table, &self.cl_table, alpha)
    }
    /// Drag coefficient at angle of attack α \[rad\].
    pub fn cd(&self, alpha: f64) -> f64 {
        Self::interp_coeff(&self.aoa_table, &self.cd_table, alpha)
    }
    /// Aerodynamic force \[N\] at velocity v \[m/s\] and angle of attack α \[rad\].
    /// Returns (lift, drag) magnitudes.
    pub fn forces(&self, v: f64, alpha: f64) -> (f64, f64) {
        let q = 0.5 * self.rho * v * v;
        let lift = self.cl(alpha) * q * self.ref_area;
        let drag = self.cd(alpha) * q * self.ref_area;
        (lift, drag)
    }
}
/// Ground effect: image-vortex lift augmentation near ground.
#[derive(Debug, Clone, Copy)]
pub struct GroundEffect {
    /// Wingspan b \[m\].
    pub wingspan: f64,
    /// Reference lift without ground effect (in free air).
    pub lift_free: f64,
    /// Aspect ratio AR.
    pub aspect_ratio: f64,
}
impl GroundEffect {
    /// Create a ground effect model.
    pub fn new(wingspan: f64, lift_free: f64, aspect_ratio: f64) -> Self {
        Self {
            wingspan,
            lift_free,
            aspect_ratio,
        }
    }
    /// Ground effect ratio as function of h/b (Cheeseman & Bennett approximation).
    /// Returns factor by which lift is multiplied (>1 near ground).
    pub fn lift_ratio(&self, height: f64) -> f64 {
        let hb = height / self.wingspan;
        if hb > 2.0 {
            return 1.0;
        }
        let k = 1.0 - 1.32 * (-(1.05 * hb * 10.0).powf(0.8)).exp();
        1.0 + (1.0 - k) / k.max(0.01)
    }
    /// Lift with ground effect at height h \[m\].
    pub fn lift_with_ground_effect(&self, height: f64) -> f64 {
        self.lift_free * self.lift_ratio(height)
    }
    /// Induced drag reduction due to ground effect.
    pub fn induced_drag_ratio(&self, height: f64) -> f64 {
        let hb = height / self.wingspan;
        if hb > 2.0 {
            return 1.0;
        }
        1.0 / self.lift_ratio(height)
    }
}
/// Lift force models including thin airfoil and Kutta-Joukowski.
#[derive(Debug, Clone, Copy)]
pub struct LiftModel {
    /// Air density \[kg/m³\].
    pub rho: f64,
    /// Flow velocity \[m/s\].
    pub velocity: f64,
    /// Reference span \[m\].
    pub span: f64,
    /// Reference chord \[m\].
    pub chord: f64,
}
impl LiftModel {
    /// Create a lift model.
    pub fn new(rho: f64, velocity: f64, span: f64, chord: f64) -> Self {
        Self {
            rho,
            velocity,
            span,
            chord,
        }
    }
    /// Thin airfoil theory: Cl = 2π α.
    pub fn thin_airfoil_cl(&self, alpha: f64) -> f64 {
        2.0 * PI * alpha
    }
    /// Lift force from thin airfoil theory \[N\].
    pub fn thin_airfoil_lift(&self, alpha: f64) -> f64 {
        let cl = self.thin_airfoil_cl(alpha);
        let q = 0.5 * self.rho * self.velocity * self.velocity;
        cl * q * self.span * self.chord
    }
    /// Kutta-Joukowski lift: L = ρ V Γ b.
    /// `circulation` = circulation \[m²/s\].
    pub fn kutta_joukowski_lift(&self, circulation: f64) -> f64 {
        self.rho * self.velocity * circulation * self.span
    }
    /// NACA 4-digit maximum camber at x/c location.
    /// `m` = max camber fraction, `p` = max camber position fraction.
    pub fn naca_camber(&self, x_over_c: f64, m: f64, p: f64) -> f64 {
        if x_over_c <= p {
            m / (p * p) * (2.0 * p * x_over_c - x_over_c * x_over_c)
        } else {
            m / ((1.0 - p) * (1.0 - p)) * (1.0 - 2.0 * p + 2.0 * p * x_over_c - x_over_c * x_over_c)
        }
    }
    /// NACA 4-digit lift coefficient using thin-airfoil ideal: Cl ≈ 2π(α + α_ZL).
    pub fn naca_lift_coefficient(&self, alpha: f64, m: f64, p: f64) -> f64 {
        let alpha_zl = -2.0 * PI * m * (1.0 - p);
        2.0 * PI * (alpha - alpha_zl)
    }
}
/// Doublet panel method for 3-D potential flow (simplified formulation).
#[derive(Debug, Clone)]
pub struct DoubletPanelMethod {
    /// Doublet strengths μ at each panel.
    pub mu: Vec<f64>,
    /// Panel areas \[m²\].
    pub areas: Vec<f64>,
    /// Number of panels.
    pub n_panels: usize,
}
impl DoubletPanelMethod {
    /// Create a doublet panel model.
    pub fn new(n_panels: usize, areas: Vec<f64>) -> Self {
        Self {
            mu: vec![0.0; n_panels],
            areas,
            n_panels,
        }
    }
    /// Set doublet strength for all panels.
    pub fn set_strengths(&mut self, mu: Vec<f64>) {
        self.mu = mu;
    }
    /// Compute velocity potential at a far-field point (simplified dipole sum).
    pub fn far_field_potential(&self, r: f64) -> f64 {
        if r < 1e-12 {
            return 0.0;
        }
        let total_mu: f64 = self
            .mu
            .iter()
            .zip(self.areas.iter())
            .map(|(m, a)| m * a)
            .sum();
        total_mu / (4.0 * PI * r * r)
    }
}
/// Prandtl lifting line theory for finite-span wings.
#[derive(Debug, Clone)]
pub struct LiftingLine {
    /// Wingspan b \[m\].
    pub span: f64,
    /// Root chord cr \[m\].
    pub root_chord: f64,
    /// Tip chord ct \[m\] (0 for triangular).
    pub tip_chord: f64,
    /// Root angle of attack \[rad\].
    pub alpha_root: f64,
    /// Air density \[kg/m³\].
    pub rho: f64,
    /// Free-stream velocity \[m/s\].
    pub v_inf: f64,
    /// Number of Fourier modes used.
    pub n_modes: usize,
}
impl LiftingLine {
    /// Create a lifting line model.
    pub fn new(
        span: f64,
        root_chord: f64,
        tip_chord: f64,
        alpha_root: f64,
        rho: f64,
        v_inf: f64,
        n_modes: usize,
    ) -> Self {
        Self {
            span,
            root_chord,
            tip_chord,
            alpha_root,
            rho,
            v_inf,
            n_modes,
        }
    }
    /// Chord at spanwise position y (linear taper).
    pub fn chord_at(&self, y: f64) -> f64 {
        let eta = 2.0 * y / self.span;
        let taper = self.tip_chord / self.root_chord;
        self.root_chord * (1.0 - (1.0 - taper) * eta.abs())
    }
    /// Elliptic lift distribution Γ(y) = Γ₀ * √(1 - (2y/b)²).
    pub fn elliptic_circulation(&self, gamma_root: f64, y: f64) -> f64 {
        let eta = 2.0 * y / self.span;
        gamma_root * (1.0 - eta * eta).max(0.0).sqrt()
    }
    /// Aspect ratio AR = b² / S.
    pub fn aspect_ratio(&self) -> f64 {
        let s = self.wing_area();
        if s < 1e-14 {
            return 0.0;
        }
        self.span * self.span / s
    }
    /// Wing area (trapezoidal).
    pub fn wing_area(&self) -> f64 {
        0.5 * (self.root_chord + self.tip_chord) * self.span
    }
    /// Lift curve slope for finite wing: a = a₀ / (1 + a₀/(π*AR)).
    /// `a0` = 2π (thin airfoil, infinite span).
    pub fn lift_curve_slope(&self) -> f64 {
        let ar = self.aspect_ratio();
        let a0 = 2.0 * PI;
        a0 / (1.0 + a0 / (PI * ar))
    }
    /// 3-D lift coefficient.
    pub fn cl_3d(&self) -> f64 {
        self.lift_curve_slope() * self.alpha_root
    }
    /// Total lift \[N\].
    pub fn total_lift(&self) -> f64 {
        let cl = self.cl_3d();
        let s = self.wing_area();
        let q = 0.5 * self.rho * self.v_inf * self.v_inf;
        cl * q * s
    }
    /// Induced drag coefficient CDi = CL² / (π * AR * e), e ≈ 1 for elliptic.
    pub fn induced_drag_coefficient(&self, oswald_e: f64) -> f64 {
        let cl = self.cl_3d();
        let ar = self.aspect_ratio();
        cl * cl / (PI * ar * oswald_e)
    }
    /// Spanwise induced downwash angle at center: αᵢ = CL / (π * AR).
    pub fn induced_angle(&self) -> f64 {
        let cl = self.cl_3d();
        let ar = self.aspect_ratio();
        if ar < 1e-14 {
            return 0.0;
        }
        cl / (PI * ar)
    }
}
/// A horseshoe vortex element for the vortex lattice method.
#[derive(Debug, Clone)]
pub struct HorseshoeVortex {
    /// Bound vortex start point A.
    pub a: [f64; 3],
    /// Bound vortex end point B.
    pub b: [f64; 3],
    /// Vortex strength Γ.
    pub gamma: f64,
    /// Control point (3/4 chord location).
    pub control_point: [f64; 3],
    /// Panel normal vector.
    pub normal: [f64; 3],
}
impl HorseshoeVortex {
    /// Create a horseshoe vortex element.
    pub fn new(
        a: [f64; 3],
        b: [f64; 3],
        gamma: f64,
        control_point: [f64; 3],
        normal: [f64; 3],
    ) -> Self {
        Self {
            a,
            b,
            gamma,
            control_point,
            normal,
        }
    }
    /// Induced velocity at point p from the bound segment only.
    pub fn bound_velocity(&self, p: [f64; 3]) -> [f64; 3] {
        biot_savart_filament(self.a, self.b, p, self.gamma)
    }
    /// Induced velocity from trailing vortex A (semi-infinite from A in x direction).
    pub fn trailing_a_velocity(&self, p: [f64; 3]) -> [f64; 3] {
        let far = [self.a[0] + 1e6, self.a[1], self.a[2]];
        biot_savart_filament(self.a, far, p, -self.gamma)
    }
    /// Induced velocity from trailing vortex B.
    pub fn trailing_b_velocity(&self, p: [f64; 3]) -> [f64; 3] {
        let far = [self.b[0] + 1e6, self.b[1], self.b[2]];
        biot_savart_filament(self.b, far, p, self.gamma)
    }
    /// Total induced velocity at p (bound + trailing).
    pub fn total_velocity(&self, p: [f64; 3]) -> [f64; 3] {
        let vb = self.bound_velocity(p);
        let va = self.trailing_a_velocity(p);
        let vbr = self.trailing_b_velocity(p);
        add3(add3(vb, va), vbr)
    }
}
