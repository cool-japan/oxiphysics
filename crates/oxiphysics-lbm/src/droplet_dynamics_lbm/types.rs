//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::{find, probit, sherwood_number, union_uf};

/// Rayleigh-Plateau instability and secondary-breakup model for LBM droplets.
///
/// Detects breakup for cylindrical jets (Rayleigh-Plateau) and spherical
/// droplets (Weber-number criterion).
#[derive(Debug, Clone)]
pub struct BreakupModel {
    /// Critical Weber number We* for secondary breakup (≈ 12 for spheres).
    pub critical_weber: f64,
    /// Surface tension σ \[lattice units\].
    pub surface_tension: f64,
    /// Carrier-phase density ρ_c \[lattice units\].
    pub carrier_density: f64,
}
impl BreakupModel {
    /// Create a new `BreakupModel`.
    pub fn new(critical_weber: f64, surface_tension: f64, carrier_density: f64) -> Self {
        Self {
            critical_weber,
            surface_tension,
            carrier_density,
        }
    }
    /// Rayleigh-Plateau growth rate σ_RP for a jet of radius R₀ and wavenumber k.
    ///
    /// σ_RP² = (σ / (ρ R₀³)) · (k R₀)² (1 − (k R₀)²)   \[inviscid\]
    pub fn rayleigh_plateau_growth_rate(&self, jet_radius: f64, wavenumber: f64) -> f64 {
        if jet_radius <= 0.0 || self.carrier_density <= 0.0 {
            return 0.0;
        }
        let kr = wavenumber * jet_radius;
        if kr >= 1.0 {
            return 0.0;
        }
        let val = (self.surface_tension / (self.carrier_density * jet_radius.powi(3)))
            * kr
            * kr
            * (1.0 - kr * kr);
        if val > 0.0 { val.sqrt() } else { 0.0 }
    }
    /// Most-unstable wavenumber k* = 1 / (√2 R₀).
    pub fn most_unstable_wavenumber(&self, jet_radius: f64) -> f64 {
        1.0 / (2.0_f64.sqrt() * jet_radius.max(1e-14))
    }
    /// Break-up time estimate: t_breakup ≈ 1 / σ_max.
    pub fn breakup_time(&self, jet_radius: f64) -> f64 {
        let k_star = self.most_unstable_wavenumber(jet_radius);
        let sigma_max = self.rayleigh_plateau_growth_rate(jet_radius, k_star);
        if sigma_max <= 0.0 {
            f64::INFINITY
        } else {
            1.0 / sigma_max
        }
    }
    /// Number of satellite droplets from Rayleigh-Plateau: N ≈ 2π / (k* · 2 R₀) ≈ π.
    pub fn satellite_count_estimate(&self, jet_length: f64, jet_radius: f64) -> f64 {
        let k_star = self.most_unstable_wavenumber(jet_radius);
        let wavelength = 2.0 * PI / k_star.max(1e-14);
        jet_length / wavelength
    }
    /// Daughter droplet radius after Rayleigh-Plateau fragmentation.
    ///
    /// Volume conservation: V_jet = N_drops * (4/3) π r_child³
    /// → r_child = (3 R₀² λ / 4)^(1/3) where λ = 2π/k*.
    pub fn child_radius_plateau(&self, jet_radius: f64) -> f64 {
        let k_star = self.most_unstable_wavenumber(jet_radius);
        let lambda = 2.0 * PI / k_star.max(1e-14);
        (3.0 * jet_radius * jet_radius * lambda / 4.0).cbrt()
    }
    /// Weber-number breakup criterion for a spherical droplet.
    pub fn weber_breakup(&self, droplet_radius: f64, relative_velocity: f64) -> bool {
        let we = self.carrier_density * relative_velocity * relative_velocity * droplet_radius
            / self.surface_tension.max(1e-30);
        we > self.critical_weber
    }
}
/// Size distribution used by the spray model.
#[derive(Debug, Clone, PartialEq)]
pub enum SizeDistribution {
    /// Log-normal distribution: ln(d) ~ N(μ_ln, σ_ln²).
    LogNormal,
    /// Rosin-Rammler distribution: F(d) = 1 - exp(-(d/X̄)^n).
    RosinRammler,
}
/// Young's equation and dynamic contact angle model (Cox-Voinov).
///
/// Implements:
/// - Young's equation: σ_SG - σ_SL = σ_LG cos θ_Y
/// - Cox-Voinov dynamic angle: θ³ ≈ θ_Y³ + 9 Ca ln(L/λ)
/// - Cassie-Baxter for rough surfaces
#[derive(Debug, Clone)]
pub struct ContactAngle {
    /// Young's equilibrium contact angle θ_Y \[rad\].
    pub theta_young: f64,
    /// Liquid-gas surface tension σ_LG \[lattice units\].
    pub sigma_lg: f64,
    /// Solid-gas surface energy σ_SG \[lattice units\].
    pub sigma_sg: f64,
    /// Solid-liquid surface energy σ_SL \[lattice units\].
    pub sigma_sl: f64,
    /// Roughness factor r for Wenzel model (r ≥ 1).
    pub roughness: f64,
    /// Solid area fraction for Cassie-Baxter model (0 < f_s ≤ 1).
    pub solid_fraction: f64,
}
impl ContactAngle {
    /// Create a `ContactAngle` from surface energies (Young's equation).
    ///
    /// Panics-safe: cos θ is clamped to \[-1, 1\].
    pub fn from_energies(sigma_lg: f64, sigma_sg: f64, sigma_sl: f64) -> Self {
        let cos_theta = ((sigma_sg - sigma_sl) / sigma_lg.max(1e-30)).clamp(-1.0, 1.0);
        let theta_young = cos_theta.acos();
        Self {
            theta_young,
            sigma_lg,
            sigma_sg,
            sigma_sl,
            roughness: 1.0,
            solid_fraction: 1.0,
        }
    }
    /// Spreading coefficient S = σ_SG - σ_SL - σ_LG.
    ///
    /// S ≥ 0 → complete wetting; S < 0 → partial wetting.
    pub fn spreading_coefficient(&self) -> f64 {
        self.sigma_sg - self.sigma_sl - self.sigma_lg
    }
    /// Work of adhesion W_A = σ_LG (1 + cos θ_Y).
    pub fn work_of_adhesion(&self) -> f64 {
        self.sigma_lg * (1.0 + self.theta_young.cos())
    }
    /// Wenzel contact angle θ_W: cos θ_W = r cos θ_Y.
    pub fn wenzel_angle(&self) -> f64 {
        let cos_w = (self.roughness * self.theta_young.cos()).clamp(-1.0, 1.0);
        cos_w.acos()
    }
    /// Cassie-Baxter contact angle θ_CB: cos θ_CB = f_s cos θ_Y − (1 − f_s).
    pub fn cassie_baxter_angle(&self) -> f64 {
        let cos_cb = (self.solid_fraction * self.theta_young.cos() - (1.0 - self.solid_fraction))
            .clamp(-1.0, 1.0);
        cos_cb.acos()
    }
    /// Dynamic contact angle θ_dyn from Cox-Voinov model.
    ///
    /// θ_dyn³ ≈ θ_Y³ + 9 Ca · ln(L / λ)
    ///
    /// where Ca = μ U / σ and ln(L/λ) is the macro-to-micro length ratio (~ln 1000).
    pub fn cox_voinov_angle(&self, capillary_number: f64, ln_ratio: f64) -> f64 {
        let theta3 = self.theta_young.powi(3) + 9.0 * capillary_number * ln_ratio;
        theta3.max(0.0).cbrt()
    }
    /// Capillary length l_c = sqrt(σ_LG / (ρ g)).
    pub fn capillary_length(&self, density: f64, gravity: f64) -> f64 {
        (self.sigma_lg / (density * gravity).max(1e-30)).sqrt()
    }
    /// Returns `true` if the surface is hydrophilic (θ_Y < 90°).
    pub fn is_hydrophilic(&self) -> bool {
        self.theta_young < PI * 0.5
    }
    /// Returns `true` if the surface is superhydrophobic (θ_Y > 150°).
    pub fn is_superhydrophobic(&self) -> bool {
        self.theta_young > 150.0_f64.to_radians()
    }
}
/// LBM color-gradient interface tracker for a two-component fluid.
///
/// Stores the red (droplet) and blue (carrier) distribution functions for a
/// 1-D test lattice of `n` nodes (D1Q3 velocity set: {-1, 0, +1}).
#[derive(Debug, Clone)]
pub struct DropletLBM {
    /// Red (droplet-phase) distribution functions, length 3n.
    pub f_red: Vec<f64>,
    /// Blue (carrier-phase) distribution functions, length 3n.
    pub f_blue: Vec<f64>,
    /// Number of lattice nodes.
    pub n: usize,
    /// Relaxation parameter ω for BGK collision.
    pub omega: f64,
    /// Surface tension parameter in the color-gradient model.
    pub surface_tension: f64,
}
impl DropletLBM {
    /// D1Q3 weights: w\[-1\]=1/6, w\[0\]=2/3, w\[+1\]=1/6.
    const WEIGHTS: [f64; 3] = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];
    /// Create a new `DropletLBM` with equilibrium initialisation.
    ///
    /// The droplet (red) phase occupies nodes `[left, right)`.
    pub fn new(n: usize, omega: f64, surface_tension: f64, left: usize, right: usize) -> Self {
        let mut f_red = vec![0.0; 3 * n];
        let mut f_blue = vec![0.0; 3 * n];
        for i in 0..n {
            let rho_r = if i >= left && i < right { 1.0 } else { 0.0 };
            let rho_b = 1.0 - rho_r;
            for k in 0..3 {
                f_red[i * 3 + k] = Self::WEIGHTS[k] * rho_r;
                f_blue[i * 3 + k] = Self::WEIGHTS[k] * rho_b;
            }
        }
        Self {
            f_red,
            f_blue,
            n,
            omega,
            surface_tension,
        }
    }
    /// Total density at node i.
    pub fn density(&self, i: usize) -> f64 {
        let base = i * 3;
        self.f_red[base]
            + self.f_red[base + 1]
            + self.f_red[base + 2]
            + self.f_blue[base]
            + self.f_blue[base + 1]
            + self.f_blue[base + 2]
    }
    /// Red phase density at node i.
    pub fn red_density(&self, i: usize) -> f64 {
        let base = i * 3;
        self.f_red[base] + self.f_red[base + 1] + self.f_red[base + 2]
    }
    /// Phase indicator field φ = (ρ_r - ρ_b) / (ρ_r + ρ_b) ∈ \[-1, 1\].
    pub fn phase_field(&self, i: usize) -> f64 {
        let rho_r = self.red_density(i);
        let rho_b = self.density(i) - rho_r;
        let total = rho_r + rho_b;
        if total < 1e-15 {
            0.0
        } else {
            (rho_r - rho_b) / total
        }
    }
    /// BGK collision + color-gradient recolouring step (simplified).
    pub fn collide(&mut self) {
        for i in 0..self.n {
            let rho = self.density(i);
            let base = i * 3;
            let f0 = self.f_red[base] + self.f_blue[base];
            let f2 = self.f_red[base + 2] + self.f_blue[base + 2];
            let u = if rho > 1e-15 { (-f0 + f2) / rho } else { 0.0 };
            let feq = |w: f64, c: f64| -> f64 {
                w * rho * (1.0 + 3.0 * c * u + 4.5 * (c * u).powi(2) - 1.5 * u * u)
            };
            let feq0 = feq(Self::WEIGHTS[0], -1.0);
            let feq1 = feq(Self::WEIGHTS[1], 0.0);
            let feq2 = feq(Self::WEIGHTS[2], 1.0);
            let total_f = [f0, self.f_red[base + 1] + self.f_blue[base + 1], f2];
            let post = [
                total_f[0] + self.omega * (feq0 - total_f[0]),
                total_f[1] + self.omega * (feq1 - total_f[1]),
                total_f[2] + self.omega * (feq2 - total_f[2]),
            ];
            let rho_r = self.red_density(i);
            let phi = if rho > 1e-15 { rho_r / rho } else { 0.0 };
            let phi_grad = self.color_gradient(i);
            let beta = self.surface_tension;
            for k in 0..3 {
                let ck = [-1.0_f64, 0.0, 1.0][k];
                let anti = beta * phi * (1.0 - phi) * phi_grad * ck;
                self.f_red[base + k] = phi * post[k] + anti;
                self.f_blue[base + k] = (1.0 - phi) * post[k] - anti;
            }
        }
    }
    /// Streaming step (periodic boundaries).
    pub fn stream(&mut self) {
        let n = self.n;
        let mut new_red = vec![0.0; 3 * n];
        let mut new_blue = vec![0.0; 3 * n];
        for i in 0..n {
            let idst_l = if i == 0 { n - 1 } else { i - 1 };
            let idst_r = if i == n - 1 { 0 } else { i + 1 };
            new_red[idst_l * 3] = self.f_red[i * 3];
            new_red[i * 3 + 1] = self.f_red[i * 3 + 1];
            new_red[idst_r * 3 + 2] = self.f_red[i * 3 + 2];
            new_blue[idst_l * 3] = self.f_blue[i * 3];
            new_blue[i * 3 + 1] = self.f_blue[i * 3 + 1];
            new_blue[idst_r * 3 + 2] = self.f_blue[i * 3 + 2];
        }
        self.f_red = new_red;
        self.f_blue = new_blue;
    }
    /// Colour gradient ∂φ/∂x at node i (centred difference, periodic).
    fn color_gradient(&self, i: usize) -> f64 {
        let ip = if i + 1 < self.n { i + 1 } else { 0 };
        let im = if i > 0 { i - 1 } else { self.n - 1 };
        (self.phase_field(ip) - self.phase_field(im)) * 0.5
    }
    /// Advance one full LBM step: collide then stream.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
    }
    /// Total red mass (sum over all nodes).
    pub fn total_red_mass(&self) -> f64 {
        (0..self.n).map(|i| self.red_density(i)).sum()
    }
}
/// A single droplet tracked in the LBM domain.
#[derive(Debug, Clone)]
pub struct Droplet {
    /// Droplet centroid \[x, y, z\] in lattice units.
    pub center: [f64; 3],
    /// Equivalent sphere radius in lattice units.
    pub radius: f64,
    /// Droplet velocity \[vx, vy, vz\] in lattice units / time step.
    pub velocity: [f64; 3],
    /// Surface tension coefficient σ (force per unit length).
    pub surface_tension: f64,
}
impl Droplet {
    /// Create a new `Droplet`.
    pub fn new(center: [f64; 3], radius: f64, velocity: [f64; 3], surface_tension: f64) -> Self {
        Self {
            center,
            radius,
            velocity,
            surface_tension,
        }
    }
    /// Droplet volume (sphere approximation).
    pub fn volume(&self) -> f64 {
        4.0 / 3.0 * PI * self.radius.powi(3)
    }
    /// Speed (magnitude of velocity).
    pub fn speed(&self) -> f64 {
        (self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2)).sqrt()
    }
    /// Weber number We = ρ u² R / σ.
    pub fn weber_number(&self, density: f64) -> f64 {
        if self.surface_tension == 0.0 {
            return f64::INFINITY;
        }
        let u2 = self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2);
        density * u2 * self.radius / self.surface_tension
    }
    /// Capillary number Ca = μ u / σ.
    pub fn capillary_number(&self, viscosity: f64) -> f64 {
        if self.surface_tension == 0.0 {
            return f64::INFINITY;
        }
        viscosity * self.speed() / self.surface_tension
    }
    /// Centre-to-centre distance to another droplet.
    pub fn distance_to(&self, other: &Droplet) -> f64 {
        let dx = self.center[0] - other.center[0];
        let dy = self.center[1] - other.center[1];
        let dz = self.center[2] - other.center[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Gap between surfaces (positive = separated, negative = overlapping).
    pub fn surface_gap(&self, other: &Droplet) -> f64 {
        self.distance_to(other) - self.radius - other.radius
    }
}
/// Spray injection model with configurable size distribution.
#[derive(Debug, Clone)]
pub struct SprayModel {
    /// Mean droplet diameter (log-normal μ_ln or Rosin-Rammler X̄).
    pub mean_diameter: f64,
    /// Standard deviation parameter (log-normal σ_ln or Rosin-Rammler n).
    pub spread_param: f64,
    /// Injection velocity magnitude \[lattice units / time step\].
    pub injection_speed: f64,
    /// Injection direction unit vector.
    pub injection_dir: [f64; 3],
    /// Number of droplets injected per time step.
    pub injection_rate: f64,
    /// Distribution type.
    pub distribution: SizeDistribution,
}
impl SprayModel {
    /// Create a new `SprayModel`.
    pub fn new(
        mean_diameter: f64,
        spread_param: f64,
        injection_speed: f64,
        injection_dir: [f64; 3],
        injection_rate: f64,
        distribution: SizeDistribution,
    ) -> Self {
        Self {
            mean_diameter,
            spread_param,
            injection_speed,
            injection_dir,
            injection_rate,
            distribution,
        }
    }
    /// Log-normal PDF evaluated at diameter d.
    ///
    /// f(d) = 1/(d σ √2π) exp(-(ln d - μ)² / (2σ²))
    pub fn log_normal_pdf(&self, d: f64) -> f64 {
        if d <= 0.0 || self.spread_param <= 0.0 {
            return 0.0;
        }
        let mu = self.mean_diameter.ln();
        let sigma = self.spread_param;
        let exponent = -((d.ln() - mu).powi(2)) / (2.0 * sigma * sigma);
        exponent.exp() / (d * sigma * (2.0 * PI).sqrt())
    }
    /// Rosin-Rammler CDF evaluated at diameter d.
    ///
    /// F(d) = 1 - exp(-(d / X̄)^n)
    pub fn rosin_rammler_cdf(&self, d: f64) -> f64 {
        if d <= 0.0 {
            return 0.0;
        }
        1.0 - (-(d / self.mean_diameter).powf(self.spread_param)).exp()
    }
    /// Mean diameter for the log-normal distribution.
    ///
    /// D_mean = exp(μ + σ²/2)
    pub fn log_normal_mean(&self) -> f64 {
        let mu = self.mean_diameter.ln();
        (mu + self.spread_param * self.spread_param * 0.5).exp()
    }
    /// Sauter mean diameter (SMD) D₃₂ for log-normal distribution.
    ///
    /// D₃₂ = D_mean * exp(5 σ²/2)
    pub fn sauter_mean_diameter(&self) -> f64 {
        let mu = self.mean_diameter.ln();
        (mu + 5.0 * self.spread_param * self.spread_param * 0.5).exp()
    }
    /// Sample injection velocity with a random component in the perpendicular plane.
    ///
    /// Uses a deterministic angle based on `index` for reproducibility.
    pub fn injection_velocity(&self, index: usize, spread_angle: f64) -> [f64; 3] {
        let angle = 2.0 * PI * (index as f64 * 0.618033988749895).fract();
        let perp_mag = self.injection_speed * spread_angle.tan();
        let [dx, dy, dz] = self.injection_dir;
        let perp = if dx.abs() < 0.9 {
            let ex = [1.0_f64, 0.0, 0.0];
            let dot = ex[0] * dx + ex[1] * dy + ex[2] * dz;
            let px = ex[0] - dot * dx;
            let py = ex[1] - dot * dy;
            let pz = ex[2] - dot * dz;
            let len = (px * px + py * py + pz * pz).sqrt().max(1e-15);
            [px / len, py / len, pz / len]
        } else {
            let ey = [0.0_f64, 1.0, 0.0];
            let dot = ey[0] * dx + ey[1] * dy + ey[2] * dz;
            let px = ey[0] - dot * dx;
            let py = ey[1] - dot * dy;
            let pz = ey[2] - dot * dz;
            let len = (px * px + py * py + pz * pz).sqrt().max(1e-15);
            [px / len, py / len, pz / len]
        };
        let (sa, ca) = angle.sin_cos();
        [
            self.injection_speed * dx
                + perp_mag * (ca * perp[0] + sa * (-dz * perp[1] + dy * perp[2])),
            self.injection_speed * dy
                + perp_mag * (ca * perp[1] + sa * (dz * perp[0] - dx * perp[2])),
            self.injection_speed * dz
                + perp_mag * (ca * perp[2] + sa * (-dy * perp[0] + dx * perp[1])),
        ]
    }
    /// Sample droplet diameter for injection parcel `index`.
    ///
    /// Uses inverse-CDF approximation (deterministic, uniform quantile from index).
    pub fn sample_diameter(&self, index: usize) -> f64 {
        let p = ((index as f64 + 0.5) * 0.618033988749895)
            .fract()
            .clamp(1e-6, 1.0 - 1e-6);
        match self.distribution {
            SizeDistribution::LogNormal => {
                let z = probit(p);
                let mu = self.mean_diameter.ln();
                (mu + self.spread_param * z).exp()
            }
            SizeDistribution::RosinRammler => {
                self.mean_diameter * (-((1.0 - p).ln())).powf(1.0 / self.spread_param)
            }
        }
    }
}
/// Connected-components labeling for droplets identified from an order-parameter
/// field (φ > threshold ⇒ droplet phase).
///
/// Uses a single-pass union-find algorithm on a 2-D grid.
#[derive(Debug, Clone)]
pub struct DropletTracking {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Threshold: nodes with φ > threshold belong to a droplet.
    pub threshold: f64,
    /// Label field: label\[k\] = droplet id (0 = background).
    pub labels: Vec<usize>,
    /// Number of distinct droplets found.
    pub n_droplets: usize,
}
impl DropletTracking {
    /// Create a `DropletTracking` from a 2-D order-parameter field.
    pub fn label_droplets(phi: &[f64], nx: usize, ny: usize, threshold: f64) -> Self {
        let n = nx * ny;
        let mut labels = vec![0usize; n];
        let mut parent: Vec<usize> = (0..=n).collect();
        let mut next_label = 1usize;
        for j in 0..ny {
            for i in 0..nx {
                let k = j * nx + i;
                if phi[k] <= threshold {
                    continue;
                }
                let mut nbr_labels: Vec<usize> = Vec::new();
                if i > 0 {
                    let kl = j * nx + i - 1;
                    if labels[kl] > 0 {
                        nbr_labels.push(find(&mut parent, labels[kl]));
                    }
                }
                if j > 0 {
                    let kt = (j - 1) * nx + i;
                    if labels[kt] > 0 {
                        nbr_labels.push(find(&mut parent, labels[kt]));
                    }
                }
                if nbr_labels.is_empty() {
                    labels[k] = next_label;
                    next_label += 1;
                } else {
                    let min_label = nbr_labels.iter().copied().min().unwrap_or(next_label);
                    labels[k] = min_label;
                    for &l in &nbr_labels {
                        union_uf(&mut parent, l, min_label);
                    }
                }
            }
        }
        let mut remap = vec![0usize; next_label];
        let mut new_id = 0usize;
        for lab in labels.iter_mut() {
            if *lab == 0 {
                continue;
            }
            let root = find(&mut parent, *lab);
            if remap[root] == 0 {
                new_id += 1;
                remap[root] = new_id;
            }
            *lab = remap[root];
        }
        Self {
            nx,
            ny,
            threshold,
            labels,
            n_droplets: new_id,
        }
    }
    /// Centroid of droplet `id` (1-indexed).
    pub fn centroid(&self, id: usize) -> Option<[f64; 2]> {
        let mut sx = 0.0_f64;
        let mut sy = 0.0_f64;
        let mut count = 0usize;
        for (k, &lab) in self.labels.iter().enumerate() {
            if lab == id {
                let ix = (k % self.nx) as f64;
                let iy = (k / self.nx) as f64;
                sx += ix;
                sy += iy;
                count += 1;
            }
        }
        if count == 0 {
            None
        } else {
            Some([sx / count as f64, sy / count as f64])
        }
    }
    /// Volume (number of nodes) of droplet `id`.
    pub fn droplet_volume(&self, id: usize) -> usize {
        self.labels.iter().filter(|&&l| l == id).count()
    }
    /// Equivalent circle radius R = sqrt(A / π) from node count.
    pub fn equivalent_radius(&self, id: usize) -> f64 {
        let area = self.droplet_volume(id) as f64;
        (area / PI).sqrt()
    }
    /// Circularity C = 4π A / P² (1 = perfect circle, < 1 = non-circular).
    ///
    /// Perimeter P is estimated by counting interface nodes (nodes with at least
    /// one background neighbour).
    pub fn circularity(&self, id: usize) -> f64 {
        let area = self.droplet_volume(id) as f64;
        if area == 0.0 {
            return 0.0;
        }
        let mut perimeter = 0usize;
        for (k, &lab) in self.labels.iter().enumerate() {
            if lab != id {
                continue;
            }
            let i = k % self.nx;
            let j = k / self.nx;
            let is_border = [
                if i > 0 {
                    self.labels[j * self.nx + i - 1] != id
                } else {
                    true
                },
                if i + 1 < self.nx {
                    self.labels[j * self.nx + i + 1] != id
                } else {
                    true
                },
                if j > 0 {
                    self.labels[(j - 1) * self.nx + i] != id
                } else {
                    true
                },
                if j + 1 < self.ny {
                    self.labels[(j + 1) * self.nx + i] != id
                } else {
                    true
                },
            ];
            if is_border.iter().any(|&b| b) {
                perimeter += 1;
            }
        }
        if perimeter == 0 {
            return 1.0;
        }
        let p = perimeter as f64;
        4.0 * PI * area / (p * p)
    }
}
/// Model for droplet coalescence: neck formation and critical radius criterion.
///
/// Based on the lubrication theory of thin-film drainage between two approaching
/// droplets and the van der Waals disjoining pressure at short range.
#[derive(Debug, Clone)]
pub struct CoalescenceModel {
    /// Surface tension σ \[lattice units\].
    pub surface_tension: f64,
    /// Critical neck radius below which coalescence is irreversible \[lattice units\].
    pub critical_neck_radius: f64,
    /// Film-drainage mobility coefficient (dimensionless).
    pub drainage_mobility: f64,
    /// Hamaker constant A for van der Waals attraction \[lattice units energy\].
    pub hamaker_constant: f64,
}
impl CoalescenceModel {
    /// Create a new `CoalescenceModel`.
    pub fn new(
        surface_tension: f64,
        critical_neck_radius: f64,
        drainage_mobility: f64,
        hamaker_constant: f64,
    ) -> Self {
        Self {
            surface_tension,
            critical_neck_radius,
            drainage_mobility,
            hamaker_constant,
        }
    }
    /// Neck growth rate from inviscid Eggers–Villermaux model.
    ///
    /// dr/dt ≈ (σ / ρ r)^(1/2) for inertial regime.
    pub fn neck_growth_rate_inertial(&self, neck_radius: f64, density: f64) -> f64 {
        if neck_radius <= 0.0 || density <= 0.0 {
            return 0.0;
        }
        (self.surface_tension / (density * neck_radius)).sqrt()
    }
    /// Neck growth rate in the viscous regime: dr/dt ≈ σ / (6 π μ ln(r₀/r)).
    pub fn neck_growth_rate_viscous(&self, neck_radius: f64, r0: f64, viscosity: f64) -> f64 {
        if neck_radius <= 0.0 || r0 <= 0.0 || viscosity <= 0.0 {
            return 0.0;
        }
        let ln_ratio = (r0 / neck_radius).ln().max(1e-10);
        self.surface_tension / (6.0 * PI * viscosity * ln_ratio)
    }
    /// van der Waals disjoining pressure between two flat films separated by gap h.
    ///
    /// Π_vdW = -A / (6π h³)
    pub fn disjoining_pressure(&self, gap: f64) -> f64 {
        if gap <= 0.0 {
            return 0.0;
        }
        -self.hamaker_constant / (6.0 * PI * gap.powi(3))
    }
    /// Returns `true` if the neck radius is below the critical value.
    pub fn neck_is_critical(&self, neck_radius: f64) -> bool {
        neck_radius <= self.critical_neck_radius
    }
    /// Merged droplet radius conserving volume (3-D spheres).
    pub fn merged_radius(&self, r1: f64, r2: f64) -> f64 {
        (r1.powi(3) + r2.powi(3)).cbrt()
    }
    /// Film drainage timescale τ_drain = h₀² / (drainage_mobility · σ).
    pub fn drainage_timescale(&self, initial_gap: f64) -> f64 {
        if self.surface_tension <= 0.0 || self.drainage_mobility <= 0.0 {
            return f64::INFINITY;
        }
        initial_gap * initial_gap / (self.drainage_mobility * self.surface_tension)
    }
}
/// Lift and drag force models for a droplet moving through a carrier fluid.
#[derive(Debug, Clone)]
pub struct LiftDragDroplet {
    /// Continuous-phase density ρ_c.
    pub fluid_density: f64,
    /// Continuous-phase dynamic viscosity μ_c.
    pub fluid_viscosity: f64,
    /// Relative velocity \[u_rel_x, u_rel_y, u_rel_z\] = u_fluid - u_droplet.
    pub rel_velocity: [f64; 3],
    /// Ambient fluid vorticity \[ω_x, ω_y, ω_z\] (for Saffman / Magnus lift).
    pub vorticity: [f64; 3],
    /// Droplet spin \[Ω_x, Ω_y, Ω_z\] (rad / time step, for Magnus lift).
    pub spin: [f64; 3],
}
impl LiftDragDroplet {
    /// Create a new `LiftDragDroplet`.
    pub fn new(
        fluid_density: f64,
        fluid_viscosity: f64,
        rel_velocity: [f64; 3],
        vorticity: [f64; 3],
        spin: [f64; 3],
    ) -> Self {
        Self {
            fluid_density,
            fluid_viscosity,
            rel_velocity,
            vorticity,
            spin,
        }
    }
    /// Particle (droplet) Reynolds number Re_p = ρ_c |u_rel| d / μ_c.
    pub fn particle_reynolds(&self, radius: f64) -> f64 {
        let u = (self.rel_velocity[0].powi(2)
            + self.rel_velocity[1].powi(2)
            + self.rel_velocity[2].powi(2))
        .sqrt();
        self.fluid_density * u * 2.0 * radius / self.fluid_viscosity.max(1e-30)
    }
    /// Schiller-Naumann drag coefficient C_D.
    ///
    /// C_D = 24/Re * (1 + 0.15 Re^0.687)  for Re ≤ 1000
    /// C_D = 0.44                           for Re > 1000
    pub fn schiller_naumann_cd(&self, radius: f64) -> f64 {
        let re = self.particle_reynolds(radius);
        if re < 1e-10 {
            return 0.0;
        }
        if re <= 1000.0 {
            24.0 / re * (1.0 + 0.15 * re.powf(0.687))
        } else {
            0.44
        }
    }
    /// Schiller-Naumann drag force vector F_D.
    ///
    /// F_D = ½ ρ_c |u_rel| u_rel π r² C_D
    pub fn drag_force(&self, radius: f64) -> [f64; 3] {
        let cd = self.schiller_naumann_cd(radius);
        let u_mag = (self.rel_velocity[0].powi(2)
            + self.rel_velocity[1].powi(2)
            + self.rel_velocity[2].powi(2))
        .sqrt();
        let area = PI * radius * radius;
        let coeff = 0.5 * self.fluid_density * u_mag * cd * area;
        [
            coeff * self.rel_velocity[0],
            coeff * self.rel_velocity[1],
            coeff * self.rel_velocity[2],
        ]
    }
    /// Saffman lift coefficient C_L (shear-induced).
    ///
    /// C_L = 1.615 * sqrt(Re_s) / Re_p  (simplified)
    ///
    /// where Re_s = ρ_c |ω| d² / μ_c is the shear Reynolds number.
    pub fn saffman_lift_coeff(&self, radius: f64) -> f64 {
        let re_p = self.particle_reynolds(radius).max(1e-10);
        let omega_mag =
            (self.vorticity[0].powi(2) + self.vorticity[1].powi(2) + self.vorticity[2].powi(2))
                .sqrt();
        let diameter = 2.0 * radius;
        let re_s =
            self.fluid_density * omega_mag * diameter * diameter / self.fluid_viscosity.max(1e-30);
        1.615 * re_s.sqrt() / re_p
    }
    /// Saffman lift force vector F_L = C_L ½ ρ_c |u_rel|² π r² * (ω̂ × û_rel).
    pub fn saffman_lift(&self, radius: f64) -> [f64; 3] {
        let cl = self.saffman_lift_coeff(radius);
        let u_mag = (self.rel_velocity[0].powi(2)
            + self.rel_velocity[1].powi(2)
            + self.rel_velocity[2].powi(2))
        .sqrt();
        let cx =
            self.vorticity[1] * self.rel_velocity[2] - self.vorticity[2] * self.rel_velocity[1];
        let cy =
            self.vorticity[2] * self.rel_velocity[0] - self.vorticity[0] * self.rel_velocity[2];
        let cz =
            self.vorticity[0] * self.rel_velocity[1] - self.vorticity[1] * self.rel_velocity[0];
        let area = PI * radius * radius;
        let coeff = cl * 0.5 * self.fluid_density * u_mag * area;
        [coeff * cx, coeff * cy, coeff * cz]
    }
    /// Magnus lift force due to droplet spin: F_M = π r³ ρ_c (Ω × u_rel).
    pub fn magnus_lift(&self, radius: f64) -> [f64; 3] {
        let cx = self.spin[1] * self.rel_velocity[2] - self.spin[2] * self.rel_velocity[1];
        let cy = self.spin[2] * self.rel_velocity[0] - self.spin[0] * self.rel_velocity[2];
        let cz = self.spin[0] * self.rel_velocity[1] - self.spin[1] * self.rel_velocity[0];
        let coeff = PI * radius.powi(3) * self.fluid_density;
        [coeff * cx, coeff * cy, coeff * cz]
    }
    /// Total force = drag + Saffman lift + Magnus lift.
    pub fn total_force(&self, radius: f64) -> [f64; 3] {
        let fd = self.drag_force(radius);
        let fs = self.saffman_lift(radius);
        let fm = self.magnus_lift(radius);
        [
            fd[0] + fs[0] + fm[0],
            fd[1] + fs[1] + fm[1],
            fd[2] + fs[2] + fm[2],
        ]
    }
}
/// Cahn-Hilliard LBM for a diffuse-interface two-phase system on a D2Q9 grid.
///
/// The order parameter φ ∈ \[-1, +1\] marks the phases (+1 = droplet, -1 = carrier).
/// The chemical potential μ = -ε² ∇²φ + f'(φ) drives the phase-field evolution.
///
/// A simplified 1-D finite-difference solver is used for the Laplacians; the
/// LBM streaming and collision are handled by the [`DropletLBM`] type.
#[derive(Debug, Clone)]
pub struct CahnHilliardLbm {
    /// Order parameter field φ_i ∈ \[-1, +1\], length `nx * ny`.
    pub phi: Vec<f64>,
    /// Chemical potential field μ_i, length `nx * ny`.
    pub mu: Vec<f64>,
    /// Grid width (columns).
    pub nx: usize,
    /// Grid height (rows).
    pub ny: usize,
    /// Lattice spacing Δx \[lattice units\].
    pub dx: f64,
    /// Interface width parameter ε \[lattice units\].
    pub epsilon: f64,
    /// Mobility M \[lattice units² / time step\].
    pub mobility: f64,
}
impl CahnHilliardLbm {
    /// Create a new `CahnHilliardLbm` on an `nx × ny` lattice with a
    /// centered circular droplet of radius `r_droplet`.
    pub fn new_circular_droplet(
        nx: usize,
        ny: usize,
        r_droplet: f64,
        epsilon: f64,
        mobility: f64,
    ) -> Self {
        let dx = 1.0_f64;
        let cx = nx as f64 * 0.5;
        let cy = ny as f64 * 0.5;
        let phi: Vec<f64> = (0..nx * ny)
            .map(|k| {
                let ix = (k % nx) as f64;
                let iy = (k / nx) as f64;
                let r = ((ix - cx).powi(2) + (iy - cy).powi(2)).sqrt();
                ((r_droplet - r) / (epsilon.max(1e-14) * 2.0_f64.sqrt())).tanh()
            })
            .collect();
        Self {
            phi,
            mu: vec![0.0; nx * ny],
            nx,
            ny,
            dx,
            epsilon,
            mobility,
        }
    }
    /// Flat index (row-major) for grid node (i, j).
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }
    /// Laplacian of φ at node (i, j) using second-order FD with periodic boundaries.
    pub fn laplacian_phi(&self, i: usize, j: usize) -> f64 {
        let ip = if i + 1 < self.nx { i + 1 } else { 0 };
        let im = if i > 0 { i - 1 } else { self.nx - 1 };
        let jp = if j + 1 < self.ny { j + 1 } else { 0 };
        let jm = if j > 0 { j - 1 } else { self.ny - 1 };
        let c = self.phi[self.idx(i, j)];
        let dx2 = self.dx * self.dx;
        (self.phi[self.idx(ip, j)]
            + self.phi[self.idx(im, j)]
            + self.phi[self.idx(i, jp)]
            + self.phi[self.idx(i, jm)]
            - 4.0 * c)
            / dx2
    }
    /// Compute the chemical potential field μ = -ε² ∇²φ + f'(φ).
    pub fn compute_mu(&mut self) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let k = self.idx(i, j);
                let phi_k = self.phi[k];
                let lap = self.laplacian_phi(i, j);
                self.mu[k] = -self.epsilon * self.epsilon * lap + phi_k * (phi_k * phi_k - 1.0);
            }
        }
    }
    /// One explicit time step of the Cahn-Hilliard equation: ∂φ/∂t = M ∇²μ.
    pub fn step(&mut self, dt: f64) {
        self.compute_mu();
        let mut new_phi = self.phi.clone();
        for j in 0..self.ny {
            for i in 0..self.nx {
                let ip = if i + 1 < self.nx { i + 1 } else { 0 };
                let im = if i > 0 { i - 1 } else { self.nx - 1 };
                let jp = if j + 1 < self.ny { j + 1 } else { 0 };
                let jm = if j > 0 { j - 1 } else { self.ny - 1 };
                let k = self.idx(i, j);
                let lap_mu = (self.mu[self.idx(ip, j)]
                    + self.mu[self.idx(im, j)]
                    + self.mu[self.idx(i, jp)]
                    + self.mu[self.idx(i, jm)]
                    - 4.0 * self.mu[k])
                    / (self.dx * self.dx);
                new_phi[k] = (self.phi[k] + dt * self.mobility * lap_mu).clamp(-1.5, 1.5);
            }
        }
        self.phi = new_phi;
    }
    /// Total order-parameter (conserved quantity): Σ_i φ_i.
    pub fn total_order_parameter(&self) -> f64 {
        self.phi.iter().sum()
    }
    /// Volume fraction of the droplet phase (φ > 0 fraction).
    pub fn droplet_volume_fraction(&self) -> f64 {
        let n_pos = self.phi.iter().filter(|&&p| p > 0.0).count();
        n_pos as f64 / self.phi.len() as f64
    }
    /// Bulk free energy F = Σ_i \[ε²/2 |∇φ|² + (1/4)(φ²-1)²\] Δx².
    pub fn bulk_free_energy(&self) -> f64 {
        let mut f = 0.0_f64;
        for j in 0..self.ny {
            for i in 0..self.nx {
                let k = self.idx(i, j);
                let phi_k = self.phi[k];
                let ip = if i + 1 < self.nx { i + 1 } else { 0 };
                let jp = if j + 1 < self.ny { j + 1 } else { 0 };
                let dphi_dx = (self.phi[self.idx(ip, j)] - phi_k) / self.dx;
                let dphi_dy = (self.phi[self.idx(i, jp)] - phi_k) / self.dx;
                let grad2 = dphi_dx * dphi_dx + dphi_dy * dphi_dy;
                f += 0.5 * self.epsilon * self.epsilon * grad2
                    + 0.25 * (phi_k * phi_k - 1.0).powi(2);
            }
        }
        f * self.dx * self.dx
    }
}
/// D² evaporation law and mass transfer for a single droplet.
///
/// Based on the classical D² law: d(D²)/dt = -K where K is the evaporation
/// constant, combined with a Sherwood-number correction.
#[derive(Debug, Clone)]
pub struct EvaporationModel {
    /// Initial droplet diameter squared D₀² \[lattice units²\].
    pub d2_initial: f64,
    /// Evaporation constant K \[lattice units² / time step\].
    pub evaporation_constant: f64,
    /// Vapor diffusion coefficient D_v \[lattice units² / time step\].
    pub vapor_diffusivity: f64,
    /// Ambient vapor mass fraction Y_inf.
    pub vapor_mass_fraction_inf: f64,
    /// Surface (saturation) vapor mass fraction Y_s.
    pub vapor_mass_fraction_surf: f64,
    /// Droplet density ρ_l.
    pub droplet_density: f64,
}
impl EvaporationModel {
    /// Create a new `EvaporationModel`.
    pub fn new(
        d2_initial: f64,
        evaporation_constant: f64,
        vapor_diffusivity: f64,
        vapor_mass_fraction_inf: f64,
        vapor_mass_fraction_surf: f64,
        droplet_density: f64,
    ) -> Self {
        Self {
            d2_initial,
            evaporation_constant,
            vapor_diffusivity,
            vapor_mass_fraction_inf,
            vapor_mass_fraction_surf,
            droplet_density,
        }
    }
    /// D² at time t using D²(t) = D₀² - K t.
    ///
    /// Returns 0 when fully evaporated.
    pub fn d2_at_time(&self, t: f64) -> f64 {
        (self.d2_initial - self.evaporation_constant * t).max(0.0)
    }
    /// Droplet diameter at time t.
    pub fn diameter_at_time(&self, t: f64) -> f64 {
        self.d2_at_time(t).sqrt()
    }
    /// Droplet lifetime (time for complete evaporation).
    pub fn lifetime(&self) -> f64 {
        if self.evaporation_constant <= 0.0 {
            return f64::INFINITY;
        }
        self.d2_initial / self.evaporation_constant
    }
    /// Spalding transfer number B_M = (Y_s - Y_inf) / (1 - Y_s).
    pub fn spalding_number(&self) -> f64 {
        let denom = 1.0 - self.vapor_mass_fraction_surf;
        if denom.abs() < 1e-15 {
            return f64::INFINITY;
        }
        (self.vapor_mass_fraction_surf - self.vapor_mass_fraction_inf) / denom
    }
    /// Evaporation rate corrected by the Sherwood number Sh.
    ///
    /// ṁ = -π D D_v ρ Sh ln(1 + B_M)
    ///
    /// Uses correlation Sh = 2 + 0.6 Re^0.5 Sc^0.33 (Ranz-Marshall).
    pub fn mass_transfer_rate(&self, diameter: f64, reynolds: f64, schmidt: f64) -> f64 {
        let sh = sherwood_number(reynolds, schmidt);
        let bm = self.spalding_number();
        if bm <= 0.0 || diameter <= 0.0 {
            return 0.0;
        }
        -PI * diameter * self.vapor_diffusivity * self.droplet_density * sh * (1.0 + bm).ln()
    }
    /// Vapor pressure at surface from Clausius-Clapeyron approximation.
    ///
    /// p_vap = p_ref * exp(-L/(R T) * (1/T - 1/T_ref))
    ///
    /// Here we use a simplified form: p_vap = exp(-l_v * (1/T - 1/t_ref)).
    pub fn vapor_pressure(temp: f64, latent_heat: f64, t_ref: f64, p_ref: f64) -> f64 {
        if temp <= 0.0 || t_ref <= 0.0 {
            return 0.0;
        }
        p_ref * (-latent_heat * (1.0 / temp - 1.0 / t_ref)).exp()
    }
}
/// Breakup models for droplets in turbulent / high-Weber-number flows.
#[derive(Debug, Clone)]
pub struct DropletBreakup {
    /// Critical Weber number We* for secondary breakup.
    pub critical_weber: f64,
    /// Atwood number At = (ρ_d - ρ_c) / (ρ_d + ρ_c) for R-T instability.
    pub atwood_number: f64,
    /// Gravity acceleration magnitude g \[lattice units\].
    pub gravity: f64,
}
impl DropletBreakup {
    /// Create a new `DropletBreakup`.
    pub fn new(critical_weber: f64, atwood_number: f64, gravity: f64) -> Self {
        Self {
            critical_weber,
            atwood_number,
            gravity,
        }
    }
    /// Weber-number breakup criterion.
    ///
    /// Returns `true` when We > We*.
    pub fn weber_breakup(&self, droplet: &Droplet, density: f64) -> bool {
        droplet.weber_number(density) > self.critical_weber
    }
    /// Rayleigh-Taylor instability growth rate.
    ///
    /// σ_RT = sqrt(At · g · k - σ k³ / (ρ_d + ρ_c))
    ///
    /// with the most-unstable wavenumber k* = sqrt(At g (ρ_d+ρ_c) / (3σ)).
    ///
    /// Here we use a simplified form using the droplet radius:
    /// k* ≈ 1/R.
    pub fn rt_growth_rate(&self, droplet: &Droplet, density_sum: f64) -> f64 {
        if droplet.surface_tension <= 0.0 || density_sum <= 0.0 {
            return 0.0;
        }
        let k = 1.0 / droplet.radius.max(1e-14);
        let term1 = self.atwood_number * self.gravity * k;
        let term2 = droplet.surface_tension * k.powi(3) / density_sum;
        let val = term1 - term2;
        if val > 0.0 { val.sqrt() } else { 0.0 }
    }
    /// Returns `true` if Rayleigh-Taylor instability drives breakup.
    ///
    /// Breakup occurs when the RT growth rate exceeds a threshold (here
    /// 0.1 in lattice units).
    pub fn rt_breakup(&self, droplet: &Droplet, density_sum: f64) -> bool {
        self.rt_growth_rate(droplet, density_sum) > 0.1
    }
    /// Daughter droplet radius after symmetric binary breakup
    /// (volume conservation): r_child = R / 2^(1/3).
    pub fn child_radius(parent_radius: f64) -> f64 {
        parent_radius / 2.0_f64.cbrt()
    }
}
/// Criterion controlling whether two neighbouring droplets merge.
#[derive(Debug, Clone)]
pub struct CoalescenceCriterion {
    /// Critical surface-to-surface gap below which coalescence is possible.
    pub critical_distance: f64,
    /// Film-drainage time constant τ_drain (time steps).
    pub film_drainage_time: f64,
    /// Hamaker constant A (van der Waals) in lattice-unit energy units.
    pub hamaker_constant: f64,
}
impl CoalescenceCriterion {
    /// Create a new `CoalescenceCriterion`.
    pub fn new(critical_distance: f64, film_drainage_time: f64, hamaker_constant: f64) -> Self {
        Self {
            critical_distance,
            film_drainage_time,
            hamaker_constant,
        }
    }
    /// Van der Waals disjoining pressure between two flat interfaces
    /// separated by gap h: Π = -A / (6π h³).
    ///
    /// Returns 0 if h ≤ 0.
    pub fn van_der_waals_pressure(&self, gap: f64) -> f64 {
        if gap <= 0.0 {
            return 0.0;
        }
        -self.hamaker_constant / (6.0 * PI * gap.powi(3))
    }
    /// Film drainage timescale scaled by capillary number.
    ///
    /// τ_eff = τ_drain * (1 + Ca)
    pub fn effective_drain_time(&self, capillary_number: f64) -> f64 {
        self.film_drainage_time * (1.0 + capillary_number)
    }
    /// Returns `true` if two droplets satisfy the coalescence criterion.
    pub fn should_coalesce(&self, a: &Droplet, b: &Droplet) -> bool {
        a.surface_gap(b) < self.critical_distance
    }
    /// Merge two droplets conserving volume and momentum.
    pub fn merge(&self, a: &Droplet, b: &Droplet) -> Droplet {
        let va = a.volume();
        let vb = b.volume();
        let vt = va + vb;
        let r_new = (3.0 * vt / (4.0 * PI)).cbrt();
        let cx = (va * a.center[0] + vb * b.center[0]) / vt;
        let cy = (va * a.center[1] + vb * b.center[1]) / vt;
        let cz = (va * a.center[2] + vb * b.center[2]) / vt;
        let vx = (va * a.velocity[0] + vb * b.velocity[0]) / vt;
        let vy = (va * a.velocity[1] + vb * b.velocity[1]) / vt;
        let vz = (va * a.velocity[2] + vb * b.velocity[2]) / vt;
        let sigma = (a.surface_tension + b.surface_tension) * 0.5;
        Droplet::new([cx, cy, cz], r_new, [vx, vy, vz], sigma)
    }
}
/// Physical parameters for an LBM droplet dynamics simulation.
///
/// Encapsulates density ratio, viscosity ratio, surface tension coefficient,
/// and the diffuse-interface width (in lattice units) used by the
/// Cahn-Hilliard and color-gradient LBM models.
#[derive(Debug, Clone)]
pub struct DropletLbmParams {
    /// Density of the droplet (heavy) phase ρ_d \[lattice units\].
    pub density_droplet: f64,
    /// Density of the carrier (light) phase ρ_c \[lattice units\].
    pub density_carrier: f64,
    /// Dynamic viscosity of the droplet phase μ_d \[lattice units\].
    pub viscosity_droplet: f64,
    /// Dynamic viscosity of the carrier phase μ_c \[lattice units\].
    pub viscosity_carrier: f64,
    /// Surface-tension coefficient σ \[force per unit length in lattice units\].
    pub surface_tension: f64,
    /// Diffuse-interface half-width W \[lattice units\] (Cahn number = W/L).
    pub interface_width: f64,
}
impl DropletLbmParams {
    /// Create a new `DropletLbmParams`.
    pub fn new(
        density_droplet: f64,
        density_carrier: f64,
        viscosity_droplet: f64,
        viscosity_carrier: f64,
        surface_tension: f64,
        interface_width: f64,
    ) -> Self {
        Self {
            density_droplet,
            density_carrier,
            viscosity_droplet,
            viscosity_carrier,
            surface_tension,
            interface_width,
        }
    }
    /// Density ratio Λ = ρ_d / ρ_c.
    pub fn density_ratio(&self) -> f64 {
        self.density_droplet / self.density_carrier.max(1e-30)
    }
    /// Viscosity ratio λ = μ_d / μ_c.
    pub fn viscosity_ratio(&self) -> f64 {
        self.viscosity_droplet / self.viscosity_carrier.max(1e-30)
    }
    /// BGK relaxation parameter for the droplet phase.
    ///
    /// τ_d = 3 μ_d / ρ_d + 0.5
    pub fn tau_droplet(&self) -> f64 {
        3.0 * self.viscosity_droplet / self.density_droplet.max(1e-30) + 0.5
    }
    /// BGK relaxation parameter for the carrier phase.
    ///
    /// τ_c = 3 μ_c / ρ_c + 0.5
    pub fn tau_carrier(&self) -> f64 {
        3.0 * self.viscosity_carrier / self.density_carrier.max(1e-30) + 0.5
    }
    /// Cahn number Cn = W / L (dimensionless interface width relative to domain L).
    pub fn cahn_number(&self, domain_length: f64) -> f64 {
        self.interface_width / domain_length.max(1e-30)
    }
    /// Capillary number Ca = μ_c U / σ.
    pub fn capillary_number(&self, velocity_scale: f64) -> f64 {
        self.viscosity_carrier * velocity_scale / self.surface_tension.max(1e-30)
    }
    /// Weber number We = ρ_c U² L / σ.
    pub fn weber_number(&self, velocity_scale: f64, length_scale: f64) -> f64 {
        self.density_carrier * velocity_scale * velocity_scale * length_scale
            / self.surface_tension.max(1e-30)
    }
}
