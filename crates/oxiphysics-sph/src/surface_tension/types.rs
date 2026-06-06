//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::kernel::CubicSplineKernel;
use crate::kernel::SphKernel;
use oxiphysics_core::math::Vec3;

/// Sessile droplet equilibrium shape (spherical cap approximation).
///
/// A droplet of given volume resting on a flat substrate in equilibrium adopts
/// a spherical cap shape characterised solely by the static contact angle θ.
///
/// The geometry uses:
/// - Contact radius: `r_c = (3V/π f(θ))^(1/3)`  where `f(θ) = (2−3cosθ+cos³θ)/sin³θ`
/// - Cap height: `h = r_c (1 − cosθ) / sinθ`
/// - Sphere radius of curvature: `R = r_c / sinθ`
#[derive(Debug, Clone)]
pub struct SessileDroplet {
    /// Droplet volume V (m³).
    pub volume: f64,
    /// Fluid density ρ (kg/m³).
    pub density: f64,
    /// Surface tension coefficient σ (N/m).
    pub sigma: f64,
    /// Static contact angle θ (radians).
    pub contact_angle: f64,
}
impl SessileDroplet {
    /// Create a sessile droplet with the given parameters.
    pub fn new(volume: f64, density: f64, sigma: f64, contact_angle: f64) -> Self {
        Self {
            volume,
            density,
            sigma,
            contact_angle,
        }
    }
    /// f(θ) = (2 − 3cosθ + cos³θ) / sin³θ  (spherical cap volume shape factor).
    fn shape_factor(&self) -> f64 {
        let ct = self.contact_angle.cos();
        let st = self.contact_angle.sin();
        if st.abs() < 1e-14 {
            return if ct > 0.0 { 1e30 } else { 0.0 };
        }
        (2.0 - 3.0 * ct + ct.powi(3)) / st.powi(3)
    }
    /// Equilibrium contact radius r_c (m).
    pub fn contact_width(&self) -> f64 {
        let f = self.shape_factor();
        if f < 1e-30 {
            return 0.0;
        }
        (3.0 * self.volume / (std::f64::consts::PI * f)).cbrt()
    }
    /// Droplet cap height h = r_c (1 − cosθ) / sinθ (m).
    pub fn height(&self) -> f64 {
        let rc = self.contact_width();
        let st = self.contact_angle.sin();
        if st.abs() < 1e-14 {
            return 0.0;
        }
        rc * (1.0 - self.contact_angle.cos()) / st
    }
    /// Sphere radius of curvature R = r_c / sinθ (m).
    pub fn sphere_radius(&self) -> f64 {
        let rc = self.contact_width();
        let st = self.contact_angle.sin();
        if st.abs() < 1e-14 {
            return f64::INFINITY;
        }
        rc / st
    }
    /// Young–Laplace pressure across the interface: ΔP = 2σ/R (Pa).
    pub fn laplace_pressure(&self) -> f64 {
        let r = self.sphere_radius();
        if r.is_infinite() || r < 1e-30 {
            return 0.0;
        }
        2.0 * self.sigma / r
    }
    /// Bond number for this droplet: Bo = ρ g R² / σ.
    pub fn bond_number(&self, g: f64) -> f64 {
        let r = self.sphere_radius();
        bond_number(self.density, g, r, self.sigma)
    }
    /// Droplet mass m = ρ V (kg).
    pub fn mass(&self) -> f64 {
        self.density * self.volume
    }
}
/// Dynamic contact angle model (Cox–Voinov / Tanner law).
///
/// The dynamic contact angle θ_d depends on the contact-line capillary number:
/// ```text
/// cos θ_d = cos θ_s − 2 A Ca^n
/// ```
/// where θ_s is the static contact angle, `A` is an empirical prefactor, and
/// `n` is an exponent (typically 1/3 for the Tanner law).
///
/// Advancing motion (Ca > 0) increases θ; receding (Ca < 0) decreases θ.
#[derive(Debug, Clone)]
pub struct DynamicContactAngle {
    /// Static contact angle θ_s (radians).
    pub theta_static: f64,
    /// Empirical prefactor A in the Cox–Voinov correction.
    pub prefactor_a: f64,
    /// Surface tension coefficient σ (N/m), needed to compute Ca.
    pub sigma: f64,
}
impl DynamicContactAngle {
    /// Create a new dynamic contact angle model.
    pub fn new(theta_static: f64, prefactor_a: f64, sigma: f64) -> Self {
        Self {
            theta_static,
            prefactor_a,
            sigma,
        }
    }
    /// Dynamic contact angle (radians) given the contact-line capillary number.
    ///
    /// `ca` = μ U_cl / σ (signed: positive = advancing, negative = receding).
    pub fn angle(&self, ca: f64) -> f64 {
        let cos_d = self.theta_static.cos() - 2.0 * self.prefactor_a * ca;
        cos_d.clamp(-1.0, 1.0).acos()
    }
    /// Advancing angle (Ca > 0): maximum expected dynamic angle.
    pub fn advancing_angle(&self, ca_magnitude: f64) -> f64 {
        self.angle(ca_magnitude.abs())
    }
    /// Receding angle (Ca < 0): minimum expected dynamic angle.
    pub fn receding_angle(&self, ca_magnitude: f64) -> f64 {
        self.angle(-ca_magnitude.abs())
    }
    /// Contact-angle hysteresis: θ_adv − θ_rec (radians).
    pub fn hysteresis(&self, ca_magnitude: f64) -> f64 {
        self.advancing_angle(ca_magnitude) - self.receding_angle(ca_magnitude)
    }
}
/// Extended CSF model that applies a renormalization correction to the
/// color-field gradient.
///
/// The standard SPH color-gradient estimator is first-order consistent only
/// on a complete kernel support.  Near interfaces or free surfaces, the
/// support is truncated and the gradient is under-estimated.  This struct
/// applies the *renormalized* (or *corrected*) gradient:
///
/// ```text
/// ñ_i = L_i⁻¹ · n_i
/// ```
///
/// where `L_i` is the first-order consistency (renormalization) matrix.
#[derive(Debug, Clone)]
pub struct CsfModel {
    /// Surface tension coefficient σ (N/m).
    pub sigma: f64,
    /// Smoothing length h.
    pub smoothing_length: f64,
}
impl CsfModel {
    /// Create a new `CsfModel`.
    pub fn new(sigma: f64, smoothing_length: f64) -> Self {
        Self {
            sigma,
            smoothing_length,
        }
    }
    /// Compute the renormalized (corrected) color-field gradient.
    ///
    /// Each raw gradient component is divided by the corresponding diagonal
    /// element of the renormalization matrix `L_i`:
    ///
    /// ```text
    /// L_i^{ab} = Σ_j (m_j / ρ_j) * (r_j - r_i)^a * ∂W_ij / ∂r^b
    /// ```
    ///
    /// For efficiency only the diagonal is used here (a simplified but robust
    /// correction adequate for near-isotropic particle distributions).
    ///
    /// # Arguments
    /// - `positions`   : particle positions.
    /// - `masses`      : particle masses.
    /// - `densities`   : particle densities.
    /// - `color`       : scalar color field c_i (e.g. from `compute_color_field`).
    /// - `neighbors`   : neighbour lists.
    ///
    /// # Returns
    /// Corrected color gradient vectors.
    pub fn compute_color_gradient_correction(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        color: &[f64],
        neighbors: &[Vec<usize>],
    ) -> Vec<Vec3> {
        let n = positions.len();
        let h = self.smoothing_length;
        let kernel = CubicSplineKernel;
        let mut raw_grad = vec![Vec3::zeros(); n];
        for i in 0..n {
            if densities[i] < 1e-14 {
                continue;
            }
            for &j in &neighbors[i] {
                let rhoj = densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_w = kernel.grad_w(r, h);
                let rhat = rij / r;
                let dc = color[j] - color[i];
                raw_grad[i] += masses[j] / rhoj * dc * grad_w * rhat;
            }
        }
        let mut l_diag = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            if densities[i] < 1e-14 {
                continue;
            }
            for &j in &neighbors[i] {
                let rhoj = densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_w = kernel.grad_w(r, h);
                let rhat = rij / r;
                let fac = masses[j] / rhoj * grad_w;
                let dr = positions[j] - positions[i];
                for a in 0..3 {
                    l_diag[i][a] += fac * dr[a] * rhat[a];
                }
            }
        }
        let mut corrected = vec![Vec3::zeros(); n];
        for i in 0..n {
            let g = raw_grad[i];
            let l = l_diag[i];
            corrected[i] = Vec3::new(
                if l[0].abs() > 1e-14 { g.x / l[0] } else { g.x },
                if l[1].abs() > 1e-14 { g.y / l[1] } else { g.y },
                if l[2].abs() > 1e-14 { g.z / l[2] } else { g.z },
            );
        }
        corrected
    }
}
/// Pairwise surface tension model (Akinci et al. 2013).
///
/// Computes inter-particle cohesion and curvature minimization forces
/// without requiring explicit interface tracking.
#[derive(Debug, Clone)]
pub struct PairwiseSurfaceTension {
    /// Surface tension coefficient γ (N/m).
    pub gamma: f64,
    /// Smoothing length h.
    pub smoothing_length: f64,
}
impl PairwiseSurfaceTension {
    /// Create a new pairwise surface tension model.
    pub fn new(gamma: f64, smoothing_length: f64) -> Self {
        Self {
            gamma,
            smoothing_length,
        }
    }
    /// Compute the cohesion kernel C(r) used by Akinci's model.
    ///
    /// `C(r) = (32 / (pi h^9)) * (h-r)^3 * r^3`  for h/2 <= r <= h
    /// `C(r) = (32 / (pi h^9)) * (2*(h-r)^3 * r^3 - h^6/64)` for 0 <= r < h/2
    pub(crate) fn cohesion_kernel(&self, r: f64) -> f64 {
        let h = self.smoothing_length;
        if r >= h || r < 0.0 {
            return 0.0;
        }
        let coeff = 32.0 / (std::f64::consts::PI * h.powi(9));
        if r >= h / 2.0 {
            coeff * (h - r).powi(3) * r.powi(3)
        } else {
            coeff * (2.0 * (h - r).powi(3) * r.powi(3) - h.powi(6) / 64.0)
        }
    }
    /// Compute pairwise surface tension forces.
    ///
    /// For each pair (i, j):
    /// `f_ij = -gamma * m_i * m_j * C(|r_ij|) * r_hat_ij / |r_ij|`
    /// plus a curvature-minimization term.
    pub fn compute_forces(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
    ) -> Vec<Vec3> {
        let n = positions.len();
        let mut forces = vec![Vec3::zeros(); n];
        for i in 0..n {
            let rhoi = densities[i];
            if rhoi < 1e-14 {
                continue;
            }
            for &j in &neighbors[i] {
                let rhoj = densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let c = self.cohesion_kernel(r);
                let rhat = rij / r;
                let k_ij = 2.0 * densities[i].max(1e-14) / (densities[i] + densities[j]).max(1e-14);
                let f_cohesion = -self.gamma * masses[i] * masses[j] * c * rhat;
                let f_curvature =
                    -self.gamma * masses[i] * (densities[i] - densities[j]) * rhat * 0.001;
                forces[i] += k_ij * (f_cohesion + f_curvature);
            }
        }
        forces
    }
}
impl PairwiseSurfaceTension {
    /// Compute Tartakovsky & Meakin (2005) cohesion forces.
    ///
    /// The cohesion potential between particles i and j is
    ///
    /// ```text
    /// F_ij = -κ * m_i * m_j / (ρ_i * ρ_j) * C(r) * r̂_ij
    /// ```
    ///
    /// where `C(r)` is the cohesion kernel and `κ` is a cohesion strength
    /// parameter (units N m⁻³).  This is distinct from the Akinci cohesion
    /// which uses `m² * C(r)` — here the density normalisation reduces
    /// clumping artefacts.
    ///
    /// # Arguments
    /// - `positions`  : particle positions.
    /// - `masses`     : particle masses.
    /// - `densities`  : particle densities.
    /// - `neighbors`  : neighbour lists.
    /// - `kappa`      : cohesion strength κ (N m⁻³).
    ///
    /// # Returns
    /// Cohesion force vectors for each particle.
    pub fn compute_cohesion_force(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
        kappa: f64,
    ) -> Vec<Vec3> {
        let n = positions.len();
        let mut forces = vec![Vec3::zeros(); n];
        for i in 0..n {
            let rhoi = densities[i];
            if rhoi < 1e-14 {
                continue;
            }
            for &j in &neighbors[i] {
                let rhoj = densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let c = self.cohesion_kernel(r);
                if c.abs() < 1e-30 {
                    continue;
                }
                let rhat = rij / r;
                let fac = -kappa * masses[i] * masses[j] / (rhoi * rhoj) * c;
                forces[i] += fac * rhat;
            }
        }
        forces
    }
}
/// Morris (2000) surface tension model.
///
/// Uses the pairwise kernel gradient to compute a surface tension force
/// that is proportional to the inter-particle color field difference.
/// This is an alternative to the CSF approach and is better suited for
/// multiphase flows with large density ratios.
#[derive(Debug, Clone)]
pub struct MorrisSurfaceTension {
    /// Surface tension coefficient σ (N/m).
    pub sigma: f64,
    /// Smoothing length h.
    pub smoothing_length: f64,
}
impl MorrisSurfaceTension {
    /// Create a new Morris surface tension model.
    pub fn new(sigma: f64, smoothing_length: f64) -> Self {
        Self {
            sigma,
            smoothing_length,
        }
    }
    /// Compute the Morris surface tension force for each particle.
    ///
    /// The force is:
    /// `f_i = sigma * sum_j (m_j / rho_j) * grad_W(r_ij) * (n_i - n_j)`
    ///
    /// where `n_i` is the (normalized) surface normal.
    pub fn compute_forces(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
    ) -> Vec<Vec3> {
        let kernel = crate::kernel::CubicSplineKernel;
        let h = self.smoothing_length;
        let n = positions.len();
        let mut normals = vec![Vec3::zeros(); n];
        for i in 0..n {
            if densities[i] < 1e-14 {
                continue;
            }
            for &j in &neighbors[i] {
                let rhoj = densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let gw = kernel.grad_w(r, h);
                let rhat = rij / r;
                normals[i] += masses[j] / rhoj * gw * rhat;
            }
        }
        let n_hats: Vec<Vec3> = normals
            .iter()
            .map(|n| {
                let mag = n.norm();
                if mag > 1e-14 { n / mag } else { Vec3::zeros() }
            })
            .collect();
        let mut curvatures = vec![0.0_f64; n];
        for i in 0..n {
            if densities[i] < 1e-14 {
                continue;
            }
            for &j in &neighbors[i] {
                let rhoj = densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let gw = kernel.grad_w(r, h);
                let rhat = rij / r;
                let dn = n_hats[j] - n_hats[i];
                curvatures[i] -= masses[j] / rhoj * dn.dot(&(gw * rhat));
            }
        }
        let mut forces = vec![Vec3::zeros(); n];
        for i in 0..n {
            forces[i] = self.sigma * curvatures[i] * normals[i];
        }
        forces
    }
}
/// Surface tension model combining CSF bulk forces with a contact-line
/// contribution for wetting on solid surfaces.
#[derive(Debug, Clone)]
pub struct SurfaceTension {
    /// Surface tension coefficient σ (N/m).
    pub sigma: f64,
    /// Smoothing length h.
    pub smoothing_length: f64,
    /// Equilibrium (static) contact angle θ_e (radians).
    pub contact_angle: f64,
}
impl SurfaceTension {
    /// Create a new `SurfaceTension` model.
    pub fn new(sigma: f64, smoothing_length: f64, contact_angle: f64) -> Self {
        Self {
            sigma,
            smoothing_length,
            contact_angle,
        }
    }
    /// Compute the contact-line force per unit length on particles near a wall.
    ///
    /// Based on the generalised Navier boundary condition, the contact-line
    /// force density on a fluid particle close to a solid wall is
    ///
    /// ```text
    /// f_cl,i = σ * cos(θ_e) * W(d_i, h) * t̂_wall
    /// ```
    ///
    /// where `d_i` is the distance from particle i to the wall, `W` is the
    /// kernel, and `t̂_wall` is the unit tangent to the wall in the plane of
    /// the contact line (here taken as the projection of the interface normal
    /// onto the wall plane).
    ///
    /// # Arguments
    /// - `positions`      : particle positions.
    /// - `wall_point`     : any point on the wall.
    /// - `wall_normal`    : outward unit normal of the wall (pointing into fluid).
    /// - `interface_normals`: interface normal at each particle (from color gradient).
    ///
    /// # Returns
    /// Contact-line force vector for each particle (zero for far particles).
    pub fn compute_contact_line_force(
        &self,
        positions: &[Vec3],
        wall_point: Vec3,
        wall_normal: Vec3,
        interface_normals: &[Vec3],
    ) -> Vec<Vec3> {
        let n = positions.len();
        let h = self.smoothing_length;
        let cos_theta = self.contact_angle.cos();
        let kernel = CubicSplineKernel;
        let mut forces = vec![Vec3::zeros(); n];
        for i in 0..n {
            let d = (positions[i] - wall_point).dot(&wall_normal);
            if d < 0.0 || d > 2.0 * h {
                continue;
            }
            let w = kernel.w(d, h);
            if w < 1e-30 {
                continue;
            }
            let ni = interface_normals[i];
            let ni_mag = ni.norm();
            if ni_mag < 1e-14 {
                continue;
            }
            let ni_hat = ni / ni_mag;
            let proj = ni_hat.dot(&wall_normal);
            let tangent = ni_hat - proj * wall_normal;
            let t_mag = tangent.norm();
            if t_mag < 1e-14 {
                continue;
            }
            let t_hat = tangent / t_mag;
            forces[i] = self.sigma * cos_theta * w * t_hat;
        }
        forces
    }
}
/// Surface tension model using color function gradient method.
///
/// The CSF (Continuum Surface Force) approach computes surface tension forces
/// by estimating the interface curvature from the color field gradient.
#[derive(Debug, Clone)]
pub struct CsfSurfaceTension {
    /// Surface tension coefficient σ (N/m).
    pub sigma: f64,
    /// Smoothing length h for kernel evaluation.
    pub smoothing_length: f64,
}
impl CsfSurfaceTension {
    /// Create a new CSF surface tension model.
    pub fn new(sigma: f64, smoothing_length: f64) -> Self {
        Self {
            sigma,
            smoothing_length,
        }
    }
    /// Compute color function gradient (surface normal proxy).
    ///
    /// For a single-phase fluid where c = 1 everywhere:
    /// `n_i = grad(c_i) = sum_j (m_j / rho_j) * grad_W(r_ij) * r_hat_ij`
    pub fn compute_color_gradient(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
    ) -> Vec<Vec3> {
        self.compute_color_gradient_with_kernel(
            positions,
            masses,
            densities,
            neighbors,
            &CubicSplineKernel,
        )
    }
    fn compute_color_gradient_with_kernel(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
    ) -> Vec<Vec3> {
        let n = positions.len();
        let h = self.smoothing_length;
        let mut gradients = vec![Vec3::zeros(); n];

        for i in 0..n {
            if densities[i] < 1e-14 {
                continue;
            }
            for &j in &neighbors[i] {
                let rhoj = densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_w = kernel.grad_w(r, h);
                let rhat = rij / r;
                gradients[i] += masses[j] / rhoj * grad_w * rhat;
            }
        }
        gradients
    }
    /// Compute curvature kappa_i = -div(n_hat) where n_hat = n / |n|.
    ///
    /// Uses the SPH divergence estimator:
    /// `kappa_i = -sum_j (m_j / rho_j) * (n_hat_j - n_hat_i) . grad_W(r_ij)`
    pub fn compute_curvature(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        normals: &[Vec3],
        neighbors: &[Vec<usize>],
    ) -> Vec<f64> {
        self.compute_curvature_with_kernel(
            positions,
            masses,
            densities,
            normals,
            neighbors,
            &CubicSplineKernel,
        )
    }
    fn compute_curvature_with_kernel(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        normals: &[Vec3],
        neighbors: &[Vec<usize>],
        kernel: &dyn SphKernel,
    ) -> Vec<f64> {
        let n = positions.len();
        let h = self.smoothing_length;
        let n_hats: Vec<Vec3> = normals
            .iter()
            .map(|v| {
                let mag = v.norm();
                if mag > 1e-14 { v / mag } else { Vec3::zeros() }
            })
            .collect();
        let mut curvatures = vec![0.0_f64; n];

        for i in 0..n {
            if densities[i] < 1e-14 {
                continue;
            }
            for &j in &neighbors[i] {
                let rhoj = densities[j];
                if rhoj < 1e-14 {
                    continue;
                }
                let rij = positions[i] - positions[j];
                let r = rij.norm();
                if r < 1e-14 {
                    continue;
                }
                let grad_w = kernel.grad_w(r, h);
                let rhat = rij / r;
                let dn = n_hats[j] - n_hats[i];
                curvatures[i] -= masses[j] / rhoj * dn.dot(&(grad_w * rhat));
            }
        }
        curvatures
    }
    /// CSF force: f_i = sigma * kappa_i * n_i (force per unit volume).
    ///
    /// Returns the surface tension force vector for each particle.
    pub fn compute_forces(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
    ) -> Vec<Vec3> {
        let normals = self.compute_color_gradient(positions, masses, densities, neighbors);
        let curvatures = self.compute_curvature(positions, masses, densities, &normals, neighbors);
        let n = positions.len();
        let mut forces = vec![Vec3::zeros(); n];

        for i in 0..n {
            forces[i] = self.sigma * curvatures[i] * normals[i];
        }
        forces
    }
    /// CSF force with a threshold on the color gradient magnitude.
    ///
    /// Only particles with `|n_i| > threshold` are considered surface
    /// particles and receive a surface tension force. This reduces noise
    /// from interior particles with near-zero gradients.
    pub fn compute_forces_with_threshold(
        &self,
        positions: &[Vec3],
        masses: &[f64],
        densities: &[f64],
        neighbors: &[Vec<usize>],
        threshold: f64,
    ) -> Vec<Vec3> {
        let normals = self.compute_color_gradient(positions, masses, densities, neighbors);
        let curvatures = self.compute_curvature(positions, masses, densities, &normals, neighbors);
        let n = positions.len();
        let mut forces = vec![Vec3::zeros(); n];

        for i in 0..n {
            if normals[i].norm() > threshold {
                forces[i] = self.sigma * curvatures[i] * normals[i];
            }
        }
        forces
    }
}
