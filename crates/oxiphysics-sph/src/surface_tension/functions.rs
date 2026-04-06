//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::kernel::SphKernel;
use oxiphysics_core::math::Vec3;

use super::types::CsfSurfaceTension;

/// Compute the color field value for each particle.
///
/// `c_i = sum_j (m_j / rho_j) * W(r_ij, h)`
///
/// The color field is 1 inside the fluid and drops to 0 outside.
/// It is useful for identifying the free surface.
#[allow(dead_code, clippy::needless_range_loop)]
pub fn compute_color_field(
    positions: &[Vec3],
    masses: &[f64],
    densities: &[f64],
    neighbors: &[Vec<usize>],
    kernel: &dyn SphKernel,
    h: f64,
) -> Vec<f64> {
    let n = positions.len();
    let mut color = vec![0.0_f64; n];
    for i in 0..n {
        let rhoi = densities[i];
        if rhoi < 1e-14 {
            continue;
        }
        color[i] += masses[i] / rhoi * kernel.w(0.0, h);
        for &j in &neighbors[i] {
            let rhoj = densities[j];
            if rhoj < 1e-14 {
                continue;
            }
            let rij = positions[i] - positions[j];
            let r = rij.norm();
            color[i] += masses[j] / rhoj * kernel.w(r, h);
        }
    }
    color
}
/// Compute the color field gradient for multi-phase flows.
///
/// In multi-phase simulation, each particle has a phase label. The color
/// field gradient points toward the interface between phases.
///
/// * `phase` — phase label for each particle (e.g. 0 or 1)
#[allow(dead_code, clippy::needless_range_loop)]
pub fn compute_multiphase_color_gradient(
    positions: &[Vec3],
    masses: &[f64],
    densities: &[f64],
    phase: &[u32],
    neighbors: &[Vec<usize>],
    kernel: &dyn SphKernel,
    h: f64,
) -> Vec<Vec3> {
    let n = positions.len();
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
            let c_diff = if phase[j] == phase[i] { 0.0 } else { 1.0 };
            let grad_w = kernel.grad_w(r, h);
            let rhat = rij / r;
            gradients[i] += c_diff * masses[j] / rhoj * grad_w * rhat;
        }
    }
    gradients
}
/// Estimate curvature from a pre-computed color gradient using the Laplacian method.
///
/// `kappa_i = -laplacian(c_i) / |grad(c_i)|`
///
/// This is a simple curvature estimator that works well for smooth interfaces.
#[allow(dead_code, clippy::needless_range_loop, clippy::too_many_arguments)]
pub fn estimate_curvature_laplacian(
    positions: &[Vec3],
    masses: &[f64],
    densities: &[f64],
    color: &[f64],
    color_grad: &[Vec3],
    neighbors: &[Vec<usize>],
    kernel: &dyn SphKernel,
    h: f64,
) -> Vec<f64> {
    let n = positions.len();
    let mut curvatures = vec![0.0_f64; n];
    for i in 0..n {
        if densities[i] < 1e-14 {
            continue;
        }
        let grad_mag = color_grad[i].norm();
        if grad_mag < 1e-14 {
            continue;
        }
        let mut laplacian_c = 0.0_f64;
        for &j in &neighbors[i] {
            let rhoj = densities[j];
            if rhoj < 1e-14 {
                continue;
            }
            let rij = positions[i] - positions[j];
            let r = rij.norm();
            let lap = kernel.laplacian_w(r, h);
            laplacian_c += masses[j] / rhoj * (color[j] - color[i]) * lap;
        }
        curvatures[i] = -laplacian_c / grad_mag;
    }
    curvatures
}
/// Compute a wetting angle correction force for particles near a wall.
///
/// Adjusts the color field gradient near wall boundaries to enforce a
/// desired contact angle `theta_w` (in radians). The wall normal points
/// from solid into fluid.
///
/// `n_wall_i = n_i + cos(theta_w) * n_wall + sin(theta_w) * t_wall`
///
/// where `t_wall` is the tangent component of n_i projected onto the wall.
#[allow(dead_code, clippy::needless_range_loop)]
pub fn apply_wetting_angle(
    normals: &mut [Vec3],
    positions: &[Vec3],
    wall_point: &Vec3,
    wall_normal: &Vec3,
    theta_w: f64,
    wall_distance_threshold: f64,
) {
    let n = normals.len();
    let nw = wall_normal / wall_normal.norm().max(1e-14);
    for i in 0..n {
        let d = (positions[i] - wall_point).dot(&nw);
        if d.abs() > wall_distance_threshold || d < 0.0 {
            continue;
        }
        let n_mag = normals[i].norm();
        if n_mag < 1e-14 {
            continue;
        }
        let n_hat = normals[i] / n_mag;
        let t_wall = n_hat - nw * n_hat.dot(&nw);
        let t_mag = t_wall.norm();
        let t_hat = if t_mag > 1e-14 {
            t_wall / t_mag
        } else {
            Vec3::zeros()
        };
        normals[i] = n_mag * (theta_w.cos() * nw + theta_w.sin() * t_hat);
    }
}
/// Identify surface particles based on color gradient magnitude.
///
/// A particle is considered a surface particle if `|grad(c_i)| > threshold`.
/// Returns a boolean mask.
#[allow(dead_code)]
pub fn detect_surface_particles(
    positions: &[Vec3],
    masses: &[f64],
    densities: &[f64],
    neighbors: &[Vec<usize>],
    h: f64,
    threshold: f64,
) -> Vec<bool> {
    let csf = CsfSurfaceTension::new(1.0, h);
    let gradients = csf.compute_color_gradient(positions, masses, densities, neighbors);
    gradients.iter().map(|g| g.norm() > threshold).collect()
}
/// Detect isolated droplets by connected component analysis.
///
/// Particles within `radius` of each other are considered connected.
/// Returns a label for each particle indicating which droplet it belongs to.
/// Label 0 is the largest component; isolated particles get unique labels.
#[allow(dead_code, clippy::needless_range_loop)]
pub fn detect_droplets(positions: &[Vec3], neighbors: &[Vec<usize>]) -> Vec<usize> {
    let n = positions.len();
    if n == 0 {
        return Vec::new();
    }
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank = vec![0_usize; n];
    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    fn union(parent: &mut [usize], rank: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra == rb {
            return;
        }
        if rank[ra] < rank[rb] {
            parent[ra] = rb;
        } else if rank[ra] > rank[rb] {
            parent[rb] = ra;
        } else {
            parent[rb] = ra;
            rank[ra] += 1;
        }
    }
    for i in 0..n {
        for &j in &neighbors[i] {
            union(&mut parent, &mut rank, i, j);
        }
    }
    for i in 0..n {
        find(&mut parent, i);
    }
    let mut component_sizes = std::collections::HashMap::new();
    for i in 0..n {
        let root = find(&mut parent, i);
        *component_sizes.entry(root).or_insert(0_usize) += 1;
    }
    let largest_root = component_sizes
        .iter()
        .max_by_key(|&(_, &count)| count)
        .map(|(&root, _)| root)
        .unwrap_or(0);
    let mut root_to_label = std::collections::HashMap::new();
    root_to_label.insert(largest_root, 0_usize);
    let mut next_label = 1_usize;
    let mut labels = vec![0_usize; n];
    for i in 0..n {
        let root = find(&mut parent, i);
        let label = *root_to_label.entry(root).or_insert_with(|| {
            let l = next_label;
            next_label += 1;
            l
        });
        labels[i] = label;
    }
    labels
}
/// Compute per-droplet center of mass.
///
/// Returns a map from droplet label to (center_of_mass, total_mass).
#[allow(dead_code)]
pub fn droplet_centers_of_mass(
    positions: &[Vec3],
    masses: &[f64],
    labels: &[usize],
) -> std::collections::HashMap<usize, (Vec3, f64)> {
    let mut acc: std::collections::HashMap<usize, (Vec3, f64)> = std::collections::HashMap::new();
    for i in 0..positions.len() {
        let entry = acc.entry(labels[i]).or_insert((Vec3::zeros(), 0.0));
        entry.0 += positions[i] * masses[i];
        entry.1 += masses[i];
    }
    for (_, (com, total_mass)) in acc.iter_mut() {
        if *total_mass > 1e-14 {
            *com /= *total_mass;
        }
    }
    acc
}
/// Compute the effective radius of a droplet assuming spherical shape.
///
/// `R = (3 * V / (4 * pi))^(1/3)` where `V = M / rho`.
#[allow(dead_code)]
pub fn droplet_effective_radius(total_mass: f64, density: f64) -> f64 {
    if density < 1e-14 || total_mass < 1e-14 {
        return 0.0;
    }
    let volume = total_mass / density;
    (3.0 * volume / (4.0 * std::f64::consts::PI)).cbrt()
}
/// Tartakovsky-Meakin (2005) diffuse interface surface tension.
///
/// Uses a pairwise interaction force derived from free energy of mixing:
/// `f_ij = -A_ij * grad_W(r_ij)`
/// where `A_ij` depends on the color field difference between particles.
#[allow(dead_code, clippy::needless_range_loop)]
pub fn tartakovsky_meakin_forces(
    positions: &[Vec3],
    masses: &[f64],
    densities: &[f64],
    phase: &[u32],
    neighbors: &[Vec<usize>],
    h: f64,
    sigma: f64,
) -> Vec<Vec3> {
    let kernel = crate::kernel::CubicSplineKernel;
    let n = positions.len();
    let mut forces = vec![Vec3::zeros(); n];
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
            let c_diff = if phase[i] != phase[j] {
                1.0_f64
            } else {
                0.0_f64
            };
            if c_diff < 1e-14 {
                continue;
            }
            let gw = kernel.grad_w(r, h);
            let rhat = rij / r;
            let a_ij = sigma * c_diff;
            forces[i] -= a_ij * masses[j] / rhoj * gw * rhat;
        }
    }
    forces
}
/// Compute the equilibrium contact angle from Young-Dupré equation.
///
/// `cos(theta) = (gamma_SG - gamma_SL) / gamma_LG`
///
/// where:
/// - `gamma_SG` = solid-gas interfacial energy (J/m^2)
/// - `gamma_SL` = solid-liquid interfacial energy (J/m^2)
/// - `gamma_LG` = liquid-gas surface tension (J/m^2)
///
/// Returns angle in radians, clamped to \[0, π\].
#[allow(dead_code)]
pub fn young_dupree_contact_angle(gamma_sg: f64, gamma_sl: f64, gamma_lg: f64) -> f64 {
    if gamma_lg < 1e-15 {
        return 0.0;
    }
    let cos_theta = (gamma_sg - gamma_sl) / gamma_lg;
    cos_theta.clamp(-1.0, 1.0).acos()
}
/// Compute the spreading coefficient S from interfacial energies.
///
/// `S = gamma_SG - gamma_SL - gamma_LG`
///
/// - S > 0: complete wetting (θ = 0)
/// - S < 0: partial wetting (θ > 0)
#[allow(dead_code)]
pub fn spreading_coefficient(gamma_sg: f64, gamma_sl: f64, gamma_lg: f64) -> f64 {
    gamma_sg - gamma_sl - gamma_lg
}
/// Compute Marangoni forces due to surface tension gradients.
///
/// The Marangoni force arises from surface tension gradients (typically
/// temperature-driven). The force per unit volume is:
/// `f_Marangoni = -dsigma/dT * grad(T)` (approximately)
///
/// Here we compute it via SPH as:
/// `f_i = -dsigma_dT * sum_j (m_j/rho_j) * (T_j - T_i) * grad_W(r_ij)`
#[allow(dead_code, clippy::too_many_arguments, clippy::needless_range_loop)]
pub fn marangoni_forces(
    positions: &[Vec3],
    masses: &[f64],
    densities: &[f64],
    temperatures: &[f64],
    neighbors: &[Vec<usize>],
    h: f64,
    dsigma_dt: f64,
) -> Vec<Vec3> {
    let kernel = crate::kernel::CubicSplineKernel;
    let n = positions.len();
    let mut forces = vec![Vec3::zeros(); n];
    for i in 0..n {
        if densities[i] < 1e-14 {
            continue;
        }
        let mut grad_t = Vec3::zeros();
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
            grad_t += masses[j] / rhoj * (temperatures[j] - temperatures[i]) * gw * rhat;
        }
        forces[i] = -dsigma_dt * grad_t;
    }
    forces
}
/// Weber number: ratio of inertial to surface tension forces.
///
/// `We = ρ U² L / σ`
///
/// - `rho` — fluid density (kg/m³)
/// - `u`   — characteristic velocity (m/s)
/// - `l`   — characteristic length (m)
/// - `sigma` — surface tension coefficient (N/m)
///
/// Returns `f64::INFINITY` when `sigma ≈ 0`.
#[allow(dead_code)]
pub fn weber_number(rho: f64, u: f64, l: f64, sigma: f64) -> f64 {
    if sigma.abs() < 1e-30 {
        return f64::INFINITY;
    }
    rho * u * u * l / sigma
}
/// Capillary number: ratio of viscous to surface tension forces.
///
/// `Ca = μ U / σ`
///
/// - `mu`    — dynamic viscosity (Pa·s)
/// - `u`     — characteristic velocity (m/s)
/// - `sigma` — surface tension coefficient (N/m)
///
/// Returns `f64::INFINITY` when `sigma ≈ 0`.
#[allow(dead_code)]
pub fn capillary_number(mu: f64, u: f64, sigma: f64) -> f64 {
    if sigma.abs() < 1e-30 {
        return f64::INFINITY;
    }
    mu * u / sigma
}
/// Ohnesorge number: ratio of viscous to inertial and surface tension forces.
///
/// `Oh = μ / sqrt(ρ σ L)`
#[allow(dead_code)]
pub fn ohnesorge_number(mu: f64, rho: f64, sigma: f64, l: f64) -> f64 {
    let denom = (rho * sigma * l).sqrt();
    if denom < 1e-30 {
        return f64::INFINITY;
    }
    mu / denom
}
/// Bond number (Eötvös number): ratio of gravitational to surface tension forces.
///
/// `Bo = ρ g L² / σ`
#[allow(dead_code)]
pub fn bond_number(rho: f64, g: f64, l: f64, sigma: f64) -> f64 {
    if sigma.abs() < 1e-30 {
        return f64::INFINITY;
    }
    rho * g * l * l / sigma
}
/// Young–Laplace pressure for a general interface with two principal radii.
///
/// `ΔP = σ (1/R₁ + 1/R₂)`
///
/// For a sphere: `R₁ = R₂ = R` → `ΔP = 2σ/R`.
/// For a cylinder: `R₁ = R`, `R₂ = ∞` → `ΔP = σ/R`.
#[allow(dead_code)]
pub fn young_laplace_pressure(sigma: f64, r1: f64, r2: f64) -> f64 {
    let k1 = if r1.abs() > 1e-30 { 1.0 / r1 } else { 0.0 };
    let k2 = if r2.abs() > 1e-30 { 1.0 / r2 } else { 0.0 };
    sigma * (k1 + k2)
}
/// Determine whether two droplets should coalesce based on their proximity.
///
/// Uses a simplified criterion: droplets coalesce when the gap between them
/// is smaller than a threshold (typically on the order of the molecular scale).
///
/// `gap = center_distance - r1 - r2`
#[allow(dead_code)]
pub fn should_coalesce(_r1: f64, _r2: f64, gap: f64, threshold: f64) -> bool {
    gap < threshold
}
/// Compute the radius of the merged droplet (volume conservation).
///
/// `R_merged = (R1^3 + R2^3)^(1/3)`
#[allow(dead_code)]
pub fn coalescence_merged_radius(r1: f64, r2: f64) -> f64 {
    (r1.powi(3) + r2.powi(3)).cbrt()
}
/// Compute the driving pressure for coalescence from the Laplace pressure.
///
/// The Laplace pressure across a spherical interface: `ΔP = 2*sigma/R`.
/// Two approaching droplets of radii r1, r2 have a net driving force:
/// `F_coalesce = 2*sigma*(1/r1 + 1/r2) * A_neck`
#[allow(dead_code)]
pub fn coalescence_driving_force(sigma: f64, r1: f64, r2: f64, neck_area: f64) -> f64 {
    if r1 < 1e-15 || r2 < 1e-15 {
        return 0.0;
    }
    sigma * (1.0 / r1 + 1.0 / r2) * neck_area
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::CubicSplineKernel;
    use crate::neighbor::SpatialHash;
    use crate::particle::{ParticleSet, SphParticle};
    use crate::surface_tension::types::*;
    use crate::wcsph::compute_density;
    /// Build a uniform 3D grid of particles within a box.
    fn uniform_box(
        nx: usize,
        ny: usize,
        nz: usize,
        spacing: f64,
        mass: f64,
    ) -> (Vec<Vec3>, Vec<f64>, Vec<f64>) {
        let mut ps = ParticleSet::new();
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let pos = Vec3::new(
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    );
                    ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                }
            }
        }
        let h = spacing * 1.3;
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        compute_density(&mut ps, &neighbors, &kernel, h);
        (ps.positions, ps.masses, ps.densities)
    }
    /// Test: uniform fluid → zero surface force (symmetric cancellation).
    ///
    /// For a fully interior particle in an infinite uniform lattice the color
    /// gradient (surface normal) cancels by symmetry, giving zero force.
    #[test]
    fn test_csf_zero_on_uniform_fluid() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(7, 7, 7, spacing, mass);
        let h = spacing * 1.3;
        let csf = CsfSurfaceTension::new(0.072, h);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let forces = csf.compute_forces(&positions, &masses, &densities, &neighbors);
        let center_idx = 3 * 7 * 7 + 3 * 7 + 3;
        let f_mag = forces[center_idx].norm();
        assert!(
            f_mag < 1e-10,
            "Interior particle force should be ~0, got {f_mag:.3e}"
        );
    }
    /// Test: particles on one side of an interface have nonzero surface force.
    #[test]
    fn test_csf_force_magnitude_positive_on_interface() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(7, 7, 4, spacing, mass);
        let h = spacing * 1.3;
        let csf = CsfSurfaceTension::new(0.072, h);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let forces = csf.compute_forces(&positions, &masses, &densities, &neighbors);
        let max_force = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_force > 1e-12,
            "At least one surface particle should have nonzero force, max = {max_force:.3e}"
        );
    }
    /// Test: color gradient at surface points toward the bulk (inward).
    ///
    /// In the SPH color-gradient formulation the gradient
    /// `n_i = sum_j (m_j/rho_j) * grad_W(r_ij) * r_hat_ij`
    /// points **toward** the bulk because grad_W < 0 (kernel is decreasing) and
    /// the neighbours of a top-surface particle are predominantly below (r_hat in +z),
    /// giving a net –z contribution.  This is the standard CSF convention: the
    /// colour-field gradient points inward, and the curvature sign is chosen so
    /// that `sigma * kappa * n` produces a restoring (inward) force.
    #[test]
    fn test_csf_normal_points_outward() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = ParticleSet::new();
        for ix in 0..7_usize {
            for iy in 0..7_usize {
                for iz in 0..4_usize {
                    let pos = Vec3::new(
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    );
                    ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                }
            }
        }
        let h = spacing * 1.3;
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        compute_density(&mut ps, &neighbors, &kernel, h);
        let csf = CsfSurfaceTension::new(0.072, h);
        let grads =
            csf.compute_color_gradient(&ps.positions, &ps.masses, &ps.densities, &neighbors);
        let top_layer: Vec<usize> = (0..ps.len())
            .filter(|&i| (ps.positions[i].z - 3.0 * spacing).abs() < 1e-10)
            .collect();
        assert!(!top_layer.is_empty(), "Should have top-layer particles");
        let nonzero_count = top_layer
            .iter()
            .filter(|&&i| grads[i].z.abs() > 1e-14)
            .count();
        let ratio = nonzero_count as f64 / top_layer.len() as f64;
        assert!(
            ratio > 0.4,
            "Expected most top-layer particles to have a nonzero z gradient, but only {:.0}% did",
            ratio * 100.0
        );
        let inward_count = top_layer.iter().filter(|&&i| grads[i].z < 0.0).count();
        let inward_ratio = inward_count as f64 / top_layer.len() as f64;
        assert!(
            inward_ratio > 0.4,
            "Expected most top-layer colour gradients to point inward (–z), but only {:.0}% did",
            inward_ratio * 100.0
        );
    }
    /// Test: doubling sigma doubles force magnitude.
    #[test]
    fn test_csf_sigma_scales_force() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(7, 7, 4, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let csf1 = CsfSurfaceTension::new(0.072, h);
        let csf2 = CsfSurfaceTension::new(0.144, h);
        let forces1 = csf1.compute_forces(&positions, &masses, &densities, &neighbors);
        let forces2 = csf2.compute_forces(&positions, &masses, &densities, &neighbors);
        let mut checked = 0;
        for i in 0..positions.len() {
            let m1 = forces1[i].norm();
            let m2 = forces2[i].norm();
            if m1 > 1e-14 {
                let ratio = m2 / m1;
                assert!(
                    (ratio - 2.0).abs() < 1e-6,
                    "Particle {i}: expected ratio 2.0 but got {ratio:.6}"
                );
                checked += 1;
            }
        }
        assert!(checked > 0, "No particles with nonzero force found");
    }
    #[test]
    fn test_pairwise_surface_tension_nonzero() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(5, 5, 3, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let pst = PairwiseSurfaceTension::new(0.072, h);
        let forces = pst.compute_forces(&positions, &masses, &densities, &neighbors);
        let max_f = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_f > 1e-14,
            "Pairwise ST should produce nonzero forces: max={max_f:.3e}"
        );
    }
    #[test]
    fn test_pairwise_cohesion_kernel_properties() {
        let h = 0.1;
        let pst = PairwiseSurfaceTension::new(1.0, h);
        let c0 = pst.cohesion_kernel(0.0);
        assert!(c0 < 0.0, "C(0) should be negative: {c0}");
        let ch = pst.cohesion_kernel(h);
        assert!(ch.abs() < 1e-14, "C(h) should be 0: {ch}");
        let c_out = pst.cohesion_kernel(h + 0.01);
        assert!(c_out.abs() < 1e-14, "C(r>h) should be 0: {c_out}");
        let c_mid = pst.cohesion_kernel(0.7 * h);
        assert!(c_mid > 0.0, "C(0.7h) should be positive: {c_mid}");
    }
    #[test]
    fn test_color_field_interior_near_one() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(7, 7, 7, spacing, mass);
        let h = spacing * 1.3;
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let color = compute_color_field(&positions, &masses, &densities, &neighbors, &kernel, h);
        let center_idx = 3 * 7 * 7 + 3 * 7 + 3;
        assert!(
            (color[center_idx] - 1.0).abs() < 0.3,
            "Interior color field should be ~1.0, got {}",
            color[center_idx]
        );
    }
    #[test]
    fn test_color_field_surface_less_than_interior() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(7, 7, 4, spacing, mass);
        let h = spacing * 1.3;
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let color = compute_color_field(&positions, &masses, &densities, &neighbors, &kernel, h);
        let center_idx = 3 * 7 * 4 + 3 * 4 + 1;
        let corner_idx = 0;
        assert!(
            color[corner_idx] < color[center_idx],
            "Corner color {} should be less than center color {}",
            color[corner_idx],
            color[center_idx]
        );
    }
    #[test]
    fn test_multiphase_gradient_at_interface() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = ParticleSet::new();
        let mut phases = Vec::new();
        for ix in 0..5_usize {
            for iy in 0..5_usize {
                for iz in 0..5_usize {
                    let pos = Vec3::new(
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    );
                    ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                    phases.push(if iz < 2 { 0_u32 } else { 1_u32 });
                }
            }
        }
        let h = spacing * 1.3;
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        compute_density(&mut ps, &neighbors, &kernel, h);
        let grad = compute_multiphase_color_gradient(
            &ps.positions,
            &ps.masses,
            &ps.densities,
            &phases,
            &neighbors,
            &kernel,
            h,
        );
        let interface_particles: Vec<usize> = (0..ps.len())
            .filter(|&i| {
                let iz = (ps.positions[i].z / spacing + 0.5) as usize;
                iz == 1 || iz == 2
            })
            .collect();
        let has_gradient = interface_particles.iter().any(|&i| grad[i].norm() > 1e-14);
        assert!(
            has_gradient,
            "Interface particles should have nonzero color gradient"
        );
    }
    #[test]
    fn test_multiphase_gradient_zero_in_bulk() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = ParticleSet::new();
        let mut phases = Vec::new();
        for ix in 0..5_usize {
            for iy in 0..5_usize {
                for iz in 0..5_usize {
                    let pos = Vec3::new(
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    );
                    ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                    phases.push(0_u32);
                }
            }
        }
        let h = spacing * 1.3;
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        compute_density(&mut ps, &neighbors, &kernel, h);
        let grad = compute_multiphase_color_gradient(
            &ps.positions,
            &ps.masses,
            &ps.densities,
            &phases,
            &neighbors,
            &kernel,
            h,
        );
        let max_grad = grad.iter().map(|g| g.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_grad < 1e-14,
            "Uniform phase should have zero gradient: {max_grad}"
        );
    }
    #[test]
    fn test_curvature_laplacian_estimator() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(7, 7, 4, spacing, mass);
        let h = spacing * 1.3;
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let csf = CsfSurfaceTension::new(0.072, h);
        let color_grad = csf.compute_color_gradient(&positions, &masses, &densities, &neighbors);
        let color = compute_color_field(&positions, &masses, &densities, &neighbors, &kernel, h);
        let curvature = estimate_curvature_laplacian(
            &positions,
            &masses,
            &densities,
            &color,
            &color_grad,
            &neighbors,
            &kernel,
            h,
        );
        let finite_count = curvature.iter().filter(|&&k| k.is_finite()).count();
        assert_eq!(
            finite_count,
            curvature.len(),
            "All curvatures should be finite"
        );
    }
    #[test]
    fn test_detect_surface_particles() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(7, 7, 4, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let surface = detect_surface_particles(&positions, &masses, &densities, &neighbors, h, 1.0);
        let surface_count = surface.iter().filter(|&&s| s).count();
        let total = positions.len();
        assert!(surface_count > 0, "Should detect some surface particles");
        assert!(
            surface_count < total,
            "Not all particles should be surface particles"
        );
    }
    #[test]
    fn test_detect_droplets_single_blob() {
        let spacing = 0.1_f64;
        let (positions, _masses, _densities) = uniform_box(4, 4, 4, spacing, 1.0);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let labels = detect_droplets(&positions, &neighbors);
        assert!(
            labels.iter().all(|&l| l == 0),
            "Single blob should have all label 0"
        );
    }
    #[test]
    fn test_detect_droplets_two_separated() {
        let spacing = 0.1_f64;
        let mut positions = Vec::new();
        for ix in 0..3_usize {
            for iy in 0..3_usize {
                for iz in 0..3_usize {
                    positions.push(Vec3::new(
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    ));
                }
            }
        }
        let blob1_size = positions.len();
        for ix in 0..3_usize {
            for iy in 0..3_usize {
                for iz in 0..3_usize {
                    positions.push(Vec3::new(
                        10.0 + ix as f64 * spacing,
                        10.0 + iy as f64 * spacing,
                        10.0 + iz as f64 * spacing,
                    ));
                }
            }
        }
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let labels = detect_droplets(&positions, &neighbors);
        let unique: std::collections::HashSet<usize> = labels.iter().copied().collect();
        assert_eq!(
            unique.len(),
            2,
            "Should detect 2 droplets, got {} ({:?})",
            unique.len(),
            unique
        );
        let _blob1_label = labels[0];
        let blob2_label = labels[blob1_size];
        assert_ne!(
            labels[0], blob2_label,
            "Two blobs should have different labels"
        );
    }
    #[test]
    fn test_droplet_centers_of_mass() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(10.0, 0.0, 0.0),
        ];
        let masses = vec![1.0, 1.0, 2.0];
        let labels = vec![0, 0, 1];
        let coms = droplet_centers_of_mass(&positions, &masses, &labels);
        let (com0, m0) = coms[&0];
        assert!((m0 - 2.0).abs() < 1e-12, "Total mass label 0: {m0}");
        assert!(
            (com0.x - 0.5).abs() < 1e-12,
            "CoM x for label 0: {}",
            com0.x
        );
        let (com1, m1) = coms[&1];
        assert!((m1 - 2.0).abs() < 1e-12, "Total mass label 1: {m1}");
        assert!(
            (com1.x - 10.0).abs() < 1e-12,
            "CoM x for label 1: {}",
            com1.x
        );
    }
    #[test]
    fn test_droplet_effective_radius() {
        let mass = 1000.0 * (4.0 / 3.0) * std::f64::consts::PI;
        let r = droplet_effective_radius(mass, 1000.0);
        assert!(
            (r - 1.0).abs() < 1e-10,
            "Effective radius should be ~1, got {r}"
        );
    }
    #[test]
    fn test_droplet_effective_radius_zero_density() {
        let r = droplet_effective_radius(1.0, 0.0);
        assert!(r.abs() < 1e-14, "Zero density should give zero radius");
    }
    #[test]
    fn test_morris_surface_tension_nonzero() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(5, 5, 3, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let morris = MorrisSurfaceTension::new(0.072, h);
        let forces = morris.compute_forces(&positions, &masses, &densities, &neighbors);
        let max_f = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_f > 1e-12,
            "Morris ST should produce nonzero forces: {max_f}"
        );
    }
    #[test]
    fn test_morris_sigma_scales_force() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(7, 7, 4, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let m1 = MorrisSurfaceTension::new(0.072, h);
        let m2 = MorrisSurfaceTension::new(0.144, h);
        let f1 = m1.compute_forces(&positions, &masses, &densities, &neighbors);
        let f2 = m2.compute_forces(&positions, &masses, &densities, &neighbors);
        let mut checked = 0;
        for i in 0..positions.len() {
            let mag1 = f1[i].norm();
            let mag2 = f2[i].norm();
            if mag1 > 1e-14 {
                let ratio = mag2 / mag1;
                assert!((ratio - 2.0).abs() < 1e-6, "ratio = {ratio}");
                checked += 1;
            }
        }
        assert!(checked > 0);
    }
    #[test]
    fn test_marangoni_zero_on_uniform_temp() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(5, 5, 5, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let temperatures = vec![300.0; positions.len()];
        let forces = marangoni_forces(
            &positions,
            &masses,
            &densities,
            &temperatures,
            &neighbors,
            h,
            -0.0001,
        );
        let max_f = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_f < 1e-10,
            "Uniform temp → zero Marangoni force: {max_f}"
        );
    }
    #[test]
    fn test_marangoni_nonzero_with_gradient() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let mut ps = ParticleSet::new();
        let mut temperatures = Vec::new();
        for ix in 0..5_usize {
            for iy in 0..5_usize {
                for iz in 0..3_usize {
                    let pos = Vec3::new(
                        ix as f64 * spacing,
                        iy as f64 * spacing,
                        iz as f64 * spacing,
                    );
                    ps.add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                    temperatures.push(300.0 + ix as f64 * 10.0);
                }
            }
        }
        let h = spacing * 1.3;
        let kernel = CubicSplineKernel;
        let neighbors = SpatialHash::find_all_neighbors(&ps.positions, 2.0 * h);
        compute_density(&mut ps, &neighbors, &kernel, h);
        let forces = marangoni_forces(
            &ps.positions,
            &ps.masses,
            &ps.densities,
            &temperatures,
            &neighbors,
            h,
            -0.0001,
        );
        let max_f = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_f > 1e-14,
            "Temperature gradient → nonzero Marangoni: {max_f}"
        );
    }
    #[test]
    fn test_young_dupree_contact_angle() {
        let theta = young_dupree_contact_angle(0.5, 0.3, 0.072);
        assert!(
            (0.0..=std::f64::consts::PI).contains(&theta),
            "theta = {theta}"
        );
        assert!(theta < 1e-6, "Complete wetting expected, theta = {theta}");
    }
    #[test]
    fn test_young_dupree_partial_wetting() {
        let theta = young_dupree_contact_angle(0.03, 0.07, 0.072);
        assert!(
            theta > std::f64::consts::FRAC_PI_2,
            "Hydrophobic: theta = {theta}"
        );
    }
    #[test]
    fn test_curvature_coalescence_criterion() {
        let r1 = 1e-3;
        let r2 = 1e-3;
        let gap = 1e-5;
        assert!(
            should_coalesce(r1, r2, gap, 1e-4),
            "Close droplets should coalesce"
        );
        assert!(
            !should_coalesce(r1, r2, 1.0, 1e-4),
            "Far droplets should not coalesce"
        );
    }
    #[test]
    fn test_merged_droplet_radius() {
        let r1 = 1e-3;
        let r2 = 1e-3;
        let r_merged = coalescence_merged_radius(r1, r2);
        let expected = (r1.powi(3) + r2.powi(3)).cbrt();
        assert!(
            (r_merged - expected).abs() < 1e-20,
            "Merged radius = {r_merged}"
        );
    }
    #[test]
    fn test_tartakovsky_meakin_diffuse() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(5, 5, 3, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let phases = vec![0_u32; positions.len()];
        let forces = tartakovsky_meakin_forces(
            &positions, &masses, &densities, &phases, &neighbors, h, 0.072,
        );
        let max_f = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_f < 1e-10,
            "Uniform phase → zero diffuse interface force: {max_f}"
        );
    }
    #[test]
    fn test_wetting_angle_modifies_normals() {
        let mut normals = vec![Vec3::new(0.0, 0.0, 1.0); 3];
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.05),
            Vec3::new(0.0, 0.0, 0.05),
            Vec3::new(0.0, 0.0, 5.0),
        ];
        let wall_point = Vec3::new(0.0, 0.0, 0.0);
        let wall_normal = Vec3::new(0.0, 0.0, 1.0);
        let theta_w = std::f64::consts::FRAC_PI_4;
        let orig = normals.clone();
        apply_wetting_angle(
            &mut normals,
            &positions,
            &wall_point,
            &wall_normal,
            theta_w,
            0.1,
        );
        let diff0 = (normals[0] - orig[0]).norm();
        assert!(
            diff0 > 1e-10,
            "Near-wall normal should be modified: diff={diff0}"
        );
        let diff2 = (normals[2] - orig[2]).norm();
        assert!(
            diff2 < 1e-10,
            "Far particle normal should be unchanged: diff={diff2}"
        );
    }
    #[test]
    fn test_weber_number_basic() {
        let we = weber_number(1000.0, 1.0, 1e-3, 0.072);
        assert!((we - 1000.0 * 1.0_f64.powi(2) * 1e-3 / 0.072).abs() < 1e-10);
    }
    #[test]
    fn test_capillary_number_basic() {
        let ca = capillary_number(1e-3, 0.1, 0.072);
        assert!((ca - 1e-3 * 0.1 / 0.072).abs() < 1e-20);
    }
    #[test]
    fn test_weber_number_zero_sigma() {
        let we = weber_number(1000.0, 1.0, 1e-3, 0.0);
        assert!(
            we.is_infinite(),
            "zero surface tension → infinite Weber number"
        );
    }
    #[test]
    fn test_dynamic_contact_angle_zero_velocity() {
        let static_deg = 70.0_f64;
        let theta_static = static_deg.to_radians();
        let dca = DynamicContactAngle::new(theta_static, 0.01, 0.072);
        let theta_d = dca.angle(0.0);
        assert!(
            (theta_d - theta_static).abs() < 1e-12,
            "zero velocity → static angle, got {theta_d}"
        );
    }
    #[test]
    fn test_dynamic_contact_angle_advancing_larger() {
        let dca = DynamicContactAngle::new(70_f64.to_radians(), 0.01, 0.072);
        let theta_adv = dca.angle(0.001);
        let theta_static = 70_f64.to_radians();
        assert!(
            theta_adv >= theta_static,
            "advancing angle should be ≥ static: {theta_adv:.4} vs {theta_static:.4}"
        );
    }
    #[test]
    fn test_dynamic_contact_angle_receding_smaller() {
        let dca = DynamicContactAngle::new(70_f64.to_radians(), 0.01, 0.072);
        let theta_rec = dca.angle(-0.001);
        let theta_static = 70_f64.to_radians();
        assert!(
            theta_rec <= theta_static,
            "receding angle should be ≤ static: {theta_rec:.4} vs {theta_static:.4}"
        );
    }
    #[test]
    fn test_sessile_droplet_contact_width_is_positive() {
        let drop = SessileDroplet::new(1e-3, 1000.0, 0.072, 70_f64.to_radians());
        let w = drop.contact_width();
        assert!(w > 0.0, "contact width must be positive, got {w}");
    }
    #[test]
    fn test_sessile_droplet_height_positive() {
        let drop = SessileDroplet::new(1e-3, 1000.0, 0.072, 70_f64.to_radians());
        let h = drop.height();
        assert!(h > 0.0, "droplet height must be positive, got {h}");
    }
    #[test]
    fn test_sessile_droplet_completely_wetting_flat() {
        let drop_flat = SessileDroplet::new(1e-3, 1000.0, 0.072, 0.01_f64.to_radians());
        let drop_normal = SessileDroplet::new(1e-3, 1000.0, 0.072, 90_f64.to_radians());
        assert!(
            drop_flat.height() < drop_normal.height(),
            "flat droplet should be shorter than 90° droplet"
        );
    }
    #[test]
    fn test_sessile_droplet_laplace_pressure() {
        let drop = SessileDroplet::new(1e-3, 1000.0, 0.072, 90_f64.to_radians());
        let dp = drop.laplace_pressure();
        assert!(dp > 0.0, "Laplace pressure must be positive");
    }
    #[test]
    fn test_sessile_droplet_bond_number() {
        let drop = SessileDroplet::new(1e-3, 1000.0, 0.072, 90_f64.to_radians());
        let bo = drop.bond_number(9.81);
        assert!(bo > 0.0, "Bond number must be positive");
    }
    #[test]
    fn test_young_laplace_sphere() {
        let dp = young_laplace_pressure(0.072, 1e-3, 1e-3);
        let expected = 2.0 * 0.072 / 1e-3;
        assert!(
            (dp - expected).abs() < 1e-10,
            "ΔP = {dp}, expected {expected}"
        );
    }
    #[test]
    fn test_young_laplace_cylinder() {
        let dp = young_laplace_pressure(0.072, 1e-3, 1e30);
        let expected = 0.072 / 1e-3;
        assert!(
            (dp - expected).abs() < 1e-5,
            "cylinder ΔP = {dp}, expected ~{expected}"
        );
    }
    #[test]
    fn test_ohnesorge_number_basic() {
        let oh = ohnesorge_number(1e-3, 1000.0, 0.072, 1e-3);
        let expected = 1e-3 / (1000.0_f64 * 0.072 * 1e-3).sqrt();
        assert!(
            (oh - expected).abs() < 1e-15,
            "Oh = {oh}, expected {expected}"
        );
    }
    #[test]
    fn test_bond_number_basic() {
        let bo = bond_number(1000.0, 9.81, 1e-3, 0.072);
        let expected = 1000.0 * 9.81 * 1e-3 * 1e-3 / 0.072;
        assert!(
            (bo - expected).abs() < 1e-12,
            "Bo = {bo}, expected {expected}"
        );
    }
    #[test]
    fn test_csf_model_correction_uniform_color_zero_grad() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(5, 5, 5, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let color = vec![1.0_f64; positions.len()];
        let model = CsfModel::new(0.072, h);
        let corrected = model
            .compute_color_gradient_correction(&positions, &masses, &densities, &color, &neighbors);
        let max_mag = corrected.iter().map(|g| g.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_mag < 1e-10,
            "Uniform color → zero corrected gradient, max_mag={max_mag}"
        );
    }
    #[test]
    fn test_csf_model_correction_returns_correct_count() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(4, 4, 4, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let color: Vec<f64> = positions.iter().map(|p| p.x).collect();
        let model = CsfModel::new(0.072, h);
        let corrected = model
            .compute_color_gradient_correction(&positions, &masses, &densities, &color, &neighbors);
        assert_eq!(
            corrected.len(),
            positions.len(),
            "Output length should match particle count"
        );
    }
    #[test]
    fn test_csf_model_correction_finite_values() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(6, 6, 4, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let color: Vec<f64> = positions.iter().map(|p| 0.5 * p.z).collect();
        let model = CsfModel::new(0.072, h);
        let corrected = model
            .compute_color_gradient_correction(&positions, &masses, &densities, &color, &neighbors);
        let all_finite = corrected
            .iter()
            .all(|g| g.x.is_finite() && g.y.is_finite() && g.z.is_finite());
        assert!(all_finite, "All corrected gradients should be finite");
    }
    #[test]
    fn test_cohesion_force_zero_kappa() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(4, 4, 4, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let pst = PairwiseSurfaceTension::new(0.072, h);
        let forces = pst.compute_cohesion_force(&positions, &masses, &densities, &neighbors, 0.0);
        let max_f = forces.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(max_f < 1e-20, "κ=0 → zero cohesion force, max_f={max_f}");
    }
    #[test]
    fn test_cohesion_force_scales_with_kappa() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(5, 5, 3, spacing, mass);
        let h = spacing * 1.3;
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let pst = PairwiseSurfaceTension::new(0.072, h);
        let f1 = pst.compute_cohesion_force(&positions, &masses, &densities, &neighbors, 1.0);
        let f2 = pst.compute_cohesion_force(&positions, &masses, &densities, &neighbors, 2.0);
        let mut checked = 0;
        for i in 0..positions.len() {
            let m1 = f1[i].norm();
            let m2 = f2[i].norm();
            if m1 > 1e-20 {
                let ratio = m2 / m1;
                assert!(
                    (ratio - 2.0).abs() < 1e-8,
                    "Force should scale 2x, got ratio={ratio}"
                );
                checked += 1;
            }
        }
        assert!(
            checked > 0,
            "Should have at least one nonzero cohesion force"
        );
    }
    #[test]
    fn test_cohesion_force_nonzero_for_nearby_particles() {
        let h = 0.3_f64;
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.1, 0.0, 0.0)];
        let masses = vec![1e-3; 2];
        let densities = vec![1000.0; 2];
        let neighbors = vec![vec![1], vec![0]];
        let pst = PairwiseSurfaceTension::new(0.072, h);
        let forces = pst.compute_cohesion_force(&positions, &masses, &densities, &neighbors, 1.0);
        let mag = forces[0].norm();
        assert!(
            mag > 1e-20,
            "Nearby particles → nonzero cohesion force, mag={mag}"
        );
    }
    #[test]
    fn test_contact_line_force_far_from_wall_is_zero() {
        let h = 0.1_f64;
        let model = SurfaceTension::new(0.072, h, std::f64::consts::FRAC_PI_4);
        let positions = vec![Vec3::new(0.0, 0.0, 5.0)];
        let wall_point = Vec3::zeros();
        let wall_normal = Vec3::new(0.0, 0.0, 1.0);
        let normals = vec![Vec3::new(1.0, 0.0, 0.0)];
        let f = model.compute_contact_line_force(&positions, wall_point, wall_normal, &normals);
        assert!(
            f[0].norm() < 1e-20,
            "Far particle → zero contact-line force"
        );
    }
    #[test]
    fn test_contact_line_force_nonzero_near_wall() {
        let h = 0.3_f64;
        let model = SurfaceTension::new(0.072, h, std::f64::consts::FRAC_PI_4);
        let positions = vec![Vec3::new(0.0, 0.0, 0.05)];
        let wall_point = Vec3::zeros();
        let wall_normal = Vec3::new(0.0, 0.0, 1.0);
        let normals = vec![Vec3::new(1.0, 0.0, 1.0).normalize()];
        let f = model.compute_contact_line_force(&positions, wall_point, wall_normal, &normals);
        let mag = f[0].norm();
        assert!(
            mag > 1e-10,
            "Near-wall particle → nonzero contact-line force, mag={mag}"
        );
    }
    #[test]
    fn test_contact_line_force_zero_interface_normal_gives_zero() {
        let h = 0.3_f64;
        let model = SurfaceTension::new(0.072, h, std::f64::consts::FRAC_PI_4);
        let positions = vec![Vec3::new(0.0, 0.0, 0.05)];
        let wall_point = Vec3::zeros();
        let wall_normal = Vec3::new(0.0, 0.0, 1.0);
        let normals = vec![Vec3::zeros()];
        let f = model.compute_contact_line_force(&positions, wall_point, wall_normal, &normals);
        assert!(
            f[0].norm() < 1e-20,
            "Zero interface normal → zero contact-line force"
        );
    }
    #[test]
    fn test_contact_line_force_scales_with_sigma() {
        let h = 0.3_f64;
        let model1 = SurfaceTension::new(0.072, h, std::f64::consts::FRAC_PI_4);
        let model2 = SurfaceTension::new(0.144, h, std::f64::consts::FRAC_PI_4);
        let positions = vec![Vec3::new(0.0, 0.0, 0.05)];
        let wall_point = Vec3::zeros();
        let wall_normal = Vec3::new(0.0, 0.0, 1.0);
        let normals = vec![Vec3::new(1.0, 0.0, 1.0).normalize()];
        let f1 = model1.compute_contact_line_force(&positions, wall_point, wall_normal, &normals);
        let f2 = model2.compute_contact_line_force(&positions, wall_point, wall_normal, &normals);
        let m1 = f1[0].norm();
        let m2 = f2[0].norm();
        if m1 > 1e-20 {
            let ratio = m2 / m1;
            assert!(
                (ratio - 2.0).abs() < 1e-10,
                "Contact-line force should scale 2x with σ, got ratio={ratio}"
            );
        }
    }
    #[test]
    fn test_csf_forces_with_threshold() {
        let spacing = 0.1_f64;
        let mass = 1000.0 * spacing.powi(3);
        let (positions, masses, densities) = uniform_box(7, 7, 4, spacing, mass);
        let h = spacing * 1.3;
        let csf = CsfSurfaceTension::new(0.072, h);
        let neighbors = SpatialHash::find_all_neighbors(&positions, 2.0 * h);
        let forces_high =
            csf.compute_forces_with_threshold(&positions, &masses, &densities, &neighbors, 1e6);
        let max_high = forces_high.iter().map(|f| f.norm()).fold(0.0_f64, f64::max);
        assert!(
            max_high < 1e-14,
            "Very high threshold should suppress all forces: {max_high}"
        );
        let forces_low =
            csf.compute_forces_with_threshold(&positions, &masses, &densities, &neighbors, 1e-30);
        let forces_no = csf.compute_forces(&positions, &masses, &densities, &neighbors);
        for i in 0..positions.len() {
            let diff = (forces_low[i] - forces_no[i]).norm();
            assert!(
                diff < 1e-14,
                "Low threshold should match no-threshold: diff={diff}"
            );
        }
    }
}
