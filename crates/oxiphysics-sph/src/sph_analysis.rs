// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH post-processing and analysis tools.
//!
//! Implements:
//!
//! - [`SphFieldInterp`]: Interpolate scalar/vector fields to query points
//! - [`SphStatistics`]: Kinetic/potential energy, momentum diagnostics
//! - [`VelocityGradient`]: SPH velocity gradient, strain rate, vorticity
//! - [`SphTurbulence`]: TKE, dissipation, Reynolds stress, structure function
//! - [`SphFreeSurface`]: Free surface detection via density criterion
//! - [`SphParticleSorting`]: Z-curve Morton code sort for cache-friendly search
//! - [`SphDiagnostics`]: Max velocity, density error, disorder parameter
//! - [`SphVtkExporter`]: Export to VTK PolyData XML format
//! - [`SphForceBalance`]: Pressure/viscous/body/surface tension force diagnostics
//! - [`SphConvergence`]: L2 error vs analytical solution, convergence rate

// ============================================================================
// Field Interpolation
// ============================================================================

/// SPH field interpolation from particles to arbitrary query points.
///
/// Uses SPH kernel-weighted interpolation:
/// φ(x) = Σⱼ (mⱼ/ρⱼ) * φⱼ * W(|x - xⱼ|, h)
pub struct SphFieldInterp {
    /// Smoothing length h (m).
    pub h: f64,
    /// Particle masses (kg).
    pub masses: Vec<f64>,
    /// Particle densities (kg/m³).
    pub densities: Vec<f64>,
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
}

impl SphFieldInterp {
    /// Create a new SPH field interpolator.
    pub fn new(h: f64, masses: Vec<f64>, densities: Vec<f64>, positions: Vec<[f64; 3]>) -> Self {
        Self {
            h,
            masses,
            densities,
            positions,
        }
    }

    /// Cubic spline kernel W(r, h).
    ///
    /// Returns normalized cubic spline value.
    pub fn cubic_spline_kernel(r: f64, h: f64) -> f64 {
        let q = r / h;
        let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
        if q <= 1.0 {
            sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        } else if q <= 2.0 {
            sigma * 0.25 * (2.0 - q).powi(3)
        } else {
            0.0
        }
    }

    /// Interpolate a scalar field at query point x.
    ///
    /// φ(x) = Σⱼ (mⱼ/ρⱼ) * φⱼ * W(|x - xⱼ|, h)
    pub fn interpolate_scalar(&self, x: [f64; 3], field: &[f64]) -> f64 {
        sph_interpolate_scalar(
            x,
            &self.positions,
            field,
            &self.masses,
            &self.densities,
            self.h,
        )
    }

    /// Interpolate a 3D vector field at query point x.
    pub fn interpolate_vector(&self, x: [f64; 3], field: &[[f64; 3]]) -> [f64; 3] {
        let fx: Vec<f64> = field.iter().map(|v| v[0]).collect();
        let fy: Vec<f64> = field.iter().map(|v| v[1]).collect();
        let fz: Vec<f64> = field.iter().map(|v| v[2]).collect();
        [
            self.interpolate_scalar(x, &fx),
            self.interpolate_scalar(x, &fy),
            self.interpolate_scalar(x, &fz),
        ]
    }

    /// Kernel-weighted number density at x.
    pub fn number_density(&self, x: [f64; 3]) -> f64 {
        self.positions
            .iter()
            .map(|&xj| {
                let r = dist3(x, xj);
                Self::cubic_spline_kernel(r, self.h)
            })
            .sum()
    }
}

// ============================================================================
// Statistics
// ============================================================================

/// SPH statistical diagnostics: energy, momentum, temperature.
#[derive(Debug, Clone, Default)]
pub struct SphStatistics {
    /// Total kinetic energy (J).
    pub kinetic_energy: f64,
    /// Total potential energy (J).
    pub potential_energy: f64,
    /// Total momentum vector (kg·m/s).
    pub momentum: [f64; 3],
    /// Angular momentum vector (kg·m²/s).
    pub angular_momentum: [f64; 3],
    /// Average temperature (K) if applicable.
    pub avg_temperature: f64,
}

impl SphStatistics {
    /// Compute statistics from particle state.
    ///
    /// - `masses`: particle masses (kg)
    /// - `velocities`: particle velocities (m/s)
    /// - `positions`: particle positions (m)
    /// - `gravity`: gravitational acceleration (m/s²)
    pub fn compute(
        masses: &[f64],
        velocities: &[[f64; 3]],
        positions: &[[f64; 3]],
        gravity: [f64; 3],
    ) -> Self {
        let mut ke = 0.0;
        let mut pe = 0.0;
        let mut mom = [0.0_f64; 3];
        let mut ang = [0.0_f64; 3];
        let n = masses.len().min(velocities.len()).min(positions.len());
        for i in 0..n {
            let m = masses[i];
            let v = velocities[i];
            let x = positions[i];
            let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            ke += 0.5 * m * v2;
            pe += -m * (gravity[0] * x[0] + gravity[1] * x[1] + gravity[2] * x[2]);
            mom[0] += m * v[0];
            mom[1] += m * v[1];
            mom[2] += m * v[2];
            // L = r × (m*v)
            ang[0] += x[1] * m * v[2] - x[2] * m * v[1];
            ang[1] += x[2] * m * v[0] - x[0] * m * v[2];
            ang[2] += x[0] * m * v[1] - x[1] * m * v[0];
        }
        Self {
            kinetic_energy: ke,
            potential_energy: pe,
            momentum: mom,
            angular_momentum: ang,
            avg_temperature: 0.0,
        }
    }

    /// Total mechanical energy E = KE + PE.
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy + self.potential_energy
    }

    /// Momentum magnitude.
    pub fn momentum_magnitude(&self) -> f64 {
        let m = self.momentum;
        (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt()
    }
}

// ============================================================================
// Velocity Gradient
// ============================================================================

/// SPH velocity gradient tensor computation.
///
/// ∂uₐ/∂xᵦ ≈ Σⱼ (mⱼ/ρⱼ) * (u^α_j - u^α_i) * ∂W/∂x^β_{ij}
pub struct VelocityGradient {
    /// Smoothing length h (m).
    pub h: f64,
}

impl VelocityGradient {
    /// Create a velocity gradient calculator.
    pub fn new(h: f64) -> Self {
        Self { h }
    }

    /// Cubic spline kernel gradient magnitude: dW/dr.
    pub fn kernel_gradient_magnitude(r: f64, h: f64) -> f64 {
        let q = r / h;
        let sigma = 1.0 / (std::f64::consts::PI * h * h * h * h);
        if q <= 1.0 {
            sigma * (-3.0 * q + 2.25 * q * q)
        } else if q <= 2.0 {
            -sigma * 0.75 * (2.0 - q).powi(2)
        } else {
            0.0
        }
    }

    /// Kernel gradient vector ∇W(rᵢⱼ, h) = dW/dr * r̂.
    pub fn kernel_gradient_vector(r_ij: [f64; 3], h: f64) -> [f64; 3] {
        let r = dist_mag(r_ij);
        if r < 1e-14 {
            return [0.0; 3];
        }
        let dw_dr = Self::kernel_gradient_magnitude(r, h);
        [
            dw_dr * r_ij[0] / r,
            dw_dr * r_ij[1] / r,
            dw_dr * r_ij[2] / r,
        ]
    }

    /// Velocity gradient tensor (3×3 row-major) at particle i.
    ///
    /// L_αβ = ∂uₐ/∂xᵦ
    pub fn gradient_tensor(
        &self,
        u_i: [f64; 3],
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        masses: &[f64],
        densities: &[f64],
        xi: [f64; 3],
    ) -> [f64; 9] {
        let mut l = [0.0_f64; 9];
        for (j, (&xj, (&vj, (&mj, &rhoj)))) in positions
            .iter()
            .zip(velocities.iter().zip(masses.iter().zip(densities.iter())))
            .enumerate()
        {
            let _ = j;
            let r_ij = [xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]];
            let grad_w = Self::kernel_gradient_vector(r_ij, self.h);
            let dv = [vj[0] - u_i[0], vj[1] - u_i[1], vj[2] - u_i[2]];
            let c = mj / rhoj;
            // L_αβ += (mⱼ/ρⱼ) * (vⱼ^α - vᵢ^α) * ∂W/∂xᵦ
            for alpha in 0..3 {
                for beta in 0..3 {
                    l[alpha * 3 + beta] += c * dv[alpha] * grad_w[beta];
                }
            }
        }
        l
    }

    /// Strain rate tensor S = (L + Lᵀ) / 2 (3×3 row-major).
    pub fn strain_rate(l: &[f64; 9]) -> [f64; 9] {
        let mut s = [0.0_f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                s[i * 3 + j] = 0.5 * (l[i * 3 + j] + l[j * 3 + i]);
            }
        }
        s
    }

    /// Rotation rate tensor Ω = (L - Lᵀ) / 2 (3×3 row-major).
    pub fn rotation_rate(l: &[f64; 9]) -> [f64; 9] {
        let mut omega = [0.0_f64; 9];
        for i in 0..3 {
            for j in 0..3 {
                omega[i * 3 + j] = 0.5 * (l[i * 3 + j] - l[j * 3 + i]);
            }
        }
        omega
    }

    /// Vorticity vector ω = \[Ω23 - Ω32, Ω31 - Ω13, Ω12 - Ω21\].
    ///
    /// ω = ∇ × u
    pub fn vorticity(l: &[f64; 9]) -> [f64; 3] {
        [l[2 * 3 + 1] - l[3 + 2], l[2] - l[2 * 3], l[3] - l[1]]
    }

    /// Second invariant of strain rate tensor: Q = -0.5 * S_ij * S_ij.
    pub fn second_invariant_s(s: &[f64; 9]) -> f64 {
        let mut q = 0.0;
        for &v in s.iter() {
            q += v * v;
        }
        -0.5 * q
    }
}

// ============================================================================
// Turbulence
// ============================================================================

/// SPH turbulence diagnostics.
///
/// Computes turbulent kinetic energy (TKE), dissipation rate,
/// Reynolds stress tensor, and second-order structure function.
pub struct SphTurbulence {
    /// Mean velocity field.
    pub mean_velocity: Vec<[f64; 3]>,
}

impl SphTurbulence {
    /// Create turbulence analyzer with given mean velocities.
    pub fn new(mean_velocity: Vec<[f64; 3]>) -> Self {
        Self { mean_velocity }
    }

    /// Velocity fluctuation u' = u - `u`.
    pub fn fluctuation(&self, idx: usize, velocity: [f64; 3]) -> [f64; 3] {
        if idx < self.mean_velocity.len() {
            let mean = self.mean_velocity[idx];
            [
                velocity[0] - mean[0],
                velocity[1] - mean[1],
                velocity[2] - mean[2],
            ]
        } else {
            velocity
        }
    }

    /// Turbulent kinetic energy per unit mass: k = 0.5 * <u'·u'>.
    pub fn tke(fluctuations: &[[f64; 3]]) -> f64 {
        let n = fluctuations.len();
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = fluctuations
            .iter()
            .map(|u| u[0] * u[0] + u[1] * u[1] + u[2] * u[2])
            .sum();
        0.5 * sum / n as f64
    }

    /// Reynolds stress tensor <u'_i u'_j> (3×3 row-major).
    pub fn reynolds_stress(fluctuations: &[[f64; 3]]) -> [f64; 9] {
        let n = fluctuations.len();
        if n == 0 {
            return [0.0; 9];
        }
        let mut r = [0.0_f64; 9];
        for u in fluctuations {
            for i in 0..3 {
                for j in 0..3 {
                    r[i * 3 + j] += u[i] * u[j];
                }
            }
        }
        r.map(|v| v / n as f64)
    }

    /// Second-order longitudinal structure function D_LL(r).
    ///
    /// D_LL(r) = <\[u(x+r) - u(x)\]²>
    ///
    /// Computes approximate value from particle pairs within (r-dr, r+dr).
    pub fn structure_function_ll(
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        r_target: f64,
        dr: f64,
    ) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        let n = positions.len().min(velocities.len());
        for i in 0..n {
            for j in (i + 1)..n {
                let r_ij = dist3(positions[i], positions[j]);
                if (r_ij - r_target).abs() < dr {
                    let r_hat = if r_ij > 1e-14 {
                        [
                            (positions[j][0] - positions[i][0]) / r_ij,
                            (positions[j][1] - positions[i][1]) / r_ij,
                            (positions[j][2] - positions[i][2]) / r_ij,
                        ]
                    } else {
                        [1.0, 0.0, 0.0]
                    };
                    let du = [
                        velocities[j][0] - velocities[i][0],
                        velocities[j][1] - velocities[i][1],
                        velocities[j][2] - velocities[i][2],
                    ];
                    let du_l = du[0] * r_hat[0] + du[1] * r_hat[1] + du[2] * r_hat[2];
                    sum += du_l * du_l;
                    count += 1;
                }
            }
        }
        if count > 0 { sum / count as f64 } else { 0.0 }
    }

    /// Turbulence dissipation rate estimate ε ≈ 15ν * <(∂u/∂x)²>.
    ///
    /// Uses isotropic approximation.
    pub fn dissipation_rate_estimate(nu: f64, strain_rate_sq: f64) -> f64 {
        2.0 * nu * strain_rate_sq
    }
}

// ============================================================================
// Free Surface Detection
// ============================================================================

/// SPH free surface detection via density criterion and color function.
///
/// A particle is flagged as free surface if its density drops below
/// a threshold fraction of the reference density.
pub struct SphFreeSurface {
    /// Reference density ρ₀ (kg/m³).
    pub rho0: f64,
    /// Detection threshold fraction (typically 0.85–0.95).
    pub threshold: f64,
    /// Smoothing length h (m).
    pub h: f64,
}

impl SphFreeSurface {
    /// Create a free surface detector.
    pub fn new(rho0: f64, threshold: f64, h: f64) -> Self {
        Self { rho0, threshold, h }
    }

    /// Check if particle i is a free surface particle.
    ///
    /// Returns true if ρᵢ < threshold * ρ₀.
    pub fn is_free_surface(&self, rho: f64) -> bool {
        rho < self.threshold * self.rho0
    }

    /// Classify all particles and return boolean vector.
    pub fn classify_particles(&self, densities: &[f64]) -> Vec<bool> {
        densities
            .iter()
            .map(|&rho| self.is_free_surface(rho))
            .collect()
    }

    /// Color function c_i: 1 for interior, 0 for free surface.
    pub fn color_function(&self, rho: f64) -> f64 {
        if self.is_free_surface(rho) { 0.0 } else { 1.0 }
    }

    /// Surface normal estimate from color function gradient.
    ///
    /// n_i = -∇c_i / |∇c_i|  (points outward from fluid).
    pub fn surface_normal_from_gradient(grad_c: [f64; 3]) -> [f64; 3] {
        let mag = dist_mag(grad_c);
        if mag < 1e-14 {
            [0.0, 0.0, 1.0]
        } else {
            [-grad_c[0] / mag, -grad_c[1] / mag, -grad_c[2] / mag]
        }
    }

    /// SPH curvature from surface normal divergence (Brackbill et al.).
    ///
    /// κ ≈ -∇ · n̂  (approximate: n̂ already normalized).
    pub fn curvature_estimate(&self, div_n: f64) -> f64 {
        -div_n
    }
}

// ============================================================================
// Particle Sorting
// ============================================================================

/// SPH particle sorting using Z-curve (Morton code) for cache-friendly access.
///
/// Morton codes interleave the bits of quantized x, y, z coordinates,
/// mapping 3D positions to a 1D ordering that preserves spatial locality.
pub struct SphParticleSorting {
    /// Domain minimum corner.
    pub domain_min: [f64; 3],
    /// Domain size.
    pub domain_size: [f64; 3],
    /// Resolution (number of cells per axis).
    pub resolution: u32,
}

impl SphParticleSorting {
    /// Create a particle sorter for the given domain.
    pub fn new(domain_min: [f64; 3], domain_size: [f64; 3], resolution: u32) -> Self {
        Self {
            domain_min,
            domain_size,
            resolution,
        }
    }

    /// Expand a 10-bit integer to 30 bits by inserting 2 zero bits between each bit.
    pub fn expand_bits(x: u32) -> u32 {
        let x = x & 0x3ff;
        let x = (x | (x << 16)) & 0x30000ff;
        let x = (x | (x << 8)) & 0x0300f00f;
        let x = (x | (x << 4)) & 0x30c30c3;
        (x | (x << 2)) & 0x9249249
    }

    /// Compute 30-bit Morton code for a 3D position.
    pub fn morton_code(&self, pos: [f64; 3]) -> u32 {
        let res = self.resolution as f64;
        let xi = ((pos[0] - self.domain_min[0]) / self.domain_size[0] * res).clamp(0.0, res - 1.0)
            as u32;
        let yi = ((pos[1] - self.domain_min[1]) / self.domain_size[1] * res).clamp(0.0, res - 1.0)
            as u32;
        let zi = ((pos[2] - self.domain_min[2]) / self.domain_size[2] * res).clamp(0.0, res - 1.0)
            as u32;
        Self::expand_bits(xi) | (Self::expand_bits(yi) << 1) | (Self::expand_bits(zi) << 2)
    }

    /// Sort particle indices by Morton code.
    ///
    /// Returns a permutation index vector sorted by Z-curve order.
    pub fn sort_indices(&self, positions: &[[f64; 3]]) -> Vec<usize> {
        let mut keys: Vec<(u32, usize)> = positions
            .iter()
            .enumerate()
            .map(|(i, &pos)| (self.morton_code(pos), i))
            .collect();
        keys.sort_unstable_by_key(|&(code, _)| code);
        keys.into_iter().map(|(_, i)| i).collect()
    }

    /// Reorder a particle array by the given permutation.
    pub fn reorder<T: Clone>(data: &[T], indices: &[usize]) -> Vec<T> {
        indices.iter().map(|&i| data[i].clone()).collect()
    }
}

// ============================================================================
// Diagnostics
// ============================================================================

/// SPH simulation diagnostics.
///
/// Computes per-timestep diagnostic quantities for monitoring
/// simulation health and convergence.
#[derive(Debug, Clone, Default)]
pub struct SphDiagnostics {
    /// Maximum particle speed (m/s).
    pub max_velocity: f64,
    /// Maximum density error |ρ - ρ₀| / ρ₀.
    pub max_density_error: f64,
    /// Average number of neighbors per particle.
    pub avg_neighbor_count: f64,
    /// Disorder parameter (deviation from ideal lattice).
    pub disorder_parameter: f64,
}

impl SphDiagnostics {
    /// Compute diagnostics from particle state.
    pub fn compute(
        velocities: &[[f64; 3]],
        densities: &[f64],
        rho0: f64,
        neighbor_counts: &[usize],
    ) -> Self {
        let max_vel = velocities
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0_f64, f64::max);

        let max_dens_err = densities
            .iter()
            .map(|&rho| (rho - rho0).abs() / rho0)
            .fold(0.0_f64, f64::max);

        let avg_nb = if neighbor_counts.is_empty() {
            0.0
        } else {
            neighbor_counts.iter().sum::<usize>() as f64 / neighbor_counts.len() as f64
        };

        // Disorder is the mean relative spread of neighbour counts about their
        // own mean (the data-derived "ideal" lattice count). A perfectly
        // regular packing has every particle at the mean count ⇒ disorder ≈ 0,
        // while an irregular neighbour-count distribution yields a positive
        // value. When the mean count is zero (no neighbours at all) there is no
        // meaningful target and the disorder is defined to be zero.
        let disorder = if avg_nb > 0.0 {
            Self::compute_disorder(neighbor_counts, avg_nb)
        } else {
            0.0
        };

        Self {
            max_velocity: max_vel,
            max_density_error: max_dens_err,
            avg_neighbor_count: avg_nb,
            disorder_parameter: disorder,
        }
    }

    /// Compute disorder parameter: average relative deviation from a target
    /// neighbour count.
    ///
    /// `disorder = mean_i |N_i − N_target| / N_target`
    ///
    /// `N_target` is the reference (ideal-lattice) neighbour count; pass the
    /// mean neighbour count to measure deviation from the configuration's own
    /// average. Returns `0.0` for an empty input or a non-positive target.
    pub fn compute_disorder(neighbor_counts: &[usize], n_target: f64) -> f64 {
        if neighbor_counts.is_empty() || n_target <= 0.0 {
            return 0.0;
        }
        let sum: f64 = neighbor_counts
            .iter()
            .map(|&n| (n as f64 - n_target).abs() / n_target)
            .sum();
        sum / neighbor_counts.len() as f64
    }
}

// ============================================================================
// VTK Exporter
// ============================================================================

/// SPH particle data exporter to VTK PolyData XML format.
///
/// Writes particle positions, velocities, densities, and pressures
/// to a VTK .vtp file for visualization in ParaView.
pub struct SphVtkExporter {
    /// Output file path prefix.
    pub path_prefix: String,
    /// Current frame number.
    pub frame: usize,
}

impl SphVtkExporter {
    /// Create a new VTK exporter.
    pub fn new(path_prefix: &str) -> Self {
        Self {
            path_prefix: path_prefix.to_string(),
            frame: 0,
        }
    }

    /// Generate VTK PolyData XML string for a particle set.
    ///
    /// Returns the XML content as a String suitable for writing to a .vtp file.
    pub fn generate_vtp(
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        densities: &[f64],
        pressures: &[f64],
        time: f64,
    ) -> String {
        let n = positions.len();
        let mut xml = String::new();
        xml.push_str("<?xml version=\"1.0\"?>\n");
        xml.push_str("<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        xml.push_str(&format!("  <!-- Time: {:.6e} -->\n", time));
        xml.push_str("  <PolyData>\n");
        xml.push_str(&format!(
            "    <Piece NumberOfPoints=\"{}\" NumberOfVerts=\"{}\">\n",
            n, n
        ));
        xml.push_str("      <Points>\n");
        xml.push_str(
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n",
        );
        for p in positions {
            xml.push_str(&format!(
                "          {:.6e} {:.6e} {:.6e}\n",
                p[0], p[1], p[2]
            ));
        }
        xml.push_str("        </DataArray>\n");
        xml.push_str("      </Points>\n");
        xml.push_str("      <Verts>\n");
        xml.push_str("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
        for i in 0..n {
            xml.push_str(&format!("          {}\n", i));
        }
        xml.push_str("        </DataArray>\n");
        xml.push_str("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
        for i in 1..=n {
            xml.push_str(&format!("          {}\n", i));
        }
        xml.push_str("        </DataArray>\n");
        xml.push_str("      </Verts>\n");
        xml.push_str("      <PointData>\n");
        // Velocity
        xml.push_str("        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n");
        for v in velocities {
            xml.push_str(&format!(
                "          {:.6e} {:.6e} {:.6e}\n",
                v[0], v[1], v[2]
            ));
        }
        xml.push_str("        </DataArray>\n");
        // Density
        xml.push_str("        <DataArray type=\"Float64\" Name=\"density\" format=\"ascii\">\n");
        for &rho in densities {
            xml.push_str(&format!("          {:.6e}\n", rho));
        }
        xml.push_str("        </DataArray>\n");
        // Pressure
        xml.push_str("        <DataArray type=\"Float64\" Name=\"pressure\" format=\"ascii\">\n");
        for &p in pressures {
            xml.push_str(&format!("          {:.6e}\n", p));
        }
        xml.push_str("        </DataArray>\n");
        xml.push_str("      </PointData>\n");
        xml.push_str("    </Piece>\n");
        xml.push_str("  </PolyData>\n");
        xml.push_str("</VTKFile>\n");
        xml
    }

    /// Get output filename for current frame.
    pub fn filename(&self) -> String {
        format!("{}_{:06}.vtp", self.path_prefix, self.frame)
    }

    /// Advance frame counter.
    pub fn next_frame(&mut self) {
        self.frame += 1;
    }
}

// ============================================================================
// Force Balance
// ============================================================================

/// SPH force balance diagnostics.
///
/// Decomposes total particle forces into pressure, viscous, body, and
/// surface tension contributions for monitoring and debugging.
#[derive(Debug, Clone, Default)]
pub struct SphForceBalance {
    /// Pressure force per particle (N).
    pub pressure_forces: Vec<[f64; 3]>,
    /// Viscous force per particle (N).
    pub viscous_forces: Vec<[f64; 3]>,
    /// Body force per particle (N).
    pub body_forces: Vec<[f64; 3]>,
    /// Surface tension force per particle (N).
    pub surface_tension_forces: Vec<[f64; 3]>,
}

impl SphForceBalance {
    /// Create empty force balance for n particles.
    pub fn new(n: usize) -> Self {
        Self {
            pressure_forces: vec![[0.0; 3]; n],
            viscous_forces: vec![[0.0; 3]; n],
            body_forces: vec![[0.0; 3]; n],
            surface_tension_forces: vec![[0.0; 3]; n],
        }
    }

    /// Total force on particle i.
    pub fn total_force(&self, i: usize) -> [f64; 3] {
        let n = self.pressure_forces.len();
        if i >= n {
            return [0.0; 3];
        }
        [
            self.pressure_forces[i][0]
                + self.viscous_forces[i][0]
                + self.body_forces[i][0]
                + self.surface_tension_forces[i][0],
            self.pressure_forces[i][1]
                + self.viscous_forces[i][1]
                + self.body_forces[i][1]
                + self.surface_tension_forces[i][1],
            self.pressure_forces[i][2]
                + self.viscous_forces[i][2]
                + self.body_forces[i][2]
                + self.surface_tension_forces[i][2],
        ]
    }

    /// L2 norm of total force vector summed over all particles.
    pub fn total_force_l2(&self) -> f64 {
        let n = self.pressure_forces.len();
        (0..n)
            .map(|i| {
                let f = self.total_force(i);
                f[0] * f[0] + f[1] * f[1] + f[2] * f[2]
            })
            .sum::<f64>()
            .sqrt()
    }

    /// Sum of pressure forces (for momentum balance check).
    pub fn sum_pressure_forces(&self) -> [f64; 3] {
        let mut s = [0.0_f64; 3];
        for f in &self.pressure_forces {
            s[0] += f[0];
            s[1] += f[1];
            s[2] += f[2];
        }
        s
    }
}

// ============================================================================
// Convergence
// ============================================================================

/// SPH convergence analysis against an analytical solution.
///
/// Computes L2 error norms and estimates convergence rates from
/// a series of resolutions.
pub struct SphConvergence {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Smoothing length h.
    pub h: f64,
}

impl SphConvergence {
    /// Create a convergence analyzer.
    pub fn new(positions: Vec<[f64; 3]>, h: f64) -> Self {
        Self { positions, h }
    }

    /// Compute L2 error norm between numerical and analytical fields.
    ///
    /// L2 = sqrt(Σ (φ_num - φ_exact)² * vol_i) / sqrt(Σ φ_exact² * vol_i)
    pub fn l2_error(numerical: &[f64], analytical: &[f64], volumes: &[f64]) -> f64 {
        let n = numerical.len().min(analytical.len()).min(volumes.len());
        let err_sq: f64 = (0..n)
            .map(|i| {
                let e = numerical[i] - analytical[i];
                e * e * volumes[i]
            })
            .sum();
        let norm_sq: f64 = (0..n)
            .map(|i| analytical[i] * analytical[i] * volumes[i])
            .sum();
        if norm_sq < 1e-30 {
            err_sq.sqrt()
        } else {
            (err_sq / norm_sq).sqrt()
        }
    }

    /// L1 error norm.
    pub fn l1_error(numerical: &[f64], analytical: &[f64]) -> f64 {
        let n = numerical.len().min(analytical.len());
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = (0..n).map(|i| (numerical[i] - analytical[i]).abs()).sum();
        let denom: f64 = (0..n).map(|i| analytical[i].abs()).sum();
        if denom < 1e-30 { sum } else { sum / denom }
    }

    /// L∞ (maximum) error.
    pub fn linf_error(numerical: &[f64], analytical: &[f64]) -> f64 {
        let n = numerical.len().min(analytical.len());
        (0..n)
            .map(|i| (numerical[i] - analytical[i]).abs())
            .fold(0.0_f64, f64::max)
    }

    /// Estimate convergence rate from two resolutions.
    ///
    /// p ≈ log(e1 / e2) / log(h1 / h2)
    pub fn convergence_rate(e1: f64, e2: f64, h1: f64, h2: f64) -> f64 {
        if e1 < 1e-30 || e2 < 1e-30 || (h1 / h2).abs() < 1e-10 {
            return 0.0;
        }
        (e1 / e2).ln() / (h1 / h2).ln()
    }
}

// ============================================================================
// Helper Functions
// ============================================================================

/// SPH interpolation of a scalar field at query point x.
///
/// φ(x) = Σⱼ (mⱼ/ρⱼ) * φⱼ * W(|x - xⱼ|, h)
pub fn sph_interpolate_scalar(
    x: [f64; 3],
    positions: &[[f64; 3]],
    field: &[f64],
    masses: &[f64],
    densities: &[f64],
    h: f64,
) -> f64 {
    let mut phi = 0.0;
    let n = positions
        .len()
        .min(field.len())
        .min(masses.len())
        .min(densities.len());
    for j in 0..n {
        let r = dist3(x, positions[j]);
        let w = cubic_spline_1d(r, h);
        phi += masses[j] / densities[j] * field[j] * w;
    }
    phi
}

/// SPH divergence of a vector field at query point x.
///
/// ∇·f(x) ≈ Σⱼ (mⱼ/ρⱼ) * f_j · ∇W_ij
pub fn sph_divergence(
    x: [f64; 3],
    positions: &[[f64; 3]],
    field: &[[f64; 3]],
    masses: &[f64],
    densities: &[f64],
    h: f64,
) -> f64 {
    let mut div = 0.0;
    let n = positions
        .len()
        .min(field.len())
        .min(masses.len())
        .min(densities.len());
    for j in 0..n {
        let r_ij = [
            x[0] - positions[j][0],
            x[1] - positions[j][1],
            x[2] - positions[j][2],
        ];
        let r = dist_mag(r_ij);
        if r < 1e-14 {
            continue;
        }
        let dw_dr = cubic_spline_grad_1d(r, h);
        let grad_w = [
            dw_dr * r_ij[0] / r,
            dw_dr * r_ij[1] / r,
            dw_dr * r_ij[2] / r,
        ];
        let c = masses[j] / densities[j];
        div += c * (field[j][0] * grad_w[0] + field[j][1] * grad_w[1] + field[j][2] * grad_w[2]);
    }
    div
}

/// SPH curl of a vector field at query point x.
///
/// (∇×f)_k = Σⱼ (mⱼ/ρⱼ) * (f_j × ∇W_ij)_k
pub fn sph_curl(
    x: [f64; 3],
    positions: &[[f64; 3]],
    field: &[[f64; 3]],
    masses: &[f64],
    densities: &[f64],
    h: f64,
) -> [f64; 3] {
    let mut curl = [0.0_f64; 3];
    let n = positions
        .len()
        .min(field.len())
        .min(masses.len())
        .min(densities.len());
    for j in 0..n {
        let r_ij = [
            x[0] - positions[j][0],
            x[1] - positions[j][1],
            x[2] - positions[j][2],
        ];
        let r = dist_mag(r_ij);
        if r < 1e-14 {
            continue;
        }
        let dw_dr = cubic_spline_grad_1d(r, h);
        let gw = [
            dw_dr * r_ij[0] / r,
            dw_dr * r_ij[1] / r,
            dw_dr * r_ij[2] / r,
        ];
        let c = masses[j] / densities[j];
        let fj = field[j];
        // cross product fj × gw
        curl[0] += c * (fj[1] * gw[2] - fj[2] * gw[1]);
        curl[1] += c * (fj[2] * gw[0] - fj[0] * gw[2]);
        curl[2] += c * (fj[0] * gw[1] - fj[1] * gw[0]);
    }
    curl
}

/// SPH Laplacian of a scalar field at particle i.
///
/// ∇²φ_i ≈ 2 * Σⱼ (mⱼ/ρⱼ) * (φⱼ - φᵢ) * (r_ij · ∇W_ij) / |r_ij|²
pub fn sph_laplacian(
    xi: [f64; 3],
    phi_i: f64,
    positions: &[[f64; 3]],
    field: &[f64],
    masses: &[f64],
    densities: &[f64],
    h: f64,
) -> f64 {
    let mut lap = 0.0;
    let n = positions
        .len()
        .min(field.len())
        .min(masses.len())
        .min(densities.len());
    for j in 0..n {
        let r_ij = [
            xi[0] - positions[j][0],
            xi[1] - positions[j][1],
            xi[2] - positions[j][2],
        ];
        let r = dist_mag(r_ij);
        if r < 1e-14 {
            continue;
        }
        let dw_dr = cubic_spline_grad_1d(r, h);
        let r_dot_grad = r_ij[0] * (dw_dr * r_ij[0] / r)
            + r_ij[1] * (dw_dr * r_ij[1] / r)
            + r_ij[2] * (dw_dr * r_ij[2] / r);
        let c = masses[j] / densities[j];
        lap += 2.0 * c * (field[j] - phi_i) * r_dot_grad / (r * r);
    }
    lap
}

// Private helpers

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn dist_mag(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn cubic_spline_1d(r: f64, h: f64) -> f64 {
    SphFieldInterp::cubic_spline_kernel(r, h)
}

fn cubic_spline_grad_1d(r: f64, h: f64) -> f64 {
    VelocityGradient::kernel_gradient_magnitude(r, h)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // SphFieldInterp tests
    #[test]
    fn test_cubic_spline_kernel_at_zero() {
        let w = SphFieldInterp::cubic_spline_kernel(0.0, 1.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_cubic_spline_kernel_beyond_support() {
        let w = SphFieldInterp::cubic_spline_kernel(3.0, 1.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_cubic_spline_kernel_positive() {
        let w = SphFieldInterp::cubic_spline_kernel(0.5, 1.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_sph_interpolate_single_particle() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let field = vec![5.0_f64];
        let masses = vec![1e-3];
        let densities = vec![1000.0];
        let interp = SphFieldInterp::new(0.1, masses.clone(), densities.clone(), positions.clone());
        let val = interp.interpolate_scalar([0.0, 0.0, 0.0], &field);
        assert!(val.is_finite());
    }

    #[test]
    fn test_sph_field_interp_number_density() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let masses = vec![1e-3];
        let densities = vec![1000.0];
        let interp = SphFieldInterp::new(0.1, masses, densities, positions);
        let nd = interp.number_density([0.0, 0.0, 0.0]);
        assert!(nd > 0.0);
    }

    // SphStatistics tests
    #[test]
    fn test_sph_statistics_ke_positive() {
        let masses = vec![1.0];
        let velocities = vec![[1.0, 0.0, 0.0]];
        let positions = vec![[0.0, 0.0, 0.0]];
        let stats = SphStatistics::compute(&masses, &velocities, &positions, [0.0, -9.81, 0.0]);
        assert!((stats.kinetic_energy - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_sph_statistics_zero_velocity() {
        let masses = vec![1.0];
        let velocities = vec![[0.0, 0.0, 0.0]];
        let positions = vec![[0.0, 0.0, 0.0]];
        let stats = SphStatistics::compute(&masses, &velocities, &positions, [0.0, -9.81, 0.0]);
        assert_eq!(stats.kinetic_energy, 0.0);
    }

    #[test]
    fn test_sph_statistics_total_energy() {
        let stats = SphStatistics {
            kinetic_energy: 10.0,
            potential_energy: 5.0,
            ..Default::default()
        };
        assert!((stats.total_energy() - 15.0).abs() < 1e-14);
    }

    #[test]
    fn test_sph_statistics_momentum_magnitude() {
        let stats = SphStatistics {
            momentum: [3.0, 4.0, 0.0],
            ..Default::default()
        };
        assert!((stats.momentum_magnitude() - 5.0).abs() < 1e-10);
    }

    // VelocityGradient tests
    #[test]
    fn test_strain_rate_symmetry() {
        let l = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let s = VelocityGradient::strain_rate(&l);
        // S should be symmetric
        assert!((s[1] - s[3]).abs() < 1e-14);
        assert!((s[2] - s[2 * 3]).abs() < 1e-14);
    }

    #[test]
    fn test_rotation_rate_antisymmetry() {
        let l = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let omega = VelocityGradient::rotation_rate(&l);
        assert!((omega[1] + omega[3]).abs() < 1e-14);
    }

    #[test]
    fn test_rotation_rate_diagonal_zero() {
        let l = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let omega = VelocityGradient::rotation_rate(&l);
        assert_eq!(omega[0], 0.0);
        assert_eq!(omega[4], 0.0);
        assert_eq!(omega[8], 0.0);
    }

    #[test]
    fn test_vorticity_identity() {
        // Simple shear: ∂u/∂y = 1, all others 0
        let l = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vort = VelocityGradient::vorticity(&l);
        assert!((vort[2] - (-1.0)).abs() < 1e-14);
    }

    // SphTurbulence tests
    #[test]
    fn test_tke_zero_fluctuations() {
        let fluct = vec![[0.0_f64; 3]; 10];
        assert_eq!(SphTurbulence::tke(&fluct), 0.0);
    }

    #[test]
    fn test_tke_unit_fluctuations() {
        let fluct = vec![[1.0, 0.0, 0.0]; 4];
        assert!((SphTurbulence::tke(&fluct) - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_reynolds_stress_symmetry() {
        let fluct = vec![[1.0, 2.0, 3.0], [0.5, 1.0, 1.5]];
        let r = SphTurbulence::reynolds_stress(&fluct);
        assert!((r[1] - r[3]).abs() < 1e-12);
    }

    #[test]
    fn test_structure_function_no_pairs() {
        let pos = vec![[0.0_f64; 3]; 3];
        let vel = vec![[1.0, 0.0, 0.0]; 3];
        let dll = SphTurbulence::structure_function_ll(&pos, &vel, 1.0, 0.1);
        assert_eq!(dll, 0.0);
    }

    // SphFreeSurface tests
    #[test]
    fn test_free_surface_detect() {
        let fs = SphFreeSurface::new(1000.0, 0.9, 0.01);
        assert!(fs.is_free_surface(800.0));
        assert!(!fs.is_free_surface(950.0));
    }

    #[test]
    fn test_free_surface_classify() {
        let fs = SphFreeSurface::new(1000.0, 0.9, 0.01);
        let result = fs.classify_particles(&[800.0, 950.0]);
        assert!(result[0]);
        assert!(!result[1]);
    }

    #[test]
    fn test_free_surface_normal_from_gradient() {
        let n = SphFreeSurface::surface_normal_from_gradient([0.0, 0.0, 1.0]);
        assert!((n[2] - (-1.0)).abs() < 1e-14);
    }

    // SphParticleSorting tests
    #[test]
    fn test_morton_expand_bits() {
        assert_eq!(SphParticleSorting::expand_bits(0), 0);
        assert_eq!(SphParticleSorting::expand_bits(1), 1);
    }

    #[test]
    fn test_morton_code_origin() {
        let sorter = SphParticleSorting::new([0.0; 3], [1.0; 3], 1024);
        let code = sorter.morton_code([0.0, 0.0, 0.0]);
        assert_eq!(code, 0);
    }

    #[test]
    fn test_sort_indices_length() {
        let sorter = SphParticleSorting::new([0.0; 3], [1.0; 3], 1024);
        let positions = vec![[0.5, 0.1, 0.2], [0.1, 0.9, 0.3], [0.8, 0.5, 0.0]];
        let indices = sorter.sort_indices(&positions);
        assert_eq!(indices.len(), 3);
    }

    #[test]
    fn test_reorder() {
        let data = vec![10, 20, 30];
        let indices = vec![2, 0, 1];
        let result = SphParticleSorting::reorder(&data, &indices);
        assert_eq!(result, vec![30, 10, 20]);
    }

    // SphDiagnostics tests
    #[test]
    fn test_diagnostics_max_velocity() {
        let vels = vec![[3.0_f64, 4.0, 0.0], [1.0, 0.0, 0.0]];
        let dens = vec![1000.0, 1000.0];
        let diag = SphDiagnostics::compute(&vels, &dens, 1000.0, &[10, 10]);
        assert!((diag.max_velocity - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_diagnostics_density_error() {
        let vels = vec![[0.0; 3]];
        let dens = vec![1100.0];
        let diag = SphDiagnostics::compute(&vels, &dens, 1000.0, &[10]);
        assert!((diag.max_density_error - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_diagnostics_disorder() {
        let nb = vec![10, 12, 8];
        let d = SphDiagnostics::compute_disorder(&nb, 10.0);
        assert!(d >= 0.0);
        // |10-10|/10 + |12-10|/10 + |8-10|/10 = 0 + 0.2 + 0.2 = 0.4, /3
        assert!(
            (d - 0.4 / 3.0).abs() < 1e-12,
            "disorder should be 0.4/3, got {d}"
        );
    }

    #[test]
    fn test_diagnostics_disorder_is_computed_not_hardcoded() {
        let vels = vec![[0.0_f64; 3]; 6];
        let dens = vec![1000.0; 6];

        // Perfectly regular neighbour-count distribution: disorder must be ~0.
        let regular = [20usize, 20, 20, 20, 20, 20];
        let diag_regular = SphDiagnostics::compute(&vels, &dens, 1000.0, &regular);
        assert!(
            diag_regular.disorder_parameter < 1e-12,
            "regular packing must have ~zero disorder, got {}",
            diag_regular.disorder_parameter
        );
        assert!((diag_regular.avg_neighbor_count - 20.0).abs() < 1e-12);

        // Irregular/disordered neighbour-count distribution: disorder must be
        // strictly positive. If `compute` still hard-coded 0.0 this fails.
        let irregular = [5usize, 35, 8, 31, 12, 29];
        let diag_irregular = SphDiagnostics::compute(&vels, &dens, 1000.0, &irregular);
        assert!(
            diag_irregular.disorder_parameter > 0.3,
            "disordered packing must have a clearly positive disorder, got {}",
            diag_irregular.disorder_parameter
        );

        // Cross-check the integrated value equals the helper applied to the
        // data-derived mean target (proves `compute` actually calls it).
        let mean = irregular.iter().sum::<usize>() as f64 / irregular.len() as f64;
        let expected = SphDiagnostics::compute_disorder(&irregular, mean);
        assert!(
            (diag_irregular.disorder_parameter - expected).abs() < 1e-12,
            "compute() disorder ({}) must match compute_disorder() ({})",
            diag_irregular.disorder_parameter,
            expected
        );
    }

    // SphVtkExporter tests
    #[test]
    fn test_vtk_generate_vtp_not_empty() {
        let pos = vec![[0.0, 0.0, 0.0]];
        let vel = vec![[1.0, 0.0, 0.0]];
        let den = vec![1000.0];
        let pre = vec![0.0];
        let xml = SphVtkExporter::generate_vtp(&pos, &vel, &den, &pre, 0.0);
        assert!(xml.contains("VTKFile"));
        assert!(xml.contains("PolyData"));
    }

    #[test]
    fn test_vtk_filename() {
        let exp = SphVtkExporter::new("output/particles");
        assert_eq!(exp.filename(), "output/particles_000000.vtp");
    }

    #[test]
    fn test_vtk_next_frame() {
        let mut exp = SphVtkExporter::new("out");
        exp.next_frame();
        assert_eq!(exp.frame, 1);
    }

    // SphForceBalance tests
    #[test]
    fn test_force_balance_total_zero() {
        let fb = SphForceBalance::new(2);
        let f = fb.total_force(0);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn test_force_balance_l2_norm() {
        let mut fb = SphForceBalance::new(1);
        fb.pressure_forces[0] = [3.0, 4.0, 0.0];
        let l2 = fb.total_force_l2();
        assert!((l2 - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_force_balance_sum_pressure() {
        let mut fb = SphForceBalance::new(2);
        fb.pressure_forces[0] = [1.0, 0.0, 0.0];
        fb.pressure_forces[1] = [2.0, 0.0, 0.0];
        let s = fb.sum_pressure_forces();
        assert!((s[0] - 3.0).abs() < 1e-14);
    }

    // SphConvergence tests
    #[test]
    fn test_l2_error_zero() {
        let num = vec![1.0, 2.0, 3.0];
        let ana = vec![1.0, 2.0, 3.0];
        let vol = vec![1.0, 1.0, 1.0];
        assert_eq!(SphConvergence::l2_error(&num, &ana, &vol), 0.0);
    }

    #[test]
    fn test_l2_error_nonzero() {
        let num = vec![1.1, 2.0, 3.0];
        let ana = vec![1.0, 2.0, 3.0];
        let vol = vec![1.0, 1.0, 1.0];
        let err = SphConvergence::l2_error(&num, &ana, &vol);
        assert!(err > 0.0);
    }

    #[test]
    fn test_l1_error_zero() {
        let num = vec![1.0; 5];
        let ana = vec![1.0; 5];
        assert_eq!(SphConvergence::l1_error(&num, &ana), 0.0);
    }

    #[test]
    fn test_linf_error() {
        let num = vec![1.0, 3.0, 5.0];
        let ana = vec![1.0, 2.0, 5.0];
        assert!((SphConvergence::linf_error(&num, &ana) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_convergence_rate_second_order() {
        // If e2 = e1/4 and h2 = h1/2 → rate should be 2
        let rate = SphConvergence::convergence_rate(0.04, 0.01, 0.2, 0.1);
        assert!((rate - 2.0).abs() < 1e-10);
    }

    // Helper function tests
    #[test]
    fn test_sph_interpolate_scalar_zero_field() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let field = vec![0.0_f64];
        let masses = vec![1e-3];
        let densities = vec![1000.0];
        let val = sph_interpolate_scalar(
            [0.0, 0.0, 0.0],
            &positions,
            &field,
            &masses,
            &densities,
            0.1,
        );
        assert_eq!(val, 0.0);
    }

    #[test]
    fn test_sph_divergence_zero_field() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let field = vec![[0.0_f64; 3]];
        let masses = vec![1e-3];
        let densities = vec![1000.0];
        let div = sph_divergence(
            [0.0, 0.0, 0.0],
            &positions,
            &field,
            &masses,
            &densities,
            0.1,
        );
        assert_eq!(div, 0.0);
    }

    #[test]
    fn test_sph_curl_zero_field() {
        let positions = vec![[0.0, 0.0, 0.0]];
        let field = vec![[0.0_f64; 3]];
        let masses = vec![1e-3];
        let densities = vec![1000.0];
        let c = sph_curl(
            [0.0, 0.0, 0.0],
            &positions,
            &field,
            &masses,
            &densities,
            0.1,
        );
        assert_eq!(c, [0.0; 3]);
    }

    #[test]
    fn test_sph_laplacian_uniform_field() {
        // For a uniform field, Laplacian should be near zero
        let positions = vec![[0.05, 0.0, 0.0], [-0.05, 0.0, 0.0]];
        let field = vec![5.0_f64, 5.0_f64];
        let masses = vec![1e-3, 1e-3];
        let densities = vec![1000.0, 1000.0];
        let lap = sph_laplacian(
            [0.0, 0.0, 0.0],
            5.0,
            &positions,
            &field,
            &masses,
            &densities,
            0.1,
        );
        assert_eq!(lap, 0.0);
    }
}
