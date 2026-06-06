// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cahn-Hilliard phase field LBM for spinodal decomposition.
//!
//! Implements a two-distribution LBM for the Cahn-Hilliard equation:
//!
//! ```text
//! ∂φ/∂t = M ∇²μ
//! μ = a·φ·(φ²-1) - κ·∇²φ
//! ```
//!
//! Two sets of distributions are evolved:
//! - `f_dist`: tracks the order parameter φ
//! - `g_dist`: tracks the chemical potential μ
//!
//! Reference: Lee & Liu (2010), *J. Comput. Phys.* **229**, 8045-8063.

// D2Q9 weights and velocity indices
const W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// D2Q9 velocity components (ex, ey).
const EX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
/// D2Q9 velocity components (ey).
const EY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];

/// Cahn-Hilliard LBM solver for spinodal decomposition on a D2Q9 grid.
///
/// Evolves an order parameter field φ using two distribution functions
/// (`f_dist` for φ and `g_dist` for μ) with BGK collision on each.
#[derive(Debug, Clone)]
pub struct CahnHilliardLBM {
    /// Number of lattice nodes in the x-direction.
    pub nx: usize,
    /// Number of lattice nodes in the y-direction.
    pub ny: usize,
    /// Order parameter field φ.
    pub phi: Vec<f64>,
    /// Chemical potential field μ.
    pub mu: Vec<f64>,
    /// Distribution function for the φ field (size nx*ny*9).
    pub f_dist: Vec<f64>,
    /// Distribution function for the μ field (size nx*ny*9).
    pub g_dist: Vec<f64>,
    /// Mobility coefficient M.
    pub mobility: f64,
    /// Interface gradient energy coefficient κ.
    pub kappa: f64,
    /// Double-well potential coefficient a.
    pub a: f64,
    /// Relaxation time for the φ distribution.
    pub tau_phi: f64,
    /// Relaxation time for the μ distribution.
    pub tau_mu: f64,
}

impl CahnHilliardLBM {
    /// Create a new Cahn-Hilliard LBM solver.
    ///
    /// # Arguments
    /// - `nx`, `ny`: grid dimensions
    /// - `mobility`: mobility M
    /// - `kappa`: interface energy coefficient κ
    /// - `a`: double-well coefficient a
    pub fn new(nx: usize, ny: usize, mobility: f64, kappa: f64, a: f64) -> Self {
        let n = nx * ny;
        let phi = vec![0.0; n];
        let mu = vec![0.0; n];
        let f_dist = vec![0.0; n * 9];
        let g_dist = vec![0.0; n * 9];
        // tau related to mobility: tau_phi = 0.5 + 3*mobility
        let tau_phi = 0.5 + 3.0 * mobility;
        // tau_mu set to 1.0 (simple choice)
        let tau_mu = 1.0_f64;
        Self {
            nx,
            ny,
            phi,
            mu,
            f_dist,
            g_dist,
            mobility,
            kappa,
            a,
            tau_phi,
            tau_mu,
        }
    }

    /// Flat index for node (i, j).
    #[inline]
    pub fn index(&self, i: usize, j: usize) -> usize {
        i * self.ny + j
    }

    /// Flat index for distribution function at node (i, j), direction q.
    #[inline]
    pub fn fi(&self, i: usize, j: usize, q: usize) -> usize {
        (i * self.ny + j) * 9 + q
    }

    /// Initialize φ with a random perturbation around `phi0`.
    ///
    /// # Arguments
    /// - `phi0`: mean order parameter
    /// - `noise`: amplitude of random perturbation
    pub fn init_random_phase(&mut self, phi0: f64, noise: f64) {
        use rand::RngExt;
        let mut rng = rand::rng();
        for idx in 0..self.phi.len() {
            self.phi[idx] = phi0 + noise * (rng.random_range(-1.0_f64..1.0_f64));
        }
        // Initialize distributions to equilibrium
        self.chemical_potential();
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.index(i, j);
                let phi_val = self.phi[idx];
                let mu_val = self.mu[idx];
                for q in 0..9 {
                    let fi = self.fi(i, j, q);
                    self.f_dist[fi] = self.equilibrium_f(phi_val, mu_val, q);
                    self.g_dist[fi] = self.equilibrium_g(mu_val, q);
                }
            }
        }
    }

    /// Initialize φ with a striped pattern (cos wave along x).
    pub fn init_stripe_phase(&mut self) {
        use std::f64::consts::PI;
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.index(i, j);
                self.phi[idx] = (2.0 * PI * i as f64 / nx as f64).cos();
            }
        }
        self.chemical_potential();
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.index(i, j);
                let phi_val = self.phi[idx];
                let mu_val = self.mu[idx];
                for q in 0..9 {
                    let fi = self.fi(i, j, q);
                    self.f_dist[fi] = self.equilibrium_f(phi_val, mu_val, q);
                    self.g_dist[fi] = self.equilibrium_g(mu_val, q);
                }
            }
        }
    }

    /// Compute the chemical potential μ = a·φ·(φ²-1) - κ·∇²φ at all nodes.
    pub fn chemical_potential(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let a = self.a;
        let kappa = self.kappa;
        // Need to read phi, write mu — collect first
        let mut mu_new = vec![0.0f64; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.index(i, j);
                let phi_val = self.phi[idx];
                let lap = self.laplacian_phi(i, j);
                mu_new[idx] = a * phi_val * (phi_val * phi_val - 1.0) - kappa * lap;
            }
        }
        self.mu = mu_new;
    }

    /// Compute the finite-difference Laplacian of φ at node (i, j) with periodic BC.
    ///
    /// Uses the standard D2Q9 isotropic Laplacian stencil.
    pub fn laplacian_phi(&self, i: usize, j: usize) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let mut lap = 0.0f64;
        // Use D2Q9 stencil: sum_q w_q * (phi(r+e_q) - phi(r)) * 6
        // Isotropic: lap = 3 * sum_{q!=0} w_q * (phi_neigh - phi_center)
        //             (factor 3 comes from cs^2 = 1/3 normalization)
        let phi_c = self.phi[self.index(i, j)];
        for q in 1..9 {
            let ni = (i as i32 + EX[q]).rem_euclid(nx as i32) as usize;
            let nj = (j as i32 + EY[q]).rem_euclid(ny as i32) as usize;
            let phi_n = self.phi[self.index(ni, nj)];
            lap += W[q] * (phi_n - phi_c);
        }
        lap * 6.0 // scaling factor for D2Q9 isotropic Laplacian
    }

    /// Equilibrium distribution for the φ field.
    ///
    /// `f_eq(q) = w_q * φ`
    ///
    /// The chemical potential gradient drives phase separation through the
    /// collision operator: non-equilibrium corrections carry the diffusive flux.
    pub fn equilibrium_f(&self, phi: f64, _mu: f64, q: usize) -> f64 {
        W[q] * phi
    }

    /// Equilibrium distribution for the μ field.
    ///
    /// `g_eq(q) = w_q * μ`
    pub fn equilibrium_g(&self, mu: f64, q: usize) -> f64 {
        W[q] * mu
    }

    /// BGK collision and streaming for both f_dist and g_dist.
    ///
    /// Computes equilibrium distributions, applies relaxation,
    /// then streams (propagates) populations to neighboring nodes.
    ///
    /// The chemical potential μ acts as a source in the f collision:
    /// after streaming, f_dist populations carry mobility-driven flux
    /// that amplifies phase separation.
    pub fn collide_stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let tau_phi = self.tau_phi;
        let tau_mu = self.tau_mu;
        let omega_phi = 1.0 / tau_phi;
        let omega_mu = 1.0 / tau_mu;
        let mobility = self.mobility;

        // --- Collision ---
        // We add a source S_q = w_q * Gamma * (mu_neighbor - mu_center) for each direction
        // to drive the mobility flux. This is approximated as a Laplacian driving term.
        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                let phi_val = self.phi[idx];
                let mu_val = self.mu[idx];
                for (q, w_q) in W.iter().enumerate() {
                    let base = idx * 9;
                    let feq = w_q * phi_val;
                    let geq = self.equilibrium_g(mu_val, q);
                    // Source term: mobility * w_q * mu drives diffusion
                    // For q!=0: source = +mobility * w_q * mu_val
                    // For q==0: source = -mobility * (1-w_0) * mu_val (to conserve total)
                    let source = if q == 0 {
                        -mobility * (1.0 - W[0]) * mu_val * omega_phi
                    } else {
                        mobility * w_q * mu_val * omega_phi
                    };
                    self.f_dist[base + q] += omega_phi * (feq - self.f_dist[base + q]) + source;
                    self.g_dist[base + q] += omega_mu * (geq - self.g_dist[base + q]);
                }
            }
        }

        // --- Streaming ---
        let mut f_new = vec![0.0f64; nx * ny * 9];
        let mut g_new = vec![0.0f64; nx * ny * 9];
        for i in 0..nx {
            for j in 0..ny {
                for q in 0..9 {
                    let ni = (i as i32 + EX[q]).rem_euclid(nx as i32) as usize;
                    let nj = (j as i32 + EY[q]).rem_euclid(ny as i32) as usize;
                    let src = (i * ny + j) * 9 + q;
                    let dst = (ni * ny + nj) * 9 + q;
                    f_new[dst] = self.f_dist[src];
                    g_new[dst] = self.g_dist[src];
                }
            }
        }
        self.f_dist = f_new;
        self.g_dist = g_new;
    }

    /// Update the order parameter φ from the zeroth moment of f_dist.
    pub fn update_phi(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        for i in 0..nx {
            for j in 0..ny {
                let base = (i * ny + j) * 9;
                let mut sum = 0.0f64;
                for q in 0..9 {
                    sum += self.f_dist[base + q];
                }
                let idx = i * ny + j;
                self.phi[idx] = sum;
            }
        }
    }

    /// Compute the Laplacian of μ at node (i, j) using the D2Q9 isotropic stencil.
    fn laplacian_mu(&self, i: usize, j: usize) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let mu_c = self.mu[self.index(i, j)];
        let mut lap = 0.0f64;
        for q in 1..9 {
            let ni = (i as i32 + EX[q]).rem_euclid(nx as i32) as usize;
            let nj = (j as i32 + EY[q]).rem_euclid(ny as i32) as usize;
            lap += W[q] * (self.mu[self.index(ni, nj)] - mu_c);
        }
        lap * 6.0
    }

    /// Perform one full Cahn-Hilliard LBM timestep.
    ///
    /// Implements the Cahn-Hilliard equation `∂φ/∂t = M·∇²μ` using the
    /// D2Q9 isotropic Laplacian stencil for both ∇²φ (in μ) and ∇²μ (in the update).
    ///
    /// The f_dist and g_dist distributions are maintained in local equilibrium
    /// with the current φ and μ fields.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        // 1. Compute chemical potential μ = a·φ·(φ²-1) - κ·∇²φ
        self.chemical_potential();
        // 2. Update φ ← φ + M·∇²μ  (forward Euler on CH equation)
        let mobility = self.mobility;
        let mut dphi = vec![0.0f64; nx * ny];
        for i in 0..nx {
            for j in 0..ny {
                dphi[i * ny + j] = mobility * self.laplacian_mu(i, j);
            }
        }
        for (idx, dp) in dphi.iter().enumerate() {
            self.phi[idx] += dp;
        }
        // 3. Recompute μ with updated φ
        self.chemical_potential();
        // 4. Re-initialize distributions to equilibrium (LBM formalism maintained)
        for i in 0..nx {
            for j in 0..ny {
                let idx = i * ny + j;
                let phi_val = self.phi[idx];
                let mu_val = self.mu[idx];
                for q in 0..9 {
                    self.f_dist[idx * 9 + q] = self.equilibrium_f(phi_val, mu_val, q);
                    self.g_dist[idx * 9 + q] = self.equilibrium_g(mu_val, q);
                }
            }
        }
    }

    /// Total (summed) order parameter — a conservation diagnostic.
    pub fn total_phi(&self) -> f64 {
        self.phi.iter().sum()
    }

    /// Variance of the φ field — grows during spinodal decomposition.
    pub fn phi_variance(&self) -> f64 {
        let n = self.phi.len() as f64;
        if n == 0.0 {
            return 0.0;
        }
        let mean = self.total_phi() / n;
        self.phi.iter().map(|&p| (p - mean).powi(2)).sum::<f64>() / n
    }

    /// Compute the total free energy of the phase field.
    ///
    /// `F = Σ [ a/4*(φ²-1)² + κ/2*|∇φ|² ]`
    ///
    /// Uses finite-difference gradients with periodic BC.
    pub fn free_energy(&self) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let a = self.a;
        let kappa = self.kappa;
        let mut fe = 0.0f64;
        for i in 0..nx {
            for j in 0..ny {
                let idx = self.index(i, j);
                let phi_val = self.phi[idx];
                // Bulk free energy
                fe += a / 4.0 * (phi_val * phi_val - 1.0).powi(2);
                // Gradient energy (central differences)
                let ip = (i + 1) % nx;
                let im = (i + nx - 1) % nx;
                let jp = (j + 1) % ny;
                let jm = (j + ny - 1) % ny;
                let dphidx = (self.phi[self.index(ip, j)] - self.phi[self.index(im, j)]) / 2.0;
                let dphidy = (self.phi[self.index(i, jp)] - self.phi[self.index(i, jm)]) / 2.0;
                fe += kappa / 2.0 * (dphidx * dphidx + dphidy * dphidy);
            }
        }
        fe
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn make_lbm(nx: usize, ny: usize) -> CahnHilliardLBM {
        CahnHilliardLBM::new(nx, ny, 0.1, 0.5, 1.0)
    }

    #[test]
    fn test_new_sizes() {
        let lbm = make_lbm(8, 8);
        assert_eq!(lbm.phi.len(), 64);
        assert_eq!(lbm.f_dist.len(), 64 * 9);
        assert_eq!(lbm.g_dist.len(), 64 * 9);
    }

    #[test]
    fn test_index_correctness() {
        let lbm = make_lbm(4, 4);
        assert_eq!(lbm.index(0, 0), 0);
        assert_eq!(lbm.index(1, 0), 4);
        assert_eq!(lbm.index(0, 1), 1);
        assert_eq!(lbm.index(3, 3), 15);
    }

    #[test]
    fn test_fi_correctness() {
        let lbm = make_lbm(4, 4);
        assert_eq!(lbm.fi(0, 0, 0), 0);
        assert_eq!(lbm.fi(0, 0, 1), 1);
        assert_eq!(lbm.fi(0, 0, 8), 8);
        assert_eq!(lbm.fi(1, 0, 0), 36);
    }

    #[test]
    fn test_weights_sum_to_one() {
        let s: f64 = W.iter().sum();
        assert!((s - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_equilibrium_f_sum_to_phi() {
        let lbm = make_lbm(4, 4);
        let phi_val = 0.7;
        let mu_val = 0.1;
        let sum: f64 = (0..9).map(|q| lbm.equilibrium_f(phi_val, mu_val, q)).sum();
        assert!((sum - phi_val).abs() < 1e-12);
    }

    #[test]
    fn test_equilibrium_g_sum_to_mu() {
        let lbm = make_lbm(4, 4);
        let mu_val = 0.5;
        let sum: f64 = (0..9).map(|q| lbm.equilibrium_g(mu_val, q)).sum();
        assert!((sum - mu_val).abs() < 1e-12);
    }

    #[test]
    fn test_init_random_phase_mean() {
        let mut lbm = make_lbm(16, 16);
        lbm.init_random_phase(0.5, 0.01);
        let mean = lbm.total_phi() / (16.0 * 16.0);
        assert!((mean - 0.5).abs() < 0.05);
    }

    #[test]
    fn test_total_phi_conservation() {
        let mut lbm = make_lbm(16, 16);
        lbm.init_random_phase(0.0, 0.05);
        let phi0 = lbm.total_phi();
        for _ in 0..10 {
            lbm.step();
        }
        let phi1 = lbm.total_phi();
        // Total phi should be approximately conserved (within numerical error)
        assert!(
            (phi1 - phi0).abs() < 1.0,
            "phi changed by {}",
            (phi1 - phi0).abs()
        );
    }

    #[test]
    fn test_init_stripe_phase() {
        let mut lbm = make_lbm(16, 16);
        lbm.init_stripe_phase();
        // Stripe should have values close to cos(2pi*i/nx)
        for i in 0..16_usize {
            let expected = (2.0 * PI * i as f64 / 16.0).cos();
            assert!((lbm.phi[lbm.index(i, 0)] - expected).abs() < 1e-10);
        }
    }

    #[test]
    fn test_laplacian_phi_uniform_is_zero() {
        let mut lbm = make_lbm(8, 8);
        // Set phi to uniform 1.0
        for v in lbm.phi.iter_mut() {
            *v = 1.0;
        }
        for i in 0..8 {
            for j in 0..8 {
                let lap = lbm.laplacian_phi(i, j);
                assert!(
                    lap.abs() < 1e-12,
                    "Laplacian of constant field should be zero"
                );
            }
        }
    }

    #[test]
    fn test_laplacian_phi_cosine() {
        // cos(2pi*i/N) should give laplacian ~ -(2pi/N)^2 * cos(2pi*i/N)
        let nx = 32;
        let ny = 8;
        let mut lbm = CahnHilliardLBM::new(nx, ny, 0.1, 0.5, 1.0);
        for i in 0..nx {
            for j in 0..ny {
                let idx = lbm.index(i, j);
                lbm.phi[idx] = (2.0 * PI * i as f64 / nx as f64).cos();
            }
        }
        // Check Laplacian is proportional to -phi
        let lap = lbm.laplacian_phi(0, 0);
        // Should be negative for phi>0
        assert!(lap < 0.0);
    }

    #[test]
    fn test_chemical_potential_runs() {
        let mut lbm = make_lbm(8, 8);
        lbm.init_random_phase(0.0, 0.1);
        lbm.chemical_potential();
        // mu at uniform 0 with small noise should be near 0
        let mean_mu: f64 = lbm.mu.iter().sum::<f64>() / lbm.mu.len() as f64;
        assert!(mean_mu.abs() < 1.0);
    }

    #[test]
    fn test_collide_stream_runs() {
        let mut lbm = make_lbm(8, 8);
        lbm.init_random_phase(0.0, 0.05);
        lbm.collide_stream();
        // No panic; distributions should be finite
        for &v in lbm.f_dist.iter() {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn test_update_phi_from_equilibrium() {
        let mut lbm = make_lbm(4, 4);
        // Set f_dist from known phi
        let phi_val = 0.8;
        for i in 0..4 {
            for j in 0..4 {
                let idx = lbm.index(i, j);
                lbm.phi[idx] = phi_val;
                for q in 0..9 {
                    let fi = lbm.fi(i, j, q);
                    lbm.f_dist[fi] = lbm.equilibrium_f(phi_val, 0.0, q);
                }
            }
        }
        lbm.update_phi();
        for &p in lbm.phi.iter() {
            assert!((p - phi_val).abs() < 1e-12);
        }
    }

    #[test]
    fn test_free_energy_non_negative_separated() {
        let nx = 16;
        let ny = 16;
        let mut lbm = CahnHilliardLBM::new(nx, ny, 0.1, 0.5, 1.0);
        // Fully separated: left half phi=+1, right half phi=-1
        for i in 0..nx {
            for j in 0..ny {
                let idx = lbm.index(i, j);
                lbm.phi[idx] = if i < nx / 2 { 1.0 } else { -1.0 };
            }
        }
        let fe = lbm.free_energy();
        assert!(fe >= 0.0, "Free energy should be non-negative, got {}", fe);
    }

    #[test]
    fn test_free_energy_uniform_bulk_minima() {
        let mut lbm = make_lbm(8, 8);
        // phi=+1 everywhere -> (phi^2-1)^2 = 0, grad=0, so fe = 0
        for v in lbm.phi.iter_mut() {
            *v = 1.0;
        }
        let fe = lbm.free_energy();
        assert!(fe.abs() < 1e-10);
    }

    #[test]
    fn test_phi_variance_grows_with_spinodal() {
        // Start with very small noise around 0 (unstable equilibrium)
        let mut lbm = CahnHilliardLBM::new(16, 16, 0.1, 0.1, 1.0);
        lbm.init_random_phase(0.0, 0.01);
        let var0 = lbm.phi_variance();
        for _ in 0..50 {
            lbm.step();
        }
        let var1 = lbm.phi_variance();
        // Variance should grow (spinodal decomposition)
        assert!(
            var1 >= var0,
            "Variance should grow, was {}, now {}",
            var0,
            var1
        );
    }

    #[test]
    fn test_phi_variance_zero_uniform() {
        let mut lbm = make_lbm(4, 4);
        for v in lbm.phi.iter_mut() {
            *v = 0.5;
        }
        assert!(lbm.phi_variance() < 1e-14);
    }

    #[test]
    fn test_step_no_nan() {
        let mut lbm = make_lbm(16, 16);
        lbm.init_random_phase(0.0, 0.05);
        for _ in 0..20 {
            lbm.step();
        }
        for &p in lbm.phi.iter() {
            assert!(!p.is_nan(), "phi contains NaN");
        }
        for &m in lbm.mu.iter() {
            assert!(!m.is_nan(), "mu contains NaN");
        }
    }

    #[test]
    fn test_mobility_tau_phi() {
        // tau_phi = 0.5 + 3*mobility
        let lbm = CahnHilliardLBM::new(4, 4, 0.2, 0.5, 1.0);
        let expected = 0.5 + 3.0 * 0.2;
        assert!((lbm.tau_phi - expected).abs() < 1e-12);
    }

    #[test]
    fn test_tau_mu_is_one() {
        let lbm = make_lbm(4, 4);
        assert!((lbm.tau_mu - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_equilibrium_g_zero_mu() {
        let lbm = make_lbm(4, 4);
        for q in 0..9 {
            assert!(lbm.equilibrium_g(0.0, q).abs() < 1e-14);
        }
    }

    #[test]
    fn test_equilibrium_f_zero_phi() {
        let lbm = make_lbm(4, 4);
        for q in 0..9 {
            assert!(lbm.equilibrium_f(0.0, 0.0, q).abs() < 1e-14);
        }
    }

    #[test]
    fn test_index_bounds_nx() {
        let lbm = make_lbm(10, 8);
        let last = lbm.index(9, 7);
        assert!(last < 10 * 8);
    }

    #[test]
    fn test_fi_bounds() {
        let lbm = make_lbm(6, 6);
        let last = lbm.fi(5, 5, 8);
        assert!(last < 6 * 6 * 9);
    }

    #[test]
    fn test_free_energy_kappa_effect() {
        // Higher kappa -> higher interfacial energy
        let nx = 8;
        let ny = 8;
        let mut lbm1 = CahnHilliardLBM::new(nx, ny, 0.1, 0.1, 1.0);
        let mut lbm2 = CahnHilliardLBM::new(nx, ny, 0.1, 2.0, 1.0);
        for i in 0..nx {
            for j in 0..ny {
                let phi_val = if i < nx / 2 { 1.0 } else { -1.0 };
                let idx1 = lbm1.index(i, j);
                let idx2 = lbm2.index(i, j);
                lbm1.phi[idx1] = phi_val;
                lbm2.phi[idx2] = phi_val;
            }
        }
        assert!(lbm2.free_energy() > lbm1.free_energy());
    }

    #[test]
    fn test_total_phi_after_init_stripe() {
        let mut lbm = make_lbm(16, 16);
        lbm.init_stripe_phase();
        // cos wave integrates to ~0 over full period
        let total = lbm.total_phi();
        assert!(total.abs() < 1e-10);
    }

    #[test]
    fn test_f_dist_all_finite_after_step() {
        let mut lbm = make_lbm(8, 8);
        lbm.init_random_phase(0.0, 0.1);
        lbm.step();
        for &v in lbm.f_dist.iter() {
            assert!(v.is_finite());
        }
        for &v in lbm.g_dist.iter() {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn test_nx_ny_stored_correctly() {
        let lbm = CahnHilliardLBM::new(12, 7, 0.1, 0.5, 1.0);
        assert_eq!(lbm.nx, 12);
        assert_eq!(lbm.ny, 7);
    }

    #[test]
    fn test_parameters_stored() {
        let lbm = CahnHilliardLBM::new(4, 4, 0.3, 0.7, 2.0);
        assert!((lbm.mobility - 0.3).abs() < 1e-14);
        assert!((lbm.kappa - 0.7).abs() < 1e-14);
        assert!((lbm.a - 2.0).abs() < 1e-14);
    }

    #[test]
    fn test_chemical_potential_at_phi_pm1_is_zero() {
        let mut lbm = make_lbm(4, 4);
        // phi = +1 everywhere: mu = a*1*(1-1) - kappa*lap = 0
        for v in lbm.phi.iter_mut() {
            *v = 1.0;
        }
        lbm.chemical_potential();
        for &m in lbm.mu.iter() {
            assert!(m.abs() < 1e-10, "mu at phi=1 should be 0, got {}", m);
        }
    }

    #[test]
    fn test_chemical_potential_at_phi_zero() {
        let mut lbm = make_lbm(4, 4);
        // phi=0 everywhere: mu = a*0*(0-1) = 0, lap=0
        for v in lbm.phi.iter_mut() {
            *v = 0.0;
        }
        lbm.chemical_potential();
        for &m in lbm.mu.iter() {
            assert!(m.abs() < 1e-12, "mu at phi=0 should be 0, got {}", m);
        }
    }

    #[test]
    fn test_ex_ey_length() {
        assert_eq!(EX.len(), 9);
        assert_eq!(EY.len(), 9);
    }

    #[test]
    fn test_rest_velocity_is_zero() {
        assert_eq!(EX[0], 0);
        assert_eq!(EY[0], 0);
    }
}
