#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Chain models: PolymerChain, WormlikeChain, FreelyJointedChain, FENEBond,
//! WLCSampler, KratkyPorodChain.

use super::{dist3, dot3, fene_force, mag3, norm3, sub3, wlc_extension};

/// Bead-spring polymer chain with Lennard-Jones and FENE bonds.
///
/// Models a linear polymer using connected Lennard-Jones beads
/// and FENE (Finitely Extensible Nonlinear Elastic) springs.
#[derive(Debug, Clone)]
pub struct PolymerChain {
    /// Number of beads.
    pub n_beads: usize,
    /// Bead positions \[x, y, z\].
    pub positions: Vec<[f64; 3]>,
    /// Bead velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Bead masses.
    pub masses: Vec<f64>,
    /// FENE maximum extension R0.
    pub fene_r0: f64,
    /// FENE spring constant.
    pub fene_k: f64,
    /// LJ epsilon.
    pub lj_eps: f64,
    /// LJ sigma.
    pub lj_sigma: f64,
    /// Persistence length in bond units.
    pub persistence_length: f64,
    /// Temperature kT.
    pub kt: f64,
    /// Time step.
    pub dt: f64,
}

impl PolymerChain {
    /// Create a new linear polymer chain.
    pub fn new(
        n_beads: usize,
        fene_r0: f64,
        fene_k: f64,
        lj_eps: f64,
        lj_sigma: f64,
        persistence_length: f64,
        kt: f64,
        dt: f64,
    ) -> Self {
        // Initialize stretched chain along x-axis
        let positions: Vec<[f64; 3]> = (0..n_beads)
            .map(|i| [i as f64 * lj_sigma, 0.0, 0.0])
            .collect();
        let velocities = vec![[0.0f64; 3]; n_beads];
        let masses = vec![1.0f64; n_beads];
        Self {
            n_beads,
            positions,
            velocities,
            masses,
            fene_r0,
            fene_k,
            lj_eps,
            lj_sigma,
            persistence_length,
            kt,
            dt,
        }
    }

    /// Compute FENE bond force between adjacent beads i and i+1.
    pub fn fene_force_pair(&self, i: usize) -> [f64; 3] {
        if i + 1 >= self.n_beads {
            return [0.0; 3];
        }
        let r = sub3(self.positions[i + 1], self.positions[i]);
        let rmag = mag3(r);
        let f = fene_force(rmag, self.fene_k, self.fene_r0);
        let rhat = norm3(r);
        [f * rhat[0], f * rhat[1], f * rhat[2]]
    }

    /// Compute LJ force between beads i and j.
    pub fn lj_force_pair(&self, i: usize, j: usize) -> [f64; 3] {
        if i == j {
            return [0.0; 3];
        }
        let r = sub3(self.positions[j], self.positions[i]);
        let rmag = mag3(r).max(1e-10);
        let eps = self.lj_eps;
        let sig = self.lj_sigma;
        let s6 = (sig / rmag).powi(6);
        let fmag = 24.0 * eps / rmag * (2.0 * s6 * s6 - s6);
        let rhat = norm3(r);
        [fmag * rhat[0], fmag * rhat[1], fmag * rhat[2]]
    }

    /// Total force on bead i.
    pub fn total_force(&self, i: usize) -> [f64; 3] {
        let mut f = [0.0f64; 3];
        // FENE bond to i-1
        if i > 0 {
            let r = sub3(self.positions[i], self.positions[i - 1]);
            let rmag = mag3(r).max(1e-10);
            let fmag = fene_force(rmag, self.fene_k, self.fene_r0);
            let rhat = norm3(r);
            f[0] -= fmag * rhat[0];
            f[1] -= fmag * rhat[1];
            f[2] -= fmag * rhat[2];
        }
        // FENE bond to i+1
        if i + 1 < self.n_beads {
            let r = sub3(self.positions[i + 1], self.positions[i]);
            let rmag = mag3(r).max(1e-10);
            let fmag = fene_force(rmag, self.fene_k, self.fene_r0);
            let rhat = norm3(r);
            f[0] += fmag * rhat[0];
            f[1] += fmag * rhat[1];
            f[2] += fmag * rhat[2];
        }
        // LJ with all non-bonded beads
        for j in 0..self.n_beads {
            if j.abs_diff(i) > 1 {
                let flj = self.lj_force_pair(i, j);
                f[0] += flj[0];
                f[1] += flj[1];
                f[2] += flj[2];
            }
        }
        f
    }

    /// Velocity Verlet step.
    pub fn step(&mut self) {
        let dt = self.dt;
        let forces: Vec<[f64; 3]> = (0..self.n_beads).map(|i| self.total_force(i)).collect();
        for i in 0..self.n_beads {
            let m = self.masses[i];
            self.velocities[i][0] += forces[i][0] / m * dt;
            self.velocities[i][1] += forces[i][1] / m * dt;
            self.velocities[i][2] += forces[i][2] / m * dt;
            self.positions[i][0] += self.velocities[i][0] * dt;
            self.positions[i][1] += self.velocities[i][1] * dt;
            self.positions[i][2] += self.velocities[i][2] * dt;
        }
    }

    /// Compute end-to-end distance.
    pub fn end_to_end(&self) -> f64 {
        if self.n_beads < 2 {
            return 0.0;
        }
        dist3(self.positions[0], self.positions[self.n_beads - 1])
    }

    /// Compute radius of gyration.
    pub fn radius_of_gyration(&self) -> f64 {
        let n = self.n_beads as f64;
        let com: [f64; 3] = self.positions.iter().fold([0.0; 3], |acc, p| {
            [acc[0] + p[0] / n, acc[1] + p[1] / n, acc[2] + p[2] / n]
        });
        let rg2: f64 = self
            .positions
            .iter()
            .map(|&p| {
                let d = sub3(p, com);
                dot3(d, d)
            })
            .sum::<f64>()
            / n;
        rg2.sqrt()
    }

    /// Mean square end-to-end distance (ideal chain: <R^2> = N*b^2).
    pub fn ideal_r2(&self) -> f64 {
        self.n_beads as f64 * self.lj_sigma * self.lj_sigma
    }
}

/// Worm-like chain (WLC) force extension model.
///
/// Implements F = kT/p * (1/(4(1-r/L)²) - 1/4 + r/L)
#[derive(Debug, Clone)]
pub struct WormlikeChain {
    /// Contour length L in nm.
    pub contour_length: f64,
    /// Persistence length p in nm.
    pub persistence_length: f64,
    /// Thermal energy kT in pN·nm.
    pub kt: f64,
}

impl WormlikeChain {
    /// Create a new WLC model.
    pub fn new(contour_length: f64, persistence_length: f64, kt: f64) -> Self {
        Self {
            contour_length,
            persistence_length,
            kt,
        }
    }

    /// Compute WLC force at extension r (in nm).
    pub fn force(&self, r: f64) -> f64 {
        wlc_extension(r, self.contour_length, self.persistence_length, self.kt)
    }

    /// Compute WLC elastic energy by integration.
    pub fn energy(&self, r: f64) -> f64 {
        let n_steps = 100usize;
        let dr = r / n_steps as f64;
        (0..n_steps).map(|i| self.force(i as f64 * dr) * dr).sum()
    }

    /// Relative extension r/L.
    pub fn relative_extension(&self, r: f64) -> f64 {
        (r / self.contour_length).clamp(0.0, 0.9999)
    }

    /// Stiffness dF/dr at extension r.
    pub fn stiffness(&self, r: f64) -> f64 {
        let x = self.relative_extension(r);
        let p = self.persistence_length;
        let l = self.contour_length;
        self.kt / p * (0.5 / ((1.0 - x) * (1.0 - x) * (1.0 - x)) + 1.0 / l) * 1.0 / l
    }
}

// ─── WLC extended model ───────────────────────────────────────────────────────

/// Worm-like chain Monte Carlo sampler.
///
/// Generates WLC conformations using discrete bond angle model
/// with bending rigidity κ = kT * persistence_length.
#[derive(Debug, Clone)]
pub struct WLCSampler {
    /// Number of segments.
    pub n_segments: usize,
    /// Segment length b.
    pub b: f64,
    /// Persistence length l_p.
    pub persistence_length: f64,
    /// Thermal energy kT.
    pub kt: f64,
}

impl WLCSampler {
    /// Create a new WLC sampler.
    pub fn new(n_segments: usize, b: f64, persistence_length: f64, kt: f64) -> Self {
        WLCSampler {
            n_segments,
            b,
            persistence_length,
            kt,
        }
    }

    /// Bending stiffness κ = kT * l_p / b.
    pub fn bending_stiffness(&self) -> f64 {
        self.kt * self.persistence_length / self.b
    }

    /// Mean cosine of bond angle: <cos θ> = exp(-b / l_p).
    pub fn mean_cos_angle(&self) -> f64 {
        (-self.b / self.persistence_length).exp()
    }

    /// End-to-end distance squared in WLC limit.
    ///
    /// `R²` = 2 * l_p * L * \[1 - (l_p/L)(1 - exp(-L/l_p))\].
    pub fn mean_r2(&self) -> f64 {
        let l = self.n_segments as f64 * self.b;
        let lp = self.persistence_length;
        2.0 * lp * l * (1.0 - (lp / l) * (1.0 - (-l / lp).exp()))
    }

    /// Compute tangent-tangent correlation <t(0)·t(s)> = exp(-s/l_p).
    pub fn tangent_correlation(&self, s: f64) -> f64 {
        (-s / self.persistence_length).exp()
    }

    /// Persistence length from given bond correlation data.
    ///
    /// Fits l_p = -b / ln(<cos θ>).
    pub fn fit_persistence_length(b: f64, mean_cos_theta: f64) -> f64 {
        if mean_cos_theta <= 0.0 || mean_cos_theta >= 1.0 {
            return f64::INFINITY;
        }
        -b / mean_cos_theta.ln()
    }

    /// Radius of gyration for WLC.
    ///
    /// Rg² = (l_p * L / 3)\[1 - (3 l_p / L)(1 - (l_p/L)(1 - exp(-L/l_p)))\].
    pub fn radius_of_gyration_sq(&self) -> f64 {
        let l = self.n_segments as f64 * self.b;
        let lp = self.persistence_length;
        let ratio = lp / l;
        (lp * l / 3.0)
            * (1.0 - 3.0 * ratio + 6.0 * ratio * ratio
                - 6.0 * ratio * ratio * ratio * (1.0 - (-1.0 / ratio).exp()))
    }
}

// ─── Freely-jointed chain ─────────────────────────────────────────────────────

/// Freely-jointed chain (FJC) model with Kuhn segments.
///
/// Each segment has length b (Kuhn length) and is independently oriented.
/// Stretching is described by the Langevin function.
#[derive(Debug, Clone)]
pub struct FreelyJointedChain {
    /// Number of Kuhn segments N.
    pub n_segments: usize,
    /// Kuhn segment length b.
    pub b: f64,
    /// Thermal energy kT.
    pub kt: f64,
}

impl FreelyJointedChain {
    /// Create a new FJC model.
    pub fn new(n_segments: usize, b: f64, kt: f64) -> Self {
        FreelyJointedChain { n_segments, b, kt }
    }

    /// Contour length L = N * b.
    pub fn contour_length(&self) -> f64 {
        self.n_segments as f64 * self.b
    }

    /// Mean-square end-to-end distance `R²` = N * b².
    pub fn mean_r2(&self) -> f64 {
        self.n_segments as f64 * self.b * self.b
    }

    /// Root-mean-square end-to-end distance.
    pub fn rms_end_to_end(&self) -> f64 {
        self.mean_r2().sqrt()
    }

    /// Langevin function L(x) = coth(x) - 1/x.
    pub fn langevin(x: f64) -> f64 {
        if x.abs() < 1e-6 {
            return x / 3.0;
        }
        1.0 / x.tanh() - 1.0 / x
    }

    /// Inverse Langevin via Padé approximation: L⁻¹(y) ≈ y(3 - y²) / (1 - y²).
    pub fn inverse_langevin(y: f64) -> f64 {
        let y = y.clamp(-0.9999, 0.9999);
        y * (3.0 - y * y) / (1.0 - y * y)
    }

    /// Extension `r` = N * b * L(f * b / kT) under force f.
    pub fn mean_extension(&self, force: f64) -> f64 {
        let x = force * self.b / self.kt;
        self.contour_length() * Self::langevin(x)
    }

    /// Force for given extension r using inverse Langevin.
    pub fn force_extension(&self, r: f64) -> f64 {
        let y = r / self.contour_length();
        let y_clamped = y.clamp(-0.9999, 0.9999);
        self.kt / self.b * Self::inverse_langevin(y_clamped)
    }

    /// Entropic spring constant k_eff = 3 kT / (N b²) (Gaussian limit).
    pub fn entropic_spring_constant(&self) -> f64 {
        3.0 * self.kt / (self.n_segments as f64 * self.b * self.b)
    }

    /// Flory exponent ν = 0.5 (ideal chain, no excluded volume).
    pub fn flory_exponent() -> f64 {
        0.5
    }
}

// ─── FENE bond extended model ─────────────────────────────────────────────────

/// FENE (Finitely Extensible Nonlinear Elastic) bond extended model.
///
/// U(r) = -K R₀²/2 * ln(1 - (r/R₀)²)
#[derive(Debug, Clone)]
pub struct FENEBond {
    /// Spring constant K.
    pub k: f64,
    /// Maximum extension R₀.
    pub r0: f64,
    /// Equilibrium bond length r_eq (for reference).
    pub r_eq: f64,
}

impl FENEBond {
    /// Create a new FENE bond.
    pub fn new(k: f64, r0: f64, r_eq: f64) -> Self {
        FENEBond { k, r0, r_eq }
    }

    /// FENE potential energy U(r).
    pub fn potential(&self, r: f64) -> f64 {
        let x = r / self.r0;
        if x >= 1.0 {
            return f64::INFINITY;
        }
        -0.5 * self.k * self.r0 * self.r0 * (1.0 - x * x).ln()
    }

    /// FENE force magnitude F(r) = K*r / (1 - (r/R₀)²).
    pub fn force_magnitude(&self, r: f64) -> f64 {
        let x = r / self.r0;
        if x >= 1.0 {
            return f64::INFINITY;
        }
        self.k * r / (1.0 - x * x)
    }

    /// Maximum force at r → R₀.
    pub fn max_force_approx(&self, epsilon: f64) -> f64 {
        self.force_magnitude(self.r0 * (1.0 - epsilon))
    }

    /// Taylor expansion of potential for small r: U ≈ (1/2) K r².
    pub fn harmonic_approx_stiffness(&self) -> f64 {
        self.k
    }

    /// Combined FENE + WCA (truncated LJ) potential (Kremer-Grest model).
    ///
    /// Returns total potential for separation r with LJ parameters eps, sig.
    pub fn fene_wca_potential(&self, r: f64, eps: f64, sig: f64) -> f64 {
        let r_cut = 2.0_f64.powf(1.0 / 6.0) * sig;
        let u_fene = if r < self.r0 {
            self.potential(r)
        } else {
            f64::INFINITY
        };
        let u_wca = if r < r_cut {
            let s6 = (sig / r).powi(6);
            4.0 * eps * (s6 * s6 - s6) + eps
        } else {
            0.0
        };
        u_fene + u_wca
    }
}

// ─── Kratky-Porod worm-like chain (extended) ─────────────────────────────────

/// Kratky-Porod worm-like chain model with exact partition function.
///
/// The Kratky-Porod chain is the continuous limit of a chain with bending
/// rigidity κ = kT * l_p. Used for semi-flexible polymers (DNA, actin, etc.)
///
/// References:
/// - Kratky, O. & Porod, G. (1949). *Rec. Trav. Chim.* 68, 1106.
/// - Marko, J.F. & Siggia, E.D. (1995). *Macromolecules* 28, 8759.
#[derive(Debug, Clone)]
pub struct KratkyPorodChain {
    /// Contour length L.
    pub contour_length: f64,
    /// Persistence length l_p.
    pub lp: f64,
    /// Thermal energy kT.
    pub kt: f64,
    /// Number of segments for discretization.
    pub n_seg: usize,
}

impl KratkyPorodChain {
    /// Create a new Kratky-Porod chain.
    pub fn new(contour_length: f64, lp: f64, kt: f64, n_seg: usize) -> Self {
        Self {
            contour_length,
            lp,
            kt,
            n_seg,
        }
    }

    /// Bending stiffness κ = kT * l_p.
    pub fn bending_stiffness(&self) -> f64 {
        self.kt * self.lp
    }

    /// Reduced length u = L / (2 l_p).
    pub fn reduced_length(&self) -> f64 {
        self.contour_length / (2.0 * self.lp)
    }

    /// Mean-square end-to-end distance.
    ///
    /// `R²` = 2 l_p L \[1 - (l_p/L)(1 - exp(-L/l_p))\]
    pub fn mean_r2(&self) -> f64 {
        let l = self.contour_length;
        let lp = self.lp;
        2.0 * lp * l * (1.0 - (lp / l) * (1.0 - (-l / lp).exp()))
    }

    /// Radius of gyration (Rg²).
    ///
    /// Rg² = l_p L / 3 * \[1 - (3l_p/L)(1 - (3l_p²/L²)(1 - exp(-L/l_p)))\]
    pub fn radius_of_gyration_sq(&self) -> f64 {
        let l = self.contour_length;
        let lp = self.lp;
        let u = l / lp;
        lp * l / 3.0 * (1.0 - 3.0 / u + 6.0 / (u * u) - 6.0 / (u * u * u) * (1.0 - (-u).exp()))
    }

    /// WLC force extension (Marko-Siggia interpolation formula).
    ///
    /// F = kT/l_p * \[1/(4(1-x)²) - 1/4 + x\]  where x = r/L.
    pub fn force_extension(&self, r: f64) -> f64 {
        let x = (r / self.contour_length).clamp(0.0, 0.9999);
        self.kt / self.lp * (0.25 / ((1.0 - x) * (1.0 - x)) - 0.25 + x)
    }

    /// Extension at given force (numerical inversion using bisection).
    pub fn extension_at_force(&self, force: f64) -> f64 {
        let l = self.contour_length;
        let mut lo = 0.0f64;
        let mut hi = 0.9999 * l;
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            if self.force_extension(mid) < force {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Tangent-tangent correlation: <t(s)·t(0)> = exp(-s/l_p).
    pub fn tangent_correlation(&self, s: f64) -> f64 {
        (-s / self.lp.max(1e-30)).exp()
    }

    /// Persistence ratio l_p / L.
    pub fn persistence_ratio(&self) -> f64 {
        self.lp / self.contour_length.max(1e-30)
    }

    /// Classify regime: rod (l_p >> L), semi-flexible (l_p ~ L), or coil (l_p << L).
    pub fn regime(&self) -> &'static str {
        let ratio = self.persistence_ratio();
        if ratio > 10.0 {
            "rod"
        } else if ratio > 0.1 {
            "semi-flexible"
        } else {
            "coil"
        }
    }

    /// Effective spring constant (entropic) at small extension.
    ///
    /// k_eff = 3 kT / (2 l_p L) (WLC linear response).
    pub fn spring_constant_linear(&self) -> f64 {
        3.0 * self.kt / (2.0 * self.lp * self.contour_length)
    }
}
