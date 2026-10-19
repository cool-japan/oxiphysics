// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cascaded (central-moment) lattice Boltzmann collision for D2Q9.
//!
//! The collision is carried out in the frame that moves with the fluid.  The
//! nine populations are mapped to the nine central moments
//!
//! ```text
//! κ_pq = Σ_i f_i (c_ix − u_x)^p (c_iy − u_y)^q ,      p, q ∈ {0, 1, 2},
//! ```
//!
//! every moment relaxes toward its equilibrium value, and the populations are
//! rebuilt **exactly** from the relaxed moments:
//!
//! 1. the relaxed central moments are shifted back to raw moments
//!    `m_pq = Σ_i f_i c_ix^p c_iy^q` with the binomial expansion of
//!    `c = (c − u) + u`,
//! 2. the raw moments are inverted with the tensor-product inverse of the
//!    one-dimensional `{−1, 0, 1}` moment problem,
//!    `f(a, b) = Σ_pq m_pq g_p(a) g_q(b)` with
//!    `g_0(c) = 1 − c²`, `g_1(c) = c / 2`, `g_2(c) = (3c² − 2) / 2`.
//!
//! Both steps are exact linear maps (the nine monomials `x^p y^q` form a basis
//! of the functions on the D2Q9 stencil), so the collision conserves mass and
//! momentum to rounding error and never needs clipping or renormalisation.
//!
//! The equilibrium central moments are those of the Maxwell–Boltzmann
//! distribution, `κ_eq = (ρ, 0, 0, ρ/3, ρ/3, 0, 0, 0, ρ/9)`.  The matching
//! populations ([`CascadedCollision2D::equilibrium`]) are the exact fixed point
//! of the collision.  They differ from the second-order polynomial equilibrium
//! of the other D2Q9 operators of this crate only by lattice-aliasing terms of
//! order `|u|³` (the central moments `κ21`, `κ12`, `κ22` of the polynomial
//! equilibrium are `−ρ u_x² u_y`, `−ρ u_x u_y²` and `ρ/9 + 3ρ u_x² u_y²`).
//!
//! Reference: Geier, Greiner & Korvink (2006), Phys. Rev. E 73, 066705.

/// D2Q9 lattice velocity x-components, ordered like the rest of `mrt.rs`.
const CX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
/// D2Q9 lattice velocity y-components.
const CY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// Orders `(p, q)` of the nine moments `[κ00, κ10, κ01, κ20, κ02, κ11, κ21, κ12, κ22]`.
const ORDERS: [(usize, usize); 9] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (2, 0),
    (0, 2),
    (1, 1),
    (2, 1),
    (1, 2),
    (2, 2),
];

/// Dual basis of the monomials `1, c, c²` on the one-dimensional stencil
/// `c ∈ {−1, 0, 1}`: the coefficient of the raw moment `m_order` in the
/// population with velocity component `c`.
#[inline]
fn dual_basis(order: usize, c: f64) -> f64 {
    match order {
        0 => 1.0 - c * c,
        1 => 0.5 * c,
        _ => 0.5 * (3.0 * c * c - 2.0),
    }
}

/// Cascaded Lattice Boltzmann (CLBM) collision for D2Q9.
///
/// Performs relaxation in the frame of central moments (co-moving with fluid),
/// which provides enhanced Galilean invariance and improved stability at high
/// Mach/Reynolds numbers.  See the [module documentation](self) for the exact
/// moment transforms.
///
/// Reference: Geier, Greiner & Korvink (2006), Phys. Rev. E 73, 066705.
#[derive(Debug, Clone)]
pub struct CascadedCollision2D {
    /// Kinematic viscosity in lattice units the operator was configured for
    /// (`rates[5] = 1 / (3 nu + 1/2)`).
    pub nu: f64,
    /// Relaxation rates of the central-moment modes, ordered like the array
    /// returned by [`CascadedCollision2D::central_moments`]:
    ///
    /// | index | mode | role |
    /// |-------|------|------|
    /// | 0, 1, 2 | `κ00`, `κ10`, `κ01` | mass and momentum, conserved: rate 0 |
    /// | 3 | `κ20 + κ02` | trace of the second moment: bulk viscosity |
    /// | 4 | `κ20 − κ02` | normal-stress difference: shear viscosity |
    /// | 5 | `κ11` | shear stress: shear viscosity |
    /// | 6, 7 | `κ21`, `κ12` | third order, non-hydrodynamic |
    /// | 8 | `κ22` | fourth order, non-hydrodynamic |
    ///
    /// Indices 3 and 4 act on the trace and on the deviator of the second-order
    /// block (not on `κ20` and `κ02` separately) so that the viscosity stays
    /// isotropic.  Stable rates lie in `(0, 2)`; the shear rates are
    /// `1 / (3ν + 1/2)`.  A non-zero rate at index 0–2 relaxes the density and
    /// the momentum toward `(ρ, ρ u)` of the arguments of
    /// [`CascadedCollision2D::collide`], which breaks conservation.
    pub rates: [f64; 9],
}

impl CascadedCollision2D {
    /// Construct CLBM from kinematic viscosity.
    ///
    /// The shear modes (`κ20 − κ02`, `κ11`) relax at `s_ν = 1 / (3ν + 1/2)`,
    /// the bulk mode (`κ20 + κ02`) and the non-hydrodynamic modes
    /// (`κ21`, `κ12`, `κ22`) at 1, and the conserved modes at 0.
    pub fn new(nu: f64) -> Self {
        Self::with_bulk_rate(nu, 1.0)
    }

    /// Construct CLBM from kinematic viscosity with an explicit relaxation
    /// rate `bulk_rate` of the bulk mode `κ20 + κ02`.
    pub fn with_bulk_rate(nu: f64, bulk_rate: f64) -> Self {
        let s_nu = 1.0 / (3.0 * nu + 0.5);
        // Conserved modes relax at 0, ghost modes at 1.
        let rates = [0.0, 0.0, 0.0, bulk_rate, s_nu, s_nu, 1.0, 1.0, 1.0];
        Self { nu, rates }
    }

    /// Kinematic viscosity this operator was configured for.
    pub fn viscosity(&self) -> f64 {
        self.nu
    }

    /// Compute central moments of order (p,q): κ_{pq} = Σ fᵢ (cix-ux)^p (ciy-uy)^q
    ///
    /// The result is ordered `[κ00, κ10, κ01, κ20, κ02, κ11, κ21, κ12, κ22]`.
    /// With `ux = uy = 0` these are the raw moments of `f`.
    pub fn central_moments(f: &[f64; 9], ux: f64, uy: f64) -> [f64; 9] {
        let mut kappa = [0.0f64; 9];
        for i in 0..9 {
            let dx = CX[i] - ux;
            let dy = CY[i] - uy;
            kappa[0] += f[i]; // κ00
            kappa[1] += f[i] * dx; // κ10
            kappa[2] += f[i] * dy; // κ01
            kappa[3] += f[i] * dx * dx; // κ20
            kappa[4] += f[i] * dy * dy; // κ02
            kappa[5] += f[i] * dx * dy; // κ11
            kappa[6] += f[i] * dx * dx * dy; // κ21
            kappa[7] += f[i] * dx * dy * dy; // κ12
            kappa[8] += f[i] * dx * dx * dy * dy; // κ22
        }
        kappa
    }

    /// Equilibrium central moments for D2Q9.
    fn equilibrium_central_moments(rho: f64) -> [f64; 9] {
        // In the central-moment frame the equilibrium takes simple form:
        // κ00=ρ, κ10=κ01=0, κ20=κ02=ρ/3, κ11=0, κ21=κ12=0, κ22=ρ/9
        [
            rho,
            0.0,
            0.0,
            rho / 3.0,
            rho / 3.0,
            0.0,
            0.0,
            0.0,
            rho / 9.0,
        ]
    }

    /// Raw moments `[m00, m10, m01, m20, m02, m11, m21, m12, m22]` of a
    /// distribution from its central moments about the velocity `(ux, uy)`
    /// (binomial expansion of `c = (c − u) + u`).
    fn raw_from_central(kappa: &[f64; 9], ux: f64, uy: f64) -> [f64; 9] {
        let [k00, k10, k01, k20, k02, k11, k21, k12, k22] = *kappa;
        let (u, v) = (ux, uy);
        [
            k00,
            k10 + u * k00,
            k01 + v * k00,
            k20 + 2.0 * u * k10 + u * u * k00,
            k02 + 2.0 * v * k01 + v * v * k00,
            k11 + u * k01 + v * k10 + u * v * k00,
            k21 + 2.0 * u * k11 + v * k20 + u * u * k01 + 2.0 * u * v * k10 + u * u * v * k00,
            k12 + u * k02 + 2.0 * v * k11 + v * v * k10 + 2.0 * u * v * k01 + u * v * v * k00,
            k22 + 2.0 * u * k12
                + 2.0 * v * k21
                + u * u * k02
                + v * v * k20
                + 4.0 * u * v * k11
                + 2.0 * u * v * v * k10
                + 2.0 * u * u * v * k01
                + u * u * v * v * k00,
        ]
    }

    /// Populations from the nine raw moments `[m00, m10, m01, m20, m02, m11, m21, m12, m22]`.
    ///
    /// Exact inverse of the raw-moment map `m_pq = Σ_i f_i c_ix^p c_iy^q`, e.g.
    /// `f0 = m00 − m20 − m02 + m22` and `f1 = (m10 + m20 − m12 − m22) / 2`.
    fn populations_from_raw(m: &[f64; 9]) -> [f64; 9] {
        let mut f = [0.0f64; 9];
        for (i, f_i) in f.iter_mut().enumerate() {
            let gx = [
                dual_basis(0, CX[i]),
                dual_basis(1, CX[i]),
                dual_basis(2, CX[i]),
            ];
            let gy = [
                dual_basis(0, CY[i]),
                dual_basis(1, CY[i]),
                dual_basis(2, CY[i]),
            ];
            *f_i = ORDERS
                .iter()
                .zip(m.iter())
                .map(|(&(p, q), &m_pq)| m_pq * gx[p] * gy[q])
                .sum();
        }
        f
    }

    /// Equilibrium populations of the cascaded scheme: the populations whose
    /// central moments about `(ux, uy)` are exactly the Maxwellian ones,
    /// `(ρ, 0, 0, ρ/3, ρ/3, 0, 0, 0, ρ/9)`.
    ///
    /// Equivalent to the product form
    /// `f(a, b) = ρ ψ(a; u_x) ψ(b; u_y)` with `ψ(0; u) = 2/3 − u²` and
    /// `ψ(±1; u) = (1/3 + u² ± u) / 2`.  It is the fixed point of
    /// [`CascadedCollision2D::collide`] for every set of rates, and equals the
    /// second-order polynomial D2Q9 equilibrium up to terms of order `|u|³`.
    pub fn equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let kappa_eq = Self::equilibrium_central_moments(rho);
        Self::populations_from_raw(&Self::raw_from_central(&kappa_eq, ux, uy))
    }

    /// Relax the central moments `kappa` toward `target` at the rates of `self`.
    fn relax_central_moments(&self, kappa: &[f64; 9], target: &[f64; 9]) -> [f64; 9] {
        let mut out = [0.0f64; 9];
        // Mass and momentum (rate 0 unless the caller asks otherwise).
        for m in 0..3 {
            out[m] = kappa[m] - self.rates[m] * (kappa[m] - target[m]);
        }
        // Second order: the trace relaxes with the bulk rate, the deviator with
        // the shear rate; this keeps the viscosity isotropic.
        let trace = kappa[3] + kappa[4];
        let deviator = kappa[3] - kappa[4];
        let trace_eq = target[3] + target[4];
        let deviator_eq = target[3] - target[4];
        let trace_new = trace - self.rates[3] * (trace - trace_eq);
        let deviator_new = deviator - self.rates[4] * (deviator - deviator_eq);
        out[3] = 0.5 * (trace_new + deviator_new);
        out[4] = 0.5 * (trace_new - deviator_new);
        // Shear stress and the higher-order moments.
        for m in 5..9 {
            out[m] = kappa[m] - self.rates[m] * (kappa[m] - target[m]);
        }
        out
    }

    /// Perform a CLBM collision step for a single cell.
    ///
    /// `rho` sets the equilibrium moments `(ρ/3, ρ/3, ρ/9)` and `(ux, uy)` is
    /// the velocity of the frame in which the moments are relaxed; pass the
    /// macroscopic density and velocity of `f`.  With the default rates the
    /// collision conserves mass and momentum for any `rho`, `ux`, `uy`.
    ///
    /// Returns post-collision distributions.
    pub fn collide(&self, f: &[f64; 9], rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let kappa = Self::central_moments(f, ux, uy);
        let kappa_eq = Self::equilibrium_central_moments(rho);
        let relaxed = self.relax_central_moments(&kappa, &kappa_eq);
        Self::populations_from_raw(&Self::raw_from_central(&relaxed, ux, uy))
    }
}
