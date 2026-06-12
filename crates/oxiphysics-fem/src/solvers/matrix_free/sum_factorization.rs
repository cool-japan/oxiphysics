// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Matrix-free sum-factorized Poisson operator on a regular hexahedral mesh.
//!
//! This module implements the Kronbichler-Kormann (2012) matrix-free evaluation
//! of the high-order Poisson stiffness operator `y = K·x` on a regular
//! Cartesian mesh of the unit cube `[0,1]^3`. The mesh is partitioned into
//! `n_elems_per_dir^3` equal cubic elements; each element carries `(p+1)^3`
//! Lagrange degrees of freedom collocated on a Gauss-Lobatto-Legendre (GLL)
//! node grid, and the global system shares DOFs across element faces, edges,
//! and corners.
//!
//! The key idea of *sum-factorization* is that the `(p+1)^3 × (p+1)^3` dense
//! element stiffness matrix never has to be formed or stored. The action of the
//! reference-element gradient (and its transpose, used for integration) is a
//! tensor product of three 1D operators, so it can be applied as a sequence of
//! 1D contractions ([`contract_axis`]). The per-element work drops from
//! `O(p^6)` to `O(p^4)`, and the only stored data are a few `(p+1) × (p+1)` 1D
//! matrices — independent of the mesh size.
//!
//! On a regular cube mesh the geometry is constant: the Jacobian of the
//! reference-to-physical map is the diagonal `J = (h/2) I` with determinant
//! `(h/2)^3` and inverse-transpose `(2/h) I`. The Laplacian therefore reduces
//! to the isotropic sum of the three squared reference gradients weighted by
//! the GLL quadrature weights, which is exactly what [`SumFactPoisson::apply`]
//! evaluates.

use super::operator::MatrixFreeOperator;
use rayon::prelude::*;

/// Evaluate the Legendre polynomial `P_n` and its derivative `P'_n` at `x`.
///
/// Uses the standard three-term recurrence
/// `(k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}`, with the derivative obtained from
/// the identity `(1 - x^2) P'_n = n (P_{n-1} - x P_n)`. The latter is singular
/// only at `x = ±1`, where the well-known closed forms `P'_n(±1) = (±1)^{n-1}
/// n(n+1)/2` are substituted to keep the evaluation numerically clean.
fn legendre_p_dp(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }
    if n == 1 {
        return (x, 1.0);
    }
    let mut p_km1 = 1.0; // P_0
    let mut p_k = x; // P_1
    for k in 1..n {
        let kf = k as f64;
        let p_kp1 = ((2.0 * kf + 1.0) * x * p_k - kf * p_km1) / (kf + 1.0);
        p_km1 = p_k;
        p_k = p_kp1;
    }
    // Now p_k = P_n, p_km1 = P_{n-1}.
    let nf = n as f64;
    let denom = 1.0 - x * x;
    let dp = if denom.abs() < 1e-14 {
        // Endpoint closed form: P'_n(±1) = (±1)^{n-1} n(n+1)/2.
        // At x = +1 the sign is +1; at x = -1 it is (-1)^{n-1}, i.e. +1 when
        // (n-1) is even and -1 otherwise.
        let mag = nf * (nf + 1.0) / 2.0;
        let positive = x > 0.0 || (n - 1).is_multiple_of(2);
        if positive { mag } else { -mag }
    } else {
        nf * (p_km1 - x * p_k) / denom
    };
    (p_k, dp)
}

/// Return the `(p+1)`-point Gauss-Lobatto-Legendre (GLL) nodes and weights on
/// the reference interval `[-1, 1]`.
///
/// Closed-form values are used for `p ∈ {1, 2, 3}`; for higher `p` the interior
/// nodes are found by Newton iteration on the roots of `P'_p` (the derivative of
/// the degree-`p` Legendre polynomial), with the two endpoints fixed at `±1`.
/// The interior weights are `w_i = 2 / (p(p+1) [P_p(x_i)]^2)` and the endpoint
/// weights are `2 / (p(p+1))`. Nodes are returned sorted ascending.
pub fn gll_nodes_weights(p: usize) -> (Vec<f64>, Vec<f64>) {
    match p {
        0 => {
            // Degenerate: a single node. Treat as the midpoint with weight 2.
            (vec![0.0], vec![2.0])
        }
        1 => (vec![-1.0, 1.0], vec![1.0, 1.0]),
        2 => (vec![-1.0, 0.0, 1.0], vec![1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0]),
        3 => {
            let r = 1.0 / 5.0_f64.sqrt();
            (
                vec![-1.0, -r, r, 1.0],
                vec![1.0 / 6.0, 5.0 / 6.0, 5.0 / 6.0, 1.0 / 6.0],
            )
        }
        _ => gll_general(p),
    }
}

/// General GLL rule for `p >= 4` via Newton iteration on the roots of `P'_p`.
///
/// This path is also a correctness cross-check: it reproduces the closed-form
/// `p ∈ {1, 2, 3}` rules to roughly `1e-12`.
fn gll_general(p: usize) -> (Vec<f64>, Vec<f64>) {
    let n = p + 1;
    let pf = p as f64;
    let mut nodes = vec![0.0f64; n];
    nodes[0] = -1.0;
    nodes[n - 1] = 1.0;

    // Interior roots of P'_p. Seed from the Chebyshev-Gauss-Lobatto points
    // cos(π i / p), which interlace the GLL nodes well, then Newton-polish.
    // P'_p has p-1 interior roots (indices 1..p in node space).
    // `i` is the global node index in `1..(n-1)`.
    for (i, slot) in nodes.iter_mut().enumerate().take(n - 1).skip(1) {
        // Seed: cos(π i / p) gives a descending sequence; we store ascending,
        // so map index i to a descending Chebyshev point and rely on the final
        // sort to order them.
        let mut x = (std::f64::consts::PI * (i as f64) / pf).cos();
        for _ in 0..100 {
            // We need the root of P'_p. Newton needs P'_p and P''_p.
            // P''_p from differentiating (1-x^2)P'_p = p(P_{p-1} - x P_p):
            //   -2x P'_p + (1-x^2) P''_p = p(P'_{p-1} - P_p - x P'_p)
            // => P''_p = [p(P'_{p-1} - P_p - x P'_p) + 2x P'_p] / (1 - x^2).
            let (pp, dpp) = legendre_p_dp(p, x);
            let (_pm1, dpm1) = legendre_p_dp(p - 1, x);
            let denom = 1.0 - x * x;
            if denom.abs() < 1e-14 {
                break;
            }
            let ddpp = (pf * (dpm1 - pp - x * dpp) + 2.0 * x * dpp) / denom;
            if ddpp.abs() < 1e-300 {
                break;
            }
            let dx = dpp / ddpp;
            x -= dx;
            if dx.abs() < 1e-14 {
                break;
            }
        }
        *slot = x;
    }

    nodes.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    // Weights: endpoints 2/(p(p+1)); interior 2/(p(p+1) [P_p(x_i)]^2).
    let end_w = 2.0 / (pf * (pf + 1.0));
    let mut weights = vec![0.0f64; n];
    for i in 0..n {
        if i == 0 || i == n - 1 {
            weights[i] = end_w;
        } else {
            let (pp, _dpp) = legendre_p_dp(p, nodes[i]);
            weights[i] = 2.0 / (pf * (pf + 1.0) * pp * pp);
        }
    }
    (nodes, weights)
}

/// Build the 1D shape-value matrix `φ` and shape-derivative matrix `dφ/dξ`,
/// both row-major `n_quad × n_nodes`.
///
/// The `n_nodes` Lagrange basis functions are nodal at `nodes` and are sampled
/// at `quad_pts`:
/// - `phi[q*n_nodes + a]  = φ_a(quad_pts[q])`
/// - `dphi[q*n_nodes + a] = φ'_a(quad_pts[q])`
///
/// The value uses the standard Lagrange product
/// `φ_a(x) = Π_{m≠a} (x - x_m)/(x_a - x_m)`. The derivative uses the *stable*
/// double-loop form
/// `φ'_a(x) = Σ_{k≠a} [ 1/(x_a - x_k) · Π_{m≠a,k} (x - x_m)/(x_a - x_m) ]`,
/// which is exact even when `x` coincides with a node (as it does for GLL
/// collocation, where `quad_pts == nodes`). No special-casing of the identity
/// is done so the helper remains reusable for a general nodal basis.
pub(crate) fn lagrange_at(nodes: &[f64], quad_pts: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n_nodes = nodes.len();
    let n_quad = quad_pts.len();
    let mut phi = vec![0.0f64; n_quad * n_nodes];
    let mut dphi = vec![0.0f64; n_quad * n_nodes];

    for q in 0..n_quad {
        let x = quad_pts[q];
        for a in 0..n_nodes {
            // Value: Π_{m≠a} (x - x_m)/(x_a - x_m).
            let mut val = 1.0;
            for m in 0..n_nodes {
                if m == a {
                    continue;
                }
                val *= (x - nodes[m]) / (nodes[a] - nodes[m]);
            }
            phi[q * n_nodes + a] = val;

            // Derivative: Σ_{k≠a} [ 1/(x_a - x_k) · Π_{m≠a,k} (x-x_m)/(x_a-x_m) ].
            let mut d = 0.0;
            for k in 0..n_nodes {
                if k == a {
                    continue;
                }
                let mut prod = 1.0 / (nodes[a] - nodes[k]);
                for m in 0..n_nodes {
                    if m == a || m == k {
                        continue;
                    }
                    prod *= (x - nodes[m]) / (nodes[a] - nodes[m]);
                }
                d += prod;
            }
            dphi[q * n_nodes + a] = d;
        }
    }
    (phi, dphi)
}

/// Contract a 1D operator along one axis of a `(np × np × np)` tensor.
///
/// `input` and `output` are length `np³`, indexed `(i*np + j)*np + k`. `op1d`
/// is row-major `np × np`: `op1d[out*np + inp]` multiplies the input index
/// `inp` to accumulate into output index `out` along `axis`
/// (`0 → i`, `1 → j`, `2 → k`). `output` is fully overwritten.
pub(crate) fn contract_axis(input: &[f64], output: &mut [f64], op1d: &[f64], np: usize, axis: usize) {
    match axis {
        0 => {
            for j in 0..np {
                for k in 0..np {
                    for out_i in 0..np {
                        let mut s = 0.0;
                        for in_i in 0..np {
                            s += op1d[out_i * np + in_i] * input[(in_i * np + j) * np + k];
                        }
                        output[(out_i * np + j) * np + k] = s;
                    }
                }
            }
        }
        1 => {
            for i in 0..np {
                for k in 0..np {
                    for out_j in 0..np {
                        let mut s = 0.0;
                        for in_j in 0..np {
                            s += op1d[out_j * np + in_j] * input[(i * np + in_j) * np + k];
                        }
                        output[(i * np + out_j) * np + k] = s;
                    }
                }
            }
        }
        _ => {
            // axis == 2 (and any out-of-range value falls back to k for safety).
            for i in 0..np {
                for j in 0..np {
                    for out_k in 0..np {
                        let mut s = 0.0;
                        for in_k in 0..np {
                            s += op1d[out_k * np + in_k] * input[(i * np + j) * np + in_k];
                        }
                        output[(i * np + j) * np + out_k] = s;
                    }
                }
            }
        }
    }
}

/// Transpose a row-major `np × np` matrix.
fn transpose(m: &[f64], np: usize) -> Vec<f64> {
    let mut t = vec![0.0f64; np * np];
    for r in 0..np {
        for c in 0..np {
            t[c * np + r] = m[r * np + c];
        }
    }
    t
}

/// Matrix-free Poisson operator on a regular unit-cube hexahedral (Q_p) mesh.
///
/// The domain is the unit cube `[0,1]^3` partitioned into
/// `n_elems_per_dir^3` cubic elements of side `h = 1/n_elems_per_dir`.
/// Each element carries `(p+1)^3` Lagrange DOFs on a Gauss-Lobatto-Legendre
/// node grid; the global mesh has `(n_elems_per_dir·p + 1)^3` DOFs with shared
/// inter-element faces/edges/corners.
///
/// The stiffness apply `y += K·x` is evaluated element-locally via
/// sum-factorization (1D tensor contractions), never forming the global `K`.
pub struct SumFactPoisson {
    /// Number of elements per Cartesian direction.
    pub n_elems_per_dir: usize,
    /// Polynomial degree `p` (supported: 1, 2, 3; general GLL for higher).
    pub p: usize,
    /// Mesh spacing `h = 1 / n_elems_per_dir`.
    h: f64,
    /// 1D shape values `φ_a(ξ_q)`, row-major `(p+1) × (p+1)` = `n_quad × n_nodes`.
    phi: Vec<f64>,
    /// 1D shape derivatives `dφ_a/dξ|_{ξ_q}`, same layout as `phi`.
    dphi: Vec<f64>,
    /// Transpose of `phi` (used for the integration / test-function step).
    phi_t: Vec<f64>,
    /// Transpose of `dphi` (used for the integration / test-function step).
    dphi_t: Vec<f64>,
    /// 1D GLL quadrature weights, length `p+1`.
    weights: Vec<f64>,
    /// Total DOF count `(n_elems_per_dir·p + 1)^3`.
    n_dofs: usize,
}

impl SumFactPoisson {
    /// Construct the operator for a `n_elems_per_dir^3` mesh of degree-`p`
    /// elements.
    ///
    /// To honor the no-panic-in-production preference, the arguments are clamped
    /// to their valid minima (`p >= 1`, `n_elems_per_dir >= 1`) rather than
    /// asserted.
    pub fn new(n_elems_per_dir: usize, p: usize) -> Self {
        let p = p.max(1);
        let n_elems_per_dir = n_elems_per_dir.max(1);
        let h = 1.0 / n_elems_per_dir as f64;

        let (nodes, weights) = gll_nodes_weights(p);
        let quad_pts = nodes.clone();
        let (phi, dphi) = lagrange_at(&nodes, &quad_pts);
        let np = p + 1;
        let phi_t = transpose(&phi, np);
        let dphi_t = transpose(&dphi, np);

        let nodes_per_dir = n_elems_per_dir * p + 1;
        let n_dofs = nodes_per_dir * nodes_per_dir * nodes_per_dir;

        Self {
            n_elems_per_dir,
            p,
            h,
            phi,
            dphi,
            phi_t,
            dphi_t,
            weights,
            n_dofs,
        }
    }

    /// Number of global DOF nodes per direction.
    #[inline]
    fn nodes_per_dir(&self) -> usize {
        self.n_elems_per_dir * self.p + 1
    }

    /// Global DOF index of local node `(i,j,k)` in element `(ex,ey,ez)`.
    ///
    /// Global node coords: `(gi,gj,gk) = (ex*p+i, ey*p+j, ez*p+k)`.
    /// Flattened: `(gi*n + gj)*n + gk` with `n = nodes_per_dir`.
    #[inline]
    fn global_dof(&self, ex: usize, ey: usize, ez: usize, i: usize, j: usize, k: usize) -> usize {
        let n = self.nodes_per_dir();
        let gi = ex * self.p + i;
        let gj = ey * self.p + j;
        let gk = ez * self.p + k;
        (gi * n + gj) * n + gk
    }

    /// Total DOF count `(n_elems_per_dir·p + 1)^3`.
    pub fn n_dofs(&self) -> usize {
        self.n_dofs
    }

    /// Polynomial degree `p`.
    pub fn p(&self) -> usize {
        self.p
    }

    /// Number of elements per Cartesian direction.
    pub fn n_elems_per_dir(&self) -> usize {
        self.n_elems_per_dir
    }

    /// Number of nodes (DOFs) along one Cartesian axis.
    pub fn dofs_per_dir(&self) -> usize {
        self.nodes_per_dir()
    }

    /// Apply the reference-element stiffness (scaled by the regular-mesh
    /// geometry) to a gathered element vector `u_e` (length `np³`), accumulating
    /// into `y_e`.
    ///
    /// This realises `y_e += K_e · u_e` where, with GLL quadrature,
    /// `K_e[m,n] = Σ_d Σ_q geo_q · G^d_m(q) · G^d_n(q)` and `G^d_m(q)` is the
    /// `d`-th reference derivative of basis `m` at quadrature point `q`. The
    /// forward step (steps 1-2) evaluates the three reference gradients of
    /// `u_e` at every quad point by sum-factorization; the integration step
    /// (step 4) applies the transpose contractions to assemble the test-function
    /// pairing.
    pub(crate) fn element_apply(&self, u_e: &[f64], y_e: &mut [f64]) {
        let np = self.p + 1;
        let n3 = np * np * np;

        // Geometry + quadrature factor (assembled per quad point below).
        //   (2/h)^2 : the diagonal entry of J^{-T} J^{-1} for the isotropic
        //             cube map J = (h/2) I.
        //   (h/2)^3 : the Jacobian determinant |J|.
        // The 3D GLL weight w_a·w_b·w_c is folded in per quad point.
        let two_over_h = 2.0 / self.h;
        let half_h = self.h / 2.0;
        let geo = two_over_h * two_over_h * half_h * half_h * half_h;

        // Scratch buffers for the two intermediate contractions.
        let mut t1 = vec![0.0f64; n3];
        let mut t2 = vec![0.0f64; n3];

        // --- Step 2: reference gradients of u_e at all quad points. ---
        // grad ξ: dphi along axis 0, phi along axes 1 & 2.
        let mut g_xi = vec![0.0f64; n3];
        contract_axis(u_e, &mut t1, &self.dphi, np, 0);
        contract_axis(&t1, &mut t2, &self.phi, np, 1);
        contract_axis(&t2, &mut g_xi, &self.phi, np, 2);

        // grad η: phi along axis 0, dphi along axis 1, phi along axis 2.
        let mut g_eta = vec![0.0f64; n3];
        contract_axis(u_e, &mut t1, &self.phi, np, 0);
        contract_axis(&t1, &mut t2, &self.dphi, np, 1);
        contract_axis(&t2, &mut g_eta, &self.phi, np, 2);

        // grad ζ: phi along axes 0 & 1, dphi along axis 2.
        let mut g_zeta = vec![0.0f64; n3];
        contract_axis(u_e, &mut t1, &self.phi, np, 0);
        contract_axis(&t1, &mut t2, &self.phi, np, 1);
        contract_axis(&t2, &mut g_zeta, &self.dphi, np, 2);

        // --- Step 3: scale each gradient component by geo · w_a·w_b·w_c. ---
        for a in 0..np {
            let wa = self.weights[a];
            for b in 0..np {
                let wab = wa * self.weights[b];
                for c in 0..np {
                    let f = geo * wab * self.weights[c];
                    let q = (a * np + b) * np + c;
                    g_xi[q] *= f;
                    g_eta[q] *= f;
                    g_zeta[q] *= f;
                }
            }
        }

        // --- Step 4: integrate back via the transpose contractions. ---
        // ξ contribution: dphi_t @ axis0, phi_t @ axes1&2.
        let mut c_acc = vec![0.0f64; n3];
        contract_axis(&g_xi, &mut t1, &self.dphi_t, np, 0);
        contract_axis(&t1, &mut t2, &self.phi_t, np, 1);
        contract_axis(&t2, &mut c_acc, &self.phi_t, np, 2);
        for idx in 0..n3 {
            y_e[idx] += c_acc[idx];
        }

        // η contribution: phi_t @ axis0, dphi_t @ axis1, phi_t @ axis2.
        contract_axis(&g_eta, &mut t1, &self.phi_t, np, 0);
        contract_axis(&t1, &mut t2, &self.dphi_t, np, 1);
        contract_axis(&t2, &mut c_acc, &self.phi_t, np, 2);
        for idx in 0..n3 {
            y_e[idx] += c_acc[idx];
        }

        // ζ contribution: phi_t @ axes0&1, dphi_t @ axis2.
        contract_axis(&g_zeta, &mut t1, &self.phi_t, np, 0);
        contract_axis(&t1, &mut t2, &self.phi_t, np, 1);
        contract_axis(&t2, &mut c_acc, &self.dphi_t, np, 2);
        for idx in 0..n3 {
            y_e[idx] += c_acc[idx];
        }
    }

    /// Return the reciprocal diagonal `1/K[i,i]` of the assembled stiffness,
    /// computed matrix-free by accumulating element diagonals. Used to build a
    /// Jacobi preconditioner without forming `K`.
    ///
    /// The diagonal assembled here is exactly the diagonal of the operator that
    /// [`SumFactPoisson::apply`] applies:
    /// `K_e[m,m] = Σ_d Σ_q geo_q · (G^d_m(q))^2`, where the reference derivative
    /// of element basis `m = (im,jm,km)` in direction `d` at quad
    /// `q = (a,b,c)` is the tensor product of the stored 1D `phi`/`dphi` rows.
    ///
    /// `dirichlet` masks boundary DOFs: where `dirichlet[i]` is true the returned
    /// value is `1.0` so the BC-handled system uses an identity row.
    pub fn jacobi_diag_inverse(&self, dirichlet: &[bool]) -> Vec<f64> {
        let np = self.p + 1;
        let nelem = self.n_elems_per_dir;
        let two_over_h = 2.0 / self.h;
        let half_h = self.h / 2.0;
        let geo = two_over_h * two_over_h * half_h * half_h * half_h;

        // Precompute the 3D quadrature weight at each quad point.
        // Assemble the global diagonal element by element (setup-time, serial).
        let mut diag = vec![0.0f64; self.n_dofs];

        let total_elems = nelem * nelem * nelem;
        for e in 0..total_elems {
            let ex = e / (nelem * nelem);
            let rem = e % (nelem * nelem);
            let ey = rem / nelem;
            let ez = rem % nelem;

            for im in 0..np {
                for jm in 0..np {
                    for km in 0..np {
                        // K_e[m,m] = Σ_d Σ_q geo · w_q · (G^d_m(q))^2.
                        let mut diag_mm = 0.0;
                        for a in 0..np {
                            let wa = self.weights[a];
                            // 1D factors for node m in each direction at quad a/b/c.
                            let phi_a_im = self.phi[a * np + im];
                            let dphi_a_im = self.dphi[a * np + im];
                            for b in 0..np {
                                let wab = wa * self.weights[b];
                                let phi_b_jm = self.phi[b * np + jm];
                                let dphi_b_jm = self.dphi[b * np + jm];
                                for c in 0..np {
                                    let w = geo * wab * self.weights[c];
                                    let phi_c_km = self.phi[c * np + km];
                                    let dphi_c_km = self.dphi[c * np + km];

                                    // G^ξ = dphi_a_im · phi_b_jm · phi_c_km
                                    let g_xi = dphi_a_im * phi_b_jm * phi_c_km;
                                    // G^η = phi_a_im · dphi_b_jm · phi_c_km
                                    let g_eta = phi_a_im * dphi_b_jm * phi_c_km;
                                    // G^ζ = phi_a_im · phi_b_jm · dphi_c_km
                                    let g_zeta = phi_a_im * phi_b_jm * dphi_c_km;

                                    diag_mm += w * (g_xi * g_xi + g_eta * g_eta + g_zeta * g_zeta);
                                }
                            }
                        }
                        let g = self.global_dof(ex, ey, ez, im, jm, km);
                        diag[g] += diag_mm;
                    }
                }
            }
        }

        // Build the reciprocal diagonal, with identity rows for Dirichlet DOFs.
        // A Dirichlet-masked DOF uses an identity row (`diag_inv = 1`); a
        // (numerically) zero diagonal is likewise floored to `1` to keep the
        // Jacobi preconditioner well defined. Both cases coincide in value.
        let mut diag_inv = vec![0.0f64; self.n_dofs];
        for i in 0..self.n_dofs {
            let is_dirichlet = dirichlet.get(i).copied().unwrap_or(false);
            diag_inv[i] = if is_dirichlet || diag[i].abs() < 1e-300 {
                1.0
            } else {
                1.0 / diag[i]
            };
        }
        diag_inv
    }
}

impl MatrixFreeOperator for SumFactPoisson {
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        // y += K·x. Parallelise over elements; combine via per-thread partial
        // vectors (fold/reduce) to avoid data races on shared DOFs.
        let np = self.p + 1;
        let nelem = self.n_elems_per_dir;
        let total_elems = nelem * nelem * nelem;
        let partial = (0..total_elems)
            .into_par_iter()
            .fold(
                || vec![0.0f64; self.n_dofs],
                |mut acc, e| {
                    let ex = e / (nelem * nelem);
                    let rem = e % (nelem * nelem);
                    let ey = rem / nelem;
                    let ez = rem % nelem;
                    // gather
                    let mut u_e = vec![0.0f64; np * np * np];
                    for i in 0..np {
                        for j in 0..np {
                            for k in 0..np {
                                let g = self.global_dof(ex, ey, ez, i, j, k);
                                u_e[(i * np + j) * np + k] = x[g];
                            }
                        }
                    }
                    let mut y_e = vec![0.0f64; np * np * np];
                    self.element_apply(&u_e, &mut y_e);
                    // scatter into the thread-local accumulator
                    for i in 0..np {
                        for j in 0..np {
                            for k in 0..np {
                                let g = self.global_dof(ex, ey, ez, i, j, k);
                                acc[g] += y_e[(i * np + j) * np + k];
                            }
                        }
                    }
                    acc
                },
            )
            .reduce(
                || vec![0.0f64; self.n_dofs],
                |mut a, b| {
                    for (ai, bi) in a.iter_mut().zip(b.iter()) {
                        *ai += bi;
                    }
                    a
                },
            );
        for (yi, pi) in y.iter_mut().zip(partial.iter()) {
            *yi += pi;
        }
    }

    fn dim(&self) -> usize {
        self.n_dofs
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gll_p2_weights() {
        let (nodes, weights) = gll_nodes_weights(2);
        let exp_n = [-1.0, 0.0, 1.0];
        let exp_w = [1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0];
        for i in 0..3 {
            assert!((nodes[i] - exp_n[i]).abs() < 1e-12, "node {i}");
            assert!((weights[i] - exp_w[i]).abs() < 1e-12, "weight {i}");
        }
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-12, "weights must sum to 2");
    }

    #[test]
    fn gll_p3_nodes() {
        let (nodes, weights) = gll_nodes_weights(3);
        let r = 1.0 / 5.0_f64.sqrt();
        let exp_n = [-1.0, -r, r, 1.0];
        for i in 0..4 {
            assert!((nodes[i] - exp_n[i]).abs() < 1e-12, "node {i}");
        }
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-12, "weights must sum to 2");
    }

    #[test]
    fn gll_general_matches_closed_form() {
        // The Newton path must reproduce the closed forms for p = 2, 3.
        for p in [2usize, 3] {
            let (cn, cw) = gll_nodes_weights(p);
            let (gn, gw) = gll_general(p);
            for i in 0..(p + 1) {
                assert!((cn[i] - gn[i]).abs() < 1e-11, "p={p} node {i}");
                assert!((cw[i] - gw[i]).abs() < 1e-11, "p={p} weight {i}");
            }
        }
    }

    #[test]
    fn gll_p4_weights_sum_to_two() {
        let (nodes, weights) = gll_nodes_weights(4);
        assert_eq!(nodes.len(), 5);
        // Symmetric nodes about 0; endpoints ±1.
        assert!((nodes[0] + 1.0).abs() < 1e-12);
        assert!((nodes[4] - 1.0).abs() < 1e-12);
        assert!(nodes[2].abs() < 1e-12, "middle node is 0");
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-11, "weights must sum to 2");
    }

    #[test]
    fn phi_is_identity_for_gll() {
        for p in [1usize, 2, 3] {
            let (nodes, _w) = gll_nodes_weights(p);
            let (phi, _dphi) = lagrange_at(&nodes, &nodes);
            let np = p + 1;
            for q in 0..np {
                for a in 0..np {
                    let expected = if q == a { 1.0 } else { 0.0 };
                    assert!(
                        (phi[q * np + a] - expected).abs() < 1e-12,
                        "p={p} phi[{q},{a}] = {}",
                        phi[q * np + a]
                    );
                }
            }
        }
    }

    #[test]
    fn partition_of_unity_dphi_rows_sum_zero() {
        // The derivative of the constant function (Σ_a φ_a ≡ 1) is zero, so each
        // row of dphi must sum to ~0.
        for p in [1usize, 2, 3, 4] {
            let (nodes, _w) = gll_nodes_weights(p);
            let (_phi, dphi) = lagrange_at(&nodes, &nodes);
            let np = p + 1;
            for q in 0..np {
                let mut s = 0.0;
                for a in 0..np {
                    s += dphi[q * np + a];
                }
                assert!(s.abs() < 1e-11, "p={p} row {q} sum {s:.3e}");
            }
        }
    }

    #[test]
    fn contract_axis_matches_naive() {
        // p = 3 → np = 4. Deterministic pseudo-random op1d and input.
        let np = 4;
        let n3 = np * np * np;
        let mut op = vec![0.0f64; np * np];
        for (idx, v) in op.iter_mut().enumerate() {
            *v = ((idx * 37 + 11) % 17) as f64 / 7.0 - 1.0;
        }
        let mut input = vec![0.0f64; n3];
        for (idx, v) in input.iter_mut().enumerate() {
            *v = ((idx * 53 + 7) % 23) as f64 / 11.0 - 1.0;
        }

        for axis in 0..3 {
            let mut out = vec![0.0f64; n3];
            contract_axis(&input, &mut out, &op, np, axis);

            // Naive reference.
            let mut refout = vec![0.0f64; n3];
            for i in 0..np {
                for j in 0..np {
                    for k in 0..np {
                        let mut s = 0.0;
                        for t in 0..np {
                            let (oi, ij, ik) = match axis {
                                0 => (t, j, k),
                                1 => (i, t, k),
                                _ => (i, j, t),
                            };
                            let out_pos = match axis {
                                0 => i,
                                1 => j,
                                _ => k,
                            };
                            s += op[out_pos * np + t] * input[(oi * np + ij) * np + ik];
                        }
                        refout[(i * np + j) * np + k] = s;
                    }
                }
            }
            for idx in 0..n3 {
                assert!(
                    (out[idx] - refout[idx]).abs() < 1e-12,
                    "axis {axis} idx {idx}: {} vs {}",
                    out[idx],
                    refout[idx]
                );
            }
        }
    }

    #[test]
    fn operator_dims_and_dofs() {
        let op = SumFactPoisson::new(2, 2);
        // nodes_per_dir = 2*2 + 1 = 5; n_dofs = 125.
        assert_eq!(op.dofs_per_dir(), 5);
        assert_eq!(op.n_dofs(), 125);
        assert_eq!(op.dim(), 125);
        assert_eq!(op.p(), 2);
        assert_eq!(op.n_elems_per_dir(), 2);
    }

    #[test]
    fn constant_field_in_nullspace() {
        // The Poisson stiffness annihilates constants: K·1 = 0 everywhere.
        let op = SumFactPoisson::new(2, 2);
        let n = op.n_dofs();
        let x = vec![1.0f64; n];
        let mut y = vec![0.0f64; n];
        op.apply(&x, &mut y);
        let max = y.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(max < 1e-10, "K·1 should vanish, max |y| = {max:.3e}");
    }

    #[test]
    fn jacobi_diag_inverse_matches_apply_diagonal() {
        // The reciprocal diagonal from jacobi_diag_inverse must equal the
        // reciprocal of K[i,i] obtained by probing apply with unit vectors.
        let op = SumFactPoisson::new(1, 2);
        let n = op.n_dofs();
        let no_bc = vec![false; n];
        let diag_inv = op.jacobi_diag_inverse(&no_bc);

        for i in 0..n {
            let mut e = vec![0.0f64; n];
            e[i] = 1.0;
            let mut col = vec![0.0f64; n];
            op.apply(&e, &mut col);
            let kii = col[i];
            assert!(kii > 0.0, "diagonal must be positive at {i}: {kii}");
            assert!(
                (diag_inv[i] - 1.0 / kii).abs() < 1e-9,
                "diag_inv[{i}] = {} vs 1/{} = {}",
                diag_inv[i],
                kii,
                1.0 / kii
            );
        }
    }

    #[test]
    fn operator_is_symmetric() {
        // K must be symmetric: xᵀ K y == yᵀ K x.
        let op = SumFactPoisson::new(2, 2);
        let n = op.n_dofs();
        let mut xv = vec![0.0f64; n];
        let mut yv = vec![0.0f64; n];
        for i in 0..n {
            xv[i] = ((i * 31 + 5) % 13) as f64 / 6.0 - 1.0;
            yv[i] = ((i * 17 + 3) % 11) as f64 / 5.0 - 1.0;
        }
        let mut kx = vec![0.0f64; n];
        let mut ky = vec![0.0f64; n];
        op.apply(&xv, &mut kx);
        op.apply(&yv, &mut ky);
        let yt_kx: f64 = yv.iter().zip(kx.iter()).map(|(a, b)| a * b).sum();
        let xt_ky: f64 = xv.iter().zip(ky.iter()).map(|(a, b)| a * b).sum();
        assert!(
            (yt_kx - xt_ky).abs() < 1e-9 * (1.0 + yt_kx.abs()),
            "symmetry: {yt_kx} vs {xt_ky}"
        );
    }
}
