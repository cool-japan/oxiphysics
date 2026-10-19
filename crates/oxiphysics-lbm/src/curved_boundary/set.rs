// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Curved wall links, lattice domains and the boundary-set driver.

use std::collections::BTreeMap;

use super::CurvedBcMethod;
use super::kernels::{
    LinkPopulations, WallLink, bouzidi_bc, bouzidi_quadratic_bc, filippova_hanel,
};
use super::lattice::CurvedLattice;
use crate::error::{Error, Result};

// ---------------------------------------------------------------------------
// CurvedWallNode
// ---------------------------------------------------------------------------

/// One wall-cut lattice link adjacent to a curved wall.
///
/// The link starts at the fluid node `fluid_idx`, points along lattice
/// direction `dir` into the solid node `wall_idx`, and is cut by the wall at
/// `x_f + q c_dir`.  The boundary scheme reconstructs the population
/// `f[fluid_idx][opposite(dir)]` that would otherwise stream in from the solid.
#[derive(Debug, Clone, PartialEq)]
pub struct CurvedWallNode {
    /// Linear index of the fluid-side node.
    pub fluid_idx: usize,
    /// Linear index of the solid-side node (`x_f + c_dir`).
    pub wall_idx: usize,
    /// Lattice direction pointing from the fluid node into the wall.
    pub dir: usize,
    /// Fraction of the link lying in the fluid, `q ∈ (0, 1]`.
    ///
    /// `q = 1` places the wall on the solid node, `q → 0` on the fluid node.
    pub q: f64,
    /// Wall velocity at the intersection point `x_f + q c_dir`.
    pub u_wall: [f64; 3],
    /// Fluid node `x_f − c_dir`, if it exists and is fluid.
    pub upstream: Option<usize>,
    /// Fluid node `x_f − 2 c_dir`, if it exists and is fluid.
    pub upstream2: Option<usize>,
}

impl CurvedWallNode {
    /// Create a stationary-wall link without upstream neighbours.
    ///
    /// # Errors
    /// Returns [`Error::ParameterOutOfBounds`] if `q ∉ (0, 1]`.
    pub fn new(fluid_idx: usize, wall_idx: usize, dir: usize, q: f64) -> Result<Self> {
        if !q.is_finite() || q <= 0.0 || q > 1.0 {
            return Err(Error::ParameterOutOfBounds {
                name: "q".to_string(),
                value: q,
                min: 0.0,
                max: 1.0,
            });
        }
        Ok(Self {
            fluid_idx,
            wall_idx,
            dir,
            q,
            u_wall: [0.0; 3],
            upstream: None,
            upstream2: None,
        })
    }

    /// Set the wall velocity at the link intersection point.
    pub fn with_wall_velocity(mut self, u_wall: [f64; 3]) -> Self {
        self.u_wall = u_wall;
        self
    }

    /// Set the upstream fluid nodes `x_f − c_dir` and `x_f − 2 c_dir`.
    pub fn with_upstream(mut self, upstream: Option<usize>, upstream2: Option<usize>) -> Self {
        self.upstream = upstream;
        self.upstream2 = upstream2;
        self
    }
}

// ---------------------------------------------------------------------------
// LatticeDomain
// ---------------------------------------------------------------------------

/// Cartesian lattice of `nx × ny × nz` nodes with spacing `grid_dx`.
///
/// Node `(ix, iy, iz)` sits at `(ix, iy, iz) · grid_dx` and has linear index
/// `ix + iy·nx + iz·nx·ny`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LatticeDomain {
    /// Nodes in x.
    pub nx: usize,
    /// Nodes in y.
    pub ny: usize,
    /// Nodes in z (1 for 2-D).
    pub nz: usize,
    /// Lattice spacing in physical units.
    pub grid_dx: f64,
    /// Periodicity flags per axis.
    pub periodic: [bool; 3],
}

impl LatticeDomain {
    /// Create a non-periodic domain.
    ///
    /// # Errors
    /// Returns an error for a zero dimension or a non-positive spacing.
    pub fn new(nx: usize, ny: usize, nz: usize, grid_dx: f64) -> Result<Self> {
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::InvalidGridDimensions3D { nx, ny, nz });
        }
        if !grid_dx.is_finite() || grid_dx <= 0.0 {
            return Err(Error::ParameterOutOfBounds {
                name: "grid_dx".to_string(),
                value: grid_dx,
                min: f64::MIN_POSITIVE,
                max: f64::MAX,
            });
        }
        Ok(Self {
            nx,
            ny,
            nz,
            grid_dx,
            periodic: [false; 3],
        })
    }

    /// Set the per-axis periodicity flags.
    pub fn with_periodic(mut self, periodic: [bool; 3]) -> Self {
        self.periodic = periodic;
        self
    }

    /// Total number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    /// Linear index of node `(ix, iy, iz)`.
    pub fn linear_index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix + iy * self.nx + iz * self.nx * self.ny
    }

    /// Integer coordinates of a linear index.
    pub fn coords(&self, idx: usize) -> [usize; 3] {
        let plane = self.nx * self.ny;
        [idx % self.nx, (idx % plane) / self.nx, idx / plane]
    }

    /// Physical position of node `(ix, iy, iz)`.
    pub fn position(&self, node: [usize; 3]) -> [f64; 3] {
        [
            node[0] as f64 * self.grid_dx,
            node[1] as f64 * self.grid_dx,
            node[2] as f64 * self.grid_dx,
        ]
    }

    /// Node reached from `node` by the integer offset `c`, honouring
    /// periodicity; `None` if it leaves a non-periodic axis.
    pub fn neighbour(&self, node: [usize; 3], c: [i32; 3]) -> Option<[usize; 3]> {
        let dims = [self.nx, self.ny, self.nz];
        let mut out = [0usize; 3];
        for d in 0..3 {
            let n = i64::try_from(dims[d]).ok()?;
            let p = i64::try_from(node[d]).ok()? + i64::from(c[d]);
            let wrapped = if self.periodic[d] {
                p.rem_euclid(n)
            } else if p < 0 || p >= n {
                return None;
            } else {
                p
            };
            out[d] = usize::try_from(wrapped).ok()?;
        }
        Some(out)
    }
}

// ---------------------------------------------------------------------------
// CurvedBoundarySet
// ---------------------------------------------------------------------------

/// A collection of curved-wall links sharing one boundary scheme.
#[derive(Debug, Clone)]
pub struct CurvedBoundarySet {
    /// All wall-cut links of this boundary.
    pub nodes: Vec<CurvedWallNode>,
    /// Numerical scheme used for the links.
    pub method: CurvedBcMethod,
    /// Uniform wall velocity assigned to links generated by
    /// [`CurvedBoundarySet::from_geometry`].  [`CurvedBoundarySet::apply`]
    /// uses the per-link [`CurvedWallNode::u_wall`].
    pub wall_velocity: [f64; 3],
}

/// Macroscopic state recovered from post-collision populations.
struct NodeState {
    rho: f64,
    u: [f64; 3],
}

impl CurvedBoundarySet {
    /// Create an empty set.
    pub fn new(method: CurvedBcMethod, wall_velocity: [f64; 3]) -> Self {
        Self {
            nodes: Vec::new(),
            method,
            wall_velocity,
        }
    }

    /// Add a wall-cut link.
    pub fn add_node(&mut self, node: CurvedWallNode) {
        self.nodes.push(node);
    }

    /// Assign the same wall velocity to the set default and to every link.
    pub fn set_uniform_wall_velocity(&mut self, u_wall: [f64; 3]) {
        self.wall_velocity = u_wall;
        for node in &mut self.nodes {
            node.u_wall = u_wall;
        }
    }

    /// Detect all wall-cut links of a geometry given by a signed-distance
    /// function (positive in the fluid, non-positive in the solid).
    ///
    /// Uses the linear Bouzidi scheme and assigns [`Self::wall_velocity`]
    /// (zero) to every link.  See [`Self::from_geometry_with`].
    pub fn from_geometry(
        sdf: &dyn Fn([f64; 3]) -> f64,
        domain: &LatticeDomain,
        lattice: &CurvedLattice,
    ) -> Self {
        let zero = |_: [f64; 3]| [0.0; 3];
        Self::from_geometry_with(sdf, &zero, domain, lattice, CurvedBcMethod::Bouzidi)
    }

    /// Detect all wall-cut links of a geometry and attach local wall velocities.
    ///
    /// Every fluid node (SDF > 0) is tested along **every** non-rest lattice
    /// direction (axis and diagonal links).  A link is cut when the target
    /// node is solid (SDF ≤ 0).  The intersection fraction `q` is located by
    /// bisection on the SDF along the link, so an exact SDF gives an exact
    /// `q`.  The wall velocity of each link is `wall_velocity(x_f + q c_i Δx)`,
    /// and the upstream fluid nodes `x_f − c_i`, `x_f − 2 c_i` are recorded
    /// for the interpolation schemes.
    pub fn from_geometry_with(
        sdf: &dyn Fn([f64; 3]) -> f64,
        wall_velocity: &dyn Fn([f64; 3]) -> [f64; 3],
        domain: &LatticeDomain,
        lattice: &CurvedLattice,
        method: CurvedBcMethod,
    ) -> Self {
        let mut set = Self::new(method, [0.0; 3]);
        let dx = domain.grid_dx;
        let is_fluid = |node: [usize; 3]| sdf(domain.position(node)) > 0.0;
        for idx in 0..domain.num_nodes() {
            let node = domain.coords(idx);
            let pos = domain.position(node);
            let phi_f = sdf(pos);
            if phi_f <= 0.0 {
                continue;
            }
            for (dir, c) in lattice.velocities().iter().enumerate() {
                if *c == [0, 0, 0] {
                    continue;
                }
                let Some(target) = domain.neighbour(node, *c) else {
                    continue;
                };
                let phi_s = sdf(domain.position(target));
                if phi_s > 0.0 {
                    continue;
                }
                let step = [
                    f64::from(c[0]) * dx,
                    f64::from(c[1]) * dx,
                    f64::from(c[2]) * dx,
                ];
                let q = link_fraction(sdf, pos, step, phi_f, phi_s);
                let hit = [
                    pos[0] + q * step[0],
                    pos[1] + q * step[1],
                    pos[2] + q * step[2],
                ];
                let back = [-c[0], -c[1], -c[2]];
                let up1 = domain.neighbour(node, back).filter(|n| is_fluid(*n));
                let up2 = up1
                    .and_then(|n| domain.neighbour(n, back))
                    .filter(|n| is_fluid(*n));
                set.nodes.push(CurvedWallNode {
                    fluid_idx: idx,
                    wall_idx: domain.linear_index(target[0], target[1], target[2]),
                    dir,
                    q,
                    u_wall: wall_velocity(hit),
                    upstream: up1.map(|n| domain.linear_index(n[0], n[1], n[2])),
                    upstream2: up2.map(|n| domain.linear_index(n[0], n[1], n[2])),
                });
            }
        }
        set
    }

    /// Apply the boundary scheme to every link.
    ///
    /// * `f_post` — post-collision, pre-streaming populations `f̃[node][dir]`
    /// * `f`      — post-streaming populations; only the unknown entries
    ///   `f[fluid_idx][opposite(dir)]` of the links are overwritten
    /// * `lattice` — velocity set matching the population layout
    /// * `omega`   — BGK relaxation frequency `1/τ`, `0 < ω < 2`
    /// * `body_force` — optional body force per node (Guo forcing).  It is used
    ///   to recover the physical velocity `u = (Σ c f̃ − F/2)/ρ` from
    ///   post-collision populations and, for [`CurvedBcMethod::MomentBased`],
    ///   to impose `Σ c f = ρ u − F/2` on the reconstructed node.
    ///
    /// Scheme fallbacks when upstream nodes are missing (e.g. in narrow gaps):
    /// quadratic BFL → linear BFL → half-way bounce-back with the Ladd term;
    /// Mei–Luo–Shyy with `q < 1/2` → half-way bounce-back with the Ladd term.
    ///
    /// # Errors
    /// Returns an error on mismatched array sizes, out-of-range link indices,
    /// `ω ∉ (0, 2)`, a non-positive density, or a singular Mei–Luo–Shyy
    /// weighting (`τ = 2` with `q < 1/2`).
    pub fn apply(
        &self,
        f_post: &[Vec<f64>],
        f: &mut [Vec<f64>],
        lattice: &CurvedLattice,
        omega: f64,
        body_force: Option<&dyn Fn(usize) -> [f64; 3]>,
    ) -> Result<()> {
        if f_post.len() != f.len() {
            return Err(Error::ComponentCountMismatch {
                expected: f_post.len(),
                found: f.len(),
                field: "post-streaming populations".to_string(),
            });
        }
        if !omega.is_finite() || omega <= 0.0 || omega >= 2.0 {
            return Err(Error::ParameterOutOfBounds {
                name: "omega".to_string(),
                value: omega,
                min: 0.0,
                max: 2.0,
            });
        }
        let force_at = |idx: usize| body_force.map_or([0.0; 3], |bf| bf(idx));
        match self.method {
            CurvedBcMethod::MomentBased => self.apply_moment_based(f_post, f, lattice, &force_at),
            method => {
                for node in &self.nodes {
                    let ctx = LinkContext::gather(node, f_post, lattice, &force_at)?;
                    let value = ctx.reconstruct(method, omega)?;
                    write_population(f, node.fluid_idx, ctx.opp, value)?;
                }
                Ok(())
            }
        }
    }

    /// Moment-based reconstruction (see [`CurvedBcMethod::MomentBased`]).
    ///
    /// For every boundary fluid node the unknown populations are first
    /// estimated with linear BFL.  The node velocity imposed by the wall is
    /// then estimated on each cut link by linear interpolation between the
    /// wall (`u_w` at distance `q`) and the upstream node (`u_ff` at distance
    /// 1): `u_t = (u_w + q u_ff)/(1 + q)`, averaged over the cut links.  The
    /// unknown populations receive the minimal `w`-weighted correction `δ`
    /// that enforces the first-moment constraint
    /// `Σ_k (c_k − u_t) f_k = −F/2` (i.e. `Σ c f = ρ u_t − F/2` with the
    /// density left free, as in Zou–He).  The constraint system is solved
    /// through the pseudo-inverse of `A W Aᵀ`, so nodes whose unknown set
    /// cannot span every momentum component receive the least-squares
    /// correction.
    fn apply_moment_based(
        &self,
        f_post: &[Vec<f64>],
        f: &mut [Vec<f64>],
        lattice: &CurvedLattice,
        force_at: &dyn Fn(usize) -> [f64; 3],
    ) -> Result<()> {
        let mut groups: BTreeMap<usize, Vec<&CurvedWallNode>> = BTreeMap::new();
        for node in &self.nodes {
            groups.entry(node.fluid_idx).or_default().push(node);
        }
        let dims = lattice.dimensions();
        for (fluid_idx, links) in groups {
            let mut unknown: Vec<usize> = Vec::with_capacity(links.len());
            let mut u_sum = [0.0; 3];
            let mut n_est = 0usize;
            for node in &links {
                let ctx = LinkContext::gather(node, f_post, lattice, force_at)?;
                let value = ctx.reconstruct(CurvedBcMethod::Bouzidi, 1.0)?;
                write_population(f, fluid_idx, ctx.opp, value)?;
                if !unknown.contains(&ctx.opp) {
                    unknown.push(ctx.opp);
                }
                if let Some(u_ff) = ctx.u_ff {
                    let q = node.q;
                    for (d, acc) in u_sum.iter_mut().enumerate() {
                        *acc += (node.u_wall[d] + q * u_ff[d]) / (1.0 + q);
                    }
                    n_est += 1;
                }
            }
            if n_est == 0 {
                continue;
            }
            let u_t = u_sum.map(|s| s / n_est as f64);
            let force = force_at(fluid_idx);
            let pops = f.get(fluid_idx).ok_or(Error::IndexOutOfBounds {
                x: fluid_idx,
                y: 0,
                nx: f_post.len(),
                ny: 1,
            })?;
            if pops.len() != lattice.num_directions() {
                return Err(Error::ComponentCountMismatch {
                    expected: lattice.num_directions(),
                    found: pops.len(),
                    field: "node populations".to_string(),
                });
            }
            // Residual of Σ (c − u_t) f = −F/2 over all populations.
            let mut residual = [0.0; 3];
            for (d, r) in residual.iter_mut().enumerate().take(dims) {
                let current: f64 = pops
                    .iter()
                    .zip(lattice.velocities())
                    .map(|(fk, c)| (f64::from(c[d]) - u_t[d]) * fk)
                    .sum();
                *r = -0.5 * force[d] - current;
            }
            // Constraint rows A[d][k] = c_k[d] − u_t[d] over unknown directions.
            let rows: Vec<[f64; 3]> = unknown
                .iter()
                .filter_map(|&k| lattice.velocity_f64(k))
                .map(|c| [c[0] - u_t[0], c[1] - u_t[1], c[2] - u_t[2]])
                .collect();
            let weights: Vec<f64> = unknown.iter().filter_map(|&k| lattice.weight(k)).collect();
            let mut m = [[0.0; 3]; 3];
            for (a_k, w_k) in rows.iter().zip(weights.iter()) {
                for a in 0..dims {
                    for b in 0..dims {
                        m[a][b] += w_k * a_k[a] * a_k[b];
                    }
                }
            }
            let lambda = symmetric_pinv_solve(&m, &residual, dims);
            let pops = f.get_mut(fluid_idx).ok_or(Error::IndexOutOfBounds {
                x: fluid_idx,
                y: 0,
                nx: f_post.len(),
                ny: 1,
            })?;
            for ((&k, a_k), w_k) in unknown.iter().zip(rows.iter()).zip(weights.iter()) {
                let delta: f64 = w_k * (0..dims).map(|a| a_k[a] * lambda[a]).sum::<f64>();
                if let Some(slot) = pops.get_mut(k) {
                    *slot += delta;
                }
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Everything one link needs from the post-collision field.
struct LinkContext {
    opp: usize,
    pops: LinkPopulations,
    link: WallLink,
    u_f: [f64; 3],
    u_ff: Option<[f64; 3]>,
    has_upstream2: bool,
}

impl LinkContext {
    fn gather(
        node: &CurvedWallNode,
        f_post: &[Vec<f64>],
        lattice: &CurvedLattice,
        force_at: &dyn Fn(usize) -> [f64; 3],
    ) -> Result<Self> {
        let nq = lattice.num_directions();
        let (c_i, w_i, opp) = match (
            lattice.velocity_f64(node.dir),
            lattice.weight(node.dir),
            lattice.opposite(node.dir),
        ) {
            (Some(c), Some(w), Some(o)) if node.dir != o => (c, w, o),
            _ => {
                return Err(Error::ParameterOutOfBounds {
                    name: "CurvedWallNode::dir".to_string(),
                    value: node.dir as f64,
                    min: 1.0,
                    max: nq.saturating_sub(1) as f64,
                });
            }
        };
        let node_pops = |idx: usize| -> Result<&[f64]> {
            let p = f_post.get(idx).ok_or(Error::IndexOutOfBounds {
                x: idx,
                y: 0,
                nx: f_post.len(),
                ny: 1,
            })?;
            if p.len() != nq {
                return Err(Error::ComponentCountMismatch {
                    expected: nq,
                    found: p.len(),
                    field: "node populations".to_string(),
                });
            }
            Ok(p.as_slice())
        };
        let state = |idx: usize| -> Result<NodeState> {
            let (rho, j) = lattice.moments(node_pops(idx)?)?;
            if !rho.is_finite() || rho <= 0.0 {
                return Err(Error::NonFiniteField {
                    field: "density".to_string(),
                    index: idx,
                    value: rho,
                });
            }
            let force = force_at(idx);
            Ok(NodeState {
                rho,
                u: [
                    (j[0] - 0.5 * force[0]) / rho,
                    (j[1] - 0.5 * force[1]) / rho,
                    (j[2] - 0.5 * force[2]) / rho,
                ],
            })
        };
        let xf = node_pops(node.fluid_idx)?;
        let fluid = state(node.fluid_idx)?;
        let mut pops = LinkPopulations {
            f_i_xf: xf[node.dir],
            f_opp_xf: xf[opp],
            ..LinkPopulations::default()
        };
        let mut u_ff = None;
        if let Some(up) = node.upstream {
            let p = node_pops(up)?;
            pops.f_i_xff = p[node.dir];
            pops.f_opp_xff = p[opp];
            u_ff = Some(state(up)?.u);
        }
        let mut has_upstream2 = false;
        if let (Some(_), Some(up2)) = (node.upstream, node.upstream2) {
            pops.f_i_xfff = node_pops(up2)?[node.dir];
            has_upstream2 = true;
        }
        let link = WallLink::new(node.q, w_i, c_i, node.u_wall, fluid.rho)?;
        Ok(Self {
            opp,
            pops,
            link,
            u_f: fluid.u,
            u_ff,
            has_upstream2,
        })
    }

    /// Reconstruct the unknown population with the given per-link scheme.
    fn reconstruct(&self, method: CurvedBcMethod, omega: f64) -> Result<f64> {
        let q = self.link.q;
        let has_up = self.u_ff.is_some();
        let linear = || {
            if q < 0.5 && !has_up {
                // No upstream node: half-way bounce-back with the Ladd term.
                self.pops.f_i_xf - self.link.moving_wall_term()
            } else {
                bouzidi_bc(&self.pops, &self.link)
            }
        };
        match method {
            CurvedBcMethod::Bouzidi | CurvedBcMethod::MomentBased => Ok(linear()),
            CurvedBcMethod::QuadraticInterpolation => {
                let complete = if q < 0.5 {
                    has_up && self.has_upstream2
                } else {
                    has_up
                };
                Ok(if complete {
                    bouzidi_quadratic_bc(&self.pops, &self.link)
                } else {
                    linear()
                })
            }
            CurvedBcMethod::FilippovaHanel => {
                if q < 0.5 && !has_up {
                    // The original Filippova–Hänel weighting (u_bf = u_f) is
                    // singular at τ = 1 and unstable near it; without an
                    // upstream node use the bounce-back fallback instead.
                    Ok(linear())
                } else {
                    filippova_hanel(self.pops.f_i_xf, &self.link, self.u_f, self.u_ff, omega)
                }
            }
        }
    }
}

/// Write `value` into `f[node][dir]` with bounds checking.
fn write_population(f: &mut [Vec<f64>], node: usize, dir: usize, value: f64) -> Result<()> {
    let len = f.len();
    let slot = f
        .get_mut(node)
        .and_then(|p| p.get_mut(dir))
        .ok_or(Error::IndexOutOfBounds {
            x: node,
            y: dir,
            nx: len,
            ny: 0,
        })?;
    *slot = value;
    Ok(())
}

/// Locate the wall crossing on the segment `pos → pos + step`.
///
/// Bisection on the SDF when the unwrapped end point is solid; otherwise
/// (periodic wrap onto a non-periodic SDF) linear interpolation of the
/// end-point values.  The result lies in `(0, 1]`.
fn link_fraction(
    sdf: &dyn Fn([f64; 3]) -> f64,
    pos: [f64; 3],
    step: [f64; 3],
    phi_f: f64,
    phi_s: f64,
) -> f64 {
    let at = |t: f64| {
        sdf([
            pos[0] + t * step[0],
            pos[1] + t * step[1],
            pos[2] + t * step[2],
        ])
    };
    let q = if at(1.0) <= 0.0 {
        let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if at(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    } else {
        phi_f / (phi_f - phi_s).max(f64::MIN_POSITIVE)
    };
    q.clamp(f64::EPSILON, 1.0)
}

/// Solve `M λ = r` for a symmetric positive semi-definite `n × n` block of `m`
/// (`n ≤ 3`) using the Moore–Penrose pseudo-inverse (cyclic Jacobi
/// eigen-decomposition, eigenvalues below `1e-12 · λ_max` discarded).
pub(super) fn symmetric_pinv_solve(m: &[[f64; 3]; 3], r: &[f64; 3], n: usize) -> [f64; 3] {
    let n = n.min(3);
    let mut a = *m;
    let mut v = [[0.0; 3]; 3];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    let pairs: Vec<(usize, usize)> = (0..n)
        .flat_map(|p| ((p + 1)..n).map(move |q| (p, q)))
        .collect();
    for _sweep in 0..50 {
        let off: f64 = pairs.iter().map(|&(p, q)| a[p][q] * a[p][q]).sum();
        if off < 1e-30 {
            break;
        }
        for &(p, q) in &pairs {
            let apq = a[p][q];
            if apq.abs() < 1e-300 {
                continue;
            }
            let theta = (a[q][q] - a[p][p]) / (2.0 * apq);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let c = 1.0 / (t * t + 1.0).sqrt();
            let s = t * c;
            // A ← A J (columns p, q)
            for row in a.iter_mut().take(n) {
                let (akp, akq) = (row[p], row[q]);
                row[p] = c * akp - s * akq;
                row[q] = s * akp + c * akq;
            }
            // A ← Jᵀ A (rows p, q)
            let (row_p, row_q) = (a[p], a[q]);
            a[p] = std::array::from_fn(|k| c * row_p[k] - s * row_q[k]);
            a[q] = std::array::from_fn(|k| s * row_p[k] + c * row_q[k]);
            // V ← V J
            for row in v.iter_mut().take(n) {
                let (vp, vq) = (row[p], row[q]);
                row[p] = c * vp - s * vq;
                row[q] = s * vp + c * vq;
            }
        }
    }
    let lambda_max = (0..n).map(|i| a[i][i].abs()).fold(0.0_f64, f64::max);
    let cutoff = 1e-12 * lambda_max.max(f64::MIN_POSITIVE);
    let mut out = [0.0; 3];
    for e in 0..n {
        let ev = a[e][e];
        if ev.abs() <= cutoff {
            continue;
        }
        let proj: f64 = (0..n).map(|k| v[k][e] * r[k]).sum::<f64>() / ev;
        for (k, o) in out.iter_mut().enumerate().take(n) {
            *o += v[k][e] * proj;
        }
    }
    out
}
