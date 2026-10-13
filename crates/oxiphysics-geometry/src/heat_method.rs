// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Heat-method geodesic distances on triangle meshes.
//!
//! This is the algorithm of Crane, Weischedel and Wardetzky, *Geodesics in Heat*
//! (2013), with the finite-element discretisation of the paper, shared by
//! [`HeatMethodSolver`](crate::discrete_geometry::HeatMethodSolver) and
//! [`heat_geodesic`](crate::geodesic::heat_geodesic):
//!
//! 1. **Heat flow.**  Integrate `du/dt = Δu` for a short time `t` from an
//!    impulse at the source vertices with backward Euler,
//!    `(M + dt·L) u⁺ = M u`, where `L` is the (positive semi-definite) cotangent
//!    stiffness matrix and `M` the lumped mixed-Voronoi vertex area matrix.  With
//!    `t` proportional to the squared edge length the result does not depend on
//!    the scale of the mesh.
//! 2. **Direction field.**  On every triangle `X = -∇u / |∇u|`, with the exact
//!    gradient of the piecewise-linear heat field,
//!    `∇u = (1 / 2A) Σᵢ uᵢ (N × eᵢ)` (`eᵢ` is the counter-clockwise edge opposite
//!    corner `i`).
//! 3. **Poisson solve.**  The distance `φ` is the piecewise-linear function whose
//!    gradient is closest to `X` in the least-squares sense: `L φ = b` with
//!    `bᵢ = Σ_f A_f X_f · ∇ψᵢ`, the integrated divergence of `X` over **both**
//!    triangles incident to each edge at `i`.
//!
//! The heat field decays geometrically away from the sources, so its far-field
//! values are tiny and must be computed with *relative* accuracy: the gradient
//! direction of a value that is merely accurate to the norm of the right-hand side
//! is noise.  The heat step therefore uses Gauss–Seidel sweeps in breadth-first
//! order from the sources, iterated until no value changes by more than a tiny
//! *relative* amount.  On Delaunay-like meshes every update only adds positive
//! terms, which keeps the full relative accuracy of the tail.  The Poisson step
//! needs only a norm-wise accurate solution and uses a Jacobi-preconditioned
//! conjugate-gradient solve.
//!
//! The heat underflows to zero beyond a few hundred mesh rings from the source; if
//! that happens before every vertex is reached the diffusion time is enlarged (see
//! [`solve`]).

use crate::error::{Error, Result};
use oxiphysics_core::sparse::{CsrMatrix, pcg_solve};

/// Tolerance of the relative change of any heat value between Gauss–Seidel sweeps.
const HEAT_SWEEP_TOLERANCE: f64 = 1e-13;
/// Gauss–Seidel sweeps (forward + backward) attempted per heat step before the
/// step falls back to a conjugate-gradient solve.
const HEAT_MAX_SWEEPS: usize = 4000;
/// Heat values below this (relative to the maximum) count as underflown.
const HEAT_UNDERFLOW: f64 = 1e-250;
/// How many times the diffusion time is enlarged when the heat underflows.
const HEAT_MAX_ENLARGEMENTS: usize = 8;
/// Factor by which the diffusion time is enlarged on underflow.
const HEAT_ENLARGEMENT: f64 = 4.0;
/// Tolerance of the preconditioned conjugate-gradient solve (relative to ‖b‖).
const PCG_TOLERANCE: f64 = 1e-10;
/// Largest relative residual accepted from the conjugate-gradient solve.
const PCG_ACCEPTED_RESIDUAL: f64 = 1e-6;

/// Output of [`solve`].
#[derive(Debug, Clone)]
pub(crate) struct HeatSolution {
    /// Solution `φ` of the Poisson problem; defined up to one additive constant
    /// per connected component.
    pub(crate) potential: Vec<f64>,
    /// Connected component of every vertex (through non-degenerate triangles).
    pub(crate) component: Vec<usize>,
    /// Number of connected components.
    pub(crate) component_count: usize,
    /// Unit direction field `X` per triangle (zero on degenerate triangles and
    /// where the heat gradient vanishes).
    pub(crate) directions: Vec<[f64; 3]>,
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Geometry of one non-degenerate triangle.
struct FaceGeometry {
    /// Unit normal (right-hand rule on the vertex order).
    normal: [f64; 3],
    /// Area.
    area: f64,
    /// Counter-clockwise edge opposite each corner: from corner `c + 1` to `c + 2`.
    opposite: [[f64; 3]; 3],
    /// Cotangent of the interior angle at each corner.
    cot: [f64; 3],
}

impl FaceGeometry {
    /// `None` for a (numerically) degenerate triangle.
    fn new(vertices: &[[f64; 3]], tri: &[usize; 3]) -> Option<Self> {
        let p = [vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]];
        let opposite = [sub(p[2], p[1]), sub(p[0], p[2]), sub(p[1], p[0])];
        let n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
        let twice_area = dot(n, n).sqrt();
        let longest_sq = opposite.iter().map(|e| dot(*e, *e)).fold(0.0_f64, f64::max);
        let usable = twice_area.is_finite() && twice_area > 1e-14 * longest_sq;
        if !usable {
            return None;
        }
        let corner_cot = |c: usize| -> f64 {
            let next = sub(p[(c + 1) % 3], p[c]);
            let prev = sub(p[(c + 2) % 3], p[c]);
            dot(next, prev) / twice_area
        };
        Some(Self {
            normal: [n[0] / twice_area, n[1] / twice_area, n[2] / twice_area],
            area: 0.5 * twice_area,
            opposite,
            cot: [corner_cot(0), corner_cot(1), corner_cot(2)],
        })
    }
}

/// Connected components of the vertices through non-degenerate triangles.
struct Components {
    /// Component index of every vertex (numbered by first appearance).
    id: Vec<usize>,
    /// Number of components.
    count: usize,
}

impl Components {
    fn label(n: usize, faces: &[[usize; 3]], geometry: &[Option<FaceGeometry>]) -> Self {
        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        let mut parent: Vec<usize> = (0..n).collect();
        for (tri, geom) in faces.iter().zip(geometry.iter()) {
            if geom.is_none() {
                continue;
            }
            for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2])] {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    parent[ra.max(rb)] = ra.min(rb);
                }
            }
        }
        let mut label_of_root: Vec<usize> = vec![usize::MAX; n];
        let mut id = vec![0usize; n];
        let mut count = 0usize;
        for v in 0..n {
            let root = find(&mut parent, v);
            if label_of_root[root] == usize::MAX {
                label_of_root[root] = count;
                count += 1;
            }
            id[v] = label_of_root[root];
        }
        Self { id, count }
    }
}

/// Lumped mixed-Voronoi vertex areas (Meyer et al. 2003).  A vertex that carries
/// no area (isolated, or only touched by degenerate triangles) gets the mean
/// area so that every diagonal entry of the heat system stays positive.
fn lumped_masses(n: usize, faces: &[[usize; 3]], geometry: &[Option<FaceGeometry>]) -> Vec<f64> {
    let mut mass = vec![0.0_f64; n];
    for (tri, geom) in faces.iter().zip(geometry.iter()) {
        let Some(g) = geom else {
            continue;
        };
        if let Some(obtuse) = (0..3).find(|&c| g.cot[c] < 0.0) {
            for (c, &v) in tri.iter().enumerate() {
                mass[v] += if c == obtuse {
                    0.5 * g.area
                } else {
                    0.25 * g.area
                };
            }
        } else {
            for (c, &v) in tri.iter().enumerate() {
                let mut voronoi = 0.0;
                for d in (0..3).filter(|&d| d != c) {
                    voronoi += dot(g.opposite[d], g.opposite[d]) * g.cot[d];
                }
                mass[v] += 0.125 * voronoi;
            }
        }
    }
    let positive: Vec<f64> = mass.iter().copied().filter(|&m| m > 0.0).collect();
    let fallback = if positive.is_empty() {
        1.0
    } else {
        positive.iter().sum::<f64>() / positive.len() as f64
    };
    for m in &mut mass {
        if m.is_nan() || *m <= 0.0 {
            *m = fallback;
        }
    }
    mass
}

/// Cotangent Laplacian as vertex adjacency with edge weights.
///
/// With `w_ij = ½ (cot α + cot β)` the stiffness matrix `L` has
/// `L_ij = -w_ij` for neighbours and `L_ii = Σ_j w_ij` (positive semi-definite).
struct LaplaceGraph {
    /// Row pointers of the neighbour lists.
    ptr: Vec<usize>,
    /// Neighbour vertex of every entry.
    idx: Vec<usize>,
    /// Cotangent weight of every entry.
    weight: Vec<f64>,
    /// `L_ii`.
    diag: Vec<f64>,
}

impl LaplaceGraph {
    fn build(n: usize, faces: &[[usize; 3]], geometry: &[Option<FaceGeometry>]) -> Self {
        let mut edges: Vec<(usize, usize, f64)> = Vec::with_capacity(faces.len() * 3);
        for (tri, geom) in faces.iter().zip(geometry.iter()) {
            let Some(g) = geom else {
                continue;
            };
            for c in 0..3 {
                let (a, b) = (tri[(c + 1) % 3], tri[(c + 2) % 3]);
                edges.push((a.min(b), a.max(b), 0.5 * g.cot[c]));
            }
        }
        // A stable sort keeps the summation order, hence the result, deterministic.
        edges.sort_by_key(|&(a, b, _)| (a, b));
        let mut merged: Vec<(usize, usize, f64)> = Vec::with_capacity(edges.len() / 2 + 1);
        for (a, b, w) in edges {
            match merged.last_mut() {
                Some(last) if last.0 == a && last.1 == b => last.2 += w,
                _ => merged.push((a, b, w)),
            }
        }
        let mut degree = vec![0usize; n];
        for &(a, b, _) in &merged {
            degree[a] += 1;
            degree[b] += 1;
        }
        let mut ptr = vec![0usize; n + 1];
        for v in 0..n {
            ptr[v + 1] = ptr[v] + degree[v];
        }
        let mut next = ptr.clone();
        let mut idx = vec![0usize; ptr[n]];
        let mut weight = vec![0.0_f64; ptr[n]];
        let mut diag = vec![0.0_f64; n];
        for &(a, b, w) in &merged {
            idx[next[a]] = b;
            weight[next[a]] = w;
            next[a] += 1;
            idx[next[b]] = a;
            weight[next[b]] = w;
            next[b] += 1;
            diag[a] += w;
            diag[b] += w;
        }
        Self {
            ptr,
            idx,
            weight,
            diag,
        }
    }

    fn vertex_count(&self) -> usize {
        self.diag.len()
    }

    /// Breadth-first order of the vertices reachable from `sources`.
    fn breadth_first_order(&self, sources: &[usize]) -> Vec<usize> {
        let n = self.vertex_count();
        let mut seen = vec![false; n];
        let mut order: Vec<usize> = Vec::with_capacity(n);
        for &s in sources {
            if !seen[s] {
                seen[s] = true;
                order.push(s);
            }
        }
        let mut head = 0;
        while head < order.len() {
            let v = order[head];
            head += 1;
            for k in self.ptr[v]..self.ptr[v + 1] {
                let w = self.idx[k];
                if !seen[w] {
                    seen[w] = true;
                    order.push(w);
                }
            }
        }
        order
    }

    /// Solve `(M + dt·L) u⁺ = rhs` in place by symmetric Gauss–Seidel sweeps in
    /// `order` until no heat value changes by more than [`HEAT_SWEEP_TOLERANCE`]
    /// (relative).  `u` holds the initial guess.
    ///
    /// For a long diffusion time the system is dominated by the Laplacian and
    /// Gauss–Seidel converges slowly; when `max_sweeps` sweeps are not enough the
    /// step falls back to a norm-wise accurate conjugate-gradient solve (the heat
    /// then spreads so far that its tail is not many orders of magnitude small).
    fn heat_step(
        &self,
        mass: &[f64],
        dt: f64,
        rhs: &[f64],
        u: &mut [f64],
        order: &[usize],
        max_sweeps: usize,
    ) -> Result<()> {
        let relax = |u: &mut [f64], v: usize, worst: &mut f64| {
            let mut acc = rhs[v];
            for k in self.ptr[v]..self.ptr[v + 1] {
                acc += dt * self.weight[k] * u[self.idx[k]];
            }
            let new = acc / (mass[v] + dt * self.diag[v]);
            let old = u[v];
            let scale = new.abs().max(old.abs());
            if scale > 1e-290 {
                *worst = worst.max((new - old).abs() / scale);
            }
            u[v] = new;
        };
        for _ in 0..max_sweeps {
            let mut worst = 0.0_f64;
            for &v in order {
                relax(u, v, &mut worst);
            }
            for &v in order.iter().rev() {
                relax(u, v, &mut worst);
            }
            if worst <= HEAT_SWEEP_TOLERANCE {
                return Ok(());
            }
        }
        let solved = self
            .scaled_system(mass, dt)
            .solve(rhs, "heat method diffusion")?;
        u.copy_from_slice(&solved);
        Ok(())
    }

    /// Backward-Euler heat flow from an impulse at `sources` for a total time `t`
    /// in `steps` equal sub-steps.  The result is scaled to a maximum of one.
    fn heat_flow(
        &self,
        mass: &[f64],
        sources: &[usize],
        t: f64,
        steps: usize,
        order: &[usize],
    ) -> Result<Vec<f64>> {
        let n = self.vertex_count();
        let dt = t / steps as f64;
        let mut u = vec![0.0_f64; n];
        for &s in sources {
            u[s] = 1.0 / mass[s];
        }
        for _ in 0..steps {
            let rhs: Vec<f64> = mass.iter().zip(u.iter()).map(|(m, h)| m * h).collect();
            self.heat_step(mass, dt, &rhs, &mut u, order, HEAT_MAX_SWEEPS)?;
            let peak = u.iter().copied().fold(0.0_f64, f64::max);
            if peak > 0.0 && peak.is_finite() {
                for h in &mut u {
                    *h /= peak;
                }
            }
        }
        Ok(u)
    }

    /// [`heat_flow`](Self::heat_flow) with the diffusion time enlarged by
    /// [`HEAT_ENLARGEMENT`] — at most `max_enlargements` times — until the heat has
    /// not underflowed at any vertex of `order` (all vertices reachable from the
    /// sources).
    fn heat_flow_reaching_all(
        &self,
        mass: &[f64],
        sources: &[usize],
        t: f64,
        steps: usize,
        order: &[usize],
        max_enlargements: usize,
    ) -> Result<Vec<f64>> {
        let mut t = t;
        let mut heat = self.heat_flow(mass, sources, t, steps, order)?;
        let mut enlargements = 0usize;
        while !order.iter().all(|&v| heat[v].abs() > HEAT_UNDERFLOW) {
            if enlargements == max_enlargements {
                return Err(Error::general(format!(
                    "heat method: the heat underflows before it reaches every vertex even \
                     for a diffusion time of {t:e}; the mesh is too large for double precision"
                )));
            }
            enlargements += 1;
            t *= HEAT_ENLARGEMENT;
            heat = self.heat_flow(mass, sources, t, steps, order)?;
        }
        Ok(heat)
    }

    /// The matrix `diag(shift) + coefficient · L`, symmetrically Jacobi-scaled.
    fn scaled_system(&self, shift: &[f64], coefficient: f64) -> ScaledSystem {
        let n = self.vertex_count();
        let diagonal: Vec<f64> = (0..n)
            .map(|v| shift[v] + coefficient * self.diag[v])
            .collect();
        let scale: Vec<f64> = diagonal
            .iter()
            .map(|&d| if d > 1e-300 { 1.0 / d.sqrt() } else { 1.0 })
            .collect();
        let mut rows = Vec::with_capacity(self.idx.len() + n);
        let mut cols = Vec::with_capacity(self.idx.len() + n);
        let mut vals = Vec::with_capacity(self.idx.len() + n);
        for v in 0..n {
            if diagonal[v] > 1e-300 {
                rows.push(v);
                cols.push(v);
                vals.push(1.0);
            }
            for k in self.ptr[v]..self.ptr[v + 1] {
                let w = self.idx[k];
                rows.push(v);
                cols.push(w);
                vals.push(-coefficient * self.weight[k] * scale[v] * scale[w]);
            }
        }
        ScaledSystem {
            matrix: CsrMatrix::from_triplets(n, n, &rows, &cols, &vals),
            scale,
        }
    }

    /// Solve the (singular) Poisson system `L φ = b`.  `b` must already sum to zero
    /// over every connected component.
    fn poisson(&self, rhs: &[f64]) -> Result<Vec<f64>> {
        let no_shift = vec![0.0_f64; self.vertex_count()];
        self.scaled_system(&no_shift, 1.0)
            .solve(rhs, "heat method Poisson solve")
    }
}

/// A symmetric system `A x = b` scaled to a unit diagonal, `Â = S A S` with
/// `S = diag(A)^-1/2`, for the conjugate-gradient routine.
///
/// Conjugate gradients on `Â` generate exactly the iterates of Jacobi-preconditioned
/// conjugate gradients on `A`, but the absolute thresholds of the routine now act
/// on a well-scaled problem whatever the size of the mesh.
struct ScaledSystem {
    /// `S A S`.
    matrix: CsrMatrix,
    /// The diagonal of `S`.
    scale: Vec<f64>,
}

impl ScaledSystem {
    /// Solve `A x = rhs` (`A` may be singular if `rhs` is consistent).
    fn solve(&self, rhs: &[f64], operation: &'static str) -> Result<Vec<f64>> {
        let n = rhs.len();
        let scaled_rhs: Vec<f64> = rhs
            .iter()
            .zip(self.scale.iter())
            .map(|(b, s)| b * s)
            .collect();
        let norm = scaled_rhs.iter().map(|b| b * b).sum::<f64>().sqrt();
        if !norm.is_finite() {
            return Err(Error::general(format!(
                "{operation}: non-finite right-hand side"
            )));
        }
        if norm == 0.0 {
            return Ok(vec![0.0; n]);
        }
        let unit_rhs: Vec<f64> = scaled_rhs.iter().map(|b| b / norm).collect();
        // Conjugate gradients terminate in at most `n` steps in exact arithmetic; the
        // factor two leaves room for rounding on long, thin meshes where the number of
        // iterations is proportional to the length of the mesh in edges.
        let max_iter = 2 * n + 1000;
        let unit_preconditioner = vec![1.0_f64; n];
        let (x_hat, iterations, _) = pcg_solve(
            &self.matrix,
            &unit_rhs,
            &unit_preconditioner,
            PCG_TOLERANCE,
            max_iter,
        );
        let applied = self.matrix.matvec(&x_hat);
        let residual = unit_rhs
            .iter()
            .zip(applied.iter())
            .map(|(b, a)| (b - a) * (b - a))
            .sum::<f64>()
            .sqrt();
        if residual.is_nan() || residual > PCG_ACCEPTED_RESIDUAL {
            return Err(Error::convergence_failure(
                operation,
                iterations as usize,
                residual,
            ));
        }
        Ok(x_hat
            .iter()
            .zip(self.scale.iter())
            .map(|(x, s)| x * norm * s)
            .collect())
    }
}

/// Reject malformed input before any indexing happens.
fn validate(vertices: &[[f64; 3]], faces: &[[usize; 3]], sources: &[usize]) -> Result<()> {
    let n = vertices.len();
    for (i, v) in vertices.iter().enumerate() {
        if !v.iter().all(|c| c.is_finite()) {
            return Err(Error::invalid_mesh(format!(
                "vertex {i} has a non-finite coordinate"
            )));
        }
    }
    for (f, tri) in faces.iter().enumerate() {
        if let Some(&bad) = tri.iter().find(|&&v| v >= n) {
            return Err(Error::invalid_mesh(format!(
                "triangle {f} references vertex {bad} but the mesh has {n} vertices"
            )));
        }
    }
    for &s in sources {
        if s >= n {
            return Err(Error::index_out_of_bounds(s, n));
        }
    }
    Ok(())
}

/// Run the heat method on the mesh `(vertices, faces)` for the given `sources`.
///
/// `diffusion_time` is the total heat-flow time `t` (it must be positive and
/// finite; a good choice is the squared mean edge length) and `steps` the number
/// of equal backward-Euler sub-steps it is split into (at least one).  When the
/// heat underflows before it reaches every vertex of a source's component — only
/// possible for meshes that are hundreds of edges across — the diffusion time is
/// enlarged by a factor four, at most eight times, until it does; if it still
/// underflows after that an error is returned.
///
/// # Errors
/// * [`Error::InvalidMesh`] for a non-finite vertex or a triangle that references a
///   vertex outside the mesh;
/// * [`Error::IndexOutOfBounds`] for a source outside the mesh;
/// * [`Error::OutOfRange`] for a non-positive or non-finite `diffusion_time`;
/// * [`Error::ConvergenceFailure`] when an iterative solve does not converge;
/// * [`Error::General`] when the heat underflows for every admissible diffusion time.
pub(crate) fn solve(
    vertices: &[[f64; 3]],
    faces: &[[usize; 3]],
    sources: &[usize],
    diffusion_time: f64,
    steps: usize,
) -> Result<HeatSolution> {
    validate(vertices, faces, sources)?;
    let usable_time = diffusion_time.is_finite() && diffusion_time > 0.0;
    if !usable_time {
        return Err(Error::out_of_range(
            "diffusion_time",
            diffusion_time,
            f64::MIN_POSITIVE,
            f64::MAX,
        ));
    }
    let n = vertices.len();
    let geometry: Vec<Option<FaceGeometry>> = faces
        .iter()
        .map(|tri| FaceGeometry::new(vertices, tri))
        .collect();
    let components = Components::label(n, faces, &geometry);
    let mut directions = vec![[0.0_f64; 3]; faces.len()];
    if n == 0 {
        return Ok(HeatSolution {
            potential: Vec::new(),
            component: Vec::new(),
            component_count: 0,
            directions,
        });
    }
    let mass = lumped_masses(n, faces, &geometry);
    let graph = LaplaceGraph::build(n, faces, &geometry);
    let mut unique_sources: Vec<usize> = sources.to_vec();
    unique_sources.sort_unstable();
    unique_sources.dedup();
    let order = graph.breadth_first_order(&unique_sources);

    // 1. Heat flow, enlarging the time while the tail underflows.
    let heat = graph.heat_flow_reaching_all(
        &mass,
        &unique_sources,
        diffusion_time,
        steps.max(1),
        &order,
        HEAT_MAX_ENLARGEMENTS,
    )?;

    // 2. Unit direction field and its integrated divergence.
    let mut rhs = vec![0.0_f64; n];
    for (f, tri) in faces.iter().enumerate() {
        let Some(g) = &geometry[f] else {
            continue;
        };
        let rotated: [[f64; 3]; 3] = std::array::from_fn(|c| cross(g.normal, g.opposite[c]));
        let mut grad = [0.0_f64; 3];
        for (&v, r) in tri.iter().zip(rotated.iter()) {
            let u = heat[v];
            grad[0] += u * r[0];
            grad[1] += u * r[1];
            grad[2] += u * r[2];
        }
        // grad is 2A times the gradient of the heat field; only its direction matters.
        let len = dot(grad, grad).sqrt();
        let usable_gradient = len > 0.0 && len.is_finite();
        if !usable_gradient {
            continue;
        }
        let x = [-grad[0] / len, -grad[1] / len, -grad[2] / len];
        directions[f] = x;
        for (&v, r) in tri.iter().zip(rotated.iter()) {
            rhs[v] += 0.5 * dot(x, *r);
        }
    }

    // 3. Poisson solve; the right-hand side must sum to zero on every component.
    let mut sums = vec![0.0_f64; components.count];
    let mut counts = vec![0usize; components.count];
    for (v, &b) in rhs.iter().enumerate() {
        sums[components.id[v]] += b;
        counts[components.id[v]] += 1;
    }
    for (v, b) in rhs.iter_mut().enumerate() {
        let c = components.id[v];
        *b -= sums[c] / counts[c] as f64;
    }
    let potential = graph.poisson(&rhs)?;
    Ok(HeatSolution {
        potential,
        component: components.id,
        component_count: components.count,
        directions,
    })
}

/// Mean and mean square of the triangle side lengths, counted once per triangle
/// side.  Triangles that reference a vertex outside the mesh are skipped; `None`
/// when no side remains.
pub(crate) fn edge_length_statistics(
    vertices: &[[f64; 3]],
    faces: &[[usize; 3]],
) -> Option<(f64, f64)> {
    let mut sum = 0.0_f64;
    let mut sum_sq = 0.0_f64;
    let mut count = 0usize;
    for tri in faces {
        if tri.iter().any(|&v| v >= vertices.len()) {
            continue;
        }
        for c in 0..3 {
            let d = sub(vertices[tri[(c + 1) % 3]], vertices[tri[c]]);
            let len_sq = dot(d, d);
            sum += len_sq.sqrt();
            sum_sq += len_sq;
            count += 1;
        }
    }
    let usable = count > 0 && sum.is_finite() && sum_sq.is_finite() && sum > 0.0;
    usable.then_some((sum / count as f64, sum_sq / count as f64))
}

/// Shift the Poisson potential into distances.
///
/// Every connected component that contains a source is shifted so that its
/// smallest potential is zero; vertices in components without a source are
/// unreachable and get `f64::INFINITY`.
pub(crate) fn distances_from_min_shift(solution: &HeatSolution, sources: &[usize]) -> Vec<f64> {
    let mut has_source = vec![false; solution.component_count];
    for &s in sources {
        has_source[solution.component[s]] = true;
    }
    let mut minimum = vec![f64::INFINITY; solution.component_count];
    for (v, &phi) in solution.potential.iter().enumerate() {
        let c = solution.component[v];
        if has_source[c] {
            minimum[c] = minimum[c].min(phi);
        }
    }
    solution
        .potential
        .iter()
        .enumerate()
        .map(|(v, &phi)| {
            let c = solution.component[v];
            if has_source[c] {
                phi - minimum[c]
            } else {
                f64::INFINITY
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_meshes::{RawMesh, flat_grid, icosphere};
    use std::time::Instant;

    fn heat_distances(mesh: &RawMesh, sources: &[usize], t: f64, steps: usize) -> Vec<f64> {
        let solution = solve(&mesh.0, &mesh.1, sources, t, steps).expect("valid mesh");
        distances_from_min_shift(&solution, sources)
    }

    fn mean_edge(mesh: &RawMesh) -> f64 {
        edge_length_statistics(&mesh.0, &mesh.1)
            .map(|(mean, _)| mean)
            .expect("mesh has edges")
    }

    fn euclid(a: [f64; 3], b: [f64; 3]) -> f64 {
        dot(sub(a, b), sub(a, b)).sqrt()
    }

    /// Great-circle distance between two points on the unit sphere.
    fn great_circle(a: [f64; 3], b: [f64; 3]) -> f64 {
        dot(a, b).clamp(-1.0, 1.0).acos()
    }

    #[test]
    fn flat_grid_distances_match_the_euclidean_distance() {
        let mesh = flat_grid(20, 20, 1.0);
        let source = 10 * 21 + 10;
        let h = mean_edge(&mesh);
        let d = heat_distances(&mesh, &[source], h * h, 1);
        assert!(d[source].abs() < 1e-9, "source distance {}", d[source]);
        let (mut sum, mut count, mut worst) = (0.0_f64, 0usize, 0.0_f64);
        for (v, &dv) in d.iter().enumerate() {
            let exact = euclid(mesh.0[v], mesh.0[source]);
            if exact < 3.0 {
                continue;
            }
            let rel = (dv - exact).abs() / exact;
            sum += rel;
            count += 1;
            worst = worst.max(rel);
        }
        // The heat method is a few-percent method: ~3 % on average and below 9 %
        // everywhere on this anisotropic right-triangle lattice.
        let mean = sum / count as f64;
        assert!(mean < 0.04, "mean relative error {mean}");
        assert!(worst < 0.10, "worst relative error {worst}");
        // Along the lattice diagonal the distance is accurate to a percent or two.
        let at = |i: usize, j: usize| j * 21 + i;
        for k in 3..9 {
            let exact = (k as f64) * 2.0_f64.sqrt();
            let diagonal = d[at(10 + k, 10 + k)];
            assert!(
                (diagonal - exact).abs() < 0.03 * exact,
                "k = {k}: {diagonal} vs {exact}"
            );
        }
    }

    #[test]
    fn sphere_antipode_distance_is_pi() {
        let mesh = icosphere(4);
        assert_eq!(mesh.0.len(), 2562);
        let h = mean_edge(&mesh);
        for steps in [1, 4] {
            let d = heat_distances(&mesh, &[0], h * h, steps);
            // Vertex 3 of the icosahedron is the antipode of vertex 0.
            let antipode = d[3];
            assert!(
                (antipode - std::f64::consts::PI).abs() < 0.03 * std::f64::consts::PI,
                "steps = {steps}: antipode distance {antipode}"
            );
            let mut worst = 0.0_f64;
            for (v, &dv) in d.iter().enumerate() {
                let exact = great_circle(mesh.0[v], mesh.0[0]);
                if exact > 0.3 {
                    worst = worst.max((dv - exact).abs() / exact);
                }
            }
            assert!(
                worst < 0.06,
                "steps = {steps}: worst relative error {worst}"
            );
        }
    }

    #[test]
    fn distances_scale_with_the_mesh() {
        let mesh = icosphere(2);
        let h = mean_edge(&mesh);
        let reference = heat_distances(&mesh, &[5], h * h, 1);
        let s = 7.5;
        let scaled: RawMesh = (
            mesh.0
                .iter()
                .map(|p| [p[0] * s, p[1] * s, p[2] * s])
                .collect(),
            mesh.1.clone(),
        );
        let big = heat_distances(&scaled, &[5], (h * s) * (h * s), 1);
        for (a, b) in reference.iter().zip(big.iter()) {
            assert!(
                (a * s - b).abs() < 1e-6 * (1.0 + b.abs()),
                "{a} * {s} vs {b}"
            );
        }
    }

    #[test]
    fn far_field_of_a_long_strip_stays_accurate() {
        // 400 cells long: the heat has to be resolved over hundreds of rings.
        let mesh = flat_grid(400, 2, 1.0);
        let h = mean_edge(&mesh);
        let d = heat_distances(&mesh, &[0], h * h, 1);
        let mut worst = 0.0_f64;
        for (v, &dv) in d.iter().enumerate() {
            let exact = euclid(mesh.0[v], mesh.0[0]);
            if exact > 5.0 {
                worst = worst.max((dv - exact).abs() / exact);
            }
        }
        assert!(worst < 0.06, "worst relative error {worst}");
        assert!(
            d[400] > 390.0 && d[400] < 410.0,
            "last vertex of the first row: {}",
            d[400]
        );
    }

    #[test]
    fn multiple_sources_give_the_distance_to_the_nearest_one() {
        let mesh = flat_grid(30, 30, 1.0);
        let sources = [0usize, 30, 31 * 30, 31 * 31 - 1];
        let h = mean_edge(&mesh);
        let d = heat_distances(&mesh, &sources, h * h, 1);
        for &s in &sources {
            assert!(d[s] < 1.5, "source {s} has distance {}", d[s]);
        }
        let centre = 15 * 31 + 15;
        let exact = sources
            .iter()
            .map(|&s| euclid(mesh.0[centre], mesh.0[s]))
            .fold(f64::INFINITY, f64::min);
        assert!(
            (d[centre] - exact).abs() < 0.1 * exact,
            "centre {} vs {exact}",
            d[centre]
        );
    }

    #[test]
    fn components_without_a_source_are_unreachable() {
        let (mut vertices, mut faces) = flat_grid(6, 6, 1.0);
        let (more_vertices, more_faces) = flat_grid(6, 6, 1.0);
        let offset = vertices.len();
        vertices.extend(more_vertices.iter().map(|p| [p[0] + 100.0, p[1], p[2]]));
        faces.extend(
            more_faces
                .iter()
                .map(|t| [t[0] + offset, t[1] + offset, t[2] + offset]),
        );
        let mesh: RawMesh = (vertices, faces);
        let h = mean_edge(&mesh);

        let single = heat_distances(&mesh, &[3], h * h, 1);
        assert!(single[..offset].iter().all(|d| d.is_finite()));
        assert!(single[offset..].iter().all(|d| d.is_infinite()));

        let both = heat_distances(&mesh, &[3, offset + 3], h * h, 1);
        assert!(both.iter().all(|d| d.is_finite()));
        assert!(both[3] < 1e-6 && both[offset + 3] < 1e-6 + 1.0);
    }

    #[test]
    fn invalid_input_is_reported() {
        let mesh = flat_grid(3, 3, 1.0);
        let n = mesh.0.len();
        match solve(&mesh.0, &mesh.1, &[n], 1.0, 1) {
            Err(Error::IndexOutOfBounds { index, len }) => assert_eq!((index, len), (n, n)),
            other => panic!("expected IndexOutOfBounds, got {other:?}"),
        }
        let mut bad_faces = mesh.1.clone();
        bad_faces[2] = [0, 1, n + 4];
        assert!(matches!(
            solve(&mesh.0, &bad_faces, &[0], 1.0, 1),
            Err(Error::InvalidMesh { .. })
        ));
        for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(matches!(
                solve(&mesh.0, &mesh.1, &[0], t, 1),
                Err(Error::OutOfRange { .. })
            ));
        }
        let mut nan_vertices = mesh.0.clone();
        nan_vertices[4][1] = f64::NAN;
        assert!(matches!(
            solve(&nan_vertices, &mesh.1, &[0], 1.0, 1),
            Err(Error::InvalidMesh { .. })
        ));
        // An empty mesh has no valid source but is otherwise fine.
        assert!(solve(&[], &[], &[], 1.0, 1).is_ok());
        assert!(matches!(
            solve(&[], &[], &[0], 1.0, 1),
            Err(Error::IndexOutOfBounds { .. })
        ));
    }

    #[test]
    fn degenerate_triangles_and_isolated_vertices_are_tolerated() {
        let (mut vertices, mut faces) = flat_grid(5, 5, 1.0);
        let isolated = vertices.len();
        vertices.push([50.0, 50.0, 0.0]);
        // A zero-area triangle on three collinear vertices.
        faces.push([0, 1, 2]);
        let mesh: RawMesh = (vertices, faces);
        let h = mean_edge(&mesh);
        let d = heat_distances(&mesh, &[14], h * h, 1);
        assert!(d[..isolated].iter().all(|x| x.is_finite()));
        assert!(d[isolated].is_infinite());
        let solution = solve(&mesh.0, &mesh.1, &[14], h * h, 1).expect("valid mesh");
        assert_eq!(solution.directions.last(), Some(&[0.0, 0.0, 0.0]));
    }

    #[test]
    fn direction_field_is_unit_length_and_points_away_from_the_source() {
        let mesh = flat_grid(16, 16, 1.0);
        let source = 8 * 17 + 8;
        let h = mean_edge(&mesh);
        let solution = solve(&mesh.0, &mesh.1, &[source], h * h, 1).expect("valid mesh");
        let mut interior = 0;
        for (f, x) in solution.directions.iter().enumerate() {
            let len = dot(*x, *x).sqrt();
            assert!((len - 1.0).abs() < 1e-9, "face {f}: |X| = {len}");
            assert!(x[2].abs() < 1e-9, "face {f}: X leaves the plane");
            let tri = mesh.1[f];
            let centroid = [
                (mesh.0[tri[0]][0] + mesh.0[tri[1]][0] + mesh.0[tri[2]][0]) / 3.0,
                (mesh.0[tri[0]][1] + mesh.0[tri[1]][1] + mesh.0[tri[2]][1]) / 3.0,
                0.0,
            ];
            let radial = sub(centroid, mesh.0[source]);
            let r = dot(radial, radial).sqrt();
            if r <= 2.0 {
                continue;
            }
            let cos = dot(*x, radial) / r;
            // Near the domain boundary the (Neumann) heat flow bends along the
            // boundary, so the direction field is only roughly radial there.
            assert!(
                cos > 0.7,
                "face {f}: X is not outward ({x:?} vs {radial:?})"
            );
            let to_boundary = centroid[0]
                .min(16.0 - centroid[0])
                .min(centroid[1])
                .min(16.0 - centroid[1]);
            if to_boundary >= 3.0 {
                interior += 1;
                assert!(cos > 0.88, "face {f}: X is not radial (cos = {cos})");
            }
        }
        assert!(
            interior > 100,
            "only {interior} interior faces were checked"
        );
    }

    #[test]
    fn long_diffusion_falls_back_to_the_conjugate_gradient_solve() {
        let (vertices, faces) = flat_grid(12, 12, 1.0);
        let geometry: Vec<Option<FaceGeometry>> = faces
            .iter()
            .map(|t| FaceGeometry::new(&vertices, t))
            .collect();
        let n = vertices.len();
        let mass = lumped_masses(n, &faces, &geometry);
        let graph = LaplaceGraph::build(n, &faces, &geometry);
        let source = 6 * 13 + 6;
        let order = graph.breadth_first_order(&[source]);
        let mut rhs = vec![0.0_f64; n];
        rhs[source] = 1.0;
        let dt = 300.0;
        let mut sweeps = vec![0.0_f64; n];
        graph
            .heat_step(&mass, dt, &rhs, &mut sweeps, &order, 100_000)
            .expect("Gauss-Seidel converges");
        // A single sweep cannot converge, which forces the conjugate-gradient path.
        let mut cg = vec![0.0_f64; n];
        graph
            .heat_step(&mass, dt, &rhs, &mut cg, &order, 1)
            .expect("conjugate-gradient fallback");
        let peak = sweeps.iter().copied().fold(0.0_f64, f64::max);
        assert!(peak > 0.0);
        for (v, (a, b)) in sweeps.iter().zip(cg.iter()).enumerate() {
            assert!(
                (a - b).abs() <= 1e-6 * peak,
                "vertex {v}: Gauss-Seidel {a} vs conjugate gradients {b}"
            );
        }
        // The whole pipeline accepts such a long diffusion time as well.
        let mesh: RawMesh = (vertices, faces);
        let d = heat_distances(&mesh, &[source], dt, 1);
        assert!(d.iter().all(|x| x.is_finite() && *x >= 0.0));
        assert!(d[source].abs() < 1e-6);
    }

    #[test]
    fn heat_that_underflows_is_given_more_time() {
        // 1500 cells long: with the default diffusion time the heat falls below
        // 1e-250 long before the far end, so the time has to be enlarged.
        let mesh = flat_grid(1500, 1, 1.0);
        let h = mean_edge(&mesh);
        let d = heat_distances(&mesh, &[0], h * h, 1);
        assert!(d.iter().all(|x| x.is_finite()));
        // The strip is effectively one-dimensional, where the method is exact.
        assert!(
            (d[1500] - 1500.0).abs() < 0.01 * 1500.0,
            "far end at {}",
            d[1500]
        );
    }

    #[test]
    fn underflow_without_enlargements_is_an_error() {
        let (vertices, faces) = flat_grid(1500, 1, 1.0);
        let geometry: Vec<Option<FaceGeometry>> = faces
            .iter()
            .map(|t| FaceGeometry::new(&vertices, t))
            .collect();
        let n = vertices.len();
        let mass = lumped_masses(n, &faces, &geometry);
        let graph = LaplaceGraph::build(n, &faces, &geometry);
        let order = graph.breadth_first_order(&[0]);
        let t = 1.3;
        assert!(matches!(
            graph.heat_flow_reaching_all(&mass, &[0], t, 1, &order, 0),
            Err(Error::General(_))
        ));
        let heat = graph
            .heat_flow_reaching_all(&mass, &[0], t, 1, &order, HEAT_MAX_ENLARGEMENTS)
            .expect("enlarging the time rescues the tail");
        assert!(order.iter().all(|&v| heat[v] > HEAT_UNDERFLOW));
    }

    #[test]
    fn solve_is_deterministic() {
        let mesh = icosphere(3);
        let h = mean_edge(&mesh);
        let a = solve(&mesh.0, &mesh.1, &[7], h * h, 2).expect("valid mesh");
        let b = solve(&mesh.0, &mesh.1, &[7], h * h, 2).expect("valid mesh");
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&a.potential), bits(&b.potential));
    }

    #[test]
    fn ten_thousand_vertices_solve_quickly() {
        let mesh = icosphere(5);
        assert_eq!(mesh.0.len(), 10242);
        let h = mean_edge(&mesh);
        let start = Instant::now();
        let d = heat_distances(&mesh, &[0], h * h, 1);
        let elapsed = start.elapsed().as_secs_f64();
        assert!(elapsed < 20.0, "10k-vertex solve took {elapsed:.2} s");
        let antipode = d[3];
        assert!(
            (antipode - std::f64::consts::PI).abs() < 0.03 * std::f64::consts::PI,
            "antipode distance {antipode}"
        );
    }

    #[test]
    fn mass_matrix_integrates_to_the_surface_area() {
        let (vertices, faces) = icosphere(3);
        let geometry: Vec<Option<FaceGeometry>> = faces
            .iter()
            .map(|t| FaceGeometry::new(&vertices, t))
            .collect();
        let mass = lumped_masses(vertices.len(), &faces, &geometry);
        let area: f64 = geometry.iter().flatten().map(|g| g.area).sum();
        let total: f64 = mass.iter().sum();
        assert!((total - area).abs() < 1e-12 * area, "{total} vs {area}");
        assert!(mass.iter().all(|&m| m > 0.0));
        // The unit sphere has area 4 pi; the polyhedron is slightly smaller.
        assert!((area / (4.0 * std::f64::consts::PI) - 1.0).abs() < 0.01);
    }

    #[test]
    fn stiffness_matrix_annihilates_constants_and_is_symmetric() {
        let (vertices, faces) = icosphere(2);
        let geometry: Vec<Option<FaceGeometry>> = faces
            .iter()
            .map(|t| FaceGeometry::new(&vertices, t))
            .collect();
        let graph = LaplaceGraph::build(vertices.len(), &faces, &geometry);
        for v in 0..vertices.len() {
            let row_sum: f64 = graph.diag[v]
                - (graph.ptr[v]..graph.ptr[v + 1])
                    .map(|k| graph.weight[k])
                    .sum::<f64>();
            assert!(row_sum.abs() < 1e-12, "row {v} sums to {row_sum}");
            assert!(graph.diag[v] > 0.0);
            for k in graph.ptr[v]..graph.ptr[v + 1] {
                let w = graph.idx[k];
                let back = (graph.ptr[w]..graph.ptr[w + 1])
                    .find(|&q| graph.idx[q] == v)
                    .expect("adjacency is symmetric");
                assert_eq!(graph.weight[k].to_bits(), graph.weight[back].to_bits());
            }
        }
    }
}
