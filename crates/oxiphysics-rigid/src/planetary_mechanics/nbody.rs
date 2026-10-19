// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Newtonian N-body gravity: pairwise summation, a Barnes-Hut octree and a
//! kick-drift-kick leapfrog integrator.
//!
//! Every interaction uses Plummer softening, i.e. a point mass `m` at
//! separation `d` pulls with `G m d / (|d|² + ε²)^(3/2)`, and the matching
//! pair potential is `-G m_i m_j / √(|d|² + ε²)`.  The force, the octree
//! multipole terms and the potential energy therefore all derive from the same
//! softened potential: the multipole terms are the Taylor expansion of the
//! softened kernel `(|r|² + ε²)^(-1/2)` about the centre of mass of a cell.

use super::{G_CONST, vec3_add, vec3_cross, vec3_dot, vec3_scale, vec3_sub};

/// Deepest octree level.
///
/// Bodies that still share a cell at this depth are merged into one
/// multi-body leaf, which bounds the recursion for (nearly) coincident bodies.
const MAX_TREE_DEPTH: usize = 64;

/// A particle in the N-body gravitational system.
#[derive(Debug, Clone)]
pub struct GravBody {
    /// Position \[m\].
    pub pos: [f64; 3],
    /// Velocity \[m/s\].
    pub vel: [f64; 3],
    /// Mass \[kg\].
    pub mass: f64,
}

impl GravBody {
    /// Create a new gravitational body.
    pub fn new(pos: [f64; 3], vel: [f64; 3], mass: f64) -> Self {
        Self { pos, vel, mass }
    }
}

/// Accelerations from the previous step, valid only for the exact state they
/// were computed from.
#[derive(Debug, Clone)]
struct AccelCache {
    positions: Vec<[f64; 3]>,
    masses: Vec<f64>,
    epsilon: f64,
    theta: f64,
    accel: Vec<[f64; 3]>,
}

impl AccelCache {
    /// Whether the cached accelerations belong to exactly this system state.
    fn matches(&self, bodies: &[GravBody], epsilon: f64, theta: f64) -> bool {
        let same_bits = |a: f64, b: f64| a.to_bits() == b.to_bits();
        self.accel.len() == bodies.len()
            && self.positions.len() == bodies.len()
            && self.masses.len() == bodies.len()
            && same_bits(self.epsilon, epsilon)
            && same_bits(self.theta, theta)
            && self.positions.iter().zip(bodies).all(|(p, b)| {
                same_bits(p[0], b.pos[0]) && same_bits(p[1], b.pos[1]) && same_bits(p[2], b.pos[2])
            })
            && self
                .masses
                .iter()
                .zip(bodies)
                .all(|(&m, b)| same_bits(m, b.mass))
    }
}

/// N-body gravitational simulation with exact pairwise or Barnes-Hut forces.
///
/// With `theta == 0` (the default of [`NBodyGravity::new`]) the forces are
/// summed over all pairs, applying Newton's third law so that linear momentum
/// is conserved to rounding.  With `theta > 0` the forces are evaluated on a
/// Barnes-Hut octree: a cell of side `s` whose centre of mass is at distance
/// `d` from the target body is replaced by its monopole plus quadrupole
/// moments when `s / d < theta` and the body lies outside the cell.  On a
/// seeded random cloud of 500 bodies, `theta = 0.5` reproduces the direct
/// accelerations with a median per-body relative error of about 0.1 % and no
/// body off by more than 0.5 % of the rms acceleration.
///
/// [`NBodyGravity::step`] advances the system with the second-order,
/// symplectic kick-drift-kick leapfrog.
#[derive(Debug, Clone)]
pub struct NBodyGravity {
    /// Bodies in the simulation.
    pub bodies: Vec<GravBody>,
    /// Plummer softening length to avoid the close-encounter singularity \[m\].
    pub epsilon: f64,
    /// Barnes-Hut opening angle (0 = direct summation).
    pub theta: f64,
    cache: Option<AccelCache>,
}

/// Plummer kernel `1 / (|dr|² + ε²)^(3/2)`; zero for a vanishing softened distance.
#[inline]
fn plummer_kernel(dr: [f64; 3], epsilon: f64) -> f64 {
    let r2 = vec3_dot(dr, dr) + epsilon * epsilon;
    if r2 > 0.0 {
        1.0 / (r2 * r2.sqrt())
    } else {
        0.0
    }
}

/// `m · d dᵀ` as `[xx, yy, zz, xy, xz, yz]`.
#[inline]
fn outer_moment(d: [f64; 3], mass: f64) -> [f64; 6] {
    [
        mass * d[0] * d[0],
        mass * d[1] * d[1],
        mass * d[2] * d[2],
        mass * d[0] * d[1],
        mass * d[0] * d[2],
        mass * d[1] * d[2],
    ]
}

/// Acceleration at `target` from a cell with the given mass, centre of mass
/// and second mass moment, using the softened monopole plus quadrupole terms.
///
/// With `R = target - com`, `ρ² = |R|² + ε²` and the kernel `f = 1/ρ`, the
/// potential of the cell is `-G (M f + ½ ∇∇f : Q)`; its negative gradient is
/// `-G M R / ρ³ + (G/2) [3 (tr Q · R + 2 Q R) / ρ⁵ - 15 (Rᵀ Q R) R / ρ⁷]`.
fn cell_acceleration(
    target: [f64; 3],
    com: [f64; 3],
    mass: f64,
    quad: &[f64; 6],
    epsilon: f64,
) -> [f64; 3] {
    let r = vec3_sub(target, com);
    let rho_sq = vec3_dot(r, r) + epsilon * epsilon;
    if rho_sq <= 0.0 {
        return [0.0; 3];
    }
    let rho = rho_sq.sqrt();
    let inv_rho3 = 1.0 / (rho_sq * rho);
    let monopole = vec3_scale(r, -G_CONST * mass * inv_rho3);

    let [qxx, qyy, qzz, qxy, qxz, qyz] = *quad;
    let trace = qxx + qyy + qzz;
    let q_r = [
        qxx * r[0] + qxy * r[1] + qxz * r[2],
        qxy * r[0] + qyy * r[1] + qyz * r[2],
        qxz * r[0] + qyz * r[1] + qzz * r[2],
    ];
    let r_q_r = vec3_dot(r, q_r);
    let inv_rho5 = inv_rho3 / rho_sq;
    let inv_rho7 = inv_rho5 / rho_sq;
    let along_r = 3.0 * trace * inv_rho5 - 15.0 * r_q_r * inv_rho7;
    let quadrupole = vec3_scale(
        vec3_add(vec3_scale(r, along_r), vec3_scale(q_r, 6.0 * inv_rho5)),
        0.5 * G_CONST,
    );
    vec3_add(monopole, quadrupole)
}

impl NBodyGravity {
    /// Create a new N-body system using direct summation.
    pub fn new(bodies: Vec<GravBody>) -> Self {
        Self {
            bodies,
            epsilon: 1e3,
            theta: 0.0,
            cache: None,
        }
    }

    /// Create with Barnes-Hut approximation enabled.
    ///
    /// `theta` is the opening angle; a value `<= 0` selects direct summation.
    pub fn with_barnes_hut(bodies: Vec<GravBody>, theta: f64) -> Self {
        Self {
            bodies,
            epsilon: 1e3,
            theta,
            cache: None,
        }
    }

    /// Exact softened accelerations of all bodies, summed over every pair.
    ///
    /// Each pair is visited once and applied to both bodies with opposite
    /// sign, so `Σ m_i a_i` vanishes to rounding.
    pub fn accelerations_direct(&self) -> Vec<[f64; 3]> {
        let mut acc = vec![[0.0_f64; 3]; self.bodies.len()];
        for (i, bi) in self.bodies.iter().enumerate() {
            for (j, bj) in self.bodies.iter().enumerate().skip(i + 1) {
                let dr = vec3_sub(bj.pos, bi.pos);
                let kernel = G_CONST * plummer_kernel(dr, self.epsilon);
                acc[i] = vec3_add(acc[i], vec3_scale(dr, kernel * bj.mass));
                acc[j] = vec3_sub(acc[j], vec3_scale(dr, kernel * bi.mass));
            }
        }
        acc
    }

    /// Accelerations of all bodies using the configured method.
    ///
    /// Barnes-Hut when `theta > 0`, otherwise (or when the positions are not
    /// finite, so no tree can be built) exact pairwise summation.
    pub fn accelerations(&self) -> Vec<[f64; 3]> {
        if self.theta > 0.0 {
            let tree = self.build_tree();
            if !tree.is_empty() {
                return (0..self.bodies.len())
                    .map(|i| self.accel_bh(&tree, i, self.theta))
                    .collect();
            }
        }
        self.accelerations_direct()
    }

    /// Build the Barnes-Hut octree for the current body positions.
    ///
    /// The result is an arena whose element 0 is the root cell (the smallest
    /// cube around all bodies); children always have larger indices than their
    /// parent.  Leaves carry one body (`body_idx`), or several merged bodies
    /// when they coincide.  Every node holds its total mass, centre of mass
    /// and second mass moment about the centre of mass.
    ///
    /// Returns an empty vector when there are no bodies or when a position is
    /// not finite.
    pub fn build_tree(&self) -> Vec<OctreeNode> {
        let Some(first) = self.bodies.first() else {
            return Vec::new();
        };
        let mut lo = first.pos;
        let mut hi = first.pos;
        for body in &self.bodies {
            for axis in 0..3 {
                let coord = body.pos[axis];
                if !coord.is_finite() {
                    return Vec::new();
                }
                lo[axis] = lo[axis].min(coord);
                hi[axis] = hi[axis].max(coord);
            }
        }
        let center = [
            0.5 * (lo[0] + hi[0]),
            0.5 * (lo[1] + hi[1]),
            0.5 * (lo[2] + hi[2]),
        ];
        let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(hi[2] - lo[2]);
        if !center.iter().all(|c| c.is_finite()) || !extent.is_finite() {
            return Vec::new();
        }
        let half_width = (0.5 * extent * (1.0 + 1e-9)).max(f64::MIN_POSITIVE);

        let mut nodes = vec![OctreeNode::new(center, half_width)];
        for index in 0..self.bodies.len() {
            self.insert_body(&mut nodes, index);
        }
        // Children are created after their parent, so a reverse sweep visits
        // every child before the node that aggregates it.
        for index in (0..nodes.len()).rev() {
            let children = nodes[index].children;
            if children.iter().all(Option::is_none) {
                continue;
            }
            let mut mass = 0.0;
            let mut moment = [0.0_f64; 3];
            for &child in children.iter().flatten() {
                let (child_mass, child_com) = (nodes[child].mass, nodes[child].com);
                mass += child_mass;
                moment = vec3_add(moment, vec3_scale(child_com, child_mass));
            }
            let com = if mass != 0.0 {
                vec3_scale(moment, 1.0 / mass)
            } else {
                nodes[index].center
            };
            // Parallel-axis theorem: shift every child's moment to the new centre of mass.
            let mut quad = [0.0_f64; 6];
            for &child in children.iter().flatten() {
                let offset = vec3_sub(nodes[child].com, com);
                let shifted = outer_moment(offset, nodes[child].mass);
                for (total, (own, shift)) in
                    quad.iter_mut().zip(nodes[child].quad.iter().zip(shifted))
                {
                    *total += own + shift;
                }
            }
            let node = &mut nodes[index];
            node.mass = mass;
            node.com = com;
            node.quad = quad;
        }
        nodes
    }

    /// Insert body `body` into the arena, subdividing occupied leaves.
    fn insert_body(&self, nodes: &mut Vec<OctreeNode>, body: usize) {
        let pos = self.bodies[body].pos;
        let mass = self.bodies[body].mass;
        let mut current = 0_usize;
        let mut depth = 0_usize;
        loop {
            if nodes[current].is_leaf() {
                let Some(resident) = nodes[current].body_idx else {
                    let node = &mut nodes[current];
                    node.body_idx = Some(body);
                    node.mass = mass;
                    node.com = pos;
                    return;
                };
                if depth >= MAX_TREE_DEPTH || self.bodies[resident].pos == pos {
                    // Merge into the resident leaf (coincident or too deep to split).
                    let node = &mut nodes[current];
                    let total = node.mass + mass;
                    if total != 0.0 {
                        let moment =
                            vec3_add(vec3_scale(node.com, node.mass), vec3_scale(pos, mass));
                        node.com = vec3_scale(moment, 1.0 / total);
                    }
                    node.mass = total;
                    return;
                }
                // Push the resident content one level down, then keep inserting.
                let (resident_mass, resident_com) = (nodes[current].mass, nodes[current].com);
                let octant = nodes[current].octant_for(self.bodies[resident].pos);
                let child = Self::add_child(nodes, current, octant);
                nodes[child].body_idx = Some(resident);
                nodes[child].mass = resident_mass;
                nodes[child].com = resident_com;
                nodes[current].body_idx = None;
            }
            let octant = nodes[current].octant_for(pos);
            current = match nodes[current].children[octant] {
                Some(existing) => existing,
                None => Self::add_child(nodes, current, octant),
            };
            depth += 1;
        }
    }

    /// Append the child cell of `parent` in `octant` and link it.
    fn add_child(nodes: &mut Vec<OctreeNode>, parent: usize, octant: usize) -> usize {
        let child = OctreeNode::new(
            nodes[parent].child_center(octant),
            0.5 * nodes[parent].half_width,
        );
        nodes.push(child);
        let index = nodes.len() - 1;
        nodes[parent].children[octant] = Some(index);
        index
    }

    /// Barnes-Hut acceleration of body `i` using a tree from [`Self::build_tree`].
    ///
    /// A cell is accepted, represented by its monopole and quadrupole moments
    /// about its centre of mass, when its side `2·half_width` is smaller than
    /// `theta` times the distance to that centre of mass and the body is
    /// outside the cell; otherwise the cell is opened.  `theta == 0`
    /// therefore reproduces direct summation.
    /// Returns zero for an invalid index or an empty tree.
    pub fn accel_bh(&self, tree: &[OctreeNode], i: usize, theta: f64) -> [f64; 3] {
        let mut acc = [0.0_f64; 3];
        if let Some(body) = self.bodies.get(i) {
            self.accumulate_cell(tree, 0, i, body.pos, theta * theta, &mut acc);
        }
        acc
    }

    /// Add the pull of the subtree rooted at `cell` on the body `target`.
    fn accumulate_cell(
        &self,
        tree: &[OctreeNode],
        cell: usize,
        target: usize,
        target_pos: [f64; 3],
        theta_sq: f64,
        acc: &mut [f64; 3],
    ) {
        let Some(node) = tree.get(cell) else {
            return;
        };
        if node.mass == 0.0 {
            return;
        }
        let leaf = node.is_leaf();
        let accept = if leaf {
            if node.body_idx == Some(target) {
                return;
            }
            true
        } else {
            let dr = vec3_sub(node.com, target_pos);
            let side = 2.0 * node.half_width;
            !node.contains(target_pos) && side * side < theta_sq * vec3_dot(dr, dr)
        };
        if accept {
            let quad = if leaf { &[0.0; 6] } else { &node.quad };
            let pull = cell_acceleration(target_pos, node.com, node.mass, quad, self.epsilon);
            *acc = vec3_add(*acc, pull);
            return;
        }
        for &child in node.children.iter().flatten() {
            self.accumulate_cell(tree, child, target, target_pos, theta_sq, acc);
        }
    }

    /// Step the simulation forward by `dt` \[s\] with kick-drift-kick leapfrog.
    ///
    /// `v += a dt/2; x += v dt; v += a' dt/2`, where `a'` is evaluated at the
    /// new positions.  The scheme is second-order accurate, time-reversible
    /// and symplectic, so with exact pairwise forces the energy error stays
    /// bounded instead of drifting (a Barnes-Hut force is only approximately
    /// conservative, so there the bound is set by the tree error).  The
    /// closing acceleration is kept and reused as the opening acceleration of
    /// the next step, so each step costs one force evaluation; the cache is
    /// discarded whenever the public fields no longer match the state it was
    /// computed for.
    pub fn step(&mut self, dt: f64) {
        if self.bodies.is_empty() {
            return;
        }
        let epsilon = self.epsilon;
        let theta = self.theta;
        let opening = match self.cache.take() {
            Some(cache) if cache.matches(&self.bodies, epsilon, theta) => cache.accel,
            _ => self.accelerations(),
        };
        for (body, acc) in self.bodies.iter_mut().zip(&opening) {
            body.vel = vec3_add(body.vel, vec3_scale(*acc, 0.5 * dt));
            body.pos = vec3_add(body.pos, vec3_scale(body.vel, dt));
        }
        let closing = self.accelerations();
        for (body, acc) in self.bodies.iter_mut().zip(&closing) {
            body.vel = vec3_add(body.vel, vec3_scale(*acc, 0.5 * dt));
        }
        self.cache = Some(AccelCache {
            positions: self.bodies.iter().map(|b| b.pos).collect(),
            masses: self.bodies.iter().map(|b| b.mass).collect(),
            epsilon,
            theta,
            accel: closing,
        });
    }

    /// Compute total kinetic energy of the system \[J\].
    pub fn kinetic_energy(&self) -> f64 {
        self.bodies
            .iter()
            .map(|b| {
                let v2 = vec3_dot(b.vel, b.vel);
                0.5 * b.mass * v2
            })
            .sum()
    }

    /// Compute total gravitational potential energy of the system \[J\].
    ///
    /// Sums the Plummer-softened pair potential `-G m_i m_j / √(r² + ε²)`
    /// over all pairs, the potential whose gradient is the force applied by
    /// [`Self::step`].
    pub fn potential_energy(&self) -> f64 {
        let mut e = 0.0;
        for (i, bi) in self.bodies.iter().enumerate() {
            for bj in self.bodies.iter().skip(i + 1) {
                let dr = vec3_sub(bj.pos, bi.pos);
                let softened = (vec3_dot(dr, dr) + self.epsilon * self.epsilon).sqrt();
                if softened > 0.0 {
                    e -= G_CONST * bi.mass * bj.mass / softened;
                }
            }
        }
        e
    }

    /// Compute total angular momentum of the system \[kg·m²/s\].
    pub fn angular_momentum(&self) -> [f64; 3] {
        let mut l = [0.0_f64; 3];
        for b in &self.bodies {
            let li = vec3_cross(b.pos, vec3_scale(b.vel, b.mass));
            l = vec3_add(l, li);
        }
        l
    }

    /// Compute center of mass position \[m\] (the origin for a massless system).
    pub fn center_of_mass(&self) -> [f64; 3] {
        let mut total_mass = 0.0;
        let mut com = [0.0_f64; 3];
        for b in &self.bodies {
            com = vec3_add(com, vec3_scale(b.pos, b.mass));
            total_mass += b.mass;
        }
        if total_mass == 0.0 {
            return [0.0; 3];
        }
        vec3_scale(com, 1.0 / total_mass)
    }
}

// ─── OctreeNode for Barnes-Hut ─────────────────────────────────────────────────

/// A node in the Barnes-Hut octree.
#[derive(Debug, Clone)]
pub struct OctreeNode {
    /// Centre of the cell \[m\].
    pub center: [f64; 3],
    /// Half-width of the cell \[m\].
    pub half_width: f64,
    /// Total mass in this cell \[kg\].
    pub mass: f64,
    /// Centre of mass of this cell \[m\].
    pub com: [f64; 3],
    /// Second mass moment about the centre of mass,
    /// `Σ m (r - com)(r - com)ᵀ` \[kg·m²\], stored as `[xx, yy, zz, xy, xz, yz]`.
    ///
    /// Zero for leaves, whose bodies sit at the centre of mass.
    pub quad: [f64; 6],
    /// Child indices (8 children) in a flat arena.
    pub children: [Option<usize>; 8],
    /// Body index if this is a leaf.
    pub body_idx: Option<usize>,
}

impl OctreeNode {
    /// Create a new empty octree node.
    pub fn new(center: [f64; 3], half_width: f64) -> Self {
        Self {
            center,
            half_width,
            mass: 0.0,
            com: [0.0; 3],
            quad: [0.0; 6],
            children: [None; 8],
            body_idx: None,
        }
    }

    /// Whether the node has no children.
    pub fn is_leaf(&self) -> bool {
        self.children.iter().all(Option::is_none)
    }

    /// Whether `pos` lies inside the cell or on its boundary.
    pub fn contains(&self, pos: [f64; 3]) -> bool {
        (0..3).all(|axis| (pos[axis] - self.center[axis]).abs() <= self.half_width)
    }

    /// Compute which of the 8 octants a position belongs to.
    pub fn octant_for(&self, pos: [f64; 3]) -> usize {
        let mut idx = 0;
        if pos[0] >= self.center[0] {
            idx |= 1;
        }
        if pos[1] >= self.center[1] {
            idx |= 2;
        }
        if pos[2] >= self.center[2] {
            idx |= 4;
        }
        idx
    }

    /// Child centre for a given octant index.
    pub fn child_center(&self, oct: usize) -> [f64; 3] {
        let h = self.half_width / 2.0;
        [
            self.center[0] + if oct & 1 != 0 { h } else { -h },
            self.center[1] + if oct & 2 != 0 { h } else { -h },
            self.center[2] + if oct & 4 != 0 { h } else { -h },
        ]
    }
}
