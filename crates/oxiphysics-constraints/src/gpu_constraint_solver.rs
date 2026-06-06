// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated constraint solver using the wgpu compute backend.
//!
//! This module provides a GPU-dispatch path for the batch Projected Gauss-Seidel
//! (PGS) constraint solver.  When a GPU device is available the constraint
//! impulse computation is uploaded as a WGSL compute shader; when no GPU is
//! present the solver falls back transparently to a CPU PGS loop.
//!
//! # Architecture
//!
//! ```text
//!  ┌──────────────────────────────────────────┐
//!  │  GpuConstraintSolver                     │
//!  │  ┌────────────────┐  ┌─────────────────┐ │
//!  │  │ GPU path        │  │ CPU fallback    │ │
//!  │  │ WgpuBackend +  │  │ CPU PGS loop    │ │
//!  │  │ WGSL kernels   │  │                 │ │
//!  │  └────────────────┘  └─────────────────┘ │
//!  └──────────────────────────────────────────┘
//! ```
//!
//! # GPU kernel overview
//!
//! The WGSL kernel `constraint_pgs_iter` runs one full PGS iteration over N
//! constraints.  Each workgroup handles 64 constraints sequentially (true
//! Gauss-Seidel within the block); across workgroups the update is Jacobi-style,
//! making this *block-PGS* with block size 64.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_constraints::gpu_constraint_solver::{
//!     GpuConstraintSolver, GpuSolverConfig, CpuConstraintData,
//! };
//!
//! let config = GpuSolverConfig { iterations: 8, omega: 1.0 };
//! let mut solver = GpuConstraintSolver::new(config);
//!
//! let data = CpuConstraintData::empty();
//! let results = solver.solve(&data);
//! println!("Solved {} constraints (GPU={})", results.n_constraints, results.used_gpu);
//! ```

use oxiphysics_gpu::compute::{WgpuBackend, WgpuBufferHandle};

// ── WGSL source ───────────────────────────────────────────────────────────────

/// WGSL source for one block-PGS iteration over N constraints.
///
/// Each workgroup of 64 threads handles a contiguous block of constraints.
/// Within the block constraints are solved sequentially (Gauss-Seidel order);
/// across blocks the velocity update is Jacobi (block-PGS).
///
/// Multiple dispatches of this shader (controlled by `GpuSolverConfig::iterations`)
/// are required to achieve convergence.
pub const WGSL_CONSTRAINT_PGS: &str = r#"
struct GpuConstraint {
    nx: f32, ny: f32, nz: f32, em: f32,
    bias: f32,
    lambda_lo: f32, lambda_hi: f32,
    body_a: u32, body_b: u32,
    rax: f32, ray: f32, raz: f32,
    rbx: f32, rby: f32, rbz: f32,
    _pad0: f32, _pad1: f32, _pad2: f32, _pad3: f32, _pad4: f32,
}

struct Uniforms {
    n_constraints: u32,
    n_bodies: u32,
    omega: f32,
    _pad: f32,
}

@group(0) @binding(0) var<storage, read>       constraints : array<GpuConstraint>;
@group(0) @binding(1) var<storage, read_write> lambda      : array<f32>;
@group(0) @binding(2) var<storage, read_write> vel_lin     : array<vec4<f32>>;
@group(0) @binding(3) var<storage, read_write> vel_ang     : array<vec4<f32>>;
@group(0) @binding(4) var<uniform>             uniforms    : Uniforms;

const WG: u32 = 64u;

@compute @workgroup_size(WG, 1, 1)
fn constraint_pgs_iter(
    @builtin(workgroup_id)           wg_id : vec3<u32>,
    @builtin(local_invocation_index) lid   : u32,
) {
    let base = wg_id.x * WG;
    if lid != 0u { return; }

    for (var ci: u32 = base; ci < min(base + WG, uniforms.n_constraints); ci = ci + 1u) {
        let c = constraints[ci];

        var vla = vec3<f32>(0.0); var wla = vec3<f32>(0.0); var inv_ma = 0.0f;
        if c.body_a != 0xFFFFFFFFu {
            vla = vel_lin[c.body_a].xyz; wla = vel_ang[c.body_a].xyz;
            inv_ma = vel_lin[c.body_a].w;
        }
        var vlb = vec3<f32>(0.0); var wlb = vec3<f32>(0.0); var inv_mb = 0.0f;
        if c.body_b != 0xFFFFFFFFu {
            vlb = vel_lin[c.body_b].xyz; wlb = vel_ang[c.body_b].xyz;
            inv_mb = vel_lin[c.body_b].w;
        }

        let n  = vec3<f32>(c.nx, c.ny, c.nz);
        let ra = vec3<f32>(c.rax, c.ray, c.raz);
        let rb = vec3<f32>(c.rbx, c.rby, c.rbz);
        let va = vla + cross(wla, ra);
        let vb = vlb + cross(wlb, rb);
        let rv = dot(n, va - vb);

        let d_lam_raw = -(rv + c.bias) * c.em * uniforms.omega;
        let old_lam   = lambda[ci];
        let new_lam   = clamp(old_lam + d_lam_raw, c.lambda_lo, c.lambda_hi);
        let d_lam     = new_lam - old_lam;
        lambda[ci]    = new_lam;

        let imp = n * d_lam;
        if c.body_a != 0xFFFFFFFFu {
            vel_lin[c.body_a] = vec4<f32>(vla + imp * inv_ma, inv_ma);
            vel_ang[c.body_a] = vec4<f32>(wla + cross(ra, imp) * inv_ma, 0.0);
        }
        if c.body_b != 0xFFFFFFFFu {
            vel_lin[c.body_b] = vec4<f32>(vlb - imp * inv_mb, inv_mb);
            vel_ang[c.body_b] = vec4<f32>(wlb - cross(rb, imp) * inv_mb, 0.0);
        }
    }
}
"#;

// ── CPU-side data structures ──────────────────────────────────────────────────

/// Single constraint in a form suitable for GPU upload (matches WGSL struct layout).
#[derive(Debug, Clone, Copy, Default)]
#[repr(C)]
pub struct GpuConstraint {
    /// Contact normal x.
    pub nx: f32,
    /// Contact normal y.
    pub ny: f32,
    /// Contact normal z.
    pub nz: f32,
    /// Effective mass (inverse of sum of mass-weighted Jacobian terms).
    pub em: f32,
    /// Bias velocity (Baumgarte stabilisation + restitution).
    pub bias: f32,
    /// Lower lambda clamp bound (0 for contact, `f32::NEG_INFINITY` for bilateral).
    pub lambda_lo: f32,
    /// Upper lambda clamp bound (`f32::MAX`).
    pub lambda_hi: f32,
    /// Body A index (`u32::MAX` for static world).
    pub body_a: u32,
    /// Body B index (`u32::MAX` for static world).
    pub body_b: u32,
    /// r-vector from body A CoM to contact point, x component.
    pub rax: f32,
    /// r-vector from body A CoM to contact point, y component.
    pub ray: f32,
    /// r-vector from body A CoM to contact point, z component.
    pub raz: f32,
    /// r-vector from body B CoM to contact point, x component.
    pub rbx: f32,
    /// r-vector from body B CoM to contact point, y component.
    pub rby: f32,
    /// r-vector from body B CoM to contact point, z component.
    pub rbz: f32,
    /// Padding byte 61–64 (WGSL 16-byte alignment requirement).
    pub _pad0: f32,
    /// Padding byte 65–68.
    pub _pad1: f32,
    /// Padding byte 69–72.
    pub _pad2: f32,
    /// Padding byte 73–76.
    pub _pad3: f32,
    /// Padding byte 77–80 — brings total to 80 bytes (80 % 16 == 0).
    pub _pad4: f32,
}

// GpuConstraint layout (repr(C)):
//   nx,ny,nz,em          =  4 × f32 = 16 bytes  (cumul  16)
//   bias,lambda_lo,lambda_hi = 3 × f32 = 12 bytes  (cumul  28)
//   body_a,body_b        =  2 × u32  =  8 bytes  (cumul  36)
//   rax,ray,raz          =  3 × f32 = 12 bytes  (cumul  48)
//   rbx,rby,rbz          =  3 × f32 = 12 bytes  (cumul  60)
//   _pad0.._pad4         =  5 × f32 = 20 bytes  (cumul  80)  ← 80 % 16 == 0 ✓
const CONSTRAINT_F64_SLOTS: usize = 20;

/// Parameters shared by [`GpuConstraint::contact`] and [`GpuConstraint::bilateral`].
#[derive(Debug, Clone, Copy)]
pub struct GpuConstraintParams {
    /// Contact normal x component.
    pub nx: f32,
    /// Contact normal y component.
    pub ny: f32,
    /// Contact normal z component.
    pub nz: f32,
    /// Effective mass (1/K).
    pub effective_mass: f32,
    /// Velocity bias (Baumgarte + restitution).
    pub bias: f32,
    /// Index of body A.
    pub body_a: u32,
    /// Index of body B.
    pub body_b: u32,
    /// r-vector from body A CoM to contact, x.
    pub rax: f32,
    /// r-vector from body A CoM to contact, y.
    pub ray: f32,
    /// r-vector from body A CoM to contact, z.
    pub raz: f32,
    /// r-vector from body B CoM to contact, x.
    pub rbx: f32,
    /// r-vector from body B CoM to contact, y.
    pub rby: f32,
    /// r-vector from body B CoM to contact, z.
    pub rbz: f32,
}

impl GpuConstraint {
    /// Create a contact constraint (lambda ≥ 0).
    pub fn contact(p: GpuConstraintParams) -> Self {
        Self {
            nx: p.nx,
            ny: p.ny,
            nz: p.nz,
            em: p.effective_mass,
            bias: p.bias,
            lambda_lo: 0.0,
            lambda_hi: f32::MAX,
            body_a: p.body_a,
            body_b: p.body_b,
            rax: p.rax,
            ray: p.ray,
            raz: p.raz,
            rbx: p.rbx,
            rby: p.rby,
            rbz: p.rbz,
            _pad0: 0.0,
            _pad1: 0.0,
            _pad2: 0.0,
            _pad3: 0.0,
            _pad4: 0.0,
        }
    }

    /// Create a bilateral (equality) constraint (lambda unconstrained).
    pub fn bilateral(p: GpuConstraintParams) -> Self {
        Self {
            lambda_lo: f32::NEG_INFINITY,
            lambda_hi: f32::MAX,
            ..Self::contact(p)
        }
    }
}

/// Body velocity entry: `[vx, vy, vz, inv_mass]`.
pub type GpuBodyVel = [f32; 4];

/// Full CPU-side constraint data ready for GPU upload or CPU fallback.
#[derive(Debug, Default)]
pub struct CpuConstraintData {
    /// Constraint descriptors.
    pub constraints: Vec<GpuConstraint>,
    /// Body linear velocities + inverse mass (indexed by body handle).
    pub vel_lin: Vec<GpuBodyVel>,
    /// Body angular velocities (indexed by body handle, w component unused).
    pub vel_ang: Vec<GpuBodyVel>,
    /// Accumulated impulses (one per constraint; zero-initialise for cold start).
    pub lambda: Vec<f32>,
}

impl CpuConstraintData {
    /// Create an empty data set.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Allocate storage for `n_constraints` constraints and `n_bodies` bodies.
    pub fn new(n_constraints: usize, n_bodies: usize) -> Self {
        Self {
            constraints: vec![GpuConstraint::default(); n_constraints],
            vel_lin: vec![[0.0; 4]; n_bodies],
            vel_ang: vec![[0.0; 4]; n_bodies],
            lambda: vec![0.0; n_constraints],
        }
    }

    /// Number of constraints.
    pub fn n_constraints(&self) -> usize {
        self.constraints.len()
    }

    /// Number of bodies.
    pub fn n_bodies(&self) -> usize {
        self.vel_lin.len()
    }
}

// ── GpuSolverConfig ───────────────────────────────────────────────────────────

/// Configuration for the GPU constraint solver.
#[derive(Debug, Clone)]
pub struct GpuSolverConfig {
    /// Number of PGS iterations per physics step (default 8).
    pub iterations: u32,
    /// SOR over-relaxation factor (1.0 = standard PGS, typically 1.0–1.5).
    pub omega: f32,
}

impl Default for GpuSolverConfig {
    fn default() -> Self {
        Self {
            iterations: 8,
            omega: 1.0,
        }
    }
}

// ── SolveResult ───────────────────────────────────────────────────────────────

/// Output from one [`GpuConstraintSolver::solve`] call.
#[derive(Debug)]
pub struct SolveResult {
    /// Number of constraints processed.
    pub n_constraints: usize,
    /// Whether the GPU path was actually used (`false` → CPU fallback).
    pub used_gpu: bool,
    /// Updated body linear velocities + inverse mass.
    pub vel_lin: Vec<GpuBodyVel>,
    /// Updated body angular velocities.
    pub vel_ang: Vec<GpuBodyVel>,
    /// Final accumulated impulses (warm-start these next frame).
    pub lambda: Vec<f32>,
}

// ── GpuConstraintSolver ───────────────────────────────────────────────────────

/// GPU-accelerated projected Gauss-Seidel constraint solver.
///
/// Attempts to use [`WgpuBackend`]; falls back to a CPU PGS loop when no GPU
/// device is available or the backend returns an error.
pub struct GpuConstraintSolver {
    /// Solver configuration.
    pub config: GpuSolverConfig,
    /// GPU backend (`None` → CPU fallback).
    backend: Option<WgpuBackend>,
    // GPU buffer handles (re-created when data size changes)
    buf_constraints: Option<WgpuBufferHandle>,
    buf_lambda: Option<WgpuBufferHandle>,
    buf_vel_lin: Option<WgpuBufferHandle>,
    buf_vel_ang: Option<WgpuBufferHandle>,
    /// Cached sizes to detect resize.
    last_nc: usize,
    last_nb: usize,
}

impl GpuConstraintSolver {
    /// Create a new solver.  Attempts GPU initialisation; falls back silently.
    pub fn new(config: GpuSolverConfig) -> Self {
        let backend = match WgpuBackend::try_new() {
            Ok(mut b) => {
                b.register_shader("constraint_pgs", WGSL_CONSTRAINT_PGS);
                Some(b)
            }
            Err(_) => None,
        };

        Self {
            config,
            backend,
            buf_constraints: None,
            buf_lambda: None,
            buf_vel_lin: None,
            buf_vel_ang: None,
            last_nc: 0,
            last_nb: 0,
        }
    }

    /// Create a solver that explicitly uses the CPU path (useful for tests).
    pub fn cpu_only(config: GpuSolverConfig) -> Self {
        Self {
            config,
            backend: None,
            buf_constraints: None,
            buf_lambda: None,
            buf_vel_lin: None,
            buf_vel_ang: None,
            last_nc: 0,
            last_nb: 0,
        }
    }

    /// `true` if a GPU backend is active.
    pub fn has_gpu(&self) -> bool {
        self.backend.is_some()
    }

    /// Solve the constraint system, choosing GPU or CPU path automatically.
    pub fn solve(&mut self, data: &CpuConstraintData) -> SolveResult {
        if self.backend.is_some() {
            self.solve_gpu(data)
        } else {
            self.solve_cpu(data)
        }
    }

    // ── GPU path ─────────────────────────────────────────────────────────────

    fn solve_gpu(&mut self, data: &CpuConstraintData) -> SolveResult {
        let nc = data.n_constraints();
        let nb = data.n_bodies();

        let backend = self
            .backend
            .as_mut()
            .expect("solve_gpu called only when backend is Some");

        // (Re-)allocate GPU buffers when data size changes.
        if nc != self.last_nc || nb != self.last_nb {
            // The WgpuBackend stub allocates Vec<f64> of length `len`.
            // We flatten each GpuConstraint into CONSTRAINT_F64_SLOTS f64 slots.
            self.buf_constraints = Some(backend.create_buffer(nc * CONSTRAINT_F64_SLOTS));
            self.buf_lambda = Some(backend.create_buffer(nc));
            self.buf_vel_lin = Some(backend.create_buffer(nb * 4));
            self.buf_vel_ang = Some(backend.create_buffer(nb * 4));
            self.last_nc = nc;
            self.last_nb = nb;
        }

        let bc = self
            .buf_constraints
            .expect("buf_constraints allocated above when size changed");
        let bl = self
            .buf_lambda
            .expect("buf_lambda allocated above when size changed");
        let bvl = self
            .buf_vel_lin
            .expect("buf_vel_lin allocated above when size changed");
        let bva = self
            .buf_vel_ang
            .expect("buf_vel_ang allocated above when size changed");

        // Upload (convert f32 → f64 for backend)
        backend.write_buffer(bc, &constraints_to_f64(&data.constraints));
        backend.write_buffer(
            bl,
            &data.lambda.iter().map(|&v| v as f64).collect::<Vec<_>>(),
        );
        backend.write_buffer(bvl, &body_vels_to_f64(&data.vel_lin));
        backend.write_buffer(bva, &body_vels_to_f64(&data.vel_ang));

        // Dispatch iterations (each dispatch = one full sweep of all constraints)
        let wg = (nc as u32).div_ceil(64);
        for _ in 0..self.config.iterations {
            backend.dispatch("constraint_pgs", &[bc, bl, bvl, bva], wg);
        }

        // Download and convert f64 → f32
        let lambda_raw = backend.read_buffer(bl);
        let vel_lin_raw = backend.read_buffer(bvl);
        let vel_ang_raw = backend.read_buffer(bva);

        let lambda = lambda_raw.iter().map(|&v| v as f32).collect();
        let vel_lin = f64_to_body_vels(&vel_lin_raw, nb);
        let vel_ang = f64_to_body_vels(&vel_ang_raw, nb);

        SolveResult {
            n_constraints: nc,
            used_gpu: true,
            vel_lin,
            vel_ang,
            lambda,
        }
    }

    // ── CPU fallback ──────────────────────────────────────────────────────────

    fn solve_cpu(&self, data: &CpuConstraintData) -> SolveResult {
        let nc = data.n_constraints();

        let mut lambda = data.lambda.clone();
        let mut vel_lin = data.vel_lin.clone();
        let mut vel_ang = data.vel_ang.clone();

        let omega = self.config.omega;

        for _iter in 0..self.config.iterations {
            for (ci, c) in data.constraints.iter().enumerate() {
                let (vla, wla, inv_ma) = gather_body(&vel_lin, &vel_ang, c.body_a);
                let (vlb, wlb, inv_mb) = gather_body(&vel_lin, &vel_ang, c.body_b);

                // Velocity at contact point: v + ω × r
                let va = add(vla, cross(wla, [c.rax, c.ray, c.raz]));
                let vb = add(vlb, cross(wlb, [c.rbx, c.rby, c.rbz]));

                let rv = dot([c.nx, c.ny, c.nz], sub(va, vb));

                let d_lam_raw = -(rv + c.bias) * c.em * omega;
                let old_lam = lambda[ci];
                let new_lam = (old_lam + d_lam_raw).clamp(c.lambda_lo, c.lambda_hi);
                let d_lam = new_lam - old_lam;
                lambda[ci] = new_lam;

                let imp = [c.nx * d_lam, c.ny * d_lam, c.nz * d_lam];

                apply_impulse(
                    &mut vel_lin,
                    &mut vel_ang,
                    c.body_a,
                    imp,
                    [c.rax, c.ray, c.raz],
                    inv_ma,
                    1.0,
                );
                apply_impulse(
                    &mut vel_lin,
                    &mut vel_ang,
                    c.body_b,
                    imp,
                    [c.rbx, c.rby, c.rbz],
                    inv_mb,
                    -1.0,
                );
            }
        }

        SolveResult {
            n_constraints: nc,
            used_gpu: false,
            vel_lin,
            vel_ang,
            lambda,
        }
    }
}

// ── small vector helpers ──────────────────────────────────────────────────────

#[inline]
fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn add(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Gather linear+angular velocity for body `idx` (u32::MAX → zeroes).
#[inline]
fn gather_body(vel_lin: &[[f32; 4]], vel_ang: &[[f32; 4]], idx: u32) -> ([f32; 3], [f32; 3], f32) {
    if idx == u32::MAX {
        ([0.0; 3], [0.0; 3], 0.0)
    } else {
        let i = idx as usize;
        (
            [vel_lin[i][0], vel_lin[i][1], vel_lin[i][2]],
            [vel_ang[i][0], vel_ang[i][1], vel_ang[i][2]],
            vel_lin[i][3],
        )
    }
}

/// Apply impulse to body `idx` (skip if u32::MAX).
/// `sign` = +1 for body A, -1 for body B.
#[inline]
fn apply_impulse(
    vel_lin: &mut [[f32; 4]],
    vel_ang: &mut [[f32; 4]],
    idx: u32,
    imp: [f32; 3],
    r: [f32; 3],
    inv_m: f32,
    sign: f32,
) {
    if idx == u32::MAX {
        return;
    }
    let i = idx as usize;
    vel_lin[i][0] += sign * imp[0] * inv_m;
    vel_lin[i][1] += sign * imp[1] * inv_m;
    vel_lin[i][2] += sign * imp[2] * inv_m;
    let torque = cross(r, imp);
    vel_ang[i][0] += sign * torque[0] * inv_m;
    vel_ang[i][1] += sign * torque[1] * inv_m;
    vel_ang[i][2] += sign * torque[2] * inv_m;
}

// ── backend upload/download helpers ──────────────────────────────────────────

/// Flatten `Vec<GpuConstraint>` to `Vec<f64>` for backend upload.
///
/// Each constraint occupies `CONSTRAINT_F64_SLOTS` slots.  The two `u32` body
/// indices are bitcast to f64 by value (safe for indices < 2^32).
fn constraints_to_f64(cs: &[GpuConstraint]) -> Vec<f64> {
    let mut out = Vec::with_capacity(cs.len() * CONSTRAINT_F64_SLOTS);
    for c in cs {
        out.push(c.nx as f64);
        out.push(c.ny as f64);
        out.push(c.nz as f64);
        out.push(c.em as f64);
        out.push(c.bias as f64);
        out.push(c.lambda_lo as f64);
        out.push(c.lambda_hi as f64);
        out.push(c.body_a as f64);
        out.push(c.body_b as f64);
        out.push(c.rax as f64);
        out.push(c.ray as f64);
        out.push(c.raz as f64);
        out.push(c.rbx as f64);
        out.push(c.rby as f64);
        out.push(c.rbz as f64);
        out.push(c._pad0 as f64);
        out.push(c._pad1 as f64);
        out.push(c._pad2 as f64);
        out.push(c._pad3 as f64);
        out.push(c._pad4 as f64);
    }
    out
}

/// Flatten `Vec<GpuBodyVel>` (`[f32; 4]`) to `Vec<f64>` for backend upload.
fn body_vels_to_f64(vels: &[[f32; 4]]) -> Vec<f64> {
    vels.iter()
        .flat_map(|v| v.iter().map(|&x| x as f64))
        .collect()
}

/// Reconstruct `Vec<[f32; 4]>` from flat `Vec<f64>` downloaded from backend.
fn f64_to_body_vels(flat: &[f64], n: usize) -> Vec<[f32; 4]> {
    (0..n)
        .map(|i| {
            let base = i * 4;
            [
                flat[base] as f32,
                flat[base + 1] as f32,
                flat[base + 2] as f32,
                flat[base + 3] as f32,
            ]
        })
        .collect()
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cpu_fallback_empty() {
        let mut solver = GpuConstraintSolver::cpu_only(GpuSolverConfig::default());
        let data = CpuConstraintData::empty();
        let r = solver.solve(&data);
        assert_eq!(r.n_constraints, 0);
        assert!(!r.used_gpu);
    }

    #[test]
    fn test_cpu_single_contact() {
        // Body 0 approaching body 1 (static) along +Y.
        // After PGS the Y velocity of body 0 should be non-negative.
        let cfg = GpuSolverConfig {
            iterations: 20,
            omega: 1.0,
        };
        let mut solver = GpuConstraintSolver::cpu_only(cfg);

        let mut data = CpuConstraintData::new(1, 2);
        // Body 0: moving in -Y at 1 m/s, mass 1 kg → inv_mass = 1
        data.vel_lin[0] = [0.0, -1.0, 0.0, 1.0];
        // Body 1: static (inv_mass = 0)
        data.vel_lin[1] = [0.0, 0.0, 0.0, 0.0];

        data.constraints[0] = GpuConstraint {
            nx: 0.0,
            ny: 1.0,
            nz: 0.0,
            em: 1.0, // 1 / (inv_ma + inv_mb + ...) = 1 / (1 + 0)
            bias: 0.0,
            lambda_lo: 0.0,
            lambda_hi: f32::MAX,
            body_a: 0,
            body_b: u32::MAX,
            rax: 0.0,
            ray: 0.0,
            raz: 0.0,
            rbx: 0.0,
            rby: 0.0,
            rbz: 0.0,
            _pad0: 0.0,
            _pad1: 0.0,
            _pad2: 0.0,
            _pad3: 0.0,
            _pad4: 0.0,
        };

        let r = solver.solve(&data);
        assert!(!r.used_gpu);
        let vy = r.vel_lin[0][1];
        assert!(vy >= -1e-4, "vy after resolution = {}", vy);
        assert!(r.lambda[0] >= 0.0);
    }

    #[test]
    fn test_cpu_bilateral_joint() {
        // Bilateral: two equal-mass bodies moving toward each other.
        // After solving, relative velocity along X should be ~0.
        let cfg = GpuSolverConfig {
            iterations: 40,
            omega: 1.0,
        };
        let mut solver = GpuConstraintSolver::cpu_only(cfg);

        let mut data = CpuConstraintData::new(1, 2);
        data.vel_lin[0] = [1.0, 0.0, 0.0, 1.0];
        data.vel_lin[1] = [-1.0, 0.0, 0.0, 1.0];

        data.constraints[0] = GpuConstraint {
            nx: 1.0,
            ny: 0.0,
            nz: 0.0,
            em: 0.5, // 1 / (1 + 1)
            bias: 0.0,
            lambda_lo: f32::NEG_INFINITY,
            lambda_hi: f32::MAX,
            body_a: 0,
            body_b: 1,
            rax: 0.0,
            ray: 0.0,
            raz: 0.0,
            rbx: 0.0,
            rby: 0.0,
            rbz: 0.0,
            _pad0: 0.0,
            _pad1: 0.0,
            _pad2: 0.0,
            _pad3: 0.0,
            _pad4: 0.0,
        };

        let r = solver.solve(&data);
        let rv = r.vel_lin[0][0] - r.vel_lin[1][0];
        assert!(rv.abs() < 1e-4, "residual rv = {}", rv);
    }

    #[test]
    fn test_constraint_struct_alignment() {
        // GpuConstraint must be a multiple of 16 bytes for WGSL storage buffer alignment.
        assert_eq!(
            std::mem::size_of::<GpuConstraint>() % 16,
            0,
            "GpuConstraint size {} is not 16-byte aligned",
            std::mem::size_of::<GpuConstraint>()
        );
    }

    #[test]
    fn test_solver_construction_no_panic() {
        let _s = GpuConstraintSolver::new(GpuSolverConfig::default());
    }
}
