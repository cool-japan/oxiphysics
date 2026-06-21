// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GPU-accelerated constraint solver using the wgpu compute backend.
//!
//! [`GpuConstraintSolver`] solves a batch of velocity-level constraints with a
//! projected Gauss-Seidel (PGS) sweep, structured around the
//! [`oxiphysics_gpu::compute::WgpuBackend`] buffer API.
//!
//! ## Honest status of the "GPU" path
//!
//! The [`WgpuBackend`] used here is currently a **CPU-emulation stub**: it
//! stores buffers as CPU-side `Vec<f64>` shadows and its `dispatch` is a no-op
//! (it does not execute WGSL).  `WgpuBackend::try_new` therefore returns `Err`,
//! so the default [`GpuConstraintSolver::new`] has no backend and `solve()`
//! runs the CPU PGS directly.  When a backend is injected via
//! [`GpuConstraintSolver::with_backend`], `solve()` stages the data in the
//! backend buffers and then runs the *same* PGS sweep on the CPU — it never
//! trusts the no-op `dispatch` to have solved anything — so the result is
//! identical to the CPU path and [`SolveResult::used_gpu`] honestly reports
//! `false` (no real GPU device did the work).
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
//! # On-device kernel reference
//!
//! The WGSL kernel [`WGSL_CONSTRAINT_PGS`] (`constraint_pgs_iter`) encodes one
//! full PGS iteration over N constraints: each workgroup handles 64 constraints
//! sequentially (Gauss-Seidel within the block) while across workgroups the
//! update is Jacobi-style (block-PGS, block size 64).  It is retained as the
//! reference for a future real-GPU backend (`WgpuBackendReal`); the current
//! CPU-emulation path computes the mathematically identical sweep in Rust.
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
///
/// This is the *on-device reference* only: the stub [`WgpuBackend`] does not
/// execute it.  [`GpuConstraintSolver`] computes the equivalent sweep on the
/// CPU (see the module-level "Honest status" section).
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

/// Projected Gauss-Seidel constraint solver structured around the
/// [`WgpuBackend`] buffer API.
///
/// In the default build `WgpuBackend::try_new` returns `Err`, so `solve()` runs
/// the CPU PGS.  An injected backend ([`Self::with_backend`]) routes `solve()`
/// through the buffer-staged path, which currently CPU-emulates the PGS sweep
/// (the stub backend cannot execute WGSL) and reports `used_gpu = false`.
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
    /// Create a new solver.
    ///
    /// `WgpuBackend::try_new` returns an honest `Err` in the default build (the
    /// real on-device backend is the separate `WgpuBackendReal` type), so
    /// `backend` is `None` and `solve()` runs the CPU PGS.  Use
    /// [`Self::with_backend`] to route the solver through the buffer-staged GPU
    /// path (for example in tests).
    pub fn new(config: GpuSolverConfig) -> Self {
        let backend = WgpuBackend::try_new().ok();

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

    /// Create a solver around an explicit [`WgpuBackend`], routing `solve()`
    /// through the buffer-staged GPU path.
    ///
    /// The default build's [`WgpuBackend`] is a CPU-emulation stub (its
    /// `dispatch` does not execute kernels), so this constructor is primarily
    /// used to exercise and regression-test the GPU code path: `solve()` stages
    /// the data in the backend buffers and computes the identical PGS result as
    /// [`Self::cpu_only`], reporting `used_gpu = false`.
    pub fn with_backend(config: GpuSolverConfig, backend: WgpuBackend) -> Self {
        Self {
            config,
            backend: Some(backend),
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

        // (Re-)allocate backend buffers when data size changes.
        if nc != self.last_nc || nb != self.last_nb {
            // The WgpuBackend stub stores each buffer as a CPU-side `Vec<f64>`
            // shadow.  Each GpuConstraint is flattened into CONSTRAINT_F64_SLOTS
            // f64 slots; body velocities use 4 slots each (`[x, y, z, inv_mass]`).
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

        // Stage the inputs in the backend buffers (convert f32 → f64).
        backend.write_buffer(bc, &constraints_to_f64(&data.constraints));
        backend.write_buffer(bl, &lambda_to_f64(&data.lambda));
        backend.write_buffer(bvl, &body_vels_to_f64(&data.vel_lin));
        backend.write_buffer(bva, &body_vels_to_f64(&data.vel_ang));

        // Honest compute.
        //
        // The wgpu-backend stub's `dispatch` is a no-op: it cannot execute WGSL
        // on the CPU.  Calling it here and reading the buffers back unchanged
        // would echo the inputs as a fabricated "solution".  Instead we run the
        // very sweep that `WGSL_CONSTRAINT_PGS` encodes, here in Rust, over the
        // data staged in the backend buffers.  The buffers carry the data in and
        // the solved state back out, so the result is identical to `solve_cpu`.
        //
        // When the real on-device backend (`WgpuBackendReal`) is wired in, this
        // block is replaced by a real `dispatch` of the WGSL kernel.
        let constraints = f64_to_constraints(&backend.read_buffer(bc), nc);
        let mut lambda = f64_to_f32(&backend.read_buffer(bl));
        let mut vel_lin = f64_to_body_vels(&backend.read_buffer(bvl), nb);
        let mut vel_ang = f64_to_body_vels(&backend.read_buffer(bva), nb);

        run_pgs(
            &constraints,
            &mut lambda,
            &mut vel_lin,
            &mut vel_ang,
            self.config.iterations,
            self.config.omega,
        );

        // Write the solved state back through the buffers and read it out, so the
        // returned values genuinely flow through the backend's memory path.
        backend.write_buffer(bl, &lambda_to_f64(&lambda));
        backend.write_buffer(bvl, &body_vels_to_f64(&vel_lin));
        backend.write_buffer(bva, &body_vels_to_f64(&vel_ang));

        let lambda_out = f64_to_f32(&backend.read_buffer(bl));
        let vel_lin_out = f64_to_body_vels(&backend.read_buffer(bvl), nb);
        let vel_ang_out = f64_to_body_vels(&backend.read_buffer(bva), nb);

        SolveResult {
            n_constraints: nc,
            // Honest: this `WgpuBackend` is a CPU-emulation stub, so no real GPU
            // device performed the work.  `is_available()` is `false` for the
            // stub and only becomes `true` once a real on-device backend is wired
            // in (at which point the dispatch above runs the WGSL kernel).
            used_gpu: backend.is_available(),
            vel_lin: vel_lin_out,
            vel_ang: vel_ang_out,
            lambda: lambda_out,
        }
    }

    // ── CPU fallback ──────────────────────────────────────────────────────────

    fn solve_cpu(&self, data: &CpuConstraintData) -> SolveResult {
        let nc = data.n_constraints();

        let mut lambda = data.lambda.clone();
        let mut vel_lin = data.vel_lin.clone();
        let mut vel_ang = data.vel_ang.clone();

        run_pgs(
            &data.constraints,
            &mut lambda,
            &mut vel_lin,
            &mut vel_ang,
            self.config.iterations,
            self.config.omega,
        );

        SolveResult {
            n_constraints: nc,
            used_gpu: false,
            vel_lin,
            vel_ang,
            lambda,
        }
    }
}

// ── projected Gauss-Seidel sweep ───────────────────────────────────────────────

/// Run `iterations` projected Gauss-Seidel sweeps over `constraints`, mutating
/// `lambda`, `vel_lin`, and `vel_ang` in place.
///
/// This is the single source of truth for the PGS arithmetic: both the CPU path
/// ([`GpuConstraintSolver::solve_cpu`]) and the buffer-staged GPU path
/// ([`GpuConstraintSolver::solve_gpu`], which CPU-emulates the no-op stub
/// `dispatch`) call it, guaranteeing identical results.  It mirrors the
/// [`WGSL_CONSTRAINT_PGS`] kernel.
fn run_pgs(
    constraints: &[GpuConstraint],
    lambda: &mut [f32],
    vel_lin: &mut [[f32; 4]],
    vel_ang: &mut [[f32; 4]],
    iterations: u32,
    omega: f32,
) {
    for _ in 0..iterations {
        for (ci, c) in constraints.iter().enumerate() {
            let (vla, wla, inv_ma) = gather_body(vel_lin, vel_ang, c.body_a);
            let (vlb, wlb, inv_mb) = gather_body(vel_lin, vel_ang, c.body_b);

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
                vel_lin,
                vel_ang,
                c.body_a,
                imp,
                [c.rax, c.ray, c.raz],
                inv_ma,
                1.0,
            );
            apply_impulse(
                vel_lin,
                vel_ang,
                c.body_b,
                imp,
                [c.rbx, c.rby, c.rbz],
                inv_mb,
                -1.0,
            );
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

/// Flatten a `lambda` slice (`[f32]`) to `Vec<f64>` for backend upload.
fn lambda_to_f64(lambda: &[f32]) -> Vec<f64> {
    lambda.iter().map(|&v| v as f64).collect()
}

/// Narrow a flat `Vec<f64>` downloaded from the backend back to `Vec<f32>`.
fn f64_to_f32(flat: &[f64]) -> Vec<f32> {
    flat.iter().map(|&v| v as f32).collect()
}

/// Reconstruct `Vec<GpuConstraint>` from the flat `Vec<f64>` produced by
/// [`constraints_to_f64`].  This is the exact inverse: because every
/// [`GpuConstraint`] field is `f32` (or a `u32` index ≤ 2^32), the `f32 → f64`
/// upload and `f64 → f32` reconstruction are lossless, so the recovered
/// constraints are bit-identical to the originals.
fn f64_to_constraints(flat: &[f64], n: usize) -> Vec<GpuConstraint> {
    (0..n)
        .map(|i| {
            let b = i * CONSTRAINT_F64_SLOTS;
            GpuConstraint {
                nx: flat[b] as f32,
                ny: flat[b + 1] as f32,
                nz: flat[b + 2] as f32,
                em: flat[b + 3] as f32,
                bias: flat[b + 4] as f32,
                lambda_lo: flat[b + 5] as f32,
                lambda_hi: flat[b + 6] as f32,
                body_a: flat[b + 7] as u32,
                body_b: flat[b + 8] as u32,
                rax: flat[b + 9] as f32,
                ray: flat[b + 10] as f32,
                raz: flat[b + 11] as f32,
                rbx: flat[b + 12] as f32,
                rby: flat[b + 13] as f32,
                rbz: flat[b + 14] as f32,
                _pad0: flat[b + 15] as f32,
                _pad1: flat[b + 16] as f32,
                _pad2: flat[b + 17] as f32,
                _pad3: flat[b + 18] as f32,
                _pad4: flat[b + 19] as f32,
            }
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

    #[test]
    fn test_gpu_path_matches_cpu() {
        // Regression guard against the silent fabrication where `solve_gpu`
        // trusted the no-op stub `dispatch` and returned the uploaded inputs as
        // a fake "solution".  We force the GPU path by injecting a CPU-emulation
        // stub backend and assert that (a) it actually changes the state and
        // (b) it matches the CPU reference solver bit-for-bit.
        let cfg = GpuSolverConfig {
            iterations: 16,
            omega: 1.2,
        };

        // A non-trivial system: 3 bodies, 2 contacts + 1 bilateral, warm-started
        // with non-zero lambda and non-zero angular velocities.
        let mut data = CpuConstraintData::new(3, 3);
        data.vel_lin[0] = [0.0, -2.0, 0.0, 1.0];
        data.vel_lin[1] = [0.5, -1.0, 0.0, 0.5];
        data.vel_lin[2] = [-0.3, 0.7, 0.2, 2.0];
        data.vel_ang[0] = [0.0, 0.0, 0.1, 0.0];
        data.vel_ang[1] = [0.2, 0.0, 0.0, 0.0];
        data.vel_ang[2] = [0.0, -0.1, 0.0, 0.0];
        data.lambda = vec![0.3, 0.0, -0.5];

        // Contact: body 0 vs static world.
        data.constraints[0] = GpuConstraint::contact(GpuConstraintParams {
            nx: 0.0,
            ny: 1.0,
            nz: 0.0,
            effective_mass: 0.8,
            bias: -0.05,
            body_a: 0,
            body_b: u32::MAX,
            rax: 0.1,
            ray: 0.0,
            raz: -0.2,
            rbx: 0.0,
            rby: 0.0,
            rbz: 0.0,
        });
        // Contact: body 1 vs body 2.
        data.constraints[1] = GpuConstraint::contact(GpuConstraintParams {
            nx: 0.6,
            ny: 0.8,
            nz: 0.0,
            effective_mass: 0.5,
            bias: 0.0,
            body_a: 1,
            body_b: 2,
            rax: 0.0,
            ray: 0.3,
            raz: 0.0,
            rbx: -0.1,
            rby: 0.0,
            rbz: 0.05,
        });
        // Bilateral: body 0 vs body 2.
        data.constraints[2] = GpuConstraint::bilateral(GpuConstraintParams {
            nx: 1.0,
            ny: 0.0,
            nz: 0.0,
            effective_mass: 0.4,
            bias: 0.0,
            body_a: 0,
            body_b: 2,
            rax: 0.0,
            ray: 0.0,
            raz: 0.0,
            rbx: 0.0,
            rby: 0.0,
            rbz: 0.0,
        });

        // Forced GPU path (injected stub backend) vs the CPU reference.
        let mut gpu_solver =
            GpuConstraintSolver::with_backend(cfg.clone(), WgpuBackend::new_stub());
        assert!(
            gpu_solver.has_gpu(),
            "with_backend must enable the GPU path"
        );
        let mut cpu_solver = GpuConstraintSolver::cpu_only(cfg.clone());

        let rg = gpu_solver.solve(&data);
        let rc = cpu_solver.solve(&data);

        // The stub backend is CPU emulation, so `used_gpu` must honestly be false.
        assert!(
            !rg.used_gpu,
            "CPU-emulation stub backend must report used_gpu = false"
        );

        // It must NOT be an echo of the inputs: the solver changed the state.
        let changed =
            rg.lambda != data.lambda || rg.vel_lin != data.vel_lin || rg.vel_ang != data.vel_ang;
        assert!(
            changed,
            "GPU path returned inputs unchanged — silent fabrication"
        );

        // And it must match the CPU reference exactly (same arithmetic, lossless
        // f32↔f64 buffer round-trip).
        assert_eq!(rg.n_constraints, rc.n_constraints);
        assert_eq!(rg.lambda.len(), rc.lambda.len());
        for (g, c) in rg.lambda.iter().zip(rc.lambda.iter()) {
            assert!((g - c).abs() <= 1e-6, "lambda mismatch: gpu={g} cpu={c}");
        }
        for (gv, cv) in rg.vel_lin.iter().zip(rc.vel_lin.iter()) {
            for (g, c) in gv.iter().zip(cv.iter()) {
                assert!((g - c).abs() <= 1e-6, "vel_lin mismatch: gpu={g} cpu={c}");
            }
        }
        for (gv, cv) in rg.vel_ang.iter().zip(rc.vel_ang.iter()) {
            for (g, c) in gv.iter().zip(cv.iter()) {
                assert!((g - c).abs() <= 1e-6, "vel_ang mismatch: gpu={g} cpu={c}");
            }
        }
    }
}
