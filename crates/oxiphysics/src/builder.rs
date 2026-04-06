// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fluent builder API for configuring physics simulations.
//!
//! This module provides a set of builder types that follow the
//! *builder pattern*, letting users compose rich simulation configurations
//! without needing to construct deeply nested structs by hand.
//!
//! # Quick-start
//!
//! ```no_run
//! use oxiphysics::builder::{SimulationBuilder, RigidBodyBuilder, BroadphaseType};
//!
//! let cfg = SimulationBuilder::new()
//!     .with_gravity([0.0, -9.81, 0.0])
//!     .with_time_step(1.0 / 60.0)
//!     .with_substeps(4)
//!     .with_collision_enabled(true)
//!     .with_solver_iterations(10, 5)
//!     .with_broadphase_type(BroadphaseType::SweepAndPrune)
//!     .build();
//!
//! let body = RigidBodyBuilder::sphere(0.5)
//!     .with_position([0.0, 5.0, 0.0])
//!     .with_density(1000.0)
//!     .with_restitution(0.3)
//!     .build();
//!
//! assert!((cfg.gravity[1] + 9.81).abs() < 1e-12);
//! assert!((body.radius - 0.5).abs() < 1e-12);
//! ```

// ---------------------------------------------------------------------------
// BroadphaseType
// ---------------------------------------------------------------------------

/// Selects the broadphase collision-detection algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BroadphaseType {
    /// Dynamic AABB tree — good for scenes with many moving objects.
    DynamicAabb,
    /// Sweep-and-prune along a primary axis — efficient for convex bodies.
    SweepAndPrune,
    /// Brute-force O(n²) check — useful for tiny scenes or debugging.
    BruteForce,
}

impl Default for BroadphaseType {
    fn default() -> Self {
        Self::SweepAndPrune
    }
}

// ---------------------------------------------------------------------------
// SimulationConfig / SimulationBuilder
// ---------------------------------------------------------------------------

/// Fully resolved configuration for a physics simulation.
///
/// Produced by [`SimulationBuilder::build`].
#[derive(Debug, Clone)]
pub struct SimulationConfig {
    /// Gravitational acceleration vector (m s⁻²).
    pub gravity: [f64; 3],
    /// Fixed simulation time-step (s).
    pub time_step: f64,
    /// Number of sub-steps per time-step.
    pub substeps: usize,
    /// Whether collision detection is active.
    pub collision_enabled: bool,
    /// Whether bodies may go to sleep when idle.
    pub sleeping_enabled: bool,
    /// Number of position-correction solver iterations per sub-step.
    pub solver_position_iterations: usize,
    /// Number of velocity-correction solver iterations per sub-step.
    pub solver_velocity_iterations: usize,
    /// Coefficient of restitution (bounciness) applied globally.
    pub restitution: f64,
    /// Coefficient of friction applied globally.
    pub friction: f64,
    /// Broadphase algorithm to use.
    pub broadphase_type: BroadphaseType,
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            time_step: 1.0 / 60.0,
            substeps: 1,
            collision_enabled: true,
            sleeping_enabled: true,
            solver_position_iterations: 10,
            solver_velocity_iterations: 5,
            restitution: 0.3,
            friction: 0.5,
            broadphase_type: BroadphaseType::SweepAndPrune,
        }
    }
}

/// Fluent builder for [`SimulationConfig`].
///
/// All methods return `self` for chaining.  Call [`SimulationBuilder::build`]
/// to produce a [`SimulationConfig`].
///
/// # Example
/// ```no_run
/// use oxiphysics::builder::{SimulationBuilder, BroadphaseType};
/// let cfg = SimulationBuilder::new()
///     .with_gravity([0.0, -9.81, 0.0])
///     .with_time_step(0.01)
///     .build();
/// assert!((cfg.time_step - 0.01).abs() < 1e-12);
/// ```
#[derive(Debug, Clone, Default)]
pub struct SimulationBuilder {
    config: SimulationConfig,
}

impl SimulationBuilder {
    /// Create a new builder pre-populated with sensible defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the global gravity vector (m s⁻²).
    pub fn with_gravity(mut self, g: [f64; 3]) -> Self {
        self.config.gravity = g;
        self
    }

    /// Set the fixed simulation time-step (s).
    ///
    /// Must be positive.
    pub fn with_time_step(mut self, dt: f64) -> Self {
        assert!(dt > 0.0, "time_step must be positive");
        self.config.time_step = dt;
        self
    }

    /// Set the number of physics sub-steps per call to `step()`.
    ///
    /// More sub-steps improve accuracy at the cost of performance.
    pub fn with_substeps(mut self, n: usize) -> Self {
        assert!(n >= 1, "substeps must be >= 1");
        self.config.substeps = n;
        self
    }

    /// Enable or disable collision detection.
    pub fn with_collision_enabled(mut self, b: bool) -> Self {
        self.config.collision_enabled = b;
        self
    }

    /// Enable or disable body sleeping.
    ///
    /// Sleeping bodies are not integrated and skip collision tests,
    /// saving CPU time.
    pub fn with_sleeping_enabled(mut self, b: bool) -> Self {
        self.config.sleeping_enabled = b;
        self
    }

    /// Set the number of position and velocity solver iterations per sub-step.
    pub fn with_solver_iterations(mut self, pos: usize, vel: usize) -> Self {
        self.config.solver_position_iterations = pos;
        self.config.solver_velocity_iterations = vel;
        self
    }

    /// Set the global coefficient of restitution (0 = inelastic, 1 = elastic).
    pub fn with_restitution(mut self, e: f64) -> Self {
        self.config.restitution = e.clamp(0.0, 1.0);
        self
    }

    /// Set the global coefficient of friction (0 = frictionless).
    pub fn with_friction(mut self, mu: f64) -> Self {
        self.config.friction = mu.max(0.0);
        self
    }

    /// Choose the broadphase algorithm.
    pub fn with_broadphase_type(mut self, t: BroadphaseType) -> Self {
        self.config.broadphase_type = t;
        self
    }

    /// Consume the builder and return the resolved [`SimulationConfig`].
    pub fn build(self) -> SimulationConfig {
        self.config
    }
}

// ---------------------------------------------------------------------------
// ShapeKind — internal helper
// ---------------------------------------------------------------------------

/// Describes the shape of a rigid body being built.
#[derive(Debug, Clone, PartialEq)]
pub enum ShapeKind {
    /// Sphere with the given radius (m).
    Sphere {
        /// Sphere radius in metres.
        radius: f64,
    },
    /// Box with half-extents along each axis (m).
    Cuboid {
        /// Half-extents \[hx, hy, hz\] in metres.
        half_extents: [f64; 3],
    },
    /// Capsule with a cylindrical shaft.
    Capsule {
        /// Capsule end-cap radius in metres.
        radius: f64,
        /// Half-height of the cylindrical shaft in metres.
        half_height: f64,
    },
}

// ---------------------------------------------------------------------------
// BodyKind — static / kinematic / dynamic
// ---------------------------------------------------------------------------

/// Motion type of a rigid body.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BodyKind {
    /// Infinitely massive body that never moves.
    Static,
    /// Moves according to an externally provided trajectory.
    Kinematic,
    /// Fully simulated body obeying Newton's laws.
    Dynamic,
}

impl Default for BodyKind {
    fn default() -> Self {
        Self::Dynamic
    }
}

// ---------------------------------------------------------------------------
// MassSpec — mass or density
// ---------------------------------------------------------------------------

/// Specifies the mass of a rigid body either directly or via density.
#[derive(Debug, Clone)]
pub enum MassSpec {
    /// Explicit total mass (kg).
    Mass(f64),
    /// Density (kg m⁻³); mass is inferred from the shape volume.
    Density(f64),
}

impl Default for MassSpec {
    fn default() -> Self {
        Self::Density(1000.0)
    }
}

// ---------------------------------------------------------------------------
// RigidBodyConfig / RigidBodyBuilder
// ---------------------------------------------------------------------------

/// Fully resolved configuration for a single rigid body.
///
/// Produced by [`RigidBodyBuilder::build`].
#[derive(Debug, Clone)]
pub struct RigidBodyConfig {
    /// Shape and its defining parameters.
    pub shape: ShapeKind,
    /// Convenience accessor — radius for sphere / capsule shapes.
    pub radius: f64,
    /// World-space position (m).
    pub position: [f64; 3],
    /// Euler angles describing the initial orientation (roll, pitch, yaw) (rad).
    pub rotation_euler: [f64; 3],
    /// Initial linear velocity (m s⁻¹).
    pub linear_velocity: [f64; 3],
    /// Initial angular velocity (rad s⁻¹).
    pub angular_velocity: [f64; 3],
    /// Mass specification (explicit or via density).
    pub mass_spec: MassSpec,
    /// Motion type.
    pub kind: BodyKind,
    /// Coefficient of restitution for this body (overrides global if ≥ 0).
    pub restitution: f64,
    /// Coefficient of friction for this body (overrides global if ≥ 0).
    pub friction: f64,
    /// Linear velocity damping coefficient.
    pub linear_damping: f64,
    /// Angular velocity damping coefficient.
    pub angular_damping: f64,
}

/// Fluent builder for a single [`RigidBodyConfig`].
///
/// Start with one of the shape factory methods:
/// - [`RigidBodyBuilder::sphere`]
/// - [`RigidBodyBuilder::cuboid`]
/// - [`RigidBodyBuilder::capsule`]
///
/// Then chain optional setters before calling [`RigidBodyBuilder::build`].
#[derive(Debug, Clone)]
pub struct RigidBodyBuilder {
    shape: ShapeKind,
    radius: f64,
    position: [f64; 3],
    rotation_euler: [f64; 3],
    linear_velocity: [f64; 3],
    angular_velocity: [f64; 3],
    mass_spec: MassSpec,
    kind: BodyKind,
    restitution: f64,
    friction: f64,
    linear_damping: f64,
    angular_damping: f64,
}

impl RigidBodyBuilder {
    // ---- shape factories ----

    /// Start building a sphere rigid body with the given radius (m).
    pub fn sphere(radius: f64) -> Self {
        Self::new_inner(ShapeKind::Sphere { radius }, radius)
    }

    /// Start building a box (cuboid) rigid body with the given half-extents (m).
    pub fn cuboid(half_extents: [f64; 3]) -> Self {
        Self::new_inner(ShapeKind::Cuboid { half_extents }, half_extents[0])
    }

    /// Start building a capsule rigid body.
    pub fn capsule(radius: f64, half_height: f64) -> Self {
        Self::new_inner(ShapeKind::Capsule { radius, half_height }, radius)
    }

    fn new_inner(shape: ShapeKind, radius: f64) -> Self {
        Self {
            shape,
            radius,
            position: [0.0; 3],
            rotation_euler: [0.0; 3],
            linear_velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            mass_spec: MassSpec::default(),
            kind: BodyKind::default(),
            restitution: -1.0, // -1 means "use global default"
            friction: -1.0,
            linear_damping: 0.0,
            angular_damping: 0.0,
        }
    }

    // ---- setters ----

    /// Set the world-space initial position (m).
    pub fn with_position(mut self, pos: [f64; 3]) -> Self {
        self.position = pos;
        self
    }

    /// Set the initial orientation as Euler angles (roll, pitch, yaw) in radians.
    pub fn with_rotation_euler(mut self, roll: f64, pitch: f64, yaw: f64) -> Self {
        self.rotation_euler = [roll, pitch, yaw];
        self
    }

    /// Set the initial linear velocity (m s⁻¹).
    pub fn with_linear_velocity(mut self, v: [f64; 3]) -> Self {
        self.linear_velocity = v;
        self
    }

    /// Set the initial angular velocity (rad s⁻¹).
    pub fn with_angular_velocity(mut self, w: [f64; 3]) -> Self {
        self.angular_velocity = w;
        self
    }

    /// Specify mass via density (kg m⁻³).  The actual mass is computed from
    /// the shape volume.
    pub fn with_density(mut self, rho: f64) -> Self {
        self.mass_spec = MassSpec::Density(rho);
        self
    }

    /// Specify mass directly (kg).
    pub fn with_mass(mut self, m: f64) -> Self {
        self.mass_spec = MassSpec::Mass(m);
        self
    }

    /// Make the body static (infinite mass, never moves).
    pub fn is_static(mut self) -> Self {
        self.kind = BodyKind::Static;
        self
    }

    /// Make the body kinematic (moves by assignment, not by forces).
    pub fn is_kinematic(mut self) -> Self {
        self.kind = BodyKind::Kinematic;
        self
    }

    /// Override the global coefficient of restitution for this body.
    pub fn with_restitution(mut self, e: f64) -> Self {
        self.restitution = e.clamp(0.0, 1.0);
        self
    }

    /// Override the global coefficient of friction for this body.
    pub fn with_friction(mut self, mu: f64) -> Self {
        self.friction = mu.max(0.0);
        self
    }

    /// Set linear and angular damping coefficients.
    pub fn with_damping(mut self, lin: f64, ang: f64) -> Self {
        self.linear_damping = lin.max(0.0);
        self.angular_damping = ang.max(0.0);
        self
    }

    /// Consume the builder and return the resolved [`RigidBodyConfig`].
    pub fn build(self) -> RigidBodyConfig {
        RigidBodyConfig {
            shape: self.shape,
            radius: self.radius,
            position: self.position,
            rotation_euler: self.rotation_euler,
            linear_velocity: self.linear_velocity,
            angular_velocity: self.angular_velocity,
            mass_spec: self.mass_spec,
            kind: self.kind,
            restitution: self.restitution,
            friction: self.friction,
            linear_damping: self.linear_damping,
            angular_damping: self.angular_damping,
        }
    }
}

// ---------------------------------------------------------------------------
// FluidConfig / FluidSimBuilder
// ---------------------------------------------------------------------------

/// Selects the LBM collision operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LbmCollisionOp {
    /// Single-relaxation-time (BGK) operator.
    Bgk,
    /// Multiple-relaxation-time (MRT) operator.
    Mrt,
    /// Two-relaxation-time (TRT) operator.
    Trt,
}

impl Default for LbmCollisionOp {
    fn default() -> Self {
        Self::Bgk
    }
}

/// Selects the dimensionality of the LBM grid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LbmDimension {
    /// Two-dimensional lattice.
    D2,
    /// Three-dimensional lattice.
    D3,
}

/// Fully resolved configuration for an LBM fluid simulation.
///
/// Produced by [`FluidSimBuilder::build`].
#[derive(Debug, Clone)]
pub struct FluidConfig {
    /// Grid resolution along X.
    pub nx: usize,
    /// Grid resolution along Y.
    pub ny: usize,
    /// Grid resolution along Z (1 for 2-D simulations).
    pub nz: usize,
    /// Physical cell size (m).
    pub dx: f64,
    /// Kinematic viscosity (m² s⁻¹).
    pub viscosity: f64,
    /// Inlet velocity components (m s⁻¹).
    pub inlet_velocity: [f64; 3],
    /// LBM collision operator.
    pub collision_op: LbmCollisionOp,
    /// Grid dimensionality.
    pub dimension: LbmDimension,
}

/// Fluent builder for an LBM [`FluidConfig`].
///
/// # Example
/// ```no_run
/// use oxiphysics::builder::FluidSimBuilder;
/// let cfg = FluidSimBuilder::lbm_2d(64, 32, 0.01)
///     .with_viscosity(1e-6)
///     .with_inlet_velocity(0.1, 0.0)
///     .with_mrt()
///     .build();
/// assert_eq!(cfg.nx, 64);
/// ```
#[derive(Debug, Clone)]
pub struct FluidSimBuilder {
    nx: usize,
    ny: usize,
    nz: usize,
    dx: f64,
    viscosity: f64,
    inlet_velocity: [f64; 3],
    collision_op: LbmCollisionOp,
    dimension: LbmDimension,
}

impl FluidSimBuilder {
    /// Start building a 2-D LBM simulation on an `nx × ny` grid with cell size `dx`.
    pub fn lbm_2d(nx: usize, ny: usize, dx: f64) -> Self {
        Self {
            nx,
            ny,
            nz: 1,
            dx,
            viscosity: 1e-6,
            inlet_velocity: [0.0; 3],
            collision_op: LbmCollisionOp::default(),
            dimension: LbmDimension::D2,
        }
    }

    /// Start building a 3-D LBM simulation on an `nx × ny × nz` grid.
    pub fn lbm_3d(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        Self {
            nx,
            ny,
            nz,
            dx,
            viscosity: 1e-6,
            inlet_velocity: [0.0; 3],
            collision_op: LbmCollisionOp::default(),
            dimension: LbmDimension::D3,
        }
    }

    /// Set kinematic viscosity (m² s⁻¹).
    pub fn with_viscosity(mut self, nu: f64) -> Self {
        self.viscosity = nu;
        self
    }

    /// Set inlet velocity components (m s⁻¹).  The z component is ignored in 2-D.
    pub fn with_inlet_velocity(mut self, vx: f64, vy: f64) -> Self {
        self.inlet_velocity[0] = vx;
        self.inlet_velocity[1] = vy;
        self
    }

    /// Use the MRT collision operator.
    pub fn with_mrt(mut self) -> Self {
        self.collision_op = LbmCollisionOp::Mrt;
        self
    }

    /// Use the TRT collision operator.
    pub fn with_trt(mut self) -> Self {
        self.collision_op = LbmCollisionOp::Trt;
        self
    }

    /// Consume the builder and return the resolved [`FluidConfig`].
    pub fn build(self) -> FluidConfig {
        FluidConfig {
            nx: self.nx,
            ny: self.ny,
            nz: self.nz,
            dx: self.dx,
            viscosity: self.viscosity,
            inlet_velocity: self.inlet_velocity,
            collision_op: self.collision_op,
            dimension: self.dimension,
        }
    }
}

// ---------------------------------------------------------------------------
// SphConfig / SphSimBuilder
// ---------------------------------------------------------------------------

/// SPH pressure solver variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SphPressureSolver {
    /// Weakly Compressible SPH (equation of state).
    WcSph,
    /// Implicit Incompressible SPH (IISPH).
    Iisph,
    /// Predictive-corrective incompressible SPH (PCISPH).
    Pcisph,
}

impl Default for SphPressureSolver {
    fn default() -> Self {
        Self::WcSph
    }
}

/// Fully resolved configuration for an SPH simulation.
///
/// Produced by [`SphSimBuilder::build`].
#[derive(Debug, Clone)]
pub struct SphConfig {
    /// Number of fluid particles.
    pub n_particles: usize,
    /// Initial inter-particle spacing (m).
    pub particle_spacing: f64,
    /// Smoothing kernel radius (m).
    pub kernel_radius: f64,
    /// Pressure solver variant.
    pub pressure_solver: SphPressureSolver,
    /// Reference fluid density (kg m⁻³).
    pub rest_density: f64,
    /// Dynamic viscosity (Pa s).
    pub viscosity: f64,
}

/// Fluent builder for an SPH [`SphConfig`].
///
/// # Example
/// ```no_run
/// use oxiphysics::builder::{SphSimBuilder, SphPressureSolver};
/// let cfg = SphSimBuilder::new(1000)
///     .with_particle_spacing(0.02)
///     .with_kernel_radius(0.04)
///     .with_pressure_solver_iisph()
///     .build();
/// assert_eq!(cfg.pressure_solver, SphPressureSolver::Iisph);
/// ```
#[derive(Debug, Clone)]
pub struct SphSimBuilder {
    n_particles: usize,
    particle_spacing: f64,
    kernel_radius: f64,
    pressure_solver: SphPressureSolver,
    rest_density: f64,
    viscosity: f64,
}

impl SphSimBuilder {
    /// Create a new SPH builder for `n_particles` particles.
    pub fn new(n_particles: usize) -> Self {
        Self {
            n_particles,
            particle_spacing: 0.02,
            kernel_radius: 0.04,
            pressure_solver: SphPressureSolver::default(),
            rest_density: 1000.0,
            viscosity: 1e-3,
        }
    }

    /// Set the initial inter-particle spacing (m).
    pub fn with_particle_spacing(mut self, h: f64) -> Self {
        self.particle_spacing = h;
        self
    }

    /// Set the SPH smoothing kernel radius (m).
    pub fn with_kernel_radius(mut self, r: f64) -> Self {
        self.kernel_radius = r;
        self
    }

    /// Use the IISPH pressure solver.
    pub fn with_pressure_solver_iisph(mut self) -> Self {
        self.pressure_solver = SphPressureSolver::Iisph;
        self
    }

    /// Use the PCISPH pressure solver.
    pub fn with_pressure_solver_pcisph(mut self) -> Self {
        self.pressure_solver = SphPressureSolver::Pcisph;
        self
    }

    /// Set the reference rest density (kg m⁻³).
    pub fn with_rest_density(mut self, rho: f64) -> Self {
        self.rest_density = rho;
        self
    }

    /// Set the dynamic viscosity (Pa s).
    pub fn with_viscosity(mut self, mu: f64) -> Self {
        self.viscosity = mu;
        self
    }

    /// Consume the builder and return the resolved [`SphConfig`].
    pub fn build(self) -> SphConfig {
        SphConfig {
            n_particles: self.n_particles,
            particle_spacing: self.particle_spacing,
            kernel_radius: self.kernel_radius,
            pressure_solver: self.pressure_solver,
            rest_density: self.rest_density,
            viscosity: self.viscosity,
        }
    }
}

// ---------------------------------------------------------------------------
// FemConfig / FemBuilder
// ---------------------------------------------------------------------------

/// Describes how a Dirichlet boundary condition is applied to an FEM node.
#[derive(Debug, Clone)]
pub struct DirichletBc {
    /// Index of the constrained node.
    pub node: usize,
    /// Degree of freedom index (0 = x, 1 = y, 2 = z).
    pub dof: usize,
    /// Prescribed displacement value (m).
    pub value: f64,
}

/// Element type in the FEM mesh.
#[derive(Debug, Clone, PartialEq)]
pub enum FemElement {
    /// Linear tetrahedral element (4 nodes).
    Tet([usize; 4]),
    /// Trilinear hexahedral element (8 nodes).
    Hex([usize; 8]),
}

/// FEM elastic material specification.
#[derive(Debug, Clone)]
pub struct FemMaterialSpec {
    /// Young's modulus (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio (dimensionless, 0 < ν < 0.5).
    pub poisson_ratio: f64,
    /// Mass density (kg m⁻³).
    pub density: f64,
}

/// Fully resolved configuration for an FEM simulation.
///
/// Produced by [`FemBuilder::build`].
#[derive(Debug, Clone)]
pub struct FemConfig {
    /// Node positions (m), stored as flat `[x, y, z]` triples.
    pub nodes: Vec<[f64; 3]>,
    /// Element connectivity list.
    pub elements: Vec<FemElement>,
    /// Dirichlet boundary conditions.
    pub dirichlet_bcs: Vec<DirichletBc>,
    /// Material specification (currently a single global material).
    pub material: Option<FemMaterialSpec>,
}

/// Fluent builder for an FEM [`FemConfig`].
///
/// # Example
/// ```no_run
/// use oxiphysics::builder::FemBuilder;
/// let cfg = FemBuilder::new()
///     .add_node([0.0, 0.0, 0.0])
///     .add_node([1.0, 0.0, 0.0])
///     .add_node([0.0, 1.0, 0.0])
///     .add_node([0.0, 0.0, 1.0])
///     .add_tet([0, 1, 2, 3])
///     .with_material_elastic(200e9, 0.3, 7800.0)
///     .with_dirichlet(0, 0, 0.0)
///     .build();
/// assert_eq!(cfg.nodes.len(), 4);
/// ```
#[derive(Debug, Clone, Default)]
pub struct FemBuilder {
    nodes: Vec<[f64; 3]>,
    elements: Vec<FemElement>,
    dirichlet_bcs: Vec<DirichletBc>,
    material: Option<FemMaterialSpec>,
}

impl FemBuilder {
    /// Create an empty FEM builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a node at the given position (m).
    pub fn add_node(mut self, pos: [f64; 3]) -> Self {
        self.nodes.push(pos);
        self
    }

    /// Add a linear tetrahedral element defined by four node indices.
    pub fn add_tet(mut self, i: [usize; 4]) -> Self {
        self.elements.push(FemElement::Tet(i));
        self
    }

    /// Add a trilinear hexahedral element defined by eight node indices.
    pub fn add_hex(mut self, i: [usize; 8]) -> Self {
        self.elements.push(FemElement::Hex(i));
        self
    }

    /// Specify a linear elastic material with Young's modulus `e_modulus` (Pa),
    /// Poisson's ratio `nu`, and density `rho` (kg m⁻³).
    pub fn with_material_elastic(mut self, e_modulus: f64, nu: f64, rho: f64) -> Self {
        self.material = Some(FemMaterialSpec {
            young_modulus: e_modulus,
            poisson_ratio: nu,
            density: rho,
        });
        self
    }

    /// Add a Dirichlet boundary condition (prescribed displacement).
    pub fn with_dirichlet(mut self, node: usize, dof: usize, value: f64) -> Self {
        self.dirichlet_bcs.push(DirichletBc { node, dof, value });
        self
    }

    /// Consume the builder and return the resolved [`FemConfig`].
    pub fn build(self) -> FemConfig {
        FemConfig {
            nodes: self.nodes,
            elements: self.elements,
            dirichlet_bcs: self.dirichlet_bcs,
            material: self.material,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- SimulationBuilder ---

    #[test]
    fn test_simulation_builder_defaults() {
        let cfg = SimulationBuilder::new().build();
        assert!((cfg.gravity[1] + 9.81).abs() < 1e-12);
        assert!((cfg.time_step - 1.0 / 60.0).abs() < 1e-10);
        assert_eq!(cfg.substeps, 1);
        assert!(cfg.collision_enabled);
        assert!(cfg.sleeping_enabled);
    }

    #[test]
    fn test_simulation_builder_gravity() {
        let cfg = SimulationBuilder::new()
            .with_gravity([0.0, -1.62, 0.0])
            .build();
        assert!((cfg.gravity[1] + 1.62).abs() < 1e-12);
    }

    #[test]
    fn test_simulation_builder_zero_gravity() {
        let cfg = SimulationBuilder::new()
            .with_gravity([0.0, 0.0, 0.0])
            .build();
        for c in cfg.gravity.iter() {
            assert!(c.abs() < 1e-12);
        }
    }

    #[test]
    fn test_simulation_builder_time_step() {
        let cfg = SimulationBuilder::new().with_time_step(0.002).build();
        assert!((cfg.time_step - 0.002).abs() < 1e-12);
    }

    #[test]
    fn test_simulation_builder_substeps() {
        let cfg = SimulationBuilder::new().with_substeps(8).build();
        assert_eq!(cfg.substeps, 8);
    }

    #[test]
    fn test_simulation_builder_collision_disabled() {
        let cfg = SimulationBuilder::new()
            .with_collision_enabled(false)
            .build();
        assert!(!cfg.collision_enabled);
    }

    #[test]
    fn test_simulation_builder_sleeping_disabled() {
        let cfg = SimulationBuilder::new()
            .with_sleeping_enabled(false)
            .build();
        assert!(!cfg.sleeping_enabled);
    }

    #[test]
    fn test_simulation_builder_solver_iterations() {
        let cfg = SimulationBuilder::new()
            .with_solver_iterations(20, 10)
            .build();
        assert_eq!(cfg.solver_position_iterations, 20);
        assert_eq!(cfg.solver_velocity_iterations, 10);
    }

    #[test]
    fn test_simulation_builder_restitution() {
        let cfg = SimulationBuilder::new().with_restitution(0.7).build();
        assert!((cfg.restitution - 0.7).abs() < 1e-12);
    }

    #[test]
    fn test_simulation_builder_restitution_clamped_high() {
        let cfg = SimulationBuilder::new().with_restitution(1.5).build();
        assert!((cfg.restitution - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_simulation_builder_restitution_clamped_low() {
        let cfg = SimulationBuilder::new().with_restitution(-0.5).build();
        assert!((cfg.restitution - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_simulation_builder_friction() {
        let cfg = SimulationBuilder::new().with_friction(0.8).build();
        assert!((cfg.friction - 0.8).abs() < 1e-12);
    }

    #[test]
    fn test_simulation_builder_friction_clamped() {
        let cfg = SimulationBuilder::new().with_friction(-1.0).build();
        assert!((cfg.friction - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_simulation_builder_broadphase_dynamic_aabb() {
        let cfg = SimulationBuilder::new()
            .with_broadphase_type(BroadphaseType::DynamicAabb)
            .build();
        assert_eq!(cfg.broadphase_type, BroadphaseType::DynamicAabb);
    }

    #[test]
    fn test_simulation_builder_broadphase_brute_force() {
        let cfg = SimulationBuilder::new()
            .with_broadphase_type(BroadphaseType::BruteForce)
            .build();
        assert_eq!(cfg.broadphase_type, BroadphaseType::BruteForce);
    }

    #[test]
    fn test_simulation_builder_chain() {
        let cfg = SimulationBuilder::new()
            .with_gravity([0.0, -9.81, 0.0])
            .with_time_step(0.01)
            .with_substeps(4)
            .with_collision_enabled(true)
            .with_sleeping_enabled(true)
            .with_solver_iterations(15, 8)
            .with_restitution(0.4)
            .with_friction(0.6)
            .with_broadphase_type(BroadphaseType::SweepAndPrune)
            .build();
        assert_eq!(cfg.substeps, 4);
        assert_eq!(cfg.solver_position_iterations, 15);
    }

    // --- RigidBodyBuilder ---

    #[test]
    fn test_rigid_body_builder_sphere() {
        let body = RigidBodyBuilder::sphere(0.5).build();
        assert!((body.radius - 0.5).abs() < 1e-12);
        assert!(matches!(body.shape, ShapeKind::Sphere { .. }));
    }

    #[test]
    fn test_rigid_body_builder_cuboid() {
        let body = RigidBodyBuilder::cuboid([1.0, 0.5, 0.25]).build();
        if let ShapeKind::Cuboid { half_extents } = body.shape {
            assert!((half_extents[0] - 1.0).abs() < 1e-12);
        } else {
            panic!("Expected Cuboid");
        }
    }

    #[test]
    fn test_rigid_body_builder_capsule() {
        let body = RigidBodyBuilder::capsule(0.3, 0.7).build();
        if let ShapeKind::Capsule { radius, half_height } = body.shape {
            assert!((radius - 0.3).abs() < 1e-12);
            assert!((half_height - 0.7).abs() < 1e-12);
        } else {
            panic!("Expected Capsule");
        }
    }

    #[test]
    fn test_rigid_body_builder_position() {
        let body = RigidBodyBuilder::sphere(1.0)
            .with_position([1.0, 2.0, 3.0])
            .build();
        assert_eq!(body.position, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_rigid_body_builder_rotation() {
        let body = RigidBodyBuilder::sphere(1.0)
            .with_rotation_euler(0.1, 0.2, 0.3)
            .build();
        assert!((body.rotation_euler[0] - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_rigid_body_builder_linear_velocity() {
        let body = RigidBodyBuilder::sphere(1.0)
            .with_linear_velocity([1.0, 0.0, 0.0])
            .build();
        assert_eq!(body.linear_velocity, [1.0, 0.0, 0.0]);
    }

    #[test]
    fn test_rigid_body_builder_angular_velocity() {
        let body = RigidBodyBuilder::sphere(1.0)
            .with_angular_velocity([0.0, 1.0, 0.0])
            .build();
        assert_eq!(body.angular_velocity, [0.0, 1.0, 0.0]);
    }

    #[test]
    fn test_rigid_body_builder_density() {
        let body = RigidBodyBuilder::sphere(1.0).with_density(500.0).build();
        assert!(matches!(body.mass_spec, MassSpec::Density(d) if (d - 500.0).abs() < 1e-12));
    }

    #[test]
    fn test_rigid_body_builder_mass() {
        let body = RigidBodyBuilder::sphere(1.0).with_mass(10.0).build();
        assert!(matches!(body.mass_spec, MassSpec::Mass(m) if (m - 10.0).abs() < 1e-12));
    }

    #[test]
    fn test_rigid_body_builder_static() {
        let body = RigidBodyBuilder::sphere(1.0).is_static().build();
        assert_eq!(body.kind, BodyKind::Static);
    }

    #[test]
    fn test_rigid_body_builder_kinematic() {
        let body = RigidBodyBuilder::sphere(1.0).is_kinematic().build();
        assert_eq!(body.kind, BodyKind::Kinematic);
    }

    #[test]
    fn test_rigid_body_builder_dynamic_default() {
        let body = RigidBodyBuilder::sphere(1.0).build();
        assert_eq!(body.kind, BodyKind::Dynamic);
    }

    #[test]
    fn test_rigid_body_builder_restitution() {
        let body = RigidBodyBuilder::sphere(1.0).with_restitution(0.9).build();
        assert!((body.restitution - 0.9).abs() < 1e-12);
    }

    #[test]
    fn test_rigid_body_builder_friction() {
        let body = RigidBodyBuilder::sphere(1.0).with_friction(1.2).build();
        assert!((body.friction - 1.2).abs() < 1e-12);
    }

    #[test]
    fn test_rigid_body_builder_damping() {
        let body = RigidBodyBuilder::sphere(1.0).with_damping(0.1, 0.05).build();
        assert!((body.linear_damping - 0.1).abs() < 1e-12);
        assert!((body.angular_damping - 0.05).abs() < 1e-12);
    }

    #[test]
    fn test_rigid_body_builder_damping_clamped() {
        let body = RigidBodyBuilder::sphere(1.0)
            .with_damping(-1.0, -2.0)
            .build();
        assert!((body.linear_damping - 0.0).abs() < 1e-12);
        assert!((body.angular_damping - 0.0).abs() < 1e-12);
    }

    // --- FluidSimBuilder ---

    #[test]
    fn test_fluid_builder_2d_dimensions() {
        let cfg = FluidSimBuilder::lbm_2d(64, 32, 0.01).build();
        assert_eq!(cfg.nx, 64);
        assert_eq!(cfg.ny, 32);
        assert_eq!(cfg.nz, 1);
    }

    #[test]
    fn test_fluid_builder_3d_dimensions() {
        let cfg = FluidSimBuilder::lbm_3d(16, 16, 16, 0.05).build();
        assert_eq!(cfg.nz, 16);
    }

    #[test]
    fn test_fluid_builder_viscosity() {
        let cfg = FluidSimBuilder::lbm_2d(32, 32, 0.01)
            .with_viscosity(1e-3)
            .build();
        assert!((cfg.viscosity - 1e-3).abs() < 1e-15);
    }

    #[test]
    fn test_fluid_builder_inlet_velocity() {
        let cfg = FluidSimBuilder::lbm_2d(32, 32, 0.01)
            .with_inlet_velocity(0.1, 0.0)
            .build();
        assert!((cfg.inlet_velocity[0] - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_fluid_builder_mrt() {
        let cfg = FluidSimBuilder::lbm_2d(32, 32, 0.01).with_mrt().build();
        assert_eq!(cfg.collision_op, LbmCollisionOp::Mrt);
    }

    #[test]
    fn test_fluid_builder_trt() {
        let cfg = FluidSimBuilder::lbm_2d(32, 32, 0.01).with_trt().build();
        assert_eq!(cfg.collision_op, LbmCollisionOp::Trt);
    }

    #[test]
    fn test_fluid_builder_default_bgk() {
        let cfg = FluidSimBuilder::lbm_2d(32, 32, 0.01).build();
        assert_eq!(cfg.collision_op, LbmCollisionOp::Bgk);
    }

    #[test]
    fn test_fluid_builder_dx() {
        let cfg = FluidSimBuilder::lbm_3d(8, 8, 8, 0.025).build();
        assert!((cfg.dx - 0.025).abs() < 1e-12);
    }

    // --- SphSimBuilder ---

    #[test]
    fn test_sph_builder_n_particles() {
        let cfg = SphSimBuilder::new(5000).build();
        assert_eq!(cfg.n_particles, 5000);
    }

    #[test]
    fn test_sph_builder_particle_spacing() {
        let cfg = SphSimBuilder::new(100).with_particle_spacing(0.01).build();
        assert!((cfg.particle_spacing - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_sph_builder_kernel_radius() {
        let cfg = SphSimBuilder::new(100).with_kernel_radius(0.04).build();
        assert!((cfg.kernel_radius - 0.04).abs() < 1e-12);
    }

    #[test]
    fn test_sph_builder_iisph() {
        let cfg = SphSimBuilder::new(100)
            .with_pressure_solver_iisph()
            .build();
        assert_eq!(cfg.pressure_solver, SphPressureSolver::Iisph);
    }

    #[test]
    fn test_sph_builder_pcisph() {
        let cfg = SphSimBuilder::new(100)
            .with_pressure_solver_pcisph()
            .build();
        assert_eq!(cfg.pressure_solver, SphPressureSolver::Pcisph);
    }

    #[test]
    fn test_sph_builder_rest_density() {
        let cfg = SphSimBuilder::new(100).with_rest_density(800.0).build();
        assert!((cfg.rest_density - 800.0).abs() < 1e-12);
    }

    #[test]
    fn test_sph_builder_viscosity() {
        let cfg = SphSimBuilder::new(100).with_viscosity(0.001).build();
        assert!((cfg.viscosity - 0.001).abs() < 1e-12);
    }

    #[test]
    fn test_sph_builder_default_wcsph() {
        let cfg = SphSimBuilder::new(100).build();
        assert_eq!(cfg.pressure_solver, SphPressureSolver::WcSph);
    }

    // --- FemBuilder ---

    #[test]
    fn test_fem_builder_nodes() {
        let cfg = FemBuilder::new()
            .add_node([0.0, 0.0, 0.0])
            .add_node([1.0, 0.0, 0.0])
            .build();
        assert_eq!(cfg.nodes.len(), 2);
    }

    #[test]
    fn test_fem_builder_tet() {
        let cfg = FemBuilder::new()
            .add_node([0.0, 0.0, 0.0])
            .add_node([1.0, 0.0, 0.0])
            .add_node([0.0, 1.0, 0.0])
            .add_node([0.0, 0.0, 1.0])
            .add_tet([0, 1, 2, 3])
            .build();
        assert_eq!(cfg.elements.len(), 1);
        assert!(matches!(cfg.elements[0], FemElement::Tet([0, 1, 2, 3])));
    }

    #[test]
    fn test_fem_builder_hex() {
        let cfg = FemBuilder::new()
            .add_hex([0, 1, 2, 3, 4, 5, 6, 7])
            .build();
        assert!(matches!(cfg.elements[0], FemElement::Hex(_)));
    }

    #[test]
    fn test_fem_builder_material() {
        let cfg = FemBuilder::new()
            .with_material_elastic(200e9, 0.3, 7800.0)
            .build();
        let mat = cfg.material.unwrap();
        assert!((mat.young_modulus - 200e9).abs() < 1.0);
        assert!((mat.poisson_ratio - 0.3).abs() < 1e-12);
        assert!((mat.density - 7800.0).abs() < 1e-12);
    }

    #[test]
    fn test_fem_builder_no_material() {
        let cfg = FemBuilder::new().build();
        assert!(cfg.material.is_none());
    }

    #[test]
    fn test_fem_builder_dirichlet() {
        let cfg = FemBuilder::new()
            .with_dirichlet(0, 1, 0.0)
            .build();
        assert_eq!(cfg.dirichlet_bcs.len(), 1);
        assert_eq!(cfg.dirichlet_bcs[0].node, 0);
        assert_eq!(cfg.dirichlet_bcs[0].dof, 1);
    }

    #[test]
    fn test_fem_builder_multiple_dirichlet() {
        let cfg = FemBuilder::new()
            .with_dirichlet(0, 0, 0.0)
            .with_dirichlet(0, 1, 0.0)
            .with_dirichlet(0, 2, 0.0)
            .build();
        assert_eq!(cfg.dirichlet_bcs.len(), 3);
    }

    #[test]
    fn test_fem_builder_empty() {
        let cfg = FemBuilder::new().build();
        assert!(cfg.nodes.is_empty());
        assert!(cfg.elements.is_empty());
    }

    // --- BroadphaseType ---

    #[test]
    fn test_broadphase_type_default() {
        let bt = BroadphaseType::default();
        assert_eq!(bt, BroadphaseType::SweepAndPrune);
    }

    #[test]
    fn test_broadphase_type_variants() {
        let variants = [
            BroadphaseType::DynamicAabb,
            BroadphaseType::SweepAndPrune,
            BroadphaseType::BruteForce,
        ];
        assert_eq!(variants.len(), 3);
    }

    // --- BodyKind ---

    #[test]
    fn test_body_kind_default() {
        let bk = BodyKind::default();
        assert_eq!(bk, BodyKind::Dynamic);
    }

    // --- MassSpec ---

    #[test]
    fn test_mass_spec_density_default() {
        let ms = MassSpec::default();
        assert!(matches!(ms, MassSpec::Density(d) if (d - 1000.0).abs() < 1e-12));
    }
}
