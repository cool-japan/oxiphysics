// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Granular sand simulation (MLS-MPM with Drucker-Prager plasticity).
//!
//! [`SandSimulation`] is a Material Point Method solver for dry, cohesionless
//! sand: moving-least-squares MPM with quadratic B-spline weights and APIC
//! affine velocities (Hu et al., 2018), Hencky elasticity, and the
//! Drucker-Prager return mapping in logarithmic-strain space of Klár et al.
//! (2016, "Drucker-Prager elastoplasticity for sand animation"). The
//! background grid is sparse (unbounded domain, memory proportional to the
//! number of particles) and a rough floor plane `y = 0` closes the domain.
//! See the private `mpm` submodule for the scheme in formulas.
//!
//! What is modelled: gravity, elastic compression, frictional shear failure
//! (a larger friction angle gives steeper stable heaps and a shorter runout
//! when a column collapses), tension cut-off (grains do not stick), and contact
//! with the floor with Coulomb friction. What is not: cohesion, Klár's volume
//! correction (the dilation of loosened sand is only accumulated in
//! [`SandParticle::plastic_strain_vol`]), pore fluid, and colliders other than
//! the floor plane.
//!
//! The friction angle `φ` calibrates the Drucker-Prager cone to the
//! Mohr-Coulomb criterion in **triaxial compression**, i.e. for the stress
//! states of a 3-D heap. In plane strain (see
//! [`SandSimulation::plane_strain_z`]) that cone is more permissive: a heap of
//! `φ = 30°` sand stands up to roughly 45° there, and for `φ ≳ 42°` the
//! idealised plane-strain stress state (`σ₂ = (σ₁ + σ₃)/2`) never reaches the
//! cone, so such sand does not fail in plane strain.
//!
//! The module also provides avalanche detection, angle-of-repose pile
//! geometry, the stand-alone stress-space [`DruckerPragerPlasticity`] helper
//! (not used by the simulation, which applies the return mapping in
//! logarithmic-strain space), and granular temperature / packing-fraction
//! utilities.

mod mpm;

use crate::{Error, Result};

// ---------------------------------------------------------------------------
// Small math helpers
// ---------------------------------------------------------------------------

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ---------------------------------------------------------------------------
// SandParticle
// ---------------------------------------------------------------------------

/// A single sand (MPM) particle.
///
/// All state the integrator needs lives in this struct (the simulation keeps no
/// hidden per-particle data), so particles may be added, removed, copied or
/// edited between steps. A particle with a non-finite state, a non-positive
/// mass or volume, or a position more than about a million cells from the origin
/// is *inert*: it is never modified and does not influence the others.
#[derive(Debug, Clone, PartialEq)]
pub struct SandParticle {
    /// World-space position \[m\].
    pub position: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Mass \[kg\].
    pub mass: f64,
    /// Reference volume \[m³\]; `mass / volume` is the bulk density of the sand.
    pub volume: f64,
    /// Elastic deformation gradient F^E stored column-major as a flat 3×3
    /// array (identity when stress free). The plastic part is discarded by the
    /// return mapping, so F^E measures the strain the grains are carrying.
    pub deformation_gradient: [f64; 9],
    /// Accumulated volumetric plastic strain: the dilation (`tr ε ≥ 0`) that the
    /// tension cut-off removed from the elastic strain; never decreases.
    pub plastic_strain_vol: f64,
    /// APIC affine velocity matrix C \[1/s\] (the local velocity gradient the
    /// particle carries), column-major; zero for a particle at rest.
    pub affine_velocity: [f64; 9],
}

impl SandParticle {
    /// Create a new `SandParticle` at rest with an identity deformation gradient.
    pub fn new(position: [f64; 3], mass: f64, volume: f64) -> Self {
        // Identity matrix stored column-major
        #[rustfmt::skip]
        let identity = [
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
        ];
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            volume,
            deformation_gradient: identity,
            plastic_strain_vol: 0.0,
            affine_velocity: [0.0; 9],
        }
    }

    /// Kinetic energy of this particle \[J\].
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.velocity, self.velocity)
    }
}

// ---------------------------------------------------------------------------
// DruckerPragerPlasticity
// ---------------------------------------------------------------------------

/// Drucker-Prager (cone) plasticity model for granular materials.
///
/// The yield criterion is:
/// `‖s‖ ≤ k_phi · p` where `p = -σ_v/3` is the mean pressure and
/// `k_phi = 2·sin(φ)/(√3·(3−sin(φ)))`.
pub struct DruckerPragerPlasticity {
    /// Friction angle φ \[radians\].
    pub friction_angle: f64,
    /// Cohesion c \[Pa\].
    pub cohesion: f64,
}

impl DruckerPragerPlasticity {
    /// Create a new `DruckerPragerPlasticity`.
    pub fn new(friction_angle: f64, cohesion: f64) -> Self {
        Self {
            friction_angle,
            cohesion,
        }
    }

    /// Drucker-Prager cone constant k_φ.
    pub fn k_phi(&self) -> f64 {
        let sin_phi = self.friction_angle.sin();
        2.0 * sin_phi / (3.0_f64.sqrt() * (3.0 - sin_phi))
    }

    /// Return `true` if the stress state (σ_v, σ_d) is *outside* the yield surface.
    ///
    /// # Arguments
    /// * `sigma_v` – Volumetric (trace) stress \[Pa\]: σ_v = σ_xx + σ_yy + σ_zz.
    /// * `sigma_d` – Deviatoric stress norm ‖s‖ \[Pa\].
    pub fn yield_criterion(&self, sigma_v: f64, sigma_d: f64) -> bool {
        // Mean pressure p = -σ_v / 3  (compression positive)
        let p = -sigma_v / 3.0;
        let limit = self.k_phi() * p + self.cohesion;
        sigma_d > limit
    }

    /// Radial return mapping: project a 6-component Voigt stress back onto the
    /// yield cone.
    ///
    /// # Arguments
    /// * `stress` – Voigt stress \[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz\] \[Pa\].
    ///
    /// Returns the projected stress in the same layout.
    pub fn plastic_return_mapping(&self, stress: [f64; 6]) -> [f64; 6] {
        let sigma_v = stress[0] + stress[1] + stress[2];
        // Deviatoric part
        let p = sigma_v / 3.0;
        let s = [
            stress[0] - p,
            stress[1] - p,
            stress[2] - p,
            stress[3],
            stress[4],
            stress[5],
        ];
        let s_norm = (s[0] * s[0]
            + s[1] * s[1]
            + s[2] * s[2]
            + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
            .sqrt();

        // If inside yield surface, return unchanged.
        if !self.yield_criterion(sigma_v, s_norm) {
            return stress;
        }

        // Compute permissible deviatoric norm.
        let mean_p = -sigma_v / 3.0;
        let limit = (self.k_phi() * mean_p + self.cohesion).max(0.0);

        if s_norm < 1e-30 {
            return stress;
        }

        let scale = limit / s_norm;
        [
            scale * s[0] + p,
            scale * s[1] + p,
            scale * s[2] + p,
            scale * s[3],
            scale * s[4],
            scale * s[5],
        ]
    }
}

// ---------------------------------------------------------------------------
// SandSimulation
// ---------------------------------------------------------------------------

/// Young's modulus of a new simulation \[Pa\] (the sand of Klár et al. 2016).
const DEFAULT_YOUNGS_MODULUS_PA: f64 = 3.537e5;
/// Poisson's ratio of a new simulation.
const DEFAULT_POISSON_RATIO: f64 = 0.3;
/// Friction angle of a new simulation \[degrees\].
const DEFAULT_FRICTION_ANGLE_DEG: f64 = 30.0;
/// Floor friction coefficient of a new simulation.
const DEFAULT_FLOOR_FRICTION: f64 = 1.0;
/// CFL number of a new simulation.
const DEFAULT_CFL: f64 = 0.3;
/// Largest particle count [`SandSimulation::fill_box`] creates in one call.
const MAX_FILL_PARTICLES: usize = 50_000_000;

/// Upper bound on the number of explicit sub-steps that one
/// [`SandSimulation::try_step`] call takes before it gives up.
pub const MAX_SUBSTEPS_PER_STEP: usize = 20_000;

/// MLS-MPM sand simulation: APIC particle-grid transfers, Hencky elasticity and
/// a Drucker-Prager return mapping (see the module documentation).
///
/// The grid is sparse and built from the particles every sub-step, so there is
/// no domain to size; the floor plane `y = 0` is the only boundary. A good
/// resolution has about two particles per cell edge: `grid_size ≈ 2 · spacing`.
#[derive(Debug, Clone)]
pub struct SandSimulation {
    /// Particle list.
    pub particles: Vec<SandParticle>,
    /// Background grid cell size \[m\].
    pub grid_size: f64,
    /// Gravity vector \[m/s²\].
    pub gravity: [f64; 3],
    /// Young's modulus E of the elastic response \[Pa\] (default `3.537e5`).
    pub youngs_modulus_pa: f64,
    /// Poisson's ratio ν, `0 ≤ ν < 0.5` (default `0.3`).
    pub poisson_ratio: f64,
    /// Friction angle φ of the sand in degrees, `0 ≤ φ < 90` (default `30`).
    ///
    /// The Drucker-Prager cone is calibrated to the Mohr-Coulomb friction angle
    /// in triaxial compression, `α = √(2/3) · 2 sin φ / (3 − sin φ)`.
    pub friction_angle_deg: f64,
    /// Coulomb friction coefficient µ ≥ 0 between the sand and the floor plane
    /// (default `1`); `f64::INFINITY` makes the floor sticky.
    pub floor_friction: f64,
    /// CFL number of the sub-step control, `0 < cfl ≤ 1` (default `0.3`).
    pub cfl: f64,
    /// Optional plane-strain slab `[z_min, z_max]` in metres (default `None`:
    /// a full 3-D simulation).
    ///
    /// The planes `z = z_min` and `z = z_max` are symmetry planes: grid nodes on
    /// or beyond them get zero `z` velocity and particles are kept inside. A slab
    /// one or two cells thick (`z_max − z_min ≈ grid_size`) filled symmetrically
    /// gives an exact two-dimensional (plane-strain) simulation at a fraction of
    /// the 3-D cost; align the planes with grid planes (multiples of
    /// `grid_size`).
    pub plane_strain_z: Option<[f64; 2]>,
}

impl SandSimulation {
    /// Create a new `SandSimulation` with the default sand (`E = 3.537e5 Pa`,
    /// `ν = 0.3`, `φ = 30°`), a rough floor (`µ = 1`) and `cfl = 0.3`.
    pub fn new(grid_size: f64, gravity: [f64; 3]) -> Self {
        Self {
            particles: Vec::new(),
            grid_size,
            gravity,
            youngs_modulus_pa: DEFAULT_YOUNGS_MODULUS_PA,
            poisson_ratio: DEFAULT_POISSON_RATIO,
            friction_angle_deg: DEFAULT_FRICTION_ANGLE_DEG,
            floor_friction: DEFAULT_FLOOR_FRICTION,
            cfl: DEFAULT_CFL,
            plane_strain_z: None,
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: SandParticle) {
        self.particles.push(p);
    }

    /// Fill the axis-aligned box `[min, max]` with particles on a regular
    /// lattice of cell centres and return how many were added.
    ///
    /// Each axis is divided into `round((max − min) / spacing)` equal cells (so
    /// the lattice fills the box exactly); every particle owns one cell, with
    /// `volume = cell volume` and `mass = density · volume`. An axis shorter than
    /// half a `spacing` yields no particles.
    ///
    /// # Errors
    /// [`Error::InvalidParameter`] if a bound is not finite, `min ≥ max` on an
    /// axis, `spacing` or `density` is not positive and finite, or the box would
    /// hold more than 50 million particles. Nothing is added then.
    pub fn fill_box(
        &mut self,
        min: [f64; 3],
        max: [f64; 3],
        spacing: f64,
        density: f64,
    ) -> Result<usize> {
        let bounds_ok = min.iter().chain(max.iter()).all(|b| b.is_finite())
            && min.iter().zip(&max).all(|(lo, hi)| lo < hi);
        if !bounds_ok {
            return Err(Error::InvalidParameter(format!(
                "fill_box needs finite bounds with min < max, got {min:?} .. {max:?}"
            )));
        }
        if !(spacing.is_finite() && spacing > 0.0 && density.is_finite() && density > 0.0) {
            return Err(Error::InvalidParameter(format!(
                "fill_box needs positive finite spacing and density, got {spacing} and {density}"
            )));
        }
        let counts: [usize; 3] =
            std::array::from_fn(|d| ((max[d] - min[d]) / spacing).round() as usize);
        let total = counts[0]
            .checked_mul(counts[1])
            .and_then(|n| n.checked_mul(counts[2]))
            .filter(|&n| n <= MAX_FILL_PARTICLES)
            .ok_or_else(|| {
                Error::InvalidParameter(format!(
                    "fill_box would create more than {MAX_FILL_PARTICLES} particles"
                ))
            })?;
        if total == 0 {
            return Ok(0);
        }
        let cell: [f64; 3] = std::array::from_fn(|d| (max[d] - min[d]) / counts[d] as f64);
        let volume = cell[0] * cell[1] * cell[2];
        let mass = density * volume;
        self.particles.reserve(total);
        for i in 0..counts[0] {
            for j in 0..counts[1] {
                for k in 0..counts[2] {
                    let position = [
                        min[0] + (i as f64 + 0.5) * cell[0],
                        min[1] + (j as f64 + 0.5) * cell[1],
                        min[2] + (k as f64 + 0.5) * cell[2],
                    ];
                    self.particles
                        .push(SandParticle::new(position, mass, volume));
                }
            }
        }
        Ok(total)
    }

    /// Validate the parameters and derive the constitutive constants.
    fn material(&self) -> Result<mpm::Material> {
        let h = self.grid_size;
        if !(h.is_finite() && h > 0.0 && (1.0 / h).is_finite()) {
            return Err(Error::InvalidParameter(format!(
                "grid_size must be positive and finite, got {h}"
            )));
        }
        if !self.gravity.iter().all(|g| g.is_finite()) {
            return Err(Error::InvalidParameter(format!(
                "gravity must be finite, got {:?}",
                self.gravity
            )));
        }
        if self.floor_friction.is_nan() || self.floor_friction < 0.0 {
            return Err(Error::InvalidParameter(format!(
                "floor_friction must be non-negative, got {}",
                self.floor_friction
            )));
        }
        if let Some([z_min, z_max]) = self.plane_strain_z
            && !(z_min.is_finite() && z_max.is_finite() && z_min < z_max)
        {
            return Err(Error::InvalidParameter(format!(
                "plane_strain_z needs finite z_min < z_max, got [{z_min}, {z_max}]"
            )));
        }
        if !(self.cfl.is_finite() && self.cfl > 0.0 && self.cfl <= 1.0) {
            return Err(Error::InvalidParameter(format!(
                "cfl must be in (0, 1], got {}",
                self.cfl
            )));
        }
        mpm::Material::new(self.youngs_modulus_pa, self.poisson_ratio, self.friction_angle_deg).ok_or_else(|| {
            Error::InvalidParameter(format!(
                "sand material out of range: youngs_modulus_pa = {} (> 0), poisson_ratio = {} (0 <= nu < 0.5), \
                 friction_angle_deg = {} (0 <= phi < 90)",
                self.youngs_modulus_pa, self.poisson_ratio, self.friction_angle_deg
            ))
        })
    }

    /// Largest stable explicit time step \[s\] of the current state:
    /// `cfl · h / (c + v_max)` with the P-wave speed `c = √((λ + 2µ) / ρ_min)` of
    /// the least dense particle and the largest particle speed `v_max`.
    /// `f64::INFINITY` if there is no active particle.
    ///
    /// # Errors
    /// [`Error::InvalidParameter`] if the grid size, gravity, material
    /// constants, floor friction or CFL number are invalid.
    pub fn stable_time_step(&self) -> Result<f64> {
        let material = self.material()?;
        Ok(mpm::stable_time_step(
            &self.particles,
            self.grid_size,
            self.cfl,
            &material,
        ))
    }

    /// Advance the simulation by `dt` \[s\], taking as many explicit MLS-MPM
    /// sub-steps as the stability limit [`SandSimulation::stable_time_step`]
    /// requires (recomputed from the current velocities before every
    /// sub-step), so any `dt` is safe.
    ///
    /// Particles that were placed below the floor (or outside the plane-strain
    /// slab) are first moved onto it. Inert particles (see [`SandParticle`]) are skipped.
    ///
    /// # Errors
    /// [`Error::InvalidParameter`] for an invalid `dt` (negative or not finite)
    /// or simulation parameter, in which case nothing changes; also if more
    /// than [`MAX_SUBSTEPS_PER_STEP`] sub-steps would be needed, in which case
    /// the sub-steps already taken stay applied.
    pub fn try_step(&mut self, dt: f64) -> Result<()> {
        let material = self.material()?;
        if !dt.is_finite() || dt < 0.0 {
            return Err(Error::InvalidParameter(format!(
                "time step must be finite and non-negative, got {dt}"
            )));
        }
        if dt == 0.0 {
            return Ok(());
        }
        let h = self.grid_size;
        let inv_h = 1.0 / h;
        for p in self
            .particles
            .iter_mut()
            .filter(|p| mpm::is_active(p, inv_h))
        {
            mpm::confine_to_domain(p, self.plane_strain_z);
        }
        let mut constants = mpm::StepConstants {
            h,
            dt,
            gravity: self.gravity,
            material,
            floor_friction: self.floor_friction,
            plane_strain_z: self.plane_strain_z,
        };
        let mut remaining = dt;
        let mut substeps = 0usize;
        while remaining > 0.0 {
            if substeps >= MAX_SUBSTEPS_PER_STEP {
                return Err(Error::InvalidParameter(format!(
                    "time step {dt} s needs more than {MAX_SUBSTEPS_PER_STEP} stability-limited sub-steps"
                )));
            }
            let stable = mpm::stable_time_step(&self.particles, h, self.cfl, &material);
            if stable.is_nan() || stable <= 0.0 {
                return Err(Error::InvalidParameter(format!(
                    "no stable time step exists for the current state ({stable})"
                )));
            }
            let pieces = (remaining / stable).ceil().max(1.0);
            constants.dt = remaining / pieces;
            mpm::substep(&mut self.particles, &constants);
            remaining -= constants.dt;
            substeps += 1;
        }
        Ok(())
    }

    /// Advance the simulation by `dt` \[s\] (see [`SandSimulation::try_step`]).
    ///
    /// This is the infallible convenience form: an invalid `dt` or invalid
    /// simulation parameters leave the particles unchanged. Use
    /// [`SandSimulation::try_step`] to learn why a step was refused.
    pub fn step(&mut self, dt: f64) {
        // The error carries no state the caller could use beyond what `try_step` reports.
        let _ = self.try_step(dt);
    }

    /// Total mass of all particles \[kg\].
    pub fn total_mass(&self) -> f64 {
        self.particles.iter().map(|p| p.mass).sum()
    }

    /// Total linear momentum `Σ m v` \[kg·m/s\].
    pub fn total_momentum(&self) -> [f64; 3] {
        self.particles
            .iter()
            .fold([0.0; 3], |acc, p| add3(acc, scale3(p.velocity, p.mass)))
    }

    /// Total kinetic energy `Σ ½ m |v|²` \[J\] (of the particle velocities).
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(SandParticle::kinetic_energy)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// AvalancheDetector
// ---------------------------------------------------------------------------

/// Detects sand particles that are at risk of avalanching.
///
/// A particle is considered unstable when the local surface slope exceeds
/// `critical_angle_deg` degrees relative to horizontal.
pub struct AvalancheDetector {
    /// Critical angle of repose \[degrees\] beyond which a particle is unstable.
    pub critical_angle_deg: f64,
    /// Reference to the particle list (indices into an external Vec).
    particles: Vec<[f64; 3]>,
}

impl AvalancheDetector {
    /// Create a new `AvalancheDetector` from a list of particle positions.
    pub fn new(critical_angle_deg: f64, positions: Vec<[f64; 3]>) -> Self {
        Self {
            critical_angle_deg,
            particles: positions,
        }
    }

    /// Return the indices of particles that are unstable.
    ///
    /// A particle at index `i` is considered unstable if its height above the
    /// lowest neighbour (within `search_radius`) corresponds to a slope angle
    /// greater than `critical_angle_deg`.  When no neighbours are found the
    /// particle is considered stable.
    pub fn detect_unstable(&self) -> Vec<usize> {
        let crit_rad = self.critical_angle_deg.to_radians();
        let tan_crit = crit_rad.tan();
        let search_radius = 2.0; // look within 2 m (normalised)

        let mut unstable = Vec::new();
        for (i, &pi) in self.particles.iter().enumerate() {
            let mut most_steep: f64 = 0.0;
            for (j, &pj) in self.particles.iter().enumerate() {
                if i == j {
                    continue;
                }
                let dx = pi[0] - pj[0];
                let dz = pi[2] - pj[2];
                let horiz = (dx * dx + dz * dz).sqrt();
                if horiz > search_radius || horiz < 1e-12 {
                    continue;
                }
                let dy = pi[1] - pj[1];
                if dy > 0.0 {
                    let slope = dy / horiz;
                    if slope > most_steep {
                        most_steep = slope;
                    }
                }
            }
            if most_steep > tan_crit {
                unstable.push(i);
            }
        }
        unstable
    }
}

// ---------------------------------------------------------------------------
// SandPile
// ---------------------------------------------------------------------------

/// Idealised conical sand pile.
pub struct SandPile {
    /// Base radius (half base-width) \[m\].
    pub base_width: f64,
    /// Angle of repose \[degrees\].
    pub angle_of_repose: f64,
}

impl SandPile {
    /// Create a new `SandPile`.
    pub fn new(base_width: f64, angle_of_repose: f64) -> Self {
        Self {
            base_width,
            angle_of_repose,
        }
    }

    /// Maximum height of the conical pile: h = r · tan(θ) \[m\].
    pub fn max_height(&self) -> f64 {
        let r = self.base_width / 2.0;
        r * self.angle_of_repose.to_radians().tan()
    }

    /// Surface height of the pile at radial distance `x` from the apex \[m\].
    ///
    /// Returns 0.0 for `x ≥ base_width / 2`.
    pub fn profile(&self, x: f64) -> f64 {
        let r = self.base_width / 2.0;
        if x.abs() >= r {
            return 0.0;
        }
        self.max_height() * (1.0 - x.abs() / r)
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Convert a Coulomb friction coefficient µ to angle of repose \[radians\].
///
/// θ = atan(µ)
pub fn angle_of_repose_from_friction(mu: f64) -> f64 {
    mu.atan()
}

/// Granular temperature of a particle assembly \[m²/s²\].
///
/// T_g = (1 / N) · Σ (v_i − ⟨v⟩)² where the sum is over all particles.
pub fn granular_temperature(particles: &[SandParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let n = particles.len() as f64;
    // Mean velocity
    let mut mean = [0.0_f64; 3];
    for p in particles {
        mean = add3(mean, p.velocity);
    }
    mean = scale3(mean, 1.0 / n);
    // Variance
    let mut var = 0.0_f64;
    for p in particles {
        let dv = sub3(p.velocity, mean);
        var += dot3(dv, dv);
    }
    var / n
}

/// Random close-packing fraction for monodisperse spheres of diameter `d`.
///
/// Returns φ ≈ 0.64 (Random Close Pack) — independent of particle diameter
/// for uniform spheres.  The `d_particle` parameter is accepted for API
/// completeness (it cancels out for uniform spheres).
pub fn packing_fraction_random(_d_particle: f64) -> f64 {
    // Empirical random close-packing fraction for monodisperse spheres.
    // The measured φ_RCP is conventionally reported as 0.6366; this is numerically
    // equivalent to 2/π and we use the constant to make the intent explicit.
    std::f64::consts::FRAC_2_PI
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // --- SandParticle ---

    #[test]
    fn test_particle_identity_deformation() {
        let p = SandParticle::new([0.0, 0.0, 0.0], 1.0, 1e-6);
        // Diagonal of identity matrix (column-major indices 0,4,8)
        assert!((p.deformation_gradient[0] - 1.0).abs() < EPS);
        assert!((p.deformation_gradient[4] - 1.0).abs() < EPS);
        assert!((p.deformation_gradient[8] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_particle_kinetic_energy_at_rest() {
        let p = SandParticle::new([0.0, 0.0, 0.0], 2.0, 1e-6);
        assert!(p.kinetic_energy().abs() < EPS);
    }

    #[test]
    fn test_particle_kinetic_energy_moving() {
        let mut p = SandParticle::new([0.0, 0.0, 0.0], 2.0, 1e-6);
        p.velocity = [3.0, 4.0, 0.0]; // speed = 5 m/s
        assert!((p.kinetic_energy() - 0.5 * 2.0 * 25.0).abs() < EPS);
    }

    #[test]
    fn test_particle_zero_plastic_strain() {
        let p = SandParticle::new([1.0, 2.0, 3.0], 0.5, 5e-7);
        assert!(p.plastic_strain_vol.abs() < EPS);
    }

    // --- DruckerPragerPlasticity ---

    #[test]
    fn test_yield_criterion_inside_cone() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 0.0);
        // Compressive pressure p = 1e4 Pa → limit = k_phi * 1e4
        // Use a deviatoric norm well below limit.
        let sigma_v = -3e4; // volumetric stress → p = 1e4
        let limit = dp.k_phi() * 1e4;
        let sigma_d = 0.5 * limit;
        assert!(
            !dp.yield_criterion(sigma_v, sigma_d),
            "Should be inside yield surface"
        );
    }

    #[test]
    fn test_yield_criterion_outside_cone() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 0.0);
        let sigma_v = -3e4;
        let limit = dp.k_phi() * 1e4;
        let sigma_d = 2.0 * limit; // outside
        assert!(
            dp.yield_criterion(sigma_v, sigma_d),
            "Should be outside yield surface"
        );
    }

    #[test]
    fn test_yield_criterion_zero_deviator() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 0.0);
        assert!(!dp.yield_criterion(-3e4, 0.0));
    }

    #[test]
    fn test_return_mapping_inside_unchanged() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 1000.0);
        let stress = [-10000.0, -10000.0, -10000.0, 0.0, 0.0, 0.0];
        let projected = dp.plastic_return_mapping(stress);
        for i in 0..6 {
            assert!((projected[i] - stress[i]).abs() < EPS);
        }
    }

    #[test]
    fn test_return_mapping_projects_onto_cone() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 0.0);
        // Large deviatoric shear stress (outside cone)
        let stress = [-30000.0, -30000.0, -30000.0, 50000.0, 0.0, 0.0];
        let projected = dp.plastic_return_mapping(stress);
        // After projection the yield criterion should be satisfied (not yielding)
        let sv = projected[0] + projected[1] + projected[2];
        let p = sv / 3.0;
        let s = [
            projected[0] - p,
            projected[1] - p,
            projected[2] - p,
            projected[3],
            projected[4],
            projected[5],
        ];
        let s_norm = (s[0] * s[0]
            + s[1] * s[1]
            + s[2] * s[2]
            + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
            .sqrt();
        assert!(
            !dp.yield_criterion(sv, s_norm),
            "After projection, should be on/inside cone"
        );
    }

    #[test]
    fn test_k_phi_increases_with_friction_angle() {
        let dp1 = DruckerPragerPlasticity::new(20_f64.to_radians(), 0.0);
        let dp2 = DruckerPragerPlasticity::new(40_f64.to_radians(), 0.0);
        assert!(dp2.k_phi() > dp1.k_phi());
    }

    #[test]
    fn test_yield_criterion_with_cohesion() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 5000.0);
        // Even with zero pressure, cohesion provides a limit.
        assert!(!dp.yield_criterion(0.0, 4999.0));
        assert!(dp.yield_criterion(0.0, 5001.0));
    }

    // --- SandSimulation ---

    #[test]
    fn test_simulation_gravity_falls() {
        let mut sim = SandSimulation::new(0.1, [0.0, -9.81, 0.0]);
        sim.add_particle(SandParticle::new([0.0, 5.0, 0.0], 1.0, 1e-6));
        let dt = 0.01;
        for _ in 0..50 {
            sim.step(dt);
        }
        assert!(
            sim.particles[0].position[1] < 5.0,
            "Particle should fall under gravity"
        );
    }

    #[test]
    fn test_simulation_ground_plane() {
        let mut sim = SandSimulation::new(0.1, [0.0, -9.81, 0.0]);
        sim.add_particle(SandParticle::new([0.0, 0.01, 0.0], 1.0, 1e-6));
        for _ in 0..1000 {
            sim.step(0.01);
        }
        assert!(
            sim.particles[0].position[1] >= 0.0,
            "Particle should not fall below ground"
        );
    }

    #[test]
    fn test_simulation_velocity_updates() {
        let mut sim = SandSimulation::new(0.1, [0.0, -9.81, 0.0]);
        sim.add_particle(SandParticle::new([0.0, 10.0, 0.0], 1.0, 1e-6));
        sim.step(0.1);
        // After one step velocity[1] should be negative (downward)
        assert!(sim.particles[0].velocity[1] < 0.0);
    }

    #[test]
    fn test_simulation_no_particles_runs() {
        let mut sim = SandSimulation::new(0.1, [0.0, -9.81, 0.0]);
        sim.step(0.01); // should not panic
    }

    // --- AvalancheDetector ---

    #[test]
    fn test_avalanche_steep_particle_unstable() {
        // A high particle directly above a low one: angle ≈ 45° → unstable if crit < 45°
        let positions = vec![
            [0.0_f64, 1.0, 0.0], // particle 0: high
            [1.0_f64, 0.0, 0.0], // particle 1: low
        ];
        let detector = AvalancheDetector::new(30.0, positions);
        let unstable = detector.detect_unstable();
        assert!(unstable.contains(&0), "Steep particle 0 should be unstable");
    }

    #[test]
    fn test_avalanche_flat_surface_stable() {
        // All particles at the same height.
        let positions = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0_f64, 0.0, 0.0],
            [2.0_f64, 0.0, 0.0],
        ];
        let detector = AvalancheDetector::new(30.0, positions);
        let unstable = detector.detect_unstable();
        assert!(
            unstable.is_empty(),
            "Flat surface should have no unstable particles"
        );
    }

    #[test]
    fn test_avalanche_single_particle_stable() {
        let positions = vec![[0.0_f64, 0.0, 0.0]];
        let detector = AvalancheDetector::new(30.0, positions);
        assert!(detector.detect_unstable().is_empty());
    }

    #[test]
    fn test_avalanche_empty_stable() {
        let detector = AvalancheDetector::new(30.0, vec![]);
        assert!(detector.detect_unstable().is_empty());
    }

    // --- SandPile ---

    #[test]
    fn test_sand_pile_max_height() {
        // base_width = 2 m, angle = 45° → r = 1, h = tan(45°) = 1 m
        let pile = SandPile::new(2.0, 45.0);
        assert!(
            (pile.max_height() - 1.0).abs() < 1e-10,
            "h={}",
            pile.max_height()
        );
    }

    #[test]
    fn test_sand_pile_profile_apex_zero() {
        let pile = SandPile::new(2.0, 45.0);
        // At x = 0 (apex): height = max_height()
        assert!((pile.profile(0.0) - pile.max_height()).abs() < EPS);
    }

    #[test]
    fn test_sand_pile_profile_base_zero() {
        let pile = SandPile::new(2.0, 45.0);
        // At x = base_radius = 1.0: height = 0
        assert!(pile.profile(1.0).abs() < EPS);
    }

    #[test]
    fn test_sand_pile_profile_outside_zero() {
        let pile = SandPile::new(2.0, 45.0);
        assert!(pile.profile(1.5).abs() < EPS);
    }

    #[test]
    fn test_sand_pile_profile_midpoint() {
        let pile = SandPile::new(2.0, 45.0);
        // At x = 0.5 (half of radius 1.0): height = max_height() / 2
        assert!((pile.profile(0.5) - pile.max_height() * 0.5).abs() < EPS);
    }

    #[test]
    fn test_sand_pile_narrow_tall() {
        // 30° angle, 1 m wide: r = 0.5, h = 0.5 * tan(30°)
        let pile = SandPile::new(1.0, 30.0);
        let expected = 0.5 * 30_f64.to_radians().tan();
        assert!((pile.max_height() - expected).abs() < 1e-10);
    }

    // --- Free functions ---

    #[test]
    fn test_angle_of_repose_from_friction_45deg() {
        // mu = 1 → atan(1) = 45°
        let angle = angle_of_repose_from_friction(1.0);
        assert!((angle - std::f64::consts::FRAC_PI_4).abs() < 1e-10);
    }

    #[test]
    fn test_angle_of_repose_increases_with_mu() {
        let a1 = angle_of_repose_from_friction(0.5);
        let a2 = angle_of_repose_from_friction(1.0);
        assert!(a2 > a1);
    }

    #[test]
    fn test_angle_of_repose_zero_friction() {
        assert!(angle_of_repose_from_friction(0.0).abs() < EPS);
    }

    #[test]
    fn test_granular_temperature_at_rest() {
        let particles = vec![
            SandParticle::new([0.0, 0.0, 0.0], 1.0, 1e-6),
            SandParticle::new([1.0, 0.0, 0.0], 1.0, 1e-6),
        ];
        assert!(granular_temperature(&particles).abs() < EPS);
    }

    #[test]
    fn test_granular_temperature_uniform_flow() {
        // All particles with same velocity → zero fluctuation temperature.
        let mut particles = vec![
            SandParticle::new([0.0, 0.0, 0.0], 1.0, 1e-6),
            SandParticle::new([1.0, 0.0, 0.0], 1.0, 1e-6),
        ];
        for p in &mut particles {
            p.velocity = [2.0, 0.0, 0.0];
        }
        assert!(granular_temperature(&particles).abs() < EPS);
    }

    #[test]
    fn test_granular_temperature_opposite_velocities() {
        let mut particles = vec![
            SandParticle::new([0.0, 0.0, 0.0], 1.0, 1e-6),
            SandParticle::new([1.0, 0.0, 0.0], 1.0, 1e-6),
        ];
        particles[0].velocity = [1.0, 0.0, 0.0];
        particles[1].velocity = [-1.0, 0.0, 0.0];
        // Mean v = 0, T_g = (1² + 1²) / 2 = 1.0
        let tg = granular_temperature(&particles);
        assert!((tg - 1.0).abs() < EPS, "tg={tg}");
    }

    #[test]
    fn test_granular_temperature_empty() {
        assert!(granular_temperature(&[]).abs() < EPS);
    }

    #[test]
    fn test_packing_fraction_random_value() {
        let phi = packing_fraction_random(0.001);
        // Random close packing ≈ 0.64
        assert!(phi > 0.60 && phi < 0.70, "phi={phi}");
    }

    #[test]
    fn test_packing_fraction_independent_of_diameter() {
        // For monodisperse spheres the packing fraction is size-independent.
        let phi1 = packing_fraction_random(0.001);
        let phi2 = packing_fraction_random(1.0);
        assert!((phi1 - phi2).abs() < EPS);
    }
}

#[cfg(test)]
mod mpm_tests;
#[cfg(test)]
mod physics_tests;
#[cfg(test)]
mod sim_tests;
#[cfg(test)]
mod test_support;
