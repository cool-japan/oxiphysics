//! Gear/mimic constraint coupling two revolute joints by a fixed ratio.
//!
//! The constraint couples bodies A→B (joint 1, axis `axis1`) with bodies C→D
//! (joint 2, axis `axis2`) by `ratio r`:
//!
//! ```text
//! velocity constraint:  (ω_B − ω_A)·a₁ − r·(ω_D − ω_C)·a₂ = 0
//! position drift:       C = θ₁ − r·θ₂
//! ```
//!
//! Backlash introduces a unilateral dead zone of half-width `backlash` (radians).
//! When `|C| < backlash` the constraint is disengaged and no impulse is applied.

use core::fmt;

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

#[inline(always)]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline(always)]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline(always)]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline(always)]
fn diag_mul(inv_i: &[f64; 3], v: [f64; 3]) -> [f64; 3] {
    [inv_i[0] * v[0], inv_i[1] * v[1], inv_i[2] * v[2]]
}

// ---------------------------------------------------------------------------
// BodyState
// ---------------------------------------------------------------------------

/// Mutable angular-velocity state and diagonal inverse inertia for one rigid body.
///
/// Used to group the per-body arguments of [`GearConstraint::solve_iteration`]
/// so that the method stays within a readable parameter count.
#[derive(Debug)]
pub struct BodyState<'a> {
    /// World-space angular velocity of the body (rad/s), mutated in-place.
    pub omega: &'a mut [f64; 3],
    /// Diagonal entries of the inverse inertia tensor `I⁻¹` (kg⁻¹·m⁻²).
    /// Set all entries to `0.0` for a kinematically fixed body.
    pub inv_inertia: &'a [f64; 3],
}

impl<'a> BodyState<'a> {
    /// Construct from mutable omega reference and immutable inverse-inertia slice.
    pub fn new(omega: &'a mut [f64; 3], inv_inertia: &'a [f64; 3]) -> Self {
        Self { omega, inv_inertia }
    }
}

// ---------------------------------------------------------------------------
// GearConstraint
// ---------------------------------------------------------------------------

/// Gear/mimic dynamic constraint coupling two revolute joints by a fixed ratio
/// with optional backlash dead zone and Baumgarte position stabilization.
///
/// # Bodies
///
/// - Joint 1: bodies **A** (parent) and **B** (child), rotation axis `axis1`
/// - Joint 2: bodies **C** (parent) and **D** (child), rotation axis `axis2`
///
/// Parent bodies (A, C) are typically the housing/ground; their inverse
/// inertias should be `[0.0; 3]` if they are kinematically fixed.
///
/// # Velocity constraint
///
/// ```text
/// J·ω = 0   where   J = [−a₁, +a₁, +r·a₂, −r·a₂]
/// ```
///
/// acting on `[ω_A, ω_B, ω_C, ω_D]`.
///
/// # Position stabilization
///
/// Baumgarte bias `β/dt · (θ₁ − r·θ₂)` is added to the velocity residual to
/// prevent position drift from accumulating over many simulation steps.
///
/// # Backlash
///
/// When `|θ₁ − r·θ₂| < backlash`, the gears are disengaged and no impulse is
/// applied. Once contact is reestablished on either boundary the constraint
/// acts as a unilateral contact, preventing further penetration.
#[derive(Debug, Clone, PartialEq)]
pub struct GearConstraint {
    /// Gear ratio r: one revolution of joint 1 drives r revolutions of joint 2.
    pub ratio: f64,
    /// Half-width of the backlash dead zone in radians. `0.0` gives a rigid
    /// constraint with no dead zone.
    pub backlash: f64,
    /// World-space unit axis of joint 1.
    pub axis1: [f64; 3],
    /// World-space unit axis of joint 2.
    pub axis2: [f64; 3],
    /// Baumgarte position stabilization coefficient (dimensionless, default `0.2`).
    pub beta: f64,
}

impl GearConstraint {
    /// Construct a new `GearConstraint`.
    ///
    /// `beta` is defaulted to `0.2`. Use direct field access to override.
    pub fn new(ratio: f64, backlash: f64, axis1: [f64; 3], axis2: [f64; 3]) -> Self {
        Self {
            ratio,
            backlash,
            axis1,
            axis2,
            beta: 0.2,
        }
    }

    /// Solve one Temporal Gauss–Seidel (TGS) iteration.
    ///
    /// Computes and applies an angular-velocity impulse λ that enforces the
    /// gear velocity constraint with Baumgarte position stabilization.
    ///
    /// # Arguments
    ///
    /// - `body_a`, `body_b` — parent/child of joint 1 (axis `axis1`)
    /// - `body_c`, `body_d` — parent/child of joint 2 (axis `axis2`)
    /// - `theta1`           — current angle of joint 1 (rad)
    /// - `theta2`           — current angle of joint 2 (rad)
    /// - `dt`               — time step (s); must be positive
    ///
    /// Each [`BodyState`] bundles the mutable angular velocity with the
    /// (immutable) diagonal inverse inertia for that body.
    ///
    /// # Returns
    ///
    /// The impulse magnitude λ applied this iteration, or `0.0` when the
    /// constraint is disengaged (inside the backlash dead zone).
    pub fn solve_iteration(
        &self,
        body_a: &mut BodyState<'_>,
        body_b: &mut BodyState<'_>,
        body_c: &mut BodyState<'_>,
        body_d: &mut BodyState<'_>,
        theta1: f64,
        theta2: f64,
        dt: f64,
    ) -> f64 {
        let r = self.ratio;
        let a1 = self.axis1;
        let a2 = self.axis2;

        // ------------------------------------------------------------------
        // Position drift
        // ------------------------------------------------------------------
        let drift = theta1 - r * theta2;

        // ------------------------------------------------------------------
        // Backlash dead zone — disengage if inside [-backlash, +backlash]
        // ------------------------------------------------------------------
        if self.backlash > 0.0 && drift.abs() < self.backlash {
            return 0.0;
        }

        // ------------------------------------------------------------------
        // Effective mass
        //   M_eff⁻¹ = a₁ᵀ I_A⁻¹ a₁ + a₁ᵀ I_B⁻¹ a₁
        //           + r²·a₂ᵀ I_C⁻¹ a₂ + r²·a₂ᵀ I_D⁻¹ a₂
        // ------------------------------------------------------------------
        let r2 = r * r;
        let inv_eff_mass = dot3(a1, diag_mul(body_a.inv_inertia, a1))
            + dot3(a1, diag_mul(body_b.inv_inertia, a1))
            + r2 * dot3(a2, diag_mul(body_c.inv_inertia, a2))
            + r2 * dot3(a2, diag_mul(body_d.inv_inertia, a2));

        // Guard against fully kinematic systems (all bodies fixed).
        if inv_eff_mass < f64::EPSILON {
            return 0.0;
        }

        let eff_mass = 1.0 / inv_eff_mass;

        // ------------------------------------------------------------------
        // Velocity error: J·ω
        //   vel_err = (ω_B − ω_A)·a₁ − r·(ω_D − ω_C)·a₂
        // ------------------------------------------------------------------
        let vel_err = dot3(sub3(*body_b.omega, *body_a.omega), a1)
            - r * dot3(sub3(*body_d.omega, *body_c.omega), a2);

        // ------------------------------------------------------------------
        // Baumgarte bias
        // ------------------------------------------------------------------
        let bias = (self.beta / dt) * drift;

        // ------------------------------------------------------------------
        // Raw impulse
        // ------------------------------------------------------------------
        let mut lambda = -(vel_err + bias) * eff_mass;

        // ------------------------------------------------------------------
        // Unilateral backlash clamping
        // ------------------------------------------------------------------
        if self.backlash > 0.0 {
            if drift > self.backlash {
                // Engaged on positive boundary — only pull back (λ ≤ 0)
                lambda = lambda.min(0.0);
            } else {
                // drift ≤ −backlash, engaged on negative boundary — only push (λ ≥ 0)
                lambda = lambda.max(0.0);
            }
        }

        // ------------------------------------------------------------------
        // Apply impulse: Δω = I⁻¹ · Jᵀ · λ
        //   ω_A += I_A⁻¹ · (−a₁) · λ
        //   ω_B += I_B⁻¹ · (+a₁) · λ
        //   ω_C += I_C⁻¹ · (+r·a₂) · λ
        //   ω_D += I_D⁻¹ · (−r·a₂) · λ
        // ------------------------------------------------------------------
        *body_a.omega = add3(*body_a.omega, diag_mul(body_a.inv_inertia, scale3(a1, -lambda)));
        *body_b.omega = add3(*body_b.omega, diag_mul(body_b.inv_inertia, scale3(a1, lambda)));
        *body_c.omega = add3(*body_c.omega, diag_mul(body_c.inv_inertia, scale3(a2, r * lambda)));
        *body_d.omega = add3(*body_d.omega, diag_mul(body_d.inv_inertia, scale3(a2, -r * lambda)));

        lambda
    }
}

impl fmt::Display for GearConstraint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "GearConstraint {{ ratio: {}, backlash: {}, beta: {} }}",
            self.ratio, self.backlash, self.beta
        )
    }
}
