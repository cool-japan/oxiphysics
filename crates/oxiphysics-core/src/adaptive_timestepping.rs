// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unified adaptive time-stepping module for multi-domain physics simulations.
//!
//! Provides CFL-based timestep constraints, Richardson extrapolation for
//! error-based control, and a [`UnifiedTimeStepper`] that coordinates
//! multiple physics domains with subcycling support.
//!
//! The stepper advances every domain by one common global step, the smallest
//! effective domain timestep, so the allowed timestep ranges `[min_dt, max_dt]`
//! of the domains must overlap; [`UnifiedTimeStepper::compute_global_dt`]
//! reports disjoint ranges as [`TimeStepError::IncompatibleDomains`] instead of
//! panicking.  Invalid configurations are rejected by
//! [`TimeDomainConfig::new`] and [`UnifiedTimeStepper::new`].
//!
//! # Example
//!
//! ```
//! use oxiphysics_core::adaptive_timestepping::*;
//!
//! # fn main() -> Result<(), TimeStepError> {
//! // Configure two domains: a fast fluid and a slow structure
//! let fluid = TimeDomainConfig::new(1e-6, 1e-3, 0.25, 0.9)?;
//! let structure = TimeDomainConfig::new(1e-5, 1e-2, 0.5, 0.9)?;
//!
//! let mut stepper = UnifiedTimeStepper::new(vec![fluid, structure])?;
//!
//! // After evaluating physics, update domain timesteps
//! stepper.update_domain_dt(0, 1e-4);
//! stepper.update_domain_dt(1, 5e-4);
//!
//! let global_dt = stepper.compute_global_dt()?;
//! assert!(global_dt > 0.0);
//!
//! // Domains whose timestep ranges do not overlap cannot share a global step.
//! let tiny = TimeDomainConfig::new(1e-6, 1e-4, 0.5, 0.9)?;
//! let coarse = TimeDomainConfig::new(1e-3, 1e-2, 0.5, 0.9)?;
//! let mut disjoint = UnifiedTimeStepper::new(vec![tiny, coarse])?;
//! assert!(matches!(
//!     disjoint.compute_global_dt(),
//!     Err(TimeStepError::IncompatibleDomains { .. })
//! ));
//! # Ok(())
//! # }
//! ```

use std::fmt;

use crate::error::PhysicsError;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Errors reported by the time-domain configuration and the
/// [`UnifiedTimeStepper`].
#[derive(Debug, Clone, PartialEq)]
pub enum TimeStepError {
    /// A configuration or state value lies outside its valid range.
    InvalidParameter {
        /// Index of the offending domain, or `None` when a standalone
        /// [`TimeDomainConfig`] is validated.
        domain: Option<usize>,
        /// Name of the offending field.
        parameter: &'static str,
        /// The rejected value.
        value: f64,
        /// What the field must satisfy.
        requirement: &'static str,
    },
    /// The timestep ranges of the domains have no value in common: the largest
    /// `min_dt` exceeds the smallest `max_dt`, so no single global step is
    /// allowed by every domain.
    IncompatibleDomains {
        /// Largest `min_dt` over the domains.
        global_min_dt: f64,
        /// Smallest `max_dt` over the domains.
        global_max_dt: f64,
    },
    /// The stepper has no domain and therefore no timestep.
    NoDomains,
}

impl fmt::Display for TimeStepError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TimeStepError::InvalidParameter {
                domain,
                parameter,
                value,
                requirement,
            } => {
                match domain {
                    Some(index) => write!(f, "domain {index}: ")?,
                    None => write!(f, "time-domain configuration: ")?,
                }
                write!(f, "{parameter} = {value} must be {requirement}")
            }
            TimeStepError::IncompatibleDomains {
                global_min_dt,
                global_max_dt,
            } => write!(
                f,
                "domain timestep ranges do not overlap: the largest min_dt \
                 ({global_min_dt}) exceeds the smallest max_dt ({global_max_dt})"
            ),
            TimeStepError::NoDomains => write!(f, "the time stepper has no domains"),
        }
    }
}

impl std::error::Error for TimeStepError {}

impl From<TimeStepError> for PhysicsError {
    fn from(error: TimeStepError) -> Self {
        match error {
            TimeStepError::InvalidParameter {
                parameter,
                value,
                requirement,
                ..
            } => PhysicsError::InvalidInput {
                field: parameter.to_string(),
                reason: format!("{value} must be {requirement}"),
            },
            TimeStepError::IncompatibleDomains { .. } => PhysicsError::InvalidInput {
                field: "domains".to_string(),
                reason: error.to_string(),
            },
            TimeStepError::NoDomains => PhysicsError::InvalidInput {
                field: "domains".to_string(),
                reason: error.to_string(),
            },
        }
    }
}

// ---------------------------------------------------------------------------
// Free functions: CFL / Courant / diffusive timestep constraints
// ---------------------------------------------------------------------------

/// Compute a CFL timestep from maximum velocity and grid spacing.
///
/// `dt = cfl_number * grid_spacing / max_velocity`
///
/// Returns `f64::MAX` when `max_velocity` is effectively zero.
pub fn cfl_timestep(max_velocity: f64, grid_spacing: f64, cfl_number: f64) -> f64 {
    let v = max_velocity.abs();
    if v < 1e-30 {
        return f64::MAX;
    }
    cfl_number * grid_spacing / v
}

/// Compute the Courant timestep for wave propagation.
///
/// `dt = cfl * dx / c` where `c` is the wave/sound speed.
///
/// Returns `f64::MAX` when `c` is effectively zero.
pub fn courant_dt(dx: f64, c: f64, cfl: f64) -> f64 {
    let c_abs = c.abs();
    if c_abs < 1e-30 {
        return f64::MAX;
    }
    cfl * dx / c_abs
}

/// Compute the diffusive (viscous) timestep constraint.
///
/// `dt = safety * dx^2 / (2 * nu)` (for 1-D; the factor 2 accounts for the
/// standard Von-Neumann stability limit).
///
/// Returns `f64::MAX` when `nu` is effectively zero.
pub fn diffusive_dt(dx: f64, nu: f64, safety: f64) -> f64 {
    let nu_abs = nu.abs();
    if nu_abs < 1e-30 {
        return f64::MAX;
    }
    safety * dx * dx / (2.0 * nu_abs)
}

// ---------------------------------------------------------------------------
// Richardson extrapolation helper
// ---------------------------------------------------------------------------

/// Compute a new timestep from an error estimate using Richardson
/// extrapolation / PI-controller logic.
///
/// `dt_new = dt_old * (tolerance / error)^(1 / (order + 1))`
///
/// The result is clamped to `[dt_old * 0.1, dt_old * 5.0]` to prevent
/// violent oscillations.
///
/// If `error` is effectively zero the timestep is multiplied by the maximum
/// growth factor (5×).
pub fn richardson_dt(dt_old: f64, error: f64, tolerance: f64, order: u32) -> f64 {
    const MIN_FACTOR: f64 = 0.1;
    const MAX_FACTOR: f64 = 5.0;

    if error < 1e-30 {
        return dt_old * MAX_FACTOR;
    }

    let exponent = 1.0 / (order as f64 + 1.0);
    let ratio = (tolerance / error).powf(exponent);
    let factor = ratio.clamp(MIN_FACTOR, MAX_FACTOR);
    dt_old * factor
}

// ---------------------------------------------------------------------------
// Per-domain configuration
// ---------------------------------------------------------------------------

/// Configuration for a single time domain.
#[derive(Debug, Clone)]
pub struct TimeDomainConfig {
    /// Minimum allowed timestep.
    pub min_dt: f64,
    /// Maximum allowed timestep.
    pub max_dt: f64,
    /// CFL factor (Courant number target).
    pub cfl_factor: f64,
    /// Safety factor applied after computing raw timestep (typically 0.8–0.95).
    pub safety_factor: f64,
}

impl TimeDomainConfig {
    /// Create a validated configuration.
    ///
    /// # Errors
    ///
    /// [`TimeStepError::InvalidParameter`] unless all values are finite with
    /// `0 < min_dt <= max_dt`, `cfl_factor > 0` and `0 < safety_factor <= 1`.
    pub fn new(
        min_dt: f64,
        max_dt: f64,
        cfl_factor: f64,
        safety_factor: f64,
    ) -> Result<Self, TimeStepError> {
        let config = Self {
            min_dt,
            max_dt,
            cfl_factor,
            safety_factor,
        };
        config.validate()?;
        Ok(config)
    }

    /// Check the invariants of the configuration: all values finite,
    /// `0 < min_dt <= max_dt`, `cfl_factor > 0` and `0 < safety_factor <= 1`.
    ///
    /// The fields are public, so a configuration can leave
    /// [`new`](Self::new) valid and be edited into an invalid one later; the
    /// stepper re-checks every domain before it relies on these invariants.
    ///
    /// # Errors
    ///
    /// [`TimeStepError::InvalidParameter`] naming the first offending field.
    pub fn validate(&self) -> Result<(), TimeStepError> {
        self.check(None)
    }

    /// [`validate`](Self::validate), tagging errors with the domain index.
    fn check(&self, domain: Option<usize>) -> Result<(), TimeStepError> {
        let reject = |parameter: &'static str, value: f64, requirement: &'static str| {
            Err(TimeStepError::InvalidParameter {
                domain,
                parameter,
                value,
                requirement,
            })
        };
        if !(self.min_dt.is_finite() && self.min_dt > 0.0) {
            return reject("min_dt", self.min_dt, "finite and positive");
        }
        if !(self.max_dt.is_finite() && self.max_dt >= self.min_dt) {
            return reject("max_dt", self.max_dt, "finite and at least min_dt");
        }
        if !(self.cfl_factor.is_finite() && self.cfl_factor > 0.0) {
            return reject("cfl_factor", self.cfl_factor, "finite and positive");
        }
        if !(self.safety_factor.is_finite()
            && self.safety_factor > 0.0
            && self.safety_factor <= 1.0)
        {
            return reject("safety_factor", self.safety_factor, "in (0, 1]");
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Per-domain state
// ---------------------------------------------------------------------------

/// Mutable state for a single time domain.
#[derive(Debug, Clone)]
pub struct TimeDomainState {
    /// Current timestep for this domain.
    pub current_dt: f64,
    /// Latest local truncation error estimate (set by the integrator).
    pub error_estimate: f64,
    /// Number of steps taken in this domain.
    pub step_count: u64,
    /// Configuration (immutable after construction, but stored here for
    /// convenient access).
    pub config: TimeDomainConfig,
}

impl TimeDomainState {
    /// Create a new domain state from its configuration.
    ///
    /// The initial `current_dt` is set to `config.max_dt * config.safety_factor`.
    pub fn from_config(config: TimeDomainConfig) -> Self {
        let initial_dt = config.max_dt * config.safety_factor;
        Self {
            current_dt: initial_dt,
            error_estimate: 0.0,
            step_count: 0,
            config,
        }
    }
}

// ---------------------------------------------------------------------------
// Step schedule entry
// ---------------------------------------------------------------------------

/// Describes how a single domain should be advanced during one global step.
#[derive(Debug, Clone, PartialEq)]
pub struct StepSchedule {
    /// Index of the domain.
    pub domain_idx: usize,
    /// Number of sub-steps this domain should take.
    pub n_substeps: usize,
    /// Timestep for each sub-step.
    pub substep_dt: f64,
}

// ---------------------------------------------------------------------------
// UnifiedTimeStepper
// ---------------------------------------------------------------------------

/// Coordinates adaptive time-stepping across multiple physics domains.
///
/// Each domain may have a different natural timestep; the stepper computes a
/// global step size (the minimum, after safety factors) and determines how
/// many sub-steps faster domains effectively need within that global window.
/// Domains whose natural timestep is larger than the global one take a single
/// step of size `global_dt`.
#[derive(Debug, Clone)]
pub struct UnifiedTimeStepper {
    /// Per-domain states.
    pub domains: Vec<TimeDomainState>,
    /// Current global timestep (minimum across domains, with safety).
    pub global_dt: f64,
    /// Accumulated simulation time.
    pub global_time: f64,
    /// How many sub-steps each domain takes per global step.
    pub subcycle_ratios: Vec<usize>,
}

impl UnifiedTimeStepper {
    /// Create a new stepper from per-domain configurations.
    ///
    /// A stepper without domains is allowed (its `global_dt` is `0`), but
    /// [`compute_global_dt`](Self::compute_global_dt) reports
    /// [`TimeStepError::NoDomains`] for it.
    ///
    /// # Errors
    ///
    /// [`TimeStepError::InvalidParameter`] (tagged with the domain index) if a
    /// configuration violates the invariants checked by
    /// [`TimeDomainConfig::validate`]; configurations built by struct literal
    /// bypass [`TimeDomainConfig::new`], so they are re-checked here.
    pub fn new(configs: Vec<TimeDomainConfig>) -> Result<Self, TimeStepError> {
        for (index, config) in configs.iter().enumerate() {
            config.check(Some(index))?;
        }
        let domains: Vec<TimeDomainState> = configs
            .into_iter()
            .map(TimeDomainState::from_config)
            .collect();
        let n = domains.len();
        let mut stepper = Self {
            domains,
            global_dt: 0.0,
            global_time: 0.0,
            subcycle_ratios: vec![1; n],
        };
        stepper.global_dt = stepper.raw_global_dt();
        stepper.recompute_subcycle_ratios();
        Ok(stepper)
    }

    // -- internal helpers ---------------------------------------------------

    /// Compute the raw global dt (minimum of domain dts with safety); `0` for a
    /// stepper without domains.
    fn raw_global_dt(&self) -> f64 {
        self.domains
            .iter()
            .map(|d| d.current_dt * d.config.safety_factor)
            .reduce(f64::min)
            .unwrap_or(0.0)
    }

    /// Recompute subcycle ratios from current domain dts.
    fn recompute_subcycle_ratios(&mut self) {
        // `subcycle_ratios` is public: keep it one-to-one with the domains.
        self.subcycle_ratios.resize(self.domains.len(), 1);
        if self.global_dt <= 0.0 {
            // Fallback: all 1.
            self.subcycle_ratios.fill(1);
            return;
        }
        for (ratio_slot, domain) in self.subcycle_ratios.iter_mut().zip(&self.domains) {
            let effective_dt = domain.current_dt * domain.config.safety_factor;
            // Ratio: how many global steps fit in one domain step.
            // If domain is slower (larger dt), it takes 1 sub-step of global_dt.
            // If domain is faster (smaller dt), it needs ceil(global_dt / domain_dt).
            *ratio_slot = if effective_dt >= self.global_dt {
                1
            } else {
                // Domain is faster than global; needs multiple sub-steps.
                ((self.global_dt / effective_dt).ceil() as usize).max(1)
            };
        }
    }

    // -- public API ---------------------------------------------------------

    /// Compute and update the global timestep (minimum effective domain
    /// timestep), restricted to the range every domain allows.
    ///
    /// The result lies in `[max_i min_dt_i, min_i max_dt_i]`.  Returns the new
    /// `global_dt`; on error the stepper is left unchanged.
    ///
    /// # Errors
    ///
    /// * [`TimeStepError::NoDomains`] if the stepper has no domain;
    /// * [`TimeStepError::InvalidParameter`] if a domain configuration (public
    ///   field, possibly edited after construction) or its `current_dt` is not
    ///   valid;
    /// * [`TimeStepError::IncompatibleDomains`] if the timestep ranges of the
    ///   domains are disjoint, i.e. the largest `min_dt` exceeds the smallest
    ///   `max_dt`.
    pub fn compute_global_dt(&mut self) -> Result<f64, TimeStepError> {
        if self.domains.is_empty() {
            return Err(TimeStepError::NoDomains);
        }
        for (index, domain) in self.domains.iter().enumerate() {
            domain.config.check(Some(index))?;
            if !(domain.current_dt.is_finite() && domain.current_dt > 0.0) {
                return Err(TimeStepError::InvalidParameter {
                    domain: Some(index),
                    parameter: "current_dt",
                    value: domain.current_dt,
                    requirement: "finite and positive",
                });
            }
        }
        let global_min = self
            .domains
            .iter()
            .map(|d| d.config.min_dt)
            .fold(0.0_f64, f64::max);
        let global_max = self
            .domains
            .iter()
            .map(|d| d.config.max_dt)
            .fold(f64::INFINITY, f64::min);
        if global_min > global_max {
            return Err(TimeStepError::IncompatibleDomains {
                global_min_dt: global_min,
                global_max_dt: global_max,
            });
        }
        // `max` then `min` instead of `f64::clamp`, which panics on min > max.
        self.global_dt = self.raw_global_dt().max(global_min).min(global_max);
        self.recompute_subcycle_ratios();
        Ok(self.global_dt)
    }

    /// Recompute the subcycle ratios based on the current `global_dt` and
    /// each domain's `current_dt`.
    pub fn compute_subcycle_ratios(&mut self) {
        self.recompute_subcycle_ratios();
    }

    /// Advance the global simulation time by one `global_dt` and increment
    /// each domain's step counter by its subcycle ratio.
    pub fn advance_global_time(&mut self) {
        self.global_time += self.global_dt;
        for (domain, &ratio) in self.domains.iter_mut().zip(&self.subcycle_ratios) {
            domain.step_count = domain.step_count.saturating_add(ratio as u64);
        }
    }

    /// Update a domain's suggested timestep.
    ///
    /// The new value is clamped to the domain's `[min_dt, max_dt]` range (an
    /// infinite suggestion means "no constraint" and becomes `max_dt`).
    ///
    /// Returns `true` if the timestep was updated; `false`, leaving the domain
    /// unchanged, if `domain_idx` is out of range or `new_dt` is NaN or not
    /// positive.
    pub fn update_domain_dt(&mut self, domain_idx: usize, new_dt: f64) -> bool {
        if new_dt.is_nan() || new_dt <= 0.0 {
            return false;
        }
        match self.domains.get_mut(domain_idx) {
            Some(domain) => {
                // `max` then `min` instead of `f64::clamp`, which panics when the
                // public bounds were edited into min_dt > max_dt.
                domain.current_dt = new_dt.max(domain.config.min_dt).min(domain.config.max_dt);
                true
            }
            None => false,
        }
    }

    /// Record an error estimate for a domain (used by Richardson /
    /// embedded-pair integrators).
    ///
    /// Returns `true` if `domain_idx` was valid, `false` otherwise.
    pub fn update_domain_error(&mut self, domain_idx: usize, error: f64) -> bool {
        if let Some(domain) = self.domains.get_mut(domain_idx) {
            domain.error_estimate = error;
            true
        } else {
            false
        }
    }

    /// Generate the step schedule describing how each domain should be
    /// advanced during the next global step.
    pub fn step_schedule(&self) -> Vec<StepSchedule> {
        self.domains
            .iter()
            .enumerate()
            .map(|(i, _domain)| {
                let n = self.subcycle_ratios.get(i).copied().unwrap_or(1);
                let sub_dt = if n > 0 {
                    self.global_dt / n as f64
                } else {
                    self.global_dt
                };
                StepSchedule {
                    domain_idx: i,
                    n_substeps: n,
                    substep_dt: sub_dt,
                }
            })
            .collect()
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -- Free-function tests ------------------------------------------------

    #[test]
    fn test_cfl_timestep_basic() {
        let dt = cfl_timestep(10.0, 1.0, 0.5);
        assert!((dt - 0.05).abs() < 1e-12, "dt = {dt}");
    }

    #[test]
    fn test_cfl_timestep_zero_velocity_returns_max() {
        let dt = cfl_timestep(0.0, 1.0, 0.5);
        assert_eq!(dt, f64::MAX);
    }

    #[test]
    fn test_courant_dt_basic() {
        let dt = courant_dt(0.1, 340.0, 0.8);
        let expected = 0.8 * 0.1 / 340.0;
        assert!((dt - expected).abs() < 1e-14, "dt = {dt}");
    }

    #[test]
    fn test_courant_dt_zero_speed() {
        assert_eq!(courant_dt(0.1, 0.0, 0.8), f64::MAX);
    }

    #[test]
    fn test_diffusive_dt_basic() {
        let dt = diffusive_dt(0.01, 1e-3, 0.9);
        let expected = 0.9 * 0.01 * 0.01 / (2.0 * 1e-3);
        assert!((dt - expected).abs() < 1e-14, "dt = {dt}");
    }

    #[test]
    fn test_diffusive_dt_zero_nu() {
        assert_eq!(diffusive_dt(0.01, 0.0, 0.9), f64::MAX);
    }

    // -- Richardson extrapolation tests ------------------------------------

    #[test]
    fn test_richardson_increases_dt_when_error_below_tolerance() {
        let dt_old = 0.01;
        let error = 1e-6;
        let tolerance = 1e-4;
        let dt_new = richardson_dt(dt_old, error, tolerance, 2);
        assert!(
            dt_new > dt_old,
            "dt_new={dt_new} should be > dt_old={dt_old}"
        );
    }

    #[test]
    fn test_richardson_decreases_dt_when_error_above_tolerance() {
        let dt_old = 0.01;
        let error = 1e-2;
        let tolerance = 1e-4;
        let dt_new = richardson_dt(dt_old, error, tolerance, 2);
        assert!(
            dt_new < dt_old,
            "dt_new={dt_new} should be < dt_old={dt_old}"
        );
    }

    #[test]
    fn test_richardson_zero_error_gives_max_growth() {
        let dt_old = 0.01;
        let dt_new = richardson_dt(dt_old, 0.0, 1e-4, 2);
        assert!((dt_new - dt_old * 5.0).abs() < 1e-14);
    }

    #[test]
    fn test_richardson_clamped_growth() {
        // Very tiny error => growth capped at 5×.
        let dt_new = richardson_dt(0.01, 1e-30, 1.0, 1);
        assert!((dt_new - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_richardson_clamped_shrink() {
        // Huge error => shrink capped at 0.1×.
        let dt_new = richardson_dt(0.01, 1e10, 1e-4, 2);
        assert!((dt_new - 0.001).abs() < 1e-10);
    }

    // -- Single-domain convergence -----------------------------------------

    #[test]
    fn test_single_domain_dt_converges() {
        let cfg = TimeDomainConfig::new(1e-6, 1e-2, 0.5, 0.9).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();

        // Simulate a controller loop: update dt based on "physics".
        for _ in 0..20 {
            let target_dt = 5e-4;
            stepper.update_domain_dt(0, target_dt);
            stepper.compute_global_dt().unwrap();
        }

        let effective = stepper.global_dt;
        let expected = 5e-4 * 0.9;
        assert!(
            (effective - expected).abs() < 1e-10,
            "effective={effective}, expected={expected}"
        );
    }

    // -- Two-domain subcycling ---------------------------------------------

    #[test]
    fn test_two_domains_subcycling() {
        let fast = TimeDomainConfig::new(1e-6, 1e-2, 0.25, 0.9).unwrap();
        let slow = TimeDomainConfig::new(1e-5, 1e-1, 0.5, 0.9).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![fast, slow]).unwrap();

        // Fast domain wants dt = 1e-4, slow domain wants dt = 1e-3.
        stepper.update_domain_dt(0, 1e-4);
        stepper.update_domain_dt(1, 1e-3);
        stepper.compute_global_dt().unwrap();

        // Global dt should be driven by the fast domain.
        let g = stepper.global_dt;
        let fast_eff = 1e-4 * 0.9;
        assert!(
            (g - fast_eff).abs() < 1e-12,
            "global_dt={g}, expected={fast_eff}"
        );

        // Fast domain: 1 sub-step; slow domain: 1 sub-step (its dt is larger).
        assert_eq!(stepper.subcycle_ratios[0], 1);
        assert_eq!(stepper.subcycle_ratios[1], 1);
    }

    #[test]
    fn test_global_dt_is_minimum_across_domains() {
        let d1 = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let d2 = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let d3 = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![d1, d2, d3]).unwrap();

        stepper.update_domain_dt(0, 0.1);
        stepper.update_domain_dt(1, 0.05);
        stepper.update_domain_dt(2, 0.2);
        let g = stepper.compute_global_dt().unwrap();

        // safety=1.0, so global = min(0.1, 0.05, 0.2) = 0.05
        assert!((g - 0.05).abs() < 1e-14, "global_dt={g}");
    }

    // -- Schedule generation -----------------------------------------------

    #[test]
    fn test_schedule_single_domain() {
        let cfg = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        stepper.update_domain_dt(0, 0.01);
        stepper.compute_global_dt().unwrap();

        let sched = stepper.step_schedule();
        assert_eq!(sched.len(), 1);
        assert_eq!(sched[0].domain_idx, 0);
        assert_eq!(sched[0].n_substeps, 1);
        assert!((sched[0].substep_dt - 0.01).abs() < 1e-14);
    }

    #[test]
    fn test_schedule_multi_domain() {
        let d1 = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let d2 = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![d1, d2]).unwrap();

        stepper.update_domain_dt(0, 0.01);
        stepper.update_domain_dt(1, 0.1);
        stepper.compute_global_dt().unwrap();

        let sched = stepper.step_schedule();
        assert_eq!(sched.len(), 2);

        // Domain 0 drives global_dt = 0.01.
        assert_eq!(sched[0].n_substeps, 1);
        assert!((sched[0].substep_dt - 0.01).abs() < 1e-14);

        // Domain 1 has larger dt, so also 1 sub-step of global_dt.
        assert_eq!(sched[1].n_substeps, 1);
        assert!((sched[1].substep_dt - 0.01).abs() < 1e-14);
    }

    // -- Advance time -------------------------------------------------------

    #[test]
    fn test_advance_global_time() {
        let cfg = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        stepper.update_domain_dt(0, 0.01);
        stepper.compute_global_dt().unwrap();

        let dt = stepper.global_dt;
        stepper.advance_global_time();
        assert!((stepper.global_time - dt).abs() < 1e-14);
        assert_eq!(stepper.domains[0].step_count, 1); // 0 from init + 1 advance
    }

    #[test]
    fn test_advance_accumulates() {
        let cfg = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        stepper.update_domain_dt(0, 0.01);
        stepper.compute_global_dt().unwrap();

        for _ in 0..100 {
            stepper.advance_global_time();
        }
        let expected_time = 100.0 * stepper.global_dt;
        assert!(
            (stepper.global_time - expected_time).abs() < 1e-10,
            "time={}, expected={}",
            stepper.global_time,
            expected_time
        );
    }

    // -- Domain error updates -----------------------------------------------

    #[test]
    fn test_update_domain_error() {
        let cfg = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        assert!(stepper.update_domain_error(0, 1e-5));
        assert!((stepper.domains[0].error_estimate - 1e-5).abs() < 1e-20);
    }

    #[test]
    fn test_update_domain_error_invalid_index() {
        let cfg = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        assert!(!stepper.update_domain_error(5, 1e-5));
    }

    #[test]
    fn test_update_domain_dt_clamps() {
        let cfg = TimeDomainConfig::new(1e-6, 1e-2, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();

        // Try to set dt above max.
        stepper.update_domain_dt(0, 1.0);
        assert!((stepper.domains[0].current_dt - 1e-2).abs() < 1e-14);

        // Try to set dt below min.
        stepper.update_domain_dt(0, 1e-10);
        assert!((stepper.domains[0].current_dt - 1e-6).abs() < 1e-14);
    }

    #[test]
    fn test_update_domain_dt_invalid_index() {
        let cfg = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        assert!(!stepper.update_domain_dt(99, 0.01));
    }

    // -- Richardson + stepper integration -----------------------------------

    #[test]
    fn test_richardson_driven_adaptation() {
        let cfg = TimeDomainConfig::new(1e-8, 1.0, 0.5, 0.9).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        stepper.update_domain_dt(0, 0.01);

        let tolerance = 1e-4;
        let order = 2_u32;

        // Simulate: error is large, dt should shrink.
        let error = 1e-2;
        let new_dt = richardson_dt(stepper.domains[0].current_dt, error, tolerance, order);
        stepper.update_domain_dt(0, new_dt);
        stepper.compute_global_dt().unwrap();
        assert!(stepper.global_dt < 0.01 * 0.9, "should have shrunk");

        // Simulate: error is tiny, dt should grow.
        let error2 = 1e-10;
        let new_dt2 = richardson_dt(stepper.domains[0].current_dt, error2, tolerance, order);
        stepper.update_domain_dt(0, new_dt2);
        stepper.compute_global_dt().unwrap();
        assert!(
            stepper.domains[0].current_dt > stepper.domains[0].config.min_dt,
            "should have grown"
        );
    }

    // -- Subcycling with a domain faster than global -----------------------

    #[test]
    fn test_subcycling_fast_domain() {
        // Construct: d0 has safety 1.0 and dt=0.01,
        //            d1 has safety 0.5, so effective dt = dt*0.5.
        // If d1.current_dt = 0.01, effective = 0.005.
        // Global dt = min(0.01, 0.005) = 0.005.
        // d0 effective = 0.01 >= 0.005 → 1 sub-step.
        // d1 effective = 0.005 >= 0.005 → 1 sub-step.
        //
        // Now if we manually set global_dt larger via a different route,
        // we can test true subcycling:
        let d0 = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let d1 = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![d0, d1]).unwrap();

        // d0 wants dt=0.1, d1 wants dt=0.025
        stepper.update_domain_dt(0, 0.1);
        stepper.update_domain_dt(1, 0.025);
        stepper.compute_global_dt().unwrap();

        // global = 0.025, d0 effective 0.1 >= 0.025 → 1 sub-step
        assert_eq!(stepper.subcycle_ratios[0], 1);
        assert_eq!(stepper.subcycle_ratios[1], 1);
        assert!((stepper.global_dt - 0.025).abs() < 1e-14);
    }

    // -- Edge cases ---------------------------------------------------------

    #[test]
    fn test_empty_domains() {
        let stepper = UnifiedTimeStepper::new(vec![]).unwrap();
        assert_eq!(stepper.domains.len(), 0);
        assert_eq!(stepper.subcycle_ratios.len(), 0);
        assert!(stepper.step_schedule().is_empty());
    }

    #[test]
    fn test_identical_domains() {
        let c1 = TimeDomainConfig::new(1e-6, 1.0, 0.5, 0.9).unwrap();
        let c2 = TimeDomainConfig::new(1e-6, 1.0, 0.5, 0.9).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![c1, c2]).unwrap();

        stepper.update_domain_dt(0, 0.01);
        stepper.update_domain_dt(1, 0.01);
        stepper.compute_global_dt().unwrap();

        assert_eq!(stepper.subcycle_ratios[0], 1);
        assert_eq!(stepper.subcycle_ratios[1], 1);
        assert!((stepper.global_dt - 0.01 * 0.9).abs() < 1e-14);
    }

    #[test]
    fn test_global_dt_respects_min_dt_bound() {
        // Both domains have min_dt = 0.001. Even if the raw safety-adjusted
        // dt would be smaller, it should be clamped.
        let c1 = TimeDomainConfig::new(0.001, 1.0, 0.5, 0.01).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![c1]).unwrap();
        stepper.update_domain_dt(0, 0.001); // safety 0.01 → effective 1e-5
        let g = stepper.compute_global_dt().unwrap();
        assert!(g >= 0.001, "global_dt={g} should be >= min_dt=0.001");
    }

    #[test]
    fn test_negative_velocity_cfl() {
        // Negative velocity should be treated as abs.
        let dt = cfl_timestep(-100.0, 1.0, 0.5);
        let expected = 0.5 * 1.0 / 100.0;
        assert!((dt - expected).abs() < 1e-14);
    }

    #[test]
    fn test_diffusive_dt_negative_nu() {
        // Negative viscosity should use abs.
        let dt = diffusive_dt(0.1, -0.01, 0.9);
        let expected = 0.9 * 0.01 / 0.02;
        assert!((dt - expected).abs() < 1e-14);
    }

    // -- Configuration validation -------------------------------------------

    #[test]
    fn test_config_new_accepts_boundary_values() {
        // min_dt == max_dt and safety == 1 are valid.
        assert!(TimeDomainConfig::new(1e-3, 1e-3, 0.5, 1.0).is_ok());
        assert!(TimeDomainConfig::new(1e-12, 1e12, 1e-9, 1e-9).is_ok());
    }

    #[test]
    fn test_config_new_rejects_invalid_values_naming_the_field() {
        let nan = f64::NAN;
        let inf = f64::INFINITY;
        let cases: [([f64; 4], &str); 12] = [
            ([0.0, 1.0, 0.5, 0.9], "min_dt"),
            ([-1e-3, 1.0, 0.5, 0.9], "min_dt"),
            ([nan, 1.0, 0.5, 0.9], "min_dt"),
            ([inf, inf, 0.5, 0.9], "min_dt"),
            ([1e-3, 1e-4, 0.5, 0.9], "max_dt"),
            ([1e-3, nan, 0.5, 0.9], "max_dt"),
            ([1e-3, inf, 0.5, 0.9], "max_dt"),
            ([1e-3, 1.0, 0.0, 0.9], "cfl_factor"),
            ([1e-3, 1.0, nan, 0.9], "cfl_factor"),
            ([1e-3, 1.0, 0.5, 0.0], "safety_factor"),
            ([1e-3, 1.0, 0.5, 1.5], "safety_factor"),
            ([1e-3, 1.0, 0.5, nan], "safety_factor"),
        ];
        for ([min_dt, max_dt, cfl, safety], field) in cases {
            match TimeDomainConfig::new(min_dt, max_dt, cfl, safety) {
                Err(TimeStepError::InvalidParameter {
                    domain, parameter, ..
                }) => {
                    assert_eq!(domain, None);
                    assert_eq!(parameter, field, "{min_dt} {max_dt} {cfl} {safety}");
                }
                other => panic!("({min_dt}, {max_dt}, {cfl}, {safety}) gave {other:?}"),
            }
        }
    }

    #[test]
    fn test_stepper_new_rechecks_struct_literal_configs() {
        // Fields are public, so a config can bypass TimeDomainConfig::new.
        let good = TimeDomainConfig::new(1e-6, 1.0, 0.5, 0.9).unwrap();
        let bad = TimeDomainConfig {
            min_dt: 1.0,
            max_dt: 1e-3,
            cfl_factor: 0.5,
            safety_factor: 0.9,
        };
        match UnifiedTimeStepper::new(vec![good, bad]) {
            Err(TimeStepError::InvalidParameter {
                domain, parameter, ..
            }) => {
                assert_eq!(domain, Some(1));
                assert_eq!(parameter, "max_dt");
            }
            other => panic!("expected an invalid-parameter error, got {other:?}"),
        }
    }

    // -- Disjoint domain ranges (the former clamp panic) --------------------

    #[test]
    fn test_disjoint_domain_ranges_return_an_error_instead_of_panicking() {
        // [1e-6, 1e-4] and [1e-3, 1e-2] have no common timestep: the largest
        // lower bound (1e-3) exceeds the smallest upper bound (1e-4), which made
        // f64::clamp panic.
        let fast = TimeDomainConfig::new(1e-6, 1e-4, 0.25, 0.9).unwrap();
        let slow = TimeDomainConfig::new(1e-3, 1e-2, 0.5, 0.9).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![fast, slow]).unwrap();
        stepper.update_domain_dt(0, 5e-5);
        stepper.update_domain_dt(1, 5e-3);
        let before = (
            stepper.global_dt,
            stepper.subcycle_ratios.clone(),
            stepper.global_time,
        );
        assert_eq!(
            stepper.compute_global_dt(),
            Err(TimeStepError::IncompatibleDomains {
                global_min_dt: 1e-3,
                global_max_dt: 1e-4,
            })
        );
        // The failed call leaves the stepper untouched.
        assert_eq!(
            (
                stepper.global_dt,
                stepper.subcycle_ratios.clone(),
                stepper.global_time
            ),
            before
        );
    }

    #[test]
    fn test_overlapping_ranges_with_different_bounds_are_clamped_to_the_overlap() {
        // The ranges overlap in [1e-3, 1e-2].  The raw minimum (5e-4 * 0.9) lies
        // below the overlap and is raised to its lower bound.
        let a = TimeDomainConfig::new(1e-5, 1e-2, 0.5, 0.9).unwrap();
        let b = TimeDomainConfig::new(1e-3, 1e-1, 0.5, 0.9).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![a, b]).unwrap();
        stepper.update_domain_dt(0, 5e-4);
        stepper.update_domain_dt(1, 5e-2);
        let g = stepper.compute_global_dt().unwrap();
        assert!((g - 1e-3).abs() < 1e-15, "global_dt = {g}");
        // A raw minimum inside the overlap is left alone.
        stepper.update_domain_dt(0, 1e-2);
        stepper.update_domain_dt(1, 1e-1);
        let g = stepper.compute_global_dt().unwrap();
        assert!((g - 9e-3).abs() < 1e-15, "global_dt = {g}");
        // Touching ranges share exactly one value.
        let low = TimeDomainConfig::new(1e-4, 1e-3, 0.5, 1.0).unwrap();
        let high = TimeDomainConfig::new(1e-3, 1e-2, 0.5, 1.0).unwrap();
        let mut touching = UnifiedTimeStepper::new(vec![low, high]).unwrap();
        touching.update_domain_dt(0, 5e-4);
        touching.update_domain_dt(1, 5e-3);
        assert_eq!(touching.compute_global_dt(), Ok(1e-3));
    }

    #[test]
    fn test_edited_public_fields_cannot_make_the_stepper_panic() {
        let cfg = TimeDomainConfig::new(1e-6, 1.0, 0.5, 0.9).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        // min_dt > max_dt after construction: update_domain_dt used to call clamp.
        stepper.domains[0].config.min_dt = 2.0;
        assert!(stepper.update_domain_dt(0, 0.5));
        assert!(stepper.domains[0].current_dt.is_finite());
        match stepper.compute_global_dt() {
            Err(TimeStepError::InvalidParameter {
                domain, parameter, ..
            }) => {
                assert_eq!(domain, Some(0));
                assert_eq!(parameter, "max_dt");
            }
            other => panic!("expected an invalid-parameter error, got {other:?}"),
        }
        // A non-finite current_dt is reported, not propagated.
        stepper.domains[0].config.min_dt = 1e-6;
        stepper.domains[0].current_dt = f64::NAN;
        match stepper.compute_global_dt() {
            Err(TimeStepError::InvalidParameter { parameter, .. }) => {
                assert_eq!(parameter, "current_dt")
            }
            other => panic!("expected an invalid-parameter error, got {other:?}"),
        }
        // Subcycle ratios and step_schedule survive a resized public vector.
        stepper.domains[0].current_dt = 0.1;
        stepper.subcycle_ratios.clear();
        assert_eq!(stepper.step_schedule()[0].n_substeps, 1);
        assert!(stepper.compute_global_dt().is_ok());
        assert_eq!(stepper.subcycle_ratios.len(), 1);
    }

    #[test]
    fn test_update_domain_dt_rejects_nan_and_non_positive_suggestions() {
        let cfg = TimeDomainConfig::new(1e-6, 1e-2, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        assert!(stepper.update_domain_dt(0, 5e-3));
        for bad in [f64::NAN, 0.0, -1e-3, f64::NEG_INFINITY] {
            assert!(!stepper.update_domain_dt(0, bad), "{bad}");
            assert_eq!(stepper.domains[0].current_dt, 5e-3, "{bad}");
        }
        // An infinite suggestion means "no constraint" and becomes max_dt.
        assert!(stepper.update_domain_dt(0, f64::INFINITY));
        assert_eq!(stepper.domains[0].current_dt, 1e-2);
    }

    #[test]
    fn test_empty_stepper_has_no_timestep() {
        let mut stepper = UnifiedTimeStepper::new(vec![]).unwrap();
        assert_eq!(stepper.global_dt, 0.0);
        assert_eq!(stepper.compute_global_dt(), Err(TimeStepError::NoDomains));
        stepper.advance_global_time();
        assert_eq!(stepper.global_time, 0.0);
    }

    #[test]
    fn test_step_counters_saturate_instead_of_overflowing() {
        let cfg = TimeDomainConfig::new(1e-6, 1.0, 0.5, 1.0).unwrap();
        let mut stepper = UnifiedTimeStepper::new(vec![cfg]).unwrap();
        stepper.domains[0].step_count = u64::MAX - 1;
        stepper.advance_global_time();
        stepper.advance_global_time();
        assert_eq!(stepper.domains[0].step_count, u64::MAX);
    }

    #[test]
    fn test_time_step_error_messages_and_physics_error_conversion() {
        let invalid = TimeDomainConfig::new(0.0, 1.0, 0.5, 0.9).unwrap_err();
        assert_eq!(
            invalid.to_string(),
            "time-domain configuration: min_dt = 0 must be finite and positive"
        );
        let tagged = TimeStepError::InvalidParameter {
            domain: Some(2),
            parameter: "safety_factor",
            value: 1.5,
            requirement: "in (0, 1]",
        };
        assert_eq!(
            tagged.to_string(),
            "domain 2: safety_factor = 1.5 must be in (0, 1]"
        );
        let incompatible = TimeStepError::IncompatibleDomains {
            global_min_dt: 1e-3,
            global_max_dt: 1e-4,
        };
        assert!(incompatible.to_string().contains("do not overlap"));
        assert_eq!(
            TimeStepError::NoDomains.to_string(),
            "the time stepper has no domains"
        );
        match PhysicsError::from(tagged) {
            PhysicsError::InvalidInput { field, reason } => {
                assert_eq!(field, "safety_factor");
                assert_eq!(reason, "1.5 must be in (0, 1]");
            }
            other => panic!("unexpected conversion {other:?}"),
        }
        assert!(matches!(
            PhysicsError::from(incompatible),
            PhysicsError::InvalidInput { .. }
        ));
        // The error type can be boxed as std::error::Error.
        let boxed: Box<dyn std::error::Error> = Box::new(TimeStepError::NoDomains);
        assert!(!boxed.to_string().is_empty());
    }
}
