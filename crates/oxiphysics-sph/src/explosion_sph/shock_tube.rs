// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Shock tube (Riemann) problem and its exact solution for an ideal gas.
//!
//! The solver follows Toro, *Riemann Solvers and Numerical Methods for Fluid
//! Dynamics* (3rd ed., 2009), chapter 4: the pressure `p*` of the star region
//! is the root of `f_L(p) + f_R(p) + (u_R − u_L) = 0`, the contact velocity is
//! `u* = ½(u_L + u_R) + ½(f_R(p*) − f_L(p*))`, and the solution at `x/t` is
//! sampled from the wave pattern (a shock or a rarefaction fan on each side of
//! the contact discontinuity, or a vacuum if the two rarefactions separate).

/// Initial conditions for one side of a shock tube.
#[derive(Debug, Clone, Copy)]
pub struct ShockTubeState {
    /// Density (kg/m^3).
    pub rho: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Velocity (m/s).
    pub velocity: f64,
}

/// Sod shock tube problem definition.
#[derive(Debug, Clone)]
pub struct SodShockTube {
    /// Left state.
    pub left: ShockTubeState,
    /// Right state.
    pub right: ShockTubeState,
    /// Ratio of specific heats.
    pub gamma: f64,
    /// Domain length (m).
    pub length: f64,
    /// Diaphragm position (m).
    pub diaphragm: f64,
}

/// Relative tolerance of the Newton iteration for the star-region pressure.
const PRESSURE_TOLERANCE: f64 = 1.0e-13;

/// Maximum number of Newton iterations for the star-region pressure.
const MAX_NEWTON_ITERATIONS: usize = 100;

/// Wave structure between the two outer waves of a Riemann problem.
#[derive(Debug, Clone, Copy, PartialEq)]
enum StarRegion {
    /// Two non-linear waves enclose a contact discontinuity moving with
    /// velocity `velocity`; both star states have pressure `pressure`.
    Contact { pressure: f64, velocity: f64 },
    /// The two rarefaction waves separate and leave a vacuum between them.
    Vacuum,
}

impl ShockTubeState {
    /// Sound speed `sqrt(γ p / ρ)` of this state.
    fn sound_speed(&self, gamma: f64) -> f64 {
        (gamma * self.pressure / self.rho).sqrt()
    }
}

/// The pressure function `f_K(p)` of side `K` and its derivative `f_K'(p)`
/// (Toro eqs. 4.6 and 4.37): the velocity change across the wave that connects
/// the state `K` to the star pressure `p`.
fn pressure_function(p: f64, state: &ShockTubeState, gamma: f64) -> (f64, f64) {
    let (rho_k, p_k) = (state.rho, state.pressure);
    let c_k = state.sound_speed(gamma);
    if p > p_k {
        // Shock: f = (p − p_K) sqrt(A_K / (p + B_K)).
        let a_k = 2.0 / ((gamma + 1.0) * rho_k);
        let b_k = (gamma - 1.0) / (gamma + 1.0) * p_k;
        let root = (a_k / (p + b_k)).sqrt();
        let f = (p - p_k) * root;
        let df = root * (1.0 - 0.5 * (p - p_k) / (p + b_k));
        (f, df)
    } else {
        // Rarefaction: f = 2 c_K/(γ−1) [(p/p_K)^((γ−1)/(2γ)) − 1].
        let ratio = p / p_k;
        let f = 2.0 * c_k / (gamma - 1.0) * (ratio.powf((gamma - 1.0) / (2.0 * gamma)) - 1.0);
        let df = ratio.powf(-(gamma + 1.0) / (2.0 * gamma)) / (rho_k * c_k);
        (f, df)
    }
}

/// Initial guess for the star pressure (Toro section 4.3.2): the
/// primitive-variable estimate when the wave pattern is mild, otherwise the
/// two-rarefaction or the two-shock approximation.
fn initial_pressure_guess(left: &ShockTubeState, right: &ShockTubeState, gamma: f64) -> f64 {
    let c_l = left.sound_speed(gamma);
    let c_r = right.sound_speed(gamma);
    let du = right.velocity - left.velocity;
    let p_min = left.pressure.min(right.pressure);
    let p_max = left.pressure.max(right.pressure);
    let p_linear =
        0.5 * (left.pressure + right.pressure) - 0.125 * du * (left.rho + right.rho) * (c_l + c_r);
    let p_linear = p_linear.max(PRESSURE_TOLERANCE * p_min);
    if p_max / p_min <= 2.0 && (p_min..=p_max).contains(&p_linear) {
        return p_linear;
    }
    if p_linear < p_min {
        // Two-rarefaction approximation.
        let z = (gamma - 1.0) / (2.0 * gamma);
        let numerator = c_l + c_r - 0.5 * (gamma - 1.0) * du;
        let denominator = c_l / left.pressure.powf(z) + c_r / right.pressure.powf(z);
        return (numerator / denominator).powf(1.0 / z);
    }
    // Two-shock approximation.
    let weight = |state: &ShockTubeState| {
        let a_k = 2.0 / ((gamma + 1.0) * state.rho);
        let b_k = (gamma - 1.0) / (gamma + 1.0) * state.pressure;
        (a_k / (p_linear + b_k)).sqrt()
    };
    let (g_l, g_r) = (weight(left), weight(right));
    ((g_l * left.pressure + g_r * right.pressure - du) / (g_l + g_r))
        .max(PRESSURE_TOLERANCE * p_min)
}

impl SodShockTube {
    /// Standard Sod problem.
    pub fn standard() -> Self {
        Self {
            left: ShockTubeState {
                rho: 1.0,
                pressure: 1.0,
                velocity: 0.0,
            },
            right: ShockTubeState {
                rho: 0.125,
                pressure: 0.1,
                velocity: 0.0,
            },
            gamma: 1.4,
            length: 1.0,
            diaphragm: 0.5,
        }
    }

    /// Whether the problem has positive finite densities and pressures on both
    /// sides and `γ > 1`, which the exact solver requires.
    fn is_physical(&self) -> bool {
        let values = [
            self.left.rho,
            self.left.pressure,
            self.left.velocity,
            self.right.rho,
            self.right.pressure,
            self.right.velocity,
            self.gamma,
        ];
        values.iter().all(|v| v.is_finite())
            && self.left.rho > 0.0
            && self.right.rho > 0.0
            && self.left.pressure > 0.0
            && self.right.pressure > 0.0
            && self.gamma > 1.0
    }

    /// Solve for the star region: the pressure root of
    /// `f_L(p) + f_R(p) + (u_R − u_L) = 0` and the contact velocity, or the
    /// vacuum case. Returns `None` for unphysical input.
    fn star_region_solution(&self) -> Option<StarRegion> {
        if !self.is_physical() {
            return None;
        }
        let (left, right, gamma) = (&self.left, &self.right, self.gamma);
        let du = right.velocity - left.velocity;
        let c_l = left.sound_speed(gamma);
        let c_r = right.sound_speed(gamma);
        // Pressure positivity condition (Toro eq. 4.40): otherwise the two
        // rarefaction waves are so strong that they generate a vacuum.
        if 2.0 * (c_l + c_r) / (gamma - 1.0) <= du {
            return Some(StarRegion::Vacuum);
        }
        let mut p = initial_pressure_guess(left, right, gamma);
        for _ in 0..MAX_NEWTON_ITERATIONS {
            let (f_l, df_l) = pressure_function(p, left, gamma);
            let (f_r, df_r) = pressure_function(p, right, gamma);
            let derivative = df_l + df_r;
            if derivative.is_nan() || derivative <= 0.0 {
                break;
            }
            let mut next = p - (f_l + f_r + du) / derivative;
            if next <= 0.0 {
                next = 0.5 * p;
            }
            let change = 2.0 * (next - p).abs() / (next + p);
            p = next;
            if change < PRESSURE_TOLERANCE {
                break;
            }
        }
        let (f_l, _) = pressure_function(p, left, gamma);
        let (f_r, _) = pressure_function(p, right, gamma);
        let velocity = 0.5 * (left.velocity + right.velocity) + 0.5 * (f_r - f_l);
        Some(StarRegion::Contact {
            pressure: p,
            velocity,
        })
    }

    /// Pressure `p*` and velocity `u*` of the star region between the two
    /// non-linear waves (the pressure and velocity of the contact
    /// discontinuity).
    ///
    /// Returns `None` if the left and right rarefaction waves generate a
    /// vacuum, or if the problem is unphysical (non-positive or non-finite
    /// densities and pressures, or `γ ≤ 1`).
    pub fn star_state(&self) -> Option<(f64, f64)> {
        match self.star_region_solution()? {
            StarRegion::Contact { pressure, velocity } => Some((pressure, velocity)),
            StarRegion::Vacuum => None,
        }
    }

    /// State inside a left rarefaction fan at the similarity coordinate
    /// `speed = (x − x₀)/t` (Toro eq. 4.56).
    fn left_fan_state(&self, speed: f64) -> (f64, f64, f64) {
        let (gm1, gp1) = (self.gamma - 1.0, self.gamma + 1.0);
        let (rho_l, u_l, p_l) = (self.left.rho, self.left.velocity, self.left.pressure);
        let c_l = self.left.sound_speed(self.gamma);
        let factor = 2.0 / gp1 + gm1 / (gp1 * c_l) * (u_l - speed);
        let rho = rho_l * factor.powf(2.0 / gm1);
        let u = 2.0 / gp1 * (c_l + 0.5 * gm1 * u_l + speed);
        let p = p_l * factor.powf(2.0 * self.gamma / gm1);
        (rho, u, p)
    }

    /// State inside a right rarefaction fan at the similarity coordinate
    /// `speed = (x − x₀)/t` (Toro eq. 4.63).
    fn right_fan_state(&self, speed: f64) -> (f64, f64, f64) {
        let (gm1, gp1) = (self.gamma - 1.0, self.gamma + 1.0);
        let (rho_r, u_r, p_r) = (self.right.rho, self.right.velocity, self.right.pressure);
        let c_r = self.right.sound_speed(self.gamma);
        let factor = 2.0 / gp1 - gm1 / (gp1 * c_r) * (u_r - speed);
        let rho = rho_r * factor.powf(2.0 / gm1);
        let u = 2.0 / gp1 * (-c_r + 0.5 * gm1 * u_r + speed);
        let p = p_r * factor.powf(2.0 * self.gamma / gm1);
        (rho, u, p)
    }

    /// Sample the solution at the similarity coordinate `speed = (x − x₀)/t`
    /// for a vacuum between two rarefaction waves (Toro section 4.6).
    fn sample_vacuum(&self, speed: f64) -> (f64, f64, f64) {
        let gm1 = self.gamma - 1.0;
        let (left, right) = (&self.left, &self.right);
        let c_l = left.sound_speed(self.gamma);
        let c_r = right.sound_speed(self.gamma);
        let head_left = left.velocity - c_l;
        let tail_left = left.velocity + 2.0 * c_l / gm1;
        let tail_right = right.velocity - 2.0 * c_r / gm1;
        let head_right = right.velocity + c_r;
        if speed <= head_left {
            (left.rho, left.velocity, left.pressure)
        } else if speed < tail_left {
            self.left_fan_state(speed)
        } else if speed <= tail_right {
            (0.0, 0.0, 0.0)
        } else if speed < head_right {
            self.right_fan_state(speed)
        } else {
            (right.rho, right.velocity, right.pressure)
        }
    }

    /// Sample the solution at the similarity coordinate `speed = (x − x₀)/t`
    /// when the star region has pressure `p_star` and velocity `u_star`.
    fn sample_star_region(&self, speed: f64, p_star: f64, u_star: f64) -> (f64, f64, f64) {
        let g = self.gamma;
        let (gm1, gp1) = (g - 1.0, g + 1.0);
        let (left, right) = (&self.left, &self.right);
        let exponent = gm1 / (2.0 * g);
        if speed <= u_star {
            // Left of the contact discontinuity.
            let (rho_l, u_l, p_l) = (left.rho, left.velocity, left.pressure);
            let c_l = left.sound_speed(g);
            let ratio = p_star / p_l;
            if p_star > p_l {
                // Left shock.
                let shock_speed = u_l - c_l * (gp1 / (2.0 * g) * ratio + exponent).sqrt();
                if speed <= shock_speed {
                    return (rho_l, u_l, p_l);
                }
                let rho_star = rho_l * (ratio + gm1 / gp1) / (gm1 / gp1 * ratio + 1.0);
                return (rho_star, u_star, p_star);
            }
            // Left rarefaction: head, fan, tail.
            let head = u_l - c_l;
            let tail = u_star - c_l * ratio.powf(exponent);
            if speed <= head {
                (rho_l, u_l, p_l)
            } else if speed >= tail {
                (rho_l * ratio.powf(1.0 / g), u_star, p_star)
            } else {
                self.left_fan_state(speed)
            }
        } else {
            // Right of the contact discontinuity.
            let (rho_r, u_r, p_r) = (right.rho, right.velocity, right.pressure);
            let c_r = right.sound_speed(g);
            let ratio = p_star / p_r;
            if p_star > p_r {
                // Right shock.
                let shock_speed = u_r + c_r * (gp1 / (2.0 * g) * ratio + exponent).sqrt();
                if speed >= shock_speed {
                    return (rho_r, u_r, p_r);
                }
                let rho_star = rho_r * (ratio + gm1 / gp1) / (gm1 / gp1 * ratio + 1.0);
                return (rho_star, u_star, p_star);
            }
            // Right rarefaction: head, fan, tail.
            let head = u_r + c_r;
            let tail = u_star + c_r * ratio.powf(exponent);
            if speed >= head {
                (rho_r, u_r, p_r)
            } else if speed <= tail {
                (rho_r * ratio.powf(1.0 / g), u_star, p_star)
            } else {
                self.right_fan_state(speed)
            }
        }
    }

    /// Compute the exact solution at time `t` and position `x`.
    ///
    /// Returns `(density, velocity, pressure)`.
    ///
    /// The solution is the exact solution of the Riemann problem with the
    /// stored left and right states (including non-zero initial velocities),
    /// the adiabatic index `gamma` and the discontinuity at `diaphragm`. For
    /// `t ≤ 0` the initial states are returned. When the two rarefaction waves
    /// generate a vacuum, the vacuum region has zero density and pressure and
    /// its (undefined) velocity is reported as zero. For unphysical input
    /// (non-positive or non-finite densities or pressures, or `γ ≤ 1`) every
    /// component is `NaN`.
    pub fn exact_solution(&self, x: f64, t: f64) -> (f64, f64, f64) {
        if t <= 0.0 {
            let side = if x < self.diaphragm {
                &self.left
            } else {
                &self.right
            };
            return (side.rho, side.velocity, side.pressure);
        }
        let speed = (x - self.diaphragm) / t;
        match self.star_region_solution() {
            Some(StarRegion::Contact { pressure, velocity }) => {
                self.sample_star_region(speed, pressure, velocity)
            }
            Some(StarRegion::Vacuum) => self.sample_vacuum(speed),
            None => (f64::NAN, f64::NAN, f64::NAN),
        }
    }

    /// Generate initial SPH particles for the shock tube.
    ///
    /// Returns `(positions, densities, pressures, velocities)` as 1-D arrays.
    pub fn generate_particles(
        &self,
        n_particles: usize,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let dx = self.length / n_particles as f64;
        let mut positions = Vec::with_capacity(n_particles);
        let mut densities = Vec::with_capacity(n_particles);
        let mut pressures = Vec::with_capacity(n_particles);
        let mut velocities = Vec::with_capacity(n_particles);

        for i in 0..n_particles {
            let x = (i as f64 + 0.5) * dx;
            positions.push(x);
            if x < self.diaphragm {
                densities.push(self.left.rho);
                pressures.push(self.left.pressure);
                velocities.push(self.left.velocity);
            } else {
                densities.push(self.right.rho);
                pressures.push(self.right.pressure);
                velocities.push(self.right.velocity);
            }
        }

        (positions, densities, pressures, velocities)
    }
}

#[cfg(test)]
mod tests;
