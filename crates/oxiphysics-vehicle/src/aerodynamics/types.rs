//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::functions::AIR_DENSITY_SEA_LEVEL;
#[allow(unused_imports)]
use super::functions::*;

/// Feedback controller that adjusts front/rear wing levels to achieve a
/// target aerodynamic balance.
#[derive(Debug, Clone)]
pub struct AeroBalanceController {
    /// Target front downforce fraction (0–1).
    pub target_front_fraction: f64,
    /// Proportional gain for wing level correction.
    pub gain: f64,
    /// Maximum allowed wing level change per update step.
    pub max_delta_level: f64,
}
impl AeroBalanceController {
    /// Compute the required change in front wing level to reduce balance error.
    ///
    /// Positive return → increase front wing level (more front downforce).
    pub fn front_wing_correction(&self, current_balance: &AeroBalance) -> f64 {
        let error = self.target_front_fraction - current_balance.front_fraction();
        let delta = self.gain * error;
        delta.clamp(-self.max_delta_level, self.max_delta_level)
    }
    /// Apply the correction to a `DownforcePackage`.
    ///
    /// Clips the result to `[1, max_level]`.
    pub fn apply_correction(&self, pkg: &mut DownforcePackage, balance: &AeroBalance) {
        let delta = self.front_wing_correction(balance);
        let new_level = (pkg.front_level as f64 + delta).round() as i32;
        pkg.front_level = new_level.clamp(1, pkg.max_level as i32) as u8;
    }
}
/// Models gradual and abrupt wing stall behaviour.
///
/// Distinguishes between leading-edge stall (abrupt) and trailing-edge stall
/// (gradual) based on the leading-edge radius.
#[derive(Debug, Clone)]
pub struct WingStallModel {
    /// Lift-curve slope below stall (per radian).
    pub lift_slope: f64,
    /// Critical angle of attack for stall onset (degrees).
    pub alpha_stall: f64,
    /// Peak lift coefficient at stall.
    pub cl_max: f64,
    /// Post-stall lift coefficient (fractional recovery relative to `cl_max`).
    pub cl_post_stall_fraction: f64,
    /// `true` for leading-edge (abrupt) stall, `false` for trailing-edge (gradual).
    pub leading_edge_stall: bool,
}
impl WingStallModel {
    /// Evaluate the lift coefficient at the given angle of attack (degrees).
    pub fn lift_coefficient(&self, alpha_deg: f64) -> f64 {
        let alpha_rad = alpha_deg.to_radians();
        let _alpha_stall_rad = self.alpha_stall.to_radians();
        if alpha_deg.abs() <= self.alpha_stall {
            self.lift_slope * alpha_rad
        } else if self.leading_edge_stall {
            let sign = alpha_deg.signum();
            let cl_post = self.cl_max * self.cl_post_stall_fraction;
            sign * cl_post
        } else {
            let sign = alpha_deg.signum();
            let excess = (alpha_deg.abs() - self.alpha_stall).to_radians();
            let decay = (excess * 0.5).cos().powi(2);
            sign * self.cl_max
                * (self.cl_post_stall_fraction + (1.0 - self.cl_post_stall_fraction) * decay)
        }
        .clamp(-self.cl_max * 1.5, self.cl_max * 1.5)
    }
    /// Drag coefficient accounting for separation drag after stall.
    ///
    /// Below stall: use induced drag.  Above stall: add separation drag.
    pub fn drag_coefficient(&self, alpha_deg: f64, span: f64, chord: f64) -> f64 {
        let cl = self.lift_coefficient(alpha_deg);
        let ar = span / chord;
        let cd_induced = cl * cl / (PI * ar * 0.9);
        let cd0 = 0.01;
        if alpha_deg.abs() > self.alpha_stall {
            let excess_rad = (alpha_deg.abs() - self.alpha_stall).to_radians();
            let cd_separation = 2.0 * excess_rad.sin().powi(2);
            cd0 + cd_induced + cd_separation
        } else {
            cd0 + cd_induced
        }
    }
    /// True if the wing is currently stalled.
    pub fn is_stalled(&self, alpha_deg: f64) -> bool {
        alpha_deg.abs() > self.alpha_stall
    }
}
/// Wing configuration package describing front and rear downforce levels.
///
/// Wing levels are integers (e.g. 1–7 in Formula 1 parlance) that map to
/// effective lift coefficients via a linear interpolation table.
#[derive(Debug, Clone)]
pub struct DownforcePackage {
    /// Front wing level (1 = minimum, `max_level` = maximum downforce).
    pub front_level: u8,
    /// Rear wing level (1 = minimum, `max_level` = maximum downforce).
    pub rear_level: u8,
    /// Maximum allowed wing level.
    pub max_level: u8,
    /// ClA at minimum wing level (most negative = most downforce at min).
    pub cl_min: f64,
    /// ClA at maximum wing level.
    pub cl_max: f64,
}
impl DownforcePackage {
    /// Effective front-wing ClA at the current `front_level`.
    pub fn front_cla(&self) -> f64 {
        self.level_to_cla(self.front_level)
    }
    /// Effective rear-wing ClA at the current `rear_level`.
    pub fn rear_cla(&self) -> f64 {
        self.level_to_cla(self.rear_level)
    }
    /// Total (front + rear) ClA.
    pub fn total_cla(&self) -> f64 {
        self.front_cla() + self.rear_cla()
    }
    fn level_to_cla(&self, level: u8) -> f64 {
        if self.max_level <= 1 {
            return self.cl_min;
        }
        let t = (level as f64 - 1.0) / (self.max_level as f64 - 1.0);
        let t = t.clamp(0.0, 1.0);
        self.cl_min + t * (self.cl_max - self.cl_min)
    }
}
/// Tabulated aerodynamic database for drag coefficient lookup.
///
/// Stores a 2-D table of `Cd` values indexed by angle-of-attack (alpha) and
/// Reynolds number (Re).  Bilinear interpolation is used for queries between
/// table entries.
#[derive(Debug, Clone)]
pub struct AeroDatabase {
    /// Alpha breakpoints (degrees), sorted ascending.
    pub alpha_table: Vec<f64>,
    /// Reynolds number breakpoints (–), sorted ascending.
    pub re_table: Vec<f64>,
    /// Cd values: row = alpha index, column = Re index.
    /// Stored as a flat vector: `cd[i * re_table.len() + j]`.
    pub cd_table: Vec<f64>,
}
impl AeroDatabase {
    /// Create a new `AeroDatabase`.
    ///
    /// # Panics
    /// Panics in debug mode if `cd_table.len() != alpha_table.len() * re_table.len()`.
    pub fn new(alpha_table: Vec<f64>, re_table: Vec<f64>, cd_table: Vec<f64>) -> Self {
        debug_assert_eq!(
            cd_table.len(),
            alpha_table.len() * re_table.len(),
            "cd_table size mismatch"
        );
        Self {
            alpha_table,
            re_table,
            cd_table,
        }
    }
    /// Build a simple single-Re database (Re dimension = 1) from alpha/Cd pairs.
    pub fn from_alpha_cd(alphas: Vec<f64>, cds: Vec<f64>) -> Self {
        debug_assert_eq!(alphas.len(), cds.len());
        Self::new(alphas, vec![1.0e6], cds)
    }
    /// Bilinear interpolation of Cd at given alpha (degrees) and Reynolds number.
    ///
    /// Clamps to the table boundaries when the query is outside the range.
    pub fn interpolate_cd(&self, alpha_deg: f64, re: f64) -> f64 {
        if self.alpha_table.is_empty() || self.re_table.is_empty() {
            return 0.0;
        }
        let (ai, at) = interp_index(&self.alpha_table, alpha_deg);
        let (ri, rt) = interp_index(&self.re_table, re);
        let na = self.alpha_table.len();
        let nr = self.re_table.len();
        let i0 = ai.min(na - 1);
        let i1 = (ai + 1).min(na - 1);
        let j0 = ri.min(nr - 1);
        let j1 = (ri + 1).min(nr - 1);
        let c00 = self.cd_table[i0 * nr + j0];
        let c01 = self.cd_table[i0 * nr + j1];
        let c10 = self.cd_table[i1 * nr + j0];
        let c11 = self.cd_table[i1 * nr + j1];
        let c0 = c00 + at * (c10 - c00);
        let c1 = c01 + at * (c11 - c01);
        c0 + rt * (c1 - c0)
    }
}
/// A radial-basis-function (RBF) CFD surrogate model.
///
/// Surrogate models replace expensive CFD runs at run-time: a small set of
/// `n_centers` sample points (from offline CFD) are stored, and predictions at
/// new query points are computed by weighted interpolation.
///
/// Only f64 arrays are used — no nalgebra dependency.
pub struct CfdSurrogateModel {
    /// Input feature vectors (each of length `feature_dim`).
    pub centers: Vec<Vec<f64>>,
    /// Scalar outputs at each centre (e.g. drag or downforce coefficient).
    pub values: Vec<f64>,
    /// RBF length scale `ε`.
    pub epsilon: f64,
}
impl CfdSurrogateModel {
    /// Create a surrogate model from sampled (feature, value) pairs.
    pub fn new(centers: Vec<Vec<f64>>, values: Vec<f64>, epsilon: f64) -> Self {
        Self {
            centers,
            values,
            epsilon,
        }
    }
    /// Predict the output at a query feature vector using Gaussian RBF
    /// interpolation.
    ///
    /// `φ(r) = exp(-(ε·r)²)`
    pub fn predict(&self, query: &[f64]) -> f64 {
        let mut weighted_sum = 0.0;
        let mut weight_sum = 0.0;
        for (center, &val) in self.centers.iter().zip(self.values.iter()) {
            let r_sq: f64 = center
                .iter()
                .zip(query.iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum();
            let r = r_sq.sqrt();
            let w = (-self.epsilon * self.epsilon * r * r).exp();
            weighted_sum += w * val;
            weight_sum += w;
        }
        if weight_sum > 1e-30 {
            weighted_sum / weight_sum
        } else {
            self.nearest_neighbour(query)
        }
    }
    /// Return the value of the nearest centre to `query`.
    fn nearest_neighbour(&self, query: &[f64]) -> f64 {
        self.centers
            .iter()
            .zip(self.values.iter())
            .map(|(center, &val)| {
                let r_sq: f64 = center
                    .iter()
                    .zip(query.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum();
                (r_sq, val)
            })
            .min_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(_, v)| v)
            .unwrap_or(0.0)
    }
    /// Number of centres.
    pub fn n_centers(&self) -> usize {
        self.centers.len()
    }
    /// Feature dimension.
    pub fn feature_dim(&self) -> usize {
        self.centers.first().map(|c| c.len()).unwrap_or(0)
    }
}
/// Lookup table of Cd and Cl as a function of vehicle speed.
pub struct SpeedAeroMap {
    /// Speed sample points in m/s (must be in ascending order).
    pub velocities: Vec<f64>,
    /// Drag coefficients at each sample speed.
    pub cds: Vec<f64>,
    /// Lift coefficients at each sample speed.
    pub cls: Vec<f64>,
}
impl SpeedAeroMap {
    /// Linearly interpolate the drag coefficient at a given speed.
    pub fn interpolate_cd(&self, v: f64) -> f64 {
        interpolate_table(&self.velocities, &self.cds, v)
    }
    /// Linearly interpolate the lift coefficient at a given speed.
    pub fn interpolate_cl(&self, v: f64) -> f64 {
        interpolate_table(&self.velocities, &self.cls, v)
    }
    /// Construct a constant (speed-independent) aero map.
    pub fn from_constant(cd: f64, cl: f64) -> Self {
        Self {
            velocities: vec![0.0, 100.0],
            cds: vec![cd, cd],
            cls: vec![cl, cl],
        }
    }
}
/// Models the aerodynamic wake interaction between a leading and a trailing
/// vehicle (e.g. dirty-air effect in motorsport).
#[derive(Debug, Clone)]
pub struct WakeInteraction {
    /// Downstream distance from the leading vehicle's rear to the trailing
    /// vehicle's front (m).
    pub following_distance: f64,
    /// Lateral offset between vehicle centre-lines (m).
    pub lateral_offset: f64,
    /// Wake decay length constant (m).  Typical value: 1.5–3.0 m.
    pub decay_length: f64,
    /// Reference frontal area of the leading vehicle (m²).
    pub leading_frontal_area: f64,
}
impl WakeInteraction {
    /// Wake velocity deficit at the trailing vehicle's nose.
    ///
    /// `ΔV/V∞ = (A_lead / (2π·L²)) · exp(-d / Λ)`
    ///
    /// where `d` is the 3-D distance and `Λ` is the decay length.
    pub fn velocity_deficit(&self, freestream_speed: f64) -> f64 {
        let d_ax = self.following_distance;
        let d_lat = self.lateral_offset;
        let d = (d_ax * d_ax + d_lat * d_lat).sqrt().max(0.01);
        let exp_decay = (-d / self.decay_length).exp();
        let deficit_fraction = self.leading_frontal_area / (2.0 * PI * d * d) * exp_decay;
        freestream_speed * deficit_fraction.min(0.5)
    }
    /// Effective dynamic pressure seen by the trailing vehicle.
    ///
    /// `q_eff = 0.5 · ρ · (V∞ − ΔV)²`
    pub fn effective_dynamic_pressure(&self, freestream_speed: f64, air_density: f64) -> f64 {
        let deficit = self.velocity_deficit(freestream_speed);
        let v_eff = (freestream_speed - deficit).max(0.0);
        0.5 * air_density * v_eff * v_eff
    }
    /// Downforce loss fraction on the trailing vehicle.
    ///
    /// Downforce ∝ v²; so the fractional loss = 1 − (v_eff / v_inf)².
    pub fn downforce_loss_fraction(&self, freestream_speed: f64) -> f64 {
        if freestream_speed < 1e-12 {
            return 0.0;
        }
        let deficit = self.velocity_deficit(freestream_speed);
        let v_eff = (freestream_speed - deficit).max(0.0);
        let ratio = v_eff / freestream_speed;
        (1.0 - ratio * ratio).max(0.0)
    }
}
/// A 2-D lookup table for downforce coefficient as a function of speed and
/// ride height.
///
/// Intended for use as a CFD surrogate model: the coefficients were computed
/// offline and stored in a grid for fast run-time evaluation.
pub struct DownforceMap {
    /// Speed sample points in m/s (ascending).
    pub speeds: Vec<f64>,
    /// Ride height sample points in m (ascending).
    pub ride_heights: Vec<f64>,
    /// ClA values at each `(speed, ride_height)` combination.
    ///
    /// Row-major: `cla[i * n_rh + j]` is for `speeds[i]`, `ride_heights[j]`.
    pub cla: Vec<f64>,
}
impl DownforceMap {
    /// Bilinear interpolation of ClA at `(speed, ride_height)`.
    pub fn interpolate(&self, speed: f64, ride_height: f64) -> f64 {
        let ns = self.speeds.len();
        let nr = self.ride_heights.len();
        if ns == 0 || nr == 0 {
            return 0.0;
        }
        let (is, ts) = bracket_index(&self.speeds, speed);
        let (ir, tr) = bracket_index(&self.ride_heights, ride_height);
        let is1 = (is + 1).min(ns - 1);
        let ir1 = (ir + 1).min(nr - 1);
        let v00 = self.cla[is * nr + ir];
        let v10 = self.cla[is1 * nr + ir];
        let v01 = self.cla[is * nr + ir1];
        let v11 = self.cla[is1 * nr + ir1];
        let va = v00 * (1.0 - ts) + v10 * ts;
        let vb = v01 * (1.0 - ts) + v11 * ts;
        va * (1.0 - tr) + vb * tr
    }
    /// Build a constant downforce map (independent of speed and ride height).
    pub fn from_constant(cla: f64) -> Self {
        Self {
            speeds: vec![0.0, 100.0],
            ride_heights: vec![0.0, 0.5],
            cla: vec![cla; 4],
        }
    }
    /// Build a simple parametric downforce map with Venturi ground-effect
    /// scaling.
    ///
    /// `cla(v, h) = cla_ref * (v / v_ref)^0.05 * (1 + k * h_ref / h)`
    ///
    /// where `k` controls the ground-effect sensitivity.
    pub fn parametric(
        cla_ref: f64,
        v_ref: f64,
        h_ref: f64,
        k: f64,
        speeds: Vec<f64>,
        ride_heights: Vec<f64>,
    ) -> Self {
        let nr = ride_heights.len();
        let mut cla = Vec::with_capacity(speeds.len() * nr);
        for &v in &speeds {
            let v_factor = if v_ref > 0.0 {
                (v / v_ref).powf(0.05)
            } else {
                1.0
            };
            for &h in &ride_heights {
                let h_factor = if h > 0.0 {
                    1.0 + k * h_ref / h
                } else {
                    1.0 + k
                };
                cla.push(cla_ref * v_factor * h_factor);
            }
        }
        Self {
            speeds,
            ride_heights,
            cla,
        }
    }
}
/// Underbody diffuser model for a racing car.
///
/// The diffuser recovers pressure after the Venturi section, generating
/// additional downforce and creating a suction region under the car.
#[derive(Debug, Clone)]
pub struct DiffuserModel {
    /// Diffuser inlet height (m).
    pub inlet_height: f64,
    /// Diffuser exit height (m).
    pub exit_height: f64,
    /// Diffuser length (m).
    pub length: f64,
    /// Diffuser width (m).
    pub width: f64,
    /// Pressure recovery efficiency (0–1).
    pub efficiency: f64,
}
impl DiffuserModel {
    /// Area ratio of the diffuser (exit / inlet).
    pub fn area_ratio(&self) -> f64 {
        if self.inlet_height < 1e-12 {
            return 1.0;
        }
        (self.exit_height * self.width) / (self.inlet_height * self.width)
    }
    /// Theoretical pressure rise coefficient `Cp = 1 − AR^{-2}`.
    pub fn ideal_cp(&self) -> f64 {
        let ar = self.area_ratio();
        if ar < 1e-10 {
            return 0.0;
        }
        (1.0 - 1.0 / (ar * ar)).max(0.0)
    }
    /// Effective pressure rise coefficient (corrected for losses).
    pub fn effective_cp(&self) -> f64 {
        self.ideal_cp() * self.efficiency
    }
    /// Downforce generated by the diffuser suction.
    ///
    /// `F = 0.5 · ρ · V² · Cp · A_diffuser`
    pub fn downforce(&self, speed: f64, air_density: f64) -> f64 {
        let q = 0.5 * air_density * speed * speed;
        let a_diff = self.length * self.width;
        q * self.effective_cp() * a_diff
    }
    /// Diffuser expansion angle in degrees.
    pub fn expansion_angle(&self) -> f64 {
        if self.length < 1e-12 {
            return 0.0;
        }
        let dh = self.exit_height - self.inlet_height;
        (dh / self.length).atan().to_degrees()
    }
}
/// Steady + gusty wind field model.
///
/// The wind is specified as a mean velocity vector `[vx, vy, vz]` (m/s in
/// world frame) plus a gust amplitude.  A simple deterministic pseudo-random
/// gust is produced by sampling a sine wave at the query position.
#[derive(Debug, Clone)]
pub struct WindModel {
    /// Mean (steady) wind velocity `[vx, vy, vz]` in m/s (world frame).
    pub mean_velocity: [f64; 3],
    /// Peak gust amplitude in m/s.
    pub gust_amplitude: f64,
    /// Spatial gust wavelength in metres.
    pub gust_wavelength: f64,
}
impl WindModel {
    /// Wind velocity at `position` in world frame (m/s).
    ///
    /// The gust component is a sine wave whose phase depends on the
    /// arc-length along the mean-wind direction.
    pub fn velocity_at(&self, position: [f64; 3]) -> [f64; 3] {
        let mean_len = vec3_len(self.mean_velocity);
        let phase = if mean_len > 1e-12 {
            let dot = self.mean_velocity[0] * position[0]
                + self.mean_velocity[1] * position[1]
                + self.mean_velocity[2] * position[2];
            dot / mean_len
        } else {
            position[0]
        };
        let gust = if self.gust_wavelength > 1e-12 {
            self.gust_amplitude * (2.0 * PI * phase / self.gust_wavelength).sin()
        } else {
            0.0
        };
        if mean_len > 1e-12 {
            let nx = self.mean_velocity[0] / mean_len;
            let ny = self.mean_velocity[1] / mean_len;
            let nz = self.mean_velocity[2] / mean_len;
            [
                self.mean_velocity[0] + gust * nx,
                self.mean_velocity[1] + gust * ny,
                self.mean_velocity[2] + gust * nz,
            ]
        } else {
            [gust, 0.0, 0.0]
        }
    }
    /// Relative wind velocity seen by a vehicle moving with `vehicle_velocity`
    /// at `position` (both in world frame, m/s).
    pub fn relative_wind(&self, position: [f64; 3], vehicle_velocity: [f64; 3]) -> [f64; 3] {
        let w = self.velocity_at(position);
        [
            w[0] - vehicle_velocity[0],
            w[1] - vehicle_velocity[1],
            w[2] - vehicle_velocity[2],
        ]
    }
}
/// A finite lifting surface (wing, diffuser, splitter).
pub struct AeroWing {
    /// Wing span b in m.
    pub span: f64,
    /// Chord length c in m.
    pub chord: f64,
    /// Geometric angle of attack α in degrees.
    pub angle_of_attack: f64,
    /// Aerofoil profile.
    pub profile: WingProfile,
}
impl AeroWing {
    /// Wing aspect ratio AR = b² / (b·c) = b / c.
    pub fn aspect_ratio(&self) -> f64 {
        self.span / self.chord
    }
    /// Planform area in m².
    fn planform_area(&self) -> f64 {
        self.span * self.chord
    }
    /// Section lift coefficient based on NACA thin-airfoil theory.
    ///
    /// Below stall (|α| ≤ 15°): `Cl = 2π·α_rad`
    /// Above stall: `Cl = Cl_max · cos(α - 15°)` where `Cl_max = 2π·(15°)`
    pub fn lift_coefficient(&self) -> f64 {
        let alpha_l0 = self.profile.zero_lift_angle();
        let alpha_eff = self.angle_of_attack - alpha_l0;
        let stall_deg = 15.0_f64;
        let alpha_rad = alpha_eff.to_radians();
        if alpha_eff.abs() <= stall_deg {
            2.0 * PI * alpha_rad
        } else {
            let cl_max = 2.0 * PI * stall_deg.to_radians();
            let delta = (alpha_eff.abs() - stall_deg).to_radians();
            let sign = if alpha_eff >= 0.0 { 1.0 } else { -1.0 };
            sign * cl_max * delta.cos()
        }
    }
    /// Section drag coefficient (profile + induced drag).
    ///
    /// `Cd = Cd0 + Cl² / (π · AR · e)` where e = 0.9 (Oswald efficiency)
    pub fn drag_coefficient(&self) -> f64 {
        let cl = self.lift_coefficient();
        let cd0 = 0.01;
        let e = 0.9;
        let ar = self.aspect_ratio();
        cd0 + cl * cl / (PI * ar * e)
    }
    /// Lift force in N.
    ///
    /// `F_l = 0.5 · ρ · v² · Cl · S`
    pub fn lift_force(&self, velocity: f64, air_density: f64) -> f64 {
        0.5 * air_density * velocity * velocity * self.lift_coefficient() * self.planform_area()
    }
    /// Drag force in N.
    ///
    /// `F_d = 0.5 · ρ · v² · Cd · S`
    pub fn drag_force(&self, velocity: f64, air_density: f64) -> f64 {
        0.5 * air_density * velocity * velocity * self.drag_coefficient() * self.planform_area()
    }
    /// Aerodynamic efficiency ratio Cl / Cd.
    pub fn efficiency_ratio(&self) -> f64 {
        let cd = self.drag_coefficient();
        if cd.abs() < f64::EPSILON {
            return 0.0;
        }
        self.lift_coefficient() / cd
    }
}
/// Simplified virtual wind tunnel for aerodynamic measurements.
pub struct WindTunnel {
    /// Free-stream flow velocity in m/s.
    pub test_velocity: f64,
    /// Air density in kg/m³.
    pub air_density: f64,
}
impl WindTunnel {
    /// Measure drag and lift forces on a body.
    ///
    /// Returns `(drag, lift)` in N.
    pub fn measure_forces(&self, body: &AerodynamicBody) -> (f64, f64) {
        let drag = body.drag_force(self.test_velocity, self.air_density);
        let lift = body.lift_force(self.test_velocity, self.air_density);
        (drag, lift)
    }
    /// Measure combined body + wing drag and lift.
    ///
    /// Returns `(drag, lift)` in N.
    pub fn measure_with_wing(&self, body: &AerodynamicBody, wing: &AeroWing) -> (f64, f64) {
        let (body_drag, body_lift) = self.measure_forces(body);
        let wing_drag = wing.drag_force(self.test_velocity, self.air_density);
        let wing_lift = wing.lift_force(self.test_velocity, self.air_density);
        (body_drag + wing_drag, body_lift + wing_lift)
    }
    /// Convert drag counts to Cd (1 drag count = 0.001 Cd).
    pub fn drag_count_to_cd(drag_count: f64) -> f64 {
        drag_count * 0.001
    }
    /// Convert Cd to drag counts.
    pub fn cd_to_drag_count(cd: f64) -> f64 {
        cd / 0.001
    }
}
/// Aerodynamic properties of a vehicle body.
pub struct AerodynamicBody {
    /// Frontal (cross-sectional) area A in m².
    pub frontal_area: f64,
    /// Drag coefficient Cd.
    pub drag_coefficient: f64,
    /// Lift coefficient Cl (negative = downforce).
    pub lift_coefficient: f64,
    /// Side-force coefficient Cs.
    pub side_force_coefficient: f64,
    /// Vehicle length in m.
    pub length: f64,
    /// Wheelbase in m.
    pub wheelbase: f64,
}
impl AerodynamicBody {
    /// Aerodynamic drag force in N.
    ///
    /// `F_d = 0.5 · ρ · v² · Cd · A`
    pub fn drag_force(&self, velocity: f64, air_density: f64) -> f64 {
        0.5 * air_density * velocity * velocity * self.drag_coefficient * self.frontal_area
    }
    /// Aerodynamic lift force in N (negative = downforce).
    ///
    /// `F_l = 0.5 · ρ · v² · Cl · A`
    pub fn lift_force(&self, velocity: f64, air_density: f64) -> f64 {
        0.5 * air_density * velocity * velocity * self.lift_coefficient * self.frontal_area
    }
    /// Reynolds number.
    ///
    /// `Re = v · L / ν`
    pub fn reynolds_number(&self, velocity: f64, kinematic_viscosity: f64) -> f64 {
        velocity * self.length / kinematic_viscosity
    }
}
/// Decomposed aerodynamic force vector.
#[derive(Debug, Clone, Default)]
pub struct AeroForces {
    /// Drag force vector `[Fx, Fy, Fz]` in N (opposite to velocity).
    pub drag: [f64; 3],
    /// Lift force vector `[Fx, Fy, Fz]` in N (along +z axis).
    pub lift: [f64; 3],
    /// Side force vector `[Fx, Fy, Fz]` in N (along ±y axis).
    pub side: [f64; 3],
}
impl AeroForces {
    /// Vector sum of all force components.
    pub fn total(&self) -> [f64; 3] {
        [
            self.drag[0] + self.lift[0] + self.side[0],
            self.drag[1] + self.lift[1] + self.side[1],
            self.drag[2] + self.lift[2] + self.side[2],
        ]
    }
}
/// Aerofoil cross-section profile.
pub enum WingProfile {
    /// Symmetric NACA 00xx profile.
    Symmetric {
        /// Maximum thickness as a fraction of chord.
        max_thickness: f64,
    },
    /// Cambered NACA 4-digit profile.
    Cambered {
        /// Maximum camber as a fraction of chord.
        camber: f64,
        /// Maximum thickness as a fraction of chord.
        max_thickness: f64,
    },
}
impl WingProfile {
    /// Zero-lift angle of attack in degrees (0 for symmetric profiles).
    pub fn zero_lift_angle(&self) -> f64 {
        match self {
            WingProfile::Symmetric { .. } => 0.0,
            WingProfile::Cambered { camber, .. } => -2.0 * camber * 180.0 / PI,
        }
    }
}
/// Component breakdown of the total aerodynamic drag.
pub struct DragBreakdown {
    /// Pressure (form) drag contribution.
    pub pressure_drag: f64,
    /// Skin-friction drag contribution.
    pub friction_drag: f64,
    /// Induced drag (lift-induced vortex drag).
    pub induced_drag: f64,
    /// Interference drag (junctions, protrusions).
    pub interference_drag: f64,
}
impl DragBreakdown {
    /// Sum of all drag components.
    pub fn total(&self) -> f64 {
        self.pressure_drag + self.friction_drag + self.induced_drag + self.interference_drag
    }
    /// Compute a simplified breakdown for a body at given speed and Reynolds number.
    ///
    /// - Skin-friction coefficient via turbulent flat-plate formula: `Cf = 0.074 / Re^0.2`
    /// - Pressure drag is estimated at 60 % of the total.
    /// - Remaining budget is friction drag; induced and interference are small corrections.
    pub fn compute(body: &AerodynamicBody, velocity: f64, re: f64) -> Self {
        let q = 0.5 * AIR_DENSITY_SEA_LEVEL * velocity * velocity;
        let total_drag = body.drag_coefficient * body.frontal_area * q;
        let cf = if re > 0.0 { 0.074 / re.powf(0.2) } else { 0.0 };
        let friction_drag = cf * body.frontal_area * q;
        let pressure_drag = 0.6 * total_drag;
        let remaining = (total_drag - pressure_drag - friction_drag).max(0.0);
        let induced_drag = 0.7 * remaining;
        let interference_drag = 0.3 * remaining;
        Self {
            pressure_drag,
            friction_drag,
            induced_drag,
            interference_drag,
        }
    }
}
/// Front/rear aerodynamic load (downforce) distribution.
pub struct AeroBalance {
    /// Downforce on the front axle in N (positive = pushes down).
    pub front_downforce: f64,
    /// Downforce on the rear axle in N (positive = pushes down).
    pub rear_downforce: f64,
}
impl AeroBalance {
    /// Total downforce in N.
    pub fn total(&self) -> f64 {
        self.front_downforce + self.rear_downforce
    }
    /// Fraction of total downforce carried by the front axle.
    pub fn front_fraction(&self) -> f64 {
        let tot = self.total();
        if tot.abs() < f64::EPSILON {
            return 0.5;
        }
        self.front_downforce / tot
    }
    /// Distance of the aerodynamic balance point from the front axle in m.
    pub fn balance_point(&self, wheelbase: f64) -> f64 {
        (1.0 - self.front_fraction()) * wheelbase
    }
    /// Compute front/rear downforce split.
    ///
    /// # Arguments
    /// * `cog_x` – centre-of-gravity x-position measured from the front axle in m
    pub fn compute(
        body: &AerodynamicBody,
        velocity: f64,
        air_density: f64,
        cog_x: f64,
        wheelbase: f64,
    ) -> Self {
        let downforce = -body.lift_force(velocity, air_density);
        let cop_x = 0.5 * wheelbase;
        let rear_fraction = if wheelbase > 0.0 {
            ((cop_x - cog_x) / wheelbase).clamp(0.0, 1.0)
        } else {
            0.5
        };
        let front_fraction = 1.0 - rear_fraction;
        Self {
            front_downforce: front_fraction * downforce,
            rear_downforce: rear_fraction * downforce,
        }
    }
}
/// Aerodynamic coefficient lookup table indexed by angle of attack `alpha`
/// and sideslip angle `beta` (both in degrees).
///
/// Coefficients `(Cd, Cl, Cs)` are stored in a flat row-major grid and
/// retrieved by bilinear interpolation.
pub struct AeroMap {
    /// Alpha (angle of attack) sample points in degrees (ascending).
    pub alphas: Vec<f64>,
    /// Beta (sideslip angle) sample points in degrees (ascending).
    pub betas: Vec<f64>,
    /// Drag coefficient grid, shape `[n_alpha × n_beta]`, row-major.
    pub cd: Vec<f64>,
    /// Lift coefficient grid, shape `[n_alpha × n_beta]`, row-major.
    pub cl: Vec<f64>,
    /// Side-force coefficient grid, shape `[n_alpha × n_beta]`, row-major.
    pub cs: Vec<f64>,
}
impl AeroMap {
    /// Interpolate `(Cd, Cl, Cs)` at the given `(alpha, beta)` angles.
    ///
    /// Returns `(Cd, Cl, Cs)`.  Values are clamped to the table boundaries.
    pub fn interpolate(&self, alpha: f64, beta: f64) -> (f64, f64, f64) {
        let cd = self.bilinear(&self.cd, alpha, beta);
        let cl = self.bilinear(&self.cl, alpha, beta);
        let cs = self.bilinear(&self.cs, alpha, beta);
        (cd, cl, cs)
    }
    /// Construct a trivial map with a single alpha/beta point (constant coefficients).
    pub fn from_constant(cd: f64, cl: f64, cs: f64) -> Self {
        Self {
            alphas: vec![0.0],
            betas: vec![0.0],
            cd: vec![cd],
            cl: vec![cl],
            cs: vec![cs],
        }
    }
    fn bilinear(&self, values: &[f64], alpha: f64, beta: f64) -> f64 {
        let na = self.alphas.len();
        let nb = self.betas.len();
        if na == 0 || nb == 0 {
            return 0.0;
        }
        let (ia, ta) = bracket_index(&self.alphas, alpha);
        let (ib, tb) = bracket_index(&self.betas, beta);
        let ia1 = (ia + 1).min(na - 1);
        let ib1 = (ib + 1).min(nb - 1);
        let v00 = values[ia * nb + ib];
        let v10 = values[ia1 * nb + ib];
        let v01 = values[ia * nb + ib1];
        let v11 = values[ia1 * nb + ib1];
        let va = v00 * (1.0 - ta) + v10 * ta;
        let vb = v01 * (1.0 - ta) + v11 * ta;
        va * (1.0 - tb) + vb * tb
    }
}
/// Advanced ground-effect model that accounts for the Venturi effect in the
/// underbody diffuser and produces separate downforce and drag corrections.
#[derive(Debug, Clone)]
pub struct GroundEffectAdvanced {
    /// Clearance from the ground surface in m.
    pub ride_height: f64,
    /// Reference ride height at which ground effect becomes significant (m).
    pub reference_height: f64,
    /// Maximum downforce amplification factor (at zero ride height limit).
    pub max_amplification: f64,
    /// Diffuser expansion angle in degrees.
    pub diffuser_angle: f64,
}
impl GroundEffectAdvanced {
    /// Create an advanced ground effect model for a racing car underbody.
    pub fn new_race_car() -> Self {
        Self {
            ride_height: 0.05,
            reference_height: 0.12,
            max_amplification: 4.0,
            diffuser_angle: 15.0,
        }
    }
    /// Venturi suction factor: `k = 1 + (h_ref / h)^2 * max_amp`.
    ///
    /// Clamped to `[1, max_amplification]`.
    pub fn venturi_suction_factor(&self) -> f64 {
        if self.ride_height <= 0.0 {
            return self.max_amplification;
        }
        let ratio = self.reference_height / self.ride_height;
        (1.0 + ratio * ratio).min(self.max_amplification)
    }
    /// Amplified downforce coefficient.
    pub fn amplified_cl(&self, base_cl: f64) -> f64 {
        base_cl * self.venturi_suction_factor()
    }
    /// Induced drag penalty from ground-effect downforce.
    ///
    /// Extra drag ≈ `0.1 * |ΔCl|` (empirical; diffuser drag).
    pub fn induced_drag_penalty(&self, base_cl: f64) -> f64 {
        let delta_cl = (self.amplified_cl(base_cl) - base_cl).abs();
        let diffuser_rad = self.diffuser_angle.to_radians();
        0.1 * delta_cl * diffuser_rad.tan()
    }
    /// Whether the vehicle is in the "ground effect stall" regime where the
    /// underbody flow separates (ride height too low).
    pub fn is_stalled(&self) -> bool {
        self.ride_height < 0.015
    }
}
/// Prandtl–Glauert compressibility correction for subsonic aerodynamics.
///
/// `Cl_c = Cl_0 / sqrt(1 - M²)`  where M is the Mach number.
///
/// Valid for M < 0.7.  Returns uncorrected coefficients beyond this.
#[derive(Debug, Clone)]
pub struct CompressibilityCorrection {
    /// Speed of sound in the current air conditions (m/s).
    pub speed_of_sound: f64,
}
impl CompressibilityCorrection {
    /// Speed of sound in dry air at sea level (340.3 m/s).
    pub fn sea_level() -> Self {
        Self {
            speed_of_sound: 340.3,
        }
    }
    /// Mach number for a given airspeed.
    pub fn mach_number(&self, speed: f64) -> f64 {
        speed / self.speed_of_sound
    }
    /// Prandtl–Glauert correction factor `1 / sqrt(1 - M²)`.
    ///
    /// Clipped to a maximum of 3.0 to avoid singularity near M=1.
    pub fn pg_factor(&self, speed: f64) -> f64 {
        let m = self.mach_number(speed);
        let m2 = (m * m).min(0.99);
        (1.0 / (1.0 - m2).sqrt()).min(3.0)
    }
    /// Corrected lift coefficient.
    pub fn corrected_cl(&self, cl0: f64, speed: f64) -> f64 {
        cl0 * self.pg_factor(speed)
    }
    /// Corrected drag coefficient (Karman–Tsien approximation).
    ///
    /// `Cd_c = Cd_0 / sqrt(1 - M²) + M² * Cl_0² / (2 * sqrt(1 - M²))`
    pub fn corrected_cd(&self, cd0: f64, cl0: f64, speed: f64) -> f64 {
        let m = self.mach_number(speed);
        let m2 = (m * m).min(0.99);
        let beta = (1.0 - m2).sqrt();
        cd0 / beta + m2 * cl0 * cl0 / (2.0 * beta)
    }
}
/// A finite lifting surface modeled by a simplified Vortex-Lattice Method.
///
/// The span is divided into `n_panels` equal strips.  Each strip is modelled
/// as a horseshoe vortex.  The circulation is solved via the Prandtl lifting-
/// line approximation (closed-form for an elliptical distribution).
#[derive(Debug, Clone)]
pub struct AeroSurface {
    /// Span of the surface (m).
    pub span: f64,
    /// Mean aerodynamic chord (m).
    pub chord: f64,
    /// Number of spanwise panels used in the VLM discretisation.
    pub n_panels: usize,
    /// Geometric angle of attack (degrees).
    pub alpha_deg: f64,
    /// Oswald span efficiency factor.
    pub oswald: f64,
}
impl AeroSurface {
    /// Create a new `AeroSurface`.
    pub fn new(span: f64, chord: f64, n_panels: usize, alpha_deg: f64, oswald: f64) -> Self {
        Self {
            span,
            chord,
            n_panels: n_panels.max(1),
            alpha_deg,
            oswald: oswald.clamp(0.1, 1.0),
        }
    }
    /// Aspect ratio AR = span / chord (for a constant-chord surface).
    pub fn aspect_ratio(&self) -> f64 {
        self.span / self.chord.max(1e-9)
    }
    /// Planform area S = span × chord.
    pub fn planform_area(&self) -> f64 {
        self.span * self.chord
    }
    /// Simplified Vortex-Lattice Method: compute the spanwise lift distribution.
    ///
    /// Uses the Prandtl lifting-line closed-form elliptic solution:
    /// `Gamma(y) = Gamma0 * sqrt(1 - (2y/b)^2)`
    ///
    /// where `Gamma0 = pi * V * cl_ref * chord / (pi + 2 * AR / e)` (approximate).
    ///
    /// Returns a vector of `(y_pos_m, local_cl)` pairs for each panel centre.
    ///
    /// # Arguments
    /// * `velocity` – freestream speed (m/s)
    pub fn vortex_lattice_method(&self, velocity: f64) -> Vec<(f64, f64)> {
        if velocity.abs() < 1e-9 {
            return vec![(0.0, 0.0); self.n_panels];
        }
        let ar = self.aspect_ratio();
        let alpha_rad = self.alpha_deg.to_radians();
        let cl_alpha_3d = 2.0 * std::f64::consts::PI / (1.0 + 2.0 / (ar * self.oswald));
        let cl_ref = cl_alpha_3d * alpha_rad;
        let n = self.n_panels;
        let b = self.span;
        let dy = b / n as f64;
        let mut distribution = Vec::with_capacity(n);
        for i in 0..n {
            let y = -b * 0.5 + (i as f64 + 0.5) * dy;
            let eta = 2.0 * y / b;
            let local_factor = (1.0 - eta * eta).max(0.0).sqrt();
            let local_cl = cl_ref * local_factor;
            distribution.push((y, local_cl));
        }
        distribution
    }
    /// Total lift coefficient from the VLM distribution.
    ///
    /// Integrates the spanwise `cl(y)` distribution over the span using the
    /// trapezoidal rule and divides by the planform area.
    pub fn total_lift_coefficient(&self, velocity: f64) -> f64 {
        let dist = self.vortex_lattice_method(velocity);
        if dist.is_empty() {
            return 0.0;
        }
        let n = dist.len() as f64;
        let dy = self.span / n;
        let sum: f64 = dist.iter().map(|(_, cl)| cl).sum();
        sum * dy / self.planform_area().max(1e-12)
    }
}
/// Compact aerodynamic body description using pre-multiplied area-coefficients.
///
/// `CdA = Cd × A`, `ClA = Cl × A`, `CsA = Cs × A` (all in m²).
/// This avoids needing to store Cd and A separately when only the product is
/// known from wind-tunnel or CFD data.
#[derive(Debug, Clone)]
pub struct AeroBody {
    /// Drag area CdA in m² (product of drag coefficient and reference area).
    pub cda: f64,
    /// Lift area ClA in m² (negative → downforce).
    pub cla: f64,
    /// Side-force area CsA in m².
    pub csa: f64,
    /// Position of the aerodynamic centre of pressure `[x, y, z]` in body
    /// frame (m).
    pub cop: [f64; 3],
}
impl AeroBody {
    /// Compute the 3-D force vector `[Fx, Fy, Fz]` for a velocity vector
    /// `velocity` (m/s in body frame, x = forward, y = left, z = up).
    ///
    /// - Drag acts opposite to the velocity direction.
    /// - Lift acts along +z (positive = upward; negative CsA → downforce).
    /// - Side-force acts along ±y.
    pub fn forces(&self, velocity: [f64; 3], air_density: f64) -> AeroForces {
        let drag_v = compute_drag(velocity, air_density, self.cda);
        let lift_v = compute_lift(velocity, air_density, self.cla, 0.0);
        let side_v = compute_side_force(velocity, air_density, self.csa);
        AeroForces {
            drag: drag_v,
            lift: lift_v,
            side: side_v,
        }
    }
}
impl AeroBody {
    /// Compute the drag polar: Cd as a function of angle-of-attack `alpha_deg`.
    ///
    /// Uses the classic parabolic drag polar:
    /// `Cd(alpha) = Cd0 + k * Cl(alpha)^2`
    ///
    /// where `k = 1 / (pi * AR * e)`, `AR` is the effective aspect ratio derived
    /// from the body geometry, and `e = 0.85` is the Oswald efficiency factor.
    ///
    /// Returns `(alpha_deg, Cd)`.
    pub fn compute_drag_polar(&self, alpha_deg: f64, aspect_ratio: f64, cd0: f64) -> (f64, f64) {
        let e = 0.85_f64;
        let cl = self.lift_coefficient_at_alpha(alpha_deg);
        let k = if aspect_ratio > 1e-6 {
            1.0 / (std::f64::consts::PI * aspect_ratio * e)
        } else {
            0.0
        };
        let cd = cd0 + k * cl * cl;
        (alpha_deg, cd)
    }
    /// Angle-of-attack lift coefficient using thin-airfoil theory with stall.
    ///
    /// `Cl = 2 * pi * alpha_rad` below stall (|alpha| <= 15 deg),
    /// cosine roll-off above.
    fn lift_coefficient_at_alpha(&self, alpha_deg: f64) -> f64 {
        let stall = 15.0_f64;
        let alpha_rad = alpha_deg.to_radians();
        if alpha_deg.abs() <= stall {
            2.0 * std::f64::consts::PI * alpha_rad
        } else {
            let cl_max = 2.0 * std::f64::consts::PI * stall.to_radians();
            let delta = (alpha_deg.abs() - stall).to_radians();
            alpha_deg.signum() * cl_max * delta.cos()
        }
    }
    /// Compute induced drag coefficient from the lift coefficient.
    ///
    /// Induced drag: `CDi = Cl^2 / (pi * AR * e)`
    ///
    /// # Arguments
    /// * `cl`           – lift coefficient at the operating condition
    /// * `aspect_ratio` – wing effective aspect ratio
    /// * `oswald`       – Oswald span efficiency factor (typical 0.7–0.95)
    ///
    /// Returns `CDi >= 0`.
    pub fn compute_induced_drag(cl: f64, aspect_ratio: f64, oswald: f64) -> f64 {
        if aspect_ratio < 1e-9 || oswald < 1e-9 {
            return 0.0;
        }
        let cdi = cl * cl / (std::f64::consts::PI * aspect_ratio * oswald);
        cdi.max(0.0)
    }
}
/// Ground-effect model based on ride height.
pub struct GroundEffect {
    /// Clearance from the ground surface in m.
    pub ride_height: f64,
}
impl GroundEffect {
    /// Enhancement factor due to ground proximity.
    ///
    /// `factor = 1 + 0.5 * min(0.3 / h, 3.0)`
    pub fn enhancement_factor(&self) -> f64 {
        1.0 + 0.5 * (0.3 / self.ride_height).min(3.0)
    }
    /// Effective lift coefficient after ground-effect amplification.
    ///
    /// For negative `base_cl` (downforce), ground effect makes it more negative.
    pub fn effective_lift_coefficient(&self, base_cl: f64) -> f64 {
        base_cl * self.enhancement_factor()
    }
}
