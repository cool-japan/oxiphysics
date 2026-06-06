// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Raycast vehicle simulation combining all vehicle subsystems.
//!
//! The [`RaycastVehicle`] simulates a vehicle by casting rays downward from
//! each wheel's suspension mount point. When a ray hits the ground, suspension
//! and tire forces are computed and applied to the chassis rigid body.

use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_rigid::RigidBody;

use crate::drivetrain::DrivetrainLegacy as Drivetrain;
use crate::steering::AckermannSteering;
use crate::suspension::SuspensionModel;
use crate::tire::TireModel;
use crate::wheel::{Wheel, WheelState};

/// A ground hit result from a suspension ray cast.
#[derive(Debug, Clone)]
pub struct GroundHit {
    /// World-space contact point.
    pub point: Vec3,
    /// World-space surface normal (pointing away from ground).
    pub normal: Vec3,
    /// Distance from ray origin to hit point.
    pub distance: Real,
    /// Surface friction coefficient.
    pub friction: Real,
}

/// Runtime state of the vehicle as a whole.
#[derive(Debug, Clone, Default)]
pub struct VehicleState {
    /// Forward speed in m/s (positive = forward).
    pub speed: Real,
    /// Lateral acceleration in g-units.
    pub lateral_g: Real,
    /// Yaw rate in rad/s.
    pub yaw_rate: Real,
    /// Current throttle input (0..1).
    pub throttle: Real,
    /// Current brake input (0..1).
    pub brake: Real,
    /// Current steering input (-1..1).
    pub steer_input: Real,
}

/// Aerodynamic drag parameters.
#[derive(Debug, Clone)]
pub struct AeroDrag {
    /// Drag coefficient * frontal area (Cd * A) in m^2.
    pub cd_area: Real,
    /// Air density in kg/m^3.
    pub air_density: Real,
}

impl Default for AeroDrag {
    fn default() -> Self {
        Self {
            cd_area: 0.8, // typical sedan
            air_density: 1.225,
        }
    }
}

impl AeroDrag {
    /// Compute aerodynamic drag force magnitude given speed.
    pub fn force(&self, speed: Real) -> Real {
        0.5 * self.air_density * self.cd_area * speed * speed
    }
}

/// Ground query callback type.
///
/// Given a ray origin and direction, returns an optional [`GroundHit`].
/// The direction vector is normalized and points downward along the suspension.
pub type GroundQuery = dyn Fn(&Vec3, &Vec3, Real) -> Option<GroundHit>;

/// Raycast vehicle simulation.
///
/// Manages wheels, drivetrain, steering, and applies computed forces
/// to a chassis rigid body each simulation step.
///
/// Wheel layout convention for a 4-wheel vehicle:
/// - Index 0: Front-Left
/// - Index 1: Front-Right
/// - Index 2: Rear-Left
/// - Index 3: Rear-Right
pub struct RaycastVehicle {
    /// Wheel configuration parameters.
    pub wheels: Vec<Wheel>,
    /// Per-wheel runtime state.
    pub wheel_states: Vec<WheelState>,
    /// Drivetrain (engine + gearbox + differentials).
    pub drivetrain: Drivetrain,
    /// Steering geometry.
    pub steering: AckermannSteering,
    /// Vehicle-level state (speed, lateral g, yaw rate).
    pub state: VehicleState,
    /// Aerodynamic drag parameters.
    pub aero: AeroDrag,
    /// Number of front-axle wheels (used for steering assignment).
    pub num_front_wheels: usize,
}

impl RaycastVehicle {
    /// Create a new raycast vehicle with the given wheels and subsystems.
    pub fn new(wheels: Vec<Wheel>, drivetrain: Drivetrain, steering: AckermannSteering) -> Self {
        let n = wheels.len();
        Self {
            wheels,
            wheel_states: vec![WheelState::default(); n],
            drivetrain,
            steering,
            state: VehicleState::default(),
            aero: AeroDrag::default(),
            num_front_wheels: 2,
        }
    }

    /// Create a standard 4-wheel vehicle with default parameters.
    pub fn new_default_4wheel() -> Self {
        let half_track = 0.75;
        let front_z = 1.25;
        let rear_z = -1.25;
        let mount_y = -0.1; // slightly below chassis center

        let make_wheel = |x: Real, z: Real| -> Wheel {
            Wheel {
                connection_point: Vec3::new(x, mount_y, z),
                ..Wheel::default()
            }
        };

        let wheels = vec![
            make_wheel(-half_track, front_z), // FL
            make_wheel(half_track, front_z),  // FR
            make_wheel(-half_track, rear_z),  // RL
            make_wheel(half_track, rear_z),   // RR
        ];

        Self::new(wheels, Drivetrain::default(), AckermannSteering::default())
    }

    /// Number of wheels.
    pub fn num_wheels(&self) -> usize {
        self.wheels.len()
    }

    /// Set the driver inputs.
    pub fn set_inputs(&mut self, throttle: Real, brake: Real, steer: Real) {
        self.state.throttle = throttle.clamp(0.0, 1.0);
        self.state.brake = brake.clamp(0.0, 1.0);
        self.state.steer_input = steer.clamp(-1.0, 1.0);
    }

    /// Update the vehicle for one simulation step.
    ///
    /// This is the main simulation entry point. It:
    /// 1. Applies steering angles to front wheels
    /// 2. Casts rays for each wheel to detect ground contact
    /// 3. Computes suspension forces
    /// 4. Computes tire slip and tire forces
    /// 5. Computes drivetrain torques
    /// 6. Applies all forces to the chassis body
    /// 7. Updates vehicle state (speed, yaw rate, etc.)
    ///
    /// # Arguments
    /// * `dt` - Time step in seconds
    /// * `chassis` - Mutable reference to the chassis rigid body
    /// * `gravity` - Gravity vector
    /// * `ground_query` - Function to query ground intersection
    /// * `suspension_model` - Suspension force model
    /// * `tire_model` - Tire force model
    pub fn update<S: SuspensionModel, T: TireModel>(
        &mut self,
        dt: Real,
        chassis: &mut RigidBody,
        gravity: &Vec3,
        ground_query: &GroundQuery,
        suspension_model: &S,
        tire_model: &T,
    ) {
        let num_wheels = self.wheels.len();
        if num_wheels == 0 {
            return;
        }

        // 1. Apply steering angles
        self.apply_steering();

        // 2-4. For each wheel: raycast, suspension, tire forces
        let mut total_force = Vec3::zeros();
        let mut total_torque = Vec3::zeros();
        let mut wheel_spin_speeds = [0.0_f64; 4];
        for (wss, ws) in wheel_spin_speeds
            .iter_mut()
            .zip(self.wheel_states.iter())
            .take(num_wheels.min(4))
        {
            *wss = ws.spin_velocity;
        }

        // 5. Compute drivetrain torques
        let drive_torques = self.drivetrain.compute_wheel_torques(
            self.state.throttle,
            self.state.brake,
            &wheel_spin_speeds,
        );

        for (i, (wheel, ws)) in self
            .wheels
            .iter()
            .zip(self.wheel_states.iter_mut())
            .enumerate()
            .take(num_wheels)
        {
            // Transform wheel mount point to world space
            let world_mount = chassis.transform.transform_point(&wheel.connection_point);
            let world_down = chassis
                .transform
                .transform_vector(&wheel.suspension_direction);
            let world_axle = chassis.transform.transform_vector(&wheel.axle_direction);

            // Apply steering rotation to axle and forward direction
            let steer_angle = ws.steering_angle;
            let (sin_s, cos_s) = steer_angle.sin_cos();

            // Forward direction is perpendicular to both axle and suspension dir
            let forward_base = world_down.cross(&world_axle);
            let forward_dir =
                forward_base * cos_s + world_axle.cross(&world_down).cross(&world_down) * sin_s;
            let forward_dir = if forward_dir.norm() > 1e-10 {
                forward_dir.normalize()
            } else {
                forward_base
            };

            // Steered axle direction
            let right_dir = world_down.cross(&forward_dir);
            let right_dir = if right_dir.norm() > 1e-10 {
                right_dir.normalize()
            } else {
                world_axle
            };

            ws.forward_direction = forward_dir;
            ws.right_direction = right_dir;

            // Raycast
            let max_ray_len =
                wheel.suspension_rest_length + wheel.max_suspension_travel + wheel.radius;
            let hit = ground_query(&world_mount, &world_down, max_ray_len);

            match hit {
                Some(ground_hit) => {
                    let was_in_contact = ws.is_in_contact;
                    ws.is_in_contact = true;
                    ws.contact_point = ground_hit.point;
                    ws.contact_normal = ground_hit.normal;

                    // Suspension length = hit distance - wheel radius
                    let susp_len = (ground_hit.distance - wheel.radius).max(0.0);
                    let susp_len =
                        susp_len.min(wheel.suspension_rest_length + wheel.max_suspension_travel);

                    // Suspension velocity (zero on first contact to avoid spike)
                    let susp_velocity = if was_in_contact {
                        let prev_len = ws.suspension_length;
                        (susp_len - prev_len) / dt.max(1e-10)
                    } else {
                        0.0
                    };
                    ws.suspension_length = susp_len;

                    // Choose damping based on compression vs extension
                    let damping = if susp_velocity < 0.0 {
                        wheel.damping_compression
                    } else {
                        wheel.damping_relaxation
                    };

                    // Compute suspension force
                    let susp_force = suspension_model.compute_force(
                        wheel.suspension_rest_length,
                        susp_len,
                        susp_velocity,
                        wheel.suspension_stiffness,
                        damping,
                    );
                    ws.suspension_force = susp_force;

                    // Apply suspension force along contact normal
                    let susp_force_vec = ground_hit.normal * susp_force;
                    let force_point = world_mount + world_down * susp_len;
                    ws.world_position = force_point;

                    // Compute wheel velocity at contact point
                    let r_contact = ground_hit.point - chassis.transform.position;
                    let contact_velocity =
                        chassis.velocity + chassis.angular_velocity.cross(&r_contact);

                    // Project velocity onto forward and lateral directions
                    let v_forward = contact_velocity.dot(&forward_dir);
                    let v_lateral = contact_velocity.dot(&right_dir);

                    // Compute slip angle: angle between wheel heading and velocity
                    ws.slip_angle = if v_forward.abs() > 0.5 {
                        (-v_lateral / v_forward.abs()).atan()
                    } else if (v_forward.abs() + v_lateral.abs()) > 0.1 {
                        (-v_lateral).atan2(v_forward.abs())
                    } else {
                        0.0
                    };

                    // Compute slip ratio
                    let wheel_speed_at_ground = ws.spin_velocity * wheel.radius;
                    let denominator = v_forward.abs().max(0.5);
                    ws.slip_ratio = (wheel_speed_at_ground - v_forward) / denominator;
                    ws.slip_ratio = ws.slip_ratio.clamp(-1.0, 1.0);

                    // Compute tire forces
                    let (fy, fx) = tire_model.compute_forces(
                        ws.slip_angle,
                        ws.slip_ratio,
                        susp_force,
                        ground_hit.friction * wheel.friction_slip,
                    );
                    ws.lateral_force = fy;
                    ws.longitudinal_force = fx;

                    // Apply forces
                    let tire_force = forward_dir * fx + right_dir * fy;
                    let contact_force = susp_force_vec + tire_force;

                    total_force += contact_force;
                    let r_force = force_point - chassis.transform.position;
                    total_torque += r_force.cross(&contact_force);

                    // Update wheel spin from drive torque
                    let wheel_inertia = 0.5 * wheel.mass * wheel.radius * wheel.radius;
                    if wheel_inertia > 1e-10 && i < 4 {
                        let drive_torque = drive_torques[i];
                        let tire_resistance = -fx * wheel.radius;
                        let net_torque = drive_torque + tire_resistance;
                        ws.spin_velocity += (net_torque / wheel_inertia) * dt;
                        ws.spin_velocity = ws.spin_velocity.max(0.0); // prevent reverse spin from braking past zero
                    }
                }
                None => {
                    ws.is_in_contact = false;
                    ws.suspension_length =
                        wheel.suspension_rest_length + wheel.max_suspension_travel;
                    ws.suspension_force = 0.0;
                    ws.lateral_force = 0.0;
                    ws.longitudinal_force = 0.0;
                    ws.world_position = world_mount
                        + world_down * (wheel.suspension_rest_length + wheel.max_suspension_travel);

                    // Wheel spin decays when not in contact
                    ws.spin_velocity *= 0.99;
                }
            }

            // Update visual rotation angle
            ws.rotation_angle += ws.spin_velocity * dt;
        }

        // 6. Apply aerodynamic drag
        let forward = chassis
            .transform
            .transform_vector(&Vec3::new(0.0, 0.0, 1.0));
        let speed = chassis.velocity.dot(&forward);
        let drag_force = -forward * self.aero.force(speed.abs()) * speed.signum();
        total_force += drag_force;

        // Apply all computed forces to chassis
        chassis.apply_force(total_force);
        chassis.apply_torque(total_torque);

        // 7. Update vehicle state
        self.state.speed = speed;
        self.state.yaw_rate = chassis.angular_velocity.dot(&Vec3::new(0.0, 1.0, 0.0));

        // Lateral g from centripetal acceleration
        let right = chassis
            .transform
            .transform_vector(&Vec3::new(1.0, 0.0, 0.0));
        let lateral_accel = chassis.velocity.dot(&right); // approximate
        let g_mag = gravity.norm();
        self.state.lateral_g = if g_mag > 1e-10 {
            (lateral_accel * self.state.yaw_rate) / g_mag
        } else {
            0.0
        };
    }

    /// Apply steering angles to the front wheels.
    fn apply_steering(&mut self) {
        if self.wheels.len() < 2 {
            return;
        }

        let (left, right) = self.steering.compute_wheel_angles(self.state.steer_input);

        // Apply to front wheels
        for i in 0..self.num_front_wheels.min(self.wheel_states.len()) {
            if i % 2 == 0 {
                self.wheel_states[i].steering_angle = left;
            } else {
                self.wheel_states[i].steering_angle = right;
            }
        }
    }

    /// Get the current speed in km/h.
    pub fn speed_kmh(&self) -> Real {
        self.state.speed * 3.6
    }

    /// Get the current speed in mph.
    pub fn speed_mph(&self) -> Real {
        self.state.speed * 2.23694
    }

    /// Get the engine RPM.
    pub fn engine_rpm(&self) -> Real {
        self.drivetrain.engine_rpm as Real
    }

    /// Get the current gear number.
    pub fn current_gear(&self) -> i32 {
        self.drivetrain.gearbox.current_gear
    }

    /// Shift up one gear.
    pub fn shift_up(&mut self) {
        self.drivetrain.gearbox.shift_up();
    }

    /// Shift down one gear.
    pub fn shift_down(&mut self) {
        self.drivetrain.gearbox.shift_down();
    }

    /// Check if any wheel is in contact with the ground.
    pub fn any_wheel_in_contact(&self) -> bool {
        self.wheel_states.iter().any(|ws| ws.is_in_contact)
    }

    /// Check if all wheels are in contact with the ground.
    pub fn all_wheels_in_contact(&self) -> bool {
        self.wheel_states.iter().all(|ws| ws.is_in_contact)
    }
}

/// A simple flat ground query for testing.
///
/// Returns a hit on a flat plane at the given Y height with friction 1.0.
pub fn flat_ground_query(ground_y: Real) -> Box<GroundQuery> {
    Box::new(
        move |origin: &Vec3, direction: &Vec3, max_dist: Real| -> Option<GroundHit> {
            // Only consider downward rays
            if direction.y >= 0.0 {
                return None;
            }
            let t = (ground_y - origin.y) / direction.y;
            if t < 0.0 || t > max_dist {
                return None;
            }
            Some(GroundHit {
                point: *origin + *direction * t,
                normal: Vec3::new(0.0, 1.0, 0.0),
                distance: t,
                friction: 1.0,
            })
        },
    )
}

// ---------------------------------------------------------------------------
// Vehicle parameter validation
// ---------------------------------------------------------------------------

/// Validation errors for vehicle parameters.
#[derive(Debug, Clone)]
pub enum VehicleValidationError {
    /// Mass must be positive.
    NonPositiveMass(Real),
    /// Wheel radius must be positive.
    NonPositiveWheelRadius {
        /// Index of the wheel with invalid radius.
        wheel_index: usize,
        /// The invalid radius value.
        radius: Real,
    },
    /// Suspension stiffness must be positive.
    NonPositiveSuspensionStiffness {
        /// Index of the wheel with invalid stiffness.
        wheel_index: usize,
        /// The invalid stiffness value.
        stiffness: Real,
    },
    /// CdA must be non-negative.
    NegativeCdArea(Real),
    /// Air density must be positive.
    NonPositiveAirDensity(Real),
    /// No wheels defined.
    NoWheels,
}

impl std::fmt::Display for VehicleValidationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            VehicleValidationError::NonPositiveMass(m) => {
                write!(f, "Mass must be positive, got {m}")
            }
            VehicleValidationError::NonPositiveWheelRadius {
                wheel_index,
                radius,
            } => {
                write!(
                    f,
                    "Wheel {wheel_index} radius must be positive, got {radius}"
                )
            }
            VehicleValidationError::NonPositiveSuspensionStiffness {
                wheel_index,
                stiffness,
            } => {
                write!(
                    f,
                    "Wheel {wheel_index} suspension stiffness must be positive, got {stiffness}"
                )
            }
            VehicleValidationError::NegativeCdArea(cd) => {
                write!(f, "CdA must be non-negative, got {cd}")
            }
            VehicleValidationError::NonPositiveAirDensity(rho) => {
                write!(f, "Air density must be positive, got {rho}")
            }
            VehicleValidationError::NoWheels => {
                write!(f, "Vehicle must have at least one wheel")
            }
        }
    }
}

/// Validate vehicle parameters and return a list of errors (empty if valid).
pub fn validate_vehicle(
    vehicle: &RaycastVehicle,
    chassis_mass: Real,
) -> Vec<VehicleValidationError> {
    let mut errors = Vec::new();

    if chassis_mass <= 0.0 {
        errors.push(VehicleValidationError::NonPositiveMass(chassis_mass));
    }

    if vehicle.wheels.is_empty() {
        errors.push(VehicleValidationError::NoWheels);
    }

    for (i, w) in vehicle.wheels.iter().enumerate() {
        if w.radius <= 0.0 {
            errors.push(VehicleValidationError::NonPositiveWheelRadius {
                wheel_index: i,
                radius: w.radius,
            });
        }
        if w.suspension_stiffness <= 0.0 {
            errors.push(VehicleValidationError::NonPositiveSuspensionStiffness {
                wheel_index: i,
                stiffness: w.suspension_stiffness,
            });
        }
    }

    if vehicle.aero.cd_area < 0.0 {
        errors.push(VehicleValidationError::NegativeCdArea(vehicle.aero.cd_area));
    }

    if vehicle.aero.air_density <= 0.0 {
        errors.push(VehicleValidationError::NonPositiveAirDensity(
            vehicle.aero.air_density,
        ));
    }

    errors
}

// ---------------------------------------------------------------------------
// Aerodynamic downforce
// ---------------------------------------------------------------------------

/// Aerodynamic downforce parameters.
#[derive(Debug, Clone)]
pub struct AeroDownforce {
    /// Lift coefficient * frontal area (Cl * A) in m^2 (negative = downforce).
    pub cl_area: Real,
    /// Air density in kg/m^3.
    pub air_density: Real,
}

impl Default for AeroDownforce {
    fn default() -> Self {
        Self {
            cl_area: -1.5, // typical formula car downforce
            air_density: 1.225,
        }
    }
}

impl AeroDownforce {
    /// Compute aerodynamic downforce magnitude given speed.
    /// Returns a positive force value when cl_area is negative (downforce).
    pub fn force(&self, speed: Real) -> Real {
        -0.5 * self.air_density * self.cl_area * speed * speed
    }
}

// ---------------------------------------------------------------------------
// Weight transfer computation
// ---------------------------------------------------------------------------

/// Weight transfer result for a vehicle under longitudinal and lateral acceleration.
#[derive(Debug, Clone, Default)]
pub struct WeightTransfer {
    /// Normal load on front-left wheel (N).
    pub front_left: Real,
    /// Normal load on front-right wheel (N).
    pub front_right: Real,
    /// Normal load on rear-left wheel (N).
    pub rear_left: Real,
    /// Normal load on rear-right wheel (N).
    pub rear_right: Real,
}

impl WeightTransfer {
    /// Total load across all four wheels.
    pub fn total_load(&self) -> Real {
        self.front_left + self.front_right + self.rear_left + self.rear_right
    }

    /// Front axle load.
    pub fn front_axle_load(&self) -> Real {
        self.front_left + self.front_right
    }

    /// Rear axle load.
    pub fn rear_axle_load(&self) -> Real {
        self.rear_left + self.rear_right
    }

    /// Left-side load.
    pub fn left_load(&self) -> Real {
        self.front_left + self.rear_left
    }

    /// Right-side load.
    pub fn right_load(&self) -> Real {
        self.front_right + self.rear_right
    }
}

/// Compute static and dynamic weight transfer for a 4-wheel vehicle.
///
/// # Arguments
/// * `mass` - Vehicle mass in kg
/// * `wheelbase` - Distance between front and rear axles (m)
/// * `track_width` - Distance between left and right wheels (m)
/// * `cg_height` - Center of gravity height (m)
/// * `front_weight_ratio` - Static front weight distribution (0..1)
/// * `long_accel` - Longitudinal acceleration (m/s^2, positive = accelerating)
/// * `lat_accel` - Lateral acceleration (m/s^2, positive = turning right)
pub fn compute_weight_transfer(
    mass: Real,
    wheelbase: Real,
    track_width: Real,
    cg_height: Real,
    front_weight_ratio: Real,
    long_accel: Real,
    lat_accel: Real,
) -> WeightTransfer {
    let g = 9.81;
    let total_weight = mass * g;

    // Static loads
    let front_static = total_weight * front_weight_ratio;
    let rear_static = total_weight * (1.0 - front_weight_ratio);

    // Longitudinal weight transfer: delta_F = m * a_x * h / L
    let long_transfer = if wheelbase.abs() > 1e-10 {
        mass * long_accel * cg_height / wheelbase
    } else {
        0.0
    };

    // Lateral weight transfer: delta_F = m * a_y * h / T
    let lat_transfer = if track_width.abs() > 1e-10 {
        mass * lat_accel * cg_height / track_width
    } else {
        0.0
    };

    // Front axle load after longitudinal transfer (braking adds to front)
    let front_total = front_static - long_transfer;
    let rear_total = rear_static + long_transfer;

    // Distribute lateral transfer proportionally to axle load
    // Positive lat_accel (right turn) -> more load on left
    let lat_front = lat_transfer * front_weight_ratio;
    let lat_rear = lat_transfer * (1.0 - front_weight_ratio);

    WeightTransfer {
        front_left: front_total * 0.5 + lat_front * 0.5,
        front_right: front_total * 0.5 - lat_front * 0.5,
        rear_left: rear_total * 0.5 + lat_rear * 0.5,
        rear_right: rear_total * 0.5 - lat_rear * 0.5,
    }
}

/// Simplified weight transfer: returns \[front, rear\] axle loads under longitudinal
/// acceleration only.
pub fn longitudinal_weight_transfer(
    mass: Real,
    wheelbase: Real,
    cg_height: Real,
    front_weight_ratio: Real,
    long_accel: Real,
) -> [Real; 2] {
    let g = 9.81;
    let total_weight = mass * g;
    let transfer = if wheelbase.abs() > 1e-10 {
        mass * long_accel * cg_height / wheelbase
    } else {
        0.0
    };
    let front = total_weight * front_weight_ratio - transfer;
    let rear = total_weight * (1.0 - front_weight_ratio) + transfer;
    [front, rear]
}

// ---------------------------------------------------------------------------
// Force balance utilities
// ---------------------------------------------------------------------------

/// Compute the net longitudinal force on the vehicle.
///
/// `engine_force` - propulsive force (N)
/// `brake_force` - braking force (N, positive = retarding)
/// `drag_force` - aerodynamic drag (N, always positive)
/// `rolling_resistance` - rolling resistance force (N, always positive)
pub fn net_longitudinal_force(
    engine_force: Real,
    brake_force: Real,
    drag_force: Real,
    rolling_resistance: Real,
) -> Real {
    engine_force - brake_force - drag_force - rolling_resistance
}

/// Compute rolling resistance force.
///
/// `normal_load` - total normal load on driven wheels (N)
/// `crr` - coefficient of rolling resistance (dimensionless, typically 0.01-0.02)
pub fn rolling_resistance_force(normal_load: Real, crr: Real) -> Real {
    normal_load * crr
}

/// Compute the maximum cornering speed given available lateral grip.
///
/// `radius` - corner radius (m)
/// `mu` - lateral friction coefficient
/// `downforce_coeff` - additional downforce as a fraction of weight (e.g., 0.5 = 50% more)
pub fn max_cornering_speed(radius: Real, mu: Real, downforce_coeff: Real) -> Real {
    let g = 9.81;
    // v^2 = mu * g * R * (1 + downforce_coeff)
    (mu * g * radius * (1.0 + downforce_coeff)).sqrt()
}

// ---------------------------------------------------------------------------
// Vehicle state integration (simple Euler)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Suspension force coupling (anti-roll bar / bump-steer)
// ---------------------------------------------------------------------------

/// Anti-roll bar parameters for one axle.
#[derive(Debug, Clone)]
pub struct AntiRollBar {
    /// Anti-roll bar stiffness (N/m).
    pub stiffness: Real,
}

impl Default for AntiRollBar {
    fn default() -> Self {
        Self { stiffness: 10000.0 }
    }
}

impl AntiRollBar {
    /// Create with a given stiffness.
    pub fn new(stiffness: Real) -> Self {
        Self { stiffness }
    }

    /// Compute the roll-couple force for the left and right wheels.
    ///
    /// `left_compression` and `right_compression` are the suspension
    /// compression displacements (rest_length - current_length).
    ///
    /// Returns `(force_left, force_right)` — positive = push wheel away from body.
    pub fn compute_forces(&self, left_compression: Real, right_compression: Real) -> (Real, Real) {
        let delta = left_compression - right_compression;
        let f = self.stiffness * delta * 0.5;
        (-f, f)
    }
}

/// Apply anti-roll bar forces to a set of wheel suspension forces for one axle.
///
/// `left_idx` and `right_idx` index into `suspension_forces`.
pub fn apply_anti_roll_bar(
    bar: &AntiRollBar,
    left_compression: Real,
    right_compression: Real,
    suspension_forces: &mut [Real],
    left_idx: usize,
    right_idx: usize,
) {
    let (fl, fr) = bar.compute_forces(left_compression, right_compression);
    if left_idx < suspension_forces.len() {
        suspension_forces[left_idx] += fl;
    }
    if right_idx < suspension_forces.len() {
        suspension_forces[right_idx] += fr;
    }
}

// ---------------------------------------------------------------------------
// Aerodynamic downforce with CoP offset
// ---------------------------------------------------------------------------

/// Full aerodynamic model including lift/drag and centre-of-pressure offset.
#[derive(Debug, Clone)]
pub struct AeroModel {
    /// Drag coefficient × frontal area (m²).
    pub cd_area: Real,
    /// Downforce coefficient × reference area (m²).  Positive = downforce.
    pub cl_area: Real,
    /// Air density (kg/m³).
    pub air_density: Real,
    /// Centre-of-pressure position along the vehicle's longitudinal axis (m,
    /// positive = forward from body origin).
    pub cop_offset: Real,
}

impl Default for AeroModel {
    fn default() -> Self {
        Self {
            cd_area: 0.8,
            cl_area: 1.2,
            air_density: 1.225,
            cop_offset: -0.2,
        }
    }
}

impl AeroModel {
    /// Aerodynamic drag force (N) for a given airspeed.
    pub fn drag_force(&self, airspeed: Real) -> Real {
        0.5 * self.air_density * self.cd_area * airspeed * airspeed
    }

    /// Aerodynamic downforce (N) for a given airspeed.
    pub fn downforce(&self, airspeed: Real) -> Real {
        0.5 * self.air_density * self.cl_area * airspeed * airspeed
    }

    /// Pitching moment (N·m) induced by downforce acting at the CoP offset.
    pub fn pitch_moment(&self, airspeed: Real) -> Real {
        self.downforce(airspeed) * self.cop_offset
    }

    /// Longitudinal load transfer induced by aerodynamic pitch moment.
    ///
    /// `wheelbase` in metres; returns `[front_delta, rear_delta]` in N.
    pub fn aero_load_transfer(&self, airspeed: Real, wheelbase: Real) -> [Real; 2] {
        if wheelbase.abs() < 1e-12 {
            return [0.0, 0.0];
        }
        let m = self.pitch_moment(airspeed);
        let delta = m / wheelbase;
        [delta, -delta]
    }
}

// ---------------------------------------------------------------------------
// Chassis flex simulation (simple torsional beam)
// ---------------------------------------------------------------------------

/// Models chassis torsional flexibility as a spring-damper between front and rear sub-frames.
#[derive(Debug, Clone)]
pub struct ChassisFlex {
    /// Torsional stiffness (N·m/rad).
    pub torsional_stiffness: Real,
    /// Torsional damping (N·m·s/rad).
    pub torsional_damping: Real,
    /// Current torsional deflection angle (rad).
    pub deflection: Real,
    /// Current torsional rate (rad/s).
    pub deflection_rate: Real,
}

impl Default for ChassisFlex {
    fn default() -> Self {
        Self {
            torsional_stiffness: 50000.0,
            torsional_damping: 2000.0,
            deflection: 0.0,
            deflection_rate: 0.0,
        }
    }
}

impl ChassisFlex {
    /// Create a new chassis flex model.
    pub fn new(stiffness: Real, damping: Real) -> Self {
        Self {
            torsional_stiffness: stiffness,
            torsional_damping: damping,
            deflection: 0.0,
            deflection_rate: 0.0,
        }
    }

    /// Apply an external torque `torque` (N·m) and integrate for `dt` seconds.
    ///
    /// Returns the updated deflection angle (rad).
    pub fn step(&mut self, external_torque: Real, moment_of_inertia: Real, dt: Real) -> Real {
        if moment_of_inertia.abs() < 1e-12 || dt < 1e-12 {
            return self.deflection;
        }
        let spring_torque = -self.torsional_stiffness * self.deflection;
        let damping_torque = -self.torsional_damping * self.deflection_rate;
        let net_torque = external_torque + spring_torque + damping_torque;
        let angular_accel = net_torque / moment_of_inertia;
        self.deflection_rate += angular_accel * dt;
        self.deflection += self.deflection_rate * dt;
        self.deflection
    }

    /// Restoring torque opposing the deflection (N·m).
    pub fn restoring_torque(&self) -> Real {
        -self.torsional_stiffness * self.deflection
    }

    /// Kinetic energy stored in the torsional mode (J).
    pub fn kinetic_energy(&self, moment_of_inertia: Real) -> Real {
        0.5 * moment_of_inertia * self.deflection_rate * self.deflection_rate
    }

    /// Potential energy stored in the torsional spring (J).
    pub fn potential_energy(&self) -> Real {
        0.5 * self.torsional_stiffness * self.deflection * self.deflection
    }
}

// ---------------------------------------------------------------------------
// Wheel load transfer with coupled aero + mechanical
// ---------------------------------------------------------------------------

/// Combined wheel load (N) for a 4-wheel vehicle including static weight,
/// mechanical weight transfer, and aerodynamic downforce distribution.
pub fn full_wheel_loads(
    mass: Real,
    wheelbase: Real,
    track_width: Real,
    cg_height: Real,
    front_weight_ratio: Real,
    long_accel: Real,
    lat_accel: Real,
    aero_model: &AeroModel,
    airspeed: Real,
) -> [Real; 4] {
    // Mechanical weight transfer
    let mech = compute_weight_transfer(
        mass,
        wheelbase,
        track_width,
        cg_height,
        front_weight_ratio,
        long_accel,
        lat_accel,
    );
    // Aerodynamic downforce distributed evenly
    let df_total = aero_model.downforce(airspeed);
    let df_per_wheel = df_total * 0.25;
    // Aero longitudinal load transfer
    let [aero_front_delta, aero_rear_delta] = aero_model.aero_load_transfer(airspeed, wheelbase);
    [
        (mech.front_left + df_per_wheel + aero_front_delta * 0.5).max(0.0),
        (mech.front_right + df_per_wheel + aero_front_delta * 0.5).max(0.0),
        (mech.rear_left + df_per_wheel + aero_rear_delta * 0.5).max(0.0),
        (mech.rear_right + df_per_wheel + aero_rear_delta * 0.5).max(0.0),
    ]
}

// ---------------------------------------------------------------------------
// Full vehicle state integration (RK4 1-D)
// ---------------------------------------------------------------------------

/// Integrate vehicle speed using a simple RK4 step given a net force function.
///
/// `f_net(speed)` returns the net longitudinal force (N) at a given speed.
/// Returns the new speed after `dt` seconds.
pub fn integrate_speed_rk4(
    speed: Real,
    mass: Real,
    dt: Real,
    f_net: impl Fn(Real) -> Real,
) -> Real {
    if mass < 1e-12 {
        return speed;
    }
    let a = |v: Real| f_net(v) / mass;
    let k1 = a(speed);
    let k2 = a(speed + 0.5 * dt * k1);
    let k3 = a(speed + 0.5 * dt * k2);
    let k4 = a(speed + dt * k3);
    let new_speed = speed + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
    new_speed.max(0.0)
}

/// Stopping distance under constant braking deceleration (m).
///
/// `mu_brake` × g is the deceleration; returns metres to stop from `initial_speed` (m/s).
pub fn braking_distance(initial_speed: Real, mu_brake: Real) -> Real {
    let g = 9.81;
    let decel = mu_brake * g;
    if decel < 1e-12 {
        return f64::INFINITY;
    }
    (initial_speed * initial_speed) / (2.0 * decel)
}

/// Compute the traction-limited acceleration for a given wheel load and friction coefficient.
///
/// Returns the maximum achievable longitudinal acceleration (m/s²).
pub fn traction_limited_acceleration(wheel_loads: &[Real], friction: Real, mass: Real) -> Real {
    if mass < 1e-12 {
        return 0.0;
    }
    let total_load: Real = wheel_loads.iter().sum();
    total_load * friction / mass
}

/// Simple vehicle state for 1D longitudinal simulation.
#[derive(Debug, Clone)]
pub struct SimpleVehicleState {
    /// Position along track (m).
    pub position: Real,
    /// Speed (m/s).
    pub speed: Real,
    /// Acceleration (m/s^2).
    pub acceleration: Real,
}

impl Default for SimpleVehicleState {
    fn default() -> Self {
        Self {
            position: 0.0,
            speed: 0.0,
            acceleration: 0.0,
        }
    }
}

impl SimpleVehicleState {
    /// Integrate one time step using Euler method.
    pub fn integrate(&mut self, dt: Real) {
        self.speed += self.acceleration * dt;
        if self.speed < 0.0 {
            self.speed = 0.0;
        }
        self.position += self.speed * dt;
    }

    /// Apply a net force to compute acceleration (F = ma).
    pub fn apply_force(&mut self, force: Real, mass: Real) {
        if mass > 1e-10 {
            self.acceleration = force / mass;
        }
    }

    /// Kinetic energy (J).
    pub fn kinetic_energy(&self, mass: Real) -> Real {
        0.5 * mass * self.speed * self.speed
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::drivetrain::DriveLayout;
    use crate::suspension::LinearSuspension;
    use crate::tire::LinearTire;

    #[test]
    fn test_vehicle_creation() {
        let v = RaycastVehicle::new_default_4wheel();
        assert_eq!(v.num_wheels(), 4);
        assert_eq!(v.wheel_states.len(), 4);
    }

    #[test]
    fn test_vehicle_inputs() {
        let mut v = RaycastVehicle::new_default_4wheel();
        v.set_inputs(0.8, 0.2, -0.5);
        assert!((v.state.throttle - 0.8).abs() < 1e-10);
        assert!((v.state.brake - 0.2).abs() < 1e-10);
        assert!((v.state.steer_input - (-0.5)).abs() < 1e-10);
    }

    #[test]
    fn test_vehicle_input_clamping() {
        let mut v = RaycastVehicle::new_default_4wheel();
        v.set_inputs(2.0, -0.5, 3.0);
        assert!((v.state.throttle - 1.0).abs() < 1e-10);
        assert!((v.state.brake - 0.0).abs() < 1e-10);
        assert!((v.state.steer_input - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_wheel_contact_detection() {
        let mut v = RaycastVehicle::new_default_4wheel();
        let mut chassis = RigidBody::new(1500.0);
        chassis.transform.position = Vec3::new(0.0, 0.6, 0.0);
        chassis.gravity_scale = 0.0;
        chassis.linear_damping = 0.0;

        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let ground = flat_ground_query(0.0);
        let suspension = LinearSuspension;
        let tire = LinearTire::default();

        v.set_inputs(0.0, 0.0, 0.0);
        v.update(0.01, &mut chassis, &gravity, &ground, &suspension, &tire);

        // All wheels should be in contact when chassis is at reasonable height
        assert!(
            v.all_wheels_in_contact(),
            "Not all wheels in contact: {:?}",
            v.wheel_states
                .iter()
                .map(|w| w.is_in_contact)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_vehicle_no_contact_high_up() {
        let mut v = RaycastVehicle::new_default_4wheel();
        let mut chassis = RigidBody::new(1500.0);
        chassis.transform.position = Vec3::new(0.0, 100.0, 0.0);
        chassis.gravity_scale = 0.0;

        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let ground = flat_ground_query(0.0);
        let suspension = LinearSuspension;
        let tire = LinearTire::default();

        v.set_inputs(0.0, 0.0, 0.0);
        v.update(0.01, &mut chassis, &gravity, &ground, &suspension, &tire);

        assert!(!v.any_wheel_in_contact());
    }

    #[test]
    fn test_vehicle_straight_line_acceleration() {
        let mut v = RaycastVehicle::new_default_4wheel();
        v.drivetrain.layout = DriveLayout::RearWheelDrive;
        v.drivetrain.gearbox.shift(1);

        let mut chassis = RigidBody::new(1500.0);
        chassis.transform.position = Vec3::new(0.0, 0.5, 0.0);
        chassis.gravity_scale = 0.0;
        chassis.linear_damping = 0.0;
        chassis.angular_damping = 0.0;

        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let ground = flat_ground_query(0.0);
        let suspension = LinearSuspension;
        let tire = LinearTire::default();

        v.set_inputs(1.0, 0.0, 0.0);

        // Run several steps
        for _ in 0..100 {
            chassis.force_accumulator = Vec3::zeros();
            chassis.torque_accumulator = Vec3::zeros();
            v.update(0.01, &mut chassis, &gravity, &ground, &suspension, &tire);
            chassis.integrate_forces(0.01, &Vec3::zeros());
            chassis.integrate_velocity(0.01);
        }

        // Vehicle should have accelerated forward
        let forward_speed = chassis.velocity.z;
        assert!(
            forward_speed > 0.1,
            "Expected forward acceleration, got speed={forward_speed}"
        );
    }

    #[test]
    fn test_vehicle_terminal_velocity_with_drag() {
        let mut v = RaycastVehicle::new_default_4wheel();
        v.drivetrain.layout = DriveLayout::RearWheelDrive;
        v.drivetrain.gearbox.shift(1);
        v.aero.cd_area = 2.0; // high drag for faster convergence

        let mut chassis = RigidBody::new(1500.0);
        chassis.transform.position = Vec3::new(0.0, 0.5, 0.0);
        chassis.gravity_scale = 0.0;
        chassis.linear_damping = 0.0;
        chassis.angular_damping = 0.0;

        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let ground = flat_ground_query(0.0);
        let suspension = LinearSuspension;
        let tire = LinearTire::default();

        v.set_inputs(1.0, 0.0, 0.0);

        let mut prev_speed = 0.0;
        let mut speed_increasing_count = 0;
        let mut speed_stable_count = 0;

        for step in 0..2000 {
            chassis.force_accumulator = Vec3::zeros();
            chassis.torque_accumulator = Vec3::zeros();
            v.update(0.01, &mut chassis, &gravity, &ground, &suspension, &tire);
            chassis.integrate_forces(0.01, &Vec3::zeros());
            chassis.integrate_velocity(0.01);

            let speed = chassis.velocity.z;
            if step > 10 {
                let delta = speed - prev_speed;
                if delta > 0.001 {
                    speed_increasing_count += 1;
                } else if delta.abs() < 0.01 {
                    speed_stable_count += 1;
                }
            }
            prev_speed = speed;
        }

        // Should have had an initial acceleration phase
        assert!(
            speed_increasing_count > 5,
            "Not enough acceleration: {speed_increasing_count}"
        );
        // And a stable/terminal velocity phase
        assert!(
            speed_stable_count > 10,
            "Not enough stable speed steps: {speed_stable_count}"
        );
    }

    #[test]
    fn test_flat_ground_query_miss() {
        let query = flat_ground_query(0.0);
        // Ray pointing up should miss
        let result = query(&Vec3::new(0.0, 1.0, 0.0), &Vec3::new(0.0, 1.0, 0.0), 10.0);
        assert!(result.is_none());
    }

    #[test]
    fn test_flat_ground_query_hit() {
        let query = flat_ground_query(0.0);
        let result = query(&Vec3::new(0.0, 1.0, 0.0), &Vec3::new(0.0, -1.0, 0.0), 10.0);
        assert!(result.is_some());
        let hit = result.unwrap();
        assert!((hit.distance - 1.0).abs() < 1e-10);
        assert!((hit.point.y - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_aero_drag_force() {
        let aero = AeroDrag::default();
        let f0 = aero.force(0.0);
        assert!((f0).abs() < 1e-10);

        let f30 = aero.force(30.0);
        // F = 0.5 * 1.225 * 0.8 * 900 = 441 N
        assert!((f30 - 441.0).abs() < 1.0, "f30={f30}");
    }

    #[test]
    fn test_speed_conversions() {
        let mut v = RaycastVehicle::new_default_4wheel();
        v.state.speed = 10.0; // 10 m/s
        assert!((v.speed_kmh() - 36.0).abs() < 0.1);
        assert!((v.speed_mph() - 22.37).abs() < 0.1);
    }

    // --- Vehicle validation tests ---

    #[test]
    fn test_validate_vehicle_valid() {
        let v = RaycastVehicle::new_default_4wheel();
        let errors = validate_vehicle(&v, 1500.0);
        assert!(
            errors.is_empty(),
            "valid vehicle should have no errors: {:?}",
            errors
        );
    }

    #[test]
    fn test_validate_vehicle_negative_mass() {
        let v = RaycastVehicle::new_default_4wheel();
        let errors = validate_vehicle(&v, -100.0);
        assert!(!errors.is_empty());
        assert!(matches!(
            errors[0],
            VehicleValidationError::NonPositiveMass(_)
        ));
    }

    #[test]
    fn test_validate_vehicle_negative_cd_area() {
        let mut v = RaycastVehicle::new_default_4wheel();
        v.aero.cd_area = -1.0;
        let errors = validate_vehicle(&v, 1500.0);
        assert!(!errors.is_empty());
    }

    #[test]
    fn test_validate_vehicle_zero_air_density() {
        let mut v = RaycastVehicle::new_default_4wheel();
        v.aero.air_density = 0.0;
        let errors = validate_vehicle(&v, 1500.0);
        assert!(!errors.is_empty());
    }

    #[test]
    fn test_validation_error_display() {
        let err = VehicleValidationError::NonPositiveMass(-5.0);
        let msg = format!("{err}");
        assert!(msg.contains("Mass"));
    }

    // --- AeroDownforce tests ---

    #[test]
    fn test_aero_downforce_zero_speed() {
        let df = AeroDownforce::default();
        assert!((df.force(0.0)).abs() < 1e-10);
    }

    #[test]
    fn test_aero_downforce_positive_at_speed() {
        let df = AeroDownforce::default();
        let f = df.force(50.0);
        assert!(f > 0.0, "downforce should be positive, got {f}");
    }

    #[test]
    fn test_aero_downforce_scales_with_speed_squared() {
        let df = AeroDownforce::default();
        let f1 = df.force(10.0);
        let f2 = df.force(20.0);
        assert!(
            (f2 / f1 - 4.0).abs() < 0.01,
            "downforce should scale with v^2"
        );
    }

    // --- Weight transfer tests ---

    #[test]
    fn test_weight_transfer_static() {
        let wt = compute_weight_transfer(1000.0, 2.5, 1.5, 0.5, 0.5, 0.0, 0.0);
        let total = wt.total_load();
        let expected = 1000.0 * 9.81;
        assert!(
            (total - expected).abs() < 1.0,
            "total={total} expected={expected}"
        );
    }

    #[test]
    fn test_weight_transfer_total_load_accessor() {
        let wt = WeightTransfer {
            front_left: 1000.0,
            front_right: 1200.0,
            rear_left: 800.0,
            rear_right: 900.0,
        };
        assert!((wt.total_load() - 3900.0).abs() < 1e-6);
        assert!((wt.front_axle_load() - 2200.0).abs() < 1e-6);
        assert!((wt.rear_axle_load() - 1700.0).abs() < 1e-6);
        assert!((wt.left_load() - 1800.0).abs() < 1e-6);
        assert!((wt.right_load() - 2100.0).abs() < 1e-6);
    }

    #[test]
    fn test_longitudinal_weight_transfer_static() {
        let [front, rear] = longitudinal_weight_transfer(1000.0, 2.5, 0.5, 0.5, 0.0);
        let total = front + rear;
        let expected = 1000.0 * 9.81;
        assert!((total - expected).abs() < 1.0);
        assert!(
            (front - rear).abs() < 1.0,
            "50/50 distribution should be equal"
        );
    }

    #[test]
    fn test_longitudinal_weight_transfer_braking() {
        // Under braking (negative accel), front gets more load
        let [front_static, _] = longitudinal_weight_transfer(1000.0, 2.5, 0.5, 0.5, 0.0);
        let [front_braking, _] = longitudinal_weight_transfer(1000.0, 2.5, 0.5, 0.5, -5.0);
        assert!(
            front_braking > front_static,
            "braking should increase front load"
        );
    }

    #[test]
    fn test_longitudinal_weight_transfer_accelerating() {
        // Under acceleration, rear gets more load
        let [_, rear_static] = longitudinal_weight_transfer(1000.0, 2.5, 0.5, 0.5, 0.0);
        let [_, rear_accel] = longitudinal_weight_transfer(1000.0, 2.5, 0.5, 0.5, 5.0);
        assert!(
            rear_accel > rear_static,
            "acceleration should increase rear load"
        );
    }

    // --- Force balance tests ---

    #[test]
    fn test_net_longitudinal_force() {
        let net = net_longitudinal_force(5000.0, 1000.0, 500.0, 200.0);
        assert!((net - 3300.0).abs() < 1e-6);
    }

    #[test]
    fn test_net_longitudinal_force_braking() {
        let net = net_longitudinal_force(0.0, 5000.0, 200.0, 100.0);
        assert!(net < 0.0, "braking should give negative net force");
    }

    #[test]
    fn test_rolling_resistance_force() {
        let f = rolling_resistance_force(5000.0, 0.015);
        assert!((f - 75.0).abs() < 1e-6);
    }

    #[test]
    fn test_rolling_resistance_zero_load() {
        let f = rolling_resistance_force(0.0, 0.015);
        assert!(f.abs() < 1e-10);
    }

    #[test]
    fn test_max_cornering_speed() {
        let v = max_cornering_speed(100.0, 1.0, 0.0);
        let expected = (9.81_f64 * 100.0).sqrt();
        assert!((v - expected).abs() < 0.01);
    }

    #[test]
    fn test_max_cornering_speed_with_downforce() {
        let v_no_df = max_cornering_speed(100.0, 1.0, 0.0);
        let v_with_df = max_cornering_speed(100.0, 1.0, 0.5);
        assert!(
            v_with_df > v_no_df,
            "downforce should increase cornering speed"
        );
    }

    // --- SimpleVehicleState tests ---

    #[test]
    fn test_simple_vehicle_state_default() {
        let s = SimpleVehicleState::default();
        assert_eq!(s.position, 0.0);
        assert_eq!(s.speed, 0.0);
        assert_eq!(s.acceleration, 0.0);
    }

    #[test]
    fn test_simple_vehicle_state_integrate() {
        let mut s = SimpleVehicleState {
            acceleration: 10.0,
            ..Default::default()
        };
        s.integrate(0.1);
        assert!((s.speed - 1.0).abs() < 1e-10);
        assert!((s.position - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_simple_vehicle_state_no_negative_speed() {
        let mut s = SimpleVehicleState {
            speed: 1.0,
            acceleration: -100.0,
            ..Default::default()
        };
        s.integrate(0.1);
        assert!(s.speed >= 0.0, "speed should not go negative");
    }

    #[test]
    fn test_simple_vehicle_state_apply_force() {
        let mut s = SimpleVehicleState::default();
        s.apply_force(1000.0, 500.0);
        assert!((s.acceleration - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_simple_vehicle_state_kinetic_energy() {
        let s = SimpleVehicleState {
            speed: 10.0,
            ..Default::default()
        };
        let ke = s.kinetic_energy(100.0);
        assert!((ke - 5000.0).abs() < 1e-6);
    }

    #[test]
    fn test_simple_vehicle_state_multi_step() {
        let mut s = SimpleVehicleState::default();
        s.apply_force(500.0, 100.0); // a = 5 m/s^2
        for _ in 0..10 {
            s.integrate(0.1);
        }
        // After 1 second: v = 5.0 m/s, position ~ 2.5 m (Euler approximation)
        assert!((s.speed - 5.0).abs() < 0.1);
        assert!(s.position > 2.0 && s.position < 3.0);
    }

    #[test]
    fn test_simple_vehicle_state_apply_force_zero_mass() {
        let mut s = SimpleVehicleState::default();
        s.apply_force(1000.0, 0.0);
        assert_eq!(s.acceleration, 0.0);
    }

    // --- AntiRollBar tests ---

    #[test]
    fn test_anti_roll_bar_zero_delta_zero_force() {
        let bar = AntiRollBar::default();
        let (fl, fr) = bar.compute_forces(0.05, 0.05);
        assert!(fl.abs() < 1e-10, "equal compression → zero roll force");
        assert!(fr.abs() < 1e-10, "equal compression → zero roll force");
    }

    #[test]
    fn test_anti_roll_bar_unequal_compression() {
        let bar = AntiRollBar::new(20000.0);
        let (fl, fr) = bar.compute_forces(0.1, 0.0);
        assert!(fl < 0.0, "left should be pushed down");
        assert!(fr > 0.0, "right should be pushed up");
        assert!(
            (fl + fr).abs() < 1e-10,
            "forces should be equal and opposite"
        );
    }

    #[test]
    fn test_apply_anti_roll_bar_modifies_forces() {
        let bar = AntiRollBar::new(10000.0);
        let mut forces = [1000.0_f64; 4];
        apply_anti_roll_bar(&bar, 0.1, 0.0, &mut forces, 0, 1);
        assert!(
            (forces[0] + forces[1] - 2000.0).abs() < 1e-6,
            "total axle load preserved"
        );
    }

    // --- AeroModel tests ---

    #[test]
    fn test_aero_model_drag_zero_speed() {
        let aero = AeroModel::default();
        assert!((aero.drag_force(0.0)).abs() < 1e-10);
    }

    #[test]
    fn test_aero_model_downforce_positive() {
        let aero = AeroModel::default();
        let f = aero.downforce(40.0);
        assert!(f > 0.0, "downforce should be positive at speed");
    }

    #[test]
    fn test_aero_model_drag_scales_v_squared() {
        let aero = AeroModel::default();
        let f1 = aero.drag_force(10.0);
        let f2 = aero.drag_force(20.0);
        assert!((f2 / f1 - 4.0).abs() < 0.01, "drag should scale v²");
    }

    #[test]
    fn test_aero_model_pitch_moment() {
        let aero = AeroModel::default();
        let m = aero.pitch_moment(30.0);
        let df = aero.downforce(30.0);
        assert!((m - df * aero.cop_offset).abs() < 1e-6);
    }

    #[test]
    fn test_aero_load_transfer_sums_to_zero() {
        let aero = AeroModel::default();
        let [front, rear] = aero.aero_load_transfer(50.0, 2.5);
        assert!(
            (front + rear).abs() < 1e-10,
            "aero load transfer must sum to zero"
        );
    }

    // --- ChassisFlex tests ---

    #[test]
    fn test_chassis_flex_initial_state() {
        let flex = ChassisFlex::default();
        assert!((flex.deflection).abs() < 1e-12);
        assert!((flex.deflection_rate).abs() < 1e-12);
    }

    #[test]
    fn test_chassis_flex_restoring_torque_zero() {
        let flex = ChassisFlex::default();
        assert!((flex.restoring_torque()).abs() < 1e-12);
    }

    #[test]
    fn test_chassis_flex_deflects_under_torque() {
        let mut flex = ChassisFlex::new(1000.0, 10.0);
        flex.step(500.0, 100.0, 0.01);
        assert!(
            flex.deflection.abs() > 0.0,
            "chassis should deflect under torque"
        );
    }

    #[test]
    fn test_chassis_flex_potential_energy_nonzero_after_deflection() {
        let mut flex = ChassisFlex::new(1000.0, 0.0);
        flex.step(500.0, 100.0, 0.1);
        let pe = flex.potential_energy();
        assert!(pe > 0.0, "deflected spring should have potential energy");
    }

    #[test]
    fn test_chassis_flex_kinetic_energy_nonzero_after_step() {
        let mut flex = ChassisFlex::new(1000.0, 0.0);
        flex.step(500.0, 100.0, 0.1);
        let ke = flex.kinetic_energy(100.0);
        assert!(ke > 0.0, "deflection rate gives kinetic energy");
    }

    #[test]
    fn test_chassis_flex_zero_inertia_no_change() {
        let mut flex = ChassisFlex::new(1000.0, 10.0);
        flex.step(500.0, 0.0, 0.1);
        assert!(
            (flex.deflection).abs() < 1e-12,
            "zero inertia should not change deflection"
        );
    }

    // --- full_wheel_loads tests ---

    #[test]
    fn test_full_wheel_loads_static_total() {
        let aero = AeroModel {
            cl_area: 0.0,
            cd_area: 0.8,
            air_density: 1.225,
            cop_offset: 0.0,
        };
        let loads = full_wheel_loads(1000.0, 2.5, 1.5, 0.5, 0.5, 0.0, 0.0, &aero, 0.0);
        let total: f64 = loads.iter().sum();
        let expected = 1000.0 * 9.81;
        assert!(
            (total - expected).abs() < 1.0,
            "total load={total}, expected={expected}"
        );
    }

    #[test]
    fn test_full_wheel_loads_all_nonnegative() {
        let aero = AeroModel::default();
        let loads = full_wheel_loads(1200.0, 2.4, 1.4, 0.5, 0.5, 0.0, 0.0, &aero, 30.0);
        for l in &loads {
            assert!(*l >= 0.0, "wheel load should be non-negative, got {l}");
        }
    }

    // --- integrate_speed_rk4 tests ---

    #[test]
    fn test_integrate_speed_rk4_constant_force() {
        let speed0 = 0.0_f64;
        let mass = 1000.0;
        let force = 2000.0; // a = 2 m/s^2
        let dt = 0.1;
        let speed1 = integrate_speed_rk4(speed0, mass, dt, |_v| force);
        assert!((speed1 - 0.2).abs() < 0.01, "speed1={speed1}");
    }

    #[test]
    fn test_integrate_speed_rk4_speed_bounded_positive() {
        let speed1 = integrate_speed_rk4(0.0, 1000.0, 0.1, |_| -10000.0);
        assert!(speed1 >= 0.0, "speed should not go negative");
    }

    #[test]
    fn test_integrate_speed_rk4_terminal_velocity() {
        // drag force = -k*v, engine = F0 → terminal at v = F0/k
        let f0 = 5000.0;
        let k = 50.0;
        let mass = 1000.0;
        let mut speed = 0.0_f64;
        for _ in 0..10000 {
            speed = integrate_speed_rk4(speed, mass, 0.01, |v| f0 - k * v);
        }
        let expected = f0 / k; // 100 m/s
        assert!(
            (speed - expected).abs() < 1.0,
            "terminal speed={speed}, expected={expected}"
        );
    }

    // --- braking_distance tests ---

    #[test]
    fn test_braking_distance_formula() {
        let d = braking_distance(20.0, 1.0); // 20 m/s, mu=1.0 → d = 400/19.62 ≈ 20.4 m
        let expected = 20.0 * 20.0 / (2.0 * 9.81);
        assert!((d - expected).abs() < 0.01, "braking distance={d}");
    }

    #[test]
    fn test_braking_distance_zero_speed() {
        assert!((braking_distance(0.0, 1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_braking_distance_higher_mu_shorter() {
        let d1 = braking_distance(30.0, 0.5);
        let d2 = braking_distance(30.0, 1.0);
        assert!(d2 < d1, "higher friction → shorter braking distance");
    }

    // --- traction_limited_acceleration ---

    #[test]
    fn test_traction_limited_acceleration_basic() {
        let loads = [2500.0_f64; 4];
        let a = traction_limited_acceleration(&loads, 1.0, 1000.0);
        let expected = 10000.0 / 1000.0;
        assert!((a - expected).abs() < 1e-6);
    }

    #[test]
    fn test_traction_limited_acceleration_zero_mass() {
        let loads = [1000.0_f64; 4];
        let a = traction_limited_acceleration(&loads, 1.0, 0.0);
        assert!((a).abs() < 1e-12);
    }
}
