// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vehicle dynamics simulation for the OxiPhysics engine.
//!
//! This crate provides a raycast-based vehicle simulation including:
//!
//! - **Wheel**: wheel configuration and runtime state
//! - **Suspension**: linear and progressive spring-damper models
//! - **Tire**: Pacejka (Magic Formula), Fiala (brush), and linear tire models
//! - **Drivetrain**: engine torque curves, gearbox, differentials
//! - **Steering**: Ackermann steering geometry
//! - **Vehicle**: full raycast vehicle combining all subsystems
//! - **Stability**: ABS, traction control, and electronic stability control
//!
//! # Usage
//!
//! ```no_run
//! use oxiphysics_vehicle::vehicle::RaycastVehicle;
//!
//! let vehicle = RaycastVehicle::new_default_4wheel();
//! assert_eq!(vehicle.num_wheels(), 4);
//! ```
#![allow(missing_docs)]
#![allow(ambiguous_glob_reexports)]
#![allow(dead_code)]

mod error;
pub use error::*;

pub mod drivetrain;
pub mod steering;
pub mod suspension;
pub mod tire;
pub mod vehicle;
pub mod wheel;

pub use drivetrain::{Differential, DriveLayout, Drivetrain, EngineCurve, Gearbox};
pub use steering::AckermannSteering;
pub use suspension::{LinearSuspension, ProgressiveSuspension, SuspensionModel};
pub use tire::{FialaTire, LinearTire, PacejkaCoeffs, PacejkaTire, TireModel};
pub use vehicle::{AeroDrag, GroundHit, RaycastVehicle, VehicleState, flat_ground_query};
pub use wheel::{Wheel, WheelConfig, WheelState};

pub mod stability;
pub use stability::{AbsController, ElectronicStabilityControl, TractionControl};

pub mod aerodynamics;
pub use aerodynamics::*;

pub mod telemetry;

pub mod path;

pub mod race_line;

pub mod sensors;

pub mod tire_thermal;
pub use tire_thermal::*;

pub mod path_planning;
pub use path_planning::*;
pub mod active_suspension;
pub mod aerodynamics_vehicle;
pub mod aircraft_dynamics;
pub mod autonomous;
pub mod autonomous_driving;
pub mod autonomous_vehicle;
pub mod brake_system;
pub mod charging_station;
pub mod chassis_dynamics;
pub mod chassis_resonance;
pub mod cooling_system;
pub mod driver_model;
pub mod electric_motor;
pub mod electric_vehicle;
pub mod energy_recovery;
pub mod fuel_cell;
pub mod fuel_system;
pub mod fuel_systems;
pub mod gearbox;
pub mod hvac;
pub mod lap_simulator;
pub mod motorcycle_dynamics;
pub mod noise_vibration;
pub mod powertrain_advanced;
pub mod powertrain_dynamics;
pub mod race_simulation;
pub mod race_strategy;
pub mod ride_quality;
pub mod suspension_analysis;
pub mod suspension_optimization;
pub mod thermal_management;
pub mod tire_wear;
pub mod track_dynamics;
pub mod vehicle_dynamics;
pub mod wind_loading;

pub mod hil;
pub use hil::{
    HilChannel, HilConfig, HilError, HilInterface, HilSignalLogger, HilTimingStats, SimHilBridge,
};

pub mod fem_chassis_cosim;
pub use fem_chassis_cosim::{
    ChassisAttachmentPoint, ChassisCosimBridge, CosimLogEntry, ModalChassisModel, bending_stress,
    populate_beam_modes,
};

pub mod gpu_multi_vehicle;
pub use gpu_multi_vehicle::{
    BatchStats, MultiVehicleBatch, VehicleInput, VehicleParams, VehicleSoaState,
};
