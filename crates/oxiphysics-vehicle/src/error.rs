// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Error types for oxiphysics-vehicle.

use thiserror::Error;

/// Main error type for the vehicle module.
#[derive(Debug, Error)]
pub enum VehicleError {
    /// Generic error with a message.
    #[error("{0}")]
    General(String),

    /// Invalid wheel index.
    #[error("invalid wheel index {index}: vehicle has {count} wheels")]
    InvalidWheelIndex {
        /// The requested index.
        index: usize,
        /// Total number of wheels.
        count: usize,
    },

    /// Invalid gear.
    #[error("invalid gear {gear}: valid range is -1..{max_gear}")]
    InvalidGear {
        /// The requested gear.
        gear: i32,
        /// Maximum forward gear number.
        max_gear: usize,
    },

    /// Invalid parameter value.
    #[error("invalid parameter '{name}': {reason}")]
    InvalidParameter {
        /// Parameter name.
        name: String,
        /// Reason the value is invalid.
        reason: String,
    },
}

/// Result type alias for vehicle operations.
pub type VehicleResult<T> = std::result::Result<T, VehicleError>;

// ──────────────────────────────────────────────────────────────────────────────
// Additional error variants for sensor and telemetry subsystems
// ──────────────────────────────────────────────────────────────────────────────

/// Errors specific to sensor subsystem operations.
#[derive(Debug, thiserror::Error)]
pub enum SensorError {
    /// Sensor configuration value is out of its valid range.
    #[error("sensor '{sensor}' parameter '{param}' out of range: value={value}")]
    OutOfRange {
        /// Sensor name or type identifier.
        sensor: String,
        /// Name of the parameter that is invalid.
        param: String,
        /// The offending value.
        value: f64,
    },
    /// Sensor returned no data (e.g. out of range, below SNR threshold).
    #[error("sensor '{sensor}' no data: {reason}")]
    NoData {
        /// Sensor name.
        sensor: String,
        /// Reason no data was produced.
        reason: String,
    },
    /// Calibration data is missing or corrupt.
    #[error("sensor '{sensor}' calibration error: {reason}")]
    CalibrationError {
        /// Sensor name.
        sensor: String,
        /// Details about the calibration failure.
        reason: String,
    },
}

/// Result type alias for sensor operations.
pub type SensorResult<T> = std::result::Result<T, SensorError>;

/// Errors specific to the telemetry recording subsystem.
#[derive(Debug, thiserror::Error)]
pub enum TelemetryError {
    /// Attempt to access a frame index that is out of bounds.
    #[error("frame index {index} out of bounds (total frames: {count})")]
    FrameIndexOutOfBounds {
        /// Requested index.
        index: usize,
        /// Total recorded frames.
        count: usize,
    },
    /// Recorder buffer is full and cannot accept more data.
    #[error("telemetry buffer full (capacity: {capacity})")]
    BufferFull {
        /// Maximum capacity of the buffer.
        capacity: usize,
    },
    /// Data export failed.
    #[error("telemetry export failed: {reason}")]
    ExportFailed {
        /// Reason for export failure.
        reason: String,
    },
    /// Channel with the given name was not found.
    #[error("telemetry channel '{name}' not found")]
    ChannelNotFound {
        /// Channel name that was requested.
        name: String,
    },
}

/// Result type alias for telemetry operations.
pub type TelemetryResult<T> = std::result::Result<T, TelemetryError>;

/// Errors relating to path and trajectory operations.
#[derive(Debug, thiserror::Error)]
pub enum PathError {
    /// Path has too few waypoints for the requested operation.
    #[error("path has {count} waypoints; need at least {required}")]
    TooFewWaypoints {
        /// Number of waypoints in the path.
        count: usize,
        /// Minimum required.
        required: usize,
    },
    /// Arc-length parameter is outside `[0, total_length]`.
    #[error("arc-length {s:.3} out of range [0, {total:.3}]")]
    ArcLengthOutOfRange {
        /// The requested arc-length.
        s: f64,
        /// Total path length.
        total: f64,
    },
    /// Spline construction failed.
    #[error("spline construction failed: {reason}")]
    SplineConstructionFailed {
        /// Detailed reason.
        reason: String,
    },
}

/// Result type alias for path operations.
pub type PathResult<T> = std::result::Result<T, PathError>;

// ──────────────────────────────────────────────────────────────────────────────
// Helper constructors
// ──────────────────────────────────────────────────────────────────────────────

impl VehicleError {
    /// Construct a [`VehicleError::General`] from any displayable value.
    #[allow(dead_code)]
    pub fn general(msg: impl std::fmt::Display) -> Self {
        VehicleError::General(msg.to_string())
    }

    /// Construct an [`VehicleError::InvalidParameter`].
    #[allow(dead_code)]
    pub fn invalid_param(name: impl Into<String>, reason: impl Into<String>) -> Self {
        VehicleError::InvalidParameter {
            name: name.into(),
            reason: reason.into(),
        }
    }

    /// Return `true` if this is a parameter validation error.
    #[allow(dead_code)]
    pub fn is_invalid_parameter(&self) -> bool {
        matches!(self, VehicleError::InvalidParameter { .. })
    }
}

impl SensorError {
    /// Construct a [`SensorError::NoData`].
    #[allow(dead_code)]
    pub fn no_data(sensor: impl Into<String>, reason: impl Into<String>) -> Self {
        SensorError::NoData {
            sensor: sensor.into(),
            reason: reason.into(),
        }
    }
}

impl PathError {
    /// Construct a [`PathError::TooFewWaypoints`].
    #[allow(dead_code)]
    pub fn too_few(count: usize, required: usize) -> Self {
        PathError::TooFewWaypoints { count, required }
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Tests
// ──────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vehicle_error_general_message() {
        let e = VehicleError::general("something went wrong");
        assert!(e.to_string().contains("something went wrong"));
    }

    #[test]
    fn vehicle_error_invalid_param_fields() {
        let e = VehicleError::invalid_param("mass", "must be positive");
        assert!(e.to_string().contains("mass"));
        assert!(e.to_string().contains("must be positive"));
    }

    #[test]
    fn vehicle_error_is_invalid_parameter() {
        let e = VehicleError::invalid_param("wheelbase", "must be > 0");
        assert!(e.is_invalid_parameter());
    }

    #[test]
    fn vehicle_error_invalid_wheel_index() {
        let e = VehicleError::InvalidWheelIndex { index: 5, count: 4 };
        assert!(e.to_string().contains("5"));
        assert!(e.to_string().contains("4"));
    }

    #[test]
    fn vehicle_error_invalid_gear() {
        let e = VehicleError::InvalidGear {
            gear: -2,
            max_gear: 7,
        };
        let msg = e.to_string();
        assert!(msg.contains("-2"));
        assert!(msg.contains("7"));
    }

    #[test]
    fn sensor_error_no_data_display() {
        let e = SensorError::no_data("radar", "target out of range");
        assert!(e.to_string().contains("radar"));
        assert!(e.to_string().contains("target out of range"));
    }

    #[test]
    fn sensor_error_out_of_range_display() {
        let e = SensorError::OutOfRange {
            sensor: "imu".into(),
            param: "accel_noise".into(),
            value: -1.0,
        };
        assert!(e.to_string().contains("imu"));
        assert!(e.to_string().contains("accel_noise"));
    }

    #[test]
    fn telemetry_error_frame_index_display() {
        let e = TelemetryError::FrameIndexOutOfBounds {
            index: 100,
            count: 50,
        };
        assert!(e.to_string().contains("100"));
        assert!(e.to_string().contains("50"));
    }

    #[test]
    fn telemetry_error_buffer_full_display() {
        let e = TelemetryError::BufferFull { capacity: 1024 };
        assert!(e.to_string().contains("1024"));
    }

    #[test]
    fn telemetry_error_channel_not_found() {
        let e = TelemetryError::ChannelNotFound {
            name: "lateral_g".into(),
        };
        assert!(e.to_string().contains("lateral_g"));
    }

    #[test]
    fn path_error_too_few_waypoints_display() {
        let e = PathError::too_few(1, 3);
        assert!(e.to_string().contains("1"));
        assert!(e.to_string().contains("3"));
    }

    #[test]
    fn path_error_arc_length_out_of_range() {
        let e = PathError::ArcLengthOutOfRange {
            s: 500.0,
            total: 300.0,
        };
        assert!(e.to_string().contains("500"));
        assert!(e.to_string().contains("300"));
    }

    #[test]
    #[allow(clippy::unnecessary_literal_unwrap)]
    fn result_aliases_compile() {
        let v: VehicleResult<i32> = Ok(42);
        assert_eq!(v.unwrap(), 42);
        let s: SensorResult<f64> = Ok(3.125);
        assert!((s.unwrap() - 3.125).abs() < 1e-12);
        let t: TelemetryResult<bool> = Ok(true);
        assert!(t.unwrap());
        let p: PathResult<()> = Ok(());
        assert!(p.is_ok());
    }
}
