//! Typed errors for URDF parsing and `ArticulatedModel` mapping.

/// Typed, non-exhaustive errors for URDF parsing and `ArticulatedModel` mapping.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum UrdfError {
    /// The XML document could not be tokenised/parsed.
    #[error("XML parse error: {0}")]
    XmlParse(String),
    /// A required attribute was missing on an element.
    #[error("Missing required attribute '{attr}' on element '{element}'")]
    MissingAttribute {
        /// Element on which the attribute was expected.
        element: String,
        /// The missing attribute name.
        attr: String,
    },
    /// An attribute value could not be interpreted.
    #[error("Invalid value '{value}' for attribute '{attr}': {reason}")]
    InvalidValue {
        /// Attribute name.
        attr: String,
        /// The offending value.
        value: String,
        /// Why it is invalid.
        reason: String,
    },
    /// More than one root link was found in the kinematic tree.
    #[error("Multiple root links found: '{0}' and '{1}'")]
    MultipleRoots(String, String),
    /// No root link was found in the URDF.
    #[error("No root link found in URDF")]
    NoRoot,
    /// A cycle was detected in the kinematic tree.
    #[error("Cycle detected in kinematic tree involving link '{0}'")]
    KinematicCycle(String),
    /// A joint referenced a link that does not exist.
    #[error("Unknown link reference '{0}' in joint '{1}'")]
    UnknownLink(String, String),
    /// The joint type is not supported by the articulated mapping.
    #[error("Unsupported joint type: {0}")]
    UnsupportedJointType(String),
    /// The inertia tensor was not positive definite.
    #[error("Inertia matrix is not positive definite for link '{0}'")]
    InvalidInertia(String),
    /// A writer error occurred.
    #[error("Writer error: {0}")]
    Write(String),
}

/// Result type for URDF operations.
pub type UrdfResult<T> = Result<T, UrdfError>;
