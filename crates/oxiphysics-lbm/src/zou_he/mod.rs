//! Auto-generated module structure
//!
//! The D3Q19 and D3Q27 Zou-He conditions (`zou_he_inlet_left_d3q19`,
//! `zou_he_top_face_d3q19`, `zou_he_velocity_inlet_d3q27`, ...) are thin
//! wrappers over the face-generic Hecht-Harting closure implemented in
//! [`face`]; the D2Q9 conditions and the outlet / sponge helpers live in
//! [`functions`] and [`functions_2`].

pub mod face;
pub mod functions;
pub mod functions_2;

// Re-export all types
pub use face::*;
pub use functions::*;
pub use functions_2::*;

#[cfg(test)]
mod tests_faces;
