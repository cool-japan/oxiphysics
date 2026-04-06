//! # AngleCV - Trait Implementations
//!
//! This module contains trait implementations for `AngleCV`.
//!
//! ## Implemented Traits
//!
//! - `CollectiveVariable`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::CollectiveVariable;
#[allow(unused_imports)]
use super::functions::*;
use super::types::AngleCV;

impl CollectiveVariable for AngleCV {
    fn value(&self, positions: &[[f64; 3]]) -> f64 {
        let ri = positions[self.atom_i];
        let rj = positions[self.atom_j];
        let rk = positions[self.atom_k];
        let u = [ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]];
        let v = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
        let u_norm = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        let v_norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if u_norm < 1e-15 || v_norm < 1e-15 {
            return 0.0;
        }
        let cos_theta = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (u_norm * v_norm);
        let cos_theta = cos_theta.clamp(-1.0, 1.0);
        cos_theta.acos()
    }
    fn gradient(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let n = positions.len();
        let mut grad = vec![[0.0f64; 3]; n];
        let ri = positions[self.atom_i];
        let rj = positions[self.atom_j];
        let rk = positions[self.atom_k];
        let u = [ri[0] - rj[0], ri[1] - rj[1], ri[2] - rj[2]];
        let v = [rk[0] - rj[0], rk[1] - rj[1], rk[2] - rj[2]];
        let u_norm = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
        let v_norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if u_norm < 1e-15 || v_norm < 1e-15 {
            return grad;
        }
        let cos_theta = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (u_norm * v_norm);
        let cos_theta = cos_theta.clamp(-1.0, 1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-15);
        let u_hat = [u[0] / u_norm, u[1] / u_norm, u[2] / u_norm];
        let v_hat = [v[0] / v_norm, v[1] / v_norm, v[2] / v_norm];
        let dtheta_du: [f64; 3] = [
            -(v_hat[0] - cos_theta * u_hat[0]) / (u_norm * sin_theta),
            -(v_hat[1] - cos_theta * u_hat[1]) / (u_norm * sin_theta),
            -(v_hat[2] - cos_theta * u_hat[2]) / (u_norm * sin_theta),
        ];
        let dtheta_dv: [f64; 3] = [
            -(u_hat[0] - cos_theta * v_hat[0]) / (v_norm * sin_theta),
            -(u_hat[1] - cos_theta * v_hat[1]) / (v_norm * sin_theta),
            -(u_hat[2] - cos_theta * v_hat[2]) / (v_norm * sin_theta),
        ];
        grad[self.atom_i] = dtheta_du;
        grad[self.atom_k] = dtheta_dv;
        grad[self.atom_j] = [
            -dtheta_du[0] - dtheta_dv[0],
            -dtheta_du[1] - dtheta_dv[1],
            -dtheta_du[2] - dtheta_dv[2],
        ];
        grad
    }
    fn name(&self) -> &str {
        &self.name
    }
}
