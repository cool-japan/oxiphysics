// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Face-generic Zou-He boundary conditions for the D3Q19 and D3Q27 lattices.
//!
//! The textbook Zou-He construction (bounce-back of the non-equilibrium part
//! of the population normal to the wall) closes the D2Q9 lattice exactly, but
//! on a three-dimensional lattice it leaves the tangential momentum of the
//! boundary node undetermined.  Hecht and Harting (J. Stat. Mech. (2010)
//! P01018) repair this by letting the unknown populations carry exactly the
//! tangential momentum that the known populations do not.  This module
//! implements that on-site reconstruction once, for every axis-aligned face of
//! a box-shaped domain and for both 3D lattices, driven entirely by the
//! lattice tables (`D3Q19_VELOCITIES`, `D3Q19_WEIGHTS`, `D3Q19_OPPOSITES` and
//! their D3Q27 counterparts).  The per-face functions of the Zou-He module
//! (`zou_he_inlet_left_d3q19`, `zou_he_top_face_d3q19`, ...) are thin wrappers
//! around the functions defined here.
//!
//! # Population groups
//!
//! For a face with inward unit normal `s e_a` (`s = +1` or `-1`, `a` the axis
//! of the face) the populations of a boundary node are split by the sign of
//! `s c_{i,a}`:
//!
//! * `T` (tangent, `s c_{i,a} = 0`) and `O` (outgoing, `s c_{i,a} = -1`) are
//!   known after streaming;
//! * `I` (incoming, `s c_{i,a} = +1`) are unknown, because they would have
//!   streamed in from outside the domain.
//!
//! The lattice opposite `ī` of every unknown `i` is a known member of `O`.  It
//! is never a mirror image of `i` across a coordinate plane.
//!
//! # Closure
//!
//! Summing the populations and their momentum along the face normal gives
//! `ρ (1 - s u_a) = Σ_T f + 2 Σ_O f`.  This fixes `ρ` when the velocity is
//! prescribed and `u_a` when the density is prescribed.  The unknowns are then
//!
//! ```text
//! f_i = f_ī + w_i Σ_d c_{i,d} (ρ u_d - M_d) / W_d
//! M_d = Σ_{j in T} c_{j,d} f_j          W_d = Σ_{j in I} w_j c_{j,d}²
//! ```
//!
//! where `d` runs over the three Cartesian axes and `M_a = 0` for the normal
//! axis.  Lattice opposites cancel in every momentum sum (`c_ī = -c_i`), so
//! each unknown only adds the momentum `w_i c_{i,d} g_d` with
//! `g_d = (ρ u_d - M_d) / W_d` beyond what its opposite removes.  The node
//! therefore carries momentum `M_d + (ρ u_d - M_d) = ρ u_d` along every axis
//! and, with the closure above, mass `ρ`: `Σ f = ρ` and `Σ f c = ρ u` hold to
//! rounding error for *any* known populations, equilibrium or not.  For an
//! equilibrium state with the prescribed `ρ` and `u` the unknowns are returned
//! unchanged.  On D3Q19 the formula is the Hecht-Harting correction: the factor
//! `w_i c_{i,t} / W_t` is `±1/2` and `M_t / 2` is their `N_t`.  On D3Q27 the
//! same weight-proportional distribution extends the correction to the corner
//! populations.
//!
//! All quantities are in lattice units (`Δx = Δt = 1`, `c_s² = 1/3`).  A
//! prescribed velocity must satisfy `|s u_a| < 1` and a prescribed density
//! must be non-zero, otherwise the division in the density / velocity closure
//! is singular and the result is not finite.

use crate::lattice::{
    D3Q19_OPPOSITES, D3Q19_VELOCITIES, D3Q19_WEIGHTS, D3Q27_OPPOSITES, D3Q27_VELOCITIES,
    D3Q27_WEIGHTS,
};

/// One of the six axis-aligned faces of a box-shaped lattice domain.
///
/// The names follow the rest of the Zou-He module: `Left` and `Right` are the
/// `x = 0` and `x = nx - 1` faces, `Bottom` and `Top` the `y = 0` and
/// `y = ny - 1` faces, `Back` and `Front` the `z = 0` and `z = nz - 1` faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ZouHeFace {
    /// The `x = 0` face; its inward normal is `+x`.
    Left,
    /// The `x = nx - 1` face; its inward normal is `-x`.
    Right,
    /// The `y = 0` face; its inward normal is `+y`.
    Bottom,
    /// The `y = ny - 1` face; its inward normal is `-y`.
    Top,
    /// The `z = 0` face; its inward normal is `+z`.
    Back,
    /// The `z = nz - 1` face; its inward normal is `-z`.
    Front,
}

impl ZouHeFace {
    /// All six faces, ordered `Left, Right, Bottom, Top, Back, Front`.
    pub const ALL: [ZouHeFace; 6] = [
        ZouHeFace::Left,
        ZouHeFace::Right,
        ZouHeFace::Bottom,
        ZouHeFace::Top,
        ZouHeFace::Back,
        ZouHeFace::Front,
    ];

    /// Cartesian axis normal to the face (`0` = x, `1` = y, `2` = z).
    pub const fn axis(self) -> usize {
        match self {
            ZouHeFace::Left | ZouHeFace::Right => 0,
            ZouHeFace::Bottom | ZouHeFace::Top => 1,
            ZouHeFace::Back | ZouHeFace::Front => 2,
        }
    }

    /// Sign of the inward-pointing face normal along [`axis`](Self::axis).
    ///
    /// It is `+1` for `Left`, `Bottom` and `Back` (the domain lies at larger
    /// coordinates) and `-1` for `Right`, `Top` and `Front`.
    pub const fn inward_sign(self) -> i32 {
        match self {
            ZouHeFace::Left | ZouHeFace::Bottom | ZouHeFace::Back => 1,
            ZouHeFace::Right | ZouHeFace::Top | ZouHeFace::Front => -1,
        }
    }

    /// The two Cartesian axes tangent to the face, in increasing order.
    pub const fn tangential_axes(self) -> [usize; 2] {
        match self.axis() {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        }
    }

    /// Position of the face in the per-lattice layout tables.
    const fn slot(self) -> usize {
        match self {
            ZouHeFace::Left => 0,
            ZouHeFace::Right => 1,
            ZouHeFace::Bottom => 2,
            ZouHeFace::Top => 3,
            ZouHeFace::Back => 4,
            ZouHeFace::Front => 5,
        }
    }
}

/// How the `Q` populations of a lattice are grouped at one face.
///
/// Each index list is stored in a fixed-size array together with the number of
/// valid entries, so the whole table can be built by a `const fn`.
#[derive(Debug, Clone, Copy)]
struct FaceLayout<const Q: usize> {
    /// Known populations tangent to the face (`s c_a = 0`).
    tangent: [usize; Q],
    n_tangent: usize,
    /// Known populations pointing out of the domain (`s c_a = -1`).
    outgoing: [usize; Q],
    n_outgoing: usize,
    /// Unknown populations pointing into the domain (`s c_a = +1`).
    incoming: [usize; Q],
    n_incoming: usize,
    /// `W_d = Σ_{i incoming} w_i c_{i,d}²` for every Cartesian axis `d`.
    incoming_moment: [f64; 3],
}

impl<const Q: usize> FaceLayout<Q> {
    /// Indices of the tangent populations.
    fn tangent(&self) -> &[usize] {
        &self.tangent[..self.n_tangent]
    }

    /// Indices of the known outgoing populations.
    fn outgoing(&self) -> &[usize] {
        &self.outgoing[..self.n_outgoing]
    }

    /// Indices of the unknown incoming populations.
    fn incoming(&self) -> &[usize] {
        &self.incoming[..self.n_incoming]
    }

    /// `Σ_T f + 2 Σ_O f`, the quantity that equals `ρ (1 - s u_a)`.
    fn known_mass(&self, f: &[f64; Q]) -> f64 {
        let tangent: f64 = self.tangent().iter().map(|&i| f[i]).sum();
        let outgoing: f64 = self.outgoing().iter().map(|&i| f[i]).sum();
        tangent + 2.0 * outgoing
    }
}

/// Classify the populations of a lattice for the face with normal axis `axis`
/// and inward sign `inward` (`+1` or `-1`).
const fn build_layout<const Q: usize>(
    velocities: &[[i32; 3]; Q],
    weights: &[f64; Q],
    axis: usize,
    inward: i32,
) -> FaceLayout<Q> {
    let mut layout = FaceLayout {
        tangent: [0; Q],
        n_tangent: 0,
        outgoing: [0; Q],
        n_outgoing: 0,
        incoming: [0; Q],
        n_incoming: 0,
        incoming_moment: [0.0; 3],
    };
    let mut i = 0;
    while i < Q {
        let along_inward_normal = velocities[i][axis] * inward;
        if along_inward_normal == 0 {
            layout.tangent[layout.n_tangent] = i;
            layout.n_tangent += 1;
        } else if along_inward_normal < 0 {
            layout.outgoing[layout.n_outgoing] = i;
            layout.n_outgoing += 1;
        } else {
            layout.incoming[layout.n_incoming] = i;
            layout.n_incoming += 1;
            let mut d = 0;
            while d < 3 {
                let c = velocities[i][d] as f64;
                layout.incoming_moment[d] += weights[i] * c * c;
                d += 1;
            }
        }
        i += 1;
    }
    layout
}

/// Lattice tables plus the per-face population layouts of a 3D lattice.
struct FaceLattice<const Q: usize> {
    velocities: &'static [[i32; 3]; Q],
    weights: &'static [f64; Q],
    opposites: &'static [usize; Q],
    /// Indexed by [`ZouHeFace::slot`].
    layouts: [FaceLayout<Q>; 6],
}

impl<const Q: usize> FaceLattice<Q> {
    const fn new(
        velocities: &'static [[i32; 3]; Q],
        weights: &'static [f64; Q],
        opposites: &'static [usize; Q],
    ) -> Self {
        FaceLattice {
            velocities,
            weights,
            opposites,
            layouts: [
                build_layout(velocities, weights, 0, 1),
                build_layout(velocities, weights, 0, -1),
                build_layout(velocities, weights, 1, 1),
                build_layout(velocities, weights, 1, -1),
                build_layout(velocities, weights, 2, 1),
                build_layout(velocities, weights, 2, -1),
            ],
        }
    }

    /// Overwrite the unknown populations of `face` so that the node carries
    /// density `rho` and velocity `u`.
    ///
    /// Only populations of the incoming group are written; they are computed
    /// from the known populations alone.
    fn reconstruct(&self, f: &mut [f64; Q], face: ZouHeFace, rho: f64, u: [f64; 3]) {
        let layout = &self.layouts[face.slot()];
        let axis = face.axis();
        // Per axis: the momentum the incoming populations must deliver on top
        // of what the known populations carry, per unit weight moment.
        let drive: [f64; 3] = std::array::from_fn(|d| {
            let carried: f64 = if d == axis {
                0.0
            } else {
                layout
                    .tangent()
                    .iter()
                    .map(|&j| f64::from(self.velocities[j][d]) * f[j])
                    .sum()
            };
            (rho * u[d] - carried) / layout.incoming_moment[d]
        });
        for &i in layout.incoming() {
            let transfer: f64 = self.velocities[i]
                .iter()
                .zip(drive.iter())
                .map(|(&c, &g)| f64::from(c) * g)
                .sum();
            f[i] = f[self.opposites[i]] + self.weights[i] * transfer;
        }
    }

    /// Velocity condition: derive `ρ`, reconstruct, return `ρ`.
    fn impose_velocity(&self, f: &mut [f64; Q], face: ZouHeFace, u: [f64; 3]) -> f64 {
        let layout = &self.layouts[face.slot()];
        let inward_speed = f64::from(face.inward_sign()) * u[face.axis()];
        let rho = layout.known_mass(f) / (1.0 - inward_speed);
        self.reconstruct(f, face, rho, u);
        rho
    }

    /// Pressure condition: derive the normal velocity, reconstruct, return it.
    fn impose_pressure(
        &self,
        f: &mut [f64; Q],
        face: ZouHeFace,
        rho: f64,
        u_tangential: [f64; 2],
    ) -> f64 {
        let layout = &self.layouts[face.slot()];
        let axis = face.axis();
        let inward_speed = 1.0 - layout.known_mass(f) / rho;
        let [first, second] = face.tangential_axes();
        let mut u = [0.0; 3];
        u[axis] = f64::from(face.inward_sign()) * inward_speed;
        u[first] = u_tangential[0];
        u[second] = u_tangential[1];
        self.reconstruct(f, face, rho, u);
        u[axis]
    }
}

static D3Q19_FACES: FaceLattice<19> =
    FaceLattice::new(&D3Q19_VELOCITIES, &D3Q19_WEIGHTS, &D3Q19_OPPOSITES);

static D3Q27_FACES: FaceLattice<27> =
    FaceLattice::new(&D3Q27_VELOCITIES, &D3Q27_WEIGHTS, &D3Q27_OPPOSITES);

/// Zou-He velocity condition on one face of a D3Q19 domain.
///
/// The five populations pointing into the domain are rebuilt from their lattice
/// opposites with the Hecht-Harting (2010) on-site closure so that the node
/// carries the prescribed velocity `u = (ux, uy, uz)`: afterwards `Σ f = ρ`
/// and `Σ f c = ρ u` hold for all three components, whatever the known
/// populations were.  The known populations are not modified and the unknown
/// ones are not read.
///
/// # Arguments
/// * `f`    – the 19 populations of the boundary node after streaming.
/// * `face` – the face the node lies on.
/// * `u`    – prescribed wall / inlet velocity (all three components); the
///   component normal to the face must satisfy `|u_n| < 1`.
///
/// # Returns
/// The density `ρ` of the boundary node.
pub fn zou_he_velocity_face_d3q19(f: &mut [f64; 19], face: ZouHeFace, u: [f64; 3]) -> f64 {
    D3Q19_FACES.impose_velocity(f, face, u)
}

/// Zou-He pressure (density) condition on one face of a D3Q19 domain.
///
/// The density `rho` is imposed, the velocity component normal to the face
/// follows from mass conservation, and the tangential velocity is the given
/// `u_tangential`.  The five populations pointing into the domain are rebuilt
/// as in [`zou_he_velocity_face_d3q19`], so afterwards `Σ f = rho` and
/// `Σ f c = rho u` hold for all three components.
///
/// # Arguments
/// * `f`            – the 19 populations of the boundary node after streaming.
/// * `face`         – the face the node lies on.
/// * `rho`          – prescribed density (pressure `p = ρ / 3`), non-zero.
/// * `u_tangential` – prescribed tangential velocity components, ordered by
///   increasing axis index (`(uy, uz)` on x-faces, `(ux, uz)` on y-faces,
///   `(ux, uy)` on z-faces); `[0.0, 0.0]` is the classical Zou-He assumption
///   of purely normal flow through the face.
///
/// # Returns
/// The velocity component along the positive coordinate axis normal to the
/// face (`ux` on x-faces, ...).  It is negative when fluid leaves through a
/// `Left` / `Bottom` / `Back` face and positive for `Right` / `Top` / `Front`.
pub fn zou_he_pressure_face_d3q19(
    f: &mut [f64; 19],
    face: ZouHeFace,
    rho: f64,
    u_tangential: [f64; 2],
) -> f64 {
    D3Q19_FACES.impose_pressure(f, face, rho, u_tangential)
}

/// Zou-He velocity condition on one face of a D3Q27 domain.
///
/// The nine populations pointing into the domain are rebuilt from their lattice
/// opposites with the same weight-proportional closure as the D3Q19 version
/// ([`zou_he_velocity_face_d3q19`]), so afterwards `Σ f = ρ` and `Σ f c = ρ u`
/// hold for all three components of the prescribed velocity `u`.
///
/// # Arguments
/// * `f`    – the 27 populations of the boundary node after streaming.
/// * `face` – the face the node lies on.
/// * `u`    – prescribed velocity; the normal component must satisfy
///   `|u_n| < 1`.
///
/// # Returns
/// The density `ρ` of the boundary node.
pub fn zou_he_velocity_face_d3q27(f: &mut [f64; 27], face: ZouHeFace, u: [f64; 3]) -> f64 {
    D3Q27_FACES.impose_velocity(f, face, u)
}

/// Zou-He pressure (density) condition on one face of a D3Q27 domain.
///
/// Counterpart of [`zou_he_pressure_face_d3q19`] for the 27-velocity lattice;
/// argument and return conventions are identical.
pub fn zou_he_pressure_face_d3q27(
    f: &mut [f64; 27],
    face: ZouHeFace,
    rho: f64,
    u_tangential: [f64; 2],
) -> f64 {
    D3Q27_FACES.impose_pressure(f, face, rho, u_tangential)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Check the population grouping of every face of `lattice`; `sizes` is the
    /// expected `(tangent, outgoing, incoming)` population count.
    fn assert_layouts<const Q: usize>(lattice: &FaceLattice<Q>, sizes: (usize, usize, usize)) {
        for face in ZouHeFace::ALL {
            let layout = &lattice.layouts[face.slot()];
            let tangent = layout.tangent();
            let outgoing = layout.outgoing();
            let incoming = layout.incoming();
            assert_eq!(
                (tangent.len(), outgoing.len(), incoming.len()),
                sizes,
                "{face:?}: group sizes"
            );
            // The three groups partition the lattice.
            let mut seen = [0usize; Q];
            for &i in tangent.iter().chain(outgoing).chain(incoming) {
                seen[i] += 1;
            }
            assert!(
                seen.iter().all(|&count| count == 1),
                "{face:?}: groups do not partition the populations"
            );
            // Membership follows the velocity component along the inward normal.
            let axis = face.axis();
            let sign = face.inward_sign();
            for &i in tangent {
                assert_eq!(lattice.velocities[i][axis], 0, "{face:?}: tangent {i}");
            }
            for &i in outgoing {
                assert_eq!(
                    lattice.velocities[i][axis] * sign,
                    -1,
                    "{face:?}: outgoing {i}"
                );
            }
            for &i in incoming {
                assert_eq!(
                    lattice.velocities[i][axis] * sign,
                    1,
                    "{face:?}: incoming {i}"
                );
                assert!(
                    outgoing.contains(&lattice.opposites[i]),
                    "{face:?}: the opposite of unknown {i} must be a known outgoing population"
                );
            }
            // Weight moments of the incoming group: 1/6 along the normal axis
            // and 1/18 along the tangential axes, for both 3D lattices.
            for (d, &moment) in layout.incoming_moment.iter().enumerate() {
                let expected = if d == axis { 1.0 / 6.0 } else { 1.0 / 18.0 };
                assert!(
                    (moment - expected).abs() < 1.0e-15,
                    "{face:?}: weight moment along axis {d} is {moment}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn d3q19_face_layouts_partition_the_populations() {
        assert_layouts(&D3Q19_FACES, (9, 5, 5));
    }

    #[test]
    fn d3q27_face_layouts_partition_the_populations() {
        assert_layouts(&D3Q27_FACES, (9, 9, 9));
    }
}
