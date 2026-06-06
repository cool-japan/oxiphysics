//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Zou-He velocity boundary condition: left inlet wall (x = 0).
///
/// Given a cell on the left boundary with prescribed x-velocity `ux_wall`
/// and zero transverse velocity, this function computes the three unknown
/// distributions `f[1]`, `f[5]`, `f[8]` (pointing into the domain from
/// the left) and updates the density stored in the array slot indirectly
/// through the return value.
///
/// # Arguments
/// * `f`       – mutable reference to the 9-element distribution array.
/// * `ux_wall` – prescribed inlet velocity in the +x direction.
///
/// # Returns
/// The computed density `ρ` at this boundary cell.
pub fn zou_he_inlet_left(f: &mut [f64; 9], ux_wall: f64) -> f64 {
    let rho = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])) / (1.0 - ux_wall);
    f[1] = f[3] + (2.0 / 3.0) * rho * ux_wall;
    f[5] = f[7] - 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * rho * ux_wall;
    f[8] = f[6] + 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * rho * ux_wall;
    rho
}
/// Zou-He pressure (density) boundary condition: right outlet wall (x = nx-1).
///
/// Given a cell on the right boundary with prescribed density `rho_wall`
/// (corresponding to the desired outlet pressure `p = cs² ρ`), this function
/// computes the three unknown distributions `f[3]`, `f[6]`, `f[7]` (pointing
/// back into the domain from the right) and returns the computed x-velocity.
///
/// # Arguments
/// * `f`        – mutable reference to the 9-element distribution array.
/// * `rho_wall` – prescribed density at the outlet.
///
/// # Returns
/// The computed outflow velocity `ux` at this boundary cell.
pub fn zou_he_outlet_right(f: &mut [f64; 9], rho_wall: f64) -> f64 {
    let ux = -1.0 + (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / rho_wall;
    f[3] = f[1] - (2.0 / 3.0) * rho_wall * ux;
    f[7] = f[5] + 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * rho_wall * ux;
    f[6] = f[8] - 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * rho_wall * ux;
    ux
}
/// Zou-He pressure inlet boundary condition: left wall (x = 0), D2Q9.
///
/// Sets the density `rho_in` at the left boundary and computes the
/// three unknown distributions `f[1]`, `f[5]`, `f[8]` from the
/// no-penetration plus zero-shear conditions (uy = 0 assumed).
///
/// # Returns
/// The computed inlet velocity `ux`.
pub fn apply_pressure_inlet_d2q9(f: &mut [f64; 9], rho_in: f64) -> f64 {
    let ux = 1.0 - (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])) / rho_in;
    f[1] = f[3] + (2.0 / 3.0) * rho_in * ux;
    f[5] = f[7] - 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * rho_in * ux;
    f[8] = f[6] + 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * rho_in * ux;
    ux
}
/// Zou-He velocity outlet boundary condition: right wall (x = nx-1), D2Q9.
///
/// Sets the x-velocity `ux_out` at the right boundary and computes the
/// three unknown distributions `f[3]`, `f[6]`, `f[7]` from the momentum
/// balance.
///
/// # Returns
/// The computed outlet density `rho`.
pub fn apply_velocity_outlet_d2q9(f: &mut [f64; 9], ux_out: f64) -> f64 {
    let rho = (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / (1.0 + ux_out);
    f[3] = f[1] - (2.0 / 3.0) * rho * ux_out;
    f[7] = f[5] + 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * rho * ux_out;
    f[6] = f[8] - 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * rho * ux_out;
    rho
}
/// Zou-He velocity inlet: left wall (x = 0) of a D3Q19 domain.
///
/// The inlet face is normal to +x.  Unknown distributions are those pointing
/// into the domain in the +x direction: 1, 7, 9, 11, 13.
/// The remaining 14 distributions are known.
///
/// # Arguments
/// * `f`     – mutable reference to the 19-element distribution array.
/// * `ux_in` – prescribed inlet velocity in +x direction.
///
/// # Returns
/// The computed density `ρ` at this boundary cell.
pub fn zou_he_inlet_left_d3q19(f: &mut [f64; 19], ux_in: f64) -> f64 {
    let sum_known = f[0]
        + f[2]
        + f[3]
        + f[4]
        + f[5]
        + f[6]
        + f[8]
        + f[10]
        + f[12]
        + f[14]
        + f[15]
        + f[16]
        + f[17]
        + f[18];
    let sum_zero_x = f[0] + f[3] + f[4] + f[5] + f[6] + f[15] + f[16] + f[17] + f[18];
    let sum_neg_x = f[2] + f[8] + f[10] + f[12] + f[14];
    let rho = (sum_zero_x + 2.0 * sum_neg_x) / (1.0 - ux_in);
    f[1] = f[2] + (2.0 / 3.0) * rho * ux_in;
    f[7] = f[10] + 0.5 * (f[4] - f[3]) + (1.0 / 12.0) * rho * ux_in;
    f[9] = f[8] - 0.5 * (f[4] - f[3]) + (1.0 / 12.0) * rho * ux_in;
    f[11] = f[14] + 0.5 * (f[6] - f[5]) + (1.0 / 12.0) * rho * ux_in;
    f[13] = f[12] - 0.5 * (f[6] - f[5]) + (1.0 / 12.0) * rho * ux_in;
    let _ = sum_known;
    rho
}
/// Zou-He pressure outlet: right wall (x = nx-1) of a D3Q19 domain.
///
/// Sets the outlet density `rho_out` (pressure) and reconstructs the
/// five unknown distributions pointing in the -x direction: `f[2]`,
/// `f[8]`, `f[10]`, `f[12]`, `f[14]`.
///
/// # Returns
/// The computed outlet velocity `ux`.
pub fn zou_he_outlet_right_d3q19(f: &mut [f64; 19], rho_out: f64) -> f64 {
    let sum_zero_x = f[0] + f[3] + f[4] + f[5] + f[6] + f[15] + f[16] + f[17] + f[18];
    let sum_pos_x = f[1] + f[7] + f[9] + f[11] + f[13];
    let ux = -1.0 + (sum_zero_x + 2.0 * sum_pos_x) / rho_out;
    f[2] = f[1] - (2.0 / 3.0) * rho_out * ux;
    f[8] = f[7] + 0.5 * (f[4] - f[3]) - (1.0 / 12.0) * rho_out * ux;
    f[10] = f[9] - 0.5 * (f[4] - f[3]) - (1.0 / 12.0) * rho_out * ux;
    f[12] = f[11] + 0.5 * (f[6] - f[5]) - (1.0 / 12.0) * rho_out * ux;
    f[14] = f[13] - 0.5 * (f[6] - f[5]) - (1.0 / 12.0) * rho_out * ux;
    ux
}
/// Zou-He pressure inlet: left wall (x = 0), D3Q19.
///
/// Sets the inlet density `rho_in` and derives the unknowns.
///
/// # Returns
/// Computed inlet velocity `ux`.
pub fn apply_pressure_inlet_d3q19(f: &mut [f64; 19], rho_in: f64) -> f64 {
    let sum_zero_x = f[0] + f[3] + f[4] + f[5] + f[6] + f[15] + f[16] + f[17] + f[18];
    let sum_neg_x = f[2] + f[8] + f[10] + f[12] + f[14];
    let ux = 1.0 - (sum_zero_x + 2.0 * sum_neg_x) / rho_in;
    f[1] = f[2] + (2.0 / 3.0) * rho_in * ux;
    f[7] = f[10] + 0.5 * (f[4] - f[3]) + (1.0 / 12.0) * rho_in * ux;
    f[9] = f[8] - 0.5 * (f[4] - f[3]) + (1.0 / 12.0) * rho_in * ux;
    f[11] = f[14] + 0.5 * (f[6] - f[5]) + (1.0 / 12.0) * rho_in * ux;
    f[13] = f[12] - 0.5 * (f[6] - f[5]) + (1.0 / 12.0) * rho_in * ux;
    ux
}
/// Zou-He velocity outlet: right wall (x = nx-1), D3Q19.
///
/// Prescribes `ux_out` and reconstructs the five unknown distributions.
///
/// # Returns
/// Computed outlet density `rho`.
pub fn apply_velocity_outlet_d3q19(f: &mut [f64; 19], ux_out: f64) -> f64 {
    let sum_zero_x = f[0] + f[3] + f[4] + f[5] + f[6] + f[15] + f[16] + f[17] + f[18];
    let sum_pos_x = f[1] + f[7] + f[9] + f[11] + f[13];
    let rho = (sum_zero_x + 2.0 * sum_pos_x) / (1.0 + ux_out);
    f[2] = f[1] - (2.0 / 3.0) * rho * ux_out;
    f[8] = f[7] + 0.5 * (f[4] - f[3]) - (1.0 / 12.0) * rho * ux_out;
    f[10] = f[9] - 0.5 * (f[4] - f[3]) - (1.0 / 12.0) * rho * ux_out;
    f[12] = f[11] + 0.5 * (f[6] - f[5]) - (1.0 / 12.0) * rho * ux_out;
    f[14] = f[13] - 0.5 * (f[6] - f[5]) - (1.0 / 12.0) * rho * ux_out;
    rho
}
/// Zou-He velocity boundary condition: top wall (y = ny-1), D2Q9.
///
/// The top face is normal to -y. Unknown distributions pointing into the
/// domain in the -y direction: `f[4]`, `f[7]`, `f[8]`.
///
/// # Arguments
/// * `f`       – mutable reference to the 9-element distribution array.
/// * `uy_wall` – prescribed velocity in the -y direction (typically negative).
///
/// # Returns
/// The computed density `ρ` at this boundary cell.
pub fn zou_he_top_wall_d2q9(f: &mut [f64; 9], ux_wall: f64, uy_wall: f64) -> f64 {
    let rho = (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6])) / (1.0 + uy_wall);
    f[4] = f[2] - (2.0 / 3.0) * rho * uy_wall;
    f[7] = f[5] + 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho * uy_wall - 0.5 * rho * ux_wall;
    f[8] = f[6] - 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho * uy_wall + 0.5 * rho * ux_wall;
    rho
}
/// Zou-He velocity boundary condition: bottom wall (y = 0), D2Q9.
///
/// Unknown distributions pointing into the domain in the +y direction:
/// `f[2]`, `f[5]`, `f[6]`.
///
/// # Returns
/// The computed density `ρ` at this boundary cell.
pub fn zou_he_bottom_wall_d2q9(f: &mut [f64; 9], ux_wall: f64, uy_wall: f64) -> f64 {
    let rho = (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8])) / (1.0 - uy_wall);
    f[2] = f[4] + (2.0 / 3.0) * rho * uy_wall;
    f[5] = f[7] - 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho * uy_wall + 0.5 * rho * ux_wall;
    f[6] = f[8] + 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho * uy_wall - 0.5 * rho * ux_wall;
    rho
}
/// Zou-He pressure BC: top wall (y = ny-1), D2Q9.
///
/// Sets the density `rho_wall` at the top boundary. Returns computed `uy`.
pub fn zou_he_pressure_top_d2q9(f: &mut [f64; 9], rho_wall: f64) -> f64 {
    let uy = -1.0 + (f[0] + f[1] + f[3] + 2.0 * (f[2] + f[5] + f[6])) / rho_wall;
    f[4] = f[2] - (2.0 / 3.0) * rho_wall * uy;
    f[7] = f[5] + 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho_wall * uy;
    f[8] = f[6] - 0.5 * (f[1] - f[3]) - (1.0 / 6.0) * rho_wall * uy;
    uy
}
/// Zou-He pressure BC: bottom wall (y = 0), D2Q9.
///
/// Sets the density `rho_wall` at the bottom boundary. Returns computed `uy`.
pub fn zou_he_pressure_bottom_d2q9(f: &mut [f64; 9], rho_wall: f64) -> f64 {
    let uy = 1.0 - (f[0] + f[1] + f[3] + 2.0 * (f[4] + f[7] + f[8])) / rho_wall;
    f[2] = f[4] + (2.0 / 3.0) * rho_wall * uy;
    f[5] = f[7] - 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho_wall * uy;
    f[6] = f[8] + 0.5 * (f[1] - f[3]) + (1.0 / 6.0) * rho_wall * uy;
    uy
}
/// Zou-He moving wall BC: top wall with prescribed horizontal velocity, D2Q9.
///
/// For lid-driven cavity: the top wall moves at `ux_wall`, `uy_wall = 0`.
/// Returns computed density.
pub fn zou_he_moving_wall_top_d2q9(f: &mut [f64; 9], ux_wall: f64) -> f64 {
    zou_he_top_wall_d2q9(f, ux_wall, 0.0)
}
/// Zou-He moving wall BC: bottom wall with prescribed horizontal velocity, D2Q9.
///
/// Returns computed density.
pub fn zou_he_moving_wall_bottom_d2q9(f: &mut [f64; 9], ux_wall: f64) -> f64 {
    zou_he_bottom_wall_d2q9(f, ux_wall, 0.0)
}
/// Corner treatment for bottom-left corner (x=0, y=0) of D2Q9.
///
/// At a corner, two boundary faces meet. The unknown distributions are
/// determined by bounce-back of the non-equilibrium part plus imposed density.
pub fn zou_he_corner_bottom_left_d2q9(f: &mut [f64; 9], rho_wall: f64) {
    f[1] = f[3];
    f[2] = f[4];
    f[5] = f[7];
    let rho_sum: f64 = f.iter().sum();
    let correction = (rho_wall - rho_sum) / 3.0;
    f[1] += correction;
    f[2] += correction;
    f[5] += correction;
}
/// Corner treatment for bottom-right corner (x=nx-1, y=0) of D2Q9.
pub fn zou_he_corner_bottom_right_d2q9(f: &mut [f64; 9], rho_wall: f64) {
    f[3] = f[1];
    f[2] = f[4];
    f[6] = f[8];
    let rho_sum: f64 = f.iter().sum();
    let correction = (rho_wall - rho_sum) / 3.0;
    f[3] += correction;
    f[2] += correction;
    f[6] += correction;
}
/// Corner treatment for top-left corner (x=0, y=ny-1) of D2Q9.
pub fn zou_he_corner_top_left_d2q9(f: &mut [f64; 9], rho_wall: f64) {
    f[1] = f[3];
    f[4] = f[2];
    f[8] = f[6];
    let rho_sum: f64 = f.iter().sum();
    let correction = (rho_wall - rho_sum) / 3.0;
    f[1] += correction;
    f[4] += correction;
    f[8] += correction;
}
/// Corner treatment for top-right corner (x=nx-1, y=ny-1) of D2Q9.
pub fn zou_he_corner_top_right_d2q9(f: &mut [f64; 9], rho_wall: f64) {
    f[3] = f[1];
    f[4] = f[2];
    f[7] = f[5];
    let rho_sum: f64 = f.iter().sum();
    let correction = (rho_wall - rho_sum) / 3.0;
    f[3] += correction;
    f[4] += correction;
    f[7] += correction;
}
/// Zou-He velocity BC: top face (y = ny-1) of a D3Q19 domain.
///
/// Unknown distributions pointing into the domain (-y direction):
/// `f[4]`, `f[9]`, `f[10]`, `f[16]`, `f[18]`.
///
/// # Returns
/// The computed density `ρ` at this boundary cell.
pub fn zou_he_top_face_d3q19(f: &mut [f64; 19], uy_wall: f64) -> f64 {
    let sum_zero_y = f[0] + f[1] + f[2] + f[5] + f[6] + f[11] + f[12] + f[13] + f[14];
    let sum_pos_y = f[3] + f[7] + f[8] + f[15] + f[17];
    let rho = (sum_zero_y + 2.0 * sum_pos_y) / (1.0 + uy_wall);
    f[4] = f[3] - (2.0 / 3.0) * rho * uy_wall;
    f[9] = f[7] + 0.5 * (f[1] - f[2]) - (1.0 / 12.0) * rho * uy_wall;
    f[10] = f[8] - 0.5 * (f[1] - f[2]) - (1.0 / 12.0) * rho * uy_wall;
    f[16] = f[15] + 0.5 * (f[5] - f[6]) - (1.0 / 12.0) * rho * uy_wall;
    f[18] = f[17] - 0.5 * (f[5] - f[6]) - (1.0 / 12.0) * rho * uy_wall;
    rho
}
/// Zou-He velocity BC: bottom face (y = 0) of a D3Q19 domain.
///
/// Unknown distributions pointing into the domain (+y direction):
/// `f[3]`, `f[7]`, `f[8]`, `f[15]`, `f[17]`.
///
/// # Returns
/// The computed density `ρ` at this boundary cell.
pub fn zou_he_bottom_face_d3q19(f: &mut [f64; 19], uy_wall: f64) -> f64 {
    let sum_zero_y = f[0] + f[1] + f[2] + f[5] + f[6] + f[11] + f[12] + f[13] + f[14];
    let sum_neg_y = f[4] + f[9] + f[10] + f[16] + f[18];
    let rho = (sum_zero_y + 2.0 * sum_neg_y) / (1.0 - uy_wall);
    f[3] = f[4] + (2.0 / 3.0) * rho * uy_wall;
    f[7] = f[10] + 0.5 * (f[2] - f[1]) + (1.0 / 12.0) * rho * uy_wall;
    f[8] = f[9] - 0.5 * (f[2] - f[1]) + (1.0 / 12.0) * rho * uy_wall;
    f[15] = f[18] + 0.5 * (f[6] - f[5]) + (1.0 / 12.0) * rho * uy_wall;
    f[17] = f[16] - 0.5 * (f[6] - f[5]) + (1.0 / 12.0) * rho * uy_wall;
    rho
}
/// Zou-He velocity BC: front face (z = nz-1) of a D3Q19 domain.
///
/// Unknown distributions pointing into the domain (-z direction):
/// `f[6]`, `f[13]`, `f[14]`, `f[17]`, `f[18]`.
///
/// # Returns
/// The computed density `ρ` at this boundary cell.
pub fn zou_he_front_face_d3q19(f: &mut [f64; 19], uz_wall: f64) -> f64 {
    let sum_zero_z = f[0] + f[1] + f[2] + f[3] + f[4] + f[7] + f[8] + f[9] + f[10];
    let sum_pos_z = f[5] + f[11] + f[12] + f[15] + f[16];
    let rho = (sum_zero_z + 2.0 * sum_pos_z) / (1.0 + uz_wall);
    f[6] = f[5] - (2.0 / 3.0) * rho * uz_wall;
    f[13] = f[11] + 0.5 * (f[1] - f[2]) - (1.0 / 12.0) * rho * uz_wall;
    f[14] = f[12] - 0.5 * (f[1] - f[2]) - (1.0 / 12.0) * rho * uz_wall;
    f[17] = f[15] + 0.5 * (f[3] - f[4]) - (1.0 / 12.0) * rho * uz_wall;
    f[18] = f[16] - 0.5 * (f[3] - f[4]) - (1.0 / 12.0) * rho * uz_wall;
    rho
}
/// Zou-He velocity BC: back face (z = 0) of a D3Q19 domain.
///
/// Unknown distributions pointing into the domain (+z direction):
/// `f[5]`, `f[11]`, `f[12]`, `f[15]`, `f[16]`.
///
/// # Returns
/// The computed density `ρ` at this boundary cell.
pub fn zou_he_back_face_d3q19(f: &mut [f64; 19], uz_wall: f64) -> f64 {
    let sum_zero_z = f[0] + f[1] + f[2] + f[3] + f[4] + f[7] + f[8] + f[9] + f[10];
    let sum_neg_z = f[6] + f[13] + f[14] + f[17] + f[18];
    let rho = (sum_zero_z + 2.0 * sum_neg_z) / (1.0 - uz_wall);
    f[5] = f[6] + (2.0 / 3.0) * rho * uz_wall;
    f[11] = f[14] + 0.5 * (f[2] - f[1]) + (1.0 / 12.0) * rho * uz_wall;
    f[12] = f[13] - 0.5 * (f[2] - f[1]) + (1.0 / 12.0) * rho * uz_wall;
    f[15] = f[18] + 0.5 * (f[4] - f[3]) + (1.0 / 12.0) * rho * uz_wall;
    f[16] = f[17] - 0.5 * (f[4] - f[3]) + (1.0 / 12.0) * rho * uz_wall;
    rho
}
/// Zou-He pressure BC: top face (y = ny-1), D3Q19. Returns computed `uy`.
pub fn zou_he_pressure_top_d3q19(f: &mut [f64; 19], rho_wall: f64) -> f64 {
    let sum_zero_y = f[0] + f[1] + f[2] + f[5] + f[6] + f[11] + f[12] + f[13] + f[14];
    let sum_pos_y = f[3] + f[7] + f[8] + f[15] + f[17];
    let uy = -1.0 + (sum_zero_y + 2.0 * sum_pos_y) / rho_wall;
    f[4] = f[3] - (2.0 / 3.0) * rho_wall * uy;
    f[9] = f[7] + 0.5 * (f[1] - f[2]) - (1.0 / 12.0) * rho_wall * uy;
    f[10] = f[8] - 0.5 * (f[1] - f[2]) - (1.0 / 12.0) * rho_wall * uy;
    f[16] = f[15] + 0.5 * (f[5] - f[6]) - (1.0 / 12.0) * rho_wall * uy;
    f[18] = f[17] - 0.5 * (f[5] - f[6]) - (1.0 / 12.0) * rho_wall * uy;
    uy
}
/// Zou-He pressure BC: bottom face (y = 0), D3Q19. Returns computed `uy`.
pub fn zou_he_pressure_bottom_d3q19(f: &mut [f64; 19], rho_wall: f64) -> f64 {
    let sum_zero_y = f[0] + f[1] + f[2] + f[5] + f[6] + f[11] + f[12] + f[13] + f[14];
    let sum_neg_y = f[4] + f[9] + f[10] + f[16] + f[18];
    let uy = 1.0 - (sum_zero_y + 2.0 * sum_neg_y) / rho_wall;
    f[3] = f[4] + (2.0 / 3.0) * rho_wall * uy;
    f[7] = f[10] + 0.5 * (f[2] - f[1]) + (1.0 / 12.0) * rho_wall * uy;
    f[8] = f[9] - 0.5 * (f[2] - f[1]) + (1.0 / 12.0) * rho_wall * uy;
    f[15] = f[18] + 0.5 * (f[6] - f[5]) + (1.0 / 12.0) * rho_wall * uy;
    f[17] = f[16] - 0.5 * (f[6] - f[5]) + (1.0 / 12.0) * rho_wall * uy;
    uy
}
/// Zou-He general inlet: left wall (x=0) with prescribed `ux_in` and `uy_in`.
///
/// Unlike `zou_he_inlet_left`, this variant allows a non-zero transverse
/// velocity component `uy_in`, useful for oblique inflow.
///
/// # Returns
/// Computed density `ρ`.
pub fn zou_he_inlet_left_general_d2q9(f: &mut [f64; 9], ux_in: f64, uy_in: f64) -> f64 {
    let rho = (f[0] + f[2] + f[4] + 2.0 * (f[3] + f[6] + f[7])) / (1.0 - ux_in);
    f[1] = f[3] + (2.0 / 3.0) * rho * ux_in;
    f[5] = f[7] - 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * rho * ux_in + 0.5 * rho * uy_in;
    f[8] = f[6] + 0.5 * (f[2] - f[4]) + (1.0 / 6.0) * rho * ux_in - 0.5 * rho * uy_in;
    rho
}
/// Zou-He general outlet: right wall (x=nx-1) with prescribed density and transverse velocity.
///
/// # Returns
/// Computed velocity `ux`.
pub fn zou_he_outlet_right_general_d2q9(f: &mut [f64; 9], rho_out: f64, uy_out: f64) -> f64 {
    let ux = -1.0 + (f[0] + f[2] + f[4] + 2.0 * (f[1] + f[5] + f[8])) / rho_out;
    f[3] = f[1] - (2.0 / 3.0) * rho_out * ux;
    f[7] = f[5] + 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * rho_out * ux - 0.5 * rho_out * uy_out;
    f[6] = f[8] - 0.5 * (f[2] - f[4]) - (1.0 / 6.0) * rho_out * ux + 0.5 * rho_out * uy_out;
    ux
}
/// Apply a full bounce-back no-slip condition on a left wall cell, D2Q9.
///
/// Reverses all distributions that would leave through x=0. This is a
/// simple half-way bounce-back (first-order accurate).
pub fn bounce_back_left_d2q9(f: &mut [f64; 9]) {
    let (f1, f3) = (f[1], f[3]);
    let (f5, f7) = (f[5], f[7]);
    let (f8, f6) = (f[8], f[6]);
    f[1] = f3;
    f[3] = f1;
    f[5] = f7;
    f[7] = f5;
    f[8] = f6;
    f[6] = f8;
}
/// Apply a full bounce-back no-slip condition on a right wall cell, D2Q9.
pub fn bounce_back_right_d2q9(f: &mut [f64; 9]) {
    let (f1, f3) = (f[1], f[3]);
    let (f5, f7) = (f[5], f[7]);
    let (f8, f6) = (f[8], f[6]);
    f[1] = f3;
    f[3] = f1;
    f[5] = f7;
    f[7] = f5;
    f[8] = f6;
    f[6] = f8;
}
/// Apply a full bounce-back no-slip condition on a bottom wall cell, D2Q9.
pub fn bounce_back_bottom_d2q9(f: &mut [f64; 9]) {
    let (f2, f4) = (f[2], f[4]);
    let (f5, f7) = (f[5], f[7]);
    let (f6, f8) = (f[6], f[8]);
    f[2] = f4;
    f[4] = f2;
    f[5] = f7;
    f[7] = f5;
    f[6] = f8;
    f[8] = f6;
}
/// Apply a full bounce-back no-slip condition on a top wall cell, D2Q9.
pub fn bounce_back_top_d2q9(f: &mut [f64; 9]) {
    let (f2, f4) = (f[2], f[4]);
    let (f5, f7) = (f[5], f[7]);
    let (f6, f8) = (f[6], f[8]);
    f[2] = f4;
    f[4] = f2;
    f[5] = f7;
    f[7] = f5;
    f[6] = f8;
    f[8] = f6;
}
/// Compute macroscopic density and velocity from a D2Q9 distribution.
///
/// Returns `(rho, ux, uy)`.
pub fn macroscopic_d2q9(f: &[f64; 9]) -> (f64, f64, f64) {
    use crate::lattice::D2Q9_VELOCITIES;
    let rho: f64 = f.iter().sum();
    let mut mx = 0.0;
    let mut my = 0.0;
    for (i, &fi) in f.iter().enumerate() {
        mx += fi * D2Q9_VELOCITIES[i][0] as f64;
        my += fi * D2Q9_VELOCITIES[i][1] as f64;
    }
    let ux = if rho > 0.0 { mx / rho } else { 0.0 };
    let uy = if rho > 0.0 { my / rho } else { 0.0 };
    (rho, ux, uy)
}
/// Compute macroscopic density and velocity from a D3Q19 distribution.
///
/// Returns `(rho, ux, uy, uz)`.
pub fn macroscopic_d3q19(f: &[f64; 19]) -> (f64, f64, f64, f64) {
    use crate::lattice::D3Q19_VELOCITIES;
    let rho: f64 = f.iter().sum();
    let mut mx = 0.0;
    let mut my = 0.0;
    let mut mz = 0.0;
    for (i, &fi) in f.iter().enumerate() {
        mx += fi * D3Q19_VELOCITIES[i][0] as f64;
        my += fi * D3Q19_VELOCITIES[i][1] as f64;
        mz += fi * D3Q19_VELOCITIES[i][2] as f64;
    }
    let ux = if rho > 0.0 { mx / rho } else { 0.0 };
    let uy = if rho > 0.0 { my / rho } else { 0.0 };
    let uz = if rho > 0.0 { mz / rho } else { 0.0 };
    (rho, ux, uy, uz)
}
/// Zou-He shear inlet: left wall with a linear shear profile (ux varies with y).
///
/// The inlet velocity is `ux_in + shear_rate * y_frac`, where `y_frac ∈ [0,1]`
/// describes the normalized position along the face height.
///
/// # Arguments
/// * `f`          – mutable reference to the 9-element distribution array.
/// * `ux_base`    – base inlet velocity (at y_frac = 0).
/// * `shear_rate` – velocity increment per unit normalized height.
/// * `y_frac`     – normalized y-position in \[0, 1\].
///
/// # Returns
/// Computed density `ρ`.
pub fn zou_he_shear_inlet_left_d2q9(
    f: &mut [f64; 9],
    ux_base: f64,
    shear_rate: f64,
    y_frac: f64,
) -> f64 {
    let ux_in = ux_base + shear_rate * y_frac;
    zou_he_inlet_left(f, ux_in)
}
/// Zou-He pressure BC: front face (z = nz-1), D3Q19. Returns computed `uz`.
pub fn zou_he_pressure_front_d3q19(f: &mut [f64; 19], rho_wall: f64) -> f64 {
    let sum_zero_z = f[0] + f[1] + f[2] + f[3] + f[4] + f[7] + f[8] + f[9] + f[10];
    let sum_pos_z = f[5] + f[11] + f[12] + f[15] + f[16];
    let uz = -1.0 + (sum_zero_z + 2.0 * sum_pos_z) / rho_wall;
    f[6] = f[5] - (2.0 / 3.0) * rho_wall * uz;
    f[13] = f[11] + 0.5 * (f[1] - f[2]) - (1.0 / 12.0) * rho_wall * uz;
    f[14] = f[12] - 0.5 * (f[1] - f[2]) - (1.0 / 12.0) * rho_wall * uz;
    f[17] = f[15] + 0.5 * (f[3] - f[4]) - (1.0 / 12.0) * rho_wall * uz;
    f[18] = f[16] - 0.5 * (f[3] - f[4]) - (1.0 / 12.0) * rho_wall * uz;
    uz
}
/// Zou-He pressure BC: back face (z = 0), D3Q19. Returns computed `uz`.
pub fn zou_he_pressure_back_d3q19(f: &mut [f64; 19], rho_wall: f64) -> f64 {
    let sum_zero_z = f[0] + f[1] + f[2] + f[3] + f[4] + f[7] + f[8] + f[9] + f[10];
    let sum_neg_z = f[6] + f[13] + f[14] + f[17] + f[18];
    let uz = 1.0 - (sum_zero_z + 2.0 * sum_neg_z) / rho_wall;
    f[5] = f[6] + (2.0 / 3.0) * rho_wall * uz;
    f[11] = f[14] + 0.5 * (f[2] - f[1]) + (1.0 / 12.0) * rho_wall * uz;
    f[12] = f[13] - 0.5 * (f[2] - f[1]) + (1.0 / 12.0) * rho_wall * uz;
    f[15] = f[18] + 0.5 * (f[4] - f[3]) + (1.0 / 12.0) * rho_wall * uz;
    f[16] = f[17] - 0.5 * (f[4] - f[3]) + (1.0 / 12.0) * rho_wall * uz;
    uz
}
/// Zou-He velocity BC: front face (z = nz-1) with prescribed `uz_wall`.
///
/// # Returns
/// Computed density `ρ`.
pub fn zou_he_velocity_front_d3q19(f: &mut [f64; 19], uz_wall: f64) -> f64 {
    zou_he_front_face_d3q19(f, uz_wall)
}
/// Zou-He velocity BC: back face (z = 0) with prescribed `uz_wall`.
///
/// # Returns
/// Computed density `ρ`.
pub fn zou_he_velocity_back_d3q19(f: &mut [f64; 19], uz_wall: f64) -> f64 {
    zou_he_back_face_d3q19(f, uz_wall)
}
/// Zou-He pressure outlet: left face (x = 0), D3Q19.
///
/// Sets the outlet density `rho_out` and reconstructs unknown distributions
/// pointing in the +x direction. Returns the computed velocity `ux`.
pub fn zou_he_pressure_outlet_left_d3q19(f: &mut [f64; 19], rho_out: f64) -> f64 {
    let sum_zero_x = f[0] + f[3] + f[4] + f[5] + f[6] + f[15] + f[16] + f[17] + f[18];
    let sum_pos_x = f[1] + f[7] + f[9] + f[11] + f[13];
    let ux = -1.0 + (sum_zero_x + 2.0 * sum_pos_x) / rho_out;
    f[2] = f[1] - (2.0 / 3.0) * rho_out * ux;
    f[8] = f[7] + 0.5 * (f[4] - f[3]) - (1.0 / 12.0) * rho_out * ux;
    f[10] = f[9] - 0.5 * (f[4] - f[3]) - (1.0 / 12.0) * rho_out * ux;
    f[12] = f[11] + 0.5 * (f[6] - f[5]) - (1.0 / 12.0) * rho_out * ux;
    f[14] = f[13] - 0.5 * (f[6] - f[5]) - (1.0 / 12.0) * rho_out * ux;
    ux
}
/// Compute the D2Q9 Maxwell-Boltzmann equilibrium distribution.
///
/// Uses the standard second-order expansion:
/// f_eq_i = w_i * rho * (1 + (c·u)/cs² + (c·u)²/(2cs⁴) - u²/(2cs²))
/// where cs² = 1/3.
///
/// # Returns
/// A 9-element equilibrium distribution array.
pub fn d2q9_equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS};
    let cs2 = 1.0 / 3.0;
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0_f64; 9];
    for i in 0..9 {
        let cx = D2Q9_VELOCITIES[i][0] as f64;
        let cy = D2Q9_VELOCITIES[i][1] as f64;
        let cu = cx * ux + cy * uy;
        feq[i] = D2Q9_WEIGHTS[i]
            * rho
            * (1.0 + cu / cs2 + cu * cu / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
    }
    feq
}
/// Compute the D3Q19 Maxwell-Boltzmann equilibrium distribution.
///
/// # Returns
/// A 19-element equilibrium distribution array.
pub fn d3q19_equilibrium(rho: f64, ux: f64, uy: f64, uz: f64) -> [f64; 19] {
    use crate::lattice::{D3Q19_VELOCITIES, D3Q19_WEIGHTS};
    let cs2 = 1.0 / 3.0;
    let u2 = ux * ux + uy * uy + uz * uz;
    let mut feq = [0.0_f64; 19];
    for i in 0..19 {
        let cx = D3Q19_VELOCITIES[i][0] as f64;
        let cy = D3Q19_VELOCITIES[i][1] as f64;
        let cz = D3Q19_VELOCITIES[i][2] as f64;
        let cu = cx * ux + cy * uy + cz * uz;
        feq[i] = D3Q19_WEIGHTS[i]
            * rho
            * (1.0 + cu / cs2 + cu * cu / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
    }
    feq
}
/// Zou-He pressure boundary condition for the left face (x = 0) of D3Q19.
///
/// Prescribes density (pressure) at the left inlet; velocity is computed from
/// the non-equilibrium bounce-back method (Zou & He 1997).
///
/// # Arguments
/// * `f`       – mutable reference to the 19-element D3Q19 distribution
/// * `rho_in`  – prescribed inlet density
///
/// # Returns
/// The computed inlet velocity `ux`.
pub fn zou_he_pressure_left_d3q19(f: &mut [f64; 19], rho_in: f64) -> f64 {
    let sum_known = f[0]
        + f[2]
        + f[3]
        + f[4]
        + f[5]
        + f[6]
        + f[8]
        + f[10]
        + f[12]
        + f[14]
        + f[15]
        + f[16]
        + f[17]
        + f[18];
    let ux = 1.0 - (sum_known + 2.0 * (f[2] + f[8] + f[10] + f[12] + f[14])) / rho_in;
    let rho_ux = rho_in * ux;
    f[1] = f[2] + (2.0 / 3.0) * rho_ux;
    f[7] = f[8] - 0.5 * (f[3] - f[4]) + (1.0 / 6.0) * rho_ux;
    f[9] = f[10] + 0.5 * (f[3] - f[4]) + (1.0 / 6.0) * rho_ux;
    f[11] = f[12] - 0.5 * (f[5] - f[6]) + (1.0 / 6.0) * rho_ux;
    f[13] = f[14] + 0.5 * (f[5] - f[6]) + (1.0 / 6.0) * rho_ux;
    ux
}
/// Zou-He pressure boundary condition for the right face (x = Nx-1) of D3Q19.
///
/// Prescribes outlet density; velocity is derived from the non-equilibrium method.
///
/// # Returns
/// The computed outlet velocity `ux`.
pub fn zou_he_pressure_right_d3q19(f: &mut [f64; 19], rho_out: f64) -> f64 {
    let sum_known = f[0]
        + f[1]
        + f[3]
        + f[4]
        + f[5]
        + f[6]
        + f[7]
        + f[9]
        + f[11]
        + f[13]
        + f[15]
        + f[16]
        + f[17]
        + f[18];
    let ux = -1.0 + (sum_known + 2.0 * (f[1] + f[7] + f[9] + f[11] + f[13])) / rho_out;
    let rho_ux = rho_out * ux;
    f[2] = f[1] - (2.0 / 3.0) * rho_ux;
    f[8] = f[7] + 0.5 * (f[3] - f[4]) - (1.0 / 6.0) * rho_ux;
    f[10] = f[9] - 0.5 * (f[3] - f[4]) - (1.0 / 6.0) * rho_ux;
    f[12] = f[11] + 0.5 * (f[5] - f[6]) - (1.0 / 6.0) * rho_ux;
    f[14] = f[13] - 0.5 * (f[5] - f[6]) - (1.0 / 6.0) * rho_ux;
    ux
}
/// D3Q27 weights for reference in boundary conditions.
pub(super) const D3Q27_W: [f64; 27] = [
    8.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    2.0 / 27.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 54.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
    1.0 / 216.0,
];
/// Compute D3Q27 equilibrium distribution.
///
/// Standard second-order quadratic expansion:
/// `feq_i = w_i ρ (1 + c·u/cs² + (c·u)²/(2cs⁴) − u²/(2cs²))`
pub fn d3q27_equilibrium(rho: f64, ux: f64, uy: f64, uz: f64) -> [f64; 27] {
    let cv: [[i32; 3]; 27] = [
        [0, 0, 0],
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
        [1, 1, 0],
        [-1, 1, 0],
        [1, -1, 0],
        [-1, -1, 0],
        [1, 0, 1],
        [-1, 0, 1],
        [1, 0, -1],
        [-1, 0, -1],
        [0, 1, 1],
        [0, -1, 1],
        [0, 1, -1],
        [0, -1, -1],
        [1, 1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [-1, -1, 1],
        [1, 1, -1],
        [-1, 1, -1],
        [1, -1, -1],
        [-1, -1, -1],
    ];
    let u_sq = ux * ux + uy * uy + uz * uz;
    let mut feq = [0.0f64; 27];
    for i in 0..27 {
        let cx = cv[i][0] as f64;
        let cy = cv[i][1] as f64;
        let cz = cv[i][2] as f64;
        let eu = cx * ux + cy * uy + cz * uz;
        feq[i] = D3Q27_W[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u_sq);
    }
    feq
}
/// Zou-He velocity inlet for D3Q27 left face (x = 0).
///
/// Prescribes the full velocity (ux, uy, uz) at the left boundary.
/// The unknown distributions pointing into the domain (+x, +x+y, +x-y,
/// +x+z, +x-z, +x+y+z, +x+y-z, +x-y+z, +x-y-z) are reconstructed
/// using the non-equilibrium bounce-back method.
///
/// # Returns
/// The computed boundary density ρ.
pub fn zou_he_velocity_inlet_d3q27(f: &mut [f64; 27], ux: f64, uy: f64, uz: f64) -> f64 {
    let sum_known = f[0]
        + f[2]
        + f[4]
        + f[6]
        + f[8]
        + f[10]
        + f[12]
        + f[14]
        + f[16]
        + f[18]
        + f[20]
        + f[22]
        + f[24]
        + f[26];
    let sum_left = 2.0 * (f[2] + f[8] + f[10] + f[12] + f[14] + f[20] + f[22] + f[24] + f[26]);
    let rho = (sum_known + sum_left) / (1.0 - ux);
    let feq = d3q27_equilibrium(rho, ux, uy, uz);
    let feq_opp = d3q27_equilibrium(rho, -ux, -uy, -uz);
    let unknowns = [1usize, 7, 9, 11, 13, 19, 21, 23, 25];
    let opposites = [2usize, 8, 10, 12, 14, 20, 22, 24, 26];
    for (&u_idx, &o_idx) in unknowns.iter().zip(opposites.iter()) {
        f[u_idx] = feq[u_idx] - feq_opp[o_idx] + f[o_idx];
    }
    rho
}
/// Zou-He pressure outlet for D3Q27 right face (x = Nx-1).
///
/// Prescribes outlet density ρ_out; velocity is computed from continuity.
///
/// # Returns
/// The computed outlet velocity ux.
pub fn zou_he_pressure_outlet_d3q27(f: &mut [f64; 27], rho_out: f64) -> f64 {
    let sum_known = f[0]
        + f[1]
        + f[3]
        + f[5]
        + f[7]
        + f[9]
        + f[11]
        + f[13]
        + f[15]
        + f[17]
        + f[19]
        + f[21]
        + f[23]
        + f[25];
    let sum_right = 2.0 * (f[1] + f[7] + f[9] + f[11] + f[13] + f[19] + f[21] + f[23] + f[25]);
    let ux = -1.0 + (sum_known + sum_right) / rho_out;
    let feq = d3q27_equilibrium(rho_out, ux, 0.0, 0.0);
    let feq_neg = d3q27_equilibrium(rho_out, -ux, 0.0, 0.0);
    let unknowns = [2usize, 8, 10, 12, 14, 20, 22, 24, 26];
    let opposites = [1usize, 7, 9, 11, 13, 19, 21, 23, 25];
    for (&u_idx, &o_idx) in unknowns.iter().zip(opposites.iter()) {
        f[u_idx] = feq[u_idx] - feq_neg[o_idx] + f[o_idx];
    }
    ux
}
/// Convective outlet boundary condition for D2Q9 right face.
///
/// Implements the advective (convective) outflow condition:
/// `∂f/∂t + uc ∂f/∂x = 0`
/// where uc is the convective velocity.  The boundary cell distribution is
/// updated by: `f_new = (f_old + uc * f_interior) / (1 + uc)`.
///
/// # Arguments
/// * `f_boundary`  – distribution at the outlet cell (x = Nx-1)
/// * `f_interior`  – distribution at the adjacent interior cell (x = Nx-2)
/// * `uc`          – convective velocity (typically local or mean ux)
///
/// # Returns
/// Updated boundary distribution.
pub fn convective_outlet_d2q9(f_boundary: &[f64; 9], f_interior: &[f64; 9], uc: f64) -> [f64; 9] {
    let mut f_new = [0.0f64; 9];
    for i in 0..9 {
        f_new[i] = (f_boundary[i] + uc * f_interior[i]) / (1.0 + uc);
    }
    f_new
}
/// Convective outlet boundary condition for D3Q19 right face.
///
/// Same as the D2Q9 version but for 19-component distributions.
///
/// # Arguments
/// * `f_boundary`  – distribution at the outlet cell
/// * `f_interior`  – distribution at the adjacent interior cell
/// * `uc`          – convective velocity
///
/// # Returns
/// Updated boundary distribution.
pub fn convective_outlet_d3q19(
    f_boundary: &[f64; 19],
    f_interior: &[f64; 19],
    uc: f64,
) -> [f64; 19] {
    let mut f_new = [0.0f64; 19];
    for i in 0..19 {
        f_new[i] = (f_boundary[i] + uc * f_interior[i]) / (1.0 + uc);
    }
    f_new
}
/// Neumann (zero-gradient) outlet boundary condition for D2Q9.
///
/// Implements a zero normal-gradient condition: f_boundary = f_interior.
/// This is the simplest open boundary and works well for fully developed flows.
///
/// # Arguments
/// * `f_boundary` – mutable distribution at outlet cell (updated in place)
/// * `f_interior` – distribution at interior cell
pub fn neumann_outlet_d2q9(f_boundary: &mut [f64; 9], f_interior: &[f64; 9]) {
    f_boundary.copy_from_slice(f_interior);
}
/// Neumann (zero-gradient) outlet boundary condition for D3Q19.
///
/// Copies the interior distribution to the boundary cell for all 19 directions.
pub fn neumann_outlet_d3q19(f_boundary: &mut [f64; 19], f_interior: &[f64; 19]) {
    f_boundary.copy_from_slice(f_interior);
}
/// Second-order extrapolation outlet for D2Q9 (non-reflective).
///
/// Uses a quadratic extrapolation scheme based on two interior nodes:
/// `f_boundary ≈ 2 f_1 − f_2`
/// where f_1 is the first interior cell and f_2 the second.
/// This reduces spurious wave reflections compared to zero-gradient.
///
/// # Arguments
/// * `f_boundary` – mutable distribution at outlet (x = Nx-1)
/// * `f_1`        – distribution at x = Nx-2
/// * `f_2`        – distribution at x = Nx-3
pub fn extrapolation_outlet_d2q9(f_boundary: &mut [f64; 9], f_1: &[f64; 9], f_2: &[f64; 9]) {
    for i in 0..9 {
        f_boundary[i] = 2.0 * f_1[i] - f_2[i];
    }
}
/// Second-order extrapolation outlet for D3Q19.
///
/// Same quadratic extrapolation as the D2Q9 version.
pub fn extrapolation_outlet_d3q19(f_boundary: &mut [f64; 19], f_1: &[f64; 19], f_2: &[f64; 19]) {
    for i in 0..19 {
        f_boundary[i] = 2.0 * f_1[i] - f_2[i];
    }
}
/// Third-order extrapolation outlet for D3Q19.
///
/// Uses three interior nodes: `f_boundary = 3f_1 - 3f_2 + f_3`.
/// Provides better accuracy at the cost of requiring one additional layer.
pub fn extrapolation3_outlet_d3q19(
    f_boundary: &mut [f64; 19],
    f_1: &[f64; 19],
    f_2: &[f64; 19],
    f_3: &[f64; 19],
) {
    for i in 0..19 {
        f_boundary[i] = 3.0 * f_1[i] - 3.0 * f_2[i] + f_3[i];
    }
}
/// Characteristic-based non-reflecting boundary condition for D2Q9 outlet.
///
/// Implements a simplified LODI (Local One-Dimensional Inviscid) approach.
/// The outgoing waves carry information out of the domain; incoming waves
/// are set to zero (fully non-reflecting) or estimated from far-field state.
///
/// Reference: Poinsot & Lele, J. Comput. Phys. 101, 104–129 (1992).
///
/// # Arguments
/// * `f`     – distribution at outlet cell (updated in place)
/// * `rho`   – current local density
/// * `ux`    – current local x-velocity
/// * `uy`    – current local y-velocity
/// * `rho0`  – target far-field density
/// * `sigma` – relaxation coefficient for incoming wave (0 = non-reflecting)
pub fn characteristic_outlet_d2q9(
    f: &mut [f64; 9],
    rho: f64,
    ux: f64,
    uy: f64,
    rho0: f64,
    sigma: f64,
) {
    let cs = (1.0_f64 / 3.0).sqrt();
    let dp = (rho - rho0) * cs * cs;
    let d_rho = -sigma * dp / (cs * cs);
    let rho_new = rho + d_rho;
    let feq = d2q9_equilibrium(rho_new, ux, uy);
    f.copy_from_slice(&feq);
}
/// Characteristic-based non-reflecting boundary condition for D3Q19 outlet.
///
/// Extends the LODI approach to 3D.  The pressure perturbation drives a
/// relaxation of the boundary density toward the target ρ₀.
///
/// # Arguments
/// * `f`     – mutable distribution at outlet
/// * `rho`   – current density
/// * `ux`    – x-velocity
/// * `uy`    – y-velocity
/// * `uz`    – z-velocity
/// * `rho0`  – target far-field density
/// * `sigma` – relaxation coefficient
pub fn characteristic_outlet_d3q19(
    f: &mut [f64; 19],
    rho: f64,
    ux: f64,
    uy: f64,
    uz: f64,
    rho0: f64,
    sigma: f64,
) {
    let cs2 = 1.0 / 3.0;
    let dp = (rho - rho0) * cs2;
    let rho_new = rho - sigma * dp / cs2;
    let feq = d3q19_equilibrium(rho_new, ux, uy, uz);
    f.copy_from_slice(&feq);
}
/// Sponge zone (buffer layer) absorption for D2Q9.
///
/// Gradually damps the distribution toward a target state within a sponge region.
/// The damping strength increases linearly from `sigma_min` at the sponge inlet
/// to `sigma_max` at the domain boundary:
/// `f_new = f + sigma(x) * (feq_target - f)`
///
/// # Arguments
/// * `f`          – mutable distribution at sponge cell
/// * `rho_target` – target density (e.g., far-field value)
/// * `ux_target`  – target x-velocity
/// * `uy_target`  – target y-velocity
/// * `sigma`      – local damping strength ∈ \[0, 1\]
pub fn sponge_zone_d2q9(
    f: &mut [f64; 9],
    rho_target: f64,
    ux_target: f64,
    uy_target: f64,
    sigma: f64,
) {
    let feq_target = d2q9_equilibrium(rho_target, ux_target, uy_target);
    for i in 0..9 {
        f[i] += sigma * (feq_target[i] - f[i]);
    }
}
/// Sponge zone absorption for D3Q19.
///
/// Same damping formula as D2Q9 but extended to 3D.
///
/// # Arguments
/// * `f`          – mutable distribution at sponge cell
/// * `rho_target` – target density
/// * `ux_target`  – target x-velocity
/// * `uy_target`  – target y-velocity
/// * `uz_target`  – target z-velocity
/// * `sigma`      – local damping strength ∈ \[0, 1\]
pub fn sponge_zone_d3q19(
    f: &mut [f64; 19],
    rho_target: f64,
    ux_target: f64,
    uy_target: f64,
    uz_target: f64,
    sigma: f64,
) {
    let feq_target = d3q19_equilibrium(rho_target, ux_target, uy_target, uz_target);
    for i in 0..19 {
        f[i] += sigma * (feq_target[i] - f[i]);
    }
}
/// Compute sponge zone damping coefficient as a function of position.
///
/// Uses a smooth cosine ramp from 0 at `x_start` to `sigma_max` at `x_end`:
/// `σ(x) = σ_max * 0.5 * (1 − cos(π (x − x_start)/(x_end − x_start)))`
///
/// # Arguments
/// * `x`         – current position
/// * `x_start`   – sponge zone start
/// * `x_end`     – sponge zone end (domain boundary)
/// * `sigma_max` – maximum damping coefficient
///
/// # Returns
/// Damping coefficient σ ∈ \[0, sigma_max\].
pub fn sponge_sigma(x: f64, x_start: f64, x_end: f64, sigma_max: f64) -> f64 {
    if x <= x_start {
        return 0.0;
    }
    if x >= x_end {
        return sigma_max;
    }
    let xi = (x - x_start) / (x_end - x_start);
    sigma_max * 0.5 * (1.0 - (std::f64::consts::PI * xi).cos())
}
/// Exponential sponge damping profile.
///
/// `σ(x) = σ_max * (exp(α ξ) − 1) / (exp(α) − 1)`
/// where ξ = (x − x_start) / (x_end − x_start) and α controls steepness.
pub fn sponge_sigma_exponential(
    x: f64,
    x_start: f64,
    x_end: f64,
    sigma_max: f64,
    alpha: f64,
) -> f64 {
    if x <= x_start {
        return 0.0;
    }
    if x >= x_end {
        return sigma_max;
    }
    let xi = (x - x_start) / (x_end - x_start);
    let denom = alpha.exp() - 1.0;
    if denom.abs() < 1e-10 {
        sigma_max * xi
    } else {
        sigma_max * ((alpha * xi).exp() - 1.0) / denom
    }
}
/// Compute macroscopic density and velocity from a D3Q27 distribution.
///
/// # Returns
/// `(rho, ux, uy, uz)` — density and three velocity components.
pub fn macroscopic_d3q27(f: &[f64; 27]) -> (f64, f64, f64, f64) {
    let cv: [[i32; 3]; 27] = [
        [0, 0, 0],
        [1, 0, 0],
        [-1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
        [1, 1, 0],
        [-1, 1, 0],
        [1, -1, 0],
        [-1, -1, 0],
        [1, 0, 1],
        [-1, 0, 1],
        [1, 0, -1],
        [-1, 0, -1],
        [0, 1, 1],
        [0, -1, 1],
        [0, 1, -1],
        [0, -1, -1],
        [1, 1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [-1, -1, 1],
        [1, 1, -1],
        [-1, 1, -1],
        [1, -1, -1],
        [-1, -1, -1],
    ];
    let mut rho = 0.0f64;
    let mut jx = 0.0f64;
    let mut jy = 0.0f64;
    let mut jz = 0.0f64;
    for i in 0..27 {
        rho += f[i];
        jx += f[i] * cv[i][0] as f64;
        jy += f[i] * cv[i][1] as f64;
        jz += f[i] * cv[i][2] as f64;
    }
    let inv_rho = if rho > 1e-14 { 1.0 / rho } else { 1.0 };
    (rho, jx * inv_rho, jy * inv_rho, jz * inv_rho)
}
/// Compute a 2D parabolic (Poiseuille) velocity profile.
///
/// Returns the x-velocity at position y ∈ \[0, H\] given peak velocity u_max:
/// `u(y) = u_max * 4y(H − y)/H²`
///
/// # Arguments
/// * `y`     – transverse coordinate
/// * `h`     – channel height
/// * `u_max` – peak centerline velocity
pub fn poiseuille_profile(y: f64, h: f64, u_max: f64) -> f64 {
    u_max * 4.0 * y * (h - y) / (h * h)
}
/// Compute an error-function (tanh) shear layer profile.
///
/// `u(y) = 0.5 * u_max * (1 + tanh((y − y_c) / δ))`
///
/// # Arguments
/// * `y`     – transverse coordinate
/// * `y_c`   – center of the shear layer
/// * `delta` – shear layer thickness
/// * `u_max` – velocity jump amplitude
pub fn tanh_shear_profile(y: f64, y_c: f64, delta: f64, u_max: f64) -> f64 {
    0.5 * u_max * (1.0 + ((y - y_c) / delta).tanh())
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS, D3Q19_VELOCITIES, D3Q19_WEIGHTS};
    fn d2q9_equilibrium_at_rest() -> [f64; 9] {
        D2Q9_WEIGHTS
    }
    fn d3q19_equilibrium_at_rest() -> [f64; 19] {
        D3Q19_WEIGHTS
    }
    #[test]
    fn test_zou_he_inlet_left_mass_conservation() {
        let mut f = d2q9_equilibrium_at_rest();
        let ux_wall = 0.05;
        let rho_computed = zou_he_inlet_left(&mut f, ux_wall);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_computed - rho_sum).abs() < 1e-12,
            "Inlet rho mismatch: returned={rho_computed}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_inlet_left_velocity() {
        let mut f = d2q9_equilibrium_at_rest();
        let ux_wall = 0.07;
        let rho = zou_he_inlet_left(&mut f, ux_wall);
        let mut mx = 0.0_f64;
        for (i, &fi) in f.iter().enumerate() {
            mx += fi * D2Q9_VELOCITIES[i][0] as f64;
        }
        let ux_actual = mx / rho;
        assert!(
            (ux_actual - ux_wall).abs() < 1e-12,
            "Inlet ux = {ux_actual}, expected {ux_wall}"
        );
    }
    #[test]
    fn test_zou_he_outlet_right_mass_conservation() {
        let mut f = d2q9_equilibrium_at_rest();
        let rho_wall = 1.0;
        zou_he_outlet_right(&mut f, rho_wall);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - rho_wall).abs() < 1e-12,
            "Outlet rho sum = {rho_sum}, expected {rho_wall}"
        );
    }
    #[test]
    fn test_pressure_inlet_d2q9_positive_velocity() {
        let mut f = d2q9_equilibrium_at_rest();
        let ux = apply_pressure_inlet_d2q9(&mut f, 1.05);
        assert!(
            ux > 0.0,
            "Pressure inlet with rho > 1 should give positive ux, got {ux}"
        );
    }
    #[test]
    fn test_velocity_outlet_d2q9_density_consistent() {
        let mut f = d2q9_equilibrium_at_rest();
        f[1] += 0.01;
        f[5] += 0.005;
        f[8] += 0.005;
        let rho_ret = apply_velocity_outlet_d2q9(&mut f, 0.02);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_ret - rho_sum).abs() < 1e-10,
            "Velocity outlet rho mismatch: returned={rho_ret}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_inlet_left_d3q19_mass_conservation() {
        let mut f = d3q19_equilibrium_at_rest();
        let ux_in = 0.03;
        let rho_computed = zou_he_inlet_left_d3q19(&mut f, ux_in);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_computed - rho_sum).abs() < 1e-12,
            "D3Q19 inlet rho mismatch: returned={rho_computed}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_inlet_left_d3q19_velocity() {
        let mut f = d3q19_equilibrium_at_rest();
        let ux_in = 0.04;
        let rho = zou_he_inlet_left_d3q19(&mut f, ux_in);
        let mut mx = 0.0_f64;
        for (i, &fi) in f.iter().enumerate() {
            mx += fi * D3Q19_VELOCITIES[i][0] as f64;
        }
        let ux_actual = mx / rho;
        assert!(
            (ux_actual - ux_in).abs() < 1e-11,
            "D3Q19 inlet ux = {ux_actual}, expected {ux_in}"
        );
    }
    #[test]
    fn test_pressure_inlet_d3q19_positive_velocity() {
        let mut f = d3q19_equilibrium_at_rest();
        let ux = apply_pressure_inlet_d3q19(&mut f, 1.05);
        assert!(
            ux > 0.0,
            "D3Q19 pressure inlet should give positive ux, got {ux}"
        );
    }
    #[test]
    fn test_velocity_outlet_d3q19_density_consistent() {
        let mut f = d3q19_equilibrium_at_rest();
        f[1] += 0.01;
        f[7] += 0.005;
        f[9] += 0.005;
        let rho_ret = apply_velocity_outlet_d3q19(&mut f, 0.02);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_ret - rho_sum).abs() < 1e-10,
            "D3Q19 velocity outlet rho mismatch: returned={rho_ret}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_d2q9_inlet_outlet_net_mass() {
        let mut f_in = d2q9_equilibrium_at_rest();
        let mut f_out = d2q9_equilibrium_at_rest();
        let rho_in = zou_he_inlet_left(&mut f_in, 0.05);
        let rho_out_wall = 1.0;
        zou_he_outlet_right(&mut f_out, rho_out_wall);
        let rho_out_actual: f64 = f_out.iter().sum();
        assert!(
            (rho_out_actual - rho_out_wall).abs() < 1e-12,
            "Outlet rho mismatch after round-trip"
        );
        assert!(
            rho_in > 0.9,
            "Inlet density should remain near 1, got {rho_in}"
        );
    }
    #[test]
    fn test_zou_he_top_wall_mass_conservation() {
        let mut f = d2q9_equilibrium_at_rest();
        let rho = zou_he_top_wall_d2q9(&mut f, 0.0, 0.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho - rho_sum).abs() < 1e-12,
            "Top wall rho mismatch: returned={rho}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_bottom_wall_mass_conservation() {
        let mut f = d2q9_equilibrium_at_rest();
        let rho = zou_he_bottom_wall_d2q9(&mut f, 0.0, 0.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho - rho_sum).abs() < 1e-12,
            "Bottom wall rho mismatch: returned={rho}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_moving_wall_top_velocity() {
        let mut f = d2q9_equilibrium_at_rest();
        let ux_wall = 0.05;
        let rho = zou_he_moving_wall_top_d2q9(&mut f, ux_wall);
        let mut mx = 0.0_f64;
        for (i, &fi) in f.iter().enumerate() {
            mx += fi * D2Q9_VELOCITIES[i][0] as f64;
        }
        let ux_actual = mx / rho;
        assert!(
            (ux_actual - ux_wall).abs() < 1e-10,
            "Moving wall ux = {ux_actual}, expected {ux_wall}"
        );
    }
    #[test]
    fn test_zou_he_pressure_top_mass() {
        let mut f = d2q9_equilibrium_at_rest();
        zou_he_pressure_top_d2q9(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-12,
            "Pressure top rho = {rho_sum}, expected 1.0"
        );
    }
    #[test]
    fn test_zou_he_pressure_bottom_mass() {
        let mut f = d2q9_equilibrium_at_rest();
        zou_he_pressure_bottom_d2q9(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-12,
            "Pressure bottom rho = {rho_sum}, expected 1.0"
        );
    }
    #[test]
    fn test_corner_bottom_left_density() {
        let mut f = d2q9_equilibrium_at_rest();
        zou_he_corner_bottom_left_d2q9(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-10,
            "Corner BL rho = {rho_sum}, expected 1.0"
        );
    }
    #[test]
    fn test_corner_top_right_density() {
        let mut f = d2q9_equilibrium_at_rest();
        zou_he_corner_top_right_d2q9(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-10,
            "Corner TR rho = {rho_sum}, expected 1.0"
        );
    }
    #[test]
    fn test_corner_bottom_right_density() {
        let mut f = d2q9_equilibrium_at_rest();
        zou_he_corner_bottom_right_d2q9(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-10,
            "Corner BR rho = {rho_sum}, expected 1.0"
        );
    }
    #[test]
    fn test_corner_top_left_density() {
        let mut f = d2q9_equilibrium_at_rest();
        zou_he_corner_top_left_d2q9(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-10,
            "Corner TL rho = {rho_sum}, expected 1.0"
        );
    }
    #[test]
    fn test_zou_he_top_face_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        let rho = zou_he_top_face_d3q19(&mut f, 0.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho - rho_sum).abs() < 1e-12,
            "D3Q19 top face rho mismatch: returned={rho}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_bottom_face_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        let rho = zou_he_bottom_face_d3q19(&mut f, 0.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho - rho_sum).abs() < 1e-12,
            "D3Q19 bottom face rho mismatch: returned={rho}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_front_face_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        let rho = zou_he_front_face_d3q19(&mut f, 0.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho - rho_sum).abs() < 1e-12,
            "D3Q19 front face rho mismatch: returned={rho}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_back_face_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        let rho = zou_he_back_face_d3q19(&mut f, 0.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho - rho_sum).abs() < 1e-12,
            "D3Q19 back face rho mismatch: returned={rho}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_pressure_top_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        zou_he_pressure_top_d3q19(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-12,
            "D3Q19 pressure top rho = {rho_sum}, expected 1.0"
        );
    }
    #[test]
    fn test_zou_he_pressure_bottom_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        zou_he_pressure_bottom_d3q19(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-12,
            "D3Q19 pressure bottom rho = {rho_sum}, expected 1.0"
        );
    }
    #[test]
    fn test_zou_he_top_face_d3q19_nonzero_velocity() {
        let mut f = d3q19_equilibrium_at_rest();
        let f_before = f;
        let rho = zou_he_top_face_d3q19(&mut f, -0.02);
        assert!(rho > 0.0);
        let changed = f
            .iter()
            .zip(f_before.iter())
            .any(|(a, b)| (a - b).abs() > 1e-14);
        assert!(
            changed,
            "Non-zero wall velocity should modify distributions"
        );
    }
    #[test]
    fn test_moving_wall_bottom_zero_vel() {
        let mut f = d2q9_equilibrium_at_rest();
        let rho = zou_he_moving_wall_bottom_d2q9(&mut f, 0.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho - rho_sum).abs() < 1e-12,
            "Moving bottom wall zero vel: rho mismatch"
        );
        assert!(
            (rho - 1.0).abs() < 1e-10,
            "Moving bottom wall zero vel: rho should be ~1.0, got {rho}"
        );
    }
    #[test]
    fn test_d2q9_equilibrium_sum_equals_rho() {
        for &rho in &[0.8, 1.0, 1.2] {
            for &ux in &[-0.05, 0.0, 0.05] {
                for &uy in &[-0.03, 0.0, 0.03] {
                    let feq = d2q9_equilibrium(rho, ux, uy);
                    let sum: f64 = feq.iter().sum();
                    assert!(
                        (sum - rho).abs() < 1e-12,
                        "d2q9_equilibrium sum={sum}, expected rho={rho}"
                    );
                }
            }
        }
    }
    #[test]
    fn test_d3q19_equilibrium_sum_equals_rho() {
        for &rho in &[0.9, 1.0, 1.1] {
            let feq = d3q19_equilibrium(rho, 0.02, -0.01, 0.03);
            let sum: f64 = feq.iter().sum();
            assert!(
                (sum - rho).abs() < 1e-12,
                "d3q19_equilibrium sum={sum}, expected rho={rho}"
            );
        }
    }
    #[test]
    fn test_macroscopic_d2q9_from_equilibrium() {
        let rho_ref = 1.05;
        let ux_ref = 0.04;
        let uy_ref = -0.02;
        let feq = d2q9_equilibrium(rho_ref, ux_ref, uy_ref);
        let (rho, ux, uy) = macroscopic_d2q9(&feq);
        assert!(
            (rho - rho_ref).abs() < 1e-12,
            "rho mismatch: {rho} vs {rho_ref}"
        );
        assert!((ux - ux_ref).abs() < 1e-12, "ux mismatch: {ux} vs {ux_ref}");
        assert!((uy - uy_ref).abs() < 1e-12, "uy mismatch: {uy} vs {uy_ref}");
    }
    #[test]
    fn test_macroscopic_d3q19_from_equilibrium() {
        let rho_ref = 0.98;
        let ux_ref = 0.03;
        let uy_ref = -0.01;
        let uz_ref = 0.02;
        let feq = d3q19_equilibrium(rho_ref, ux_ref, uy_ref, uz_ref);
        let (rho, ux, uy, uz) = macroscopic_d3q19(&feq);
        assert!((rho - rho_ref).abs() < 1e-12, "rho mismatch: {rho}");
        assert!((ux - ux_ref).abs() < 1e-12, "ux mismatch: {ux}");
        assert!((uy - uy_ref).abs() < 1e-12, "uy mismatch: {uy}");
        assert!((uz - uz_ref).abs() < 1e-12, "uz mismatch: {uz}");
    }
    #[test]
    fn test_zou_he_inlet_left_general_mass() {
        let mut f = d2q9_equilibrium_at_rest();
        let rho = zou_he_inlet_left_general_d2q9(&mut f, 0.05, 0.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho - rho_sum).abs() < 1e-12,
            "General inlet rho mismatch: returned={rho}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_inlet_left_general_ux() {
        let mut f = d2q9_equilibrium_at_rest();
        let ux_in = 0.06;
        let rho = zou_he_inlet_left_general_d2q9(&mut f, ux_in, 0.0);
        let mx: f64 = f
            .iter()
            .enumerate()
            .map(|(i, fi)| fi * D2Q9_VELOCITIES[i][0] as f64)
            .sum();
        let ux = mx / rho;
        assert!(
            (ux - ux_in).abs() < 1e-11,
            "General inlet ux={ux}, expected {ux_in}"
        );
    }
    #[test]
    fn test_zou_he_outlet_right_general_mass() {
        let mut f = d2q9_equilibrium_at_rest();
        f[1] += 0.01;
        let rho_out = 1.0;
        zou_he_outlet_right_general_d2q9(&mut f, rho_out, 0.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - rho_out).abs() < 1e-10,
            "General outlet rho_sum={rho_sum}, expected {rho_out}"
        );
    }
    #[test]
    fn test_bounce_back_left_density_conserved() {
        let mut f = [0.1, 0.2, 0.15, 0.12, 0.08, 0.07, 0.09, 0.1, 0.09];
        let rho_before: f64 = f.iter().sum();
        bounce_back_left_d2q9(&mut f);
        let rho_after: f64 = f.iter().sum();
        assert!(
            (rho_before - rho_after).abs() < 1e-14,
            "BB left density changed"
        );
    }
    #[test]
    fn test_bounce_back_bottom_density_conserved() {
        let mut f = [0.11, 0.13, 0.14, 0.12, 0.10, 0.09, 0.08, 0.12, 0.11];
        let rho_before: f64 = f.iter().sum();
        bounce_back_bottom_d2q9(&mut f);
        let rho_after: f64 = f.iter().sum();
        assert!(
            (rho_before - rho_after).abs() < 1e-14,
            "BB bottom density changed"
        );
    }
    #[test]
    fn test_bounce_back_top_involutory() {
        let f_orig = [0.1, 0.2, 0.15, 0.12, 0.08, 0.07, 0.09, 0.1, 0.09];
        let mut f = f_orig;
        bounce_back_top_d2q9(&mut f);
        bounce_back_top_d2q9(&mut f);
        for i in 0..9 {
            assert!(
                (f[i] - f_orig[i]).abs() < 1e-15,
                "BB top not involutory at {i}"
            );
        }
    }
    #[test]
    fn test_bounce_back_right_involutory() {
        let f_orig = [0.12, 0.18, 0.14, 0.11, 0.09, 0.08, 0.1, 0.11, 0.07];
        let mut f = f_orig;
        bounce_back_right_d2q9(&mut f);
        bounce_back_right_d2q9(&mut f);
        for i in 0..9 {
            assert!(
                (f[i] - f_orig[i]).abs() < 1e-15,
                "BB right not involutory at {i}"
            );
        }
    }
    #[test]
    fn test_shear_inlet_at_zero_frac() {
        let ux_base = 0.03;
        let shear = 0.05;
        let mut f = d2q9_equilibrium_at_rest();
        let rho = zou_he_shear_inlet_left_d2q9(&mut f, ux_base, shear, 0.0);
        let mx: f64 = f
            .iter()
            .enumerate()
            .map(|(i, fi)| fi * D2Q9_VELOCITIES[i][0] as f64)
            .sum();
        let ux = mx / rho;
        assert!(
            (ux - ux_base).abs() < 1e-11,
            "Shear inlet at y_frac=0: ux={ux}, expected {ux_base}"
        );
    }
    #[test]
    fn test_shear_inlet_at_one_frac() {
        let ux_base = 0.02;
        let shear = 0.04;
        let mut f = d2q9_equilibrium_at_rest();
        let rho = zou_he_shear_inlet_left_d2q9(&mut f, ux_base, shear, 1.0);
        let mx: f64 = f
            .iter()
            .enumerate()
            .map(|(i, fi)| fi * D2Q9_VELOCITIES[i][0] as f64)
            .sum();
        let ux = mx / rho;
        let ux_expected = ux_base + shear;
        assert!(
            (ux - ux_expected).abs() < 1e-11,
            "Shear inlet at y_frac=1: ux={ux}, expected {ux_expected}"
        );
    }
    #[test]
    fn test_zou_he_pressure_front_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        zou_he_pressure_front_d3q19(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-12,
            "Pressure front rho={rho_sum}"
        );
    }
    #[test]
    fn test_zou_he_pressure_back_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        zou_he_pressure_back_d3q19(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!((rho_sum - 1.0).abs() < 1e-12, "Pressure back rho={rho_sum}");
    }
    #[test]
    fn test_zou_he_velocity_front_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        let rho = zou_he_velocity_front_d3q19(&mut f, 0.02);
        let rho_sum: f64 = f.iter().sum();
        assert!((rho - rho_sum).abs() < 1e-12, "Velocity front rho mismatch");
    }
    #[test]
    fn test_zou_he_velocity_back_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        let rho = zou_he_velocity_back_d3q19(&mut f, 0.02);
        let rho_sum: f64 = f.iter().sum();
        assert!((rho - rho_sum).abs() < 1e-12, "Velocity back rho mismatch");
    }
    #[test]
    fn test_zou_he_pressure_outlet_left_d3q19_mass() {
        let mut f = d3q19_equilibrium_at_rest();
        zou_he_pressure_outlet_left_d3q19(&mut f, 1.0);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_sum - 1.0).abs() < 1e-12,
            "Pressure outlet left rho={rho_sum}"
        );
    }
    #[test]
    fn test_pressure_inlet_d3q19_ux_positive_for_high_rho() {
        let mut f = d3q19_equilibrium_at_rest();
        let ux = apply_pressure_inlet_d3q19(&mut f, 1.1);
        assert!(ux > 0.0, "Expected positive ux for rho_in=1.1, got {ux}");
    }
    #[test]
    fn test_velocity_outlet_d3q19_rho_consistent_with_prescribed_ux() {
        let mut f = d3q19_equilibrium_at_rest();
        f[1] += 0.02;
        f[7] += 0.01;
        let rho_ret = apply_velocity_outlet_d3q19(&mut f, 0.03);
        let rho_sum: f64 = f.iter().sum();
        assert!(
            (rho_ret - rho_sum).abs() < 1e-10,
            "Velocity outlet D3Q19 rho mismatch: returned={rho_ret}, sum={rho_sum}"
        );
    }
    #[test]
    fn test_pressure_front_back_d3q19_density() {
        let rho_set = 1.02;
        let mut f_front = d3q19_equilibrium_at_rest();
        let mut f_back = d3q19_equilibrium_at_rest();
        zou_he_pressure_front_d3q19(&mut f_front, rho_set);
        zou_he_pressure_back_d3q19(&mut f_back, rho_set);
        let rho_front: f64 = f_front.iter().sum();
        let rho_back: f64 = f_back.iter().sum();
        assert!(
            (rho_front - rho_set).abs() < 1e-12,
            "Pressure front rho={rho_front}"
        );
        assert!(
            (rho_back - rho_set).abs() < 1e-12,
            "Pressure back  rho={rho_back}"
        );
    }
    #[test]
    fn test_macroscopic_d2q9_zero_velocity_at_rest() {
        let f = d2q9_equilibrium_at_rest();
        let (rho, ux, uy) = macroscopic_d2q9(&f);
        assert!(
            (rho - 1.0).abs() < 1e-14,
            "rho at rest should be 1.0: {rho}"
        );
        assert!(ux.abs() < 1e-14, "ux at rest should be 0: {ux}");
        assert!(uy.abs() < 1e-14, "uy at rest should be 0: {uy}");
    }
    #[test]
    fn test_macroscopic_d3q19_zero_velocity_at_rest() {
        let f = d3q19_equilibrium_at_rest();
        let (rho, ux, uy, uz) = macroscopic_d3q19(&f);
        assert!(
            (rho - 1.0).abs() < 1e-14,
            "D3Q19 rho at rest should be 1.0: {rho}"
        );
        assert!(ux.abs() < 1e-14, "D3Q19 ux at rest: {ux}");
        assert!(uy.abs() < 1e-14, "D3Q19 uy at rest: {uy}");
        assert!(uz.abs() < 1e-14, "D3Q19 uz at rest: {uz}");
    }
    #[test]
    fn test_d2q9_equilibrium_zero_vel_equals_weights() {
        let feq = d2q9_equilibrium(1.0, 0.0, 0.0);
        let weights = d2q9_equilibrium_at_rest();
        for i in 0..9 {
            assert!(
                (feq[i] - weights[i]).abs() < 1e-14,
                "feq[{i}]={} != weight[{i}]={}",
                feq[i],
                weights[i]
            );
        }
    }
    #[test]
    fn test_d3q19_equilibrium_zero_vel_equals_weights() {
        let feq = d3q19_equilibrium(1.0, 0.0, 0.0, 0.0);
        let weights = d3q19_equilibrium_at_rest();
        for i in 0..19 {
            assert!(
                (feq[i] - weights[i]).abs() < 1e-14,
                "feq[{i}]={} != weight[{i}]={}",
                feq[i],
                weights[i]
            );
        }
    }
    #[test]
    fn test_zou_he_inlet_left_known_distributions_unchanged() {
        let mut f = d2q9_equilibrium_at_rest();
        let f_orig = f;
        zou_he_inlet_left(&mut f, 0.05);
        for &i in &[0usize, 2, 3, 4, 6, 7] {
            assert!(
                (f[i] - f_orig[i]).abs() < 1e-15,
                "Known distribution f[{i}] was modified by inlet BC"
            );
        }
    }
    #[test]
    fn test_zou_he_outlet_right_known_distributions_unchanged() {
        let mut f = d2q9_equilibrium_at_rest();
        let f_orig = f;
        zou_he_outlet_right(&mut f, 1.0);
        for &i in &[0usize, 1, 2, 4, 5, 8] {
            assert!(
                (f[i] - f_orig[i]).abs() < 1e-15,
                "Known distribution f[{i}] was modified by outlet BC"
            );
        }
    }
    #[test]
    fn test_all_pressure_faces_d3q19_positive_rho() {
        let rho_set = 1.0;
        for apply in [
            zou_he_pressure_top_d3q19 as fn(&mut [f64; 19], f64) -> f64,
            zou_he_pressure_bottom_d3q19,
            zou_he_pressure_front_d3q19,
            zou_he_pressure_back_d3q19,
        ] {
            let mut f = d3q19_equilibrium_at_rest();
            apply(&mut f, rho_set);
            let sum: f64 = f.iter().sum();
            assert!(
                (sum - rho_set).abs() < 1e-12,
                "D3Q19 pressure face BC: rho={sum}, expected {rho_set}"
            );
        }
    }
    #[test]
    fn test_zou_he_inlet_left_d3q19_known_unchanged() {
        let mut f = d3q19_equilibrium_at_rest();
        let f_orig = f;
        zou_he_inlet_left_d3q19(&mut f, 0.03);
        let known = [0, 2, 3, 4, 5, 6, 8, 10, 12, 14, 15, 16, 17, 18];
        for &i in &known {
            assert!(
                (f[i] - f_orig[i]).abs() < 1e-15,
                "Known D3Q19 inlet f[{i}] was modified"
            );
        }
    }
    #[test]
    fn test_general_outlet_nonzero_uy_modifies_unknowns() {
        let mut f1 = d2q9_equilibrium_at_rest();
        let mut f2 = d2q9_equilibrium_at_rest();
        zou_he_outlet_right_general_d2q9(&mut f1, 1.0, 0.0);
        zou_he_outlet_right_general_d2q9(&mut f2, 1.0, 0.02);
        let differs = (f1[6] - f2[6]).abs() > 1e-14 || (f1[7] - f2[7]).abs() > 1e-14;
        assert!(differs, "Non-zero uy_out should change f[6] or f[7]");
    }
    #[test]
    fn test_d2q9_equilibrium_non_negative_small_vel() {
        let feq = d2q9_equilibrium(1.0, 0.05, 0.02);
        for (i, &fi) in feq.iter().enumerate() {
            assert!(fi >= 0.0, "feq[{i}]={fi} is negative for small velocity");
        }
    }
    #[test]
    fn test_d3q19_equilibrium_non_negative_small_vel() {
        let feq = d3q19_equilibrium(1.0, 0.04, 0.01, -0.02);
        for (i, &fi) in feq.iter().enumerate() {
            assert!(fi >= 0.0, "D3Q19 feq[{i}]={fi} is negative");
        }
    }
    #[test]
    fn test_inlet_outlet_round_trip_d3q19() {
        let mut f_in = d3q19_equilibrium_at_rest();
        let rho_in = zou_he_inlet_left_d3q19(&mut f_in, 0.03);
        let mut f_out = d3q19_equilibrium_at_rest();
        zou_he_outlet_right_d3q19(&mut f_out, 1.0);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_in - 1.0).abs() < 0.05,
            "Inlet rho should be near 1: {rho_in}"
        );
        assert!(
            (rho_out - 1.0).abs() < 1e-12,
            "Outlet rho should equal prescribed 1: {rho_out}"
        );
    }
}
