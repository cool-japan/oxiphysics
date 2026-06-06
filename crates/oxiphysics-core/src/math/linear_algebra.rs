//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

pub use nalgebra::{Matrix3, Matrix4, Point3, Unit, UnitQuaternion, Vector3, Vector4};

use super::geometry::{quat_exp, quat_log};
use super::types::*;

/// Scalar type used throughout the engine (f64 by default).
pub type Real = f64;
/// 3D vector type alias.
pub type Vec3 = Vector3<Real>;
/// 4D vector type alias.
pub type Vec4 = Vector4<Real>;
/// 3x3 matrix type alias.
pub type Mat3 = Matrix3<Real>;
/// 4x4 matrix type alias.
pub type Mat4 = Matrix4<Real>;
/// Quaternion type alias.
pub type Quat = UnitQuaternion<Real>;
/// Create a quaternion from an axis-angle representation.
///
/// `axis` must be a unit vector. `angle` is in radians.
pub fn quat_from_axis_angle(axis: &Vec3, angle: Real) -> Quat {
    UnitQuaternion::from_axis_angle(&Unit::new_normalize(*axis), angle)
}
/// Convert a quaternion to Euler angles (roll, pitch, yaw) in radians.
///
/// Uses the ZYX (yaw-pitch-roll) convention.
/// Returns `(roll, pitch, yaw)`.
pub fn quat_to_euler(q: &Quat) -> (Real, Real, Real) {
    let (roll, pitch, yaw) = q.euler_angles();
    (roll, pitch, yaw)
}
/// Create a quaternion from Euler angles (roll, pitch, yaw) in radians.
///
/// Uses the ZYX convention: rotation = Rz(yaw) * Ry(pitch) * Rx(roll).
pub fn quat_from_euler(roll: Real, pitch: Real, yaw: Real) -> Quat {
    UnitQuaternion::from_euler_angles(roll, pitch, yaw)
}
/// Spherical linear interpolation between two quaternions.
pub fn quat_slerp(a: &Quat, b: &Quat, t: Real) -> Quat {
    a.slerp(b, t)
}
/// Compute the determinant of a 3x3 matrix.
pub fn mat3_determinant(m: &Mat3) -> Real {
    m.determinant()
}
/// Compute the inverse of a 3x3 matrix.
///
/// Returns `None` if the matrix is singular.
pub fn mat3_inverse(m: &Mat3) -> Option<Mat3> {
    m.try_inverse()
}
/// Compute real eigenvalues of a symmetric 3x3 matrix.
///
/// Returns eigenvalues sorted in ascending order. For non-symmetric matrices
/// the result is an approximation.
pub fn mat3_eigenvalues_symmetric(m: &Mat3) -> [Real; 3] {
    let eigen = m.symmetric_eigen();
    let mut vals = [
        eigen.eigenvalues[0],
        eigen.eigenvalues[1],
        eigen.eigenvalues[2],
    ];
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    vals
}
/// Compute the trace of a 3x3 matrix.
pub fn mat3_trace(m: &Mat3) -> Real {
    m.trace()
}
/// Compute the Frobenius norm of a 3x3 matrix.
pub fn mat3_frobenius_norm(m: &Mat3) -> Real {
    let sum: f64 = (0..3)
        .flat_map(|i| (0..3).map(move |j| m[(i, j)] * m[(i, j)]))
        .sum();
    sum.sqrt()
}
/// Build a skew-symmetric matrix from a vector (for cross product as matrix multiply).
pub fn skew_symmetric(v: &Vec3) -> Mat3 {
    Mat3::new(0.0, -v.z, v.y, v.z, 0.0, -v.x, -v.y, v.x, 0.0)
}
/// Create a perspective projection matrix (right-handed, depth \[0, 1\]).
///
/// * `fov_y` - Vertical field of view in radians.
/// * `aspect` - Width / height ratio.
/// * `near` - Near clipping plane distance (positive).
/// * `far` - Far clipping plane distance (positive).
pub fn perspective(fov_y: Real, aspect: Real, near: Real, far: Real) -> Mat4 {
    let f = 1.0 / (fov_y / 2.0).tan();
    let nf = 1.0 / (near - far);
    Mat4::new(
        f / aspect,
        0.0,
        0.0,
        0.0,
        0.0,
        f,
        0.0,
        0.0,
        0.0,
        0.0,
        (far + near) * nf,
        2.0 * far * near * nf,
        0.0,
        0.0,
        -1.0,
        0.0,
    )
}
/// Create a look-at view matrix (right-handed).
///
/// * `eye` - Camera position.
/// * `target` - Look-at target position.
/// * `up` - World up direction.
pub fn look_at(eye: &Vec3, target: &Vec3, up: &Vec3) -> Mat4 {
    let f = (target - eye).normalize();
    let s = f.cross(up).normalize();
    let u = s.cross(&f);
    Mat4::new(
        s.x,
        s.y,
        s.z,
        -s.dot(eye),
        u.x,
        u.y,
        u.z,
        -u.dot(eye),
        -f.x,
        -f.y,
        -f.z,
        f.dot(eye),
        0.0,
        0.0,
        0.0,
        1.0,
    )
}
/// Create an orthographic projection matrix (right-handed).
pub fn orthographic(
    left: Real,
    right: Real,
    bottom: Real,
    top: Real,
    near: Real,
    far: Real,
) -> Mat4 {
    let rml = right - left;
    let tmb = top - bottom;
    let fmn = far - near;
    Mat4::new(
        2.0 / rml,
        0.0,
        0.0,
        -(right + left) / rml,
        0.0,
        2.0 / tmb,
        0.0,
        -(top + bottom) / tmb,
        0.0,
        0.0,
        -2.0 / fmn,
        -(far + near) / fmn,
        0.0,
        0.0,
        0.0,
        1.0,
    )
}
/// Create a Vec4 from components.
pub fn vec4(x: Real, y: Real, z: Real, w: Real) -> Vec4 {
    Vec4::new(x, y, z, w)
}
/// Homogeneous point (w=1).
pub fn vec4_point(v: &Vec3) -> Vec4 {
    Vec4::new(v.x, v.y, v.z, 1.0)
}
/// Homogeneous direction (w=0).
pub fn vec4_direction(v: &Vec3) -> Vec4 {
    Vec4::new(v.x, v.y, v.z, 0.0)
}
/// Project a Vec4 back to Vec3 by dividing by w.
///
/// Returns `None` if w is near zero.
pub fn vec4_to_vec3(v: &Vec4) -> Option<Vec3> {
    if v.w.abs() < 1e-10 {
        return None;
    }
    Some(Vec3::new(v.x / v.w, v.y / v.w, v.z / v.w))
}
/// Compute the intersection line of two planes.
///
/// Returns `None` if planes are parallel.
/// On success, returns `(point_on_line, direction)`.
pub fn plane_plane_intersection(a: &Plane, b: &Plane) -> Option<(Vec3, Vec3)> {
    let dir = a.normal.cross(&b.normal);
    let len_sq = dir.norm_squared();
    if len_sq < 1e-10 {
        return None;
    }
    let p = (b.normal.cross(&dir) * a.distance + dir.cross(&a.normal) * b.distance) / len_sq;
    Some((p, dir.normalize()))
}
/// Compute the intersection point of three planes.
///
/// Returns `None` if any two planes are parallel or all three meet in a line.
pub fn three_plane_intersection(a: &Plane, b: &Plane, c: &Plane) -> Option<Vec3> {
    let denom = a.normal.dot(&b.normal.cross(&c.normal));
    if denom.abs() < 1e-10 {
        return None;
    }
    let p = (b.normal.cross(&c.normal) * a.distance
        + c.normal.cross(&a.normal) * b.distance
        + a.normal.cross(&b.normal) * c.distance)
        / denom;
    Some(p)
}
/// Compute the polar decomposition of a 3×3 matrix `M = R * S` where `R` is
/// orthogonal (rotation) and `S` is symmetric positive semi-definite.
///
/// Uses iterative polar decomposition: `R_{k+1} = 0.5 * (R_k + (R_k^{-T}))`.
/// Returns `(R, S)` after convergence or `max_iter` iterations.
///
/// Returns `None` if the matrix is singular.
pub fn polar_decomposition(m: &Mat3, max_iter: usize) -> Option<(Mat3, Mat3)> {
    let mut r = *m;
    for _ in 0..max_iter {
        let r_inv_t = r.try_inverse()?.transpose();
        let r_new = (r + r_inv_t) * 0.5;
        let diff = (r_new - r).norm();
        r = r_new;
        if diff < 1e-12 {
            break;
        }
    }
    let s = r.transpose() * m;
    Some((r, s))
}
/// Compute eigenvalues and eigenvectors of a real symmetric 3×3 matrix using
/// the Jacobi iteration method.
///
/// Returns `(eigenvalues, eigenvectors)` where each column of `eigenvectors` is
/// an eigenvector corresponding to the eigenvalue at the same index.
/// Eigenvalues are *not* sorted.
pub fn symmetric_eigen3(m: &Mat3) -> ([Real; 3], Mat3) {
    let mut a = *m;
    let mut v = Mat3::identity();
    for _ in 0..100 {
        let mut max_val = 0.0_f64;
        let mut p = 0;
        let mut q = 1;
        for i in 0..3 {
            for j in (i + 1)..3 {
                if a[(i, j)].abs() > max_val {
                    max_val = a[(i, j)].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < 1e-12 {
            break;
        }
        let theta = (a[(q, q)] - a[(p, p)]) / (2.0 * a[(p, q)]);
        let t = if theta >= 0.0 {
            1.0 / (theta + (1.0 + theta * theta).sqrt())
        } else {
            1.0 / (theta - (1.0 + theta * theta).sqrt())
        };
        let cos = 1.0 / (1.0 + t * t).sqrt();
        let sin = t * cos;
        let a_pp = a[(p, p)];
        let a_qq = a[(q, q)];
        let a_pq = a[(p, q)];
        a[(p, p)] = cos * cos * a_pp - 2.0 * sin * cos * a_pq + sin * sin * a_qq;
        a[(q, q)] = sin * sin * a_pp + 2.0 * sin * cos * a_pq + cos * cos * a_qq;
        a[(p, q)] = 0.0;
        a[(q, p)] = 0.0;
        for r in 0..3 {
            if r != p && r != q {
                let a_rp = a[(r, p)];
                let a_rq = a[(r, q)];
                a[(r, p)] = cos * a_rp - sin * a_rq;
                a[(p, r)] = a[(r, p)];
                a[(r, q)] = sin * a_rp + cos * a_rq;
                a[(q, r)] = a[(r, q)];
            }
        }
        for r in 0..3 {
            let v_rp = v[(r, p)];
            let v_rq = v[(r, q)];
            v[(r, p)] = cos * v_rp - sin * v_rq;
            v[(r, q)] = sin * v_rp + cos * v_rq;
        }
    }
    ([a[(0, 0)], a[(1, 1)], a[(2, 2)]], v)
}
/// Evaluate real spherical harmonic Y_0^0 (degree 0, order 0).
///
/// Y_0^0 = 1 / (2 * sqrt(π))
pub fn sh_y00() -> Real {
    0.5 / std::f64::consts::PI.sqrt()
}
/// Evaluate real spherical harmonic Y_1^{-1} (degree 1, order -1).
///
/// Y_1^{-1}(θ,φ) = sqrt(3/(4π)) * sin(θ) * sin(φ) = sqrt(3/(4π)) * y/r
pub fn sh_y1m1(dir: &Vec3) -> Real {
    let len = dir.norm();
    if len < 1e-12 {
        return 0.0;
    }
    let y = dir.y / len;
    (3.0 / (4.0 * std::f64::consts::PI)).sqrt() * y
}
/// Evaluate real spherical harmonic Y_1^0 (degree 1, order 0).
///
/// Y_1^0(θ,φ) = sqrt(3/(4π)) * cos(θ) = sqrt(3/(4π)) * z/r
pub fn sh_y10(dir: &Vec3) -> Real {
    let len = dir.norm();
    if len < 1e-12 {
        return 0.0;
    }
    let z = dir.z / len;
    (3.0 / (4.0 * std::f64::consts::PI)).sqrt() * z
}
/// Evaluate real spherical harmonic Y_1^1 (degree 1, order 1).
///
/// Y_1^1(θ,φ) = sqrt(3/(4π)) * sin(θ) * cos(φ) = sqrt(3/(4π)) * x/r
pub fn sh_y11(dir: &Vec3) -> Real {
    let len = dir.norm();
    if len < 1e-12 {
        return 0.0;
    }
    let x = dir.x / len;
    (3.0 / (4.0 * std::f64::consts::PI)).sqrt() * x
}
/// Project a radiance function sampled at `n` uniformly-spread directions onto
/// the 4 L0+L1 SH coefficients \[c00, c1m1, c10, c11\].
///
/// `radiance_fn` takes a unit direction and returns the sampled value.
/// Directions are constructed from a uniform icosahedron-inspired distribution.
pub fn sh_project_l1(radiance_fn: impl Fn(&Vec3) -> Real, n: usize) -> [Real; 4] {
    let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut coeffs = [0.0_f64; 4];
    let weight = 4.0 * std::f64::consts::PI / n as Real;
    for i in 0..n {
        let theta = (1.0 - 2.0 * (i as Real + 0.5) / n as Real).acos();
        let phi = 2.0 * std::f64::consts::PI * (i as Real / golden).fract();
        let dir = Vec3::new(
            theta.sin() * phi.cos(),
            theta.sin() * phi.sin(),
            theta.cos(),
        );
        let rad = radiance_fn(&dir);
        coeffs[0] += rad * sh_y00() * weight;
        coeffs[1] += rad * sh_y1m1(&dir) * weight;
        coeffs[2] += rad * sh_y10(&dir) * weight;
        coeffs[3] += rad * sh_y11(&dir) * weight;
    }
    coeffs
}
/// Compute the derivative of `f` at `x` using dual numbers.
///
/// `f` must be a function `Dual -> Dual`.
pub fn dual_differentiate(f: impl Fn(Dual) -> Dual, x: Real) -> Real {
    f(Dual::variable(x)).du
}
/// Hermite spline interpolation between two points `p0` and `p1` with
/// tangents `m0` and `m1`.  Parameter `t` runs from 0 (at p0) to 1 (at p1).
pub fn hermite_interpolate(p0: &Vec3, m0: &Vec3, p1: &Vec3, m1: &Vec3, t: Real) -> Vec3 {
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    p0 * h00 + m0 * h10 + p1 * h01 + m1 * h11
}
/// Catmull-Rom spline interpolation.
///
/// Interpolates between `p1` and `p2` using the four control points
/// `p0, p1, p2, p3`.  `t` ∈ \[0, 1\].
pub fn catmull_rom(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3, t: Real) -> Vec3 {
    let m1 = (p2 - p0) * 0.5;
    let m2 = (p3 - p1) * 0.5;
    hermite_interpolate(p1, &m1, p2, &m2, t)
}
/// Cubic Bézier curve.
///
/// Evaluates at parameter `t` ∈ \[0, 1\] given control points `p0..p3`.
pub fn bezier_cubic(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3, t: Real) -> Vec3 {
    let mt = 1.0 - t;
    p0 * (mt * mt * mt) + p1 * (3.0 * mt * mt * t) + p2 * (3.0 * mt * t * t) + p3 * (t * t * t)
}
/// B-spline basis functions for a clamped uniform cubic B-spline (degree 3).
///
/// `i` is the span index, `t` ∈ \[0, 1\] within the span.
/// Returns the four non-zero basis function values `[N_{i-3}, N_{i-2}, N_{i-1}, N_i]`.
pub fn bspline_basis3(t: Real) -> [Real; 4] {
    let t2 = t * t;
    let t3 = t2 * t;
    let mt = 1.0 - t;
    [
        mt * mt * mt / 6.0,
        (3.0 * t3 - 6.0 * t2 + 4.0) / 6.0,
        (-3.0 * t3 + 3.0 * t2 + 3.0 * t + 1.0) / 6.0,
        t3 / 6.0,
    ]
}
/// Transpose a 3×3 matrix stored as `[[f64;3\];3]` (row-major).
pub fn mat3_transpose(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}
/// Multiply a 3×3 matrix by a 3-vector: `m * v`.
pub fn mat3_mul_vec3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
/// Multiply two 3×3 matrices: `a * b`.
pub fn mat3_mul_mat3(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0_f64; 3]; 3];
    for (i, ri) in r.iter_mut().enumerate() {
        for (j, rij) in ri.iter_mut().enumerate() {
            for k in 0..3 {
                *rij += a[i][k] * b[k][j];
            }
        }
    }
    r
}
/// Determinant of a 3×3 matrix.
pub fn mat3_det(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Inverse of a 3×3 matrix.  Returns `None` if singular (|det| < 1e-14).
pub fn mat3_arr_inverse(m: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let det = mat3_det(m);
    if det.abs() < 1e-14 {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ])
}
/// Build a 3×3 matrix from three column vectors.
pub fn mat3_from_cols(c0: [f64; 3], c1: [f64; 3], c2: [f64; 3]) -> [[f64; 3]; 3] {
    [
        [c0[0], c1[0], c2[0]],
        [c0[1], c1[1], c2[1]],
        [c0[2], c1[2], c2[2]],
    ]
}
/// Outer product of two 3-vectors: `a ⊗ b` (a 3×3 matrix).
pub fn mat3_outer_product(a: [f64; 3], b: [f64; 3]) -> [[f64; 3]; 3] {
    [
        [a[0] * b[0], a[0] * b[1], a[0] * b[2]],
        [a[1] * b[0], a[1] * b[1], a[1] * b[2]],
        [a[2] * b[0], a[2] * b[1], a[2] * b[2]],
    ]
}
/// Create a unit quaternion `[x, y, z, w]` from an axis-angle pair.
///
/// `axis` should be a unit vector; `angle` is in radians.
pub fn quat_arr_from_axis_angle(axis: [f64; 3], angle: f64) -> [f64; 4] {
    let half = angle * 0.5;
    let s = half.sin();
    [axis[0] * s, axis[1] * s, axis[2] * s, half.cos()]
}
/// Multiply two quaternions `p * q` (Hamilton product).  Format: `[x, y, z, w]`.
pub fn quat_multiply(p: [f64; 4], q: [f64; 4]) -> [f64; 4] {
    let [px, py, pz, pw] = p;
    let [qx, qy, qz, qw] = q;
    [
        pw * qx + px * qw + py * qz - pz * qy,
        pw * qy - px * qz + py * qw + pz * qx,
        pw * qz + px * qy - py * qx + pz * qw,
        pw * qw - px * qx - py * qy - pz * qz,
    ]
}
/// Convert a unit quaternion `[x, y, z, w]` to a 3×3 rotation matrix (row-major).
pub fn quat_to_mat3(q: [f64; 4]) -> [[f64; 3]; 3] {
    let [x, y, z, w] = q;
    let x2 = x * x;
    let y2 = y * y;
    let z2 = z * z;
    let xy = x * y;
    let xz = x * z;
    let yz = y * z;
    let wx = w * x;
    let wy = w * y;
    let wz = w * z;
    [
        [1.0 - 2.0 * (y2 + z2), 2.0 * (xy - wz), 2.0 * (xz + wy)],
        [2.0 * (xy + wz), 1.0 - 2.0 * (x2 + z2), 2.0 * (yz - wx)],
        [2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (x2 + y2)],
    ]
}
/// Spherical linear interpolation between two unit quaternions (format: `[x,y,z,w]`).
pub fn quat_arr_slerp(p: [f64; 4], q: [f64; 4], t: f64) -> [f64; 4] {
    let dot = p[0] * q[0] + p[1] * q[1] + p[2] * q[2] + p[3] * q[3];
    let (q2, dot2) = if dot < 0.0 {
        ([-q[0], -q[1], -q[2], -q[3]], -dot)
    } else {
        (q, dot)
    };
    let dot2 = dot2.clamp(-1.0, 1.0);
    if dot2 > 0.9995 {
        let r = [
            p[0] + t * (q2[0] - p[0]),
            p[1] + t * (q2[1] - p[1]),
            p[2] + t * (q2[2] - p[2]),
            p[3] + t * (q2[3] - p[3]),
        ];
        let norm = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2] + r[3] * r[3]).sqrt();
        [r[0] / norm, r[1] / norm, r[2] / norm, r[3] / norm]
    } else {
        let theta = dot2.acos();
        let sin_theta = theta.sin();
        let s0 = ((1.0 - t) * theta).sin() / sin_theta;
        let s1 = (t * theta).sin() / sin_theta;
        [
            s0 * p[0] + s1 * q2[0],
            s0 * p[1] + s1 * q2[1],
            s0 * p[2] + s1 * q2[2],
            s0 * p[3] + s1 * q2[3],
        ]
    }
}
/// Project vector `a` onto `onto`.  Returns the zero vector if `onto` is zero.
pub fn vec3_project(a: [f64; 3], onto: [f64; 3]) -> [f64; 3] {
    let denom = onto[0] * onto[0] + onto[1] * onto[1] + onto[2] * onto[2];
    if denom < 1e-30 {
        return [0.0; 3];
    }
    let s = (a[0] * onto[0] + a[1] * onto[1] + a[2] * onto[2]) / denom;
    [onto[0] * s, onto[1] * s, onto[2] * s]
}
/// Reflect vector `v` about unit normal `n`: `v - 2*(v·n)*n`.
pub fn vec3_reflect(v: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    let dot2 = 2.0 * (v[0] * n[0] + v[1] * n[1] + v[2] * n[2]);
    [v[0] - dot2 * n[0], v[1] - dot2 * n[1], v[2] - dot2 * n[2]]
}
/// Linear interpolation between two 3-vectors: `a + t*(b-a)`.
pub fn vec3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + t * (b[0] - a[0]),
        a[1] + t * (b[1] - a[1]),
        a[2] + t * (b[2] - a[2]),
    ]
}
/// Create a plane `[nx, ny, nz, d]` from a point `p` and unit normal `n`.
///
/// The plane equation is `n·x = d` (equivalently `n·x - d = 0`).
pub fn plane_from_point_normal(p: [f64; 3], n: [f64; 3]) -> [f64; 4] {
    let d = n[0] * p[0] + n[1] * p[1] + n[2] * p[2];
    [n[0], n[1], n[2], d]
}
/// Signed distance from `point` to the plane `[nx, ny, nz, d]`.
///
/// Positive on the side the normal points toward.
pub fn plane_signed_dist(plane: [f64; 4], point: [f64; 3]) -> f64 {
    plane[0] * point[0] + plane[1] * point[1] + plane[2] * point[2] - plane[3]
}
/// Arc length of a cubic Bézier curve approximated with `n` segments.
pub fn bezier_arc_length(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3, n: usize) -> Real {
    let mut len = 0.0;
    let mut prev = *p0;
    for i in 1..=n {
        let t = i as Real / n as Real;
        let pt = bezier_cubic(p0, p1, p2, p3, t);
        len += (pt - prev).norm();
        prev = pt;
    }
    len
}
/// Cross product of two 3-vectors: `a × b`.
pub fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Scalar triple product: `a · (b × c)`.
///
/// Equals the signed volume of the parallelepiped spanned by `a`, `b`, `c`.
pub fn vec3_triple_product(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let bc = vec3_cross(b, c);
    a[0] * bc[0] + a[1] * bc[1] + a[2] * bc[2]
}
/// Rotate vector `v` around unit `axis` by `angle` radians (Rodrigues' formula).
///
/// `axis` must be a unit vector; the result has the same magnitude as `v`.
pub fn vec3_rotate_by_angle(v: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let dot = v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2];
    let cross = vec3_cross(axis, v);
    [
        cos_a * v[0] + sin_a * cross[0] + (1.0 - cos_a) * dot * axis[0],
        cos_a * v[1] + sin_a * cross[1] + (1.0 - cos_a) * dot * axis[1],
        cos_a * v[2] + sin_a * cross[2] + (1.0 - cos_a) * dot * axis[2],
    ]
}
/// Adjugate (classical adjoint) of a 3×3 matrix: `adj(M) = det(M) * M^{-1}`.
///
/// The adjugate exists even for singular matrices.
pub fn mat3_adjugate(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let c = |r: usize, c_idx: usize| -> f64 {
        let rows: Vec<usize> = (0..3).filter(|&x| x != r).collect();
        let cols: Vec<usize> = (0..3).filter(|&x| x != c_idx).collect();
        let det2 =
            m[rows[0]][cols[0]] * m[rows[1]][cols[1]] - m[rows[0]][cols[1]] * m[rows[1]][cols[0]];
        let sign = if (r + c_idx).is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        sign * det2
    };
    [
        [c(0, 0), c(1, 0), c(2, 0)],
        [c(0, 1), c(1, 1), c(2, 1)],
        [c(0, 2), c(1, 2), c(2, 2)],
    ]
}
/// Verify the Cayley-Hamilton theorem for a 3×3 matrix M.
///
/// Computes `p(M) = M^3 - tr(M)*M^2 + ((tr(M)^2 - tr(M^2))/2)*M - det(M)*I`
/// and returns the result (should be the zero matrix by Cayley-Hamilton).
pub fn mat3_cayley_hamilton(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let i3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let tr_m = m[0][0] + m[1][1] + m[2][2];
    let m2 = mat3_mul_mat3(m, m);
    let m3 = mat3_mul_mat3(m2, m);
    let tr_m2 = m2[0][0] + m2[1][1] + m2[2][2];
    let c2 = (tr_m * tr_m - tr_m2) * 0.5;
    let det = mat3_det(m);
    let mut result = [[0.0_f64; 3]; 3];
    for (i, (ri, ((m3i, m2i), (mi, i3i)))) in result
        .iter_mut()
        .zip(m3.iter().zip(m2.iter()).zip(m.iter().zip(i3.iter())))
        .enumerate()
    {
        for (j, (rij, ((&m3ij, &m2ij), (&mij, &i3ij)))) in ri
            .iter_mut()
            .zip(m3i.iter().zip(m2i.iter()).zip(mi.iter().zip(i3i.iter())))
            .enumerate()
        {
            let _ = (i, j);
            *rij = m3ij - tr_m * m2ij + c2 * mij - det * i3ij;
        }
    }
    result
}
/// Quaternion logarithm for a unit quaternion `q = [x, y, z, w]`.
///
/// `log(q) = [θ * n̂, 0]` where `θ = acos(w)` and `n̂ = (x,y,z)/sin(θ)`.
/// For the identity quaternion the result is the zero vector.
pub fn quat_arr_log(q: [f64; 4]) -> [f64; 4] {
    let [x, y, z, w] = q;
    let w_clamped = w.clamp(-1.0, 1.0);
    let theta = w_clamped.acos();
    let sin_theta = theta.sin();
    if sin_theta.abs() < 1e-10 {
        return [0.0, 0.0, 0.0, 0.0];
    }
    let s = theta / sin_theta;
    [x * s, y * s, z * s, 0.0]
}
/// Quaternion exponential for a pure quaternion `v = [x, y, z, 0]`.
///
/// `exp(v) = [sin(|v|)*v̂, cos(|v|)]`.  If `v = 0` the identity is returned.
pub fn quat_arr_exp(v: [f64; 4]) -> [f64; 4] {
    let [x, y, z, _] = v;
    let theta = (x * x + y * y + z * z).sqrt();
    if theta < 1e-15 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    let s = theta.sin() / theta;
    [x * s, y * s, z * s, theta.cos()]
}
/// Convert a unit quaternion `[x, y, z, w]` to `(axis, angle)`.
///
/// Returns `([0,1,0], 0)` for the identity quaternion.
pub fn quat_to_axis_angle(q: [f64; 4]) -> ([f64; 3], f64) {
    let [x, y, z, w] = q;
    let w_c = w.clamp(-1.0, 1.0);
    let angle = 2.0 * w_c.acos();
    let sin_half = (1.0 - w_c * w_c).sqrt();
    if sin_half < 1e-10 {
        return ([0.0, 1.0, 0.0], 0.0);
    }
    ([x / sin_half, y / sin_half, z / sin_half], angle)
}
/// SQUAD (Spherical Quadrangle interpolation) for smooth quaternion paths.
///
/// Interpolates between `q1` and `q2` at parameter `t ∈ [0,1]` using the
/// surrounding control quaternions `s1` and `s2`.
pub fn quat_squad(q1: [f64; 4], q2: [f64; 4], s1: [f64; 4], s2: [f64; 4], t: f64) -> [f64; 4] {
    let slerp_q = quat_arr_slerp(q1, q2, t);
    let slerp_s = quat_arr_slerp(s1, s2, t);
    quat_arr_slerp(slerp_q, slerp_s, 2.0 * t * (1.0 - t))
}
#[cfg(test)]
mod proptest_tests {

    use crate::Vec3;

    use crate::math::types::{Plane, Ray};

    use proptest::prelude::*;
    fn coord() -> impl Strategy<Value = f64> {
        -1e3_f64..1e3_f64
    }
    fn pos_coord() -> impl Strategy<Value = f64> {
        0.01_f64..100.0_f64
    }
    fn vec3() -> impl Strategy<Value = Vec3> {
        (coord(), coord(), coord()).prop_map(|(x, y, z)| Vec3::new(x, y, z))
    }
    fn nonzero_vec3() -> impl Strategy<Value = Vec3> {
        (coord(), coord(), coord())
            .prop_filter("norm must be nonzero", |(x, y, z)| {
                (x * x + y * y + z * z).sqrt() > 1e-6
            })
            .prop_map(|(x, y, z)| Vec3::new(x, y, z))
    }
    proptest! {
        #[test] fn prop_vec3_dot_commutative(a in vec3(), b in vec3()) { let ab = a.dot(&
        b); let ba = b.dot(& a); prop_assert!((ab - ba).abs() < 1e-10,
        "dot not commutative: {} vs {}", ab, ba); } #[test] fn
        prop_vec3_cross_anti_commutative(a in vec3(), b in vec3()) { let ab = a.cross(&
        b); let ba = b.cross(& a); let diff = (ab + ba).norm(); prop_assert!(diff < 1e-6,
        "cross not anti-commutative: diff norm={}", diff); } #[test] fn
        prop_vec3_triangle_inequality(a in vec3(), b in vec3()) { let sum_norm = (a + b)
        .norm(); let norm_sum = a.norm() + b.norm(); prop_assert!(sum_norm <= norm_sum +
        1e-10, "triangle inequality violated: |a+b|={} > |a|+|b|={}", sum_norm,
        norm_sum); } #[test] fn prop_vec3_normalize_unit(v in nonzero_vec3()) { let n = v
        .normalize(); let len = n.norm(); prop_assert!((len - 1.0).abs() < 1e-10,
        "normalized vector has norm {} != 1", len); } #[test] fn
        prop_ray_point_at_parametric(ox in coord(), oy in coord(), oz in coord(), dx in
        coord(), dy in coord(), dz in coord(), t in 0.0_f64..10.0_f64,) { let origin =
        Vec3::new(ox, oy, oz); let direction = Vec3::new(dx, dy, dz); let ray =
        Ray::new(origin, direction); let p = ray.point_at(t); let expected = origin +
        direction * t; let diff = (p - expected).norm(); prop_assert!(diff < 1e-10,
        "ray.point_at mismatch: diff={}", diff); } #[test] fn
        prop_plane_signed_distance_positive_side(nx in pos_coord(), ny in pos_coord(), nz
        in pos_coord(), d in coord(), s in 0.01_f64..10.0_f64,) { let n = Vec3::new(nx,
        ny, nz).normalize(); let plane = Plane::new(n, d); let point = n * (d + s); let
        dist = plane.signed_distance(& point); prop_assert!((dist - s).abs() < 1e-6,
        "expected signed distance {}, got {}", s, dist); }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::Mat3;
    use crate::Vec3;
    use crate::differential_geometry::mat3_det;
    use crate::differential_geometry::mat3_mul_vec3;
    use crate::differential_geometry::mat3_transpose;
    use crate::math::Vec4;
    use crate::math::bezier_arc_length;
    use crate::math::bezier_cubic;
    use crate::math::bspline_basis3;
    use crate::math::dual_differentiate;
    use crate::math::hermite_interpolate;
    use crate::math::look_at;
    use crate::math::mat3_adjugate;
    use crate::math::mat3_arr_inverse;
    use crate::math::mat3_cayley_hamilton;
    use crate::math::mat3_determinant;
    use crate::math::mat3_eigenvalues_symmetric;
    use crate::math::mat3_inverse;
    use crate::math::mat3_mul_mat3;
    use crate::math::mat3_outer_product;
    use crate::math::orthographic;
    use crate::math::perspective;
    use crate::math::plane_from_point_normal;
    use crate::math::plane_plane_intersection;
    use crate::math::plane_signed_dist;
    use crate::math::polar_decomposition;
    use crate::math::quat_arr_exp;
    use crate::math::quat_arr_from_axis_angle;
    use crate::math::quat_arr_log;
    use crate::math::quat_arr_slerp;
    use crate::math::quat_from_axis_angle;
    use crate::math::quat_from_euler;
    use crate::math::quat_multiply;
    use crate::math::quat_slerp;
    use crate::math::quat_to_axis_angle;
    use crate::math::quat_to_euler;
    use crate::math::quat_to_mat3;
    use crate::math::sh_project_l1;
    use crate::math::sh_y00;
    use crate::math::sh_y10;
    use crate::math::skew_symmetric;
    use crate::math::symmetric_eigen3;
    use crate::math::three_plane_intersection;
    use crate::math::types::{Aabb, Dual, Frustum, Plane, Ray};
    use crate::math::vec3_cross;
    use crate::math::vec3_lerp;
    use crate::math::vec3_project;
    use crate::math::vec3_reflect;
    use crate::math::vec3_rotate_by_angle;
    use crate::math::vec3_triple_product;
    use crate::math::vec4_direction;
    use crate::math::vec4_point;
    use crate::math::vec4_to_vec3;
    #[test]
    fn test_ray_point_at() {
        let ray = Ray::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let p = ray.point_at(3.0);
        assert!((p.x - 3.0).abs() < 1e-10);
        assert!(p.y.abs() < 1e-10);
    }
    #[test]
    fn test_plane_signed_distance() {
        let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        assert!((plane.signed_distance(&Vec3::new(0.0, 5.0, 0.0)) - 5.0).abs() < 1e-10);
        assert!((plane.signed_distance(&Vec3::new(0.0, -3.0, 0.0)) + 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_ray_plane_intersection() {
        let ray = Ray::new(Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
        let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        let t = ray.intersect_plane(&plane).expect("should intersect");
        assert!((t - 1.0).abs() < 1e-10, "expected t=1.0, got {}", t);
        let hit = ray.point_at(t);
        assert!(hit.y.abs() < 1e-10);
    }
    #[test]
    fn test_ray_plane_miss() {
        let ray = Ray::new(Vec3::new(0.0, 1.0, 0.0), Vec3::new(1.0, 0.0, 0.0));
        let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        assert!(ray.intersect_plane(&plane).is_none());
    }
    #[test]
    fn test_plane_project() {
        let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        let point = Vec3::new(3.0, 5.0, 7.0);
        let d = plane.signed_distance(&point);
        let projected = point - plane.normal * d;
        assert!(plane.signed_distance(&projected).abs() < 1e-10);
    }
    #[test]
    fn test_quat_from_axis_angle_identity() {
        let q = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.0);
        let (roll, pitch, yaw) = quat_to_euler(&q);
        assert!(roll.abs() < 1e-10);
        assert!(pitch.abs() < 1e-10);
        assert!(yaw.abs() < 1e-10);
    }
    #[test]
    fn test_quat_from_axis_angle_90_degrees() {
        let q = quat_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_2);
        let v = q.transform_vector(&Vec3::new(1.0, 0.0, 0.0));
        assert!((v.x).abs() < 1e-10, "x={}", v.x);
        assert!((v.y - 1.0).abs() < 1e-10, "y={}", v.y);
        assert!((v.z).abs() < 1e-10, "z={}", v.z);
    }
    #[test]
    fn test_quat_euler_roundtrip() {
        let roll = 0.3;
        let pitch = 0.5;
        let yaw = 0.7;
        let q = quat_from_euler(roll, pitch, yaw);
        let (r2, p2, y2) = quat_to_euler(&q);
        assert!((r2 - roll).abs() < 1e-10, "roll: {} vs {}", r2, roll);
        assert!((p2 - pitch).abs() < 1e-10, "pitch: {} vs {}", p2, pitch);
        assert!((y2 - yaw).abs() < 1e-10, "yaw: {} vs {}", y2, yaw);
    }
    #[test]
    fn test_quat_slerp_endpoints() {
        let a = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.0);
        let b = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 1.0);
        let q0 = quat_slerp(&a, &b, 0.0);
        let q1 = quat_slerp(&a, &b, 1.0);
        assert!(q0.angle_to(&a) < 1e-10);
        assert!(q1.angle_to(&b) < 1e-10);
    }
    #[test]
    fn test_quat_slerp_midpoint() {
        let a = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.0);
        let b = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 1.0);
        let mid = quat_slerp(&a, &b, 0.5);
        let expected = quat_from_axis_angle(&Vec3::new(0.0, 1.0, 0.0), 0.5);
        assert!(mid.angle_to(&expected) < 1e-10);
    }
    #[test]
    fn test_mat3_determinant_identity() {
        let m = Mat3::identity();
        assert!((mat3_determinant(&m) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat3_determinant_scaled() {
        let m = Mat3::identity() * 2.0;
        assert!((mat3_determinant(&m) - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat3_inverse_identity() {
        let m = Mat3::identity();
        let inv = mat3_inverse(&m).expect("identity should be invertible");
        let diff = (inv - Mat3::identity()).norm();
        assert!(diff < 1e-10);
    }
    #[test]
    fn test_mat3_inverse_product_is_identity() {
        let m = Mat3::new(1.0, 2.0, 3.0, 0.0, 1.0, 4.0, 5.0, 6.0, 0.0);
        let inv = mat3_inverse(&m).expect("should be invertible");
        let product = m * inv;
        let diff = (product - Mat3::identity()).norm();
        assert!(diff < 1e-8, "M * M^-1 should be identity, diff={}", diff);
    }
    #[test]
    fn test_mat3_singular_no_inverse() {
        let m = Mat3::new(1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 1.0, 1.0, 1.0);
        assert!(mat3_inverse(&m).is_none());
    }
    #[test]
    fn test_mat3_eigenvalues_diagonal() {
        let m = Mat3::new(3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0);
        let eigs = mat3_eigenvalues_symmetric(&m);
        assert!((eigs[0] - 1.0).abs() < 1e-10);
        assert!((eigs[1] - 2.0).abs() < 1e-10);
        assert!((eigs[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat3_trace() {
        let m = Mat3::new(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0);
        assert!((mat3_trace(&m) - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat3_frobenius_norm_identity() {
        let m = Mat3::identity();
        assert!((mat3_frobenius_norm(&m) - 3.0_f64.sqrt()).abs() < 1e-10);
    }
    #[test]
    fn test_skew_symmetric_cross_product() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        let cross_direct = a.cross(&b);
        let cross_matrix = skew_symmetric(&a) * b;
        let diff = (cross_direct - cross_matrix).norm();
        assert!(
            diff < 1e-10,
            "skew_symmetric cross product mismatch: diff={}",
            diff
        );
    }
    #[test]
    fn test_perspective_nonzero() {
        let m = perspective(std::f64::consts::FRAC_PI_4, 16.0 / 9.0, 0.1, 100.0);
        assert!(m.norm() > 0.0);
        assert!((m[(3, 0)]).abs() < 1e-10);
        assert!((m[(3, 1)]).abs() < 1e-10);
        assert!((m[(3, 2)] + 1.0).abs() < 1e-10);
        assert!((m[(3, 3)]).abs() < 1e-10);
    }
    #[test]
    fn test_look_at_origin() {
        let eye = Vec3::new(0.0, 0.0, 5.0);
        let target = Vec3::new(0.0, 0.0, 0.0);
        let up = Vec3::new(0.0, 1.0, 0.0);
        let m = look_at(&eye, &target, &up);
        let p = m * vec4_point(&Vec3::zeros());
        let v3 = vec4_to_vec3(&p).unwrap();
        assert!((v3.x).abs() < 1e-8, "x={}", v3.x);
        assert!((v3.y).abs() < 1e-8, "y={}", v3.y);
        assert!((v3.z + 5.0).abs() < 1e-8, "z={}", v3.z);
    }
    #[test]
    fn test_orthographic_identity_like() {
        let m = orthographic(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
        let p = m * vec4_point(&Vec3::zeros());
        let v3 = vec4_to_vec3(&p).unwrap();
        assert!(v3.norm() < 1e-8);
    }
    #[test]
    fn test_vec4_point_w_is_one() {
        let p = vec4_point(&Vec3::new(1.0, 2.0, 3.0));
        assert!((p.w - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_vec4_direction_w_is_zero() {
        let d = vec4_direction(&Vec3::new(1.0, 0.0, 0.0));
        assert!(d.w.abs() < 1e-10);
    }
    #[test]
    fn test_vec4_to_vec3_perspective_divide() {
        let v = Vec4::new(4.0, 6.0, 8.0, 2.0);
        let v3 = vec4_to_vec3(&v).unwrap();
        assert!((v3.x - 2.0).abs() < 1e-10);
        assert!((v3.y - 3.0).abs() < 1e-10);
        assert!((v3.z - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_vec4_to_vec3_w_zero_returns_none() {
        let v = Vec4::new(1.0, 2.0, 3.0, 0.0);
        assert!(vec4_to_vec3(&v).is_none());
    }
    #[test]
    fn test_plane_from_points() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        let c = Vec3::new(0.0, 1.0, 0.0);
        let plane = Plane::from_points(&a, &b, &c).unwrap();
        assert!((plane.normal.z.abs() - 1.0).abs() < 1e-10);
        assert!(plane.signed_distance(&a).abs() < 1e-10);
        assert!(plane.signed_distance(&b).abs() < 1e-10);
        assert!(plane.signed_distance(&c).abs() < 1e-10);
    }
    #[test]
    fn test_plane_from_collinear_returns_none() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        let c = Vec3::new(2.0, 0.0, 0.0);
        assert!(Plane::from_points(&a, &b, &c).is_none());
    }
    #[test]
    fn test_plane_project_point() {
        let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        let point = Vec3::new(3.0, 5.0, 7.0);
        let projected = plane.project_point(&point);
        assert!(plane.signed_distance(&projected).abs() < 1e-10);
        assert!((projected.x - 3.0).abs() < 1e-10);
        assert!((projected.z - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_plane_reflect_point() {
        let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        let point = Vec3::new(3.0, 5.0, 7.0);
        let reflected = plane.reflect_point(&point);
        assert!((reflected.x - 3.0).abs() < 1e-10);
        assert!((reflected.y + 5.0).abs() < 1e-10);
        assert!((reflected.z - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_plane_classify_point() {
        let plane = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        assert_eq!(plane.classify_point(&Vec3::new(0.0, 1.0, 0.0), 1e-6), 1);
        assert_eq!(plane.classify_point(&Vec3::new(0.0, -1.0, 0.0), 1e-6), -1);
        assert_eq!(plane.classify_point(&Vec3::new(0.0, 0.0, 0.0), 1e-6), 0);
    }
    #[test]
    fn test_plane_plane_intersection() {
        let a = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        let b = Plane::new(Vec3::new(1.0, 0.0, 0.0), 0.0);
        let (point, dir) = plane_plane_intersection(&a, &b).unwrap();
        assert!(a.signed_distance(&point).abs() < 1e-8);
        assert!(b.signed_distance(&point).abs() < 1e-8);
        assert!((dir.z.abs() - 1.0).abs() < 1e-8, "dir={:?}", dir);
    }
    #[test]
    fn test_plane_plane_parallel_returns_none() {
        let a = Plane::new(Vec3::new(0.0, 1.0, 0.0), 0.0);
        let b = Plane::new(Vec3::new(0.0, 1.0, 0.0), 5.0);
        assert!(plane_plane_intersection(&a, &b).is_none());
    }
    #[test]
    fn test_three_plane_intersection() {
        let a = Plane::new(Vec3::new(1.0, 0.0, 0.0), 1.0);
        let b = Plane::new(Vec3::new(0.0, 1.0, 0.0), 2.0);
        let c = Plane::new(Vec3::new(0.0, 0.0, 1.0), 3.0);
        let p = three_plane_intersection(&a, &b, &c).unwrap();
        assert!((p.x - 1.0).abs() < 1e-8);
        assert!((p.y - 2.0).abs() < 1e-8);
        assert!((p.z - 3.0).abs() < 1e-8);
    }
    #[test]
    fn test_polar_decomposition_identity() {
        let m = Mat3::identity();
        let (r, s) = polar_decomposition(&m, 100).expect("identity is invertible");
        let diff_r = (r - Mat3::identity()).norm();
        let diff_s = (s - Mat3::identity()).norm();
        assert!(diff_r < 1e-8, "R should be identity: diff={}", diff_r);
        assert!(diff_s < 1e-8, "S should be identity: diff={}", diff_s);
    }
    #[test]
    fn test_polar_decomposition_rotation_only() {
        let q = quat_from_axis_angle(&Vec3::new(0.0, 0.0, 1.0), std::f64::consts::FRAC_PI_4);
        let m = *q.to_rotation_matrix().matrix();
        let (r, s) = polar_decomposition(&m, 100).expect("rotation is invertible");
        let product = r * s;
        let diff = (product - m).norm();
        assert!(diff < 1e-6, "R*S should equal M: diff={}", diff);
        let s_diff = (s - Mat3::identity()).norm();
        assert!(s_diff < 1e-6, "S should be near identity: diff={}", s_diff);
    }
    #[test]
    fn test_polar_decomposition_product_equals_m() {
        let m = Mat3::new(1.0, 0.5, 0.0, 0.0, 1.0, 0.3, 0.0, 0.0, 1.0);
        let (r, s) = polar_decomposition(&m, 100).expect("M should be invertible");
        let product = r * s;
        let diff = (product - m).norm();
        assert!(diff < 1e-6, "R*S should equal M: diff={}", diff);
    }
    #[test]
    fn test_symmetric_eigen3_diagonal() {
        let m = Mat3::new(3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0);
        let (eigs, _) = symmetric_eigen3(&m);
        let mut sorted = eigs;
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((sorted[0] - 1.0).abs() < 1e-8, "eig0={}", sorted[0]);
        assert!((sorted[1] - 2.0).abs() < 1e-8, "eig1={}", sorted[1]);
        assert!((sorted[2] - 3.0).abs() < 1e-8, "eig2={}", sorted[2]);
    }
    #[test]
    fn test_symmetric_eigen3_eigenvectors_orthonormal() {
        let m = Mat3::new(4.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 2.0);
        let (_, v) = symmetric_eigen3(&m);
        for i in 0..3 {
            let col_i = v.column(i);
            let norm = col_i.norm();
            assert!((norm - 1.0).abs() < 1e-6, "col {} norm={}", i, norm);
            for j in (i + 1)..3 {
                let col_j = v.column(j);
                let dot = col_i.dot(&col_j);
                assert!(
                    dot.abs() < 1e-6,
                    "cols {} {} not orthogonal: dot={}",
                    i,
                    j,
                    dot
                );
            }
        }
    }
    #[test]
    fn test_symmetric_eigen3_reconstruction() {
        let m = Mat3::new(2.0, 1.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 5.0);
        let (eigs, v) = symmetric_eigen3(&m);
        let diag = Mat3::from_diagonal(&nalgebra::Vector3::new(eigs[0], eigs[1], eigs[2]));
        let reconstructed = v * diag * v.transpose();
        let diff = (reconstructed - m).norm();
        assert!(diff < 1e-6, "reconstruction diff={}", diff);
    }
    #[test]
    fn test_sh_y00_constant() {
        let expected = 0.5 / std::f64::consts::PI.sqrt();
        assert!((sh_y00() - expected).abs() < 1e-12);
    }
    #[test]
    fn test_sh_y1_orthogonal_to_constant() {
        let n = 1000;
        let golden = (1.0 + 5.0_f64.sqrt()) / 2.0;
        let mut sum_y10 = 0.0_f64;
        for i in 0..n {
            let theta = (1.0 - 2.0 * (i as f64 + 0.5) / n as f64).acos();
            let phi = 2.0 * std::f64::consts::PI * ((i as f64) / golden).fract();
            let dir = Vec3::new(
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            );
            sum_y10 += sh_y10(&dir);
        }
        let avg = sum_y10 / n as f64;
        assert!(avg.abs() < 0.01, "Y_1^0 average should be ~0, got {}", avg);
    }
    #[test]
    fn test_sh_project_l1_uniform_radiance() {
        let coeffs = sh_project_l1(|_| 1.0, 1000);
        let expected_c00 = 2.0 * std::f64::consts::PI.sqrt();
        assert!((coeffs[0] - expected_c00).abs() < 0.1, "c00={}", coeffs[0]);
        for &c in &coeffs[1..] {
            assert!(c.abs() < 0.1, "L1 coeff should be ~0 for uniform: {}", c);
        }
    }
    #[test]
    fn test_dual_add() {
        let a = Dual::new(3.0, 1.0);
        let b = Dual::new(2.0, 0.5);
        let c = a + b;
        assert!((c.re - 5.0).abs() < 1e-12);
        assert!((c.du - 1.5).abs() < 1e-12);
    }
    #[test]
    fn test_dual_mul() {
        let x = Dual::variable(3.0);
        let y = x * x;
        assert!((y.re - 9.0).abs() < 1e-12);
        assert!((y.du - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_dual_differentiate_polynomial() {
        let deriv = dual_differentiate(|x| x.powi(3) - Dual::constant(2.0) * x, 2.0);
        assert!((deriv - 10.0).abs() < 1e-10, "deriv={}", deriv);
    }
    #[test]
    fn test_dual_sin_derivative() {
        let x = std::f64::consts::FRAC_PI_4;
        let deriv = dual_differentiate(|d| d.sin(), x);
        assert!((deriv - x.cos()).abs() < 1e-12, "deriv={}", deriv);
    }
    #[test]
    fn test_dual_exp_derivative() {
        let deriv = dual_differentiate(|d| d.exp(), 1.0);
        assert!(
            (deriv - std::f64::consts::E).abs() < 1e-10,
            "deriv={}",
            deriv
        );
    }
    #[test]
    fn test_dual_sqrt_derivative() {
        let deriv = dual_differentiate(|d| d.sqrt(), 4.0);
        assert!((deriv - 0.25).abs() < 1e-10, "deriv={}", deriv);
    }
    #[test]
    fn test_hermite_interpolate_endpoints() {
        let p0 = Vec3::new(0.0, 0.0, 0.0);
        let p1 = Vec3::new(1.0, 0.0, 0.0);
        let m0 = Vec3::new(0.0, 1.0, 0.0);
        let m1 = Vec3::new(0.0, 1.0, 0.0);
        let at_0 = hermite_interpolate(&p0, &m0, &p1, &m1, 0.0);
        let at_1 = hermite_interpolate(&p0, &m0, &p1, &m1, 1.0);
        assert!((at_0 - p0).norm() < 1e-12, "t=0 should be p0");
        assert!((at_1 - p1).norm() < 1e-12, "t=1 should be p1");
    }
    #[test]
    fn test_bezier_cubic_endpoints() {
        let p0 = Vec3::new(0.0, 0.0, 0.0);
        let p1 = Vec3::new(1.0, 2.0, 0.0);
        let p2 = Vec3::new(2.0, 2.0, 0.0);
        let p3 = Vec3::new(3.0, 0.0, 0.0);
        let at_0 = bezier_cubic(&p0, &p1, &p2, &p3, 0.0);
        let at_1 = bezier_cubic(&p0, &p1, &p2, &p3, 1.0);
        assert!((at_0 - p0).norm() < 1e-12);
        assert!((at_1 - p3).norm() < 1e-12);
    }
    #[test]
    fn test_catmull_rom_interpolates_endpoints() {
        let p0 = Vec3::new(-1.0, 0.0, 0.0);
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(1.0, 0.0, 0.0);
        let p3 = Vec3::new(2.0, 0.0, 0.0);
        let at_0 = catmull_rom(&p0, &p1, &p2, &p3, 0.0);
        let at_1 = catmull_rom(&p0, &p1, &p2, &p3, 1.0);
        assert!((at_0 - p1).norm() < 1e-12);
        assert!((at_1 - p2).norm() < 1e-12);
    }
    #[test]
    fn test_bspline_basis3_partition_of_unity() {
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let b = bspline_basis3(t);
            let sum: f64 = b.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "t={}: sum={}", t, sum);
        }
    }
    #[test]
    fn test_bezier_arc_length_straight_line() {
        let p0 = Vec3::new(0.0, 0.0, 0.0);
        let p1 = Vec3::new(1.0, 0.0, 0.0);
        let p2 = Vec3::new(2.0, 0.0, 0.0);
        let p3 = Vec3::new(3.0, 0.0, 0.0);
        let len = bezier_arc_length(&p0, &p1, &p2, &p3, 100);
        assert!((len - 3.0).abs() < 1e-6, "arc length={}", len);
    }
    #[test]
    fn test_frustum_contains_origin() {
        let vp = perspective(std::f64::consts::FRAC_PI_4, 1.0, 0.1, 100.0);
        let frustum = Frustum::from_view_projection(&vp);
        assert_eq!(frustum.planes.len(), 6);
        for plane in &frustum.planes {
            assert!((plane.normal.norm() - 1.0).abs() < 1e-6);
        }
    }
    #[test]
    fn test_frustum_sphere_outside() {
        let vp = perspective(std::f64::consts::FRAC_PI_4, 1.0, 0.1, 100.0);
        let frustum = Frustum::from_view_projection(&vp);
        let result = frustum.intersects_sphere(&Vec3::new(1000.0, 0.0, 0.0), 1.0);
        let _ = result;
    }
    #[test]
    fn test_mat3_transpose_roundtrip() {
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let t = mat3_transpose(m);
        let tt = mat3_transpose(t);
        for i in 0..3 {
            for j in 0..3 {
                assert!((m[i][j] - tt[i][j]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_mat3_mul_vec3_identity() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let v = [3.0, 4.0, 5.0];
        let result = mat3_mul_vec3(identity, v);
        for i in 0..3 {
            assert!((result[i] - v[i]).abs() < 1e-12);
        }
    }
    #[test]
    fn test_mat3_arr_inverse() {
        let m = [[1.0, 2.0, 0.0], [0.0, 1.0, 3.0], [0.0, 0.0, 1.0]];
        let inv = mat3_arr_inverse(m).expect("should be invertible");
        let prod = mat3_mul_mat3(m, inv);
        for (i, row) in prod.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < 1e-10, "prod[{i}][{j}]={}", val);
            }
        }
    }
    #[test]
    fn test_mat3_outer_product() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let op = mat3_outer_product(a, b);
        assert!((op[0][1] - 1.0).abs() < 1e-12);
        assert!(op[0][0].abs() < 1e-12);
        assert!(op[1][1].abs() < 1e-12);
    }
    #[test]
    fn test_quat_arr_from_axis_angle_to_mat3_roundtrip() {
        let q = quat_arr_from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let m = quat_to_mat3(q);
        let v = mat3_mul_vec3(m, [1.0, 0.0, 0.0]);
        assert!(v[0].abs() < 1e-10, "x={}", v[0]);
        assert!((v[1] - 1.0).abs() < 1e-10, "y={}", v[1]);
        assert!(v[2].abs() < 1e-10, "z={}", v[2]);
    }
    #[test]
    fn test_quat_multiply_identity() {
        let id = [0.0, 0.0, 0.0, 1.0];
        let q = quat_arr_from_axis_angle([0.0, 1.0, 0.0], 0.5);
        let result = quat_multiply(id, q);
        for i in 0..4 {
            assert!((result[i] - q[i]).abs() < 1e-12, "component {i}");
        }
    }
    #[test]
    fn test_quat_arr_slerp_at_t0_and_t1() {
        let p = quat_arr_from_axis_angle([0.0, 1.0, 0.0], 0.0);
        let q = quat_arr_from_axis_angle([0.0, 1.0, 0.0], 1.0);
        let s0 = quat_arr_slerp(p, q, 0.0);
        let s1 = quat_arr_slerp(p, q, 1.0);
        for i in 0..4 {
            assert!(
                (s0[i] - p[i]).abs() < 1e-10,
                "t=0 component {i}: got {}, expected {}",
                s0[i],
                p[i]
            );
            assert!(
                (s1[i] - q[i]).abs() < 1e-10,
                "t=1 component {i}: got {}, expected {}",
                s1[i],
                q[i]
            );
        }
    }
    #[test]
    fn test_vec3_reflect() {
        let v = [1.0, -1.0, 0.0];
        let n = [0.0, 1.0, 0.0];
        let r = vec3_reflect(v, n);
        assert!((r[0] - 1.0).abs() < 1e-12);
        assert!((r[1] - 1.0).abs() < 1e-12);
        assert!(r[2].abs() < 1e-12);
    }
    #[test]
    fn test_vec3_project() {
        let a = [3.0, 4.0, 0.0];
        let onto = [1.0, 0.0, 0.0];
        let proj = vec3_project(a, onto);
        assert!((proj[0] - 3.0).abs() < 1e-12);
        assert!(proj[1].abs() < 1e-12);
        assert!(proj[2].abs() < 1e-12);
    }
    #[test]
    fn test_vec3_lerp() {
        let a = [0.0, 0.0, 0.0];
        let b = [2.0, 4.0, 6.0];
        let mid = vec3_lerp(a, b, 0.5);
        assert!((mid[0] - 1.0).abs() < 1e-12);
        assert!((mid[1] - 2.0).abs() < 1e-12);
        assert!((mid[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_plane_from_point_normal_and_signed_dist() {
        let p = [0.0, 5.0, 0.0];
        let n = [0.0, 1.0, 0.0];
        let plane = plane_from_point_normal(p, n);
        let dist_above = plane_signed_dist(plane, [0.0, 8.0, 0.0]);
        let dist_below = plane_signed_dist(plane, [0.0, 3.0, 0.0]);
        assert!((dist_above - 3.0).abs() < 1e-12);
        assert!((dist_below + 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_cross() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [0.0_f64, 1.0, 0.0];
        let c = vec3_cross(a, b);
        assert!(c[0].abs() < 1e-12);
        assert!(c[1].abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_triple_product() {
        let a = [1.0_f64, 0.0, 0.0];
        let b = [0.0_f64, 1.0, 0.0];
        let c = [0.0_f64, 0.0, 1.0];
        assert!((vec3_triple_product(a, b, c) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_rotate_by_angle_around_z() {
        let v = [1.0_f64, 0.0, 0.0];
        let axis = [0.0_f64, 0.0, 1.0];
        let r = vec3_rotate_by_angle(v, axis, std::f64::consts::FRAC_PI_2);
        assert!(r[0].abs() < 1e-12, "x={}", r[0]);
        assert!((r[1] - 1.0).abs() < 1e-12, "y={}", r[1]);
        assert!(r[2].abs() < 1e-12, "z={}", r[2]);
    }
    #[test]
    fn test_vec3_rotate_zero_angle() {
        let v = [3.0_f64, 4.0, 5.0];
        let axis = [0.0_f64, 1.0, 0.0];
        let r = vec3_rotate_by_angle(v, axis, 0.0);
        for i in 0..3 {
            assert!((r[i] - v[i]).abs() < 1e-12, "component {i}");
        }
    }
    #[test]
    fn test_mat3_adjugate_identity() {
        let m = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let adj = mat3_adjugate(m);
        for (i, row) in adj.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!((val - exp).abs() < 1e-12, "adj[{i}][{j}]={}", val);
            }
        }
    }
    #[test]
    fn test_mat3_adjugate_relation() {
        let m = [[1.0_f64, 2.0, 0.0], [0.0, 1.0, 3.0], [0.0, 0.0, 1.0]];
        let adj = mat3_adjugate(m);
        let det = mat3_det(m);
        let prod = mat3_mul_mat3(m, adj);
        for (i, row) in prod.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let exp = if i == j { det } else { 0.0 };
                assert!(
                    (val - exp).abs() < 1e-10,
                    "prod[{i}][{j}]={}  exp={}",
                    val,
                    exp
                );
            }
        }
    }
    #[test]
    fn test_cayley_hamilton() {
        let m = [[1.0_f64, 2.0, 3.0], [0.0, 4.0, 5.0], [0.0, 0.0, 6.0]];
        let zero = mat3_cayley_hamilton(m);
        for (i, row) in zero.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(val.abs() < 1e-8, "CH[{i}][{j}]={}", val);
            }
        }
    }
    #[test]
    fn test_quat_log_exp_roundtrip() {
        let q = quat_arr_from_axis_angle([0.0, 1.0, 0.0], 0.8);
        let log_q = quat_arr_log(q);
        let exp_log_q = quat_arr_exp(log_q);
        for i in 0..4 {
            assert!((exp_log_q[i] - q[i]).abs() < 1e-10, "component {i}");
        }
    }
    #[test]
    fn test_quat_log_identity() {
        let id = [0.0_f64, 0.0, 0.0, 1.0];
        let log_id = quat_arr_log(id);
        for &v in &log_id {
            assert!(v.abs() < 1e-12, "log(identity) should be zero: {}", v);
        }
    }
    #[test]
    fn test_quat_to_axis_angle_roundtrip() {
        let angle = 1.2_f64;
        let axis = [0.0_f64, 1.0, 0.0];
        let q = quat_arr_from_axis_angle(axis, angle);
        let (out_axis, out_angle) = quat_to_axis_angle(q);
        assert!(
            (out_angle - angle).abs() < 1e-10,
            "angle mismatch: {} vs {}",
            out_angle,
            angle
        );
        for i in 0..3 {
            assert!(
                (out_axis[i] - axis[i]).abs() < 1e-10,
                "axis[{i}]: {} vs {}",
                out_axis[i],
                axis[i]
            );
        }
    }
    #[test]
    fn test_aabb_contains() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        assert!(aabb.contains([1.0, 1.0, 1.0]));
        assert!(!aabb.contains([3.0, 1.0, 1.0]));
    }
    #[test]
    fn test_aabb_union() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([-1.0, -1.0, -1.0], [0.5, 0.5, 0.5]);
        let u = a.union(&b);
        for i in 0..3 {
            assert!((u.min[i] - (-1.0)).abs() < 1e-12);
        }
        assert!((u.max[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_aabb_intersect() {
        let a = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let b = Aabb::new([1.0, 1.0, 1.0], [3.0, 3.0, 3.0]);
        let inter = a.intersect(&b).expect("should overlap");
        for i in 0..3 {
            assert!((inter.min[i] - 1.0).abs() < 1e-12);
            assert!((inter.max[i] - 2.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_aabb_no_intersect() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = Aabb::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(a.intersect(&b).is_none());
    }
    #[test]
    fn test_aabb_expand() {
        let a = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let exp = a.expand(0.5);
        for i in 0..3 {
            assert!((exp.min[i] - (-0.5)).abs() < 1e-12);
            assert!((exp.max[i] - 1.5).abs() < 1e-12);
        }
    }
    #[test]
    fn test_aabb_intersect_ray() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        let hit = aabb.intersect_ray([-1.0, 1.0, 1.0], [1.0, 0.0, 0.0]);
        assert!(hit.is_some(), "ray should hit AABB");
        let (t0, _t1) = hit.unwrap();
        assert!((t0 - 1.0).abs() < 1e-10, "t0={}", t0);
    }
    #[test]
    fn test_aabb_ray_miss() {
        let aabb = Aabb::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let hit = aabb.intersect_ray([0.0, 5.0, 0.5], [1.0, 0.0, 0.0]);
        assert!(hit.is_none(), "ray should miss AABB");
    }
}
/// Outer product of two Vec3 vectors: `a ⊗ b` as a nalgebra `Mat3`.
///
/// The (i,j) element equals `a[i] * b[j]`.
pub fn vec3_outer_product(a: &Vec3, b: &Vec3) -> Mat3 {
    Mat3::new(
        a.x * b.x,
        a.x * b.y,
        a.x * b.z,
        a.y * b.x,
        a.y * b.y,
        a.y * b.z,
        a.z * b.x,
        a.z * b.y,
        a.z * b.z,
    )
}
/// Build the skew-symmetric (cross-product) matrix for `v`.
///
/// Equivalent to `skew_symmetric` but with an explicit name.
/// `cross_matrix(v) * w == v.cross(&w)`.
pub fn cross_matrix(v: &Vec3) -> Mat3 {
    skew_symmetric(v)
}
/// Gram-Schmidt orthogonalization of three linearly-independent vectors.
///
/// Returns an orthonormal basis `(e0, e1, e2)` where:
/// - `e0` is `v0` normalized,
/// - `e1` is `v1` minus its projection onto `e0`, normalized,
/// - `e2` is `v2` minus its projections onto `e0` and `e1`, normalized.
///
/// Returns `None` if any intermediate vector becomes near-zero.
pub fn gram_schmidt(v0: &Vec3, v1: &Vec3, v2: &Vec3) -> Option<(Vec3, Vec3, Vec3)> {
    let e0_len = v0.norm();
    if e0_len < 1e-12 {
        return None;
    }
    let e0 = v0 / e0_len;
    let u1 = v1 - e0 * e0.dot(v1);
    let u1_len = u1.norm();
    if u1_len < 1e-12 {
        return None;
    }
    let e1 = u1 / u1_len;
    let u2 = v2 - e0 * e0.dot(v2) - e1 * e1.dot(v2);
    let u2_len = u2.norm();
    if u2_len < 1e-12 {
        return None;
    }
    let e2 = u2 / u2_len;
    Some((e0, e1, e2))
}
/// Convert Cartesian coordinates `(x, y, z)` to spherical `(r, theta, phi)`.
///
/// - `r` = radial distance ≥ 0
/// - `theta` = polar angle ∈ `[0, π]` (from +Z axis)
/// - `phi` = azimuthal angle ∈ `(-π, π]` (in XY plane from +X axis)
pub fn cartesian_to_spherical(v: &Vec3) -> (Real, Real, Real) {
    let r = v.norm();
    if r < 1e-300 {
        return (0.0, 0.0, 0.0);
    }
    let theta = (v.z / r).clamp(-1.0, 1.0).acos();
    let phi = v.y.atan2(v.x);
    (r, theta, phi)
}
/// Convert spherical coordinates `(r, theta, phi)` to Cartesian `(x, y, z)`.
///
/// - `r` – radial distance
/// - `theta` – polar angle from +Z
/// - `phi` – azimuthal angle from +X in XY plane
pub fn spherical_to_cartesian(r: Real, theta: Real, phi: Real) -> Vec3 {
    Vec3::new(
        r * theta.sin() * phi.cos(),
        r * theta.sin() * phi.sin(),
        r * theta.cos(),
    )
}
/// Convert Cartesian coordinates to cylindrical `(rho, phi, z)`.
///
/// - `rho` = radial distance in XY plane
/// - `phi` = azimuthal angle ∈ `(-π, π]`
/// - `z`   = height
pub fn cartesian_to_cylindrical(v: &Vec3) -> (Real, Real, Real) {
    let rho = (v.x * v.x + v.y * v.y).sqrt();
    let phi = v.y.atan2(v.x);
    (rho, phi, v.z)
}
/// Convert cylindrical coordinates `(rho, phi, z)` to Cartesian.
pub fn cylindrical_to_cartesian(rho: Real, phi: Real, z: Real) -> Vec3 {
    Vec3::new(rho * phi.cos(), rho * phi.sin(), z)
}
/// Rodrigues rotation: rotate vector `v` around unit `axis` by `angle` radians.
///
/// Equivalent to `vec3_rotate_by_angle` but operates on `Vec3` directly.
/// `axis` must be a unit vector.
pub fn rodrigues_rotate(v: &Vec3, axis: &Vec3, angle: Real) -> Vec3 {
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    v * cos_a + axis.cross(v) * sin_a + axis * axis.dot(v) * (1.0 - cos_a)
}
/// Matrix exponential of a skew-symmetric matrix (rotation matrix).
///
/// Computes `exp(S)` where `S` is the skew-symmetric matrix corresponding
/// to angular velocity `omega` via Rodrigues' formula:
///
/// `exp(S) = I + sin(θ)/θ * S + (1 - cos(θ))/θ² * S²`
///
/// where `θ = |omega|`.  For `θ ≈ 0` returns the identity.
pub fn mat3_exp_skew(omega: &Vec3) -> Mat3 {
    let theta = omega.norm();
    if theta < 1e-12 {
        return Mat3::identity();
    }
    let axis = omega / theta;
    let s = skew_symmetric(&axis);
    let s2 = s * s;
    Mat3::identity() + s * theta.sin() + s2 * (1.0 - theta.cos())
}
/// Matrix logarithm of a rotation matrix `R`.
///
/// Returns the skew-symmetric matrix `S` such that `exp(S) = R`.
/// The returned matrix encodes the axis-angle `theta * n̂` in its entries.
/// Returns the zero matrix for the identity rotation.
pub fn mat3_log_rotation(r: &Mat3) -> Mat3 {
    let trace = r.trace();
    let cos_theta = ((trace - 1.0) / 2.0).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();
    if theta.abs() < 1e-12 {
        return Mat3::zeros();
    }
    (r - r.transpose()) * (theta / (2.0 * theta.sin()))
}
/// Normalized linear interpolation (nlerp) between two unit quaternions.
///
/// Faster but less accurate than slerp; the result is normalized to stay
/// on the unit sphere.
pub fn quat_nlerp(a: &Quat, b: &Quat, t: Real) -> Quat {
    let ai = a.into_inner();
    let bi = b.into_inner();
    let bi = if ai.dot(&bi) < 0.0 { -bi } else { bi };
    let lerped = ai * (1.0 - t) + bi * t;
    UnitQuaternion::new_normalize(lerped)
}
/// Build the inner control quaternion `s_i` needed for SQUAD interpolation.
///
/// Given three consecutive key quaternions `q_prev`, `q_curr`, `q_next`,
/// returns `s_i = q_curr * exp( -(log(q_curr^{-1} q_next) + log(q_curr^{-1} q_prev)) / 4 )`.
pub fn quat_squad_control(q_prev: &Quat, q_curr: &Quat, q_next: &Quat) -> Quat {
    let qi_inv = q_curr.inverse();
    let log_next = quat_log(&(qi_inv * q_next));
    let log_prev = quat_log(&(qi_inv * q_prev));
    let sum = log_next + log_prev;
    let half = sum * (-0.25);
    q_curr * quat_exp(&half)
}
/// Geodesic (angular) distance between two unit quaternions.
///
/// Returns the minimal angle in `[0, π]` needed to rotate from `a` to `b`.
pub fn quat_geodesic_distance(a: &Quat, b: &Quat) -> Real {
    a.angle_to(b)
}
/// Project `v` onto the unit direction `onto_unit`.
///
/// Returns the component of `v` parallel to `onto_unit`.
/// `onto_unit` must be a unit vector.
pub fn vec3_project_onto(v: &Vec3, onto_unit: &Vec3) -> Vec3 {
    onto_unit * v.dot(onto_unit)
}
/// Reject `v` from the unit direction `onto_unit` (perpendicular component).
///
/// Returns `v - project_onto(v, onto_unit)`.
/// `onto_unit` must be a unit vector.
pub fn vec3_reject_from(v: &Vec3, onto_unit: &Vec3) -> Vec3 {
    v - vec3_project_onto(v, onto_unit)
}
/// Reflect `v` about the unit normal `n`: `v - 2*(v·n)*n`.
pub fn vec3_reflect_about(v: &Vec3, n: &Vec3) -> Vec3 {
    v - n * (2.0 * v.dot(n))
}
/// Refract `v` about the unit normal `n` with relative index of refraction `eta`.
///
/// `v` must be normalized, `n` must be a unit normal pointing away from the
/// surface on the same side as `v`.  Returns `None` on total internal
/// reflection.
pub fn vec3_refract(v: &Vec3, n: &Vec3, eta: Real) -> Option<Vec3> {
    let cos_i = -(v.dot(n));
    let sin2_t = eta * eta * (1.0 - cos_i * cos_i);
    if sin2_t > 1.0 {
        return None;
    }
    let cos_t = (1.0 - sin2_t).sqrt();
    Some(v * eta + n * (eta * cos_i - cos_t))
}
/// Angle in radians between two non-zero vectors `a` and `b`.
///
/// Returns a value in `[0, π]`.
pub fn vec3_angle_between(a: &Vec3, b: &Vec3) -> Real {
    let denom = a.norm() * b.norm();
    if denom < 1e-300 {
        return 0.0;
    }
    (a.dot(b) / denom).clamp(-1.0, 1.0).acos()
}
/// Rotate `v` by quaternion `q`.
pub fn vec3_rotate_by_quat(v: &Vec3, q: &Quat) -> Vec3 {
    q.transform_vector(v)
}
/// Build a rotation matrix from an `axis` (unit vector) and `angle` (radians).
///
/// Uses the Rodrigues rotation formula.
pub fn mat3_from_axis_angle(axis: &Vec3, angle: Real) -> Mat3 {
    let c = angle.cos();
    let s = angle.sin();
    let t = 1.0 - c;
    let [x, y, z] = [axis.x, axis.y, axis.z];
    Mat3::new(
        t * x * x + c,
        t * x * y - s * z,
        t * x * z + s * y,
        t * x * y + s * z,
        t * y * y + c,
        t * y * z - s * x,
        t * x * z - s * y,
        t * y * z + s * x,
        t * z * z + c,
    )
}
