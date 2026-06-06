//! Internal helper math functions for medical visualization.

#[inline]
pub(crate) fn lerp(a: f32, b: f32, t: f32) -> f32 {
    a + (b - a) * t
}

#[inline]
pub(crate) fn clamp_f32(x: f32, lo: f32, hi: f32) -> f32 {
    x.max(lo).min(hi)
}

#[inline]
pub(crate) fn norm3f(v: [f32; 3]) -> f32 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

#[inline]
pub(crate) fn normalize3f(v: [f32; 3]) -> [f32; 3] {
    let n = norm3f(v);
    if n < 1e-10 {
        [0.0; 3]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

#[inline]
pub(crate) fn dot3f(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
