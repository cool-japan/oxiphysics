//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

/// FFT-based spectral differentiation on a periodic 1D domain.
///
/// Implements differentiation by multiplying Fourier coefficients
/// by the wavenumber ik.  Uses O(N²) DFT internally (no external FFT crate).
pub struct FftDiff1D {
    /// Number of grid points (ideally a power of 2).
    pub n: usize,
    /// Grid spacing dx = L / n.
    pub dx: f64,
}
impl FftDiff1D {
    /// Creates a new spectral differentiator.
    pub fn new(n: usize, dx: f64) -> Self {
        Self { n, dx }
    }
    /// Computes the DFT of `x` (complex output as flat `(re, im)` pairs).
    fn dft(&self, x: &[f64]) -> Vec<(f64, f64)> {
        let n = self.n;
        let mut out = vec![(0.0_f64, 0.0_f64); n];
        for (k, ok) in out.iter_mut().enumerate() {
            let mut re = 0.0;
            let mut im = 0.0;
            for (j, &xj) in x.iter().enumerate() {
                let angle = -2.0 * PI * (k as f64) * (j as f64) / n as f64;
                re += xj * angle.cos();
                im += xj * angle.sin();
            }
            *ok = (re, im);
        }
        out
    }
    /// Inverse DFT; takes `(re, im)` pairs and returns real part.
    fn idft(&self, xk: &[(f64, f64)]) -> Vec<f64> {
        let n = self.n;
        let mut out = vec![0.0_f64; n];
        for (j, oj) in out.iter_mut().enumerate() {
            let mut re = 0.0;
            for (k, &xkk) in xk.iter().enumerate() {
                let angle = 2.0 * PI * (k as f64) * (j as f64) / n as f64;
                re += xkk.0 * angle.cos() - xkk.1 * angle.sin();
            }
            *oj = re / n as f64;
        }
        out
    }
    /// Computes the first derivative du/dx via spectral multiplication by ik.
    ///
    /// Assumes periodic boundary conditions on \[0, n*dx\].
    pub fn differentiate(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let l = n as f64 * self.dx;
        let uk = self.dft(u);
        let mut duk: Vec<(f64, f64)> = (0..n)
            .map(|k| {
                let kk = if k <= n / 2 {
                    k as f64
                } else {
                    k as f64 - n as f64
                };
                let omega = 2.0 * PI * kk / l;
                (-uk[k].1 * omega, uk[k].0 * omega)
            })
            .collect();
        if n.is_multiple_of(2) {
            duk[n / 2] = (0.0, 0.0);
        }
        self.idft(&duk)
    }
    /// Computes the second derivative d²u/dx² via spectral multiplication by -k².
    pub fn differentiate2(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n;
        let l = n as f64 * self.dx;
        let uk = self.dft(u);
        let d2uk: Vec<(f64, f64)> = (0..n)
            .map(|k| {
                let kk = if k <= n / 2 {
                    k as f64
                } else {
                    k as f64 - n as f64
                };
                let omega2 = -(2.0 * PI * kk / l).powi(2);
                (uk[k].0 * omega2, uk[k].1 * omega2)
            })
            .collect();
        self.idft(&d2uk)
    }
}
/// 3D finite-difference operator set.
///
/// All fields are flat z-then-y-then-x (i + j*nx + k*nx*ny) vectors.
pub struct FiniteDiffOps3D {
    /// Grid points in x.
    pub nx: usize,
    /// Grid points in y.
    pub ny: usize,
    /// Grid points in z.
    pub nz: usize,
    /// Uniform grid spacing.
    pub dx: f64,
}
impl FiniteDiffOps3D {
    /// Creates a new 3D operator set (uniform spacing `dx`).
    pub fn new(nx: usize, ny: usize, nz: usize, dx: f64) -> Self {
        Self { nx, ny, nz, dx }
    }
    fn idx3(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.ny * self.nx + j * self.nx + i
    }
    /// 7-point 3D Laplacian ∇²u.
    ///
    /// Returns zero at all boundary planes.
    pub fn laplacian(&self, u: &[f64]) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let dx2 = self.dx * self.dx;
        let mut lap = vec![0.0_f64; nx * ny * nz];
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let c = self.idx3(i, j, k);
                    lap[c] = (u[c - 1] - 2.0 * u[c] + u[c + 1]) / dx2
                        + (u[c - nx] - 2.0 * u[c] + u[c + nx]) / dx2
                        + (u[c - nx * ny] - 2.0 * u[c] + u[c + nx * ny]) / dx2;
                }
            }
        }
        lap
    }
    /// Gradient of a scalar field in 3D.
    ///
    /// Returns (∂u/∂x, ∂u/∂y, ∂u/∂z) each of size nx*ny*nz.
    pub fn gradient(&self, u: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let two_dx = 2.0 * self.dx;
        let mut gx = vec![0.0_f64; nx * ny * nz];
        let mut gy = vec![0.0_f64; nx * ny * nz];
        let mut gz = vec![0.0_f64; nx * ny * nz];
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let c = self.idx3(i, j, k);
                    gx[c] = (u[c + 1] - u[c - 1]) / two_dx;
                    gy[c] = (u[c + nx] - u[c - nx]) / two_dx;
                    gz[c] = (u[c + nx * ny] - u[c - nx * ny]) / two_dx;
                }
            }
        }
        (gx, gy, gz)
    }
    /// Curl of a 3D vector field (fx, fy, fz).
    ///
    /// Returns the three components (curl_x, curl_y, curl_z).
    pub fn curl(&self, fx: &[f64], fy: &[f64], fz: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let two_dx = 2.0 * self.dx;
        let nn = nx * ny * nz;
        let mut cx = vec![0.0_f64; nn];
        let mut cy = vec![0.0_f64; nn];
        let mut cz = vec![0.0_f64; nn];
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let c = self.idx3(i, j, k);
                    cx[c] = (fz[c + nx] - fz[c - nx]) / two_dx
                        - (fy[c + nx * ny] - fy[c - nx * ny]) / two_dx;
                    cy[c] = (fx[c + nx * ny] - fx[c - nx * ny]) / two_dx
                        - (fz[c + 1] - fz[c - 1]) / two_dx;
                    cz[c] = (fy[c + 1] - fy[c - 1]) / two_dx - (fx[c + nx] - fx[c - nx]) / two_dx;
                }
            }
        }
        (cx, cy, cz)
    }
}
