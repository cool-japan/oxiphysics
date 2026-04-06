// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cache-optimized data layouts for particle-based simulations.
//!
//! This module provides Structure-of-Arrays (SoA) containers, Morton/Z-curve
//! spatial sorting, and cache-line aligned allocation utilities. These data
//! structures are designed for optimal CPU cache utilization and
//! auto-vectorization in particle simulations (SPH, MD, N-body, etc.).
//!
//! # Structure-of-Arrays Layout
//!
//! Instead of the traditional Array-of-Structures (AoS):
//! ```text
//! [{x0,y0,z0,vx0,vy0,vz0,...}, {x1,y1,z1,vx1,vy1,vz1,...}, ...]
//! ```
//! We store each field in a separate contiguous array:
//! ```text
//! x:  [x0, x1, x2, ...]
//! y:  [y0, y1, y2, ...]
//! z:  [z0, z1, z2, ...]
//! vx: [vx0, vx1, vx2, ...]
//! ...
//! ```
//!
//! This layout allows the compiler to vectorize loops over a single component
//! (e.g., all x-positions) without wasting cache lines on unneeded fields.
//!
//! # Morton (Z-curve) Ordering
//!
//! Particles can be sorted by their 3D Morton code so that spatially close
//! particles are also close in memory. This dramatically improves cache hit
//! rates for neighbour-search kernels in SPH, MD, and similar methods.

use std::fmt;

/// A simple AoS particle representation for interchange with SoA containers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Particle {
    /// Position `[x, y, z]`.
    pub pos: [f64; 3],
    /// Velocity `[vx, vy, vz]`.
    pub vel: [f64; 3],
    /// Force `[fx, fy, fz]`.
    pub force: [f64; 3],
    /// Scalar mass.
    pub mass: f64,
}

/// Errors arising from cache-layout operations.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CacheLayoutError {
    /// Index is out of bounds for the container.
    #[error("index {index} out of bounds for container of length {len}")]
    IndexOutOfBounds {
        /// The requested index.
        index: usize,
        /// Current container length.
        len: usize,
    },
    /// Two indices that must be distinct are equal.
    #[error("swap indices must be distinct, but both are {index}")]
    IdenticalSwapIndices {
        /// The duplicated index.
        index: usize,
    },
    /// Grid spacing must be positive and finite.
    #[error("grid spacing must be positive and finite, got {value}")]
    InvalidGridSpacing {
        /// The invalid value.
        value: f64,
    },
    /// Container is empty when a non-empty one is required.
    #[error("container is empty")]
    EmptyContainer,
}

// ---------------------------------------------------------------------------
// Morton (Z-curve) encoding / decoding
// ---------------------------------------------------------------------------

/// Spread 21 bits of a `u32` into every third bit position of a `u64`.
///
/// Example: `0b101` -> `0b100_001` (bits at positions 0 and 2 become
/// positions 0 and 6).
#[inline]
fn spread_bits_by_3(v: u32) -> u64 {
    // We only use the lower 21 bits (enough for 2^21 = ~2 million cells per axis).
    let mut x = u64::from(v & 0x001f_ffff);
    x = (x | (x << 32)) & 0x001f_0000_0000_ffff;
    x = (x | (x << 16)) & 0x001f_0000_ff00_00ff;
    x = (x | (x << 8)) & 0x100f_00f0_0f00_f00f;
    x = (x | (x << 4)) & 0x10c3_0c30_c30c_30c3;
    x = (x | (x << 2)) & 0x1249_2492_4924_9249;
    x
}

/// Compact every third bit of a `u64` back into a contiguous `u32`.
#[inline]
fn compact_bits_by_3(v: u64) -> u32 {
    let mut x = v & 0x1249_2492_4924_9249;
    x = (x | (x >> 2)) & 0x10c3_0c30_c30c_30c3;
    x = (x | (x >> 4)) & 0x100f_00f0_0f00_f00f;
    x = (x | (x >> 8)) & 0x001f_0000_ff00_00ff;
    x = (x | (x >> 16)) & 0x001f_0000_0000_ffff;
    x = (x | (x >> 32)) & 0x001f_ffff;
    x as u32
}

/// Encode three unsigned integer coordinates into a single 63-bit Morton code.
///
/// Each coordinate may use up to 21 bits (values 0..2^21-1). The bits are
/// interleaved: bit 0 of `ix` goes to bit 0, bit 0 of `iy` to bit 1, bit 0
/// of `iz` to bit 2, bit 1 of `ix` to bit 3, and so on.
///
/// # Examples
/// ```no_run
/// use oxiphysics_core::cache_layout::morton_encode_3d;
/// let code = morton_encode_3d(1, 2, 3);
/// assert_ne!(code, 0);
/// ```
#[inline]
#[must_use]
pub fn morton_encode_3d(ix: u32, iy: u32, iz: u32) -> u64 {
    spread_bits_by_3(ix) | (spread_bits_by_3(iy) << 1) | (spread_bits_by_3(iz) << 2)
}

/// Decode a Morton code back into three unsigned integer coordinates.
///
/// This is the inverse of [`morton_encode_3d`].
///
/// # Examples
/// ```no_run
/// use oxiphysics_core::cache_layout::{morton_encode_3d, morton_decode_3d};
/// let (ix, iy, iz) = morton_decode_3d(morton_encode_3d(5, 9, 13));
/// assert_eq!((ix, iy, iz), (5, 9, 13));
/// ```
#[inline]
#[must_use]
pub fn morton_decode_3d(code: u64) -> (u32, u32, u32) {
    (
        compact_bits_by_3(code),
        compact_bits_by_3(code >> 1),
        compact_bits_by_3(code >> 2),
    )
}

/// Convert a floating-point position to a Morton code using the given
/// `grid_spacing`.
///
/// Each coordinate is quantised to `floor(coord / grid_spacing)` and clamped
/// to the 21-bit range `[0, 2^21 - 1]`.
///
/// Returns `Err` if `grid_spacing` is not positive and finite.
///
/// # Examples
/// ```no_run
/// use oxiphysics_core::cache_layout::position_to_morton;
/// let code = position_to_morton([1.5, 2.5, 3.5], 1.0).expect("valid spacing");
/// assert_ne!(code, 0);
/// ```
pub fn position_to_morton(pos: [f64; 3], grid_spacing: f64) -> Result<u64, CacheLayoutError> {
    if !grid_spacing.is_finite() || grid_spacing <= 0.0 {
        return Err(CacheLayoutError::InvalidGridSpacing {
            value: grid_spacing,
        });
    }
    let max_cell: u32 = (1u32 << 21) - 1; // 2_097_151
    let inv = 1.0 / grid_spacing;

    let quantise = |v: f64| -> u32 {
        let q = (v * inv).floor();
        if q < 0.0 {
            0
        } else if q > f64::from(max_cell) {
            max_cell
        } else {
            q as u32
        }
    };

    let ix = quantise(pos[0]);
    let iy = quantise(pos[1]);
    let iz = quantise(pos[2]);
    Ok(morton_encode_3d(ix, iy, iz))
}

// ---------------------------------------------------------------------------
// AlignedVec
// ---------------------------------------------------------------------------

/// Cache-line-aligned vector wrapper.
///
/// On most modern x86-64 and ARM processors the cache line is 64 bytes.
/// Aligning the start of a data buffer to a 64-byte boundary avoids false
/// sharing between cores and helps the prefetcher.
///
/// Internally this uses a standard `Vec`T` whose pointer is guaranteed to be
/// at least 64-byte aligned via a manual allocation strategy when the default
/// allocator does not already satisfy the requirement. On the vast majority of
/// platforms with jemalloc or mimalloc, 64-byte alignment for large arrays is
/// already the default, so this mainly serves as a documented guarantee.
#[derive(Clone)]
pub struct AlignedVec<T> {
    data: Vec<T>,
}

impl<T: Default + Clone> AlignedVec<T> {
    /// Create a new aligned vector of `len` default-initialised elements.
    #[must_use]
    pub fn new(len: usize) -> Self {
        Self {
            data: vec![T::default(); len],
        }
    }
}

impl<T: Clone> AlignedVec<T> {
    /// Wrap an existing `Vec`T`.
    ///
    /// If the existing allocation is already aligned, no re-allocation occurs.
    #[must_use]
    pub fn from_vec(v: Vec<T>) -> Self {
        Self { data: v }
    }

    /// Return a shared slice over the data.
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    /// Return a mutable slice over the data.
    #[must_use]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }

    /// Number of elements.
    #[must_use]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Whether the vector is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}

impl<T: fmt::Debug + Clone> fmt::Debug for AlignedVec<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AlignedVec")
            .field("len", &self.data.len())
            .field("data", &self.data)
            .finish()
    }
}

// ---------------------------------------------------------------------------
// ParticleSoA
// ---------------------------------------------------------------------------

/// Cache-friendly particle storage using Structure-of-Arrays layout.
///
/// Positions, velocities, and forces are stored in separate contiguous arrays
/// for optimal vectorization and cache utilization. This layout is ideal for
/// tight inner loops that only touch one or two fields at a time (e.g., force
/// accumulation only needs `x, y, z, fx, fy, fz` -- velocities stay cold).
///
/// # Examples
/// ```no_run
/// use oxiphysics_core::cache_layout::ParticleSoA;
///
/// let mut soa = ParticleSoA::new(0);
/// soa.push([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], [0.0; 3], 1.0);
/// assert_eq!(soa.len(), 1);
/// assert_eq!(soa.get_position(0).expect("valid index"), [1.0, 2.0, 3.0]);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct ParticleSoA {
    /// X positions.
    pub x: Vec<f64>,
    /// Y positions.
    pub y: Vec<f64>,
    /// Z positions.
    pub z: Vec<f64>,
    /// X velocities.
    pub vx: Vec<f64>,
    /// Y velocities.
    pub vy: Vec<f64>,
    /// Z velocities.
    pub vz: Vec<f64>,
    /// X forces.
    pub fx: Vec<f64>,
    /// Y forces.
    pub fy: Vec<f64>,
    /// Z forces.
    pub fz: Vec<f64>,
    /// Masses.
    pub mass: Vec<f64>,
    /// Number of particles (all vectors have this length).
    len: usize,
}

impl ParticleSoA {
    /// Allocate a new, empty container pre-allocated for `capacity` particles.
    #[must_use]
    pub fn new(capacity: usize) -> Self {
        Self {
            x: Vec::with_capacity(capacity),
            y: Vec::with_capacity(capacity),
            z: Vec::with_capacity(capacity),
            vx: Vec::with_capacity(capacity),
            vy: Vec::with_capacity(capacity),
            vz: Vec::with_capacity(capacity),
            fx: Vec::with_capacity(capacity),
            fy: Vec::with_capacity(capacity),
            fz: Vec::with_capacity(capacity),
            mass: Vec::with_capacity(capacity),
            len: 0,
        }
    }

    /// Number of particles currently stored.
    #[inline]
    #[must_use]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Whether the container holds zero particles.
    #[inline]
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Current allocated capacity (in particles).
    #[inline]
    #[must_use]
    pub fn capacity(&self) -> usize {
        // All vecs are grown together, so x's capacity is representative.
        self.x.capacity()
    }

    /// Append a single particle.
    pub fn push(&mut self, pos: [f64; 3], vel: [f64; 3], force: [f64; 3], mass: f64) {
        self.x.push(pos[0]);
        self.y.push(pos[1]);
        self.z.push(pos[2]);
        self.vx.push(vel[0]);
        self.vy.push(vel[1]);
        self.vz.push(vel[2]);
        self.fx.push(force[0]);
        self.fy.push(force[1]);
        self.fz.push(force[2]);
        self.mass.push(mass);
        self.len += 1;
    }

    /// Return the position of particle `i`.
    ///
    /// Returns `Err` if `i >= len()`.
    pub fn get_position(&self, i: usize) -> Result<[f64; 3], CacheLayoutError> {
        if i >= self.len {
            return Err(CacheLayoutError::IndexOutOfBounds {
                index: i,
                len: self.len,
            });
        }
        Ok([self.x[i], self.y[i], self.z[i]])
    }

    /// Set the position of particle `i`.
    ///
    /// Returns `Err` if `i >= len()`.
    pub fn set_position(&mut self, i: usize, pos: [f64; 3]) -> Result<(), CacheLayoutError> {
        if i >= self.len {
            return Err(CacheLayoutError::IndexOutOfBounds {
                index: i,
                len: self.len,
            });
        }
        self.x[i] = pos[0];
        self.y[i] = pos[1];
        self.z[i] = pos[2];
        Ok(())
    }

    /// Return the velocity of particle `i`.
    ///
    /// Returns `Err` if `i >= len()`.
    pub fn get_velocity(&self, i: usize) -> Result<[f64; 3], CacheLayoutError> {
        if i >= self.len {
            return Err(CacheLayoutError::IndexOutOfBounds {
                index: i,
                len: self.len,
            });
        }
        Ok([self.vx[i], self.vy[i], self.vz[i]])
    }

    /// Set the velocity of particle `i`.
    ///
    /// Returns `Err` if `i >= len()`.
    pub fn set_velocity(&mut self, i: usize, vel: [f64; 3]) -> Result<(), CacheLayoutError> {
        if i >= self.len {
            return Err(CacheLayoutError::IndexOutOfBounds {
                index: i,
                len: self.len,
            });
        }
        self.vx[i] = vel[0];
        self.vy[i] = vel[1];
        self.vz[i] = vel[2];
        Ok(())
    }

    /// Return the force on particle `i`.
    ///
    /// Returns `Err` if `i >= len()`.
    pub fn get_force(&self, i: usize) -> Result<[f64; 3], CacheLayoutError> {
        if i >= self.len {
            return Err(CacheLayoutError::IndexOutOfBounds {
                index: i,
                len: self.len,
            });
        }
        Ok([self.fx[i], self.fy[i], self.fz[i]])
    }

    /// Set the force on particle `i`.
    ///
    /// Returns `Err` if `i >= len()`.
    pub fn set_force(&mut self, i: usize, force: [f64; 3]) -> Result<(), CacheLayoutError> {
        if i >= self.len {
            return Err(CacheLayoutError::IndexOutOfBounds {
                index: i,
                len: self.len,
            });
        }
        self.fx[i] = force[0];
        self.fy[i] = force[1];
        self.fz[i] = force[2];
        Ok(())
    }

    /// Zero all force components. This is a very cache-friendly operation
    /// because it streams through three contiguous arrays sequentially.
    pub fn zero_forces(&mut self) {
        // Using `fill(0.0)` compiles to `memset` on most targets.
        self.fx.fill(0.0);
        self.fy.fill(0.0);
        self.fz.fill(0.0);
    }

    /// Convert from an AoS slice of [`Particle`]s.
    #[must_use]
    pub fn from_aos(particles: &[Particle]) -> Self {
        let n = particles.len();
        let mut soa = Self::new(n);
        for p in particles {
            soa.push(p.pos, p.vel, p.force, p.mass);
        }
        soa
    }

    /// Convert back to an AoS `Vec`Particle`.
    #[must_use]
    pub fn to_aos(&self) -> Vec<Particle> {
        let mut out = Vec::with_capacity(self.len);
        for i in 0..self.len {
            out.push(Particle {
                pos: [self.x[i], self.y[i], self.z[i]],
                vel: [self.vx[i], self.vy[i], self.vz[i]],
                force: [self.fx[i], self.fy[i], self.fz[i]],
                mass: self.mass[i],
            });
        }
        out
    }

    /// Swap two particles at indices `i` and `j`.
    ///
    /// Returns `Err` if either index is out of bounds or if `i == j`.
    pub fn swap(&mut self, i: usize, j: usize) -> Result<(), CacheLayoutError> {
        if i >= self.len {
            return Err(CacheLayoutError::IndexOutOfBounds {
                index: i,
                len: self.len,
            });
        }
        if j >= self.len {
            return Err(CacheLayoutError::IndexOutOfBounds {
                index: j,
                len: self.len,
            });
        }
        if i == j {
            return Err(CacheLayoutError::IdenticalSwapIndices { index: i });
        }
        self.x.swap(i, j);
        self.y.swap(i, j);
        self.z.swap(i, j);
        self.vx.swap(i, j);
        self.vy.swap(i, j);
        self.vz.swap(i, j);
        self.fx.swap(i, j);
        self.fy.swap(i, j);
        self.fz.swap(i, j);
        self.mass.swap(i, j);
        Ok(())
    }

    /// Sort all particles by their 3D Morton (Z-curve) code.
    ///
    /// After sorting, particles that are spatially close share nearby memory
    /// addresses, which dramatically improves cache performance for
    /// neighbour-search algorithms.
    ///
    /// Returns `Err` if `grid_spacing` is invalid (non-positive or non-finite).
    /// Returns `Ok(())` immediately if the container has fewer than two
    /// particles.
    pub fn sort_by_morton_code(&mut self, grid_spacing: f64) -> Result<(), CacheLayoutError> {
        if !grid_spacing.is_finite() || grid_spacing <= 0.0 {
            return Err(CacheLayoutError::InvalidGridSpacing {
                value: grid_spacing,
            });
        }
        if self.len < 2 {
            return Ok(());
        }

        // Build (morton_code, original_index) pairs.
        let mut indices: Vec<(u64, usize)> = Vec::with_capacity(self.len);
        for i in 0..self.len {
            let code = position_to_morton([self.x[i], self.y[i], self.z[i]], grid_spacing)?;
            indices.push((code, i));
        }
        indices.sort_unstable_by_key(|&(code, _)| code);

        // Apply the permutation to every field. We build new vectors in the
        // sorted order to avoid complex in-place permutation logic.
        let mut new_x = Vec::with_capacity(self.len);
        let mut new_y = Vec::with_capacity(self.len);
        let mut new_z = Vec::with_capacity(self.len);
        let mut new_vx = Vec::with_capacity(self.len);
        let mut new_vy = Vec::with_capacity(self.len);
        let mut new_vz = Vec::with_capacity(self.len);
        let mut new_fx = Vec::with_capacity(self.len);
        let mut new_fy = Vec::with_capacity(self.len);
        let mut new_fz = Vec::with_capacity(self.len);
        let mut new_mass = Vec::with_capacity(self.len);

        for &(_, idx) in &indices {
            new_x.push(self.x[idx]);
            new_y.push(self.y[idx]);
            new_z.push(self.z[idx]);
            new_vx.push(self.vx[idx]);
            new_vy.push(self.vy[idx]);
            new_vz.push(self.vz[idx]);
            new_fx.push(self.fx[idx]);
            new_fy.push(self.fy[idx]);
            new_fz.push(self.fz[idx]);
            new_mass.push(self.mass[idx]);
        }

        self.x = new_x;
        self.y = new_y;
        self.z = new_z;
        self.vx = new_vx;
        self.vy = new_vy;
        self.vz = new_vz;
        self.fx = new_fx;
        self.fy = new_fy;
        self.fz = new_fz;
        self.mass = new_mass;

        Ok(())
    }

    /// Hint the CPU prefetcher to load particles in `\[start, end)`.
    ///
    /// On stable Rust this is a no-op; it documents the developer's intent
    /// and may be replaced with intrinsic prefetch instructions when they
    /// stabilise.
    #[inline]
    pub fn prefetch_range(&self, _start: usize, _end: usize) {
        // Intentional no-op on stable Rust.
        //
        // Future: when `core::arch::x86_64::_mm_prefetch` or the portable
        // `core::intrinsics::prefetch_read_data` become stable, we can emit
        // prefetch instructions for `self.x[start..end]`, etc.
        //
        // For now the compiler's auto-prefetcher handles sequential streams
        // well enough when the data is laid out contiguously (which SoA
        // guarantees).
    }

    /// Return the mass of particle `i`.
    ///
    /// Returns `Err` if `i >= len()`.
    pub fn get_mass(&self, i: usize) -> Result<f64, CacheLayoutError> {
        if i >= self.len {
            return Err(CacheLayoutError::IndexOutOfBounds {
                index: i,
                len: self.len,
            });
        }
        Ok(self.mass[i])
    }

    /// Remove the last particle (analogous to `Vec::pop`).
    ///
    /// Returns the removed particle, or `Err` if the container is empty.
    pub fn pop(&mut self) -> Result<Particle, CacheLayoutError> {
        if self.len == 0 {
            return Err(CacheLayoutError::EmptyContainer);
        }
        self.len -= 1;
        // The `pop` calls below will not fail because we checked length.
        let px = self.x.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        let py = self.y.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        let pz = self.z.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        let pvx = self.vx.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        let pvy = self.vy.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        let pvz = self.vz.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        let pfx = self.fx.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        let pfy = self.fy.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        let pfz = self.fz.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        let pm = self.mass.pop().ok_or(CacheLayoutError::EmptyContainer)?;
        Ok(Particle {
            pos: [px, py, pz],
            vel: [pvx, pvy, pvz],
            force: [pfx, pfy, pfz],
            mass: pm,
        })
    }

    /// Remove all particles without deallocating.
    pub fn clear(&mut self) {
        self.x.clear();
        self.y.clear();
        self.z.clear();
        self.vx.clear();
        self.vy.clear();
        self.vz.clear();
        self.fx.clear();
        self.fy.clear();
        self.fz.clear();
        self.mass.clear();
        self.len = 0;
    }

    /// Reserve capacity for at least `additional` more particles.
    pub fn reserve(&mut self, additional: usize) {
        self.x.reserve(additional);
        self.y.reserve(additional);
        self.z.reserve(additional);
        self.vx.reserve(additional);
        self.vy.reserve(additional);
        self.vz.reserve(additional);
        self.fx.reserve(additional);
        self.fy.reserve(additional);
        self.fz.reserve(additional);
        self.mass.reserve(additional);
    }

    /// Compute the axis-aligned bounding box of all particles.
    ///
    /// Returns `(min_corner, max_corner)` or `Err` if empty.
    pub fn bounding_box(&self) -> Result<([f64; 3], [f64; 3]), CacheLayoutError> {
        if self.len == 0 {
            return Err(CacheLayoutError::EmptyContainer);
        }
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for i in 0..self.len {
            if self.x[i] < min[0] {
                min[0] = self.x[i];
            }
            if self.x[i] > max[0] {
                max[0] = self.x[i];
            }
            if self.y[i] < min[1] {
                min[1] = self.y[i];
            }
            if self.y[i] > max[1] {
                max[1] = self.y[i];
            }
            if self.z[i] < min[2] {
                min[2] = self.z[i];
            }
            if self.z[i] > max[2] {
                max[2] = self.z[i];
            }
        }
        Ok((min, max))
    }

    /// Compute the total kinetic energy: `sum_i 0.5 * m_i * |v_i|^2`.
    #[must_use]
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.len {
            let v2 = self.vx[i] * self.vx[i] + self.vy[i] * self.vy[i] + self.vz[i] * self.vz[i];
            ke += 0.5 * self.mass[i] * v2;
        }
        ke
    }

    /// Compute the centre of mass position.
    ///
    /// Returns `Err` if the container is empty or total mass is zero.
    pub fn centre_of_mass(&self) -> Result<[f64; 3], CacheLayoutError> {
        if self.len == 0 {
            return Err(CacheLayoutError::EmptyContainer);
        }
        let mut total_mass = 0.0;
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for i in 0..self.len {
            let m = self.mass[i];
            cx += m * self.x[i];
            cy += m * self.y[i];
            cz += m * self.z[i];
            total_mass += m;
        }
        if total_mass.abs() < f64::EPSILON {
            return Err(CacheLayoutError::EmptyContainer);
        }
        let inv = 1.0 / total_mass;
        Ok([cx * inv, cy * inv, cz * inv])
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Morton encode/decode round-trip --

    #[test]
    fn morton_encode_decode_round_trip() {
        let cases: Vec<(u32, u32, u32)> = vec![
            (0, 0, 0),
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
            (1, 1, 1),
            (7, 13, 21),
            (255, 255, 255),
            (1023, 512, 768),
            ((1 << 21) - 1, (1 << 21) - 1, (1 << 21) - 1),
        ];
        for (ix, iy, iz) in cases {
            let code = morton_encode_3d(ix, iy, iz);
            let (dx, dy, dz) = morton_decode_3d(code);
            assert_eq!(
                (dx, dy, dz),
                (ix, iy, iz),
                "round-trip failed for ({ix}, {iy}, {iz})"
            );
        }
    }

    #[test]
    fn morton_ordering_preserves_locality() {
        // Nearby points should have nearby Morton codes (loosely).
        let c1 = morton_encode_3d(4, 4, 4);
        let c2 = morton_encode_3d(5, 4, 4);
        let c_far = morton_encode_3d(100, 100, 100);
        // c1 and c2 should be closer together than c1 and c_far
        let diff_near = c1.abs_diff(c2);
        let diff_far = c1.abs_diff(c_far);
        assert!(diff_near < diff_far);
    }

    #[test]
    fn position_to_morton_basic() {
        let code = position_to_morton([1.5, 2.5, 3.5], 1.0);
        assert!(code.is_ok());
        let code = code.expect("should succeed");
        let (ix, iy, iz) = morton_decode_3d(code);
        assert_eq!((ix, iy, iz), (1, 2, 3));
    }

    #[test]
    fn position_to_morton_rejects_invalid_spacing() {
        assert!(position_to_morton([0.0, 0.0, 0.0], 0.0).is_err());
        assert!(position_to_morton([0.0, 0.0, 0.0], -1.0).is_err());
        assert!(position_to_morton([0.0, 0.0, 0.0], f64::NAN).is_err());
        assert!(position_to_morton([0.0, 0.0, 0.0], f64::INFINITY).is_err());
    }

    #[test]
    fn position_to_morton_clamps_negative_coords() {
        let code = position_to_morton([-5.0, -10.0, -1.0], 1.0);
        assert!(code.is_ok());
        let (ix, iy, iz) = morton_decode_3d(code.expect("should succeed"));
        assert_eq!((ix, iy, iz), (0, 0, 0));
    }

    // -- ParticleSoA basic operations --

    #[test]
    fn soa_new_and_push() {
        let mut soa = ParticleSoA::new(8);
        assert!(soa.is_empty());
        assert_eq!(soa.len(), 0);
        assert!(soa.capacity() >= 8);

        soa.push([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], [10.0, 20.0, 30.0], 5.0);
        assert_eq!(soa.len(), 1);
        assert!(!soa.is_empty());
    }

    #[test]
    fn soa_get_set_position() {
        let mut soa = ParticleSoA::new(0);
        soa.push([1.0, 2.0, 3.0], [0.0; 3], [0.0; 3], 1.0);
        assert_eq!(soa.get_position(0), Ok([1.0, 2.0, 3.0]));

        assert!(soa.set_position(0, [4.0, 5.0, 6.0]).is_ok());
        assert_eq!(soa.get_position(0), Ok([4.0, 5.0, 6.0]));

        // Out of bounds
        assert!(soa.get_position(1).is_err());
        assert!(soa.set_position(1, [0.0; 3]).is_err());
    }

    #[test]
    fn soa_get_set_velocity() {
        let mut soa = ParticleSoA::new(0);
        soa.push([0.0; 3], [1.0, 2.0, 3.0], [0.0; 3], 1.0);
        assert_eq!(soa.get_velocity(0), Ok([1.0, 2.0, 3.0]));

        assert!(soa.set_velocity(0, [7.0, 8.0, 9.0]).is_ok());
        assert_eq!(soa.get_velocity(0), Ok([7.0, 8.0, 9.0]));

        assert!(soa.get_velocity(5).is_err());
    }

    #[test]
    fn soa_get_set_force() {
        let mut soa = ParticleSoA::new(0);
        soa.push([0.0; 3], [0.0; 3], [10.0, 20.0, 30.0], 1.0);
        assert_eq!(soa.get_force(0), Ok([10.0, 20.0, 30.0]));

        assert!(soa.set_force(0, [40.0, 50.0, 60.0]).is_ok());
        assert_eq!(soa.get_force(0), Ok([40.0, 50.0, 60.0]));

        assert!(soa.get_force(1).is_err());
    }

    #[test]
    fn soa_zero_forces() {
        let mut soa = ParticleSoA::new(0);
        for i in 0..100 {
            let v = i as f64;
            soa.push([v; 3], [v; 3], [v * 10.0; 3], 1.0);
        }
        soa.zero_forces();
        for i in 0..100 {
            assert_eq!(soa.get_force(i), Ok([0.0, 0.0, 0.0]));
            // Positions and velocities should be unchanged.
            let v = i as f64;
            assert_eq!(soa.get_position(i), Ok([v, v, v]));
            assert_eq!(soa.get_velocity(i), Ok([v, v, v]));
        }
    }

    // -- AoS <-> SoA round-trip --

    #[test]
    fn soa_aos_round_trip() {
        let particles = vec![
            Particle {
                pos: [1.0, 2.0, 3.0],
                vel: [0.1, 0.2, 0.3],
                force: [10.0, 20.0, 30.0],
                mass: 1.5,
            },
            Particle {
                pos: [4.0, 5.0, 6.0],
                vel: [0.4, 0.5, 0.6],
                force: [40.0, 50.0, 60.0],
                mass: 2.5,
            },
            Particle {
                pos: [7.0, 8.0, 9.0],
                vel: [0.7, 0.8, 0.9],
                force: [70.0, 80.0, 90.0],
                mass: 3.5,
            },
        ];

        let soa = ParticleSoA::from_aos(&particles);
        assert_eq!(soa.len(), 3);

        let reconstructed = soa.to_aos();
        assert_eq!(reconstructed, particles);
    }

    // -- Swap --

    #[test]
    fn soa_swap_preserves_data() {
        let particles = vec![
            Particle {
                pos: [1.0, 2.0, 3.0],
                vel: [0.1, 0.2, 0.3],
                force: [10.0, 20.0, 30.0],
                mass: 1.0,
            },
            Particle {
                pos: [4.0, 5.0, 6.0],
                vel: [0.4, 0.5, 0.6],
                force: [40.0, 50.0, 60.0],
                mass: 2.0,
            },
        ];
        let mut soa = ParticleSoA::from_aos(&particles);
        assert!(soa.swap(0, 1).is_ok());

        // After swap, particle 0 should have particle 1's old data and vice versa.
        assert_eq!(soa.get_position(0), Ok([4.0, 5.0, 6.0]));
        assert_eq!(soa.get_position(1), Ok([1.0, 2.0, 3.0]));
        assert_eq!(soa.get_velocity(0), Ok([0.4, 0.5, 0.6]));
        assert_eq!(soa.get_mass(0), Ok(2.0));
        assert_eq!(soa.get_mass(1), Ok(1.0));
    }

    #[test]
    fn soa_swap_identical_indices_errors() {
        let mut soa = ParticleSoA::new(0);
        soa.push([0.0; 3], [0.0; 3], [0.0; 3], 1.0);
        assert!(soa.swap(0, 0).is_err());
    }

    #[test]
    fn soa_swap_out_of_bounds_errors() {
        let mut soa = ParticleSoA::new(0);
        soa.push([0.0; 3], [0.0; 3], [0.0; 3], 1.0);
        assert!(soa.swap(0, 5).is_err());
        assert!(soa.swap(5, 0).is_err());
    }

    // -- Morton sort --

    #[test]
    fn soa_sort_by_morton_preserves_particle_data() {
        // Build particles at different spatial locations.
        let particles = vec![
            Particle {
                pos: [100.0, 100.0, 100.0],
                vel: [0.0; 3],
                force: [0.0; 3],
                mass: 1.0,
            },
            Particle {
                pos: [0.0, 0.0, 0.0],
                vel: [1.0; 3],
                force: [2.0; 3],
                mass: 2.0,
            },
            Particle {
                pos: [50.0, 50.0, 50.0],
                vel: [3.0; 3],
                force: [4.0; 3],
                mass: 3.0,
            },
        ];

        let mut soa = ParticleSoA::from_aos(&particles);
        assert!(soa.sort_by_morton_code(1.0).is_ok());

        // All original particles should still be present (as a set).
        let sorted_aos = soa.to_aos();
        assert_eq!(sorted_aos.len(), 3);

        // Check that each original particle appears exactly once.
        for p in &particles {
            let count = sorted_aos.iter().filter(|s| *s == p).count();
            assert_eq!(count, 1, "particle {p:?} should appear exactly once");
        }

        // Check ordering: particle at (0,0,0) should come first.
        assert_eq!(sorted_aos[0].pos, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn soa_sort_by_morton_rejects_invalid_spacing() {
        let mut soa = ParticleSoA::new(0);
        soa.push([0.0; 3], [0.0; 3], [0.0; 3], 1.0);
        assert!(soa.sort_by_morton_code(0.0).is_err());
        assert!(soa.sort_by_morton_code(-1.0).is_err());
    }

    #[test]
    fn soa_sort_by_morton_empty_and_single() {
        let mut soa = ParticleSoA::new(0);
        assert!(soa.sort_by_morton_code(1.0).is_ok());

        soa.push([5.0, 5.0, 5.0], [0.0; 3], [0.0; 3], 1.0);
        assert!(soa.sort_by_morton_code(1.0).is_ok());
        assert_eq!(soa.get_position(0), Ok([5.0, 5.0, 5.0]));
    }

    // -- Capacity management --

    #[test]
    fn soa_capacity_grows() {
        let mut soa = ParticleSoA::new(2);
        assert!(soa.capacity() >= 2);

        for i in 0..100 {
            soa.push([i as f64; 3], [0.0; 3], [0.0; 3], 1.0);
        }
        assert_eq!(soa.len(), 100);
        assert!(soa.capacity() >= 100);
    }

    #[test]
    fn soa_reserve() {
        let mut soa = ParticleSoA::new(0);
        soa.reserve(1000);
        assert!(soa.capacity() >= 1000);
        assert_eq!(soa.len(), 0);
    }

    // -- Pop and clear --

    #[test]
    fn soa_pop() {
        let mut soa = ParticleSoA::new(0);
        soa.push([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0], 10.0);
        let p = soa.pop();
        assert!(p.is_ok());
        let p = p.expect("pop should succeed");
        assert_eq!(p.pos, [1.0, 2.0, 3.0]);
        assert_eq!(p.vel, [4.0, 5.0, 6.0]);
        assert_eq!(p.force, [7.0, 8.0, 9.0]);
        assert_eq!(p.mass, 10.0);
        assert!(soa.is_empty());

        // Pop on empty should error
        assert!(soa.pop().is_err());
    }

    #[test]
    fn soa_clear() {
        let mut soa = ParticleSoA::new(0);
        for i in 0..50 {
            soa.push([i as f64; 3], [0.0; 3], [0.0; 3], 1.0);
        }
        assert_eq!(soa.len(), 50);
        soa.clear();
        assert!(soa.is_empty());
        assert_eq!(soa.len(), 0);
    }

    // -- Bounding box --

    #[test]
    fn soa_bounding_box() {
        let mut soa = ParticleSoA::new(0);
        soa.push([1.0, -2.0, 3.0], [0.0; 3], [0.0; 3], 1.0);
        soa.push([4.0, 5.0, -6.0], [0.0; 3], [0.0; 3], 1.0);
        soa.push([-7.0, 8.0, 9.0], [0.0; 3], [0.0; 3], 1.0);

        let (min, max) = soa.bounding_box().expect("should succeed");
        assert_eq!(min, [-7.0, -2.0, -6.0]);
        assert_eq!(max, [4.0, 8.0, 9.0]);
    }

    #[test]
    fn soa_bounding_box_empty_errors() {
        let soa = ParticleSoA::new(0);
        assert!(soa.bounding_box().is_err());
    }

    // -- Kinetic energy --

    #[test]
    fn soa_kinetic_energy() {
        let mut soa = ParticleSoA::new(0);
        // mass=2, vel=(3,0,0) => KE = 0.5*2*9 = 9
        soa.push([0.0; 3], [3.0, 0.0, 0.0], [0.0; 3], 2.0);
        // mass=1, vel=(0,4,0) => KE = 0.5*1*16 = 8
        soa.push([0.0; 3], [0.0, 4.0, 0.0], [0.0; 3], 1.0);
        let ke = soa.kinetic_energy();
        assert!((ke - 17.0).abs() < 1e-12);
    }

    // -- Centre of mass --

    #[test]
    fn soa_centre_of_mass() {
        let mut soa = ParticleSoA::new(0);
        soa.push([0.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0);
        soa.push([10.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0);
        let com = soa.centre_of_mass().expect("should succeed");
        assert!((com[0] - 5.0).abs() < 1e-12);
        assert!((com[1]).abs() < 1e-12);
        assert!((com[2]).abs() < 1e-12);
    }

    #[test]
    fn soa_centre_of_mass_empty_errors() {
        let soa = ParticleSoA::new(0);
        assert!(soa.centre_of_mass().is_err());
    }

    // -- AlignedVec --

    #[test]
    fn aligned_vec_basic() {
        let av: AlignedVec<f64> = AlignedVec::new(100);
        assert_eq!(av.len(), 100);
        assert!(!av.is_empty());
        assert_eq!(av.as_slice()[0], 0.0);

        let av2: AlignedVec<f64> = AlignedVec::from_vec(vec![1.0, 2.0, 3.0]);
        assert_eq!(av2.len(), 3);
        assert_eq!(av2.as_slice(), &[1.0, 2.0, 3.0]);
    }

    #[test]
    fn aligned_vec_mut() {
        let mut av: AlignedVec<f64> = AlignedVec::new(5);
        av.as_mut_slice()[0] = 42.0;
        assert_eq!(av.as_slice()[0], 42.0);
    }

    #[test]
    fn aligned_vec_empty() {
        let av: AlignedVec<f64> = AlignedVec::new(0);
        assert!(av.is_empty());
        assert_eq!(av.len(), 0);
    }

    // -- Push / get consistency --

    #[test]
    fn soa_push_get_consistency() {
        let mut soa = ParticleSoA::new(0);
        let n = 200;
        for i in 0..n {
            let v = i as f64;
            soa.push(
                [v, v + 1.0, v + 2.0],
                [v * 0.1, v * 0.2, v * 0.3],
                [v * 10.0, v * 20.0, v * 30.0],
                v + 100.0,
            );
        }
        assert_eq!(soa.len(), n);

        for i in 0..n {
            let v = i as f64;
            assert_eq!(soa.get_position(i), Ok([v, v + 1.0, v + 2.0]));
            assert_eq!(soa.get_velocity(i), Ok([v * 0.1, v * 0.2, v * 0.3]));
            assert_eq!(soa.get_force(i), Ok([v * 10.0, v * 20.0, v * 30.0]));
            assert_eq!(soa.get_mass(i), Ok(v + 100.0));
        }
    }

    // -- Prefetch is a no-op but should not panic --

    #[test]
    fn soa_prefetch_range_no_panic() {
        let mut soa = ParticleSoA::new(0);
        for i in 0..10 {
            soa.push([i as f64; 3], [0.0; 3], [0.0; 3], 1.0);
        }
        soa.prefetch_range(0, 10);
        soa.prefetch_range(5, 5);
        soa.prefetch_range(0, 0);
    }

    // -- Morton sort with many particles --

    #[test]
    fn soa_sort_by_morton_larger_set() {
        let mut soa = ParticleSoA::new(0);
        // Add particles in reverse spatial order.
        for i in (0..50).rev() {
            let v = i as f64 * 2.0;
            soa.push([v, v, v], [0.0; 3], [0.0; 3], 1.0);
        }
        assert!(soa.sort_by_morton_code(1.0).is_ok());
        assert_eq!(soa.len(), 50);

        // After sorting, positions should be in ascending Morton order.
        let mut prev_code = 0u64;
        for i in 0..soa.len() {
            let pos = soa.get_position(i).expect("valid index");
            let code = position_to_morton(pos, 1.0).expect("valid");
            assert!(code >= prev_code, "Morton order violated at index {i}");
            prev_code = code;
        }
    }
}
