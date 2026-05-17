//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

use super::types::ShaderMetadata;

/// WGSL compute shader for SPH density computation.
pub const SPH_DENSITY_WGSL: &str = r#"
// Binding layout:
// @group(0) @binding(0) var<storage, read> positions: array<vec4<f32>>;
// @group(0) @binding(1) var<storage, read> masses: array<f32>;
// @group(0) @binding(2) var<storage, read_write> densities: array<f32>;
// @group(0) @binding(3) var<uniform> params: SphParams;

pub(super) struct SphParams {
    n_particles: u32,
    h: f32,
    h2: f32,
    h3: f32,
}

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let i = id.x;
    if (i >= params.n_particles) { return; }

    var density: f32 = 0.0;
    let pi = positions[i].xyz;

    for (var j: u32 = 0u; j < params.n_particles; j++) {
        let r = pi - positions[j].xyz;
        let r2 = dot(r, r);
        if (r2 < params.h2) {
            let q = 1.0 - r2 / params.h2;
            density += masses[j] * q * q * q;
        }
    }

    densities[i] = density * 315.0 / (64.0 * 3.14159265 * params.h3);
}
"#;
/// WGSL shader for particle integration (velocity Verlet half-step).
pub const INTEGRATE_WGSL: &str = r#"
pub(super) struct Particle { pos: vec4<f32>, vel: vec4<f32>, }
pub(super) struct IntegParams { dt: f32, gravity_y: f32, n: u32, _pad: u32, }

@group(0) @binding(0) var<storage, read_write> particles: array<Particle>;
@group(0) @binding(1) var<storage, read> forces: array<vec4<f32>>;
@group(0) @binding(2) var<uniform> params: IntegParams;

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let i = id.x;
    if (i >= params.n) { return; }

    var p = particles[i];
    let f = forces[i].xyz + vec3<f32>(0.0, params.gravity_y, 0.0);
    p.vel = vec4<f32>(p.vel.xyz + f * params.dt, 0.0);
    p.pos = vec4<f32>(p.pos.xyz + p.vel.xyz * params.dt, 1.0);
    particles[i] = p;
}
"#;
/// WGSL shader for LBM BGK collision (D2Q9).
pub const LBM_BGK_D2Q9_WGSL: &str = r#"
// D2Q9 LBM BGK collision kernel
// Layout: f[node * 9 + direction]
pub(super) struct LbmParams { nx: u32, ny: u32, tau: f32, _pad: u32, }

@group(0) @binding(0) var<storage, read_write> f: array<f32>;
@group(0) @binding(1) var<storage, read_write> rho: array<f32>;
@group(0) @binding(2) var<storage, read_write> ux: array<f32>;
@group(0) @binding(3) var<storage, read_write> uy: array<f32>;
@group(0) @binding(4) var<uniform> params: LbmParams;

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let idx = id.x;
    if (idx >= params.nx * params.ny) { return; }

    // Compute macroscopic values
    var r: f32 = 0.0;
    var u: f32 = 0.0;
    var v: f32 = 0.0;
    let w = array<f32, 9>(4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0);
    let ex = array<f32, 9>(0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0);
    let ey = array<f32, 9>(0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0);

    for (var i: u32 = 0u; i < 9u; i++) {
        let fi = f[idx * 9u + i];
        r += fi;
        u += ex[i] * fi;
        v += ey[i] * fi;
    }
    if (r > 0.0) { u /= r; v /= r; }

    rho[idx] = r; ux[idx] = u; uy[idx] = v;

    // BGK collision
    for (var i: u32 = 0u; i < 9u; i++) {
        let eu = ex[i]*u + ey[i]*v;
        let feq = w[i] * r * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*(u*u+v*v));
        f[idx * 9u + i] += (feq - f[idx * 9u + i]) / params.tau;
    }
}
"#;
/// WGSL compute shader for parallel cell-list construction.
pub const CELL_LIST_WGSL: &str = r#"
pub(super) struct CellParams {
    nx:      u32,
    ny:      u32,
    nz:      u32,
    cell_size: f32,
    n_atoms: u32,
    box_x:   f32,
    box_y:   f32,
    box_z:   f32,
}

@group(0) @binding(0) var<storage, read>       positions:    array<vec4<f32>>;
@group(0) @binding(1) var<storage, read_write> cell_indices: array<u32>;
@group(0) @binding(2) var<uniform>             params:       CellParams;

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let atom_id = id.x;
    if (atom_id >= params.n_atoms) { return; }

    let pos = positions[atom_id].xyz;

    let ix = clamp(u32(pos.x / params.cell_size), 0u, params.nx - 1u);
    let iy = clamp(u32(pos.y / params.cell_size), 0u, params.ny - 1u);
    let iz = clamp(u32(pos.z / params.cell_size), 0u, params.nz - 1u);

    cell_indices[atom_id] = iz * params.ny * params.nx
                          + iy * params.nx
                          + ix;
}
"#;
/// WGSL compute shader that computes a sphere SDF on the GPU.
pub const SDF_COMPUTE_WGSL: &str = r#"
pub(super) struct SdfSphereParams {
    nx:       u32,
    ny:       u32,
    nz:       u32,
    dx:       f32,
    origin_x: f32,
    origin_y: f32,
    origin_z: f32,
    cx:       f32,
    cy:       f32,
    cz:       f32,
    radius:   f32,
    _pad:     u32,
}

@group(0) @binding(0) var<storage, read_write> sdf_values: array<f32>;
@group(0) @binding(1) var<uniform>             params:     SdfSphereParams;

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let total = params.nx * params.ny * params.nz;
    let idx   = id.x;
    if (idx >= total) { return; }

    let i = idx / (params.ny * params.nz);
    let j = (idx / params.nz) % params.ny;
    let k = idx % params.nz;

    let px = params.origin_x + (f32(i) + 0.5) * params.dx;
    let py = params.origin_y + (f32(j) + 0.5) * params.dx;
    let pz = params.origin_z + (f32(k) + 0.5) * params.dx;

    let dist = distance(vec3<f32>(px, py, pz),
                        vec3<f32>(params.cx, params.cy, params.cz));
    sdf_values[idx] = dist - params.radius;
}
"#;
/// WGSL compute shader for LBM streaming step (D3Q19).
pub const LBM_STREAMING_SHADER: &str = r#"
// LBM D3Q19 streaming step
// Distributes f values from source to destination according to lattice velocities.
pub(super) struct StreamParams { nx: u32, ny: u32, nz: u32, _pad: u32, }

@group(0) @binding(0) var<storage, read>       f_src:  array<f32>;
@group(0) @binding(1) var<storage, read_write> f_dst:  array<f32>;
@group(0) @binding(2) var<uniform>             params: StreamParams;

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let nx = params.nx;
    let ny = params.ny;
    let nz = params.nz;
    let total = nx * ny * nz;
    let idx = id.x;
    if (idx >= total) { return; }

    let iz = idx / (ny * nx);
    let iy = (idx / nx) % ny;
    let ix = idx % nx;

    // D3Q19 lattice velocities (abbreviated: direction 0 = rest)
    let ex = array<i32, 19>( 0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0);
    let ey = array<i32, 19>( 0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1);
    let ez = array<i32, 19>( 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1);

    for (var q: u32 = 0u; q < 19u; q++) {
        let sx = (i32(ix) - ex[q] + i32(nx)) % i32(nx);
        let sy = (i32(iy) - ey[q] + i32(ny)) % i32(ny);
        let sz = (i32(iz) - ez[q] + i32(nz)) % i32(nz);
        let src_idx = u32(sz) * ny * nx + u32(sy) * nx + u32(sx);
        f_dst[idx * 19u + q] = f_src[src_idx * 19u + q];
    }
}
"#;
/// WGSL compute shader for rigid body semi-implicit Euler integration.
pub const RIGID_INTEGRATE_SHADER: &str = r#"
// Semi-implicit Euler integration for rigid bodies.
pub(super) struct RigidBody {
    pos:       vec4<f32>,  // xyz = position, w = mass
    vel:       vec4<f32>,  // xyz = linear velocity, w unused
    ang_vel:   vec4<f32>,  // xyz = angular velocity, w unused
    force:     vec4<f32>,  // xyz = accumulated force, w unused
    torque:    vec4<f32>,  // xyz = accumulated torque, w unused
}
pub(super) struct IntegRigidParams { dt: f32, n: u32, _pad0: u32, _pad1: u32, }

@group(0) @binding(0) var<storage, read_write> bodies: array<RigidBody>;
@group(0) @binding(1) var<uniform>             params: IntegRigidParams;

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let i = id.x;
    if (i >= params.n) { return; }

    var b = bodies[i];
    let mass = b.pos.w;
    if (mass <= 0.0) { return; }

    // Semi-implicit Euler: update velocity first, then position.
    let acc = b.force.xyz / mass;
    b.vel = vec4<f32>(b.vel.xyz + acc * params.dt, 0.0);
    b.pos = vec4<f32>(b.pos.xyz + b.vel.xyz * params.dt, mass);

    // Clear forces.
    b.force  = vec4<f32>(0.0, 0.0, 0.0, 0.0);
    b.torque = vec4<f32>(0.0, 0.0, 0.0, 0.0);
    bodies[i] = b;
}
"#;
/// WGSL compute shader for bitonic sort used in SAP broadphase.
pub const BROADPHASE_SORT_SHADER: &str = r#"
// Bitonic sort pass for SAP broadphase.
// Each invocation compares and optionally swaps one pair of elements.
pub(super) struct SortParams { n: u32, step: u32, stage: u32, _pad: u32, }

@group(0) @binding(0) var<storage, read_write> keys:   array<f32>;
@group(0) @binding(1) var<storage, read_write> values: array<u32>;
@group(0) @binding(2) var<uniform>             params: SortParams;

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let idx = id.x;
    let n   = params.n;
    if (idx >= n / 2u) { return; }

    let step  = params.step;
    let stage = params.stage;

    let j   = idx & (step - 1u);
    let i   = (idx - j) * 2u + j;
    let i2  = i + step;

    if (i2 >= n) { return; }

    // Ascending sort when the block bit is 0.
    let ascending = ((i / (stage * 2u)) & 1u) == 0u;
    let swap = (keys[i] > keys[i2]) == ascending;

    if (swap) {
        let tmp_k   = keys[i];   keys[i]   = keys[i2];   keys[i2]   = tmp_k;
        let tmp_v   = values[i]; values[i] = values[i2]; values[i2] = tmp_v;
    }
}
"#;
/// WGSL compute shader for SPH pressure force computation.
pub const SPH_FORCE_WGSL: &str = r#"
pub(super) struct SphForceParams {
    n_particles: u32,
    h: f32,
    mu: f32,
    _pad: u32,
}

@group(0) @binding(0) var<storage, read>       positions:  array<vec4<f32>>;
@group(0) @binding(1) var<storage, read>       velocities: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read>       densities:  array<f32>;
@group(0) @binding(3) var<storage, read>       pressures:  array<f32>;
@group(0) @binding(4) var<storage, read>       masses:     array<f32>;
@group(0) @binding(5) var<storage, read_write> forces:     array<vec4<f32>>;
@group(0) @binding(6) var<uniform>             params:     SphForceParams;

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let i = id.x;
    if (i >= params.n_particles) { return; }

    let pi = positions[i].xyz;
    let vi = velocities[i].xyz;
    var force: vec3<f32> = vec3<f32>(0.0, 0.0, 0.0);
    let h2 = params.h * params.h;

    for (var j: u32 = 0u; j < params.n_particles; j++) {
        if (j == i) { continue; }
        let r_vec = pi - positions[j].xyz;
        let r2 = dot(r_vec, r_vec);
        let r = sqrt(r2);
        if (r < 0.0001 || r >= params.h) { continue; }

        // Spiky gradient
        let h_r = params.h - r;
        let grad_mag = -45.0 / (3.14159265 * pow(params.h, 6.0)) * h_r * h_r;
        let grad = r_vec / r * grad_mag;

        // Pressure force
        let p_term = -masses[j] * (pressures[i] + pressures[j]) / (2.0 * densities[j]);
        force += p_term * grad;

        // Viscosity
        let v_lap = 45.0 / (3.14159265 * pow(params.h, 6.0)) * h_r;
        force += params.mu * masses[j] * (velocities[j].xyz - vi) / densities[j] * v_lap;
    }

    forces[i] = vec4<f32>(force, 0.0);
}
"#;
/// WGSL compute shader for boundary condition enforcement.
pub const BOUNDARY_ENFORCE_WGSL: &str = r#"
pub(super) struct BoundaryParams {
    n: u32,
    box_min_x: f32,
    box_min_y: f32,
    box_min_z: f32,
    box_max_x: f32,
    box_max_y: f32,
    box_max_z: f32,
    restitution: f32,
}

@group(0) @binding(0) var<storage, read_write> positions:  array<vec4<f32>>;
@group(0) @binding(1) var<storage, read_write> velocities: array<vec4<f32>>;
@group(0) @binding(2) var<uniform>             params:     BoundaryParams;

@compute @workgroup_size(64)
pub(super) fn main(@builtin(global_invocation_id) id: vec3<u32>) {
    let i = id.x;
    if (i >= params.n) { return; }

    var pos = positions[i].xyz;
    var vel = velocities[i].xyz;
    let e = params.restitution;

    if (pos.x < params.box_min_x) { pos.x = params.box_min_x; vel.x = -vel.x * e; }
    if (pos.x > params.box_max_x) { pos.x = params.box_max_x; vel.x = -vel.x * e; }
    if (pos.y < params.box_min_y) { pos.y = params.box_min_y; vel.y = -vel.y * e; }
    if (pos.y > params.box_max_y) { pos.y = params.box_max_y; vel.y = -vel.y * e; }
    if (pos.z < params.box_min_z) { pos.z = params.box_min_z; vel.z = -vel.z * e; }
    if (pos.z > params.box_max_z) { pos.z = params.box_max_z; vel.z = -vel.z * e; }

    positions[i]  = vec4<f32>(pos, positions[i].w);
    velocities[i] = vec4<f32>(vel, 0.0);
}
"#;
/// Resolves `#include "filename"` directives in shader source.
///
/// This is a simple preprocessor-like mechanism for composing shaders
/// from reusable fragments.
pub fn resolve_includes(source: &str, includes: &HashMap<&str, &str>) -> String {
    let mut result = String::new();
    for line in source.lines() {
        let trimmed = line.trim();
        if trimmed.starts_with("#include") {
            if let Some(start) = trimmed.find('"')
                && let Some(end) = trimmed[start + 1..].find('"')
            {
                let filename = &trimmed[start + 1..start + 1 + end];
                if let Some(content) = includes.get(filename) {
                    result.push_str(content);
                    result.push('\n');
                    continue;
                }
            }
            result.push_str("// UNRESOLVED: ");
            result.push_str(line);
            result.push('\n');
        } else {
            result.push_str(line);
            result.push('\n');
        }
    }
    result
}
/// Validate that a WGSL source has basic structural correctness.
///
/// Checks:
/// 1. The source contains at least one `fn` keyword.
/// 2. Braces are balanced (same number of `{` and `}`).
pub fn validate_wgsl_structure(source: &str) -> bool {
    if !source.contains("fn") {
        return false;
    }
    let open: usize = source.chars().filter(|&c| c == '{').count();
    let close: usize = source.chars().filter(|&c| c == '}').count();
    open == close && open > 0
}
/// Validate that all built-in shader sources contain the `@compute` entry point marker.
pub fn validate_shader_sources() -> bool {
    SPH_DENSITY_WGSL.contains("@compute")
        && INTEGRATE_WGSL.contains("@compute")
        && LBM_BGK_D2Q9_WGSL.contains("@compute")
        && CELL_LIST_WGSL.contains("@compute")
        && SDF_COMPUTE_WGSL.contains("@compute")
        && LBM_STREAMING_SHADER.contains("@compute")
        && RIGID_INTEGRATE_SHADER.contains("@compute")
        && BROADPHASE_SORT_SHADER.contains("@compute")
        && SPH_FORCE_WGSL.contains("@compute")
        && BOUNDARY_ENFORCE_WGSL.contains("@compute")
}
/// Count the total number of binding annotations across all built-in shaders.
pub fn count_total_bindings() -> usize {
    let all_shaders = [
        SPH_DENSITY_WGSL,
        INTEGRATE_WGSL,
        LBM_BGK_D2Q9_WGSL,
        CELL_LIST_WGSL,
        SDF_COMPUTE_WGSL,
        LBM_STREAMING_SHADER,
        RIGID_INTEGRATE_SHADER,
        BROADPHASE_SORT_SHADER,
        SPH_FORCE_WGSL,
        BOUNDARY_ENFORCE_WGSL,
    ];
    all_shaders
        .iter()
        .map(|s| s.matches("@binding(").count())
        .sum()
}
/// Simple 64-bit hash of a string (FNV-1a variant).
pub(super) fn simple_hash(s: &str) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in s.bytes() {
        h ^= byte as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}
/// Produce a mock SPIR-V byte sequence from WGSL source.
///
/// In a real pipeline this would invoke naga or spirv-cross.  Here we
/// produce a deterministic stub that starts with the SPIR-V magic number
/// and encodes basic information from the source.
pub fn mock_compile_to_spirv(source: &str, entry_point: &str) -> Vec<u8> {
    let mut out = vec![0x07u8, 0x23, 0x02, 0x03];
    out.extend_from_slice(&[0x00, 0x01, 0x05, 0x00]);
    let src_len = source.len() as u32;
    out.extend_from_slice(&src_len.to_le_bytes());
    let ep_hash = simple_hash(entry_point) as u32;
    out.extend_from_slice(&ep_hash.to_le_bytes());
    let bindings = source.matches("@binding(").count() as u32;
    out.extend_from_slice(&bindings.to_le_bytes());
    while !out.len().is_multiple_of(4) {
        out.push(0x00);
    }
    out
}
/// Parse `@workgroup_size(x)` or `@workgroup_size(x, y, z)` from WGSL.
pub(super) fn parse_workgroup_size(source: &str) -> [u32; 3] {
    if let Some(pos) = source.find("@workgroup_size(") {
        let rest = &source[pos + 16..];
        if let Some(end) = rest.find(')') {
            let inner = &rest[..end];
            let parts: Vec<u32> = inner
                .split(',')
                .map(|s| s.trim().parse::<u32>().unwrap_or(1))
                .collect();
            return match parts.len() {
                1 => [parts[0], 1, 1],
                2 => [parts[0], parts[1], 1],
                3 => [parts[0], parts[1], parts[2]],
                _ => [1, 1, 1],
            };
        }
    }
    [1, 1, 1]
}
/// Validate that a [`ShaderMetadata`] record is internally consistent.
///
/// Checks:
/// 1. Workgroup size is non-zero in every dimension.
/// 2. Total threads-per-workgroup does not exceed 1024 (common GPU limit).
/// 3. `bind_group_count` is <= 4 (common GPU limit).
/// 4. `entry_point` is non-empty.
pub fn validate_shader_metadata(meta: &ShaderMetadata) -> Result<(), String> {
    if meta.workgroup_size[0] == 0 || meta.workgroup_size[1] == 0 || meta.workgroup_size[2] == 0 {
        return Err("workgroup_size must be non-zero in every dimension".to_string());
    }
    let total = meta.workgroup_size[0]
        .saturating_mul(meta.workgroup_size[1])
        .saturating_mul(meta.workgroup_size[2]);
    if total > 1024 {
        return Err(format!(
            "total threads-per-workgroup {} exceeds 1024",
            total
        ));
    }
    if meta.bind_group_count > 4 {
        return Err(format!(
            "bind_group_count {} exceeds maximum of 4",
            meta.bind_group_count
        ));
    }
    if meta.entry_point.is_empty() {
        return Err("entry_point must not be empty".to_string());
    }
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::shaders::AddressMode;
    use crate::shaders::BindGroupLayout;
    use crate::shaders::ComputeShaderDesc;
    use crate::shaders::DescriptorSetLayout;
    use crate::shaders::DescriptorType;
    use crate::shaders::FilterMode;
    use crate::shaders::PushConstantRange;
    use crate::shaders::RenderPassDesc;
    use crate::shaders::SamplerDesc;
    use crate::shaders::ShaderCache;
    use crate::shaders::ShaderCompilationPipeline;
    use crate::shaders::ShaderHotReloadManager;
    use crate::shaders::ShaderRegistry;
    use crate::shaders::ShaderStage;
    use crate::shaders::ShaderTemplate;
    use crate::shaders::SpecializationMap;
    use crate::shaders::SpirVModule;
    use crate::shaders::TextureFormat;
    use crate::shaders::UniformBufferDesc;
    #[test]
    fn test_shader_sources_valid() {
        assert!(validate_shader_sources());
    }
    #[test]
    fn test_shader_sph_contains_density() {
        assert!(SPH_DENSITY_WGSL.contains("density"));
    }
    #[test]
    fn test_shader_lbm_contains_bgk() {
        assert!(LBM_BGK_D2Q9_WGSL.contains("feq"));
    }
    /// ShaderTemplate replaces placeholders correctly.
    #[test]
    fn test_template_instantiation() {
        let tmpl = ShaderTemplate::new(
            "@compute @workgroup_size({WORKGROUP_SIZE})\nfn main() { var x: {BUFFER_TYPE}; }",
        );
        let mut params = HashMap::new();
        params.insert("WORKGROUP_SIZE", "128");
        params.insert("BUFFER_TYPE", "f32");
        let result = tmpl.instantiate(&params);
        assert!(
            result.contains("@workgroup_size(128)"),
            "Workgroup size not substituted: {result}"
        );
        assert!(
            result.contains("var x: f32;"),
            "Buffer type not substituted: {result}"
        );
        assert!(
            !result.contains("{WORKGROUP_SIZE}"),
            "Placeholder not removed"
        );
    }
    /// ShaderTemplate leaves unknown placeholders unchanged.
    #[test]
    fn test_template_unknown_placeholder() {
        let tmpl = ShaderTemplate::new("fn main() { var x: {UNKNOWN}; }");
        let params = HashMap::new();
        let result = tmpl.instantiate(&params);
        assert!(
            result.contains("{UNKNOWN}"),
            "Unknown placeholder should be preserved"
        );
    }
    /// ShaderRegistry registers and retrieves shaders by name.
    #[test]
    fn test_shader_registry_register_get() {
        let mut reg = ShaderRegistry::new();
        assert!(reg.is_empty());
        let desc = ComputeShaderDesc::new("main", [64, 1, 1], LBM_STREAMING_SHADER);
        reg.register("lbm_stream", desc);
        assert_eq!(reg.len(), 1);
        let retrieved = reg.get("lbm_stream").expect("Shader should be present");
        assert_eq!(retrieved.entry_point, "main");
        assert_eq!(retrieved.workgroup_size, [64, 1, 1]);
        assert!(retrieved.source.contains("@compute"));
    }
    /// ShaderRegistry::with_builtins contains expected shaders.
    #[test]
    fn test_shader_registry_builtins() {
        let reg = ShaderRegistry::with_builtins();
        assert!(reg.get("sph_density").is_some());
        assert!(reg.get("lbm_streaming").is_some());
        assert!(reg.get("rigid_integrate").is_some());
        assert!(reg.get("broadphase_sort").is_some());
        assert!(reg.get("sph_force").is_some());
        assert!(reg.get("boundary_enforce").is_some());
        assert!(reg.get("integrate").is_some());
        assert!(reg.get("lbm_bgk_d2q9").is_some());
    }
    /// validate_wgsl_structure returns true for balanced, valid WGSL.
    #[test]
    fn test_validate_wgsl_structure_valid() {
        let src = "@compute @workgroup_size(64)\nfn main() { let x = 1u; }";
        assert!(
            validate_wgsl_structure(src),
            "Valid WGSL should pass validation"
        );
    }
    /// validate_wgsl_structure returns false for unbalanced braces.
    #[test]
    fn test_validate_wgsl_structure_unbalanced() {
        let src = "fn main() { let x = 1u;";
        assert!(
            !validate_wgsl_structure(src),
            "Unbalanced WGSL should fail validation"
        );
    }
    /// All built-in shaders pass the structural validator.
    #[test]
    fn test_all_builtins_pass_structural_validation() {
        let shaders = [
            SPH_DENSITY_WGSL,
            INTEGRATE_WGSL,
            LBM_BGK_D2Q9_WGSL,
            CELL_LIST_WGSL,
            SDF_COMPUTE_WGSL,
            LBM_STREAMING_SHADER,
            RIGID_INTEGRATE_SHADER,
            BROADPHASE_SORT_SHADER,
            SPH_FORCE_WGSL,
            BOUNDARY_ENFORCE_WGSL,
        ];
        for (i, src) in shaders.iter().enumerate() {
            assert!(
                validate_wgsl_structure(src),
                "Shader {i} failed structural validation"
            );
        }
    }
    #[test]
    fn test_template_placeholders() {
        let tmpl = ShaderTemplate::new("fn main() { var x: {TYPE}; var y: {SIZE}; }");
        let placeholders = tmpl.placeholders();
        assert!(placeholders.contains(&"TYPE".to_string()));
        assert!(placeholders.contains(&"SIZE".to_string()));
    }
    #[test]
    fn test_template_all_placeholders_provided() {
        let tmpl = ShaderTemplate::new("fn main() { var x: {TYPE}; }");
        let mut params = HashMap::new();
        assert!(!tmpl.all_placeholders_provided(&params));
        params.insert("TYPE", "f32");
        assert!(tmpl.all_placeholders_provided(&params));
    }
    #[test]
    fn test_specialization_map() {
        let mut sm = SpecializationMap::new();
        assert!(sm.is_empty());
        sm.define("WORKGROUP_SIZE", "64", "Threads per workgroup");
        sm.define("MAX_NEIGHBORS", "128", "Max neighbors per particle");
        assert_eq!(sm.len(), 2);
        assert_eq!(sm.get("WORKGROUP_SIZE"), Some("64"));
        sm.set("WORKGROUP_SIZE", "128");
        assert_eq!(sm.get("WORKGROUP_SIZE"), Some("128"));
        assert_eq!(sm.get("MAX_NEIGHBORS"), Some("128"));
    }
    #[test]
    fn test_specialization_map_apply() {
        let mut sm = SpecializationMap::new();
        sm.define("WG", "64", "workgroup");
        sm.set("WG", "128");
        let source = "const WG = 64;\nfn main() {}";
        let result = sm.apply(source);
        assert!(
            result.contains("const WG = 128;"),
            "specialization not applied: {result}"
        );
    }
    #[test]
    fn test_resolve_includes() {
        let mut includes = HashMap::new();
        includes.insert("common.wgsl", "// common code\nfn helper() { }");
        let source = "#include \"common.wgsl\"\nfn main() { }";
        let resolved = resolve_includes(source, &includes);
        assert!(
            resolved.contains("common code"),
            "include not resolved: {resolved}"
        );
        assert!(resolved.contains("fn main()"));
    }
    #[test]
    fn test_resolve_includes_unresolved() {
        let includes = HashMap::new();
        let source = "#include \"missing.wgsl\"\nfn main() { }";
        let resolved = resolve_includes(source, &includes);
        assert!(
            resolved.contains("UNRESOLVED"),
            "unresolved include should be marked"
        );
    }
    #[test]
    fn test_shader_cache() {
        let mut cache = ShaderCache::new();
        assert!(cache.is_empty());
        let source = cache.get_or_insert("test", || "fn main() { }".to_string());
        assert_eq!(source, "fn main() { }");
        assert!(cache.contains("test"));
        assert_eq!(cache.len(), 1);
        cache.remove("test");
        assert!(!cache.contains("test"));
    }
    #[test]
    fn test_shader_cache_returns_cached() {
        let mut cache = ShaderCache::new();
        cache.get_or_insert("test", || "first".to_string());
        let second = cache.get_or_insert("test", || "second".to_string());
        assert_eq!(second, "first", "cache should return first insertion");
    }
    #[test]
    fn test_compute_shader_desc_binding_count() {
        let desc = ComputeShaderDesc::new("main", [64, 1, 1], SPH_DENSITY_WGSL);
        assert!(desc.binding_count() > 0);
        assert_eq!(desc.threads_per_workgroup(), 64);
    }
    #[test]
    fn test_registry_unregister() {
        let mut reg = ShaderRegistry::new();
        let desc = ComputeShaderDesc::new("main", [64, 1, 1], "fn main() { }");
        reg.register("test", desc);
        assert!(reg.contains("test"));
        reg.unregister("test");
        assert!(!reg.contains("test"));
    }
    #[test]
    fn test_registry_names() {
        let mut reg = ShaderRegistry::new();
        reg.register(
            "a",
            ComputeShaderDesc::new("main", [64, 1, 1], "fn main() { }"),
        );
        reg.register(
            "b",
            ComputeShaderDesc::new("main", [64, 1, 1], "fn main() { }"),
        );
        let names: Vec<&str> = reg.names().collect();
        assert_eq!(names.len(), 2);
        assert!(names.contains(&"a"));
        assert!(names.contains(&"b"));
    }
    #[test]
    fn test_compilation_pipeline() {
        let mut pipeline = ShaderCompilationPipeline::new();
        let source = "@compute @workgroup_size(64)\nfn main() { let x = 1u; }";
        let result = pipeline.compile("test", source, None).unwrap();
        assert!(result.contains("fn main()"));
        assert_eq!(pipeline.cache_size(), 1);
        let result2 = pipeline.compile("test", source, None).unwrap();
        assert_eq!(result, result2);
    }
    #[test]
    fn test_compilation_pipeline_with_includes() {
        let mut pipeline = ShaderCompilationPipeline::new();
        pipeline.add_include("common.wgsl", "// included code");
        let source = "#include \"common.wgsl\"\n@compute @workgroup_size(64)\nfn main() { }";
        let result = pipeline.compile("test_inc", source, None).unwrap();
        assert!(result.contains("included code"));
    }
    #[test]
    fn test_compilation_pipeline_invalid_source() {
        let mut pipeline = ShaderCompilationPipeline::new();
        let source = "no functions here";
        let result = pipeline.compile("bad", source, None);
        assert!(result.is_err());
    }
    #[test]
    fn test_count_total_bindings() {
        let total = count_total_bindings();
        assert!(
            total > 10,
            "should have many bindings across all shaders, got {total}"
        );
    }
    #[test]
    fn test_sph_force_shader_valid() {
        assert!(validate_wgsl_structure(SPH_FORCE_WGSL));
        assert!(SPH_FORCE_WGSL.contains("@compute"));
        assert!(SPH_FORCE_WGSL.contains("forces"));
    }
    #[test]
    fn test_boundary_enforce_shader_valid() {
        assert!(validate_wgsl_structure(BOUNDARY_ENFORCE_WGSL));
        assert!(BOUNDARY_ENFORCE_WGSL.contains("@compute"));
        assert!(BOUNDARY_ENFORCE_WGSL.contains("restitution"));
    }
    #[test]
    fn test_registry_with_builtins_count() {
        let reg = ShaderRegistry::with_builtins();
        assert!(
            reg.len() >= 8,
            "should have at least 8 built-in shaders, got {}",
            reg.len()
        );
    }
    #[test]
    fn test_uniform_buffer_binding_desc() {
        let desc = UniformBufferDesc::new("SceneParams", 0, 0, 256);
        assert_eq!(desc.name, "SceneParams");
        assert_eq!(desc.group, 0);
        assert_eq!(desc.binding, 0);
        assert_eq!(desc.size_bytes, 256);
    }
    #[test]
    fn test_uniform_buffer_binding_layout() {
        let mut layout = BindGroupLayout::new();
        layout.add_uniform("CameraUbo", 0, 0, 64);
        layout.add_uniform("LightUbo", 0, 1, 128);
        assert_eq!(layout.binding_count(), 2);
    }
    #[test]
    fn test_push_constant_range() {
        let range = PushConstantRange::new(0, 128, ShaderStage::Compute);
        assert_eq!(range.offset, 0);
        assert_eq!(range.size, 128);
        assert_eq!(range.stage, ShaderStage::Compute);
    }
    #[test]
    fn test_push_constants_exceed_limit_detection() {
        let range = PushConstantRange::new(0, 256, ShaderStage::Compute);
        assert!(range.size > 128, "Detected large push constant range");
    }
    #[test]
    fn test_descriptor_set_layout_storage() {
        let mut dsl = DescriptorSetLayout::new(0);
        dsl.add_storage_buffer(0, ShaderStage::Compute, false);
        dsl.add_storage_buffer(1, ShaderStage::Compute, true);
        assert_eq!(dsl.bindings.len(), 2);
        assert!(!dsl.bindings[0].read_only);
        assert!(dsl.bindings[1].read_only);
    }
    #[test]
    fn test_descriptor_set_layout_uniform() {
        let mut dsl = DescriptorSetLayout::new(0);
        dsl.add_uniform_buffer(2, ShaderStage::Vertex);
        assert_eq!(dsl.bindings.len(), 1);
        assert_eq!(dsl.bindings[0].binding, 2);
        assert_eq!(
            dsl.bindings[0].descriptor_type,
            DescriptorType::UniformBuffer
        );
    }
    #[test]
    fn test_sampler_descriptor_linear() {
        let sd = SamplerDesc::linear();
        assert_eq!(sd.filter_min, FilterMode::Linear);
        assert_eq!(sd.filter_mag, FilterMode::Linear);
        assert_eq!(sd.address_mode, AddressMode::ClampToEdge);
    }
    #[test]
    fn test_sampler_descriptor_nearest() {
        let sd = SamplerDesc::nearest();
        assert_eq!(sd.filter_min, FilterMode::Nearest);
        assert_eq!(sd.filter_mag, FilterMode::Nearest);
    }
    #[test]
    fn test_sampler_descriptor_repeat() {
        let sd = SamplerDesc {
            address_mode: AddressMode::Repeat,
            ..SamplerDesc::linear()
        };
        assert_eq!(sd.address_mode, AddressMode::Repeat);
    }
    #[test]
    fn test_render_pass_desc_color_attachment() {
        let rp = RenderPassDesc::new_simple_color();
        assert_eq!(rp.color_attachments.len(), 1);
        assert_eq!(rp.color_attachments[0].format, TextureFormat::Rgba8Unorm);
    }
    #[test]
    fn test_render_pass_desc_depth_attachment() {
        let rp = RenderPassDesc::new_with_depth();
        assert!(rp.depth_attachment.is_some());
        let depth = rp.depth_attachment.unwrap();
        assert_eq!(depth.format, TextureFormat::Depth32Float);
    }
    #[test]
    fn test_render_pass_desc_no_depth() {
        let rp = RenderPassDesc::new_simple_color();
        assert!(rp.depth_attachment.is_none());
    }
    #[test]
    fn test_hot_reload_manager_watch() {
        let mut mgr = ShaderHotReloadManager::new();
        mgr.watch("particles.wgsl", "fn main() { }");
        assert!(mgr.is_watched("particles.wgsl"));
    }
    #[test]
    fn test_hot_reload_manager_update() {
        let mut mgr = ShaderHotReloadManager::new();
        mgr.watch("test.wgsl", "fn main() { let x = 1u; }");
        let changed = mgr.update("test.wgsl", "@compute\nfn main() { let x = 2u; }");
        assert!(changed);
        assert!(mgr.get_source("test.wgsl").unwrap().contains("2u"));
    }
    #[test]
    fn test_hot_reload_manager_no_change() {
        let src = "fn main() { }";
        let mut mgr = ShaderHotReloadManager::new();
        mgr.watch("test.wgsl", src);
        let changed = mgr.update("test.wgsl", src);
        assert!(!changed, "Identical content should not trigger a reload");
    }
    #[test]
    fn test_hot_reload_manager_unwatch() {
        let mut mgr = ShaderHotReloadManager::new();
        mgr.watch("a.wgsl", "fn main() {}");
        mgr.unwatch("a.wgsl");
        assert!(!mgr.is_watched("a.wgsl"));
    }
    #[test]
    fn test_spirv_compilation_mock_produces_bytes() {
        let src = "@compute @workgroup_size(64)\nfn main() { let x = 1u; }";
        let spirv = mock_compile_to_spirv(src, "main");
        assert!(
            !spirv.is_empty(),
            "Mock SPIR-V should produce non-empty bytes"
        );
        assert_eq!(spirv[0], 0x07, "SPIR-V first byte mismatch");
    }
    #[test]
    fn test_spirv_compilation_deterministic() {
        let src = "@compute @workgroup_size(64)\nfn main() { }";
        let s1 = mock_compile_to_spirv(src, "main");
        let s2 = mock_compile_to_spirv(src, "main");
        assert_eq!(s1, s2, "Compilation should be deterministic");
    }
    #[test]
    fn test_spirv_module_reflection() {
        let src = "@group(0) @binding(0) var<storage, read_write> buf: array<f32>;\n@compute @workgroup_size(64)\nfn main() { }";
        let module = SpirVModule::from_wgsl(src);
        assert!(module.entry_points.contains(&"main".to_string()));
        assert!(module.binding_count >= 1);
    }
    #[test]
    fn test_spirv_module_workgroup_size() {
        let src = "@compute @workgroup_size(128)\nfn compute_main() { }";
        let module = SpirVModule::from_wgsl(src);
        assert_eq!(module.workgroup_size[0], 128);
        assert!(module.entry_points.contains(&"compute_main".to_string()));
    }
    #[test]
    fn test_bind_group_layout_serialize() {
        let mut layout = BindGroupLayout::new();
        layout.add_uniform("Ubo", 0, 0, 64);
        layout.add_storage("Sbo", 0, 1, false);
        let desc = layout.to_wgsl_snippet();
        assert!(desc.contains("@binding(0)"));
        assert!(desc.contains("@binding(1)"));
    }
    #[test]
    fn test_bind_group_layout_empty() {
        let layout = BindGroupLayout::new();
        assert_eq!(layout.binding_count(), 0);
        assert!(layout.is_empty());
    }
}
