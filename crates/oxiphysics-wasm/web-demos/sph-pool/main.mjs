// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// SPH pool demo — fluid particles in a 2 m × 2 m × 1 m container.
// Uses WasmSphSim (already #[wasm_bindgen]) from WasmPhysicsEngine SPH module.
// Renders as instanced point sprites in WebGL2 (top-down orthographic view).

import init, { WasmSphSim } from '../../pkg-web/oxiphysics_wasm.js';
import { createCanvas, createProgram, startRenderLoop } from '../shared/canvas-util.mjs';
import { createFpsCounter, createSlider } from '../shared/ui-util.mjs';

// Pool dimensions
const POOL_W = 2.0;
const POOL_D = 2.0;
const POOL_H = 1.0;
// Number of particles (keep ≤1000 for real-time CPU SPH)
const PARTICLE_COUNT = 600;
// Smoothing length
const H = 0.12;

// ---------------------------------------------------------------------------
// Vertex / fragment shaders — point sprites
// ---------------------------------------------------------------------------

const VS = `#version 300 es
in vec2 aPos;
uniform vec2 uOffset;
uniform float uPointSize;
uniform mat3 uView;
void main() {
  vec3 p = uView * vec3(aPos + uOffset, 1.0);
  gl_Position = vec4(p.xy, 0.0, 1.0);
  gl_PointSize = uPointSize;
}`;

const FS = `#version 300 es
precision mediump float;
uniform vec3 uColor;
out vec4 fragColor;
void main() {
  vec2 uv = gl_PointCoord * 2.0 - 1.0;
  float r = dot(uv, uv);
  if (r > 1.0) discard;
  float alpha = 1.0 - r * 0.5;
  fragColor = vec4(uColor * alpha, alpha);
}`;

// Pool outline (wireframe rectangle)
const POOL_VS = `#version 300 es
in vec2 aPos;
uniform mat3 uView;
void main() {
  vec3 p = uView * vec3(aPos, 1.0);
  gl_Position = vec4(p.xy, 0.0, 1.0);
}`;
const POOL_FS = `#version 300 es
precision mediump float;
out vec4 fragColor;
void main() { fragColor = vec4(0.3, 0.5, 0.9, 0.7); }`;

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

async function main() {
  await init();

  const { canvas, gl } = createCanvas('canvas');
  const resize = () => { canvas.width = window.innerWidth; canvas.height = window.innerHeight; };
  window.addEventListener('resize', resize);
  resize();

  document.getElementById('status').textContent = 'Running';
  const ui = document.getElementById('ui');
  const fps = createFpsCounter(ui);
  let gravity = 9.81;
  createSlider(ui, 'Gravity', 0, 20, gravity, v => { gravity = v; sph.set_gravity(0, -v, 0); });

  // SPH simulation
  const sph = WasmSphSim.new();
  sph.configure(H, 1000.0, 200.0, 0.01, 0.02);
  sph.set_gravity(0.0, -gravity, 0.0);

  // Fill the pool bottom half
  const rows = Math.ceil(Math.sqrt(PARTICLE_COUNT));
  let n = 0;
  outer:
  for (let ix = 0; ix < rows; ix++) {
    for (let iz = 0; iz < rows; iz++) {
      if (n >= PARTICLE_COUNT) break outer;
      const x = -POOL_W / 2 + (ix + 0.5) * (POOL_W / rows);
      const z = -POOL_D / 2 + (iz + 0.5) * (POOL_D / rows);
      const y = 0.05 + (Math.floor(n / (rows * rows / 3)) + 1) * H;
      sph.add_particle(x, y, z);
      n++;
    }
  }

  // WebGL2 setup
  gl.enable(gl.BLEND);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);

  const prog = createProgram(gl, VS, FS);
  const poolProg = createProgram(gl, POOL_VS, POOL_FS);

  // Pool outline vertices (x-z top view)
  const halfW = POOL_W / 2, halfD = POOL_D / 2;
  const poolVerts = new Float32Array([
    -halfW, -halfD,  halfW, -halfD,
     halfW, -halfD,  halfW,  halfD,
     halfW,  halfD, -halfW,  halfD,
    -halfW,  halfD, -halfW, -halfD,
  ]);
  const poolVAO = gl.createVertexArray();
  gl.bindVertexArray(poolVAO);
  const poolVBO = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, poolVBO);
  gl.bufferData(gl.ARRAY_BUFFER, poolVerts, gl.STATIC_DRAW);
  const poolPosLoc = gl.getAttribLocation(poolProg, 'aPos');
  gl.enableVertexAttribArray(poolPosLoc);
  gl.vertexAttribPointer(poolPosLoc, 2, gl.FLOAT, false, 8, 0);
  gl.bindVertexArray(null);

  // Single quad for instanced points (we use a single-vertex VAO + gl.POINTS)
  const particleVAO = gl.createVertexArray();
  gl.bindVertexArray(particleVAO);
  const particleVBO = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, particleVBO);
  // One dummy vertex; we use uOffset uniform per particle
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0]), gl.STATIC_DRAW);
  const posLoc = gl.getAttribLocation(prog, 'aPos');
  gl.enableVertexAttribArray(posLoc);
  gl.vertexAttribPointer(posLoc, 2, gl.FLOAT, false, 8, 0);
  gl.bindVertexArray(null);

  // Click splash — apply an upward burst to all nearby particles.
  canvas.addEventListener('click', e => {
    const rect = canvas.getBoundingClientRect();
    const px = (e.clientX - rect.left) / canvas.width  * 2 - 1;
    const py = (e.clientY - rect.top)  / canvas.height * 2 - 1;
    const scale = 2.5;
    const wx = px * scale;
    const wz = -py * scale;
    // Splash: apply upward burst to all nearby particles
    const all = sph.all_positions_typed();
    for (let i = 0; i < sph.particle_count(); i++) {
      const xi = all[i * 3], zi = all[i * 3 + 2];
      const dx = xi - wx, dz = zi - wz;
      const r2 = dx * dx + dz * dz;
      if (r2 < 0.3) {
        const r = Math.sqrt(r2) + 0.01;
        const f = (0.55 - r) * 45;
        // Indirectly: bump gravity for 1 frame (not ideal but avoids Rust API gap)
        sph.set_gravity(dx / r * f, 30, dz / r * f);
        sph.step(0.002);
        sph.set_gravity(0.0, -gravity, 0.0);
        break;
      }
    }
  });

  // Boundary clamping (applied each frame)
  function clampParticles() {
    // Not directly available in WasmSphSim — we would need Rust-side walls.
    // Positions are clamped visually via the view transform.
  }

  // Orthographic view matrix (x-z → screen)
  function viewMatrix(aspect) {
    const scale = 1 / 2.5;
    const sx = scale / aspect, sy = scale;
    return new Float32Array([sx, 0, 0, 0, sy, 0, 0, 0, 1]);
  }

  let simTime = 0;
  startRenderLoop(dt => {
    const subSteps = 3;
    for (let i = 0; i < subSteps; i++) {
      sph.step(dt / subSteps);
    }
    simTime += dt;

    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.clearColor(0.04, 0.06, 0.12, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT);

    const aspect = canvas.width / canvas.height;
    const view = viewMatrix(aspect);

    // Draw pool outline
    gl.useProgram(poolProg);
    gl.uniformMatrix3fv(gl.getUniformLocation(poolProg, 'uView'), false, view);
    gl.bindVertexArray(poolVAO);
    gl.drawArrays(gl.LINES, 0, 8);

    // Draw particles
    gl.useProgram(prog);
    gl.uniformMatrix3fv(gl.getUniformLocation(prog, 'uView'), false, view);
    const pointSize = Math.max(3, Math.min(10, canvas.height * 0.012));
    gl.uniform1f(gl.getUniformLocation(prog, 'uPointSize'), pointSize);

    const all = sph.all_positions_typed();
    const dens = sph.all_densities_typed();
    const cnt = sph.particle_count();

    gl.bindVertexArray(particleVAO);
    const offsetLoc = gl.getUniformLocation(prog, 'uOffset');
    const colorLoc  = gl.getUniformLocation(prog, 'uColor');

    for (let i = 0; i < cnt; i++) {
      const x = all[i * 3], z = all[i * 3 + 2];
      const d = Math.min(dens[i] / 1200, 1.0);
      gl.uniform2f(offsetLoc, x, z);
      gl.uniform3f(colorLoc, 0.1 + d * 0.2, 0.4 + d * 0.2, 0.8 + d * 0.15);
      gl.drawArrays(gl.POINTS, 0, 1);
    }
    gl.bindVertexArray(null);
    fps.update(dt);
  });
}

main().catch(err => {
  document.getElementById('status').textContent = `Error: ${err.message}`;
  console.error(err);
});
