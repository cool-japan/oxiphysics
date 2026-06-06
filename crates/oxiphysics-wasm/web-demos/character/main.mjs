// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// Character controller demo — kinematic capsule on a sine-wave heightfield.
// Uses WasmCharacterControllerJs (WASM) + WebGL2 terrain rendering.

import init, { WasmCharacterControllerJs } from '../../pkg-web/oxiphysics_wasm.js';
import { createCanvas, createProgram, Camera, startRenderLoop } from '../shared/canvas-util.mjs';
import { createFpsCounter, createSlider } from '../shared/ui-util.mjs';

// ---------------------------------------------------------------------------
// Heightfield geometry
// ---------------------------------------------------------------------------

const GRID_W = 64;
const GRID_H = 64;
const CELL_SIZE = 0.5;

function sampleHeight(x, z) {
  return Math.sin(x * 0.5) * Math.cos(z * 0.5) * 1.5;
}

function buildTerrainMesh() {
  const verts = [];
  const indices = [];
  for (let iz = 0; iz <= GRID_H; iz++) {
    for (let ix = 0; ix <= GRID_W; ix++) {
      const x = (ix - GRID_W / 2) * CELL_SIZE;
      const z = (iz - GRID_H / 2) * CELL_SIZE;
      const y = sampleHeight(x, z);
      verts.push(x, y, z);
    }
  }
  const stride = GRID_W + 1;
  for (let iz = 0; iz < GRID_H; iz++) {
    for (let ix = 0; ix < GRID_W; ix++) {
      const a = iz * stride + ix;
      const b = a + 1;
      const c = a + stride;
      const d = c + 1;
      indices.push(a, c, b, b, c, d);
    }
  }
  return { verts: new Float32Array(verts), indices: new Uint32Array(indices) };
}

// ---------------------------------------------------------------------------
// Vertex / fragment shaders
// ---------------------------------------------------------------------------

const VS = `#version 300 es
in vec3 aPos;
uniform mat4 uMVP;
out float vHeight;
void main() {
  vHeight = aPos.y;
  gl_Position = uMVP * vec4(aPos, 1.0);
}`;

const FS = `#version 300 es
precision mediump float;
in float vHeight;
out vec4 fragColor;
void main() {
  float t = clamp(vHeight * 0.3 + 0.5, 0.0, 1.0);
  vec3 col = mix(vec3(0.15, 0.35, 0.15), vec3(0.7, 0.65, 0.55), t);
  fragColor = vec4(col, 1.0);
}`;

const CHAR_VS = `#version 300 es
in vec3 aPos;
uniform mat4 uMVP;
uniform vec3 uCharPos;
void main() {
  gl_Position = uMVP * vec4(aPos + uCharPos, 1.0);
}`;

const CHAR_FS = `#version 300 es
precision mediump float;
out vec4 fragColor;
void main() { fragColor = vec4(1.0, 0.4, 0.2, 1.0); }`;

// ---------------------------------------------------------------------------
// Character capsule mesh (very simple — 16-sided cylinder)
// ---------------------------------------------------------------------------

function buildCapsuleMesh(radius, halfHeight, segs) {
  const verts = [];
  const idx = [];
  // Simple: just build a vertical cylinder for visual representation
  for (let i = 0; i <= segs; i++) {
    const a = (i / segs) * Math.PI * 2;
    const x = Math.cos(a) * radius;
    const z = Math.sin(a) * radius;
    verts.push(x, -halfHeight - radius, z);
    verts.push(x,  halfHeight + radius, z);
  }
  for (let i = 0; i < segs; i++) {
    const b = i * 2, t = b + 1;
    const bn = b + 2, tn = t + 2;
    idx.push(b, t, tn, b, tn, bn);
  }
  return { verts: new Float32Array(verts), indices: new Uint32Array(idx) };
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

async function main() {
  await init();

  const { canvas, gl } = createCanvas('canvas');
  canvas.width = window.innerWidth;
  canvas.height = window.innerHeight;

  const status = document.getElementById('status');
  status.textContent = 'Running';

  const ui = document.getElementById('ui');
  const fps = createFpsCounter(ui);
  let gravity = 9.81;
  let speed = 4.0;
  createSlider(ui, 'Gravity', 0, 20, gravity, v => { gravity = v; });
  createSlider(ui, 'Speed', 1, 12, speed, v => { speed = v; });

  // Character controller (radius 0.4 m, half-height 0.9 m, spawn above terrain)
  const spawnY = sampleHeight(0, 0) + 2.5;
  const ctrl = WasmCharacterControllerJs.new(0.0, spawnY, 0.0, 0.4, 0.9);

  // Terrain
  const terrain = buildTerrainMesh();
  const terrainVAO = gl.createVertexArray();
  gl.bindVertexArray(terrainVAO);
  const terrainVBO = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, terrainVBO);
  gl.bufferData(gl.ARRAY_BUFFER, terrain.verts, gl.STATIC_DRAW);
  const terrainEBO = gl.createBuffer();
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, terrainEBO);
  gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, terrain.indices, gl.STATIC_DRAW);
  const terrProg = createProgram(gl, VS, FS);
  gl.enableVertexAttribArray(gl.getAttribLocation(terrProg, 'aPos'));
  gl.vertexAttribPointer(gl.getAttribLocation(terrProg, 'aPos'), 3, gl.FLOAT, false, 12, 0);
  gl.bindVertexArray(null);

  // Capsule mesh
  const cap = buildCapsuleMesh(0.4, 0.9, 16);
  const capVAO = gl.createVertexArray();
  gl.bindVertexArray(capVAO);
  const capVBO = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, capVBO);
  gl.bufferData(gl.ARRAY_BUFFER, cap.verts, gl.STATIC_DRAW);
  const capEBO = gl.createBuffer();
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, capEBO);
  gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, cap.indices, gl.STATIC_DRAW);
  const capProg = createProgram(gl, CHAR_VS, CHAR_FS);
  gl.enableVertexAttribArray(gl.getAttribLocation(capProg, 'aPos'));
  gl.vertexAttribPointer(gl.getAttribLocation(capProg, 'aPos'), 3, gl.FLOAT, false, 12, 0);
  gl.bindVertexArray(null);

  // Camera
  const cam = new Camera({ position: [0, spawnY + 4, 8], pitch: -0.4 });
  let keys = {};
  window.addEventListener('keydown', e => { keys[e.code] = true; });
  window.addEventListener('keyup',   e => { keys[e.code] = false; });

  // Mouse look
  let dragging = false;
  canvas.addEventListener('mousedown', () => { dragging = true; canvas.requestPointerLock?.(); });
  canvas.addEventListener('mouseup',   () => { dragging = false; });
  canvas.addEventListener('mousemove', e => {
    if (dragging || document.pointerLockElement === canvas) {
      cam.yaw   += e.movementX * 0.003;
      cam.pitch  = Math.max(-1.4, Math.min(1.4, cam.pitch - e.movementY * 0.003));
    }
  });
  document.addEventListener('pointerlockchange', () => {});

  startRenderLoop(dt => {
    const pos = ctrl.get_position();

    // Input
    let moveX = 0, moveZ = 0;
    if (keys['KeyW'] || keys['ArrowUp'])    moveZ -= 1;
    if (keys['KeyS'] || keys['ArrowDown'])  moveZ += 1;
    if (keys['KeyA'] || keys['ArrowLeft'])  moveX -= 1;
    if (keys['KeyD'] || keys['ArrowRight']) moveX += 1;
    if (keys['Space'] && ctrl.is_grounded()) ctrl.jump(5.0);

    // Rotate move direction by camera yaw
    const cy = Math.cos(cam.yaw), sy = Math.sin(cam.yaw);
    const dx = (cy * moveX - sy * moveZ) * speed * dt;
    const dz = (sy * moveX + cy * moveZ) * speed * dt;

    ctrl.apply_gravity(gravity, dt);
    const vy = ctrl.get_velocity()[1];

    // Simple ground check against heightfield
    const charPos = ctrl.get_position();
    const groundY = sampleHeight(charPos[0], charPos[2]);
    const hits = [];
    if (charPos[1] - 0.4 < groundY + 0.05) {
      // Provide a ground sweep hit
      hits.push(1.0, 0.0, 1.0, 0.0); // toi, nx, ny, nz
    }

    ctrl.move_and_slide(dx, vy * dt, dz, hits);

    // Clamp to terrain
    const newPos = ctrl.get_position();
    if (newPos[1] < groundY + 0.4) {
      ctrl.set_position(newPos[0], groundY + 0.4, newPos[2]);
      ctrl.set_velocity(ctrl.get_velocity()[0], 0, ctrl.get_velocity()[2]);
    }

    // Camera follow
    cam.position[0] = newPos[0] + Math.sin(cam.yaw) * 6;
    cam.position[1] = newPos[1] + 3;
    cam.position[2] = newPos[2] + Math.cos(cam.yaw) * 6;

    // Render
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.clearColor(0.12, 0.13, 0.18, 1.0);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    const mvp = cam.mvp(canvas.width / canvas.height);
    const mvpLoc = gl.getUniformLocation(terrProg, 'uMVP');

    // Terrain
    gl.useProgram(terrProg);
    gl.uniformMatrix4fv(mvpLoc, false, mvp);
    gl.bindVertexArray(terrainVAO);
    gl.drawElements(gl.TRIANGLES, terrain.indices.length, gl.UNSIGNED_INT, 0);

    // Character capsule
    const final = ctrl.get_position();
    gl.useProgram(capProg);
    gl.uniformMatrix4fv(gl.getUniformLocation(capProg, 'uMVP'), false, mvp);
    gl.uniform3f(gl.getUniformLocation(capProg, 'uCharPos'), final[0], final[1], final[2]);
    gl.bindVertexArray(capVAO);
    gl.drawElements(gl.TRIANGLES, cap.indices.length, gl.UNSIGNED_INT, 0);
    gl.bindVertexArray(null);

    fps.update(dt);
  });
}

main().catch(err => {
  document.getElementById('status').textContent = `Error: ${err.message}`;
  console.error(err);
});
