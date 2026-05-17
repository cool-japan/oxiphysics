// Copyright (C) 2026 COOLJAPAN OU (Team Kitasan)
// SPDX-License-Identifier: Apache-2.0
//
// NavMesh click-to-move demo — 10×10 grid navmesh with A* pathfinding.
// Click anywhere on the floor to set a goal; the cursor animates along the path.

import init, { WasmNavMeshJs } from '../../pkg-web/oxiphysics_wasm.js';
import { createFpsCounter } from '../shared/ui-util.mjs';

const GRID_ROWS = 10;
const GRID_COLS = 10;
const CELL = 1.0;
const WORLD_W = GRID_COLS * CELL;
const WORLD_H = GRID_ROWS * CELL;
const AGENT_SPEED = 3.0;  // world units per second

async function main() {
  await init();

  const canvas = document.getElementById('canvas');
  const ctx = canvas.getContext('2d');
  if (!ctx) throw new Error('Canvas 2D context unavailable');

  const ui = document.getElementById('ui');
  document.getElementById('status').textContent = 'Running';
  const fps = createFpsCounter(ui);

  function resize() {
    canvas.width  = window.innerWidth;
    canvas.height = window.innerHeight;
  }
  window.addEventListener('resize', resize);
  resize();

  // Build navmesh
  const nav = WasmNavMeshJs.grid(GRID_ROWS, GRID_COLS, CELL);

  // Agent state
  let agentPos = [0.5, 0.0, 0.5];  // world coords
  let path = [];       // array of [x, z] waypoints
  let pathIdx = 0;

  // World ↔ screen transforms (top-down, Y-up world displayed as 2D x-z)
  function worldToScreen(wx, wz) {
    const margin = 40;
    const scaleX = (canvas.width  - margin * 2) / WORLD_W;
    const scaleZ = (canvas.height - margin * 2) / WORLD_H;
    const scale = Math.min(scaleX, scaleZ);
    const ox = (canvas.width  - WORLD_W * scale) / 2;
    const oz = (canvas.height - WORLD_H * scale) / 2;
    return [ox + wx * scale, oz + wz * scale];
  }

  function screenToWorld(px, py) {
    const margin = 40;
    const scaleX = (canvas.width  - margin * 2) / WORLD_W;
    const scaleZ = (canvas.height - margin * 2) / WORLD_H;
    const scale = Math.min(scaleX, scaleZ);
    const ox = (canvas.width  - WORLD_W * scale) / 2;
    const oz = (canvas.height - WORLD_H * scale) / 2;
    return [(px - ox) / scale, (py - oz) / scale];
  }

  canvas.addEventListener('click', e => {
    const rect = canvas.getBoundingClientRect();
    const [wx, wz] = screenToWorld(e.clientX - rect.left, e.clientY - rect.top);
    // Clamp to navmesh bounds
    const gx = Math.max(0.1, Math.min(WORLD_W - 0.1, wx));
    const gz = Math.max(0.1, Math.min(WORLD_H - 0.1, wz));

    const flat = nav.find_path_flat(
      agentPos[0], 0.0, agentPos[2],
      gx, 0.0, gz,
      0.2,
    );
    if (flat.length >= 6) {
      path = [];
      for (let i = 0; i < flat.length; i += 3) {
        path.push([flat[i], flat[i + 2]]);  // x and z
      }
      pathIdx = 1;  // skip start (already at agentPos)
    }
  });

  // Draw helpers
  function drawGrid(ctx) {
    ctx.strokeStyle = 'rgba(40,80,40,0.6)';
    ctx.lineWidth = 0.5;
    for (let r = 0; r <= GRID_ROWS; r++) {
      const [x0, y0] = worldToScreen(0, r * CELL);
      const [x1, y1] = worldToScreen(WORLD_W, r * CELL);
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    }
    for (let c = 0; c <= GRID_COLS; c++) {
      const [x0, y0] = worldToScreen(c * CELL, 0);
      const [x1, y1] = worldToScreen(c * CELL, WORLD_H);
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    }
  }

  function drawPath(ctx) {
    if (path.length < 2) return;
    ctx.strokeStyle = 'rgba(100,200,255,0.7)';
    ctx.lineWidth = 2;
    ctx.setLineDash([6, 4]);
    ctx.beginPath();
    const [sx, sy] = worldToScreen(agentPos[0], agentPos[2]);
    ctx.moveTo(sx, sy);
    for (let i = pathIdx; i < path.length; i++) {
      const [px, pz] = worldToScreen(path[i][0], path[i][1]);
      ctx.lineTo(px, pz);
    }
    ctx.stroke();
    ctx.setLineDash([]);

    // Draw waypoints
    for (let i = pathIdx; i < path.length; i++) {
      const [px, pz] = worldToScreen(path[i][0], path[i][1]);
      const isGoal = i === path.length - 1;
      ctx.fillStyle = isGoal ? '#ffcc00' : 'rgba(100,200,255,0.5)';
      ctx.beginPath();
      ctx.arc(px, pz, isGoal ? 7 : 3.5, 0, Math.PI * 2);
      ctx.fill();
    }
  }

  function drawAgent(ctx) {
    const [sx, sy] = worldToScreen(agentPos[0], agentPos[2]);
    ctx.fillStyle = '#ff5533';
    ctx.strokeStyle = '#fff';
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    ctx.arc(sx, sy, 9, 0, Math.PI * 2);
    ctx.fill();
    ctx.stroke();
  }

  let last = performance.now();
  function frame(now) {
    const dt = Math.min((now - last) / 1000, 0.05);
    last = now;

    // Move agent along path
    if (path.length > 0 && pathIdx < path.length) {
      const target = path[pathIdx];
      const dx = target[0] - agentPos[0];
      const dz = target[1] - agentPos[2];
      const dist = Math.sqrt(dx * dx + dz * dz);
      const step = AGENT_SPEED * dt;
      if (dist <= step) {
        agentPos[0] = target[0];
        agentPos[2] = target[1];
        pathIdx++;
      } else {
        agentPos[0] += (dx / dist) * step;
        agentPos[2] += (dz / dist) * step;
      }
    }

    // Render
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Background
    ctx.fillStyle = '#0d1410';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // Floor fill
    const [fx0, fz0] = worldToScreen(0, 0);
    const [fx1, fz1] = worldToScreen(WORLD_W, WORLD_H);
    ctx.fillStyle = 'rgba(20, 50, 25, 0.7)';
    ctx.fillRect(fx0, fz0, fx1 - fx0, fz1 - fz0);

    drawGrid(ctx);
    drawPath(ctx);
    drawAgent(ctx);

    fps.update(dt);
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
}

main().catch(err => {
  document.getElementById('status').textContent = `Error: ${err.message}`;
  console.error(err);
});
