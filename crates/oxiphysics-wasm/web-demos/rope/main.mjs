// Copyright (C) 2026 COOLJAPAN OU (Team Kitasan)
// SPDX-License-Identifier: Apache-2.0
//
// Rope ragdoll demo — 32-segment Verlet rope, top pinned, mouse-drag impulse.
// Renders as a line strip with tapered width and colour-by-velocity.

import init, { WasmRopeJs } from '../../pkg-web/oxiphysics_wasm.js';
import { createFpsCounter, createSlider } from '../shared/ui-util.mjs';

const LINK_COUNT = 32;
const SEG_LEN = 0.18;
const ANCHOR_X = 0.0;
const ANCHOR_Y = 5.5;

async function main() {
  await init();

  const canvas = document.getElementById('canvas');
  const ctx = canvas.getContext('2d');
  if (!ctx) throw new Error('Canvas 2D not available');

  const ui = document.getElementById('ui');
  document.getElementById('status').textContent = 'Running';
  const fps = createFpsCounter(ui);

  // Build rope: anchor top at world (0, 5.5)
  const rope = WasmRopeJs.new(ANCHOR_X, ANCHOR_Y, 0.0, LINK_COUNT, SEG_LEN, 0.05);
  if (!rope) { document.getElementById('status').textContent = 'Rope init failed'; return; }
  rope.pin_head(ANCHOR_X, ANCHOR_Y, 0.0);

  // Damping slider — wired directly to Rust via set_damping().
  createSlider(ui, 'Damping', 0.90, 1.0, rope.get_damping(), v => rope.set_damping(v));

  // Mouse drag state
  let mouse = null;
  let prevMouse = null;
  canvas.addEventListener('mousedown', e => { mouse = screenToWorld(e); prevMouse = mouse; });
  canvas.addEventListener('mousemove', e => {
    if (mouse) {
      const cur = screenToWorld(e);
      const dx = cur[0] - prevMouse[0], dy = cur[1] - prevMouse[1];
      // Apply impulse to the link nearest the mouse
      const closest = findClosestLink(rope, cur);
      rope.apply_impulse(closest, dx * 8, dy * 8, 0.0);
      prevMouse = cur;
    }
  });
  canvas.addEventListener('mouseup',   () => { mouse = null; prevMouse = null; });
  canvas.addEventListener('mouseleave',() => { mouse = null; prevMouse = null; });

  // World-to-screen helpers (Y-up → Y-down)
  function worldToScreen(wx, wy) {
    const cx = canvas.width  / 2, cy = canvas.height / 2;
    const scale = canvas.height / 14;
    return [cx + wx * scale, cy - wy * scale];
  }
  function screenToWorld(e) {
    const rect = canvas.getBoundingClientRect();
    const px = e.clientX - rect.left, py = e.clientY - rect.top;
    const cx = canvas.width / 2, cy = canvas.height / 2;
    const scale = canvas.height / 14;
    return [(px - cx) / scale, -(py - cy) / scale];
  }
  function findClosestLink(r, pos) {
    const all = r.all_positions_flat();
    let best = 0, bestD = Infinity;
    for (let i = 0; i < r.link_count(); i++) {
      const dx = all[i * 3] - pos[0], dy = all[i * 3 + 1] - pos[1];
      const d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; best = i; }
    }
    return best;
  }

  function resize() {
    canvas.width  = window.innerWidth;
    canvas.height = window.innerHeight;
  }
  window.addEventListener('resize', resize);
  resize();

  let last = performance.now();
  function frame(now) {
    const dt = Math.min((now - last) / 1000, 0.05);
    last = now;

    // Step physics
    rope.step(dt);

    const positions = rope.all_positions_flat();
    const count = rope.link_count();

    // Draw
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // Background gradient
    const bg = ctx.createLinearGradient(0, 0, 0, canvas.height);
    bg.addColorStop(0, '#0d0d22');
    bg.addColorStop(1, '#0d1422');
    ctx.fillStyle = bg;
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    // Rope
    ctx.save();
    for (let i = 0; i < count - 1; i++) {
      const [ax, ay] = worldToScreen(positions[i * 3], positions[i * 3 + 1]);
      const [bx, by] = worldToScreen(positions[(i + 1) * 3], positions[(i + 1) * 3 + 1]);
      const t = i / (count - 1);
      const r = Math.round(80 + t * 175);
      const g = Math.round(200 - t * 100);
      const w = 3 + (1 - t) * 4;
      ctx.strokeStyle = `rgb(${r},${g},120)`;
      ctx.lineWidth = w;
      ctx.lineCap = 'round';
      ctx.beginPath();
      ctx.moveTo(ax, ay);
      ctx.lineTo(bx, by);
      ctx.stroke();
    }

    // Anchor dot
    const [ancX, ancY] = worldToScreen(ANCHOR_X, ANCHOR_Y);
    ctx.fillStyle = '#fff';
    ctx.beginPath();
    ctx.arc(ancX, ancY, 5, 0, Math.PI * 2);
    ctx.fill();

    ctx.restore();
    fps.update(dt);
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
}

main().catch(err => {
  document.getElementById('status').textContent = `Error: ${err.message}`;
  console.error(err);
});
