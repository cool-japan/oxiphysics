// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// canvas-util.mjs — WebGL2 canvas setup and rendering helpers.

/**
 * Acquire a WebGL2 context from a canvas element.
 * @param {string} id  Canvas element id.
 * @returns {{ canvas: HTMLCanvasElement, gl: WebGL2RenderingContext }}
 */
export function createCanvas(id) {
  const canvas = document.getElementById(id);
  if (!canvas) throw new Error(`Canvas element #${id} not found`);
  const gl = canvas.getContext('webgl2');
  if (!gl) throw new Error('WebGL2 is not available in this browser');
  gl.enable(gl.DEPTH_TEST);
  gl.enable(gl.CULL_FACE);
  return { canvas, gl };
}

/**
 * Compile a GLSL shader.
 * @param {WebGL2RenderingContext} gl
 * @param {number} type  gl.VERTEX_SHADER or gl.FRAGMENT_SHADER
 * @param {string} src
 * @returns {WebGLShader}
 */
function compileShader(gl, type, src) {
  const shader = gl.createShader(type);
  gl.shaderSource(shader, src);
  gl.compileShader(shader);
  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
    const info = gl.getShaderInfoLog(shader);
    gl.deleteShader(shader);
    throw new Error(`Shader compile error:\n${info}`);
  }
  return shader;
}

/**
 * Create and link a WebGL2 program from vertex and fragment GLSL sources.
 * @param {WebGL2RenderingContext} gl
 * @param {string} vsSource  Vertex shader source.
 * @param {string} fsSource  Fragment shader source.
 * @returns {WebGLProgram}
 */
export function createProgram(gl, vsSource, fsSource) {
  const vs = compileShader(gl, gl.VERTEX_SHADER, vsSource);
  const fs = compileShader(gl, gl.FRAGMENT_SHADER, fsSource);
  const prog = gl.createProgram();
  gl.attachShader(prog, vs);
  gl.attachShader(prog, fs);
  gl.linkProgram(prog);
  gl.deleteShader(vs);
  gl.deleteShader(fs);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) {
    const info = gl.getProgramInfoLog(prog);
    gl.deleteProgram(prog);
    throw new Error(`Program link error:\n${info}`);
  }
  return prog;
}

// ---------------------------------------------------------------------------
// Simple perspective camera with yaw / pitch orbit
// ---------------------------------------------------------------------------

/**
 * A minimal perspective camera supporting mouse-look and keyboard walk.
 *
 * Produces a column-major MVP matrix (Float32Array, 16 elements) suitable for
 * passing to a `uniform mat4` in GLSL.
 */
export class Camera {
  constructor({
    fovDeg = 60,
    near = 0.1,
    far = 1000,
    position = [0, 5, 10],
    yaw = 0,
    pitch = -0.3,
  } = {}) {
    this.fovDeg = fovDeg;
    this.near = near;
    this.far = far;
    this.position = Float64Array.from(position);
    this.yaw = yaw;
    this.pitch = pitch;
  }

  /** Compute the view direction unit vector. */
  get forward() {
    const cy = Math.cos(this.yaw), sy = Math.sin(this.yaw);
    const cp = Math.cos(this.pitch), sp = Math.sin(this.pitch);
    return [-sy * cp, sp, -cy * cp];
  }

  /** Move the camera in local space. */
  translate(dx, dy, dz) {
    const [fx, fy, fz] = this.forward;
    // right = forward × world_up
    const rx = fz, ry = 0, rz = -fx;
    this.position[0] += fx * dz + rx * dx;
    this.position[1] += fy * dz + dy;
    this.position[2] += fz * dz + rz * dx;
  }

  /**
   * Produce the combined projection × view matrix as a `Float32Array(16)`.
   * @param {number} aspect  Canvas width / height.
   */
  mvp(aspect) {
    const f = 1 / Math.tan((this.fovDeg * Math.PI) / 360);
    const rangeInv = 1 / (this.near - this.far);
    const proj = new Float64Array([
      f / aspect, 0, 0, 0,
      0, f, 0, 0,
      0, 0, (this.near + this.far) * rangeInv, -1,
      0, 0, this.near * this.far * rangeInv * 2, 0,
    ]);

    const [fx, fy, fz] = this.forward;
    // right = normalize(forward × up)
    const up = [0, 1, 0];
    let rx = fy * up[2] - fz * up[1];
    let ry = fz * up[0] - fx * up[2];
    let rz = fx * up[1] - fy * up[0];
    const rlen = Math.sqrt(rx * rx + ry * ry + rz * rz) || 1;
    rx /= rlen; ry /= rlen; rz /= rlen;
    // up' = right × forward
    const ux = ry * fz - rz * fy;
    const uy = rz * fx - rx * fz;
    const uz = rx * fy - ry * fx;
    const px = this.position[0], py = this.position[1], pz = this.position[2];
    const view = new Float64Array([
      rx, ux, -fx, 0,
      ry, uy, -fy, 0,
      rz, uz, -fz, 0,
      -(rx * px + ry * py + rz * pz),
      -(ux * px + uy * py + uz * pz),
        (fx * px + fy * py + fz * pz),
      1,
    ]);

    // Multiply proj × view (column-major)
    const out = new Float32Array(16);
    for (let col = 0; col < 4; col++) {
      for (let row = 0; row < 4; row++) {
        let sum = 0;
        for (let k = 0; k < 4; k++) {
          sum += proj[k * 4 + row] * view[col * 4 + k];
        }
        out[col * 4 + row] = sum;
      }
    }
    return out;
  }
}

// ---------------------------------------------------------------------------
// RAF render loop
// ---------------------------------------------------------------------------

/**
 * Start a requestAnimationFrame loop, calling `fn(dtSeconds)` each frame.
 * @param {function(number): void} fn
 */
export function startRenderLoop(fn) {
  let last = performance.now();
  function frame(now) {
    const dt = Math.min((now - last) / 1000, 0.1);
    last = now;
    fn(dt);
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
}
