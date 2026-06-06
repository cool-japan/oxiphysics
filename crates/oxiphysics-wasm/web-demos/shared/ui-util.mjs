// Copyright (C) 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//
// ui-util.mjs — Lightweight HUD utility helpers.

/**
 * Create an FPS counter that updates a DOM element.
 *
 * @param {HTMLElement} container  Element to append the counter into.
 * @returns {{ update: function(number): void }}  Call `update(dt)` each frame.
 */
export function createFpsCounter(container) {
  const el = document.createElement('p');
  el.style.cssText = 'margin:0;font-size:12px;color:#8f8';
  el.textContent = 'FPS: --';
  container.appendChild(el);

  let frames = 0;
  let elapsed = 0;

  return {
    update(dt) {
      frames++;
      elapsed += dt;
      if (elapsed >= 0.5) {
        el.textContent = `FPS: ${Math.round(frames / elapsed)}`;
        frames = 0;
        elapsed = 0;
      }
    },
  };
}

/**
 * Create a labelled range slider.
 *
 * @param {HTMLElement} container   Parent element.
 * @param {string}      label       Text label displayed next to the slider.
 * @param {number}      min         Minimum value.
 * @param {number}      max         Maximum value.
 * @param {number}      value       Initial value.
 * @param {function(number): void} onChange  Called whenever the value changes.
 * @returns {{ getValue: function(): number }}
 */
export function createSlider(container, label, min, max, value, onChange) {
  const wrap = document.createElement('div');
  wrap.style.cssText = 'margin:2px 0;font-size:12px;color:#ccc;display:flex;gap:6px;align-items:center';

  const lbl = document.createElement('span');
  lbl.textContent = label;
  lbl.style.minWidth = '90px';

  const input = document.createElement('input');
  input.type = 'range';
  input.min = String(min);
  input.max = String(max);
  input.step = String((max - min) / 200);
  input.value = String(value);
  input.style.width = '120px';

  const val = document.createElement('span');
  val.textContent = value.toFixed(2);
  val.style.minWidth = '40px';

  input.addEventListener('input', () => {
    const v = parseFloat(input.value);
    val.textContent = v.toFixed(2);
    onChange(v);
  });

  wrap.appendChild(lbl);
  wrap.appendChild(input);
  wrap.appendChild(val);
  container.appendChild(wrap);

  return {
    getValue() {
      return parseFloat(input.value);
    },
  };
}
