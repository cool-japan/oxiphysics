"""
Numpy bulk-array bridge tests for oxiphysics Python bindings.

These tests verify the three high-volume numpy array accessors added in
Block 6:
  - ValueNoise3D.sample_grid_to_numpy
  - DebugDrawSession.vertex_buffer_to_numpy
  - TelemetrySession.time_series_to_numpy

Run via::

    python -m pytest python/tests/test_numpy_bridge.py -v

after ``maturin develop --release``.
"""

import pytest

try:
    import oxiphysics
    import numpy as np

    HAS_OXIPHYSICS = True
except ImportError:
    HAS_OXIPHYSICS = False

pytestmark = pytest.mark.skipif(
    not HAS_OXIPHYSICS, reason="oxiphysics extension module not built"
)


# ---------------------------------------------------------------------------
# ValueNoise3D.sample_grid_to_numpy
# ---------------------------------------------------------------------------


def test_value_noise_to_numpy_shape_dtype():
    """sample_grid_to_numpy returns correct shape (nz, ny, nx) and dtype float64."""
    noise = oxiphysics.ValueNoise3D(seed=42)
    arr = noise.sample_grid_to_numpy(
        origin=(0.0, 0.0, 0.0),
        step=(0.1, 0.1, 0.1),
        shape=(4, 8, 16),
    )
    assert arr.shape == (4, 8, 16), f"Expected (4, 8, 16), got {arr.shape}"
    assert arr.dtype == np.float64, f"Expected float64, got {arr.dtype}"


def test_value_noise_to_numpy_values_in_range():
    """All sampled values are in [-1.0, 1.0] as documented."""
    noise = oxiphysics.ValueNoise3D(seed=0)
    arr = noise.sample_grid_to_numpy(
        origin=(0.0, 0.0, 0.0),
        step=(0.05, 0.05, 0.05),
        shape=(5, 5, 5),
    )
    assert arr.min() >= -1.0, f"Min value {arr.min()} < -1.0"
    assert arr.max() <= 1.0, f"Max value {arr.max()} > 1.0"


def test_value_noise_to_numpy_is_3d():
    """Array is exactly 3-dimensional."""
    noise = oxiphysics.ValueNoise3D(seed=7)
    arr = noise.sample_grid_to_numpy(
        origin=(1.0, 2.0, 3.0),
        step=(0.2, 0.2, 0.2),
        shape=(2, 3, 4),
    )
    assert arr.ndim == 3


def test_value_noise_to_numpy_consistent_with_sample():
    """Grid values match per-point sample() at corresponding coordinates."""
    noise = oxiphysics.ValueNoise3D(seed=99)
    origin = (0.0, 0.0, 0.0)
    step = (0.1, 0.1, 0.1)
    shape = (2, 2, 2)
    arr = noise.sample_grid_to_numpy(origin=origin, step=step, shape=shape)
    # Spot-check the [iz, iy, ix] = [1, 0, 1] cell
    iz, iy, ix = 1, 0, 1
    x = origin[0] + ix * step[0]
    y = origin[1] + iy * step[1]
    z = origin[2] + iz * step[2]
    expected = noise.sample(x, y, z)
    assert abs(arr[iz, iy, ix] - expected) < 1e-12


# ---------------------------------------------------------------------------
# DebugDrawSession.vertex_buffer_to_numpy
# ---------------------------------------------------------------------------


def test_debug_draw_vertex_buffer_to_numpy_shape():
    """vertex_buffer_to_numpy returns a (n, 3) float64 array."""
    session = oxiphysics.DebugDrawSession()
    session.begin_step(0)
    session.add_line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.0])
    session.add_sphere([0.5, 0.5, 0.5], 0.25, [0.0, 1.0, 0.0, 1.0])
    arr = session.vertex_buffer_to_numpy()
    assert arr.ndim == 2, f"Expected 2-D array, got ndim={arr.ndim}"
    assert arr.shape[1] == 3, f"Expected shape[-1]==3, got {arr.shape}"
    assert arr.dtype == np.float64, f"Expected float64, got {arr.dtype}"


def test_debug_draw_vertex_buffer_to_numpy_empty():
    """Empty session returns a (0, 3) array, not an error."""
    session = oxiphysics.DebugDrawSession()
    arr = session.vertex_buffer_to_numpy()
    assert arr.shape == (0, 3), f"Expected (0, 3), got {arr.shape}"
    assert arr.dtype == np.float64


def test_debug_draw_vertex_buffer_to_numpy_count():
    """One command → exactly one vertex row."""
    session = oxiphysics.DebugDrawSession()
    session.begin_step(0)
    session.add_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0])
    arr = session.vertex_buffer_to_numpy()
    assert arr.shape[0] == 1, f"Expected 1 vertex row, got {arr.shape[0]}"


def test_debug_draw_vertex_buffer_to_numpy_values():
    """Vertex coordinates match the positions passed to add_line."""
    session = oxiphysics.DebugDrawSession()
    session.begin_step(0)
    session.add_line([3.0, 4.0, 5.0], [6.0, 7.0, 8.0], [1.0, 1.0, 0.0, 1.0])
    arr = session.vertex_buffer_to_numpy()
    # The Line command reports its start vertex
    np.testing.assert_allclose(arr[0], [3.0, 4.0, 5.0])


# ---------------------------------------------------------------------------
# TelemetrySession.time_series_to_numpy
# ---------------------------------------------------------------------------


def _push_stats(session, n: int = 100):
    """Push n stats entries to a TelemetrySession."""
    for i in range(n):
        session.push(
            step=i,
            dt=0.016,
            body_count=10,
            sleeping_count=0,
            contact_count=i % 5,
            island_count=2,
            solve_iterations=8,
            broad_phase_pairs=20,
            kinetic_energy=float(i),
            elapsed_ms=1.0 + 0.01 * i,
        )


def test_telemetry_time_series_length_match():
    """time_series_to_numpy returns array with length == sample count."""
    session = oxiphysics.TelemetrySession(100)
    _push_stats(session, 100)
    arr = session.time_series_to_numpy("kinetic_energy")
    assert arr.shape == (100,), f"Expected (100,), got {arr.shape}"
    assert arr.dtype == np.float64


def test_telemetry_time_series_values_kinetic_energy():
    """Extracted kinetic_energy values match what was pushed."""
    session = oxiphysics.TelemetrySession(50)
    _push_stats(session, 50)
    arr = session.time_series_to_numpy("kinetic_energy")
    expected = np.arange(50, dtype=np.float64)
    np.testing.assert_allclose(arr, expected)


def test_telemetry_time_series_dt():
    """Extracting the 'dt' field returns the correct constant values."""
    session = oxiphysics.TelemetrySession(10)
    _push_stats(session, 10)
    arr = session.time_series_to_numpy("dt")
    assert arr.shape == (10,)
    np.testing.assert_allclose(arr, 0.016, atol=1e-12)


def test_telemetry_time_series_step():
    """Extracting 'step' returns monotonically increasing indices."""
    session = oxiphysics.TelemetrySession(20)
    _push_stats(session, 20)
    arr = session.time_series_to_numpy("step")
    expected = np.arange(20, dtype=np.float64)
    np.testing.assert_allclose(arr, expected)


def test_telemetry_time_series_unknown_field():
    """Unknown field name raises KeyError."""
    session = oxiphysics.TelemetrySession(10)
    _push_stats(session, 5)
    with pytest.raises(KeyError):
        session.time_series_to_numpy("nonexistent_field")


def test_telemetry_time_series_empty_session():
    """Empty session returns a zero-length array, not an error."""
    session = oxiphysics.TelemetrySession(100)
    arr = session.time_series_to_numpy("kinetic_energy")
    assert arr.shape == (0,)
    assert arr.dtype == np.float64


def test_telemetry_time_series_all_fields_accepted():
    """All documented field names are accepted without error."""
    session = oxiphysics.TelemetrySession(10)
    _push_stats(session, 3)
    valid_fields = [
        "step",
        "dt",
        "body_count",
        "sleeping_count",
        "contact_count",
        "island_count",
        "solve_iterations",
        "broad_phase_pairs",
        "kinetic_energy",
        "elapsed_ms",
    ]
    for field in valid_fields:
        arr = session.time_series_to_numpy(field)
        assert arr.shape == (3,), f"Field '{field}': expected shape (3,), got {arr.shape}"
        assert arr.dtype == np.float64, f"Field '{field}': expected float64"
