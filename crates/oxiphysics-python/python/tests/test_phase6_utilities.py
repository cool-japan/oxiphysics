"""Tests for Phase 6.2 Physics Utility classes.

Covers: DebugDrawSession, DrawList, BuoyancyWorld, Scheduler, SpatialGrid,
LodSystem, ValueNoise3D, FractalNoise, TelemetrySession, ContactCache.

Uses ``hasattr`` guards so individual tests skip gracefully if a class is
not yet bound.

Run after ``maturin develop --release`` via::

    python -m pytest python/tests/test_phase6_utilities.py -v
"""

import json
import pytest

try:
    import oxiphysics as ox
    HAS_OXIPHYSICS = True
except ImportError:
    HAS_OXIPHYSICS = False

pytestmark = pytest.mark.skipif(not HAS_OXIPHYSICS, reason="oxiphysics not built")


# ---------------------------------------------------------------------------
# DebugDrawSession
# ---------------------------------------------------------------------------

def test_debug_draw_session_creates():
    if not hasattr(ox, "DebugDrawSession"):
        pytest.skip("DebugDrawSession not yet bound")
    session = ox.DebugDrawSession()
    assert session is not None


def test_debug_draw_session_begin_step_and_step():
    if not hasattr(ox, "DebugDrawSession"):
        pytest.skip("DebugDrawSession not yet bound")
    session = ox.DebugDrawSession()
    session.begin_step(0)
    assert session.step() == 0


def test_debug_draw_session_add_line():
    if not hasattr(ox, "DebugDrawSession"):
        pytest.skip("DebugDrawSession not yet bound")
    session = ox.DebugDrawSession()
    session.begin_step(1)
    session.add_line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.0])


def test_debug_draw_session_add_sphere():
    if not hasattr(ox, "DebugDrawSession"):
        pytest.skip("DebugDrawSession not yet bound")
    session = ox.DebugDrawSession()
    session.begin_step(0)
    session.add_sphere([0.5, 0.5, 0.5], 0.25, [0.0, 1.0, 0.0, 1.0])


def test_debug_draw_session_end_step_json():
    if not hasattr(ox, "DebugDrawSession"):
        pytest.skip("DebugDrawSession not yet bound")
    session = ox.DebugDrawSession()
    session.begin_step(0)
    session.add_aabb([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0])
    result = session.end_step_json()
    assert result is not None
    assert isinstance(json.loads(result), (dict, list))


# ---------------------------------------------------------------------------
# DrawList
# ---------------------------------------------------------------------------

def test_draw_list_creates():
    if not hasattr(ox, "DrawList"):
        pytest.skip("DrawList not yet bound")
    dl = ox.DrawList()
    assert dl is not None


def test_draw_list_is_empty():
    if not hasattr(ox, "DrawList"):
        pytest.skip("DrawList not yet bound")
    dl = ox.DrawList()
    assert dl.is_empty()
    assert dl.len() == 0


def test_draw_list_add_line():
    if not hasattr(ox, "DrawList"):
        pytest.skip("DrawList not yet bound")
    dl = ox.DrawList()
    dl.add_line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.0])
    assert not dl.is_empty()
    assert dl.len() == 1


def test_draw_list_to_json():
    if not hasattr(ox, "DrawList"):
        pytest.skip("DrawList not yet bound")
    dl = ox.DrawList()
    dl.add_sphere([0.0, 0.0, 0.0], 1.0, [0.0, 0.0, 1.0, 1.0])
    result = dl.to_json()
    assert result is not None
    assert isinstance(json.loads(result), (dict, list))


def test_draw_list_clear():
    if not hasattr(ox, "DrawList"):
        pytest.skip("DrawList not yet bound")
    dl = ox.DrawList()
    dl.add_line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 1.0])
    dl.clear()
    assert dl.is_empty()


# ---------------------------------------------------------------------------
# BuoyancyWorld
# ---------------------------------------------------------------------------

def test_buoyancy_world_creates():
    if not hasattr(ox, "BuoyancyWorld"):
        pytest.skip("BuoyancyWorld not yet bound")
    world = ox.BuoyancyWorld(-9.81)
    assert world is not None


def test_buoyancy_world_add_fluid():
    if not hasattr(ox, "BuoyancyWorld"):
        pytest.skip("BuoyancyWorld not yet bound")
    world = ox.BuoyancyWorld(-9.81)
    # add_fluid(min [x,y,z], max [x,y,z], density, surface_y) -> u32
    fluid_id = world.add_fluid([-10.0, -2.0, -10.0], [10.0, 0.0, 10.0], 1000.0, 0.0)
    assert isinstance(fluid_id, int)


def test_buoyancy_world_fluid_count():
    if not hasattr(ox, "BuoyancyWorld"):
        pytest.skip("BuoyancyWorld not yet bound")
    world = ox.BuoyancyWorld(-9.81)
    assert world.fluid_count() == 0
    world.add_fluid([-10.0, -2.0, -10.0], [10.0, 0.0, 10.0], 1000.0, 0.0)
    assert world.fluid_count() == 1


def test_buoyancy_world_remove_fluid():
    if not hasattr(ox, "BuoyancyWorld"):
        pytest.skip("BuoyancyWorld not yet bound")
    world = ox.BuoyancyWorld(-9.81)
    fid = world.add_fluid([-10.0, -2.0, -10.0], [10.0, 0.0, 10.0], 1000.0, 0.0)
    world.remove_fluid(fid)
    assert world.fluid_count() == 0


# ---------------------------------------------------------------------------
# Scheduler
# ---------------------------------------------------------------------------

def test_scheduler_creates():
    if not hasattr(ox, "Scheduler"):
        pytest.skip("Scheduler not yet bound")
    sched = ox.Scheduler(4, 16)
    assert sched is not None


def test_scheduler_island_count_starts_zero():
    if not hasattr(ox, "Scheduler"):
        pytest.skip("Scheduler not yet bound")
    sched = ox.Scheduler(4, 16)
    assert sched.island_count() == 0


def test_scheduler_add_island():
    if not hasattr(ox, "Scheduler"):
        pytest.skip("Scheduler not yet bound")
    sched = ox.Scheduler(4, 16)
    iid = sched.add_island("Normal", 10)
    assert iid >= 0


def test_scheduler_schedule_json():
    if not hasattr(ox, "Scheduler"):
        pytest.skip("Scheduler not yet bound")
    sched = ox.Scheduler(4, 16)
    sched.add_island("Normal", 5)
    result = sched.schedule_json(1.0 / 60.0)
    assert result is not None
    assert isinstance(json.loads(result), (dict, list))


def test_scheduler_auto_sleep_step():
    if not hasattr(ox, "Scheduler"):
        pytest.skip("Scheduler not yet bound")
    sched = ox.Scheduler(4, 16)
    sched.add_island("Normal", 5)
    sleeping = sched.auto_sleep_step()
    assert isinstance(sleeping, list)


# ---------------------------------------------------------------------------
# SpatialGrid
# ---------------------------------------------------------------------------

def test_spatial_grid_creates():
    if not hasattr(ox, "SpatialGrid"):
        pytest.skip("SpatialGrid not yet bound")
    grid = ox.SpatialGrid(1.0)
    assert grid is not None


def test_spatial_grid_empty_initially():
    if not hasattr(ox, "SpatialGrid"):
        pytest.skip("SpatialGrid not yet bound")
    grid = ox.SpatialGrid(1.0)
    assert grid.is_empty()
    assert grid.len() == 0


def test_spatial_grid_insert_and_len():
    if not hasattr(ox, "SpatialGrid"):
        pytest.skip("SpatialGrid not yet bound")
    grid = ox.SpatialGrid(1.0)
    grid.insert(1, [0.0, 0.0, 0.0])
    assert grid.len() == 1
    assert grid.contains(1)


def test_spatial_grid_query_radius():
    if not hasattr(ox, "SpatialGrid"):
        pytest.skip("SpatialGrid not yet bound")
    grid = ox.SpatialGrid(1.0)
    grid.insert(1, [0.0, 0.0, 0.0])
    grid.insert(2, [0.5, 0.0, 0.0])
    grid.insert(3, [10.0, 0.0, 0.0])
    results = grid.query_radius([0.0, 0.0, 0.0], 1.0)
    assert isinstance(results, list)
    assert len(results) >= 2


def test_spatial_grid_remove():
    if not hasattr(ox, "SpatialGrid"):
        pytest.skip("SpatialGrid not yet bound")
    grid = ox.SpatialGrid(1.0)
    grid.insert(42, [1.0, 2.0, 3.0])
    grid.remove(42)
    assert not grid.contains(42)


def test_spatial_grid_k_nearest():
    if not hasattr(ox, "SpatialGrid"):
        pytest.skip("SpatialGrid not yet bound")
    grid = ox.SpatialGrid(1.0)
    for i in range(5):
        grid.insert(i, [float(i), 0.0, 0.0])
    results = grid.k_nearest([0.0, 0.0, 0.0], 3)
    assert isinstance(results, list)
    assert len(results) <= 3


# ---------------------------------------------------------------------------
# LodSystem
# ---------------------------------------------------------------------------

def test_lod_system_creates():
    if not hasattr(ox, "LodSystem"):
        pytest.skip("LodSystem not yet bound")
    lod = ox.LodSystem([0.0, 0.0, 0.0])
    assert lod is not None


def test_lod_system_add_body():
    if not hasattr(ox, "LodSystem"):
        pytest.skip("LodSystem not yet bound")
    lod = ox.LodSystem([0.0, 0.0, 0.0])
    bid = lod.add_body([1.0, 0.0, 0.0], 1.0)
    assert isinstance(bid, int)


def test_lod_system_body_count():
    if not hasattr(ox, "LodSystem"):
        pytest.skip("LodSystem not yet bound")
    lod = ox.LodSystem([0.0, 0.0, 0.0])
    assert lod.body_count() == 0
    lod.add_body([0.0, 0.0, 0.0], 1.0)
    assert lod.body_count() == 1


def test_lod_system_update_json():
    if not hasattr(ox, "LodSystem"):
        pytest.skip("LodSystem not yet bound")
    lod = ox.LodSystem([0.0, 0.0, 0.0])
    lod.add_body([1.0, 0.0, 0.0], 1.0)
    result = lod.update_json()
    assert result is not None
    assert isinstance(json.loads(result), (dict, list))


def test_lod_system_tier_counts_returns_tuple():
    if not hasattr(ox, "LodSystem"):
        pytest.skip("LodSystem not yet bound")
    lod = ox.LodSystem([0.0, 0.0, 0.0])
    counts = lod.tier_counts()
    assert isinstance(counts, tuple)
    assert len(counts) == 4


# ---------------------------------------------------------------------------
# ValueNoise3D
# ---------------------------------------------------------------------------

def test_value_noise_3d_creates():
    if not hasattr(ox, "ValueNoise3D"):
        pytest.skip("ValueNoise3D not yet bound")
    noise = ox.ValueNoise3D(42)
    assert noise is not None


def test_value_noise_3d_seed():
    if not hasattr(ox, "ValueNoise3D"):
        pytest.skip("ValueNoise3D not yet bound")
    noise = ox.ValueNoise3D(123)
    assert noise.seed() == 123


def test_value_noise_3d_sample_in_range():
    if not hasattr(ox, "ValueNoise3D"):
        pytest.skip("ValueNoise3D not yet bound")
    noise = ox.ValueNoise3D(0)
    for x in [0.0, 0.5, 1.0, 1.5]:
        val = noise.sample(x, 0.0, 0.0)
        assert -1.0 <= val <= 1.0, f"Value {val} at x={x} out of [-1,1]"


def test_value_noise_3d_same_seed_reproducible():
    if not hasattr(ox, "ValueNoise3D"):
        pytest.skip("ValueNoise3D not yet bound")
    n1 = ox.ValueNoise3D(7)
    n2 = ox.ValueNoise3D(7)
    assert n1.sample(1.0, 2.0, 3.0) == n2.sample(1.0, 2.0, 3.0)


# ---------------------------------------------------------------------------
# FractalNoise
# ---------------------------------------------------------------------------

def test_fractal_noise_creates():
    if not hasattr(ox, "FractalNoise"):
        pytest.skip("FractalNoise not yet bound")
    fn_ = ox.FractalNoise(42, 4)
    assert fn_ is not None


def test_fractal_noise_seed():
    if not hasattr(ox, "FractalNoise"):
        pytest.skip("FractalNoise not yet bound")
    fn_ = ox.FractalNoise(99, 4)
    assert fn_.seed() == 99


def test_fractal_noise_fbm():
    if not hasattr(ox, "FractalNoise"):
        pytest.skip("FractalNoise not yet bound")
    fn_ = ox.FractalNoise(0, 4)
    val = fn_.fbm(1.0, 2.0, 3.0)
    assert isinstance(val, float)


def test_fractal_noise_turbulence():
    if not hasattr(ox, "FractalNoise"):
        pytest.skip("FractalNoise not yet bound")
    fn_ = ox.FractalNoise(0, 4)
    val = fn_.turbulence(1.0, 2.0, 3.0)
    assert isinstance(val, float)


def test_fractal_noise_ridged():
    if not hasattr(ox, "FractalNoise"):
        pytest.skip("FractalNoise not yet bound")
    fn_ = ox.FractalNoise(0, 4)
    val = fn_.ridged(1.0, 2.0, 3.0)
    assert isinstance(val, float)


# ---------------------------------------------------------------------------
# TelemetrySession
# ---------------------------------------------------------------------------

def test_telemetry_session_creates():
    if not hasattr(ox, "TelemetrySession"):
        pytest.skip("TelemetrySession not yet bound")
    session = ox.TelemetrySession(100)
    assert session is not None


def test_telemetry_session_empty_initially():
    if not hasattr(ox, "TelemetrySession"):
        pytest.skip("TelemetrySession not yet bound")
    session = ox.TelemetrySession(50)
    assert session.is_empty()
    assert session.len() == 0


def test_telemetry_session_push_and_len():
    if not hasattr(ox, "TelemetrySession"):
        pytest.skip("TelemetrySession not yet bound")
    session = ox.TelemetrySession(10)
    session.push(
        step=0, dt=0.016, body_count=5, sleeping_count=0, contact_count=2,
        island_count=1, solve_iterations=4, broad_phase_pairs=10,
        kinetic_energy=42.0, elapsed_ms=1.0,
    )
    assert session.len() == 1


def test_telemetry_session_latest_json():
    if not hasattr(ox, "TelemetrySession"):
        pytest.skip("TelemetrySession not yet bound")
    session = ox.TelemetrySession(10)
    session.push(
        step=0, dt=0.016, body_count=5, sleeping_count=0, contact_count=2,
        island_count=1, solve_iterations=4, broad_phase_pairs=10,
        kinetic_energy=42.0, elapsed_ms=1.0,
    )
    result = session.latest_json()
    assert result is not None
    data = json.loads(result)
    assert isinstance(data, dict)


def test_telemetry_session_clear():
    if not hasattr(ox, "TelemetrySession"):
        pytest.skip("TelemetrySession not yet bound")
    session = ox.TelemetrySession(10)
    session.push(
        step=0, dt=0.016, body_count=5, sleeping_count=0, contact_count=2,
        island_count=1, solve_iterations=4, broad_phase_pairs=10,
        kinetic_energy=42.0, elapsed_ms=1.0,
    )
    session.clear()
    assert session.is_empty()


def test_telemetry_session_to_csv():
    if not hasattr(ox, "TelemetrySession"):
        pytest.skip("TelemetrySession not yet bound")
    session = ox.TelemetrySession(10)
    session.push(
        step=0, dt=0.016, body_count=5, sleeping_count=0, contact_count=2,
        island_count=1, solve_iterations=4, broad_phase_pairs=10,
        kinetic_energy=42.0, elapsed_ms=1.0,
    )
    csv = session.to_csv()
    assert isinstance(csv, str)
    assert "kinetic_energy" in csv or "step" in csv


# ---------------------------------------------------------------------------
# ContactCache
# ---------------------------------------------------------------------------

def test_contact_cache_creates():
    if not hasattr(ox, "ContactCache"):
        pytest.skip("ContactCache not yet bound")
    cache = ox.ContactCache(3)
    assert cache is not None


def test_contact_cache_empty_initially():
    if not hasattr(ox, "ContactCache"):
        pytest.skip("ContactCache not yet bound")
    cache = ox.ContactCache(3)
    assert cache.is_empty()
    assert cache.entry_count() == 0


def test_contact_cache_update_pair():
    if not hasattr(ox, "ContactCache"):
        pytest.skip("ContactCache not yet bound")
    cache = ox.ContactCache(3)
    cache.begin_step()
    # update_pair(a, b, points)
    # points: list of (pos_a [x,y,z], pos_b [x,y,z], normal [x,y,z], depth)
    points = [([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.01)]
    cache.update_pair(1, 2, points)
    assert not cache.is_empty()


def test_contact_cache_evict_stale():
    if not hasattr(ox, "ContactCache"):
        pytest.skip("ContactCache not yet bound")
    cache = ox.ContactCache(1)
    cache.begin_step()
    points = [([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.01)]
    cache.update_pair(1, 2, points)
    # Advance past TTL without updating
    cache.begin_step()
    cache.begin_step()
    evicted = cache.evict_stale()
    assert isinstance(evicted, int)


def test_contact_cache_clear():
    if not hasattr(ox, "ContactCache"):
        pytest.skip("ContactCache not yet bound")
    cache = ox.ContactCache(3)
    cache.begin_step()
    points = [([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.01)]
    cache.update_pair(1, 2, points)
    cache.clear()
    assert cache.is_empty()


def test_contact_cache_pairs_above_threshold_json():
    if not hasattr(ox, "ContactCache"):
        pytest.skip("ContactCache not yet bound")
    cache = ox.ContactCache(3)
    cache.begin_step()
    points = [([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.1)]
    cache.update_pair(1, 2, points)
    result = cache.pairs_above_impulse_threshold_json(0.05)
    assert result is not None
    assert isinstance(json.loads(result), (dict, list))
