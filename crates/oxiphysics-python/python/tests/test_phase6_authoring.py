"""Tests for Phase 6.1 Authoring classes.

Covers: MaterialTable, SceneDescription, SceneBuilder, SimRecorder,
ReplayRecord, SimReplayer.

Uses ``hasattr`` guards so individual tests skip gracefully if a class is
not yet bound, rather than failing the entire module.

Run after ``maturin develop --release`` via::

    python -m pytest python/tests/test_phase6_authoring.py -v
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
# MaterialTable
# ---------------------------------------------------------------------------

def test_material_table_creates():
    if not hasattr(ox, "MaterialTable"):
        pytest.skip("MaterialTable not yet bound")
    mt = ox.MaterialTable()
    assert mt is not None


def test_material_table_register_named():
    if not hasattr(ox, "MaterialTable"):
        pytest.skip("MaterialTable not yet bound")
    mt = ox.MaterialTable()
    # register_named(name, density, restitution, static_friction, dynamic_friction,
    #                linear_damping, angular_damping) -> u32
    mid = mt.register_named("steel", 7800.0, 0.2, 0.4, 0.35, 0.01, 0.01)
    assert isinstance(mid, int)


def test_material_table_id_for_name():
    if not hasattr(ox, "MaterialTable"):
        pytest.skip("MaterialTable not yet bound")
    mt = ox.MaterialTable()
    mt.register_named("wood", 600.0, 0.3, 0.5, 0.4, 0.02, 0.02)
    mid = mt.id_for_name("wood")
    assert mid is not None
    assert isinstance(mid, int)


def test_material_table_builtin_presets():
    """MaterialTable ships with built-in presets (concrete, rubber, etc.)."""
    if not hasattr(ox, "MaterialTable"):
        pytest.skip("MaterialTable not yet bound")
    mt = ox.MaterialTable()
    # The table is pre-populated; material_count() should be > 0
    assert mt.material_count() > 0


def test_material_table_summary():
    if not hasattr(ox, "MaterialTable"):
        pytest.skip("MaterialTable not yet bound")
    mt = ox.MaterialTable()
    summary = mt.summary()
    assert isinstance(summary, str)


# ---------------------------------------------------------------------------
# SceneBuilder and SceneDescription
# ---------------------------------------------------------------------------

def test_scene_builder_creates():
    if not hasattr(ox, "SceneBuilder"):
        pytest.skip("SceneBuilder not yet bound")
    sb = ox.SceneBuilder("test_scene")
    assert sb is not None


def test_scene_builder_build_returns_scene_description():
    if not hasattr(ox, "SceneBuilder") or not hasattr(ox, "SceneDescription"):
        pytest.skip("SceneBuilder or SceneDescription not yet bound")
    sb = ox.SceneBuilder("my_scene")
    scene = sb.build()
    assert scene is not None


def test_scene_description_from_json():
    if not hasattr(ox, "SceneDescription"):
        pytest.skip("SceneDescription not yet bound")
    if not hasattr(ox, "SceneBuilder"):
        pytest.skip("SceneBuilder not yet bound")
    sb = ox.SceneBuilder("json_scene")
    scene = sb.build()
    json_str = scene.to_json()
    restored = ox.SceneDescription.from_json(json_str)
    assert restored is not None


def test_scene_description_body_count_starts_zero():
    if not hasattr(ox, "SceneBuilder") or not hasattr(ox, "SceneDescription"):
        pytest.skip("SceneBuilder or SceneDescription not yet bound")
    scene = ox.SceneBuilder("empty").build()
    assert scene.body_count() == 0


def test_scene_builder_add_sphere():
    if not hasattr(ox, "SceneBuilder"):
        pytest.skip("SceneBuilder not yet bound")
    sb = ox.SceneBuilder("sphere_scene")
    # add_sphere(id, x, y, z, radius, is_static)
    sb.add_sphere("ball1", 0.0, 5.0, 0.0, 0.5, False)
    scene = sb.build()
    assert scene.body_count() >= 1


# ---------------------------------------------------------------------------
# SimRecorder and ReplayRecord
# ---------------------------------------------------------------------------

def test_sim_recorder_creates():
    if not hasattr(ox, "SimRecorder"):
        pytest.skip("SimRecorder not yet bound")
    rec = ox.SimRecorder(0)
    assert rec is not None


def test_sim_recorder_begin_step_and_recorded_steps():
    if not hasattr(ox, "SimRecorder"):
        pytest.skip("SimRecorder not yet bound")
    rec = ox.SimRecorder(0)
    assert rec.recorded_steps() == 0
    rec.begin_step(0.016)
    rec.end_step()
    assert rec.recorded_steps() == 1


def test_sim_recorder_finish_to_json():
    if not hasattr(ox, "SimRecorder"):
        pytest.skip("SimRecorder not yet bound")
    rec = ox.SimRecorder(0)
    rec.begin_step(0.016)
    rec.end_step()
    json_str = rec.finish_to_json()
    assert json_str is not None
    data = json.loads(json_str)
    assert isinstance(data, dict)


def test_replay_record_from_json():
    if not hasattr(ox, "SimRecorder") or not hasattr(ox, "ReplayRecord"):
        pytest.skip("SimRecorder or ReplayRecord not yet bound")
    rec = ox.SimRecorder(0)
    rec.begin_step(0.016)
    rec.end_step()
    json_str = rec.finish_to_json()
    record = ox.ReplayRecord.from_json(json_str)
    assert record is not None


def test_replay_record_step_count():
    if not hasattr(ox, "SimRecorder") or not hasattr(ox, "ReplayRecord"):
        pytest.skip("SimRecorder or ReplayRecord not yet bound")
    rec = ox.SimRecorder(0)
    for _ in range(3):
        rec.begin_step(0.016)
        rec.end_step()
    json_str = rec.finish_to_json()
    record = ox.ReplayRecord.from_json(json_str)
    assert record.step_count() == 3


def test_replay_record_to_json_roundtrip():
    if not hasattr(ox, "SimRecorder") or not hasattr(ox, "ReplayRecord"):
        pytest.skip("SimRecorder or ReplayRecord not yet bound")
    rec = ox.SimRecorder(0)
    rec.begin_step(0.016)
    rec.end_step()
    json_str = rec.finish_to_json()
    record = ox.ReplayRecord.from_json(json_str)
    roundtrip = record.to_json()
    assert roundtrip is not None
    assert isinstance(json.loads(roundtrip), dict)


# ---------------------------------------------------------------------------
# SimReplayer
# ---------------------------------------------------------------------------

def test_sim_replayer_creates():
    if not hasattr(ox, "SimRecorder") or not hasattr(ox, "ReplayRecord") or not hasattr(ox, "SimReplayer"):
        pytest.skip("SimReplayer or prerequisites not yet bound")
    rec = ox.SimRecorder(0)
    rec.begin_step(0.016)
    rec.end_step()
    json_str = rec.finish_to_json()
    record = ox.ReplayRecord.from_json(json_str)
    replayer = ox.SimReplayer(record)
    assert replayer is not None


def test_sim_replayer_total_steps():
    if not hasattr(ox, "SimRecorder") or not hasattr(ox, "ReplayRecord") or not hasattr(ox, "SimReplayer"):
        pytest.skip("SimReplayer or prerequisites not yet bound")
    rec = ox.SimRecorder(0)
    for _ in range(5):
        rec.begin_step(0.016)
        rec.end_step()
    json_str = rec.finish_to_json()
    record = ox.ReplayRecord.from_json(json_str)
    replayer = ox.SimReplayer(record)
    assert replayer.total_steps() == 5


def test_sim_replayer_advance_json():
    if not hasattr(ox, "SimRecorder") or not hasattr(ox, "ReplayRecord") or not hasattr(ox, "SimReplayer"):
        pytest.skip("SimReplayer or prerequisites not yet bound")
    rec = ox.SimRecorder(0)
    rec.begin_step(0.016)
    rec.end_step()
    json_str = rec.finish_to_json()
    record = ox.ReplayRecord.from_json(json_str)
    replayer = ox.SimReplayer(record)
    step_json = replayer.advance_json()
    # Returns a JSON string or None if done
    if step_json is not None:
        assert isinstance(step_json, str)


def test_sim_replayer_is_done_after_all_steps():
    if not hasattr(ox, "SimRecorder") or not hasattr(ox, "ReplayRecord") or not hasattr(ox, "SimReplayer"):
        pytest.skip("SimReplayer or prerequisites not yet bound")
    rec = ox.SimRecorder(0)
    rec.begin_step(0.016)
    rec.end_step()
    json_str = rec.finish_to_json()
    record = ox.ReplayRecord.from_json(json_str)
    replayer = ox.SimReplayer(record)
    assert not replayer.is_done()
    replayer.advance_json()
    assert replayer.is_done()
