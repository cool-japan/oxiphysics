# Copyright 2026 COOLJAPAN OU (Team KitaSan)
# SPDX-License-Identifier: Apache-2.0

"""test_pyi_consistency.py — stub/runtime consistency checks.

For each public class declared in oxiphysics/__init__.pyi:
  1. Assert the class exists in the runtime `oxiphysics` module (or skip if
     the extension is not installed).
  2. Assert every method declared in the stub is present on the runtime class.

These tests do NOT import oxiphysics directly at module level so that the
entire file can be collected and the tests skipped gracefully when the
extension is not built (e.g. in a pure stub-checking CI step).
"""

from __future__ import annotations

import ast
import importlib
import inspect
import pathlib
import types
from typing import Dict, List, Set

import pytest

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_STUB_PATH = pathlib.Path(__file__).parent.parent / "oxiphysics" / "__init__.pyi"


def _parse_stub_classes() -> Dict[str, List[str]]:
    """Return {class_name: [method_names]} from the .pyi stub file."""
    source = _STUB_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    result: Dict[str, List[str]] = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef):
            methods: List[str] = []
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    methods.append(item.name)
            result[node.name] = methods
    return result


def _parse_stub_functions() -> List[str]:
    """Return top-level function names declared in the .pyi stub."""
    source = _STUB_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    return [
        node.name
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.col_offset == 0
    ]


# Parse once at collection time — the .pyi file must always exist.
_STUB_CLASSES: Dict[str, List[str]] = _parse_stub_classes()
_STUB_FUNCTIONS: List[str] = _parse_stub_functions()

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def oxiphysics_module():
    """Import the compiled oxiphysics extension, skip if not available."""
    try:
        mod = importlib.import_module("oxiphysics")
    except ImportError:
        pytest.skip("oxiphysics extension not installed — stub syntax-only run")
    return mod


# ---------------------------------------------------------------------------
# Test: stub classes exist in runtime module
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("class_name", sorted(_STUB_CLASSES.keys()))
def test_stub_classes_exist_in_runtime(oxiphysics_module: types.ModuleType, class_name: str) -> None:
    """Every class declared in the .pyi stub must exist in the runtime module."""
    assert hasattr(oxiphysics_module, class_name), (
        f"Class '{class_name}' declared in __init__.pyi "
        f"but not found in the oxiphysics runtime module."
    )
    runtime_obj = getattr(oxiphysics_module, class_name)
    assert isinstance(runtime_obj, type), (
        f"'{class_name}' exists in the module but is not a class (got {type(runtime_obj)!r})."
    )


# ---------------------------------------------------------------------------
# Test: stub methods exist on runtime classes
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "class_name,method_name",
    [
        (cls, meth)
        for cls, methods in sorted(_STUB_CLASSES.items())
        for meth in methods
    ],
)
def test_stub_methods_exist_on_classes(
    oxiphysics_module: types.ModuleType,
    class_name: str,
    method_name: str,
) -> None:
    """Every method declared in the .pyi stub must exist on the runtime class."""
    if not hasattr(oxiphysics_module, class_name):
        pytest.skip(f"Class '{class_name}' not in runtime — skipping method check.")

    runtime_cls = getattr(oxiphysics_module, class_name)

    # Accept dunder __init__ always (PyO3 always has it).
    if method_name == "__init__":
        return

    assert hasattr(runtime_cls, method_name), (
        f"{class_name}.{method_name}() declared in __init__.pyi "
        f"but not found on the runtime class."
    )


# ---------------------------------------------------------------------------
# Test: module-level functions exist in runtime module
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("func_name", sorted(set(_STUB_FUNCTIONS)))
def test_stub_module_functions_exist_in_runtime(
    oxiphysics_module: types.ModuleType,
    func_name: str,
) -> None:
    """Every module-level function in the stub must exist in the runtime module."""
    assert hasattr(oxiphysics_module, func_name), (
        f"Function '{func_name}' declared in __init__.pyi "
        f"but not found in the oxiphysics runtime module."
    )


# ---------------------------------------------------------------------------
# Test: .pyi stub file parses without errors (syntax check, no runtime needed)
# ---------------------------------------------------------------------------


def test_stub_file_valid_python_syntax() -> None:
    """The .pyi stub file must be valid Python (parseable by ast.parse)."""
    source = _STUB_PATH.read_text(encoding="utf-8")
    try:
        ast.parse(source)
    except SyntaxError as exc:
        pytest.fail(f"Syntax error in __init__.pyi: {exc}")


# ---------------------------------------------------------------------------
# Test: stub class set is non-trivially populated
# ---------------------------------------------------------------------------


def test_stub_covers_phase6_classes() -> None:
    """A representative sample of Phase-6 class names must appear in the stub."""
    required: Set[str] = {
        # Phase 6.1
        "AnimationPlayer",
        "Vec3Track",
        "QuatTrack",
        "MaterialTable",
        "SceneDescription",
        "SceneBuilder",
        "SimRecorder",
        "ReplayRecord",
        "SimReplayer",
        "BodySnapshot",
        "WorldSnapshot",
        "SnapshotManager",
        # Phase 6.2
        "DrawList",
        "DebugDrawSession",
        "ContactCache",
        "BuoyancyWorld",
        "Scheduler",
        "SpatialGrid",
        "LodSystem",
        "ValueNoise3D",
        "FractalNoise",
        "SpringFollower",
        "SpringFollower3",
        "TelemetrySession",
        # Phase 6.3
        "CharacterController",
        "Rope",
        "IkChain",
        "XpbdSolver",
        "AeroSystem",
        "NavMesh",
        "RollbackBuffer",
        "EventBus",
        "ProfilerSession",
        "TriggerWorld",
        "ForceFieldAabbRegion",
        "ForceFieldSystem",
        "QueryWorld",
    }
    stub_names = set(_STUB_CLASSES.keys())
    missing = required - stub_names
    assert not missing, (
        f"The following Phase-6 classes are absent from the stub: {sorted(missing)}"
    )
