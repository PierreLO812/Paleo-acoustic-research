from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "models"
OUTPUT = ROOT / "output"
PROMPTS = ROOT / "prompts"
RUNS_STATE = ROOT / "pipeline" / ".runs"
