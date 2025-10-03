"""Generate demo barycentric samples from DE200 constants.

This script integrates a basic Newtonian N-body model seeded with the
initial state vectors and GM values packaged in the DE200 header. The
resulting samples are exported to ``data/de200_demo_positions.json`` for
use by the browser demo.

The integration is intentionally lightweight. For high-accuracy work use
a full ephemeris evaluation pipeline.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from math import sqrt
from pathlib import Path
from typing import Dict, List

import numpy as np
from jplephem.ascii import parse_header

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
HEADER_PATH = DATA_DIR / "header.200"
OUTPUT_PATH = DATA_DIR / "de200_demo_positions.json"

@dataclass
class BodyState:
    name: str
    position: np.ndarray
    velocity: np.ndarray
    gm: float


def load_constants() -> Dict[str, float]:
    if not HEADER_PATH.exists():
        raise FileNotFoundError(
            "Missing data/header.200. Run `python DEMO/fetch_ephemeris.py` first."
        )
    with HEADER_PATH.open() as fh:
        header = parse_header(fh)
    return dict(zip(header["names"], header["values"]))


def build_bodies(constants: Dict[str, float]) -> List[BodyState]:
    bodies = [
        ("Mercury", "1"),
        ("Venus", "2"),
        ("Earth-Moon Barycenter", "B"),
        ("Mars", "4"),
        ("Jupiter", "5"),
        ("Saturn", "6"),
        ("Uranus", "7"),
        ("Neptune", "8"),
        ("Pluto", "9"),
        ("Sun", "S"),
    ]
    gm_lookup = {
        "Mercury": constants["GM1"],
        "Venus": constants["GM2"],
        "Earth-Moon Barycenter": constants["GMB"],
        "Mars": constants["GM4"],
        "Jupiter": constants["GM5"],
        "Saturn": constants["GM6"],
        "Uranus": constants["GM7"],
        "Neptune": constants["GM8"],
        "Pluto": constants["GM9"],
        "Sun": constants["GMS"],
    }

    states: List[BodyState] = []
    for name, token in bodies:
        position = np.array(
            [constants[f"X{token}"], constants[f"Y{token}"], constants[f"Z{token}"]],
            dtype=float,
        )
        velocity = np.array(
            [constants[f"XD{token}"], constants[f"YD{token}"], constants[f"ZD{token}"]],
            dtype=float,
        )
        states.append(BodyState(name, position, velocity, gm_lookup[name]))
    return states


def integrate(states: List[BodyState], constants: Dict[str, float]) -> Dict[str, object]:
    n = len(states)
    pos = np.vstack([state.position for state in states])
    vel = np.vstack([state.velocity for state in states])
    mu = np.array([state.gm for state in states])

    step_days = 1.0
    total_steps = 720
    output_stride = 15
    samples: List[Dict[str, object]] = []

    def accelerations(positions: np.ndarray) -> np.ndarray:
        acc = np.zeros_like(positions)
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                delta = positions[j] - positions[i]
                r2 = float(delta.dot(delta))
                if r2 == 0.0:
                    continue
                acc[i] += mu[j] * delta / (r2 * sqrt(r2))
        return acc

    jd0 = constants["JDEPOC"]
    au_km = constants["AU"]

    current_pos = pos.copy()
    current_vel = vel.copy()
    current_acc = accelerations(current_pos)

    for step_index in range(total_steps + 1):
        if step_index % output_stride == 0:
            jd = jd0 + step_index * step_days
            samples.append(
                {
                    "julian_date": jd,
                    "positions_km": {
                        state.name: (current_pos[i] * au_km).tolist()
                        for i, state in enumerate(states)
                    },
                }
            )
        next_pos = (
            current_pos
            + current_vel * step_days
            + 0.5 * current_acc * step_days**2
        )
        next_acc = accelerations(next_pos)
        next_vel = current_vel + 0.5 * (current_acc + next_acc) * step_days
        current_pos, current_vel, current_acc = next_pos, next_vel, next_acc

    return {
        "metadata": {
            "description": "Newtonian integration seeded by DE200 constants",
            "start_julian_date": jd0,
            "step_days": step_days,
            "output_stride_days": output_stride,
            "au_km": au_km,
        },
        "bodies": [state.name for state in states],
        "samples": samples,
    }


def main() -> None:
    constants = load_constants()
    states = build_bodies(constants)
    dataset = integrate(states, constants)
    OUTPUT_PATH.write_text(json.dumps(dataset, indent=2))
    print(f"Wrote {OUTPUT_PATH} with {len(dataset['samples'])} frames")


if __name__ == "__main__":
    main()
