# SigmaDrift

A biomechanically-grounded mouse movement algorithm that outperforms WindMouse across every metric that matters for human-like trajectory generation.

Built for my PhD research on novel mouse movement humanization techniques.

## What it does

SigmaDrift generates point-to-point mouse trajectories using six interacting components from computational motor control research:

- **Sigma-lognormal velocity primitives** — asymmetric bell-shaped speed profiles from Plamondon's Kinematic Theory
- **Two-phase surge architecture** — ballistic stroke (~93% of distance) followed by 0-2 corrective sub-movements
- **Ornstein-Uhlenbeck lateral drift** — mean-reverting stochastic hand drift
- **Signal-dependent noise** — motor noise scales with command magnitude (Harris-Wolpert), making Fitts' Law emerge naturally
- **Speed-modulated physiological tremor** — 8-12 Hz tremor suppressed during fast ballistic movement
- **Gamma-distributed inter-sample timing** — non-constant polling intervals matching real hardware behavior

## Results vs WindMouse

Same distance (~630px), same target width (20px):

| Metric | SigmaDrift | WindMouse | Real Human |
|---|---|---|---|
| Movement Time | 827 ms | 499 ms | ~750-850 ms |
| Fitts' Compliance | Yes (~8%) | No | Yes |
| Sub-Movements | 2 | 15 | 1-3 |
| Path Efficiency | 0.985 | 0.973 | 0.95-0.99 |
| Velocity Profile | Bell-shaped | Jagged | Bell-shaped |

## Usage

Header-only C++20, zero dependencies.

```cpp
#include "motor_synergy.h"

auto path = motor_synergy::generate(start_x, start_y, target_x, target_y);

for (auto& pt : path) {
    // pt.x, pt.y = position
    // pt.t = timestamp in ms
}
```

Custom configuration:

```cpp
motor_synergy::config cfg;
cfg.target_width = 16.0;
cfg.overshoot_prob = 0.20;

auto path = motor_synergy::generate(x0, y0, x1, y1, cfg);
auto m = motor_synergy::compute_metrics(path, x1, y1, cfg.target_width, dist);
```

## Harness

The included Win32 visualization harness (`main.cpp`) provides side-by-side comparison:

- **Space** — generate SigmaDrift trajectory (animated)
- **W** — generate WindMouse trajectory
- **R** — record your own mouse movement
- **S** — export trajectories to CSV
- **+/-** — adjust target width

Bottom panel shows velocity profile graphs for all trajectories overlaid.

## Building

Visual Studio 2022 with C++20. Open `SigmaDrift.slnx`, build x64 Release.

## References

- Plamondon — Kinematic Theory of Rapid Human Movements
- Flash & Hogan (1985) — Minimum-jerk model
- Harris & Wolpert (1998) — Signal-dependent noise
- Muller et al. (2017) — Control-theoretic models of pointing
- Acien et al. (2022) — BeCAPTCHA-Mouse
- Liu et al. (2024) — DMTG diffusion-based generation
