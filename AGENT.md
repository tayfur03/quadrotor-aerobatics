# AGENT.md

## Purpose

This repository is a MATLAB quadrotor simulation and planning workspace. It combines:

- Differential-flatness trajectory generation and INDI-based control
- Terrain-aware and radar-aware path planning
- Radar threat-map and line-of-sight modeling
- Script demos for missions, terrain masking, planner benchmarks, FMM/FSM/BIT*, and AIT* experiments

Treat this as a script-driven MATLAB project rather than a packaged library.

## Repo Layout

- `demos/`: primary entry points and visual benchmark/demo scripts
- `tests/`: lightweight validation and smoke-test scripts
- `motion_planner/`: planners, corridor generation, path simplification, trajectory smoothing
- `terrain/`: DEM loading, terrain maps/meshes, LOS checks, ray casting
- `radar/`: radar models, threat maps, FMM/FSM heuristics, SkyMap visualization
- `flight/`: dynamics, INDI control, differential flatness, aerobatic helpers
- `math/`: quaternion and rotation helpers
- `safety/`: CBF safety-filter prototype
- `DEM/`: local terrain assets used by selected demos
- `project/`: shared path bootstrap
- `docs/`: repo notes and audits

## Primary Entry Points

Common scripts to run from the repo root in MATLAB:

- `run('demos/demo_indi_6dof_ned.m')`
- `run('demos/demo_mission_1.m')`
- `run('demos/demo_mission_2.m')`
- `run('demos/demo_mission_3.m')`
- `run('demos/demo_aerobatics.m')`
- `run('demos/demo_radar_binary_masking.m')`
- `run('demos/demo_radar_fmm_bit.m')`
- `run('demos/demo_fsm_bit_star_direct.m')`
- `run('demos/demo_ait_star_direct.m')`
- `run('demos/demo_radar_planner_benchmark.m')`
- `run('demos/final_mission_demo.m')`

Smaller validation scripts:

- `run('tests/test_trajectory_smoother.m')`
- `run('tests/test_mission2_traj.m')`
- `run('tests/test_ait_star_direct_smoke.m')`

Prefer the smallest relevant script for validation before running the heavier demos.

## Current Trajectory Flow

The active planner-to-trajectory flow is:

- planner path from `rrt_star_radar`, `fsm_bit_star_planner_direct`, `ait_star_planner_direct`, or related BIT* code
- optional geometric cleanup with `simplify_path` or `smooth_path_geometric`
- trajectory generation with `trajectory_smoother`
- optional downstream flight/control visualization

The old sliding-window trajectory generator has been removed. New work should use `trajectory_smoother` and `optimize_time_allocation`.

## MATLAB Conventions In This Repo

- Top-level scripts usually begin with `clear; clc; close all;`.
- Most scripts should call `setup_project_paths` or explicitly add the folders they use.
- Do not assume a fresh MATLAB session already has project paths loaded.
- Preserve script-first workflows unless intentionally refactoring the full call chain.

## Coordinate And Sign Conventions

- The project uses NED coordinates.
- Position is typically `[N; E; D]`.
- Several planner/demo flows use negative Z/D values for altitude above terrain.
- Terrain height, AGL, MSL, and planner altitude signs are easy to break silently.

Before changing planner or demo math, trace the sign convention end-to-end.

## Dependencies And Environment Assumptions

Expected MATLAB functionality includes:

- Core MATLAB scripting, plotting, interpolation, and numeric routines
- `minsnappolytraj` for minimum-snap trajectory generation
- Image Processing Toolbox for `bwdist` in corridor and binary-mask flows
- Statistics and Machine Learning Toolbox for `createns`, `knnsearch`, and `rangesearch`
- Parallel Computing Toolbox is optional; some threat-map and planner code can use `parfor` but should fall back to serial mode

Some demos depend on local DEM assets under `DEM/`. Preserve synthetic-terrain fallback behavior where it exists.

## High-Risk Areas

Use extra care when editing:

- `motion_planner/rrt_star_radar.m`
  - performance-sensitive
  - many interacting planner options
  - KD-tree helpers, rewiring logic, hard radar constraints, terrain queries
- `motion_planner/bit_star_planner.m`
  - depends on the `compute_stealth_corridor` output contract
  - has different behavior in weighted-risk vs strict-no-risk modes
- `motion_planner/trajectory_smoother.m`
  - depends on `minsnappolytraj` and `optimize_time_allocation`
  - numerical conditioning can degrade with too many waypoints or long segment times
- Demo scripts that convert between terrain height, AGL targets, and negative-down planner state

Small logic changes in these files can affect both correctness and runtime.

## Editing Guidance For Agents

- Keep changes narrowly scoped.
- Preserve NED and altitude sign conventions.
- Preserve graceful degradation when assets or toolboxes are unavailable.
- Avoid reintroducing removed experimental helpers unless a current demo or test needs them.
- Do not treat older generated notes as the sole source of truth for the latest radar/FMM/FSM/BIT*/AIT* flow.

## Validation Guidance

After planner or trajectory changes, validate with the smallest relevant script first:

- Smoothing/continuity changes: `run('tests/test_trajectory_smoother.m')`
- Mission trajectory changes: `run('tests/test_mission2_traj.m')`
- AIT* direct planner changes: `run('tests/test_ait_star_direct_smoke.m')`
- FSM/BIT* direct planner changes: `run('demos/demo_fsm_bit_star_direct.m')`
- Broader radar planner behavior: `run('demos/demo_radar_planner_benchmark.m')`

If a full demo is too expensive to run, document what was not verified and why.
