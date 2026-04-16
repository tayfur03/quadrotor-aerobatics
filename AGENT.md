# AGENT.md

## Purpose

This repository is a MATLAB quadrotor simulation and planning workspace. It combines:

- Differential-flatness trajectory generation and INDI-based control
- Terrain-aware and radar-aware path planning
- Demo scripts for missions, aerobatics, terrain masking, and FMM/BIT* radar planning

Treat this as a script-driven MATLAB project rather than a packaged library.

## Repo Layout

- `demos/`: primary entry points and interactive benchmark/demo scripts
- `motion_planner/`: planners, corridor generation, smoothing, waypoint logic
- `terrain/`: DEM loading, terrain meshes/maps, LOS checks, ray casting
- `radar/`: radar models, threat maps, SkyMap visualization
- `DEM/`: local terrain assets used by some demos
- `backup_original/`: preserved older copies of original files
- `CLAUDE.md`: useful project notes, but partially stale relative to newer radar/FMM/BIT* work

## Primary Entry Points

Common scripts to run from the repo root in MATLAB:

- `run('demos/demo_indi_6dof_ned.m')`
- `run('demos/demo_mission_1.m')`
- `run('demos/demo_mission_2.m')`
- `run('demos/demo_mission_3.m')`
- `run('demos/demo_aerobatics.m')`
- `run('demos/demo_radar_fmm_bit.m')`
- `run('demos/demo_radar_planner_benchmark.m')`
- `run('demos/final_mission_demo.m')`

Smaller validation scripts:

- `run('tests/test_trajectory_smoother.m')`
- `run('tests/test_trajectory_upgrades.m')`

Prefer the smallest relevant script for validation before running the heavier demos.

## MATLAB Conventions In This Repo

- Top-level scripts usually begin with `clear; clc; close all;`.
- Scripts manage their own paths with:
  - `addpath('terrain')`
  - `addpath('radar')`
  - `addpath('motion_planner')`
- Do not assume the MATLAB session already has the right paths loaded.
- Preserve script-first workflows unless you are intentionally refactoring the full call chain.

## Coordinate And Sign Conventions

- The project uses NED coordinates.
- Position is typically `[N; E; D]`.
- Several planner/demo flows use negative Z to represent altitude in world/planner state.
- Terrain height, AGL, MSL, and planner altitude signs are easy to break silently.

Before changing planner or demo math, trace the sign convention end-to-end.

## Dependencies And Environment Assumptions

Expected MATLAB functionality includes:

- Core MATLAB scripting, plotting, and interpolation
- Image Processing Toolbox for `bwdist` in corridor generation and `demo_radar_fmm_bit`
- Statistics and Machine Learning Toolbox functions such as `createns`, `knnsearch`, and `rangesearch` in radar-aware planning
- Parallel Computing Toolbox is optional; some threat-map code can use `parfor` but falls back to serial mode

Some demos depend on local DEM assets under `DEM/`. Newer scripts often fall back to synthetic terrain if a DEM is missing. Preserve that behavior.

## High-Risk Areas

Use extra care when editing:

- `motion_planner/rrt_star_radar.m`
  - performance-sensitive
  - many interacting planner options
  - KD-tree helpers, rewiring logic, hard radar constraints, terrain queries
- `motion_planner/bit_star_planner.m`
  - depends on the `compute_stealth_corridor` output contract
  - has different behavior in weighted-risk vs strict-no-risk modes
- Demo scripts that convert between terrain height, AGL targets, and negative-down planner state

Small logic changes in these files can affect both correctness and runtime.

## Editing Guidance For Agents

- Keep changes narrowly scoped.
- Preserve NED and altitude sign conventions.
- Preserve explicit `addpath(...)` usage unless you update every affected script.
- Prefer graceful degradation when assets or toolboxes are unavailable.
- Do not treat `CLAUDE.md` as the sole source of truth for the latest radar/FMM/BIT* flow.

## Validation Guidance

After planner or trajectory changes, validate with the smallest relevant script first:

- Smoothing/continuity changes: `test_trajectory_smoother`
- Trajectory upgrade or planner-integration changes: `test_trajectory_upgrades`
- Radar/FMM/BIT* changes: `demo_radar_fmm_bit`
- End-to-end integrated behavior: `final_mission_demo`

If a full demo is too expensive to run, at minimum document what you could not verify and why.
