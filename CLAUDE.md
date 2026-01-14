# Quadrotor Simulation - CLAUDE.md

## Project Overview

MATLAB-based quadrotor simulation implementing differential flatness-based control from **Tal & Karaman (2021)**: "Minimum Snap Trajectory Generation and Control for Quadrotors". The simulation features a 6-DOF quadrotor with cascaded control architecture using incremental nonlinear dynamic inversion (INDI).

## Directory Structure

```
test_rig/
├── demo_indi_6dof_ned.m      # Main simulation entry point (original trajectory tracking)
├── demo_mission_1.m          # Point-to-point navigation (no obstacles)
├── demo_mission_2.m          # Static obstacle avoidance with RRT*
├── demo_mission_3.m          # Pop-up threat avoidance with replanning
├── demo_aerobatics.m         # [NEW] Aerobatic maneuvers demonstration
├── aerobatic_maneuvers.m     # [NEW] Generate waypoints for aerobatic maneuvers
├── flat_outputs_demo.m        # Trajectory generation (heart/eight/roulette curves)
├── yaw_planner.m              # [UPDATED] Yaw planning with multiple modes
├── attitude_pd.m              # Quaternion PD attitude controller
├── flatness.m                 # Differential flatness feedforward (Eq. 14-15)
├── motor_inversion.m          # Motor command mapping with saturation
├── quad_params_indi.m         # Quadrotor physical/control parameters
├── quat_*.m (5 files)         # Quaternion utilities (mul, conj, normalize, to_R, R_to)
├── motion_planner/            # Motion planning module
│   ├── collision_checker.m    # Real-time collision detection
│   ├── obstacle_manager.m     # Static/dynamic obstacle handling (class)
│   ├── rrt_planner.m          # [UPDATED] RRT and RRT* path planning
│   ├── trajectory_smoother.m  # Uses MATLAB minsnappolytraj
│   └── waypoint_manager.m     # Mission sequencing and replanning (class)
├── papers/                    # Reference papers
│   ├── Tal et al. - 2023 - Aerobatic Trajectory Generation...
│   └── (other Tal & Karaman papers)
├── backup_original/           # Backup of original codebase
└── CLAUDE.md                  # This file
```

## Coordinate System

**NED (North-East-Down) Frame:**
- Position: `x = [N; E; D]` (meters)
- Velocity: `v = [vN; vE; vD]` (m/s)
- Attitude: quaternion `q = [qw; qx; qy; qz]` (body→inertial)
- Angular velocity: `omega = [p; q; r]` (body frame, rad/s)

## Control Architecture

```
Reference Trajectory → Position Control → Attitude Command → Attitude Control → Motor Inversion → Dynamics
     (x,v,a,j,s)         (Outer Loop)       (Incremental)      (Inner Loop)       (Actuator)      (6-DOF)
```

### Key Control Gains (from Table II of paper)
- Position: `Kx = diag([18, 18, 13.5])`
- Velocity: `Kv = diag([7.8, 7.8, 5.9])`
- Attitude: `K_xi = diag([175, 175, 82])`

## Simulation Parameters

- **Timestep:** dt = 0.002s (500 Hz)
- **Default Duration:** Tf = 16.0s
- **Mass:** 1.0 kg
- **Motor Bandwidth:** 50 Hz (τ_m = 0.02s)
- **Accelerometer Filter:** 30 Hz cutoff

## Current Trajectory System

The current system uses `flat_outputs_demo.m` to generate pre-defined smooth curves:
- **"heart"**: Parametric heart shape
- **"eight"**: Figure-8 lemniscate
- **"roulette"**: Complex multi-harmonic curve

All curves compute derivatives up to **snap** (4th derivative) for differential flatness.

---

## Motion Planner Module (NEW)

### Purpose
Add autonomous waypoint navigation and obstacle avoidance to the simulation, enabling:
1. Point-to-point navigation without obstacles
2. Navigation around known obstacles
3. Reactive avoidance of pop-up threats during flight

### Mission Types

| Mission | Description | Planner Mode |
|---------|-------------|--------------|
| **Task 1** | Go to target point, no obstacles | Direct waypoint with smooth trajectory |
| **Task 2** | Navigate around static obstacles | RRT* path planning (optimal) |
| **Task 3** | Avoid pop-up threats mid-flight | RRT replanning (fast) |

### RRT vs RRT* Algorithm Selection

The `rrt_planner.m` supports both algorithms via the `algorithm` parameter:

```matlab
% In waypoint_manager or directly:
params.algorithm = 'rrt';      % Fast, first valid path
params.algorithm = 'rrt_star'; % Optimal path with rewiring
params.rewire_radius = 1.5;    % RRT* neighborhood radius
```

- **RRT**: Fast planning, suitable for reactive replanning (demo_mission_3)
- **RRT***: Near-optimal paths, suitable for pre-planned routes (demo_mission_2)

### Key Components

#### 1. `rrt_planner.m`
- Rapidly-exploring Random Tree algorithm
- Inputs: start position, goal position, obstacle list, bounds
- Outputs: waypoint path avoiding obstacles
- Parameters: step size, max iterations, goal bias

#### 2. `waypoint_manager.m`
- Sequences waypoints for mission execution
- Tracks current waypoint and progress
- Triggers replanning when needed

#### 3. `obstacle_manager.m`
- Maintains list of obstacles (spheres/boxes)
- Supports static obstacles and dynamic "pop-up" threats
- API: `add_obstacle()`, `remove_obstacle()`, `get_obstacles()`

#### 4. `trajectory_smoother.m`
- Converts discrete waypoints to smooth trajectory
- Generates minimum-snap polynomials between waypoints
- Ensures C³ continuity for differential flatness

#### 5. `collision_checker.m`
- Real-time collision detection during simulation
- Returns true/false and distance to nearest obstacle
- Used for triggering replanning

### Integration Points

The motion planner integrates with existing code at:
1. **`reference_flat_outputs()`**: Replace fixed curves with planner-generated trajectories
2. **Main loop in `demo_indi_6dof_ned.m`**: Add collision checking and replanning hooks
3. **New demo scripts**: `demo_mission_*.m` for each task type

### Usage Example

```matlab
% Run demo scripts directly from MATLAB command window:
demo_mission_1   % Point-to-point: [0,0,0] -> [10,5,-3]
demo_mission_2   % Static obstacles: RRT planning around 3 spheres
demo_mission_3   % Pop-up threat: Replanning mid-flight at t=4s

% Or use the motion planner components directly:
addpath('motion_planner');

% Create waypoint manager
wm = waypoint_manager();
wm.set_mission([0;0;-2], [10;5;-3], obstacles, bounds);
[success, path] = wm.plan_path();
traj = wm.generate_trajectory(12.0);

% Get reference at time t
ref = wm.get_reference(t);  % Returns pos, vel, acc, jerk, snap

% For pop-up threats, use obstacle_manager
obs_mgr = obstacle_manager();
obs_mgr.add_popup_threat('sphere', [5;0;-2], struct('radius', 1.5), 4.0);
[new_threat, ~] = obs_mgr.update(current_time);
if new_threat
    [needs_replan, ~] = wm.check_replan_needed(current_pos, obs_mgr.get_active_obstacles());
    if needs_replan
        wm.replan(current_pos, current_time, obs_mgr.get_active_obstacles());
    end
end
```

---

## Aerobatic Maneuvers Module (NEW)

Based on **Tal et al. (2023)**: "Aerobatic Trajectory Generation for a VTOL Fixed-Wing Aircraft Using Differential Flatness".

### Yaw Planner Modes

The `yaw_planner.m` supports multiple yaw reference strategies:

| Mode | Description | Use Case |
|------|-------------|----------|
| `constant` | Fixed yaw angle | Hovering, simple flight |
| `constant_rate` | Linear yaw change (ψ = ψ₀ + ω·t) | Barrel rolls, rolling maneuvers |
| `tangent` | Align with horizontal velocity | 2D coordinated flight |
| `coordinated` | Full 3D velocity alignment | General aerobatics |
| `knife_edge` | Body Y vertical (ψ ± 90° from velocity) | Knife-edge flight |
| `poi` | Track a point of interest | Cinematic/inspection |
| `velocity_aligned` | Body X along 3D velocity | High-speed flight |

### Aerobatic Maneuvers

The `aerobatic_maneuvers.m` function generates waypoints compatible with `minsnappolytraj`:

| Maneuver | Description | Recommended Yaw Mode |
|----------|-------------|---------------------|
| `vertical_loop` | Loop in N-D plane | `coordinated` |
| `barrel_roll` | Helical roll with yaw rotation | `constant_rate` |
| `immelmann` | Half loop + half roll | `coordinated` |
| `split_s` | Half roll + half loop | `coordinated` |
| `climbing_turn` | 270° banked turn with climb | `coordinated` |
| `knife_edge` | Transition to knife-edge flight | `knife_edge` |
| `hover_flip` | Aggressive flip from hover | `constant` |
| `figure_eight` | Horizontal figure-8 | `tangent` |

### Usage Example

```matlab
% Generate barrel roll maneuver
params = struct('center', [0;0;-5], 'radius', 3, 'height_gain', 6, 'V', 8);
[waypoints, timePoints, velBC, accBC, info] = aerobatic_maneuvers('barrel_roll', params);

% Create smooth trajectory with minsnappolytraj
numSamples = 500;
[xref, vref, aref, jref, sref] = minsnappolytraj(waypoints, timePoints, numSamples, ...
    'VelocityBoundaryCondition', velBC, ...
    'AccelerationBoundaryCondition', accBC);

% Use recommended yaw mode
[psi, psi_dot, psi_ddot] = yaw_planner(t, info.yaw_mode, state_ref, info.yaw_params);
```

---

## Development Notes

### Running the Simulation
```matlab
% In MATLAB, from test_rig directory:
demo_indi_6dof_ned   % Original trajectory tracking demo
demo_mission_1       % Point-to-point mission
demo_mission_2       % Static obstacle avoidance (RRT*)
demo_mission_3       % Pop-up threat avoidance (RRT)
demo_aerobatics      % Aerobatic maneuvers demonstration
```

### Backup & Recovery
Original files backed up in `backup_original/`. To restore:
```matlab
copyfile('backup_original/*.m', '.')
```

### Key Equations Reference
- **Position Control:** Eq. 17 (Tal & Karaman)
- **Attitude Increment:** Eq. 22-26
- **Angular Velocity FF:** Eq. 14
- **Angular Acceleration FF:** Eq. 15

### Testing Checklist
- [ ] Verify original demos still work after changes
- [ ] Task 1: Drone reaches target within tolerance
- [ ] Task 2: Path avoids all static obstacles (RRT*)
- [ ] Task 3: Drone replans and avoids pop-up threat (RRT)
- [ ] Trajectory smoothness maintained (no discontinuities)
- [ ] Control effort within motor limits
- [ ] Aerobatic maneuvers track correctly
- [ ] Yaw modes produce expected behavior
