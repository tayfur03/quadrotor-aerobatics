function manifests = run_hard_budget_dataset_presets(preset)
%RUN_HARD_BUDGET_DATASET_PRESETS Run hard-only dataset generation presets.
%
% Usage:
%   run_hard_budget_dataset_presets("aggressive")
%   run_hard_budget_dataset_presets("more_aggressive")
%   run_hard_budget_dataset_presets("both")

if nargin < 1 || isempty(preset)
    preset = "aggressive";
end
preset = lower(string(preset));

switch preset
    case "aggressive"
        manifests = generate_radar_uav_dataset(make_hard_budget_aggressive_cfg());
    case "more_aggressive"
        manifests = generate_radar_uav_dataset(make_hard_budget_more_aggressive_cfg());
    case "both"
        manifests = struct();
        manifests.aggressive = generate_radar_uav_dataset(make_hard_budget_aggressive_cfg());
        manifests.more_aggressive = generate_radar_uav_dataset(make_hard_budget_more_aggressive_cfg());
    otherwise
        error('run_hard_budget_dataset_presets:UnknownPreset', ...
            'Unknown preset "%s". Use "aggressive", "more_aggressive", or "both".', preset);
end
end

function cfg = make_common_hard_budget_cfg(output_name)
cfg = struct();
cfg.output_dir = fullfile(pwd, output_name);

cfg.total_maps = 20;
cfg.K = 4;

cfg.scenario_sampling_mode = 'radar_relevant';
cfg.scenario_mix = [0.0, 0.0, 1.0]; % hard only
cfg.enforce_scenario_quota = true;

cfg.hard_straight_visible_range = [0.50, 1.00];
cfg.max_planned_visible_ratio = 0.55;
cfg.max_planned_mean_radar_cost = 0.45;
cfg.min_radar_improvement_ratio = 0.10;

cfg.enforce_visibility_acceptance = true;
cfg.visibility_margin_m = 0.1;
cfg.max_planned_visibility_violations = [16, 20, 24];

cfg.planner_max_time = 4;
cfg.planner_max_batches = 18;
cfg.planner_batch_size = 160;
cfg.planner_max_nodes = 3500;
cfg.planner_post_solution_batches = 2;
cfg.trajectory_batch_size = 8;
cfg.trajectory_attempt_oversample_factor = 8;

cfg.use_parallel = true;
cfg.use_parallel_cost_field = true;
cfg.use_parallel_trajectories = true;
cfg.use_gpu = true;
cfg.use_gpu_cost_field = true;
cfg.gpu_device_index = 1;
cfg.gpu_gather_each_layer = true;
cfg.use_parallel_ait_reverse = false;
cfg.use_parallel_ait_edge_checks = false;
cfg.forward_expand_chunk = 32;
cfg.graph_rebuild_interval = 3;
cfg.reverse_full_rebuild_interval = 3;

cfg.skip_existing = false;
cfg.write_partial_maps = true;
cfg.checkpoint_manifest = true;
end

function cfg = make_hard_budget_aggressive_cfg()
cfg = make_common_hard_budget_cfg('exports_hard_budget_aggressive');

cfg.scenario_planner_time_scale = [1.0, 1.8, 4.0];
cfg.scenario_planner_batch_scale = [1.0, 1.6, 2.6];
cfg.scenario_planner_node_scale = [1.0, 2.0, 4.0];
cfg.scenario_planner_post_solution_scale = [1.0, 1.5, 2.5];

cfg.max_path_attempts = max(80 * cfg.K, 1000);
cfg.radar_relevant_sampling_attempts = 1000;
end

function cfg = make_hard_budget_more_aggressive_cfg()
cfg = make_common_hard_budget_cfg('exports_hard_budget_more_aggressive');

cfg.scenario_planner_time_scale = [1.0, 2.5, 5.0];
cfg.scenario_planner_batch_scale = [1.0, 2.2, 4.0];
cfg.scenario_planner_node_scale = [1.0, 3.0, 6.0];
cfg.scenario_planner_post_solution_scale = [1.0, 2.0, 4.0];

cfg.max_path_attempts = max(120 * cfg.K, 1500);
cfg.radar_relevant_sampling_attempts = 1500;
end
