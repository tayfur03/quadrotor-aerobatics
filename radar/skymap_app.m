classdef skymap_app < handle
%SKYMAP_APP Interactive uifigure app for sliding binary threat visualization.
%   Panels:
%     - Top-down binary threat map (sliding window)
%     - Vertical binary slice along selected heading mode
%   Class colors:
%     1 -> below terrain, 2 -> hidden/safe, 3 -> visible

    properties (Access = private)
        cfg
        tm
        threat
        path
        start_pos
        goal_pos
        visibility_threshold
        n_wp
        current_idx

        class_cmap
        terrain_h_grid
        slice_alt_vec

        UIFigure
        TopAxes
        VertAxes

        WPSlider
        WPLabel
        PlayButton
        ModeDropDown
        WindowSpinner
        SliceSpinner

        TopImg
        TopCurr
        TopHeading
        VertImg
        VertTerrain
        VertCurr

        PlayTimer
        IsPlaying = false
    end

    methods
        function obj = skymap_app(cfg, tm, threat, path_rrt, start_pos, goal_pos, visibility_threshold)
            if nargin < 7 || isempty(visibility_threshold)
                visibility_threshold = 0.5;
            end
            if nargin < 1 || isempty(cfg)
                cfg = struct();
            end

            obj.cfg = obj.apply_defaults(cfg);
            obj.tm = tm;
            obj.threat = threat;
            obj.path = path_rrt;
            obj.start_pos = start_pos;
            obj.goal_pos = goal_pos;
            obj.visibility_threshold = visibility_threshold;
            obj.n_wp = size(path_rrt, 2);
            obj.current_idx = 1;
            obj.class_cmap = [0.35 0.35 0.35; 0.10 0.60 0.10; 0.80 0.10 0.10];

            if obj.n_wp < 2
                warning('skymap_app:shortPath', 'Path has fewer than 2 points. App not created.');
                return;
            end

            [N_grid, E_grid] = meshgrid(threat.N_vec, threat.E_vec);
            terrain_h = tm.get_height(N_grid(:), E_grid(:));
            obj.terrain_h_grid = reshape(terrain_h, size(N_grid));
            obj.slice_alt_vec = linspace(threat.alt_vec(1), threat.alt_vec(end), obj.cfg.skymap_slice_n_alt);

            obj.build_ui();
            obj.update_view(1);
        end

        function delete(obj)
            obj.stop_timer();
        end
    end

    methods (Access = private)
        function cfg = apply_defaults(~, cfg)
            cfg.skymap_half_window = get_cfg(cfg, 'skymap_half_window', 1800);
            cfg.skymap_slice_half_length = get_cfg(cfg, 'skymap_slice_half_length', 2500);
            cfg.skymap_slice_n_horiz = get_cfg(cfg, 'skymap_slice_n_horiz', 181);
            cfg.skymap_slice_n_alt = get_cfg(cfg, 'skymap_slice_n_alt', 120);
            cfg.skymap_play_period = get_cfg(cfg, 'skymap_play_period', 0.3);
            cfg.skymap_heading_mode = get_cfg(cfg, 'skymap_heading_mode', 'Path Heading');
            cfg.skymap_fixed_heading_deg = get_cfg(cfg, 'skymap_fixed_heading_deg', 0);
        end

        function build_ui(obj)
            obj.UIFigure = uifigure('Name', 'SkyMap App', 'Position', [140, 90, 1360, 820]);
            obj.UIFigure.CloseRequestFcn = @(src, evt) obj.on_close(src, evt);
            obj.UIFigure.WindowScrollWheelFcn = @(~, evt) obj.on_scroll(evt);

            gl = uigridlayout(obj.UIFigure, [2, 2]);
            gl.RowHeight = {'1x', 88};
            gl.ColumnWidth = {'1x', '1x'};
            gl.Padding = [8 8 8 8];
            gl.RowSpacing = 8;
            gl.ColumnSpacing = 8;

            obj.TopAxes = uiaxes(gl);
            obj.TopAxes.Layout.Row = 1;
            obj.TopAxes.Layout.Column = 1;
            hold(obj.TopAxes, 'on');
            obj.TopAxes.YDir = 'normal';
            axis(obj.TopAxes, 'equal');
            xlabel(obj.TopAxes, 'North [m]');
            ylabel(obj.TopAxes, 'East [m]');
            title(obj.TopAxes, 'Top-Down Binary Threat');
            colormap(obj.TopAxes, obj.class_cmap);
            caxis(obj.TopAxes, [1 3]);
            colorbar(obj.TopAxes, 'Ticks', [1 2 3], 'TickLabels', {'Below Terrain', 'Hidden', 'Visible'});

            obj.VertAxes = uiaxes(gl);
            obj.VertAxes.Layout.Row = 1;
            obj.VertAxes.Layout.Column = 2;
            hold(obj.VertAxes, 'on');
            obj.VertAxes.YDir = 'normal';
            xlabel(obj.VertAxes, 'Along-Slice Distance [m]');
            ylabel(obj.VertAxes, 'Altitude [m]');
            title(obj.VertAxes, 'Vertical Binary Slice');
            colormap(obj.VertAxes, obj.class_cmap);
            caxis(obj.VertAxes, [1 3]);
            colorbar(obj.VertAxes, 'Ticks', [1 2 3], 'TickLabels', {'Below Terrain', 'Hidden', 'Visible'});
            grid(obj.VertAxes, 'on');

            cp = uipanel(gl, 'Title', 'Controls');
            cp.Layout.Row = 2;
            cp.Layout.Column = [1 2];
            cgl = uigridlayout(cp, [2, 8]);
            cgl.RowHeight = {22, 34};
            cgl.ColumnWidth = {70, '1x', 100, 130, 120, 120, 95, 120};
            cgl.Padding = [6 6 6 6];

            lbl_wp = uilabel(cgl, 'Text', 'Waypoint', 'HorizontalAlignment', 'left');
            lbl_wp.Layout.Row = 1;
            lbl_wp.Layout.Column = 1;

            obj.WPSlider = uislider(cgl, ...
                'Limits', [1 obj.n_wp], ...
                'Value', 1, ...
                'MajorTicks', [], ...
                'ValueChangingFcn', @(src, evt) obj.on_slider_changing(evt), ...
                'ValueChangedFcn', @(src, evt) obj.on_slider_changed(evt));
            obj.WPSlider.Layout.Row = [1 2];
            obj.WPSlider.Layout.Column = 2;

            obj.WPLabel = uilabel(cgl, 'Text', sprintf('WP %d/%d', 1, obj.n_wp), ...
                'HorizontalAlignment', 'left');
            obj.WPLabel.Layout.Row = 1;
            obj.WPLabel.Layout.Column = 3;

            obj.PlayButton = uibutton(cgl, 'Text', 'Play', ...
                'ButtonPushedFcn', @(~, ~) obj.on_play_toggle());
            obj.PlayButton.Layout.Row = 2;
            obj.PlayButton.Layout.Column = 3;

            lbl_mode = uilabel(cgl, 'Text', 'Slice Mode');
            lbl_mode.Layout.Row = 1;
            lbl_mode.Layout.Column = 4;
            obj.ModeDropDown = uidropdown(cgl, ...
                'Items', {'Path Heading', 'Nearest Radar Bearing', 'Fixed Heading'}, ...
                'Value', obj.cfg.skymap_heading_mode, ...
                'ValueChangedFcn', @(~, ~) obj.update_view(obj.current_idx));
            obj.ModeDropDown.Layout.Row = 2;
            obj.ModeDropDown.Layout.Column = 4;

            lbl_win = uilabel(cgl, 'Text', 'Window [m]');
            lbl_win.Layout.Row = 1;
            lbl_win.Layout.Column = 5;
            obj.WindowSpinner = uispinner(cgl, ...
                'Limits', [200, 8000], 'Step', 100, ...
                'Value', obj.cfg.skymap_half_window, ...
                'ValueChangedFcn', @(~, ~) obj.update_view(obj.current_idx));
            obj.WindowSpinner.Layout.Row = 2;
            obj.WindowSpinner.Layout.Column = 5;

            lbl_slice = uilabel(cgl, 'Text', 'Slice Half [m]');
            lbl_slice.Layout.Row = 1;
            lbl_slice.Layout.Column = 6;
            obj.SliceSpinner = uispinner(cgl, ...
                'Limits', [200, 10000], 'Step', 100, ...
                'Value', obj.cfg.skymap_slice_half_length, ...
                'ValueChangedFcn', @(~, ~) obj.update_view(obj.current_idx));
            obj.SliceSpinner.Layout.Row = 2;
            obj.SliceSpinner.Layout.Column = 6;

            btn_prev = uibutton(cgl, 'Text', 'Prev', ...
                'ButtonPushedFcn', @(~, ~) obj.step_wp(-1));
            btn_prev.Layout.Row = 2;
            btn_prev.Layout.Column = 7;
            btn_next = uibutton(cgl, 'Text', 'Next', ...
                'ButtonPushedFcn', @(~, ~) obj.step_wp(1));
            btn_next.Layout.Row = 2;
            btn_next.Layout.Column = 8;

            % Initial plots
            obj.TopImg = imagesc(obj.TopAxes, obj.threat.N_vec, obj.threat.E_vec, ones(size(obj.terrain_h_grid)));
            plot(obj.TopAxes, obj.path(1, :), obj.path(2, :), 'w-', 'LineWidth', 2.2);
            plot(obj.TopAxes, obj.start_pos(1), obj.start_pos(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
            plot(obj.TopAxes, obj.goal_pos(1), obj.goal_pos(2), 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
            for r = 1:length(obj.threat.radars)
                rr = obj.threat.radars{r};
                plot(obj.TopAxes, rr.position(1), rr.position(2), 'r^', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            end
            obj.TopCurr = plot(obj.TopAxes, obj.path(1,1), obj.path(2,1), 'yo', 'MarkerSize', 9, 'MarkerFaceColor', 'y');
            obj.TopHeading = plot(obj.TopAxes, [obj.path(1,1), obj.path(1,1)], [obj.path(2,1), obj.path(2,1)], 'y-', 'LineWidth', 2.0);

            s0 = linspace(-obj.cfg.skymap_slice_half_length, obj.cfg.skymap_slice_half_length, obj.cfg.skymap_slice_n_horiz);
            obj.VertImg = imagesc(obj.VertAxes, s0, obj.slice_alt_vec, 2 * ones(length(obj.slice_alt_vec), length(s0)));
            obj.VertTerrain = plot(obj.VertAxes, s0, zeros(size(s0)), 'k-', 'LineWidth', 1.8);
            obj.VertCurr = plot(obj.VertAxes, 0, -obj.path(3,1), 'wo', 'MarkerSize', 7, 'MarkerFaceColor', 'y');
        end

        function on_slider_changing(obj, evt)
            obj.update_view(round(evt.Value));
        end

        function on_slider_changed(obj, evt)
            obj.update_view(round(evt.Value));
        end

        function on_scroll(obj, evt)
            obj.step_wp(sign(evt.VerticalScrollCount));
        end

        function step_wp(obj, step)
            if obj.n_wp < 2
                return;
            end
            idx = obj.current_idx + step;
            idx = max(1, min(obj.n_wp, idx));
            obj.update_view(idx);
        end

        function on_play_toggle(obj)
            if obj.IsPlaying
                obj.stop_timer();
                obj.PlayButton.Text = 'Play';
                obj.IsPlaying = false;
            else
                obj.start_timer();
                obj.PlayButton.Text = 'Pause';
                obj.IsPlaying = true;
            end
        end

        function start_timer(obj)
            obj.stop_timer();
            obj.PlayTimer = timer( ...
                'ExecutionMode', 'fixedRate', ...
                'Period', obj.cfg.skymap_play_period, ...
                'BusyMode', 'drop', ...
                'TimerFcn', @(~, ~) obj.on_timer_tick());
            start(obj.PlayTimer);
        end

        function stop_timer(obj)
            if ~isempty(obj.PlayTimer) && isvalid(obj.PlayTimer)
                stop(obj.PlayTimer);
                delete(obj.PlayTimer);
            end
            obj.PlayTimer = [];
        end

        function on_timer_tick(obj)
            if isempty(obj.UIFigure) || ~isvalid(obj.UIFigure)
                obj.stop_timer();
                return;
            end
            if obj.current_idx >= obj.n_wp
                obj.on_play_toggle();
                return;
            end
            obj.update_view(obj.current_idx + 1);
        end

        function on_close(obj, src, ~)
            obj.stop_timer();
            delete(src);
        end

        function update_view(obj, idx)
            if isempty(obj.UIFigure) || ~isvalid(obj.UIFigure)
                return;
            end
            idx = max(1, min(obj.n_wp, round(idx)));
            obj.current_idx = idx;

            p = obj.path(:, idx);
            N_curr = p(1);
            E_curr = p(2);
            alt_curr = -p(3);
            heading = obj.heading_at(idx);

            [class_map, ~] = obj.threat.get_binary_horizontal_slice( ...
                alt_curr, obj.visibility_threshold, obj.tm, obj.terrain_h_grid);
            obj.TopImg.CData = class_map;
            obj.TopCurr.XData = N_curr;
            obj.TopCurr.YData = E_curr;

            heading_len = max(120, 0.18 * obj.WindowSpinner.Value);
            obj.TopHeading.XData = [N_curr, N_curr + heading_len * cos(heading)];
            obj.TopHeading.YData = [E_curr, E_curr + heading_len * sin(heading)];

            N_min = max(obj.threat.N_vec(1), N_curr - obj.WindowSpinner.Value);
            N_max = min(obj.threat.N_vec(end), N_curr + obj.WindowSpinner.Value);
            E_min = max(obj.threat.E_vec(1), E_curr - obj.WindowSpinner.Value);
            E_max = min(obj.threat.E_vec(end), E_curr + obj.WindowSpinner.Value);
            xlim(obj.TopAxes, [N_min, N_max]);
            ylim(obj.TopAxes, [E_min, E_max]);
            title(obj.TopAxes, sprintf('Top-Down Threat | WP %d/%d | Alt %.1f m', idx, obj.n_wp, alt_curr));

            [class_vert, s_vec, alt_vec, terrain_line] = obj.threat.get_binary_vertical_slice( ...
                [N_curr; E_curr], heading, obj.SliceSpinner.Value, obj.cfg.skymap_slice_n_horiz, ...
                obj.slice_alt_vec, obj.visibility_threshold, obj.tm);
            obj.VertImg.XData = s_vec;
            obj.VertImg.YData = alt_vec;
            obj.VertImg.CData = class_vert;
            obj.VertTerrain.XData = s_vec;
            obj.VertTerrain.YData = terrain_line;
            obj.VertCurr.XData = 0;
            obj.VertCurr.YData = alt_curr;
            ylim(obj.VertAxes, [alt_vec(1), alt_vec(end)]);
            title(obj.VertAxes, sprintf('Vertical Slice | %s | %.1f deg', obj.ModeDropDown.Value, rad2deg(heading)));

            obj.WPSlider.Value = idx;
            obj.WPLabel.Text = sprintf('WP %d/%d', idx, obj.n_wp);
        end

        function heading = heading_at(obj, idx)
            mode = string(obj.ModeDropDown.Value);
            switch mode
                case "Nearest Radar Bearing"
                    p = obj.path(:, idx);
                    d_best = inf;
                    heading = 0;
                    for r = 1:length(obj.threat.radars)
                        rp = obj.threat.radars{r}.position;
                        d = hypot(rp(1) - p(1), rp(2) - p(2));
                        if d < d_best
                            d_best = d;
                            heading = atan2(rp(2) - p(2), rp(1) - p(1));
                        end
                    end
                case "Fixed Heading"
                    heading = deg2rad(obj.cfg.skymap_fixed_heading_deg);
                otherwise
                    if idx < obj.n_wp
                        d = obj.path(1:2, idx + 1) - obj.path(1:2, idx);
                    else
                        d = obj.path(1:2, idx) - obj.path(1:2, idx - 1);
                    end
                    if norm(d) < 1e-9
                        heading = 0;
                    else
                        heading = atan2(d(2), d(1));
                    end
            end
        end
    end
end

function val = get_cfg(cfg, name, default)
if isfield(cfg, name)
    val = cfg.(name);
else
    val = default;
end
end
