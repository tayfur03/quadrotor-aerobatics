classdef skymap_viewer < handle
%SKYMAP_VIEWER Sliding SkyMap-style GUI for binary radar threat inspection.
%   Shows:
%     - Top-down binary threat map (sliding window)
%     - Vertical binary slice along current heading
%   Color classes:
%     1 -> below terrain
%     2 -> hidden/safe
%     3 -> visible

    properties (Access = private)
        cfg
        tm
        threat
        path
        start_pos
        goal_pos
        visibility_threshold
        n_wp

        fig
        ax_map
        ax_vert
        slider
        txt

        h_map_img
        h_map_curr
        h_map_heading
        h_vert_img
        h_vert_terrain
        h_vert_curr

        terrain_h_grid
        class_cmap
        slice_alt_vec
    end

    methods
        function obj = skymap_viewer(cfg, tm, threat, path_rrt, start_pos, goal_pos, visibility_threshold)
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
            obj.class_cmap = [0.35 0.35 0.35; 0.10 0.60 0.10; 0.80 0.10 0.10];

            if obj.n_wp < 2
                warning('skymap_viewer:shortPath', 'Path has fewer than 2 points. GUI not created.');
                return;
            end

            [N_grid, E_grid] = meshgrid(threat.N_vec, threat.E_vec);
            terrain_h = tm.get_height(N_grid(:), E_grid(:));
            obj.terrain_h_grid = reshape(terrain_h, size(N_grid));
            obj.slice_alt_vec = linspace(threat.alt_vec(1), threat.alt_vec(end), obj.cfg.skymap_slice_n_alt);

            obj.build_ui();
            obj.update(1);
        end

        function update(obj, idx)
            if isempty(obj.fig) || ~isgraphics(obj.fig, 'figure')
                return;
            end
            idx = max(1, min(obj.n_wp, round(idx)));

            p = obj.path(:, idx);
            N_curr = p(1);
            E_curr = p(2);
            alt_curr = -p(3);
            heading = obj.heading_at(idx);

            [class_map, ~] = obj.threat.get_binary_horizontal_slice( ...
                alt_curr, obj.visibility_threshold, obj.tm, obj.terrain_h_grid);
            set(obj.h_map_img, 'CData', class_map);
            set(obj.h_map_curr, 'XData', N_curr, 'YData', E_curr);

            heading_len = max(120, 0.18 * obj.cfg.skymap_half_window);
            set(obj.h_map_heading, ...
                'XData', [N_curr, N_curr + heading_len * cos(heading)], ...
                'YData', [E_curr, E_curr + heading_len * sin(heading)]);

            N_min = max(obj.threat.N_vec(1), N_curr - obj.cfg.skymap_half_window);
            N_max = min(obj.threat.N_vec(end), N_curr + obj.cfg.skymap_half_window);
            E_min = max(obj.threat.E_vec(1), E_curr - obj.cfg.skymap_half_window);
            E_max = min(obj.threat.E_vec(end), E_curr + obj.cfg.skymap_half_window);
            xlim(obj.ax_map, [N_min, N_max]);
            ylim(obj.ax_map, [E_min, E_max]);
            title(obj.ax_map, sprintf('Top-Down Threat | WP %d/%d | Alt %.1f m', idx, obj.n_wp, alt_curr));

            [class_vert, s_vec, alt_vec, terrain_line] = obj.threat.get_binary_vertical_slice( ...
                [N_curr; E_curr], heading, obj.cfg.skymap_slice_half_length, ...
                obj.cfg.skymap_slice_n_horiz, obj.slice_alt_vec, obj.visibility_threshold, obj.tm);
            set(obj.h_vert_img, 'XData', s_vec, 'YData', alt_vec, 'CData', class_vert);
            set(obj.h_vert_terrain, 'XData', s_vec, 'YData', terrain_line);
            set(obj.h_vert_curr, 'XData', 0, 'YData', alt_curr);
            ylim(obj.ax_vert, [alt_vec(1), alt_vec(end)]);
            title(obj.ax_vert, sprintf('Vertical Slice | Heading %.1f deg', rad2deg(heading)));

            if isgraphics(obj.slider)
                set(obj.slider, 'Value', idx);
            end
            if isgraphics(obj.txt)
                set(obj.txt, 'String', sprintf('WP %d/%d', idx, obj.n_wp));
            end
        end
    end

    methods (Access = private)
        function cfg = apply_defaults(~, cfg)
            cfg.skymap_half_window = get_cfg(cfg, 'skymap_half_window', 1800);
            cfg.skymap_slice_half_length = get_cfg(cfg, 'skymap_slice_half_length', 2500);
            cfg.skymap_slice_n_horiz = get_cfg(cfg, 'skymap_slice_n_horiz', 181);
            cfg.skymap_slice_n_alt = get_cfg(cfg, 'skymap_slice_n_alt', 120);
        end

        function build_ui(obj)
            obj.fig = figure('Name', 'SkyMap Sliding Threat View', ...
                'Position', [140, 90, 1320, 760], 'Color', 'w');
            tl = tiledlayout(obj.fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
            obj.ax_map = nexttile(tl, 1);
            obj.ax_vert = nexttile(tl, 2);

            obj.h_map_img = imagesc(obj.ax_map, obj.threat.N_vec, obj.threat.E_vec, ones(size(obj.terrain_h_grid)));
            set(obj.ax_map, 'YDir', 'normal');
            axis(obj.ax_map, 'equal');
            hold(obj.ax_map, 'on');
            colormap(obj.ax_map, obj.class_cmap);
            caxis(obj.ax_map, [1 3]);
            colorbar(obj.ax_map, 'Ticks', [1 2 3], 'TickLabels', {'Below Terrain', 'Hidden', 'Visible'});
            plot(obj.ax_map, obj.path(1, :), obj.path(2, :), 'w-', 'LineWidth', 2.3);
            plot(obj.ax_map, obj.start_pos(1), obj.start_pos(2), 'go', 'MarkerSize', 9, 'MarkerFaceColor', 'g');
            plot(obj.ax_map, obj.goal_pos(1), obj.goal_pos(2), 'bs', 'MarkerSize', 9, 'MarkerFaceColor', 'b');
            for r = 1:length(obj.threat.radars)
                rr = obj.threat.radars{r};
                plot(obj.ax_map, rr.position(1), rr.position(2), 'r^', 'MarkerSize', 9, 'MarkerFaceColor', 'r');
            end
            obj.h_map_curr = plot(obj.ax_map, obj.path(1, 1), obj.path(2, 1), 'yo', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
            obj.h_map_heading = plot(obj.ax_map, [obj.path(1, 1), obj.path(1, 1)], [obj.path(2, 1), obj.path(2, 1)], ...
                'y-', 'LineWidth', 2.2);
            xlabel(obj.ax_map, 'North [m]');
            ylabel(obj.ax_map, 'East [m]');
            title(obj.ax_map, 'Top-Down Binary Threat (Sliding Window)');

            s_vec = linspace(-obj.cfg.skymap_slice_half_length, obj.cfg.skymap_slice_half_length, obj.cfg.skymap_slice_n_horiz);
            obj.h_vert_img = imagesc(obj.ax_vert, s_vec, obj.slice_alt_vec, 2 * ones(length(obj.slice_alt_vec), length(s_vec)));
            set(obj.ax_vert, 'YDir', 'normal');
            hold(obj.ax_vert, 'on');
            colormap(obj.ax_vert, obj.class_cmap);
            caxis(obj.ax_vert, [1 3]);
            colorbar(obj.ax_vert, 'Ticks', [1 2 3], 'TickLabels', {'Below Terrain', 'Hidden', 'Visible'});
            obj.h_vert_terrain = plot(obj.ax_vert, s_vec, zeros(size(s_vec)), 'k-', 'LineWidth', 1.8);
            obj.h_vert_curr = plot(obj.ax_vert, 0, -obj.path(3, 1), 'wo', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
            xlabel(obj.ax_vert, 'Along-Heading Distance [m]');
            ylabel(obj.ax_vert, 'Altitude [m]');
            title(obj.ax_vert, 'Vertical Slice Along Current Heading');
            grid(obj.ax_vert, 'on');

            slider_step = [1 / (obj.n_wp - 1), min(20 / (obj.n_wp - 1), 1)];
            obj.slider = uicontrol(obj.fig, 'Style', 'slider', 'Units', 'normalized', ...
                'Position', [0.08 0.02 0.75 0.045], ...
                'Min', 1, 'Max', obj.n_wp, 'Value', 1, ...
                'SliderStep', slider_step, ...
                'Callback', @(src, ~) obj.on_slider(src));
            obj.txt = uicontrol(obj.fig, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.845 0.02 0.145 0.045], ...
                'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
                'String', 'WP 1');

            set(obj.fig, 'WindowScrollWheelFcn', @(~, evt) obj.on_scroll(evt));
        end

        function on_slider(obj, src)
            idx = round(get(src, 'Value'));
            obj.update(idx);
        end

        function on_scroll(obj, evt)
            if isempty(obj.slider) || ~isgraphics(obj.slider)
                return;
            end
            curr = round(get(obj.slider, 'Value'));
            idx = curr + sign(evt.VerticalScrollCount);
            idx = max(1, min(obj.n_wp, idx));
            obj.update(idx);
        end

        function heading = heading_at(obj, idx)
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

function val = get_cfg(cfg, name, default)
if isfield(cfg, name)
    val = cfg.(name);
else
    val = default;
end
end
