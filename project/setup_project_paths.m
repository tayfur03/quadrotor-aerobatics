function repo_root = setup_project_paths()
%SETUP_PROJECT_PATHS Add project folders to the MATLAB path and normalize CWD.

this_file = mfilename('fullpath');
project_dir = fileparts(this_file);
repo_root = fileparts(project_dir);

cd(repo_root);

path_dirs = {
    repo_root
    fullfile(repo_root, 'project')
    fullfile(repo_root, 'demos')
    fullfile(repo_root, 'tests')
    fullfile(repo_root, 'flight')
    fullfile(repo_root, 'math')
    fullfile(repo_root, 'motion_planner')
    fullfile(repo_root, 'radar')
    fullfile(repo_root, 'safety')
    fullfile(repo_root, 'terrain')
};

for i = 1:numel(path_dirs)
    dir_i = path_dirs{i};
    if exist(dir_i, 'dir')
        addpath(dir_i);
    end
end
end
