function [trajectory] = LoadTrajectoryCAlpha(sim_name, n)

trajectory_file_prefix = sprintf('d:/Files/Workspace/Sosnick/%s_trajectories/', sim_name);
TrajectoryFilepath = @(n)sprintf('%s%s_%d.vtf.mat', trajectory_file_prefix, sim_name, n);
LoadTrajectoryFile = @(n)load(TrajectoryFilepath(n));

model_filepath = sprintf('d:/Files/Workspace/Sosnick/%s_mapped_atoms.mat', sim_name);
LoadModelFile = @()load(model_filepath);

model = LoadModelFile();
which_atoms_are_calpha = strcmp('CA', {model.mapped_atoms.name});

trajectory = LoadTrajectoryFile(n);
trajectory.frames = trajectory.frames(:, which_atoms_are_calpha, :); % Keep only C-alpha atoms

end
