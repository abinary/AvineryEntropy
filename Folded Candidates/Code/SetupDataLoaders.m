function [loaders] = SetupDataLoaders(sim_name)

loaders = struct();

trajectories_root = 'd:/Files/Workspace/Sosnick/';
data_root = 'c:/Files/The Lab/Projects/Bridging Physical and Information Entropy/Simulation/Sosnick/Data/20190308/';
results_root = 'c:/Files/The Lab/Projects/Bridging Physical and Information Entropy/Simulation/Sosnick/results201903/';
structures_root = 'c:/Files/The Lab/Projects/Bridging Physical and Information Entropy/Sosnick/Structures';

model_filepath = sprintf('%s%s_mapped_atoms.mat', structures_root, sim_name);
loaders.LoadModelFile = @()load(model_filepath);

pdb_filepath = sprintf('%s%s.pdb', structures_root, sim_name);
loaders.LoadPdbFile = @()pdbread(pdb_filepath);

trajectory_file_prefix = sprintf('%s%s_trajectories/', trajectories_root, sim_name);
loaders.TrajectoryFilepath = @(n)sprintf('%s%s_%d.vtf.mat', trajectory_file_prefix, sim_name, n);
loaders.LoadTrajectoryFile = @(n)load(loaders.TrajectoryFilepath(n));
loaders.LoadTrajectoryFileCAlpha = @LoadTrajectoryFileCAlpha;

loaders.SCPotentialFilepath = @(n)sprintf('%s%s/%s.%d.sc_potential.pkl.mat', data_root, sim_name, sim_name, n);
loaders.LoadSCPotentialFile = @(n)load(loaders.SCPotentialFilepath(n));

loaders.SCEntropyFilepath = @(n)sprintf('%s%s/%s.%d.sc_entropy.pkl.mat', data_root, sim_name, sim_name, n);
loaders.LoadSCEntropyFile = @(n)load(loaders.SCEntropyFilepath(n));

loaders.RamaFilepath = @(n)sprintf('%s%s.%d.rama.pkl.mat', data_root, sim_name, n);
loaders.LoadRamaFile = @(n)load(loaders.RamaFilepath(n));

loaders.rmsd_folder = sprintf('%s%s_trajectories/', root, sim_name);
loaders.RmsdFile = @(n)sprintf('%s%s_%d - rmsd.mat', results_root, loaders.rmsd_folder, sim_name, n);

loaders.LoadRmsd = @LoadRmsd;
loaders.LoadEnergy = @LoadEnergy;

    % Load a trajectory (+model) and remove all non-CA atoms
    function [trajectory] = LoadTrajectoryFileCAlpha(n)
        model = loaders.LoadModelFile();
        which_atoms_are_calpha = strcmp('CA', {model.mapped_atoms.name});
        
        trajectory = loaders.LoadTrajectoryFile(n);
        trajectory.frames = trajectory.frames(:, which_atoms_are_calpha, :); % Keep only C-alpha atoms
    end

    function [pe] = LoadEnergy(n, num_of_frames_to_skip)
        pe = load(sprintf('%s.%d.pe.pkl.mat', sim_name, n));
        pe = pe.data';
        
        if (nargin >= 2)
            pe(1:num_of_frames_to_skip) = []; % throw out the first 2000 frames
        end
    end

    function [rmsd] = LoadRmsd(n, num_of_frames_to_skip)
        rmsd = load(loaders.RmsdFile(n));
        rmsd = rmsd.rmsd;
        
        if (nargin >= 2)
            % throw out the first X frames. RMSD was calculated for frames saved every other one (hence the /2)
            rmsd(1:(num_of_frames_to_skip/2)) = [];
        end
    end
end
