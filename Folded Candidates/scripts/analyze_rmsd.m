
prefix = 'NuG2';

structures_root = '../structures/';
trajectories_root = 'd:\Files\Workspace\Sosnick\dataset 2 - 20190308\';

pdb = pdbread(sprintf('%s%s.pdb', structures_root, prefix));
atoms = pdb.Model.Atom;

load(sprintf('%s%s_mapped_atoms.mat', structures_root, prefix));

pdb_atoms_xyz = [[atoms.X]' [atoms.Y]' [atoms.Z]'];

atoms_index_in_pdb = [mapped_atoms.index_in_pdb]';
which_atoms_are_in_pdb = [0 ~= atoms_index_in_pdb];
pdb_ref_xyz = pdb_atoms_xyz(atoms_index_in_pdb(which_atoms_are_in_pdb), :);

which_atoms_are_calpha = strcmp('CA', {mapped_atoms.name});
which_atoms_are_calpha = which_atoms_are_calpha(which_atoms_are_in_pdb);

pdb_ref_xyz_ca = pdb_ref_xyz(which_atoms_are_calpha, :);

%%

windowLength = 500 / 2; % These frames were saved every 2nd one (see email from Tobin)
windowStride = 50 / 2;
%windowLength = 1000;
%windowStride = 50;

for sim_num = 0:27
    sim_num
    filename = sprintf('%s%s.vtf/%s.%d.vtf.mat', trajectories_root, prefix, prefix, sim_num)
    load(filename);
    
    frames = frames(:, which_atoms_are_in_pdb, :);
    frames_ca = frames(:, which_atoms_are_calpha, :);
    
    rmsd = arrayfun(@(i)...
        CalculateRmsd(squeeze(frames_ca(i, :, :)), pdb_ref_xyz_ca), ...
        1:size(frames_ca, 1));
    
    save(sprintf('../results201903/%s/%s_%d - rmsd.mat', prefix, prefix, sim_num), 'rmsd');

    
    if (0)
        %%
        % TODO: Sliding rmsd vs. rmsd of sliding-averaged configurations
        [rmsdMovingAvg, slidingFrameMeanIndex] = SlidingAverage(rmsd, windowLength, windowStride);
        
        figure(1);
        histogram(rmsdMovingAvg)
        ylabel('Å');
    end

    if (0)
        %%
        framesIndexes = 1:size(frames, 1);
        movingAvgKernel = ones(1, windowLength) * (1/windowLength);
        
        slidingFrameMeanIndex = conv(framesIndexes, movingAvgKernel, 'valid');
        slidingFrameMeanIndex = slidingFrameMeanIndex(1:windowStride:end);
        
        framesCaMovingAvg = convn(frames_ca, reshape(movingAvgKernel, [numel(movingAvgKernel), 1, 1]), 'valid');
        framesCaMovingAvg = framesCaMovingAvg(1:windowStride:end, :, :);
        
        rmsdMovingAvg2 = arrayfun(@(i)...
            CalculateRmsd(squeeze(framesCaMovingAvg(i, :, :)), pdb_ref_xyz_ca), ...
            1:size(framesCaMovingAvg, 1));
        
        figure(1);
        histogram(rmsdMovingAvg2)
        ylabel('Å');
    end


end

%%
figure(1);
histogram(rmsd, 200)
ylabel('Å');




