

kmedoids = true;

%results_file = 'folded-analysis 20191130-201509.results BBL-136.mat';
%results_file = 'folded-analysis 20191130-201205.results BBL-139.mat';
%results_file = 'folded-analysis 20191130-202059.results A3D-222.mat';
%results_file = 'folded-analysis 20191130-202422.results NuG2-136.mat';
results_file = 'folded-analysis 20191205-184607.results NuG2-136 proper linkage.mat';
%results_file = 'folded-analysis 20191207-232500.results NuG2-136 linkage complete.mat';

%load('folded-analysis 20191112-180440.results.mat');
%load('folded-analysis NuG2 20191112-184714.results.mat');
%load('folded-analysis 20191128-104901.results.mat');
%load('folded-analysis BBL 20191112-182531.results.mat');

[~, filename, ~] = fileparts(results_file);

if (kmedoids)
    output_folder = [filename ' - kmedoids pdbs'];
else
    output_folder = [filename ' - linkage pdbs'];
end

if (~exist(output_folder, 'dir'))
    mkdir(output_folder);
end

load(results_file);

%%
figure(2);
dendrogram(results.Linkage.Hierarchy);

%%
[folder, name, ext] = fileparts(results.TrajectoryFiles{1});
where_dot = strfind(name, '.');

sim_name = name(1:where_dot(1)-1);

%%
%sim_name = 'nug2';

pdb_filepath = ['c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Sosnick\Structures\' ...
    sprintf('%s.pdb', sim_name)];
pdb = pdbread(pdb_filepath);

ref_ca_xyz = load(['c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Sosnick\Structures\' ...
    sprintf('%s.ref_ca_xyz.mat', sim_name)]);
ref_ca_xyz = ref_ca_xyz.ref_ca_xyz;

load(['c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Sosnick\Structures\' ...
    sprintf('%s_mapped_atoms.mat', sim_name)]);

vtf_filepath = @(i)['d:\Files\Workspace\Sosnick\dataset 2 - 20190308\' ...
    sprintf('%s.vtf\\%s.%d.vtf.mat', sim_name, sim_name, i)];

ca_filepath = @(i)['d:\Files\Workspace\Sosnick\dataset 2 - 20190308\' ...
    sprintf('%s.vtf\\%s.%d.CA.mat', sim_name, sim_name, i)];


%%
if (kmedoids)
    representatives = results.KMedoids.RepresentativeSimAndFrameIndex;
else
    representatives = results.Linkage.RepresentativeSimAndFrameIndex;
end

representatives_indices = accumarray(representatives(:, 1), representatives(:, 2), [], @(x){x});

if (size(representatives, 2) > 2)
    representatives_rmsd = accumarray(representatives(:, 1), representatives(:, 3), [], @(x){x});
else
    representatives_rmsd = [];
end

w = [mapped_atoms.in_pdb];

for sim = 1:numel(representatives_indices)
    if (isempty(representatives_indices{sim}))
        continue;
    end
    
    [folder, name, ext] = fileparts(results.TrajectoryFiles{sim});
    parsed = regexp(name, '.*?\.(?<idx>\d+)\.', 'names');
    idx = str2double(parsed.idx);
    vtf = load(vtf_filepath(idx));
    ca = load(ca_filepath(idx));
    
    frame_indexes = representatives_indices{sim};
    for i = 1:numel(frame_indexes)
        xyz_ca = squeeze(ca.frames(frame_indexes(i), :, :));
        xyz = squeeze(vtf.frames(frame_indexes(i), :, :));
        
        [~, ~, transform] = procrustes(ref_ca_xyz, xyz_ca, 'reflection', false, 'scaling', false);
        
        xyz = xyz * transform.T + transform.c(1, :);

        for j = find(w) % Set mapped atoms (as indicated by "w")
            pdb.Model.Atom(mapped_atoms(j).index_in_pdb).X = xyz(j, 1);
            pdb.Model.Atom(mapped_atoms(j).index_in_pdb).Y = xyz(j, 2);
            pdb.Model.Atom(mapped_atoms(j).index_in_pdb).Z = xyz(j, 3);
        end
        
        if (~isempty(representatives_rmsd))
            pdbwrite(sprintf('%s\\%s.%d.%d (RMSD %0.2f).pdb', output_folder, sim_name, idx, frame_indexes(i), representatives_rmsd{sim}(i)), pdb);
        else
            pdbwrite(sprintf('%s\\%s.%d.%d.pdb', output_folder, sim_name, idx, frame_indexes(i)), pdb);
        end
        
    end
end
