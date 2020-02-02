

kmedoids = true;

%results_file = 'folded-analysis 20191130-201509.results BBL-136.mat';
results_file = 'folded-analysis 20191130-201205.results BBL-139.mat';
%results_file = 'folded-analysis 20191130-202059.results A3D-222.mat';
%results_file = 'folded-analysis 20191130-202422.results NuG2-136.mat';
%results_file = 'folded-analysis 20191208-174919.results NuG2-136 with additional window data.mat';
%results_file = 'folded-analysis 20191205-184607.results NuG2-136 proper linkage.mat';
%results_file = 'folded-analysis 20191216-133055.results.mat'; % NuG2, with SC contributions
%results_file = 'folded-analysis 20191217-102206.results.mat'; % BBL, with SC contributions


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

candidate_frames = [];
candidate_rmsd = [];

candidate_num = 0;
for sim = 1:numel(representatives_indices)
    if (isempty(representatives_indices{sim}))
        continue;
    end
    
    fprintf('Loading sim %d/%d \n', sim, numel(results.TrajectoryFiles));
    
    [folder, name, ext] = fileparts(results.TrajectoryFiles{sim});
    parsed = regexp(name, '.*?\.(?<idx>\d+)\.', 'names');
    idx = str2double(parsed.idx);
    vtf = load(vtf_filepath(idx));
    ca = load(ca_filepath(idx));
    
    frame_indexes = representatives_indices{sim};
    for i = 1:numel(frame_indexes)
        candidate_num = candidate_num + 1;
        xyz_ca = squeeze(ca.frames(frame_indexes(i), :, :));
        candidate_frames(end + 1, :, :) = xyz_ca;
        candidate_rmsd(end + 1) = representatives_rmsd{sim}(i);

        if (1)
            xyz = squeeze(vtf.frames(frame_indexes(i), :, :));
            
            [~, ~, transform] = procrustes(ref_ca_xyz, xyz_ca, 'reflection', false, 'scaling', false);
            
            xyz = xyz * transform.T + transform.c(1, :);
            
            for j = find(w) % Set mapped atoms (as indicated by "w")
                pdb.Model.Atom(mapped_atoms(j).index_in_pdb).X = xyz(j, 1);
                pdb.Model.Atom(mapped_atoms(j).index_in_pdb).Y = xyz(j, 2);
                pdb.Model.Atom(mapped_atoms(j).index_in_pdb).Z = xyz(j, 3);
            end
            
            if (~isempty(representatives_rmsd))
                pdbwrite(sprintf('%s\\[%02d] %s.%d.%d (RMSD %0.2f).pdb', output_folder, candidate_num, sim_name, idx, frame_indexes(i), representatives_rmsd{sim}(i)), pdb);
            else
                pdbwrite(sprintf('%s\\[%02d] %s.%d.%d.pdb', output_folder, candidate_num, sim_name, idx, frame_indexes(i)), pdb);
            end
        end
    end
end

if (1)
    %candidate_frames = reshape(ref_ca_xyz, [1 size(ref_ca_xyz)]);
    %candidate_frames = cat(1, reshape(ref_ca_xyz, [1 size(ref_ca_xyz)]), candidate_frames);
    candidate_frames(end + 1, :, :) = ref_ca_xyz;
end

%% Calculate RMSD between candidates and all frames
candidates_rmsd = [];
for sim = 1:numel(results.TrajectoryFiles)
    fprintf('Loading sim %d/%d \n', sim, numel(results.TrajectoryFiles));
    ca = load(ca_filepath(sim - 1));
    
    fprintf('Calculating RMSD ... \n');
    
    for i = 1:ceil(size(ca.frames, 1) / 10000)
        start_idx = (i - 1) * 10000;
        len = min(size(ca.frames, 1) - start_idx, 10000);
        rmsd = AllPairsRmsd2(candidate_frames, ca.frames(start_idx + [1:len], :, :), 1, 2, 3);
        candidates_rmsd = [candidates_rmsd, rmsd];
    end
end

save('candidates_rmsd.mat', 'candidates_rmsd');

%%
S_list = [];
std_list = [];
p_list = [];

figure(1);
clf;
cla;

edges = 0:0.25:35;
centers = 0.5 * (edges(1:end-1) + edges(2:end));

hold on;
for i = 1:size(candidates_rmsd, 1)
    n = histcounts(candidates_rmsd(i, :), edges);
    p = n ./ sum(n);
    std_list(i) = std(candidates_rmsd(i, :));
    S_list(i) = -sum(p .* log2(p + 1e-14));
    plot(centers, p * 1e2 + 3 * (i - 1), 'LineWidth', 2);
    plot(centers, centers .* 0 + 3 * (i - 1), '--k', 'LineWidth', 1);
    
    p_list(end + 1, :) = prctile(candidates_rmsd(i, :), [0.01 0.02 0.05 0.1 0.3 0.5 0.7 0.9 0.95 0.99]);
end
hold off;

SetMyDefaultFigureSettings();

ax = gca();
ax.YTick = [0];

xlabel('RMSD (Å)');
ylabel('Percentage (%)');

xlim([0 20]);
%%
figure(2);
plot3(candidate_rmsd, S_list(2:end), std_list(2:end), '*b');
hold on;
plot3(0, S_list(1), std_list(1), 'or');
hold off;
%plot(candidate_rmsd, S_list, '*');
%plot(candidate_rmsd, std_list, '*');
%plot(std_list, S_list, '*');
%plot(std_list, S_list, '*');

xlabel('RMSD');
ylabel('Information Entropy');
zlabel('Std');

%%
[~, order] = sort(candidate_rmsd);
figure(3);
plot(candidate_rmsd(order), p_list(order + 1, :), '-x');
hold on;
plot(p_list(1:end, :) .* 0, p_list(1, :), 'o');
hold off;

figure(4);
c = corrcoef(candidates_rmsd');
%c(c < 0.5) = 0;
imagesc(c); 
colorbar();

%adjacencyMatrix = 1 - corrcoef(candidates_rmsd');
%imagesc(adjacencyMatrix);
axis xy

SetMyDefaultFigureSettings();

%%
adjacencyMatrix = 1 - corrcoef(candidates_rmsd');

%Y = cmdscale(adjacencyMatrix);
%Y = Y(:, 1:2);

figure(5);
%plot(Y(:, 1), Y(:, 2), 'x');

G = graph(adjacencyMatrix(1:10, 1:10));
G = rmedge(G, find(G.Edges.Weight > 0.5));

%figure(5);
%p = plot(G, 'EdgeColor', double(G.Edges.Weight > 0.5) .* [1 1 1]);
p = plot(G, 'EdgeColor', 'k', 'LineWidth', 2, 'Layout', 'force');
p.NodeFontSize = 16;
p.NodeFontWeight = 'bold';
p.EdgeAlpha = 0.5;

data_table = [candidate_rmsd(:) results.CandidateWindows(results.KMedoids.MedoidIndexes, :)];

for i = 1:numel(p.NodeLabel)
    p.NodeLabel{i} = sprintf('%d (%0.1f)', i, candidate_rmsd(i));
end

w = G.Edges.Weight;

Scale = @(x, a, b)(x - a) ./ (b - a);

% p.EdgeColor = (w < 0.1) .* Scale(w, 0.1, 0.0) .* [1 0 0] + ...
%     (w < 0.2) .* Scale(w, 0.2, 0.1) .* [0 0 1];

p.EdgeColor = (w < 0.1) .* [1 0 0] + ...
    (w < 0.2 & w > 0.1) .* [0 0 1];
% 
1

f = gcf();
f.KeyPressFcn = @GraphKeyPressFcn;
f.UserData = {p, data_table};
w = 2 - G.Edges.Weight;
%G.Edges.LWidths = (w - min(w)) / (max(w) - min(w)) * 5 + 1e-3;
%G.Edges.LWidths = double(G.Edges.Weight > 0.5);
%p.LineWidth = G.Edges.LWidths;

%[T, pred] = minspantree(G);
%p = plot(T, 'Layout', 'layered');
%highlight(p, T, 'LineWidth', 3, 'EdgeColor', 'g');

SetMyDefaultFigureSettings();
