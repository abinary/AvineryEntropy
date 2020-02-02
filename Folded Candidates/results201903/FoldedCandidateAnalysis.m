function [results] = FoldedCandidateAnalysis(trajectory_files, candidate_windows, window_length, options)

if (nargin == 0)
    options = struct();
    %options.ChooseEvenlyAcrossSims = true;
    %options.NumOfCandidates = 1000;
    options.NumOfFinalists = 10;
    options.MinimalClusterSize = 10;
    options.WindowRepresentativeSelectionMethod = 2; % 1/2/3 - first/middle/last frame, 4 - largest std for all-pairs RMSD
    
    results = options;
    return;
    
    %load('FoldedCandidateAnalysis.debug.mat');
end

results = struct();
results.TrajectoryFiles = trajectory_files;
results.CandidateWindows = candidate_windows;
results.WindowLength = window_length;
results.Options = options;

if (0)
    %%
    save('FoldedCandidateAnalysis.debug.mat');
    return;
end

%%
addpath('../Code');
addpath('../../Simulation/Common Code');

%%
[~, order] = sort(candidate_windows(:, 1)); % Sort by simulation
candidate_windows = candidate_windows(order, :);

window_start_by_sim = accumarray(candidate_windows(:, 1), candidate_windows(:, 2), [numel(trajectory_files) 1], @(x) {x});

%% Load windows

representative_frames_from_all_windows = [];
representative_frame_indexes_from_all_windows = [];

% For debug
if (1)
    [folder, name, ext] = fileparts(results.TrajectoryFiles{1});
    where_dot = strfind(name, '.');
    
    sim_name = name(1:where_dot(1)-1);
    
    %sim_name = 'nug2';
    rmsd = [];
    all_rmsd = {};
    ref_ca_xyz = [];
    ref_ca_xyz = load(['c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Sosnick\Structures\' sprintf('%s.ref_ca_xyz.mat', sim_name)]);
    ref_ca_xyz = ref_ca_xyz.ref_ca_xyz;
end

fprintf('Loading frames ... \n');
for sim = 1:numel(window_start_by_sim)
    if (isempty(window_start_by_sim{sim}))
        continue;
    end
    
    fprintf('Loading from sim %d/%d \n', sim, numel(window_start_by_sim));
    
    frames = load(trajectory_files{sim});
    frames = frames.frames;
    
    % for debug
    if (1)
        rmsd = arrayfun(@(i)procrustes_rmsd(squeeze(frames(i, :, :)), ref_ca_xyz), 1:size(frames, 1));
        all_rmsd{end+1} = rmsd;
    end
    
    window_starts = window_start_by_sim{sim};
    
    frame_indexes = arrayfun(@(s)s - 1 + [1:window_length]', window_starts, 'UniformOutput', false);
    %frame_indexes = vertcat(frame_indexes{:});
    %frame_indexes = unique(frame_indexes);
    
    frames_per_window = cellfun(@(i)frames(i, :, :), frame_indexes, 'UniformOutput', false);
    
    % For debug
    if (1)
        rmsd_per_window = cellfun(@(i)rmsd(i), frame_indexes, 'UniformOutput', false);
    end
    
    for window_idx = 1:numel(frames_per_window)
        [rep_frame, rep_frame_index] = ChooseRepresentativeFromWindow(frames_per_window{window_idx});
        rep_frame_sim_index = frame_indexes{window_idx}(rep_frame_index);
        
        representative_frames_from_all_windows = ...
            [representative_frames_from_all_windows; rep_frame];
        
        representative_frame_indexes_from_all_windows = ...
            [representative_frame_indexes_from_all_windows; [sim rep_frame_sim_index rmsd_per_window{window_idx}(rep_frame_index)]];
    end
end

if (0)
    figure(3);
    histogram(representative_frame_indexes_from_all_windows(:, 3));
    save('FoldedCandidateAnalysis.checkpoint1.debug.mat');
end

%load('FoldedCandidateAnalysis.checkpoint1.debug.mat');
%%
frames = representative_frames_from_all_windows;

fprintf('Calculating all-pairs RMSD (this takes a short while) ... \n');

t = tic();
%ij = AllPairs(size(frames, 1), 1);
%all_pairs_rmsd = arrayfun(@(k)procrustes_rmsd(squeeze(frames(ij(k, 1), :, :)), squeeze(frames(ij(k, 2), :, :))), 1:size(ij, 1));
[all_pairs_rmsd, ij] = AllPairsRmsd(frames);
toc(t)
%figure(10); plot(all_pairs_rmsd, rmsd2, '.')

%%
rmsd_matrix = zeros(size(frames, 1));
%rmsd_matrix(sub2ind(size(rmsd_matrix), ij(:, 1), ij(:, 2))) = all_pairs_rmsd;
%rmsd_matrix(sub2ind(size(rmsd_matrix), ij(:, 2), ij(:, 1))) = all_pairs_rmsd;

rmsd_matrix = squareform(all_pairs_rmsd(:));

figure(1); imagesc(rmsd_matrix); colorbar
axis xy;

%%
%save('FoldedCandidateAnalysis.checkpoint2.debug.mat');

results.RMSD_matrix = rmsd_matrix;
results.RMSD_with_ref = all_rmsd;

%%

[idx,~,sumd,D,midx,info] = kmedoids([1:size(frames, 1)]', ...
    options.NumOfFinalists, 'Distance', @(i, j)rmsd_matrix(i,j));

results.KMedoids = struct();
results.KMedoids.Assignment = idx;
results.KMedoids.ClusterSize = accumarray(idx, 1, [], @(x)numel(x));

results.KMedoids.RMSDtoClusters = D;
results.KMedoids.MedoidIndexes = midx;
results.KMedoids.Info = info;
results.KMedoids.ClusterCorrelation = corrcoef(D);
results.KMedoids.MeanInClusterRMSD = sumd ./ (results.KMedoids.ClusterSize - 1);
results.KMedoids.MedoidSimAndFrameIndex = representative_frame_indexes_from_all_windows(midx, 1:2);

clusters = accumarray(idx, [1:numel(idx)], [], @(x){x});
results.KMedoids.Clusters = cellfun(@(c)representative_frame_indexes_from_all_windows(c, 1:2),...
    clusters, 'UniformOutput', false);

representatives = cellfun(@(c)IndexOfMax(std(rmsd_matrix(c, c))), clusters);

results.KMedoids.RepresentativeSimAndFrameIndex = representative_frame_indexes_from_all_windows(midx, :);

%figure(2); imagesc(D)
%figure(3); imagesc(corrcoef(D))

if (1)
    %%
    results.Linkage = struct();
    
    % Has to be a row vector for "linkage"
    all_pairs_rmsd = reshape(all_pairs_rmsd, [1 numel(all_pairs_rmsd)]);
    
    %Z = linkage(all_pairs_rmsd, 'single');
    %Z = linkage(all_pairs_rmsd, 'ward');
    Z = linkage(all_pairs_rmsd, 'complete');
    T = cluster(Z, 'maxclust', options.NumOfFinalists, 'Criterion', 'distance');
    cluster_size = accumarray(T, 1, [], @(x)numel(x));
    
    
    % Calculate the cluster index (matching entries in "Z") for each
    % conformation
    m = size(frames, 1);
    cluster_idx = [1:m]';
    
    for i = 1:size(Z, 1)
        if (all(accumarray(T, cluster_idx, [], @(x)all(x == x(1)))))
            break;
        end
        
        cluster_idx(cluster_idx == Z(i, 1)) = m + i;
        cluster_idx(cluster_idx == Z(i, 2)) = m + i;
    end
    
    if (0)
        % TODO: Ansorb small clusters into the level above in the hierarchy
        if (options.MinimalClusterSize > 0)
            for i = numel(cluster_size):-1:1
                if (cluster_size(i) < options.MinimalClusterSize)
                    % Join clusters
                    
                    cidx = cluster_idx(find(T == i, 1));
                    
                    higher_cluster_idx = m + find(any(Z(:, 1:2) == cidx, 2));
                end
            end
        end
    end
    
    
    results.Linkage.Hierarchy = Z;
    results.Linkage.Assignment = T;
    results.Linkage.ClusterSize = cluster_size;

    % Choose representative conformations
    if (1)
        clusters = accumarray(T, [1:numel(T)], [], @(x){x});
        
        % for each cluster, choose as a representative the conformation
        % resulting in widest (std) distribution of RMSD to it within the
        % cluster (that is - the most discerning conformation)
        
        representatives = cellfun(@(c)c(IndexOfMax(std(rmsd_matrix(c, c)))), clusters);
        results.Linkage.RepresentativeSimAndFrameIndex = ...
            representative_frame_indexes_from_all_windows(representatives, :);
    end
    
    %c = cophenet(Z,Y)
    
    %figure(4);
    %dendrogram(Z);
    %T = cluster(Z, 'cutoff', 1.1)
    %cutoff = median([Z(end-9,3) Z(end-8,3)]);
    %dendrogram(Z,'ColorThreshold',cutoff);
    %accumarray(T, representative_frame_indexes_from_all_windows(:, 3), [], @(x)min(x))
end

%%

    function [frame, frame_index] = ChooseRepresentativeFromWindow(window_frames)
        switch (options.WindowRepresentativeSelectionMethod)
            case 1
                frame_index = 1; % first frames
                
            case 2
                frame_index = floor(0.5 * size(window_frames, 1)); % middle frame
                
            case 3
                frame_index = size(window_frames, 1); % last frame
        end
        
        if (options.WindowRepresentativeSelectionMethod == 4)
            %%
            t = tic();
            [all_pairs_rmsd, ij] = AllPairsRmsd(window_frames);
            rmsd_matrix = zeros(size(window_frames, 1));    
            rmsd_matrix(sub2ind(size(rmsd_matrix), ij(:, 1), ij(:, 2))) = all_pairs_rmsd;
            rmsd_matrix(sub2ind(size(rmsd_matrix), ij(:, 2), ij(:, 1))) = all_pairs_rmsd;
            
            frame_index = IndexOfMax(std(rmsd_matrix));
            toc(t);
        end
        
        frame = window_frames(frame_index, :, :);
    end

    function [i] = IndexOfMax(x)
        [~, i] = max(x);
    end

end