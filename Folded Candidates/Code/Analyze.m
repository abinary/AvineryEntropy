function [] = Analyze(src)

if (~isempty(src.UserData))
    addpath('d:\Files\Workspace\Sosnick\');
    
    dcm = datacursormode(src);
    
    state = src.UserData;
    %index = get(event_obj, 'DataIndex'); % An undocumented property?
    index = state.SelectedIndex;
    
    sim_num = state.Data.SimNums(index);
    start_frame_index = state.Data.StartFrames(index);
    frame_indexes = start_frame_index - 1 + [1:state.WindowLength];
    trajectory_frame_indexes = frame_indexes(2:2:end) / 2; % Frames in the trajectory are actually saved every other one
    
    if (length(state.Trajectories) < sim_num || isempty(state.Trajectories{sim_num}))
        fprintf('Loading trajectories for sim# %d\n', sim_num);
        trajectory = state.Loaders.LoadTrajectoryFileCAlpha(sim_num);
        state.Trajectories{sim_num} = trajectory;
    else
        trajectory = state.Trajectories{sim_num};
    end
    
    rmsd = state.Loaders.LoadRmsd(sim_num);
    
    %%
    frames_in_window = trajectory.frames(trajectory_frame_indexes, :, :);
    
    %%
    [I, J] = meshgrid(1:size(frames_in_window, 1));
    which = J > I;
    I = I(which);
    J = J(which);
    
    rmsd_within_window = arrayfun(@(i)...
        CalculateRmsd(squeeze(frames_in_window(I(i), :, :)), squeeze(frames_in_window(J(i), :, :))), ...
        1:numel(I));
    
    %%
    rmsd_image = zeros(size(frames_in_window, 1));
    
    rmsd_image(sub2ind(size(rmsd_image), I, J)) = rmsd_within_window;
    rmsd_image = rmsd_image + rmsd_image';
    
    figure(3);
    imagesc(rmsd_image);
    axis xy;
    colorbar;
    title('RMSDs within selected window');
    
    %% Calculate entropy of each RMSD histogram
    rmsd_binned = floor(rmsd_image);
    [n] = histcounts2([1:size(rmsd_image, 1)]' + 0 .* rmsd_binned, rmsd_binned, -0.5:size(rmsd_image, 1)+0.5, 0.5:max(rmsd_binned(:))+0.5);
    n = n ./ size(rmsd_image, 1);

    s = -n .* log2(n + 1e-10);
    s = sum(s, 2);

    figure(2);
    imagesc(n');
    %bar(s);
    plot(n(s < 1.8, :)');
    
    
    %%
    figure(2);
    histogram(rmsd_within_window, 100);
    title('RMSD between all frames within window');
    xlim([0 15]);
    
    %%
    rmsd_within_window_matrix = zeros(size(frames_in_window, 1));
    rmsd_within_window_matrix(which) = rmsd_within_window;
    rmsd_within_window_matrix = rmsd_within_window_matrix + rmsd_within_window_matrix';
    
    if (0)
        figure(4);
        imagesc(rmsd_within_window_matrix);
        axis xy;
        
        figure(5);
        bar(mean(rmsd_within_window_matrix, 1));
    end
    
    %%
    [~, representative_frame_index] = min(mean(rmsd_within_window_matrix, 1));
    representative_frame_index = trajectory_frame_indexes(representative_frame_index);
    representative_frame = squeeze(trajectory.frames(representative_frame_index, :, :));
    
    if (0)
        rmsd2 = arrayfun(@(i)...
            CalculateRmsd(representative_frame, squeeze(trajectory.frames(i, :, :))), ...
            1:size(trajectory.frames, 1));
    else
        if (numel(state.Trajectories) < 27)
            state.Trajectories{27} = [];
        end
        
        for i = 1:27
            if (isempty(state.Trajectories{i}))
                fprintf('Loading trajectory #%d \n', i);
                state.Trajectories{i} = state.Loaders.LoadTrajectoryFileCAlpha(i);
            else
                1;
            end
        end

        if (1)
            display('Calculating RMSD for all frames');
            
            all_frames = arrayfun(@(i)state.Trajectories{i}.frames, 1:numel(state.Trajectories), 'UniformOutput', false);
            all_frames = vertcat(all_frames{:});
            
            rmsd2 = arrayfun(@(i)...
                CalculateRmsd(representative_frame, squeeze(all_frames(i, :, :))), ...
                1:size(all_frames, 1));
        else
            rmsd2 = [];
        end
    end
    
    display('Plotting now...');
    
    %%
    f = figure(6);
    f.Position = [354   141   896   420];
    
    edges = linspace(0, 30, 301);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    
    n = histcounts(rmsd2, edges);
    h = area(centers, n);
    h.FaceAlpha = 0.5;
    %histogram(rmsd2, 200);
    title('RMSD of all frames to selected frame');
    xlim([0 25]);
    
    %
    hold on;
    plot(centers, n, '-b', 'LineWidth', 2);

    
    n = histcounts(state.Data.Rmsd, edges);
    %histogram(state.Data.Rmsd, 200)
    
    h = area(centers, n);
    h.FaceAlpha = 0.2;
    h.FaceColor = [0 0 0];
    %h.EdgeColor = [1 0 0];
    p = plot(centers, n, '-k', 'LineWidth', 2);
    hold off;
    
    xlim([0 25]);
    ylim([0 3.5e4]);
    xlabel('RMSD', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('#', 'FontSize', 16, 'FontWeight', 'bold');
    title('RMSD of all frames vs. PDB');
    legend({'Selected frame', '', 'PDB', ''});

    
    %%
    figure(7);
    %plot(rmsd, rmsd2, '.');
    %h = histogram2(rmsd, rmsd2, 100);
    h = histogram2(state.Data.Rmsd, rmsd2, 300);
    h.FaceColor = 'flat';
    xlabel('RMSD(pdb)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('RMSD(picked frame)', 'FontSize', 16, 'FontWeight', 'bold');
    title('RMSD(pdb) vs. RMSD(picked frame)');
    view(0, 90);
    
    hold on;
    plot3([0 20], [0 20], 1e3 * [1 1], '--k', 'LineWidth', 2);
    hold off;
    
    
    %% RMSD to mean conformation
    if (0)
        mean_conformation = squeeze(median(frames_in_window, 1));
        rmsd_within_window_to_mean = arrayfun(@(i)...
            CalculateRmsd(mean_conformation, squeeze(frames_in_window(i, :, :))), ...
            1:size(frames_in_window, 1));
        
        figure(3);
        histogram(rmsd_within_window_to_mean, 20);
        title('RMSD within window to average conformation');
        xlim([0 15]);
        
    end
    
    
end
