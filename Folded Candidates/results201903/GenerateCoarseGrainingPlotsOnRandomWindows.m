function [] = GenerateCoarseGrainingPlotsOnRandomWindows(cg_levels, file_list, num_windows, window_lengths, seed)

VerySmallNumber = 1e-10;

if (nargin == 5)
    rng(seed); % Fix the lottery (for debug)
end

use_clustering = num_windows < 0; num_windows = abs(num_windows);

file_indexes = randi(numel(file_list), [num_windows 1]);
num_in_each = accumarray(file_indexes, 1);

window_frames = {};

max_window_length = max(window_lengths);



for i = 1:numel(num_in_each)
    if (num_in_each(i) == 0)
        continue;
    end
    
    [frames, nd] = LoadTrajectoryFileAsDihedralAngles(file_list{i});
    
    % Automatically detect limits of data
    if (~exist('L', 'var') || isempty(L))
        [L, U] = BestMatchLimits(frames, nd);
    end
    
    max_window_start = size(frames, 1) - max_window_length;
    
    if (0) % for debug
        %%
        x = frames(:, :, 2); x = x(:);
        [n, edges] = histcounts(x);
        centers = 0.5 * (edges(1:end-1) + edges(2:end));
        n = n ./ max(n);
        
        figure(1);
        bar(centers, n);
        
        [edges] = CalculateEqualProbabilityEdges(x(:), 5);
        [edges] = CalculateEqualProbabilityEdges(x(:), 5);
        hold on;
        for e = edges(2:end-1)
            plot([e e], [0 1], '--g');
        end
        hold off;
        
        %%
        histogram(x, edges);
    end
    
    for j = 1:num_in_each(i)
        windows_start_pos = randi(max_window_start);
        window_frames{end + 1} = frames(windows_start_pos - 1 + [1:max_window_length], :, :);
        fprintf('File %03d Window start: %06d \r\n', i, windows_start_pos);
    end
end



addpath('c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\Analysis Scripts');
addpath('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\Common Code');

% Load .Net engines
NET.addAssembly('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\PhysiCasino\PhysiCasino\bin\x64\Release\PhysiCasino.exe');


analysis = AnalysisPlan();

current_pool = gcp('nocreate');
if (0 && ~isempty(current_pool)) % code disabled
    analysis.UseParFor = true;
    
    PrepWorkersPool();
end

if (use_clustering)
    analysis.BranchParameter('CoarseGrain_Level', cg_levels);
else
    analysis.BranchParameter('CoarseGrain_NumOfStates', cg_levels);
end

%analysis.Parameters.WindowLength = 1;
analysis.BranchParameter('WindowLength', window_lengths);

L = zeros(size(frames, 2), 1) + L;
U = zeros(size(frames, 2), 1) + U;

analysis.Parameters.CoarseGrain_RangeMin = (L - VerySmallNumber);
analysis.Parameters.CoarseGrain_RangeMax = (U + VerySmallNumber);
analysis.Parameters.CoarseGrain_NumOfBits = 8;

% Calibrate with 10,000 configurations, to roughly get the asymptotic value
analysis.Parameters.MinimalNumConfigurationsForCalibration = 1e3;

analysis.AddAdvancedConversion(@TruncateWindow, 'truncate window');

if (use_clustering)
    analysis.AddAdvancedConversion(@ClusterAssignmentCoarseGraining);
else
    analysis.AddConversion(@LinearizeConfigurations);
end

analysis.AddAdvancedAnalysis(@EstimateEntropyForCoarseGrainedConfigurations);

% Cause calibration first
if (analysis.UseParFor)
    results = analysis.Analyze(window_frames(1));
end

% Branch on window lengths
results = analysis.Analyze(window_frames);

1;
%%
save('coarse_graining_opt_results.mat', 'cg_levels', 'file_list', ...
    'num_windows', 'window_lengths', 'results');

%%
was_calculated = false;
if (isfield(results(1), 'EntropyEstimateUnadjusted'))
    x = arrayfun(@(x)x.EntropyEstimateUnadjusted, results);
    
    if (1 && isfield(results(1), 'EntropyBinWidthContribution'))
        fprintf('Using the "EntropyBinWidthContribution" field \n');
        b = arrayfun(@(x)x.DOF * x.EntropyBinWidthContribution, results);
        was_calculated = true;
    elseif (isfield(results(1), 'EntropyBaselinePerDOF'))
        fprintf('Using the "EntropyBaselinePerDOF" field \n');
        b = arrayfun(@(x)x.DOF * x.EntropyBaselinePerDOF, results);
        was_calculated = true;
    elseif(0)
        fprintf('Using UNadjusted entropy estimates - log(n_s) + logV \n');
        logV = arrayfun(@(x)sum(log(x.Parameters.CoarseGrain_RangeMax - x.Parameters.CoarseGrain_RangeMin)), results);
        b = arrayfun(@(x)- x.DOF * log(x.Parameters.CoarseGrain_NumOfStates), results);
        b = b + logV;
        was_calculated = true;
    elseif(1)
        fprintf('Using UNadjusted entropy estimates - log(n_s) \n');
        b = arrayfun(@(x)- x.DOF * log(x.Parameters.CoarseGrain_NumOfStates), results);
        was_calculated = true;
    end
end

if (~was_calculated)    
    fprintf('Using adjusted entropy estimates + logV \n');
    x = arrayfun(@(x)x.EntropyEstimate, results);
    logV = arrayfun(@(x)sum(log(x.Parameters.CoarseGrain_RangeMax - x.Parameters.CoarseGrain_RangeMin)), results);
    b = logV;
    was_calculated = true;
end

%b = 0;
x = x + b;
%x = arrayfun(@(x)x.EntropyEstimate, results);

f = figure(1);

w = ceil(sqrt(numel(window_lengths)));
h = ceil(numel(window_lengths) / w);

for i = 1:size(x, 3)
    subplot(h, w, i);
    p = plot(cg_levels, x(:, :, i), 'LineWidth', 2);
    ax = gca();
    ax.FontSize = 12;
    ax.FontWeight = 'bold';
    
    xlabel('# CG levels', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('$S_\mathrm{A}$ ($K_\mathrm{B}$)', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
    title(sprintf('Window Length: %d', window_lengths(i)), 'FontSize', 16, 'FontWeight', 'bold');
end

if (0)
    %%
    x = arrayfun(@(x)x.EntropyEstimateUnadjusted, results);
    b = arrayfun(@(x)- x.DOF * log(x.Parameters.CoarseGrain_NumOfStates), results);
    x = x + b;
    mean_estimate = squeeze(mean(x, 1));
    %data = [cg_levels' mean_estimate];
    data = [cg_levels' mean_estimate' x'];
    open('data'); % for copying to OriginPro
end

1;

if (0)
    SetMyDefaultFigureSettings();
    
    fig_scale_factor = 1.5;
    
    f = gcf();
    ax = gca();
    ax.Units = 'centimeters';
    ax.Position = [2 1.5 2*[4.3 4.3]] * fig_scale_factor;
    
    f.Color = [1 1 1];
    f.Units = 'centimeters';
    %f.Position = [f.Position(1:2) ax.Position(3:4) + [2.5 2.5] * fig_scale_factor];
    f.Position = [[5 8] ax.Position(3:4) + [2.5 2.5] * fig_scale_factor];
end

if (ndims(x) == 2)
    %%
    [x_min, I] = min(x, [], 2);
    I = cg_levels(I);
    
    K1 = 1 + mod(I - 1, 16);
    K2 = ceil(I / 16);
    %(K2-1)*16 + (K1-1) + 1
    
    hold on;
    plot(I, x_min, '*r', 'LineWidth', 3);
    hold off;
    
end

if (1)
    %%
    figure(4);
    histogram2(K2, K1, [1:17], [1:17]);

    %xl = xlabel('\phi coarse graining level');
    %yl = ylabel('\psi coarse graining level');
    xl = xlabel('\phi');
    yl = ylabel('\psi');
    zl = zlabel('Count');

    SetMyDefaultFigureSettings();
    %axis equal;
    
    ax = gca();
    ax.XTick = [2:2:16];
    ax.YTick = [2:2:16];
    
    %ax.CameraPosition = [89.5451  -89.4075   51.2129];
    %ax.CameraPosition = [106.7870  -85.7641   27.3336];
    %ax.CameraPosition = [110.8714  -79.8699   31.2027];
    %ax.CameraTarget = [9.0000    9.0000    6.5000];

    xl.Rotation = -13;
    xl.Units = 'normalized';
    xl.Position = [0.1 0.07 0];

    yl.Rotation = 8;
    yl.Units = 'normalized';
    yl.Position = [0.85 0.03 0];
    
    zl.Units = 'normalized';
    zl.Position = [ -0.075    0.5082         0];

    ax.Position = [0.12 0.13 0.8 0.85];
    
    % 

    %%
    K1 = 1 + mod(cg_levels - 1, 16);
    K2 = ceil(cg_levels / 16);
    num_levels = K1 .* K2;
    %x = arrayfun(@(x)x.EntropyEstimateUnadjusted, results);
    b = -results(1).DOF * log(num_levels);
    y = arrayfun(@(x)x.EntropyEstimateUnadjusted, results) + b;
    %y = arrayfun(@(x)x.EntropyBinWidthContribution, results);
    
    y = x;
    y = permute(reshape(y, [size(y, 1) 16 16]), [2 3 1]);
    figure(3);
    imagesc(mean(y, 3));
    %imagesc(min(y, [], 3));
    %imagesc(max(y-mean(y, 3), [], 3));
    axis xy
    hold on;
    %plot(K2, K1, '*r', 'LineWidth', 3);
    hold off;
    
    xlabel('\phi coarse graining level');
    ylabel('\psi coarse graining level');

    SetMyDefaultFigureSettings();
    %axis equal;
    
    ax = gca();
    ax.XTick = [2:2:16];
    ax.YTick = [2:2:16];
    
    cb = colorbar();
    cb.Ticks = [-150:25:-50];
    
    %text(17.6, 15.8, '($k_\mathrm{B}$)', 'FontSize', 18, 'FontWeight', 'bold', 'interpreter', 'latex');
    text(17.6, 15.8, '(k_B)', 'FontSize', 18, 'FontWeight', 'bold');
    %text(17.6, 17.2, '(k_B)', 'FontSize', 18, 'FontWeight', 'bold');
    
end
%%
%SaveFig('CG optimization regular 1-35.jpg');

f = figure(2);

x = x - mean(x, 1);

w = ceil(sqrt(numel(window_lengths)));
h = ceil(numel(window_lengths) / w);

for i = 1:size(x, 3)
    subplot(h, w, i);
    p = plot(cg_levels, x(:, :, i), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    ax = gca();
    ax.FontSize = 12;
    ax.FontWeight = 'bold';
    
    hold on;
    plot(cg_levels, std(x(:, :, i), 1, 1), '-g', 'LineWidth', 3);
    hold off;
    
    xlabel('# CG levels', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('STD of $S_\mathrm{A}$ ($K_\mathrm{B}$)', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
    title(sprintf('Window Length: %d', window_lengths(i)), 'FontSize', 16, 'FontWeight', 'bold');
end
1;

%%

    function [window] = TruncateWindow(analysis, window)
        window = window(1:analysis.Parameters.WindowLength, :, :);
    end

    function [] = GenerateEdges(analysis, data)
        if (isfield(analysis.Parameters, 'CoarseGrain_Edges'))
            return;
        end
        
        cg = analysis.Parameters.CoarseGrain_NumOfStates;
        
        x = frames(:, :, 1); x = x(:);
        [edges1] = CalculateEqualProbabilityEdges(x(:), cg);
        
        x = frames(:, :, 2); x = x(:);
        [edges2] = CalculateEqualProbabilityEdges(x(:), cg);
        
        analysis.Parameters.CoarseGrain_Edges = [edges1; edges2];
    end

    function [data] = ClusterAssignmentCoarseGraining(analysis, data)
        cg = analysis.Parameters.CoarseGrain_Level;
        
        if (~isfield(analysis.Parameters, 'AssignmentFunc'))
            options = GenerateDihedralAngleClusters();
            
            options.ShouldAdjustNearestNeighborsDistance = true;
            %options.Method = 'PAM';%options.MethodsList{1};
            %options.Method = 'linkage';
            %options.SubMethod = 'complete';
            %options.Method = 'decomp2';
            options.Method = 'dummy';
            
            options.DecompCondDir = -1;
            options.DecompMaxAspectRatio = 1;
            
            options.HistogramGraining = 101;%2^10;
            options.RelativeProbabilityThreshold = 1e-3;
            
            [assignment_func, Z, p] = GenerateDihedralAngleClusters(file_list, cg, options);
            Z = Z';

            if (1)
                num_bins = max(Z(:));
                entropy_bin_width_contributions = zeros(1, num_bins);
                
                %pixel_area = (4 * pi^2) / numel(Z);
                pixel_area = 1.0 / numel(Z);
                
                for z_idx = 1:num_bins
                    w = nnz(Z == z_idx);
                    assert(w > 0);
                    
                    entropy_bin_width_contributions(z_idx) = log(w * pixel_area);
                    %entropy_bin_width_contributions(z_idx) = log(w);
                end
            end
            
            entropy_baseline = dot(p, entropy_bin_width_contributions);
            
            analysis.Parameters.AssignmentFunc = assignment_func;
            analysis.Parameters.AssignmentProbabilities = p;
            analysis.Parameters.EntropyBaselinePerDOF = entropy_baseline;
            analysis.Parameters.EntropyBinWidthFactors = entropy_bin_width_contributions;
            analysis.Parameters.CoarseGrain_NumOfStates = max(Z(:));
            
            if (0)
                %%
                figure(3);
                imagesc(Z); axis xy;
                colormap(jet(256));
                drawnow();
            end
        else
            assignment_func = analysis.Parameters.AssignmentFunc;
        end
        
        num_bins = analysis.Parameters.CoarseGrain_NumOfStates;
        data = assignment_func(data);
        assignments_stats = accumarray(data(:), 1, [num_bins 1], @sum);
        
        data = data - 1; % shift indexes to 0
        
        analysis.Results.CoarseGrain_AssignmentStats = assignments_stats';
        
        p = assignments_stats ./ numel(data);
        analysis.Results.EntropyBinWidthContribution = dot(p, analysis.Parameters.EntropyBinWidthFactors);
        analysis.Results.EntropyBaselinePerDOF = analysis.Parameters.EntropyBaselinePerDOF;
        analysis.Results.CoarseGrain_NumOfStates = analysis.Parameters.CoarseGrain_NumOfStates;
    end

end
