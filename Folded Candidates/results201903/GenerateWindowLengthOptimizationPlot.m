function [] = GenerateWindowLengthOptimizationPlot(cg_levels, file_list, num_windows, window_lengths, seed)

VerySmallNumber = 1e-10;

if (nargin == 5)
    rng(seed); % Fix the lottery (for debug)
end

use_clustering = false;

if (num_windows < 0 || cg_levels < 0)
    num_windows = abs(num_windows);
    cg_levels = abs(cg_levels);
    use_clustering = true;
end

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
if (~isempty(current_pool))
    analysis.UseParFor = true;
    
    PrepWorkersPool();
end

if (use_clustering)
    analysis.Parameters.CoarseGrain_Level = cg_levels;
else
    analysis.Parameters.CoarseGrain_NumOfStates = cg_levels;
end

%analysis.BranchParameter('CoarseGrain_NumOfStates', cg_levels);
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
results = analysis.Analyze(window_frames(1));

results = analysis.Analyze(window_frames);

save('window_length_opt_results.mat', 'max_window_length', 'cg_levels', 'file_list', ...
    'num_windows', 'window_lengths', 'results');

%%

if (0)
    %%
    load('window_length_opt_results.mat'); % for debug
end

ok = false;
if (isfield(results(1), 'EntropyEstimateUnadjusted'))
    x = arrayfun(@(x)x.EntropyEstimateUnadjusted, results);
    
    if (1 && isfield(results(1), 'EntropyBinWidthContribution'))
        fprintf('Using the "EntropyBinWidthContribution" field \n');
        b = arrayfun(@(x)x.DOF * x.EntropyBinWidthContribution, results);
        ok = true;
    elseif (isfield(results(1), 'EntropyBaselinePerDOF'))
        fprintf('Using the "EntropyBaselinePerDOF" field \n');
        b = arrayfun(@(x)x.DOF * x.EntropyBaselinePerDOF, results);
        ok = true;
    elseif(1)
        fprintf('Using UNadjusted entropy estimates - log(n_s) + logV \n');
        %logV = arrayfun(@(x)sum(log(x.Parameters.CoarseGrain_RangeMax - x.Parameters.CoarseGrain_RangeMin)), results);
        b = arrayfun(@(x)- x.DOF * log(x.Parameters.CoarseGrain_NumOfStates), results);
        %b = b + logV;
        ok = true;
    end
end

if(~ok)
    fprintf('Using adjusted entropy estimates + logV \n');
    x = arrayfun(@(x)x.EntropyEstimate, results);
    logV = arrayfun(@(x)sum(log(x.Parameters.CoarseGrain_RangeMax - x.Parameters.CoarseGrain_RangeMin)), results);
    b = logV;
end

%b = 0;
x = x + b;
%x = squeeze(x);
y = x - min(x, [], 2);

f = figure();
subplot(1, 2, 1);
p = plot(window_lengths, y, 'LineWidth', 2);
ax = gca();
ax.FontSize = 12;
ax.FontWeight = 'bold';

xlabel('Window Length (frames)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$S_\mathrm{A}$ ($K_\mathrm{B}$)', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
title(sprintf('%d windows', num_windows), 'FontSize', 16, 'FontWeight', 'bold');

%y_kb = y;
%[~, sel] = sort(rand([1 size(y, 1)])); % randomly select windows
%sel = sel(1:100);
%y_kb = y_kb(sel, :)';

subplot(1, 2, 2);
%y = y ./ max(y, [], 2);
y = y ./ max(y(:));
p = plot(window_lengths, y, 'Color', 0.5 * [1 1 1], 'LineWidth', 2);

hold on;
mean_plot = plot(window_lengths, mean(y, 1), 'Color', 'b', 'LineWidth', 2);

x = window_lengths(:);
y = mean(y, 1);
y = y(:);

model = 'a*exp(-x/b) + c';
fo = fitoptions(model);
fo.Lower = [0 1 0];
fo.Upper = [2 max_window_length * 0.1 1];
fo.StartPoint = [1 100 0];
f1 = fit(x, y, model, fo);

model = 'a*exp(-x/b) + c*exp(-x/d) + e';
fo = fitoptions(model);
fo.Lower = [0 1 0 1 0];
fo.Upper = [2 max_window_length 2 max_window_length 1];
fo.StartPoint = [1 max_window_length * 0.1 0 max_window_length * 0.5 0];
f2 = fit(x, y, model, fo);

f1p = plot(x, feval(f1, x), ':', 'Color', 'g', 'LineWidth', 2);
f2p = plot(x, feval(f2, x), '--', 'Color', 'm', 'LineWidth', 2);

hold off;

ax = gca();
ax.FontSize = 12;
ax.FontWeight = 'bold';

xlabel('Window Length (frames)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Normalized $S_\mathrm{A}$', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'latex');
title('Normalized', 'FontSize', 16, 'FontWeight', 'bold');

legend([mean_plot f1p f2p], ...
    {'Mean convergence', ...
    sprintf('single exp + const (t = %2.0f)', f1.b), ...
    sprintf('2 exp + const (t1 = %2.0f , t2 = %2.0f)', f2.b, f2.d)}, 'Interpreter', 'tex');

1;

%%

    function [window] = TruncateWindow(analysis, window)
        window = window(1:analysis.Parameters.WindowLength, :, :);
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
