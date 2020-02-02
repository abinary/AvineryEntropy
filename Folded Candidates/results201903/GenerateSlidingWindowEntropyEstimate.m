function [entropy_estimates, results] = GenerateSlidingWindowEntropyEstimate(file_list, are_files_continuous, cg_levels, window_length, stride, use_worker_pool)

if (~exist('use_worker_pool', 'var'))
    use_worker_pool = false;
end

VerySmallNumber = 1e-10;

if (isempty(file_list))
    return;
end

if (isempty(are_files_continuous))
    are_files_continuous = false;
end

if (numel(cg_levels) ~= 1)
    error('GenerateSlidingWindowEntropyEstimate must be given a single coarse-graining level value');
end

if (numel(window_length) ~= 1)
    error('GenerateSlidingWindowEntropyEstimate must be given a single windows length value');
end

use_clustering = false;
if (cg_levels < 0)
    cg_levels = abs(cg_levels);
    use_clustering = true;
end

addpath('c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\Analysis Scripts');
addpath('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\Common Code');

% Load .Net engines
NET.addAssembly('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\PhysiCasino\PhysiCasino\bin\x64\Release\PhysiCasino.exe');

if (use_worker_pool)
    % pctRunOnAll
    PrepWorkersPool();
end

if (are_files_continuous)
    num_frames_kept_between_files = window_length - stride;
    results  = [];
else
    num_frames_kept_between_files = 0;
    results  = {};
end

wb = waitbar(0, 'Processing ...', 'Name', 'Sliding entropy estimate progress', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1); delete(gcbf)', ...
    'CloseRequestFcn', 'delete(gcbf)');
setappdata(wb, 'canceling', 0);

for file_index = 1:numel(file_list)
    waitbar(file_index/numel(file_list), wb, sprintf('Processing file %d/%d', file_index, numel(file_list)));
    
    % Check for clicked Cancel button
    if (~ishandle(wb) || getappdata(wb, 'canceling'))
        break
    end
    
    % If the files contain a continuous simulation, keep the last frames
    % and add the loaded file
    if (num_frames_kept_between_files ~= 0 && file_index > 1)
        frames = frames(end-num_frames_kept_between_files+1:end, :, :);
        frames = [frames; LoadTrajectoryFileAsDihedralAngles(file_list{file_index})];
    else
        [frames, nd] = LoadTrajectoryFileAsDihedralAngles(file_list{file_index});
    end
    
    % Automatically detect limits of data
    if (~exist('L', 'var') || isempty(L))
        [L, U] = BestMatchLimits(frames, nd);
    end
    
    if (file_index == 1)
        analysis = AnalysisPlan();
        analysis.UseParFor = use_worker_pool;
        
        % Set coarse-graining
        if (use_clustering)
            analysis.Parameters.CoarseGrain_Level = cg_levels;
        else
            analysis.Parameters.CoarseGrain_NumOfStates = cg_levels;
        end
        
        %analysis.Parameters.CoarseGrain_NumOfStates = cg_levels; % Coarse graining levels and range, for the entropy estimate
        analysis.Parameters.CoarseGrain_RangeMin = L;
        analysis.Parameters.CoarseGrain_RangeMax = U;
        analysis.Parameters.CoarseGrain_NumOfBits = 8;

        % Calibrate with 10,000 configurations, to roughly get the asymptotic value
        analysis.Parameters.MinimalNumConfigurationsForCalibration = 1e4;
        %analysis.Parameters.MinimalLengthForCalibration = 1e7;
        
        % Convert starting indexes into batches of frames
        % This has to be added every time, instead of using the same "analysis"
        % instance, because the referenced "frames" variable changes content
        %         analysis.AddAdvancedConversion(...
        %             @(a, s)a.Parameters.Frames(s - 1 + [1:window_length], :, :), 'GetWindow');
        %analysis.AddConversion(@(s)frames(s - 1 + [1:window_length], :, :));
        
        if (use_clustering)
            analysis.AddAdvancedConversion(@ClusterAssignmentCoarseGraining);
        else
            analysis.AddConversion(@LinearizeConfigurations);
        end
        
        % Add entropy estimate
        analysis.AddAdvancedAnalysis(@EstimateEntropyForCoarseGrainedConfigurations);
        
        % If parallel, force calibration first
        if (analysis.UseParFor)
            analysis.UseParFor = false;
            analysis.Analyze(1);
            analysis.UseParFor = true;
        end
        
        calibrated_analysis = analysis;
    end
    
    % Create the analysis manager and copy calibrated parameters
    analysis = AnalysisPlan();
    analysis.UseParFor = use_worker_pool;
    analysis.Parameters = calibrated_analysis.Parameters;
    
    % To be saved into the results
    analysis.Parameters.FilePath = file_list{file_index};
    analysis.Parameters.FileIndex = file_index;
    
    %analysis.Parameters.Frames = frames;
    
    % Convert starting indexes into batches of frames
    % This has to be added every time, instead of using the same "analysis"
    % instance, because the referenced "frames" variable changes content
    %         analysis.AddAdvancedConversion(...
    %             @(a, s)a.Parameters.Frames(s - 1 + [1:window_length], :, :), 'GetWindow');
    %analysis.AddConversion(@(s)frames(s - 1 + [1:window_length], :, :));
    
    if (use_clustering)
        analysis.AddAdvancedConversion(@ClusterAssignmentCoarseGraining);
    else
        analysis.AddConversion(@LinearizeConfigurations);
    end
    
    analysis.AddAdvancedConversion(@CoarseGrain);

    % Add entropy estimate
    %analysis.AddAdvancedAnalysis(@EstimateEntropyForCoarseGrainedConfigurations);
    analysis.AddAdvancedAnalysis(@SlidingCompressedLength);
    analysis.AddAdvancedAnalysis(@ScaleToEntropy);
    
    % Run on all windows
    window_start_indexes = 1:stride:(size(frames, 1) - window_length + 1);
    window_start_indexes = arrayfun(@(x)x, window_start_indexes, 'UniformOutput', false); % Causes analysis on each index separately
    %results_for_single_file = analysis.Analyze(window_start_indexes);
    results_for_single_file = analysis.Analyze(frames);
    
    % Aggregate results
    if (are_files_continuous)
        results = [results; results_for_single_file];
    else
        results{end + 1} = results_for_single_file;
    end
end

delete(wb);

% Extract the entropy estimates from the results structure
if (are_files_continuous)
    entropy_estimates = arrayfun(@(x)x.EntropyEstimate, results);
else
    entropy_estimates = cellfun(@(c)c.EntropyEstimate, results, 'UniformOutput', false);
end

    function [cgConfigurations] = CoarseGrain(analysisObj, configurations)
        numOfStates = double(analysisObj.Parameters.CoarseGrain_NumOfStates);
        
        if (isinteger(configurations)) % Skip coarse graining if it's already coarse grained
            cgConfigurations = configurations;
        elseif (isfield(analysisObj.Parameters, 'CoarseGrain_RangeMin'))
            %%
            [cgConfigurations, range] = CoarseGrainInGivenRange(configurations, ...
                analysisObj.Parameters.CoarseGrain_RangeMin, analysisObj.Parameters.CoarseGrain_RangeMax,...
                numOfStates, false);
        end
    end

    function [data] = WriteToByteArray(cgConfigurations)
        ms = System.IO.MemoryStream();
        bs = PhysiCasino.BitWriter(ms);
        PhysiCasino.FloatingPointProcessor.WriteBinaryIntegers(bs, cgConfigurations, 8);
        data = ms.ToArray();
    end

    function [] = SlidingCompressedLength(analysisObj, data)
        data = WriteToByteArray(data);
        dof = size(frames, 2) * size(frames, 3);
        rawSize = dof * window_length;%double(data.Length);
        %c = PhysiCasino.CompressionProcessor.CompressAndReportSize(data, 26);
        c = PhysiCasino.CompressionProcessor.SlidingCompressedSize(data, stride * dof, window_length * dof, 26);
        compressedSize = double(c);
        incompressibility = compressedSize / rawSize;
        
        analysisObj.Results.SizeRaw = rawSize;
        analysisObj.Results.SizeCompressed = compressedSize;
        analysisObj.Results.Incompressibility = incompressibility;
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
                
                pixel_area = (4 * pi^2) / numel(Z);
                
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
