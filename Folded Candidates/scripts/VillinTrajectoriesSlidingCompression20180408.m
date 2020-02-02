
% Protein folding kinetics and thermodynamics from atomistic simulation
% (2012)
% doi: 10.1073/pnas.1201811109
% http://www.pnas.org/content/109/44/17845.full

cd('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\Analysis Scripts');

addpath('../Common Code');
addpath('../Common Code/splinefit');

% Load .Net engines
NET.addAssembly('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\PhysiCasino\PhysiCasino\bin\x64\Release\PhysiCasino.exe');
%NET.addAssembly('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\PhysiCasino\DcdSharp\bin\x64\Release\DcdSharp.exe');

outputRoot = '../PhysiCasino Output/';
timeStampString = datestr(now, 'yyyymmdd-hhMMss');

%%
%simulationPrefix = 'pnas2012-1yrf-WT-F10L-0';
simulationPrefix = 'pnas2012-2f4k-360K';
simulationPrefix = [simulationPrefix '-protein'];

windowLength = 2000;
windowStride = 0.1 * windowLength;

VerySmallNumber = 1e-10;

dt = 200;

%%
dataRoot = 'c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\Data\Villin trajectories\';
%dataRoot = 'd:\Files\Workspace\Villin\';
dataFolderUpper = [dataRoot 'DESRES-Trajectory_' simulationPrefix '\' ];
dataFolder = [dataRoot 'DESRES-Trajectory_' simulationPrefix '\' simulationPrefix '\'];

%%
dcdTable = readtable([dataFolder simulationPrefix '_times.csv']);
try
    dcdTimes = double(dcdTable.Var1);
    dcdFilenames = dcdTable.Var2;
catch
    %dcdTimes = double(dcdTable.StartTime_ps_);
    [dcdFilenames, firstOfEach, dcdTimeAssociation] = unique(dcdTable.File);
    dcdFilenames = cellfun(@(f)f(max(find(f=='/'))+1:end), dcdFilenames, 'UniformOutput', false);
    dcdTimes = {};
    startTimes = dcdTable.StartTime_ps_(firstOfEach);
    
    for i = 1:numel(dcdFilenames)
        which = (dcdTimeAssociation == i);
        numOfFrames = dcdTable.Count(which);
        timeSteps = cellfun(@(x)str2double(x), dcdTable.TimeBetweenFrames_ps_(which));
        
        frameTimes = arrayfun(@(i)ones(numOfFrames(i), 1)*timeSteps(i), 1:numel(numOfFrames), 'UniformOutput', false);
        frameTimes = vertcat(frameTimes{:});
        dcdTimes{i} = startTimes(i) + cumsum([0; frameTimes(2:end)]);
    end
end

states = LoadFilePattern([dataFolderUpper 'States_*.dat'], '-ascii');
%energies = LoadFilePattern([dataFolderUpper 'Energies_*.dat'], '-ascii');

shouldProcessParallel = true;

if (shouldProcessParallel)
    % pctRunOnAll
    PrepWorkersPool();
end



%%

collectedResults = {};
collectedFrames = [];
collectedFramesTimes = [];

for fileIndex = 1:numel(dcdFilenames)
    %%
    fprintf('File: %d/%d\n', fileIndex, numel(dcdFilenames));
    
    data = load([dataFolder dcdFilenames{fileIndex} '.angles.mat']);
    frames = data.angleFrames;
    clear data;
    
    size(collectedFrames, 1)
    
    % Read new frames
    numFramesRead = size(frames, 1);
    
    if (isnumeric(dcdTimes))
        currentSimulationTime = dcdTimes(fileIndex);
        frameTimes = currentSimulationTime + dt*[0:numFramesRead-1]';
        currentSimulationTime = currentSimulationTime + dt*numFramesRead;
    else
        frameTimes = dcdTimes{fileIndex};
        currentSimulationTime = frameTimes(end);
    end
    
    
%     whichFramesAreFolded = interp1(t0 + states(:, 1) * dt, 1-states(:, 2), frameTimes);
%     whichFramesAreFolded = (whichFramesAreFolded == 1);
%     
%     if (0) % folded/unfolded derived from our partitioning
%         whichFramesAreFolded = interp1(windowMeanTimes, double(~aboveTheLine), frameTimes);
%         whichFramesAreFolded = (whichFramesAreFolded == 1);
%     end
    1;
    
    %% Remove frames sampled too frequently
    which = true;
    while (any(which))
        deltaTimes = diff(frameTimes);
        which = (deltaTimes(1:end-1) > 199) & (deltaTimes(2:end) <= 199);
        which = 1 + find(which);
        
        frameTimes(which) = [];
        frames(which, :, :) = [];
    end
    
    %% Collect frames
    collectedFrames = [collectedFrames; frames]; % Accoring to correlation times, >30-50 should be good, so 200 just in case
    collectedFramesTimes = [collectedFramesTimes; frameTimes];
    
    % Process if reached 1M frames or this is the last file
    if (size(collectedFrames, 1) >= 2e5 || fileIndex == numel(dcdFilenames))
        1;
        
        %%
        analysis = AnalysisPlan();
        analysis.UseParFor = shouldProcessParallel;
        analysis.UseProgressMonitor = false;
        
        %analysis.BranchParameter('WindowLength', windowLengthList);
        %analysis.AddAdvancedConversion(@(obj, i)collectedFrames(s - 1 + [1:obj.Parameters.WindowLength], :, :));
        analysis.AddConversion(@(s)collectedFrames(s - 1 + [1:windowLength], :, :));
        
        % Order coordinates by XYXYXY...
        analysis.AddConversion(@(c)permute(c, [1 3 2]));
        
        % Linearize configuration for the following coarse graining
        analysis.AddConversion(@LinearizeConfigurations);
        
        %
        analysis.Parameters.CoarseGrain_NumOfStates = 24;
        %analysis.Parameters.CoarseGrain_kStd = 4.5;
        analysis.Parameters.TightCoordinatePacking = false;
        
        numOfStatesList = [20:25];
        %analysis.BranchParameter('CoarseGrain_NumOfStates', numOfStatesList);
        
        
        analysis.Parameters.CoarseGrain_EmbeddingAlphabetFactor = 1;
        analysis.Parameters.CoarseGrain_EmbeddingAlphabetMode = 5;
        
        %embeddingAlphabetFactorList = [0 2.^[0:12]];
        embeddingAlphabetFactorList = [0 1];
        %analysis.BranchParameter('CoarseGrain_EmbeddingAlphabetFactor', embeddingAlphabetFactorList);
        %analysis.Parameters.CoarseGrain_NumOfBits = 16;
        %analysis.BranchParameter('CoarseGrain_NumOfBits', [8:16]);
        
        %analysis.BranchParameter('CoarseGrain_EmbeddingAlphabetMode', [1 3]);
        
        analysis.Parameters.CoarseGrain_RangeMin = (-pi - VerySmallNumber) * ones(1, 66);
        analysis.Parameters.CoarseGrain_RangeMax = ( pi + VerySmallNumber) * ones(1, 66);
        
        analysis.Parameters.CompressionEngine = 1; % LZMA
        %analysis.Parameters.CompressionEngine = 2; % my 16 bit level compression, ignores the number of used bits
        %analysis.BranchParameter('CompressionEngine', [1 4 3]);
        
        analysis.AddAnalysis(@EstimateEntropyForCoarseGrainedConfigurations);
        %analysis.AddAnalysis(@EstimateEntropyOfUncorrelatedVariables);
        
        %analysis.Parameters.OutputStyle = 2;
        %analysis.AddAnalysis(@EstimateEntropyForCoarseGrainedConfigurationsOld);
        
        windowStartList = 1 + [0:windowStride:(size(collectedFrames, 1) - windowLength)];
        results = analysis.Analyze(arrayfun(@(x)x, windowStartList, 'UniformOutput', false));
        collectedResults{end + 1} = results;
        
        % Keep last frames to continue sliding the window
        %collectedFrames = collectedFrames((windowStride+1):windowLength, :, :);
        numOfFramesToRemove = FloorTo(size(collectedFrames, 1) - windowLength + windowStride, windowStride);
        
        collectedFrames = collectedFrames(numOfFramesToRemove+1:end, :, :);
        1;
    end
    
    1;
end



%%


timeStampString = datestr(now, 'yyyymmdd-hhMMss');

clear frames collectedFrames analysis;
save([timeStampString sprintf(' - Villin_results %s.mat', simulationPrefix)]);
%save([timeStampString sprintf(' - Villin_results.mat')], 'numOfStatesList', 'windowStartList', 'results', 'windowLength', 'windowStride', 'states', 'simulationPrefix', 'outputRoot', 'dataRoot', 'dataFolder');

%
%%
x = arrayfun(@(results)results.EntropyEstimate, results);
%x = min(x, [], 2);

figure(1);
%plot(numOfStatesList, x, 'LineWidth', 2);
bar(windowStartList, x - mean(x));



a = gca();
a.FontSize = 14;


xlabel('Time', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('$S_{A}$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');


%return;
%%
results = vertcat(collectedResults{:});
x = arrayfun(@(results)results.EntropyEstimate, results);
x = min(x, [], 2);

movingAvgKernel = ones(1, windowLength) * (1/windowLength);
slidingFrameMeanTime = conv(collectedFramesTimes, movingAvgKernel, 'valid');
slidingFrameMeanTime = slidingFrameMeanTime(1:windowStride:end);

sortedValues = sort(x);
x05 = sortedValues(round(end*0.05));
x95 = sortedValues(round(end*0.95));

figure(1);
%plot(numOfStatesList, x, 'LineWidth', 2);
%bar(slidingFrameMeanTime, x - mean(x));
%bar(slidingFrameMeanTime, x - 0.5 * [x05 + x95]);
bar(slidingFrameMeanTime, x - 47.8);

stateTimes = collectedFramesTimes(1) + states(:, 1) * dt;
hold on;
%plot(stateTimes, states(:, 2) * 50 - 15);
plot(stateTimes, states(:, 2) * 60 - 30);
hold off;


a = gca();
a.FontSize = 14;


xlabel('Time', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('$S_{A}$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');

%xlim(1.0e+07 * [         0    3.2107]);

%% Prepare data for origin plots (timeline)

results = vertcat(collectedResults{:});
entropy = arrayfun(@(results)results.EntropyEstimate, results);
entropy = min(entropy, [], 2);

sortedValues = sort(entropy);
x05 = sortedValues(round(end*0.05));
x95 = sortedValues(round(end*0.95));
shiftValue = 0.5 * [x05 + x95];

movingAvgKernel = ones(1, windowLength) * (1/windowLength);
slidingFrameMeanTime = conv(collectedFramesTimes, movingAvgKernel, 'valid');
slidingFrameMeanTime = slidingFrameMeanTime(1:windowStride:end);

% Resample avg state
stateTimes = collectedFramesTimes(1) + states(:, 1) * dt;
%slidingStatesMeanTime = conv(stateTimes', movingAvgKernel, 'valid');
%slidingAvgState = conv(states(:, 2)', movingAvgKernel, 'valid');
resampledStates = interp1(stateTimes, states(:, 2), slidingFrameMeanTime);
%resampledStates = interp1(slidingStatesMeanTime, slidingAvgState, slidingFrameMeanTime);
%resampledStates = round(resampledStates);

m = [slidingFrameMeanTime(:) entropy-shiftValue 1-resampledStates(:) resampledStates(:)];
open('m');


%% Prepare data for origin plots (scatter)

% Columns are: Time , Total energy , Potential energy , Kinetic energy , Thermostat
energies = LoadFilePattern([dataFolderUpper 'Energies_*.dat'], '-ascii');


%%
m = regexp(simulationPrefix, '-(?<T>\d+)K-', 'names');
T = str2double(m.T);

kcalpermolInJoule = 6.9477e-21; % J
boltzmannConst = 1.38064852e-23; % J/K

kT_in_J = boltzmannConst * T; % Notice the temperature
kT_in_kcalpermol = kT_in_J/kcalpermolInJoule;

% Columns are: Time , Total energy , Potential energy , Kinetic energy , Thermostat
enthalpy_kcal_per_mol = energies(:, 3); % potential energy
sortedValues = sort(enthalpy_kcal_per_mol);
x05 = sortedValues(round(end*0.05));
x95 = sortedValues(round(end*0.95));
enthalpyShiftValue = 0.5 * [x05 + x95] + 5; % shift for 360K
%enthalpyShiftValue = 0.5 * [x05 + x95]; % shift for 370K
%enthalpyShiftValue = -3.9724e+04; % "0.5 * [x05 + x95]" for 360K
%enthalpyShiftValue = 0.5 * [x05 + x95] - 5; % shift for 380K

enthalpy_kcal_per_mol = enthalpy_kcal_per_mol - enthalpyShiftValue;

movingAvgKernel = ones(1, windowLength) * (1/windowLength);

slidingFrameMeanTime = conv(collectedFramesTimes, movingAvgKernel, 'valid');
slidingFrameMeanTime = slidingFrameMeanTime(1:windowStride:end);

slidingEnergyMeanTime = conv(energies(:, 1)', movingAvgKernel, 'valid');
slidingAvgEnergy = conv(enthalpy_kcal_per_mol', movingAvgKernel, 'valid');

resampledEnergy = interp1(slidingEnergyMeanTime, slidingAvgEnergy, slidingFrameMeanTime);

results = vertcat(collectedResults{:});

entropy = arrayfun(@(results)results.EntropyEstimate, results);
entropy = min(entropy, [], 2);

sortedValues = sort(entropy);
x05 = sortedValues(round(end*0.05));
x95 = sortedValues(round(end*0.95));
entropyShiftValue = 0.5 * [x05 + x95];
entropy = entropy - entropyShiftValue;

stateTimes = collectedFramesTimes(1) + states(:, 1) * dt;
resampledStates = interp1(stateTimes, states(:, 2), slidingFrameMeanTime);
resampledStates = (resampledStates == 0);

rng(11);

figure(2);
%histogram(enthalpy_kcal_per_mol);
%bar(slidingEnergyMeanTime, slidingAvgEnergy);
%plot(resampledEnergy, entropy, '-k');
plot(resampledEnergy(resampledStates), entropy(resampledStates), '.r');
hold on;
plot(resampledEnergy(~resampledStates), entropy(~resampledStates), '.b');

if (1)
    t = slidingFrameMeanTime(resampledStates);
    x = resampledEnergy(resampledStates);
    y = entropy(resampledStates);
    
    xy0 = [-16.1597  -23.9552];
    wh = [9.0750    3.6695];
    xy0 = xy0 + 0.5 * wh;
    
    which = ((x(:)-xy0(1))/wh(1)).^2 + ((y(:)-xy0(2))/wh(2)).^2 <= 0.25;
    which = find(which);
    
    [~, order] = sort(rand(size(which)));
    %which = which(order(1:20));
    
    t1 = t(which);
    plot(x(which), y(which), 'ok');
    
    
    xy0 = [-16.2202  -16.5398];
    wh = [9.0750    3.6695];
    xy0 = xy0 + 0.5 * wh;
    
    which = ((x(:)-xy0(1))/wh(1)).^2 + ((y(:)-xy0(2))/wh(2)).^2 <= 0.25;
    which = find(which);
    
    [~, order] = sort(rand(size(which)));
    %which = which(order(1:20));
    
    t2 = t(which);
    plot(x(which), y(which), 'ok');
    
     %save('Villin 360K folded times at two entropy clusters.mat', 't1', 't2');
end

if (0)
    
    t = slidingFrameMeanTime(resampledStates);
    x = resampledEnergy(resampledStates);
    y = entropy(resampledStates);
    
    xy0 = [-15.1315  -26.1722];
    wh = [0.3048   19.0356];
    
    %
    xy = bsxfun(@times, [(x(:)-xy0(1)) (y(:)-xy0(2))], 1 ./ wh);
    %figure(10);
    %plot(xy(:, 1), xy(:, 2));
    which = (xy >= 0) & (xy <= 1);
    which = which(:, 1) & which(:, 2);
    which = find(which);
    
    % shuffle
    [~, order] = sort(rand(size(which)));
    which = which(order);
    
    [N,edges,bin] = histcounts(y(which), 20);
    % Choose one of each bin
    [C,ia,ic] = unique(bin);
    
    which = which(ia);
    
    % Order according to entropy
    [~, order] = sort(y(which));
    which = which(order);
    
    t1 = t(which);
    
    plot(x(which), y(which), 'ok');
    
    
    %figure(10);
    %centers = 0.5 * [edges(2:end) + edges(1:end-1)];
    %bar(centers, N);
    
    %figure(11);
    %histogram(y(which), 100)
    
    
    %save('Villin 360K times with constant energy -15.mat', 't1');
    
    
end


if (0)
    %%
    
    t = slidingFrameMeanTime(resampledStates);
    x = resampledEnergy(resampledStates);
    y = entropy(resampledStates);
    
    xy0 = [-15.1315  -26.1722];
    wh = [0.3048   19.0356];
    
    %
    xy = bsxfun(@times, [(x(:)-xy0(1)) (y(:)-xy0(2))], 1 ./ wh);
    %figure(10);
    %plot(xy(:, 1), xy(:, 2));
    which = (xy >= 0) & (xy <= 1);
    which = which(:, 1) & which(:, 2);
    which = find(which);
    
    % shuffle
    [~, order] = sort(rand(size(which)));
    which = which(order);
    
    [N,edges,bin] = histcounts(y(which), 20);
    % Choose one of each bin
    [C,ia,ic] = unique(bin);
    
    which = which(ia);
    
    % Order according to entropy
    [~, order] = sort(y(which));
    which = which(order);
    
    t1 = t(which);
    
    plot(x(which), y(which), '*k');
    
    
    %figure(10);
    %centers = 0.5 * [edges(2:end) + edges(1:end-1)];
    %bar(centers, N);
    
    %figure(11);
    %histogram(y(which), 100)
    
    
    save('Villin 360K times with constant energy -15.mat', 't1');
    
    
end

hold off;

m1 = [resampledEnergy(resampledStates)', entropy(resampledStates)];
open('m1');

m2 = [resampledEnergy(~resampledStates)', entropy(~resampledStates)];
open('m2');


%%
figure(4);
histogram(entropy, 100);

figure(5);
histogram(resampledEnergy, 100);

%%

m = [slidingFrameMeanTime(:) entropy-shiftValue 1-resampledStates(:) resampledStates(:)];
open('m');



%% Prepare data for origin plots (scatter histograms)

Xedges = linspace(-30, 35, 66);
Yedges = linspace(-30, 30, 61);
Xedges = linspace(-30, 35, 45);
Yedges = linspace(-30, 30, 45);

histogram_x = 0.5 * (Xedges(1:end-1) + Xedges(2:end));
histogram_y = 0.5 * (Yedges(1:end-1) + Yedges(2:end));

[h] = histcounts2(resampledEnergy(:), entropy, Xedges, Yedges);
h = h';

figure(4);
clf;
%surface(h);
imagesc(histogram_x, histogram_y, h)
axis xy;

xlabel('Energy (kcal/mol)');
ylabel('S / KB');

m = zeros(size(h) + 1);
m(2:end, 2:end) = h ./ sum(h(:));
m(1, 2:end) = histogram_x;
m(2:end, 1) = histogram_y;

