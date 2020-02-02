
% Protein folding kinetics and thermodynamics from atomistic simulation
% (2012)
% doi: 10.1073/pnas.1201811109
% http://www.pnas.org/content/109/44/17845.full

addpath('c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\Analysis Scripts');
addpath('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\Common Code');

% Load .Net engines
NET.addAssembly('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\PhysiCasino\PhysiCasino\bin\x64\Release\PhysiCasino.exe');
NET.addAssembly('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\PhysiCasino\DcdSharp\bin\x64\Release\DcdSharp.exe');

outputRoot = '../PhysiCasino Output/';
timeStampString = datestr(now, 'yyyymmdd-hhMMss');

%%
VerySmallNumber = 1e-10;

%%
pdb = pdbread('a3d.pdb');
atomNames = {pdb.Model.Atom.AtomName};
atomResidueId = [pdb.Model.Atom.resSeq];

%%
for n = 1:27
%n = 2;
rama = load(['a3d.' num2str(n) '.rama.pkl.mat']);
pe = load(['a3d.' num2str(n) '.pe.pkl.mat']);
pe = pe.data';
frames = double(rama.data);

% throw out the first 2000 frames
frames = frames(2000:end, :, :);
pe = pe(2000:end);

clear rama;


%%
numOfStatesOptimizationMode = false;

shouldProcessParallel = true;

if (shouldProcessParallel)
    % pctRunOnAll
    PrepWorkersPool();
end


%%
windowLength = 1000;
windowStride = 50;%0.1 * windowLength;

analysis = AnalysisPlan();
analysis.UseParFor = shouldProcessParallel;

analysis.Parameters.CoarseGrain_NumOfStates = 14;
%analysis.Parameters.CoarseGrain_kStd = 4.5;

if (numOfStatesOptimizationMode)
    %numOfStatesList = [2:4:250];
    numOfStatesList = [2:60];
    analysis.BranchParameter('CoarseGrain_NumOfStates', numOfStatesList);
end

analysis.Parameters.CoarseGrain_RangeMin = (-pi - VerySmallNumber) * ones([1 size(frames, 2)*size(frames, 3)]);
analysis.Parameters.CoarseGrain_RangeMax = ( pi + VerySmallNumber) * ones([1 size(frames, 2)*size(frames, 3)]);

if (~numOfStatesOptimizationMode)
    analysis.AddConversion(@(s)frames(s - 1 + [1:windowLength], :, :));
end

if (1)
    % Order coordinates by XYXYXY...
    analysis.AddConversion(@(c)permute(c, [1 3 2]));
end

analysis.AddConversion(@LinearizeConfigurations);

analysis.AddAnalysis(@EstimateEntropyForCoarseGrainedConfigurations);

if (numOfStatesOptimizationMode)
    results = analysis.Analyze(frames(1e3 + [1:windowLength], :, :));
else
    windowStartList = 1 + [0:windowStride:(size(frames, 1) - windowLength)];
    results = analysis.Analyze(arrayfun(@(x)x, windowStartList, 'UniformOutput', false));
end

%clear frames collectedFrames analysis;
if (numOfStatesOptimizationMode)
    save([timeStampString ' states optimization.mat'], 'results');
else
    save([timeStampString sprintf(' a3d_%d sliding window (%d frames, %d stride).mat', n, windowLength, windowStride)], 'results');
end

if (numOfStatesOptimizationMode)
    %%
    x = arrayfun(@(results)results.EntropyEstimate, results);
    %x = min(x, [], 2);
    
    figure(1);
    plot(numOfStatesList, x, 'LineWidth', 2);
    
    
    a = gca();
    a.FontSize = 14;
    
    
    xlabel('num of states', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('$S_{A}$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
else
    %%
    x = arrayfun(@(results)results.EntropyEstimate, results);
    
    %x(isinf(x)) = 1e5;
    %x(isnan(x)) = -1;
    
    figure(1);
    bar(windowStartList, x - mean(x));
    
    
    a = gca();
    a.FontSize = 14;
    
    
    xlabel('Window start', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
    ylabel('$S_{A}$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
    
end

if (~numOfStatesOptimizationMode)
    %%
    %w = windowStartList > 100;
    w = windowStartList > 0;
    
    energy = pe;
    
    framesIndexes = 1:size(frames, 1);
    movingAvgKernel = ones(1, windowLength) * (1/windowLength);
    slidingFrameMeanIndex = conv(framesIndexes, movingAvgKernel, 'valid');
    slidingFrameMeanIndex = slidingFrameMeanIndex(1:windowStride:end);
    
    slidingEnergyMean = conv(energy, movingAvgKernel', 'valid');
    slidingEnergyMean = slidingEnergyMean(1:windowStride:end);
    
    x = arrayfun(@(results)results.EntropyEstimate, results);
    
    figure(2);
    %plot(slidingEnergyMean(w), x(w) - mean(x(w)), 'x');
    plot3(slidingEnergyMean(w), x(w) - mean(x(w)), windowStartList(w), '.');
    xlabel('H', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('S', 'FontSize', 16, 'FontWeight', 'bold');
    
    view(0, 90);
    
    dcm = datacursormode(gcf());
    dcm.UpdateFcn = @MyDataCursorText;

end

title(sprintf('n = %d', n));
pause(1);

end % for n = ...

%%
return;



%
%%
x = arrayfun(@(results)results.EntropyEstimate, results);
%x = min(x, [], 2);

figure(1);
plot(numOfStatesList, x, 'LineWidth', 2);


a = gca();
a.FontSize = 14;


xlabel('num of states', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
ylabel('$S_{A}$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');


%%
x = arrayfun(@(results)results.EntropyEstimate, results);

w = windowStartList > 100;

figure(1);
bar(windowStartList(w), x(w) - mean(x(w)));

%%
%w = windowStartList > 100;
w = windowStartList > 0;

energy = pe;

framesIndexes = 1:size(frames, 1);
movingAvgKernel = ones(1, windowLength) * (1/windowLength);
slidingFrameMeanIndex = conv(framesIndexes, movingAvgKernel, 'valid');
slidingFrameMeanIndex = slidingFrameMeanIndex(1:windowStride:end);

slidingEnergyMean = conv(energy, movingAvgKernel', 'valid');
slidingEnergyMean = slidingEnergyMean(1:windowStride:end);

x = arrayfun(@(results)results.EntropyEstimate, results);

figure(2);
plot(slidingEnergyMean(w), x(w) - mean(x(w)), 'x');
xlabel('H', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('S', 'FontSize', 16, 'FontWeight', 'bold');

dcmObj = datacursormode(gcf());
%dcmObj.UpdateFcn = @(~,event_obj)num2str(slidingRmsdMean(event_obj.DataIndex));
dcmObj.UpdateFcn = @(~,event_obj)PlotWindow(event_obj, frames(windowStride*event_obj.DataIndex - 1 + [1:windowLength], :, :), slidingRmsdMean);

%%
framesIndexes = 1:size(frames, 1);
movingAvgKernel = ones(1, windowLength) * (1/windowLength);
slidingFrameMeanIndex = conv(framesIndexes, movingAvgKernel, 'valid');
slidingFrameMeanIndex = slidingFrameMeanIndex(1:windowStride:end);

slidingRmsdMean = conv(rmsd, movingAvgKernel', 'valid');
slidingRmsdMean = slidingRmsdMean(1:windowStride:end);

x = arrayfun(@(results)results.EntropyEstimate, results);

figure(2);
plot(slidingRmsdMean(w), x(w) - mean(x(w)), 'x');
%plot(slidingRmsdMean, slidingEnergyMean - 0.5*(x - mean(x)), 'x');

xlabel('rmsd', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('S', 'FontSize', 16, 'FontWeight', 'bold');

figure(3);
plot(slidingRmsdMean(w), slidingEnergyMean(w), 'x');

xlabel('rmsd', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('H', 'FontSize', 16, 'FontWeight', 'bold');

%%
figure(3);
histogram(slidingRmsdMean(w), 50)
xlabel('rmsd', 'FontSize', 16, 'FontWeight', 'bold');

%%
figure(3);
histogram(slidingEnergyMean(w), 50)
xlabel('H', 'FontSize', 16, 'FontWeight', 'bold');

%%
figure(3);
histogram(x(w) - mean(x(w)), 50)
xlabel('S', 'FontSize', 16, 'FontWeight', 'bold');

%%
x = arrayfun(@(results)results.EntropyEstimate, results);

figure(2);
plot3(slidingRmsdMean, slidingEnergyMean, x - mean(x), 'o', 'LineWidth', 2, 'MarkerSize', 6);
xlabel('rmsd');
ylabel('H');
zlabel('S');

dcmObj = datacursormode(gcf());
%dcmObj.UpdateFcn = @(~,event_obj)num2str(event_obj.DataIndex);
dcmObj.UpdateFcn = @(~,event_obj)PlotWindow(event_obj, frames0(windowStride*event_obj.DataIndex - 1 + [1:windowLength], :, :));


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

