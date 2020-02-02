
k = 1;
switch (k)
    case 1
        s100 = load('20180527-164110 hyp4 sliding window (100 frames, 50 stride).mat');
        s500 = load('20180527-164439 hyp4 sliding window (500 frames, 50 stride).mat');
        s1000 = load('20180527-164538 hyp4 sliding window (1000 frames, 50 stride).mat');
        s2000 = load('20180527-164725 hyp4 sliding window (2000 frames, 50 stride).mat');
        n = 4;
        windowLengthList = [100 500 1000 2000];
        pe = load('hyp.4.pe.pkl.mat'); pe = pe.data;
        
    case 2
        s100 = load('20180527-165608 hyp3 sliding window (100 frames, 50 stride).mat');
        s500 = load('20180527-165649 hyp3 sliding window (500 frames, 50 stride).mat');
        s1000 = load('20180527-165747 hyp3 sliding window (1000 frames, 50 stride).mat');
        s2000 = load('20180527-165916 hyp3 sliding window (2000 frames, 50 stride).mat');
        n = 3;
        windowLengthList = [100 500 1000 2000];
        pe = load('hyp.3.pe.pkl.mat'); pe = pe.data;
        
    case 3
        s = [];
        s = [s load('20180527-165608 hyp3 sliding window (100 frames, 50 stride).mat')];
        s = [s load('20180528-140322 hyp3 sliding window (200 frames, 50 stride).mat')];
        s = [s load('20180528-140353 hyp3 sliding window (300 frames, 50 stride).mat')];
        s = [s load('20180528-140542 hyp3 sliding window (400 frames, 50 stride).mat')];
        s = [s load('20180527-165649 hyp3 sliding window (500 frames, 50 stride).mat')];
        s = [s load('20180528-141222 hyp3 sliding window (600 frames, 50 stride).mat')];
        n = 3;
        windowLengthList = [100:100:600];
        pe = load('hyp.3.pe.pkl.mat'); pe = pe.data;
end


pe = load(['hyp.' num2str(n) '.pe.pkl.mat']);
pe = pe.data';

figure(2);
clf;
cla;

for i = 1:numel(s)
    results = s(i).results;
    
    windowLength = windowLengthList(i);
    windowStride = 50;
    windowStartList = ([1:numel(results)]-1)*windowStride + 1;
    
    w = windowStartList > 0;
    
    energy = pe;
    
    movingAvgKernel = ones(1, windowLength) * (1/windowLength);
    
    slidingEnergyMean = conv(energy, movingAvgKernel', 'valid');
    slidingEnergyMean = slidingEnergyMean(1:windowStride:end);
    
    x = arrayfun(@(results)results.EntropyEstimate, results);
    
    if (i ~= 1)
        hold on;
    end
    
    plot(slidingEnergyMean(w), x(w) - mean(x(w)), 'x');
    %plot3(slidingEnergyMean(w), x(w) - mean(x(w)), windowStartList, '.');
    
    if (i ~= 1)
        hold off;
    end
end


xlabel('H', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('S', 'FontSize', 16, 'FontWeight', 'bold');
legend(arrayfun(@(n)sprintf('%d frames', n), windowLengthList, 'UniformOutput', false), 'FontSize', 16, 'FontWeight', 'bold');

view(0, 90);
    