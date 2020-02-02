function output_txt = PolyChargedDataCursorFunc(obj, event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).
% Get the position of the point being clicked on 
pos = get(event_obj,'Position');
I = get(event_obj, 'DataIndex');

filepath = event_obj.Target.UserData;

%%
entropy_data = load(filepath);

[folder, filename] = fileparts(filepath);

m = regexp(filename, 'conf\w* (\d+)', 'tokens');
m = m{1};
if (iscell(m))
    m = m{1};
end
conf = str2double(m);

m = regexp(filename, '(?:Entropy for )?(.+\.angles)', 'tokens');
m = m{1};
if (iscell(m))
    m = m{1};
end

frames = load([folder filesep m '.mat']);
frames = frames.frames;

frames(:, 1) = 0;
frames = cat(3, cos(frames), sin(frames));
frames = cat(2, frames(:, 1, :) .* 0, frames);
frames = cumsum(frames, 2);

start_frame = (I - 1) * entropy_data.parameters.stride;
frame_indexes = start_frame + [1:entropy_data.parameters.window_length];

%frame_indexes = frame_indexes(1:10:end); % dilute

frames = frames(frame_indexes, :, :);

%%
%ij = AllPairs(numel(frame_indexes), 1);

frames1 = frames;
frames1(:, :, 3) = frames1(:, :, 2) .* 0;
%[all_pairs_rmsd, ij] = AllPairsRmsd(permute(frames1, [2 3 1]));
[rmsd, ~, T, adjusted_frames] = CalculateAllRmsd(frames1(1, :, :), frames1, 1, 2, 3);

frames = permute(adjusted_frames(:, 1:2, :), [3 1 2]);

%%
figure(10);
%p = plot(frames(:, :, 1), frames(:, :, 2), '-', 'Color', [0 0 0 0.05]);
plot(frames(1, :, 1), frames(1, :, 2), ':g', 'LineWidth', 3);
hold on;
%plot(frames(frame_indexes, [1 20 40], 1), frames(frame_indexes, [1 20 40], 2), '+r', 'Color', [1 0 0]);
%plot(frames(frame_indexes, [60 80 99], 1), frames(frame_indexes, [60 80 99], 2), 'ob', 'Color', [0 0 1]);

switch (conf)
    case 6
        % conf 6
        p = [1 100];
        n = [49 54];
        
    case 7
        % conf 7
        p = [1 40 100];
        n = [20 60 80];
        
    case 8
        % conf 8
        p = [1 5 10];
        n = [90 95 100];
        
    case 9
        % conf 7
        p = [5 40 100];
        n = [20 60 80];
        
    case 10
        % conf 8
        p = [5 10 15];
        n = [90 95 100];
end

if (0)
plot(frames(1, :, 1), frames(1, :, 2), ':g', 'LineWidth', 3);
plot(frames(end, :, 1), frames(end, :, 2), ':k', 'LineWidth', 3);
end

plot(frames(1, p, 1), frames(1, p, 2), 'og', 'LineWidth', 4);
plot(frames(1, n, 1), frames(1, n, 2), 'og', 'LineWidth', 4);
plot(frames(1, p, 1), frames(1, p, 2), '+r', 'LineWidth', 4);
plot(frames(1, n, 1), frames(1, n, 2), 'xb', 'LineWidth', 4);

if (0)
plot(frames(end, p, 1), frames(end, p, 2), 'ok', 'LineWidth', 4);
plot(frames(end, n, 1), frames(end, n, 2), 'ok', 'LineWidth', 4);
plot(frames(end, p, 1), frames(end, p, 2), '+r', 'LineWidth', 4);
plot(frames(end, n, 1), frames(end, n, 2), 'xb', 'LineWidth', 4);
end
hold off;

SetMyDefaultFigureSettings();
ax = gca();
ax.Box = 'off'
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

set(gcf, 'Position', [ 530   225   800   600]);

output_txt = sprintf('H = %0.1f, S = %0.1f', pos(1), pos(2));
