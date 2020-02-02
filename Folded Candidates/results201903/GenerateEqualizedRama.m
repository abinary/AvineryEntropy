sim_name = 'nug2';

rama_filepath = @(i)['c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Sosnick\Data\20190308\' ...
    sprintf('%s\\%s.%d.rama.pkl.mat', sim_name, sim_name, i)];

out_filepath = @(i)['c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Sosnick\Data\20190308\' ...
    sprintf('%s\\%s.%d.eql.mat', sim_name, sim_name, i)];

all_frames = [];
for i = 0:27
    rama_filepath(i)
    load(rama_filepath(i));    
    all_frames = [all_frames; data];
end

all_frames = permute(all_frames, [2 3 1]);

%%
x = all_frames(:, 1, :);
y = all_frames(:, 2, :);

x = x(:);
y = y(:);

if (1)
    %%
    dihedral_histograms = {};
    
    for k = 1:256
        k
        k1 = 1 + mod(k - 1, 16);
        k2 = ceil(k / 16);

        edges1 = linspace(-pi, pi, k1 + 1);
        edges2 = linspace(-pi, pi, k2 + 1);
        
        edges1(1) = -4;
        edges1(end) = 4;
        edges2(1) = -4;
        edges2(end) = 4;
        
        x = all_frames(:, 1, :);
        y = all_frames(:, 2, :);
        
        x = x(:);
        y = y(:);
        
        x = discretize(x, edges1);
        y = discretize(y, edges2);
        
        x = int32(round(x));
        y = int32(round(y));

        n = accumarray([x y], 1, [k1 k2], @sum);
        dihedral_histograms{k1, k2} = n;
    end
end

%kmedoids([x y], 10, 'Distance', @(p1, p2)norm([ShortestAngleDiff(p1(1), p2(1)) ShortestAngleDiff(p1(2), p2(2))]));

%figure(1);

k = 2^7;
edges = linspace(-pi, pi, k + 1);
centers = 0.5 * (edges(1:end-1) + edges(2:end));
edges(1) = -4;
edges(end) = 4;

x = discretize(x, edges);
y = discretize(y, edges);

x = int32(round(x));
y = int32(round(y));


n = accumarray([x y], 1, [k k], @sum);

%%
figure(1);
clf;
imagesc(centers, centers, n');

SetMyDefaultFigureSettings();

fig_scale_factor = 1.5;

f = gcf();
ax = gca();
ax.Units = 'centimeters';
ax.Position = [1.7 1.7 2*[4.3 4.3]] * fig_scale_factor;

f.Color = [1 1 1];
f.Units = 'centimeters';
f.Position = [f.Position(1:2) ax.Position(3:4) + [2.5 2.5] * fig_scale_factor];

xlabel('\phi (rad)');
ylabel('\psi (rad)');

ax.XTick = ax.YTick;

hold on;
for i = 1:11
    plot([-pi pi], [1 1] * 2*pi*i/12 - pi, '--g', 'LineWidth', 2);
    plot([1 1] * 2*pi*i/12 - pi, [-pi pi], '--g', 'LineWidth', 2);
end
hold off;

axis equal xy;
%SetPlotRectRatio(1.0);
%SaveFig('Rama Histogram.jpg');

%% 1st boxing algorithm
boxes = [1 1 k k];

BoxSum = @(box)sum(sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1)));

max_aspect_ratio = 4;
cond_dir = 0;
    
for split_idx = 1:50
    
    % Choose the box containing the max intensity
    b = arrayfun(@(b)BoxSum(boxes(b, :)), 1:size(boxes, 1));
    %b = arrayfun(@(b)max(max(n(boxes(b, 1):boxes(b, 3)-1, boxes(b, 2):boxes(b, 4)-1))), 1:size(boxes, 1));
    
    b((boxes(:, 3) == 1) & (boxes(:, 4) == 1)) = inf; % ignore boxes which cannot be split
    [~, b] = max(b);
    
    % split either horizontally or vertically, whichever is more equalized
    
    box = boxes(b, :);
    dv = [];
    dh = [];
    clear bv1 bv2 bh1 bh2;
    
    if (box(4) >= 2 && box(3)/box(4) <= max_aspect_ratio)
        bv1 = [box(1:3) box(4) / 2];
        bv2 = [box(1) box(2)+box(4)/2 box(3) box(4)/2];
        i1 = BoxSum(bv1);
        i2 = BoxSum(bv2);
        dv = abs(i1 - i2);
    else
        1;
    end
    
    if (box(3) >= 2 && box(4)/box(3) <= max_aspect_ratio)
        bh1 = [box(1:2) box(3)/2 box(4)];
        bh2 = [box(1)+box(3)/2 box(2) box(3)/2 box(4)];
        i1 = BoxSum(bh1);
        i2 = BoxSum(bh2);
        dh = abs(i1 - i2);
    else
        1;
    end
    
    % Choose the split resulting in more equal boxes
    if (isempty(dv) || (~isempty(dh) && (cond_dir * dh < cond_dir * dv)))
        boxes = [boxes(1:b-1, :); boxes(b+1:end, :); bh1; bh2];
    else
        boxes = [boxes(1:b-1, :); boxes(b+1:end, :); bv1; bv2];
    end
    
    1;
end

figure(1)
imagesc(centers, centers, n);
axis xy;
%imagesc(n);
%axis ij;

hold on;
for b = 1:size(boxes, 1)
    
    plot(...
        centers([boxes(b, 2) boxes(b, 2)+boxes(b, 4)-1 boxes(b, 2)+boxes(b, 4)-1 boxes(b, 2) boxes(b, 2)]), ...
        centers([boxes(b, 1) boxes(b, 1) boxes(b, 1)+boxes(b, 3)-1 boxes(b, 1)+boxes(b, 3)-1 boxes(b, 1)]), '--', 'LineWidth', 2);
end
hold off;

%% 2nd boxing algorithm
boxes = [1 1 k k];

BoxSum = @(box)sum(sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1)));

%BoxVSum = @(box)sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1), 1);
BoxVCumSum = @(box)cumsum(sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1), 2), 1);
%BoxHSum = @(box)sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1), 2);
BoxHCumSum = @(box)cumsum(sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1), 1), 2);

max_aspect_ratio = 4;
cond_dir = 0;
    
clc;
for split_idx = 1:15
    % Choose the box containing the max intensity
    b = arrayfun(@(b)BoxSum(boxes(b, :)), 1:size(boxes, 1));
    %b = arrayfun(@(b)max(max(n(boxes(b, 1):boxes(b, 3)-1, boxes(b, 2):boxes(b, 4)-1))), 1:size(boxes, 1));
    
    b((boxes(:, 3) == 1) & (boxes(:, 4) == 1)) = 0; % ignore boxes which cannot be split
    [~, b] = max(b);
    
    % split either horizontally or vertically, whichever is more equalized
    
    box = boxes(b, :);
    
    % decide on split direction
    if (box(3)/box(4) > 1)
        if (box(3) == 2)
            idx = 1;
        else
            bs = BoxSum(box);
            bcs = BoxVCumSum(box);
            
            idx = round(interp1(bcs ./ bs, 0.5:box(3), 0.5));
            
            if (idx < 2)
                idx = 2;
            elseif (idx > box(3)-1)
                idx = box(3)-1;
            end
        end
        
        b1 = [box(1:2) idx box(4)];
        b2 = [box(1)+idx box(2) box(3)-idx box(4)];
    else
        if (box(4) == 2)
            idx = 1;
        else
            bs = BoxSum(box);
            bcs = BoxHCumSum(box);
            idx = round(interp1(bcs ./ bs, 0.5:box(4), 0.5));
            
            if (idx < 2)
                idx = 2;
            elseif (idx > box(4)-1)
                idx = box(4)-1;
            end
        end
        
        b1 = [box(1:3) idx];
        b2 = [box(1) box(2)+idx box(3) box(4)-idx];
    end
    
    boxes = [boxes(1:b-1, :); boxes(b+1:end, :); b1; b2];
end

figure(1)
imagesc(centers, centers, n);
axis xy;
%imagesc(n);
%axis ij;

hold on;
for b = 1:size(boxes, 1)
    if (boxes(b, 3) == 1 && boxes(b, 4) == 1)
        plot(centers(boxes(b, 2)), centers(boxes(b, 1)), '*');
    else
        plot(...
            centers([boxes(b, 2) boxes(b, 2)+boxes(b, 4)-1 boxes(b, 2)+boxes(b, 4)-1 boxes(b, 2) boxes(b, 2)]), ...
            centers([boxes(b, 1) boxes(b, 1) boxes(b, 1)+boxes(b, 3)-1 boxes(b, 1)+boxes(b, 3)-1 boxes(b, 1)]), '--', 'LineWidth', 2);
    end
end
hold off;

% assignment matrix
Z = zeros(size(n));
for b = 1:size(boxes, 1)
    Z(boxes(b, 1)+[1:boxes(b, 3)]-1, boxes(b, 2)+[1:boxes(b, 4)]-1) = b;
end

figure(2);
imagesc(Z);
axis xy;

if (0)
    %% sanity check
    Z = zeros(size(n));
    for b = 1:size(boxes, 1)
        Z(boxes(b, 1)+[1:boxes(b, 3)]-1, boxes(b, 2)+[1:boxes(b, 4)]-1) = ...
            Z(boxes(b, 1)+[1:boxes(b, 3)]-1, boxes(b, 2)+[1:boxes(b, 4)]-1) + 1;
    end
    
    figure(2);
    imagesc(Z);
    axis xy;
end

%%
p = n / sum(n(:));
k = 15;
[row, col] = find(p > 0.0001);
w = n(sub2ind(size(n), row, col));
[X, Y] = meshgrid(centers, centers);
XY = [X(sub2ind(size(X), row, col)) Y(sub2ind(size(X), row, col))];

[i, j] = AllPairs(size(XY, 1), 1);
D = sqrt(sum([ShortestAngleDiff(XY(i, 1), XY(j, 1)) ShortestAngleDiff(XY(i, 2), XY(j, 2))] .^ 2, 2));
D = squareform(D);

% TODO: Distance of non-neighboring pixels is as usual, but for neighboring
% pixels make it dependent on intensity (n)
% A - area of pixel
% typical distance within points in the pixel: (A / n) ^ 0.5
% between two pixels (approximation): 0.5 * ((A1 / n1) ^ 0.5 + (A2 / n2) ^ 0.5)

D_sqr_int = [ShortestPeriodicDiff(row(i), row(j), size(n, 1)) ShortestPeriodicDiff(col(i), col(j), size(n, 2))] .^ 2;
D_sqr_int = D_sqr_int(:, 1) + D_sqr_int(:, 2);

which_neighbors = (abs(row(i) - row(j)) <= 1) & (abs(col(i) - col(j)) <= 1);

D_sqr_int(which_neighbors) = 0.5 * (...
    n(sub2ind(size(n), row(i(which_neighbors)), col(i(which_neighbors)))) .^ -0.5 + ...
    n(sub2ind(size(n), row(j(which_neighbors)), col(j(which_neighbors)))) .^ -0.5);

D_sqr_int = squareform(D_sqr_int);

assignment = PAM_weighted(XY, w, D, k, 100);

figure(1);
imagesc(centers, centers, p .* (p > 0.0001));
axis xy;
hold on;
for i = 1:max(assignment)
    plot(XY(assignment == i, 1), XY(assignment == i, 2), 'o');
end
hold off;

%%
if (0)
    [~, order] = sort(w, 'descend');
    D_sqr_int = D_sqr_int(order, order);
    XY = XY(order, :);
end

%%
% TODO: Recursively split the domain. Always split the highest (total)
% intensity first


%%
%Z = linkage(D, 'average');
Z = linkage(D_sqr_int, 'average');
%Z = linkage_weighted(D, w);
assignment = cluster(Z, 'maxclust', 11);


figure(1);
imagesc(centers, centers, p .* (p > 0.0001));
axis xy;

hold on;
for i = 1:max(assignment)
    plot(XY(assignment == i, 1), XY(assignment == i, 2), 'o');
end
hold off;

%%

assignment = zeros(size(n));
ignore = (p < 0.0001);

n_vec = n(:);
[i_vec, j_vec] = ind2sub(size(n), 1:numel(n));
ij = [i_vec(:), j_vec(:)];
clear i_vec j_vec;

ij(ignore(:), :) = []; % Remove the pixels to be ignored
n_vec(ignore(:)) = []; 

is_assigned = false(size(n_vec));

[~, m] = max(n_vec);
is_assigned(m) = true;

assignment(ij(m, 1), ij(m, 2)) = 1;

for iteration = 1:10
    [~, m] = max(n_vec .* ~is_assigned);
    is_assigned(m) = true;

    assignment(ij(m, 1), ij(m, 2)) = iteration;
end

kernel = [0 1 0; 1 1 1; 0 1 0];

%assignment = conv2(assignment, kernel, 'same');

figure(2);
imagesc(assignment - ignore);
colorbar;
axis xy;



%%
x = all_frames(:, 1, :);
figure(1);
histogram(x)

[edges] = CalculateEqualProbabilityEdges(x, 5, 51);

hold on;
for i = 2:numel(edges)-1
    plot(edges(i) * [1 1], [0 15e5], '--g');
end
hold off;

%%
%[~, x] = CalculateEqualProbabilityEdges(all_frames(:, 1, :), 11, 51);
%[~, y] = CalculateEqualProbabilityEdges(all_frames(:, 2, :), 11, 51);
%figure(2);
%histogram(x(:))

%figure(3);
%histogram(y(:))

%%
e1 = Equalize(all_frames(:, 1, :), 51);
e2 = Equalize(all_frames(:, 2, :), 51);

%%
for i = 0:27
    rama_filepath(i)
    load(rama_filepath(i));    
    
    data(:, :, 1) = (e1.CoordinateToPosition(data(:, :, 1)) - 0.5) * (2 * pi());
    data(:, :, 2) = (e2.CoordinateToPosition(data(:, :, 2)) - 0.5) * (2 * pi());

    save(out_filepath(i), 'data');    
end


%%
clear e;

for i = 1:size(all_frames, 1)
    i
    e(i, 1) = Equalize(all_frames(i, 1, :), 51);
    e(i, 2) = Equalize(all_frames(i, 2, :), 51);
end

%%
out2_filepath = @(i)['c:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Sosnick\Data\20190308\' ...
    sprintf('%s\\%s.%d.eql2.mat', sim_name, sim_name, i)];

for fi = 0:27
    rama_filepath(fi)
    load(rama_filepath(fi));    
    
    for i = 1:size(data, 2)
        data(:, i, 1) = (e(i, 1).CoordinateToPosition(data(:, i, 1)) - 0.5) * (2 * pi());
        data(:, i, 2) = (e(i, 2).CoordinateToPosition(data(:, i, 2)) - 0.5) * (2 * pi());
    end

    save(out2_filepath(fi), 'data');    
end
