function [assignment_func, Z, p] = GenerateDihedralAngleClusters(file_list, cg_levels, options)

persistent cache;

should_plot = false;

if (nargin == 0)
    options = struct();
    
    options.HistogramGraining = 2^6;
    options.RelativeProbabilityThreshold = 0.005; % relative to the max
    
    options.SparseContributionPerDOFThreshold = 0.01; % in k_B per DOF
    
    options.Method = 'linkage';
    options.SubMethod = 'average';
    
    options.MethodsList = {'linkage', 'PAM', 'decomp1', 'decomp2'};
    
    options.LinkageSubMethodsList = {'average', 'single', 'complete', 'median', 'ward', 'weighted'};
    
    options.DecompCondDir = 1;
    options.DecompMaxAspectRatio = 2;
    
    options.MedoidsMaxIterations = 100;
    options.ShouldAdjustNearestNeighborsDistance = true;
    
    assignment_func = options;
    return;
end

if (nargin < 3 || isempty(options))
    options = GenerateDihedralAngleClusters();
end

assignment_func = [];

histogram_graining = options.HistogramGraining;
threshold = options.RelativeProbabilityThreshold;

is_cached = (~isempty(cache));

if (~is_cached && exist('GenerateDihedralAngleClusters.cache.mat', 'file'))
    fprintf('<GenerateDihedralAngleClusters> Loading cache ...\n');
    cache = load('GenerateDihedralAngleClusters.cache.mat');
    is_cached = true;
end

if (is_cached) % validate cache match
    is_cached = (options.HistogramGraining == cache.options.HistogramGraining);
    is_cached = is_cached && (numel(file_list) == numel(cache.file_list));
    
    if (is_cached)
        for i = 1:numel(file_list)
            if (~strcmp(file_list{i}, cache.file_list{i}))
                is_cached = false;
                break;
            end
        end
    end
end

if (is_cached)
    fprintf('<GenerateDihedralAngleClusters> Using cache \n');
    %x = cache.x;
    %y = cache.y;
    n = cache.n;
    edges = cache.edges;
    centers = cache.centers;
else
    all_frames = [];
    data = [];
    for i = 1:numel(file_list)
        file_list{i}
        
        data = load(file_list{i});
        
        if (isfield(data, 'data'))
            data = data.data;
        elseif (isfield(data, 'frames'))
            data = data.frames;
        end
        
        all_frames = [all_frames; data];
    end
    
    all_frames = permute(all_frames, [2 3 1]);
    
    
    %%
    x = all_frames(:, 1, :);
    y = all_frames(:, 2, :);
    
    clear all_frames; % save on memory
    
    x = x(:);
    y = y(:);
    
    edges = linspace(-pi, pi, histogram_graining + 1);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    edges(1) = -4;
    edges(end) = 4;
    
    x = (x + pi) .* (histogram_graining / (2 * pi));
    x = floor(x) + 1;
    x(x < 1) = 1;
    x(x > histogram_graining) = histogram_graining;
    
    y = (y + pi) .* (histogram_graining / (2 * pi));
    y = floor(y) + 1;
    y(y < 1) = 1;
    y(y > histogram_graining) = histogram_graining;
    
    %x = discretize(x, edges);
    %y = discretize(y, edges);
    
    x = int32(x);
    y = int32(y);
    
    n = accumarray([x y], 1, [histogram_graining histogram_graining], @sum);
    
    cache.n = n;
    cache.edges = edges;
    cache.centers = centers;
    cache.options = options;
    cache.file_list = file_list;
    
    save('GenerateDihedralAngleClusters.cache.mat', 'n', 'edges', 'centers', 'options', 'file_list');
end

%%
if (should_plot)
    figure(1);
    imagesc(centers, centers, n);
    axis xy;
end

if (strcmp(options.Method, 'linkage') || strcmp(options.Method, 'PAM'))
    %%
    p = n / max(n(:));
    [row, col] = find(p > threshold);
    n_vec = n(sub2ind(size(n), row, col)); % weights
    
    [X, Y] = meshgrid(centers, centers);
    XY = [X(sub2ind(size(X), row, col)) Y(sub2ind(size(X), row, col))];
    
    [i, j] = AllPairs(size(row, 1), 1);
    %D = sqrt(sum([ShortestAngleDiff(XY(i, 1), XY(j, 1)) ShortestAngleDiff(XY(i, 2), XY(j, 2))] .^ 2, 2));
    %D = squareform(D);
    
    D_sqr_int = [ShortestPeriodicDiff(row(i), row(j), size(n, 1)) ShortestPeriodicDiff(col(i), col(j), size(n, 2))] .^ 2;
    D_sqr_int = sqrt(D_sqr_int(:, 1) + D_sqr_int(:, 2));
    
    if (options.ShouldAdjustNearestNeighborsDistance)
        % Distance of non-neighboring pixels is as usual, but for neighboring
        % pixels make it dependent on intensity (n)
        % A - area of pixel
        % typical distance within points in the pixel: (A / n) ^ 0.5
        % between two pixels (approximation): 0.5 * ((A1 / n1) ^ 0.5 + (A2 / n2) ^ 0.5)
        
        which_neighbors = (abs(row(i) - row(j)) <= 1) & (abs(col(i) - col(j)) <= 1);
        
        D_sqr_int(which_neighbors) = 0.5 * (...
            n(sub2ind(size(n), row(i(which_neighbors)), col(i(which_neighbors)))) .^ -0.5 + ...
            n(sub2ind(size(n), row(j(which_neighbors)), col(j(which_neighbors)))) .^ -0.5);
    end
    
    D_sqr_int = squareform(D_sqr_int);
end

edges1 = edges;
edges2 = edges;

switch (options.Method)
    case 'linkage'
        %%
        %Z = linkage(D, 'average');
        L = linkage(D_sqr_int, options.SubMethod);
        assignment = cluster(L, 'maxclust', cg_levels);
        
        if (should_plot)
            figure(2);
            imagesc(centers, centers, p .* (p > threshold));
            axis xy;
        end
        1;
        
    case {'PAM', 'medoids'}
        assignment = PAM_weighted([row col], n_vec, D_sqr_int, cg_levels, options.MedoidsMaxIterations);
        
        if (should_plot)
            figure(3);
            imagesc(centers, centers, p .* (p > threshold));
            axis xy;
        end
        1;
        
    case 'dummy'
        assignment = [];
        [Z, edges1, edges2] = DummyDecomp(cg_levels);
        
    case 'decomp0'
        assignment = [];
        Z = BoxDecomp0(cg_levels);
        
    case 'decomp1'
        assignment = [];
        Z = BoxDecomp1(cg_levels);
        1;
        
    case 'decomp2'
        assignment = [];
        Z = BoxDecomp2(cg_levels);
        1;
        
    otherwise
        1;
end

if (should_plot)
    hold on;
end

cluster_convex_hull = cell([cg_levels 1]);
cluster_ref_points = cell([cg_levels 1]);

for i = 1:max(assignment)
    x = XY(assignment == i, 1);
    y = XY(assignment == i, 2);
    
    mx1 = median(x);
    my1 = median(y);
    
    mx2 = median(ShiftAngle(x, -mx1));
    my2 = median(ShiftAngle(y, -my1));
    
    mx = mx1 + mx2;
    my = my1 + my2;
    
    %mx = mx1;
    %my = my1;
    %mx = 0; my = 0;
    
    if (numel(x) <= 10)
        k = [1:numel(x)]';
    else
        x = x + rand(size(x)) * 1e-6; % solves a co-linearity issue with convhull
        y = y + rand(size(y)) * 1e-6;
        %[k,V] = convhull(ShiftAngle(x, -mx), ShiftAngle(y, -my));
        [k,V] = boundary(ShiftAngle(x, -mx), ShiftAngle(y, -my));
    end
    
    cluster_convex_hull{i} = [x(k), y(k)];
    %cluster_ref_points{i} = [x(k(1:end-1)), y(k(1:end-1))];
    cluster_ref_points{i} = [x(:) y(:)];
    
    if (should_plot)
        %plot(ShiftAngle(x, -mx), ShiftAngle(y, -my), 'o');
        plot(x, y, 'o');
        
        %hold on;
        %plot(ShiftAngle(x(k), -mx), ShiftAngle(y(k), -my), '--');
        plot(x(k), y(k), 'yx');
        
        if (0)
            xy = [x(k), y(k)];
            
            d = diff(xy, 1);
            w = abs(d) > pi;
            
            [split_rows, split_cols] = find(w);
            
            segments = [[1; split_rows] , [split_rows + 1; size(xy, 1)]];
            
            for seg = 1:size(segments, 1)
                k = segments(seg, 1):segments(seg, 2);
                plot(xy(k, 1), xy(k, 2), '--');
            end
        end
        
        %hold off;
        %xlim([-pi pi]);
        %ylim([-pi pi]);
    end
end

if (should_plot)
    hold off;
end

if (exist('Z', 'var') && ~isempty(Z))
    assignment_func = @AssignBoxDecomp;
else
    
    [X, Y] = meshgrid(centers, centers);
    Z = reshape(Assign(reshape([X(:) Y(:)], [numel(X) 1 2])), size(X));
    clear X Y;
    assignment_func = @AssignBoxDecomp;
    %assignment_func = @Assign;
end

num_bins = max(Z(:));
p = zeros(1, num_bins);

for z_idx = 1:num_bins
    p(z_idx) = sum(sum(n(Z == z_idx))) / sum(n(:));
end

clear n x y;


    function [x] = ShiftAngle(x, delta)
        x = mod(x + delta + pi, 2*pi) - pi;
    end

    function [point_assignment] = AssignBoxDecomp(points)
        point_assignment = int32(cat(3, discretize(points(:, :, 1), edges1), discretize(points(:, :, 2), edges2)));
        point_assignment = Z(sub2ind(size(Z), point_assignment(:, :, 1), point_assignment(:, :, 2)));
        point_assignment = uint32(point_assignment);
    end

    function [point_assignment] = Assign(points)
        sz = size(points);
        points = reshape(points, [sz(1)*sz(2), 2]);
        
        all_d = zeros(size(points, 1), cg_levels);
        
        for k = 1:cg_levels
            ref = cluster_ref_points{k};
            
            if (isempty(ref))
                all_d(:, k) = inf;
            else
                ref = shiftdim(ref, 1);
                ref = reshape(ref, [1 size(ref)]);
                d = min(sqrt(sum(ShortestAngleDiff(points, ref) .^ 2, 2)), [], 3);
                %d = max(sqrt(sum(ShortestAngleDiff(points, ref) .^ 2, 2)), [], 3);
                %d = mean(sqrt(sum(ShortestAngleDiff(points, ref) .^ 2, 2)), 3);
                
                all_d(:, k) = d;
            end
        end
        
        [~, point_assignment] = min(all_d, [], 2);
        point_assignment = uint32(reshape(point_assignment, sz(1:2)));
    end

    function [Z, edges1, edges2] = DummyDecomp(k)
        k1 = 1 + mod(k - 1, 16);
        k2 = ceil(k / 16);
        
        edges1 = linspace(-pi, pi, k1 + 1);
        edges2 = linspace(-pi, pi, k2 + 1);
        
        edges1(1) = -4;
        edges1(end) = 4;
        edges2(1) = -4;
        edges2(end) = 4;
        
        load('dihedral_histograms.mat');
        n = dihedral_histograms{k1, k2};
        p = n ./ sum(n(:));
        
        if (numel(n) <= 1)
            Z = 1;
        else
            if (0)
                sp = sort(p(:));
                csp = cumsum(sp);
                sparse_entropy_estimate = csp .* log([1:numel(csp)]');
                [threshold] = find(sparse_entropy_estimate >= (options.SparseContributionPerDOFThreshold), 1);
                threshold = sp(threshold);
            else
                threshold = 0;
            end
            
            w = ~(p <= threshold);
            %w = true(size(n)); % do not use a threshold
            
            Z = zeros(size(n));
            
            if (0) % it seems like, for some reason, assigning numers in order of probability yields higher values
                [~, order] = sort(p(w), 'descend');
                order(order) = 1:nnz(w); % set each to its order index
                Z(w) = order;
            else
                Z(w) = 1:nnz(w);
            end
            
            Z(~w) = nnz(w) + 1;
        end
        
        %[I, J] = meshgrid(1:size(n, 1));
        
        %d1 = floor(size(n, 1) / k1);
        %d2 = floor(size(n, 1) / k2);
        
        %Z = ceil(I / d1) + (ceil(J / d2) - 1) * k1;
        %Z(Z > k) = k;
    end

    function [Z] = BoxDecomp0(k)
        boxes = [1 1 size(n)];
        
        BoxSum = @(box)box(3)*box(4);
        
        for split_idx = 1:k-1
            
            % Choose the box containing the max intensity
            b = arrayfun(@(b)BoxSum(boxes(b, :)), 1:size(boxes, 1));
            %b = arrayfun(@(b)max(max(n(boxes(b, 1):boxes(b, 3)-1, boxes(b, 2):boxes(b, 4)-1))), 1:size(boxes, 1));
            
            b((boxes(:, 3) == 1) & (boxes(:, 4) == 1)) = 0; % ignore boxes which cannot be split
            [~, b] = max(b);
            
            % split either horizontally or vertically, whichever is more equalized
            
            box = boxes(b, :);
            dv = [];
            dh = [];
            clear bv1 bv2 bh1 bh2;
            
            if (box(3)/box(4) < 1)
                bv1 = [box(1:3) box(4) / 2];
                bv2 = [box(1) box(2)+box(4)/2 box(3) box(4)/2];
                boxes = [boxes(1:b-1, :); boxes(b+1:end, :); bv1; bv2];
            else
                bh1 = [box(1:2) box(3)/2 box(4)];
                bh2 = [box(1)+box(3)/2 box(2) box(3)/2 box(4)];
                boxes = [boxes(1:b-1, :); boxes(b+1:end, :); bh1; bh2];
            end
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
        
        % assignment matrix
        Z = zeros(size(n));
        for b = 1:size(boxes, 1)
            Z(boxes(b, 1)+[1:boxes(b, 3)]-1, boxes(b, 2)+[1:boxes(b, 4)]-1) = b;
        end
    end

    function [Z] = BoxDecomp1(k)
        boxes = [1 1 size(n)];
        
        BoxSum = @(box)sum(sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1)));
        
        max_aspect_ratio = options.DecompMaxAspectRatio;
        cond_dir = options.DecompCondDir;
        
        for split_idx = 1:k-1
            
            % Choose the box containing the max intensity
            b = arrayfun(@(b)BoxSum(boxes(b, :)), 1:size(boxes, 1));
            %b = arrayfun(@(b)max(max(n(boxes(b, 1):boxes(b, 3)-1, boxes(b, 2):boxes(b, 4)-1))), 1:size(boxes, 1));
            
            b((boxes(:, 3) == 1) & (boxes(:, 4) == 1)) = 0; % ignore boxes which cannot be split
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
        
        % assignment matrix
        Z = zeros(size(n));
        for b = 1:size(boxes, 1)
            Z(boxes(b, 1)+[1:boxes(b, 3)]-1, boxes(b, 2)+[1:boxes(b, 4)]-1) = b;
        end
    end

    function [Z] = BoxDecomp2(k) % unequal splits
        
        is_cached = false;
        
        if (isfield(cache, 'decomp2_n'))
            is_cached = all(size(n) == size(cache.decomp2_n));
        end
        
        is_cached = is_cached && all(n(:) == cache.decomp2_n(:));
        
        if (~is_cached)
            cache.decomp2_box_sets = {};
            cache.decomp2_n = n;
        end
        
        boxes = [1 1 size(n)];
        
        which_are_cached = ~cellfun(@isempty, cache.decomp2_box_sets);
        cached_k = find(which_are_cached, 1, 'last');
        
        if (cached_k >= k)
            cached_k = k;
        end
        
        if (cached_k >= 2)
            boxes = cache.decomp2_box_sets{cached_k};
        end
        
        BoxSum = @(box)sum(sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1)));
        BoxVCumSum = @(box)cumsum(sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1), 2), 1);
        BoxHCumSum = @(box)cumsum(sum(n(box(1)+[1:box(3)]-1, box(2)+[1:box(4)]-1), 1), 2);
        
        for split_idx = size(boxes, 1):k-1
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
                    
                    %idx = round(interp1(bcs ./ bs, 0.5:box(3), 0.5));
                    idx = find(bcs ./ bs > 0.5, 1);
                    
                    if (idx == numel(bcs))
                        idx = idx - 1;
                    end
                    
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
                    %idx = round(interp1(bcs ./ bs, 0.5:box(4), 0.5));
                    idx = find(bcs ./ bs > 0.5, 1);
                    
                    if (idx == numel(bcs))
                        idx = idx - 1;
                    end
                    
                    if (idx < 2)
                        idx = 2;
                    elseif (idx > box(4)-1)
                        idx = box(4)-1;
                    end
                end
                
                b1 = [box(1:3) idx];
                b2 = [box(1) box(2)+idx box(3) box(4)-idx];
            end
            
            if (any(isnan(b1)) || any(isnan(b2)))
                1;
            end
            
            boxes = [boxes(1:b-1, :); boxes(b+1:end, :); b1; b2];
            
            cache.decomp2_box_sets{size(boxes, 1)} = boxes;
        end
        
        if (0)
            %%
            figure(1)
            imagesc(centers, centers, n);
            axis xy;
            imagesc(n);
            axis ij;
            
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
        end
        
        % assignment matrix
        Z = zeros(size(n));
        for b = 1:size(boxes, 1)
            Z(boxes(b, 1)+[1:boxes(b, 3)]-1, boxes(b, 2)+[1:boxes(b, 4)]-1) = b;
        end
        
        %         figure(2);
        %         imagesc(centers, centers, Z);
        %         axis xy;
    end
end
