function [rmsd, ij] = AllPairsRmsd(X)

% Code is based on the "procrustes" matlab function

persistent kernel pkernel;

if (isempty(pkernel))
    pkernel = parallel.gpu.CUDAKernel('procrustes.ptx', 'procrustes.cu');
    %kernel = parallel.gpu.CUDAKernel('svd.ptx', 'svd.cu');
end

if (size(X, 3) == 3)
    X = permute(X, [2 3 1]); % frame index should be last
end

num_items = size(X, 3);

% Assuming a 3d array is given

%save('AllPairsRmsd.debug.mat');
%load('AllPairsRmsd.debug.mat');

%%
%X = gpuArray(X);

% Center at the origin.
mu = mean(X, 1);
X0 = X - mu;

ssqX = sum(X0 .^ 2, 1);
ssqX = sum(ssqX, 2);

normX = sqrt(ssqX); % == sqrt(trace(X0*X0'))

% Scale to equal (unit) norm.
X1 = X0 ./ normX;

%%
ij = AllPairs(num_items, 1);

%rmsd = arrayfun(@(k)CalculatePairRmsd(k), 1:size(ij, 1));
%return;

rmsd = zeros(size(ij, 1), 1);

% Setup the data
if (1)
    x = permute(X1, [3 2 1]);
    y = X1;
    A_all = reshape(x, [size(x, 1)*size(x, 2), size(x, 3)]) * reshape(y, [size(y, 1) , size(y, 2)*size(y, 3)]);
    clear x y;
    A_all = reshape(A_all, [size(A_all, 1)/3 3 3 size(A_all, 2)/3]);
    A_all = permute(A_all, [1 4 3 2]);
    A_all = reshape(A_all, [numel(A_all) / 9 3 3]);
end

% Compute transformation matrices (using SVD, as in "procrustes")
if (1)
    n = int32(size(A_all, 1));

    %%
    if (0)
        kernel.ThreadBlockSize = kernel.MaxThreadsPerBlock;
        kernel.GridSize = [ceil(n / kernel.MaxThreadsPerBlock) 1 1];
        
        %[~, T_all] = feval(kernel, gpuArray(single(A_all)), zeros(21, n), n);
        [~, T_all] = feval(kernel, gpuArray(single(permute(A_all, [1 2 3]))), zeros(n, 21), n);
        T_all = reshape(gather(T_all), size(T_all));
        L_all = reshape(T_all(:, 13:21), [n 3 3]);
        M_all = reshape(T_all(:, 1:9), [n 3 3]);
        
        L_all = permute(L_all, [3 2 1]);
        M_all = permute(M_all, [3 2 1]);
        %squeeze(T_all(2, :, :))
    end

    if (1)
        pkernel.ThreadBlockSize = pkernel.MaxThreadsPerBlock;
        pkernel.GridSize = [ceil(n / pkernel.MaxThreadsPerBlock) 1 1];
        [~, T_all] = feval(pkernel, gpuArray(single(A_all)), zeros(n, 9), n);
        T_all = reshape(gather(T_all), [n 3 3]);
        T_all = permute(T_all, [2 3 1]);
    end
end

for k = 1:size(ij, 1)
    i = ij(k, 1);
    j = ij(k, 2);
    

    if (0)
        Y = X1(:, :, j);
        %A = squeeze(X0(i, :, :))' * Y;
        A = X1(:, :, i)' * Y;
        [L, ~, M] = svd(A);
        
        %L = L_all(:, :, (i - 1) * num_items + j);
        %M = M_all(:, :, (i - 1) * num_items + j);
        T = M * L';
        
        haveReflection = (det(T) < 0);
        % If we don't have what was asked for ...
        if (haveReflection)
            % ... then either force a reflection, or undo one.
            M(:,end) = -M(:,end);
            %D(end,end) = -D(end,end);
            T = M * L';
        end
    end
    
    T = T_all(:, :, (i - 1) * num_items + j);
    
%     if (abs(abs(T) - abs(T2)) > 1e-5)
%         1;
%     end
    
    % The minimized unstandardized distance D(X0,b*Y0*T) is
    % ||X0||^2 + b^2*||Y0||^2 - 2*b*trace(T*X0'*Y0)
    %traceTA = sum(diag(D)); % == trace(sqrtm(A'*A)) when doReflection is 'best'
    
    % The standardized distance between X and Y*T+c.
    %d = 1 + ssqX(j)/ssqX(i) - 2*traceTA*normX(j)/normX(i);
    
    % Restore original scale, rotate to match and shift to same origin as
    % the compared coordinates
    Z = X0(:, :, j) * T;
    x = X0(:, :, i);
    
    rmsd(k) = sqrt(sum((x(:) - Z(:)) .^ 2) ./ size(X0, 1));
end

1;
end
