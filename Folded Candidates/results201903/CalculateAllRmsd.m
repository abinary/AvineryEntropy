function [rmsd, ij, T_all, adjusted] = CalculateAllRmsd(ref, Y, frame_idx, dof_idx, coord_idx)

% Assuming a 3d array is given

%save('AllPairsRmsd.debug.mat');
%load('AllPairsRmsd.debug.mat');

persistent kernel pkernel;

if (isempty(pkernel))
    pkernel = parallel.gpu.CUDAKernel('procrustes.ptx', 'procrustes.cu');
    %kernel = parallel.gpu.CUDAKernel('svd.ptx', 'svd.cu');
end

if (nargin < 5)
    sz = size(ref);
    
    [~, frame_idx] = max(sz);
    [~, coord_idx] = min(sz);
    
    dof_idx = [1 2 3];
    dof_idx([frame_idx coord_idx]) = [];
end

if (1)
    ref = permute(ref, [dof_idx coord_idx frame_idx]); % frame index should be last
    Y = permute(Y, [dof_idx coord_idx frame_idx]); % assume "Y" is structured the same
    
    frame_idx = 3;
    coord_idx = 2;
    dof_idx = 1;
end


%%
%X = gpuArray(X);

% Center at the origin.
muX = mean(ref, dof_idx);
X0 = ref - muX;

muY = mean(Y, dof_idx);
Y0 = Y - muY;

ssqX = sum(X0 .^ 2, dof_idx);
ssqX = sum(ssqX, coord_idx);

ssqY = sum(Y0 .^ 2, dof_idx);
ssqY = sum(ssqY, coord_idx);

normX = sqrt(ssqX); % == sqrt(trace(X0*X0'))
normY = sqrt(ssqY); % == sqrt(trace(X0*X0'))

% Scale to equal (unit) norm.
X1 = X0 ./ normX;
Y1 = Y0 ./ normY;

%%
[i, j] = meshgrid([1:size(X0, frame_idx)], [1:size(Y0, frame_idx)]);
ij = cat(3, i, j);
ij = reshape(ij, [size(X0, frame_idx)*size(Y0, frame_idx), 2]);

rmsd = zeros(size(ij, 1), 1);

adjusted = [];
if (nargout >= 4)
    adjusted = zeros(size(X0, dof_idx), 3, size(ij, 1));
end

% Setup the data
if (1)
    x = permute(X1, [frame_idx coord_idx dof_idx]);
    y = Y1;
    A_all = reshape(x, [size(x, 1)*size(x, 2), size(x, 3)]) * reshape(y, [size(y, 1) , size(y, 2)*size(y, 3)]);
    clear x y;
    A_all = reshape(A_all, [size(A_all, 1)/3 3 3 size(A_all, 2)/3]);
    A_all = permute(A_all, [1 4 3 2]);
    A_all = reshape(A_all, [numel(A_all) / 9 3 3]);
end

n = int32(size(A_all, 1));

if (1)
    pkernel.ThreadBlockSize = pkernel.MaxThreadsPerBlock;
    pkernel.GridSize = [ceil(double(n) / double(pkernel.MaxThreadsPerBlock)) 1 1];
    [~, T_all] = feval(pkernel, gpuArray(single(A_all)), zeros(n, 9), n);
    T_all = gather(T_all);
    T_all = reshape(T_all, [size(ref, frame_idx) size(Y, frame_idx) 3 3]);
    T_all = permute(T_all, [3 4 1 2]);
end

for k = 1:size(ij, 1)
    i = ij(k, 1);
    j = ij(k, 2);
    
    if (0)
        % The optimum rotation matrix of Y.
        y = squeeze(Y1(:, :, j));
        A = squeeze(X1(:, :, i))' * y;
        %A = squeeze(A_all(i, :, :, j));
        [L, ~, M] = svd(A);
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
    
    %T = T_all(:, :, (i - 1) * size(Y, 3) + j);
    %T = T_all(:, :, (j - 1) * size(X, frame_idx) + i);
    T = T_all(:, :, i, j)';

    % The minimized unstandardized distance D(X0,b*Y0*T) is
    % ||X0||^2 + b^2*||Y0||^2 - 2*b*trace(T*X0'*Y0)
    %traceTA = sum(diag(D)); % == trace(sqrtm(A'*A)) when doReflection is 'best'
    
    % The standardized distance between X and Y*T+c.
    %d = 1 + ssqX(j)/ssqX(i) - 2*traceTA*normX(j)/normX(i);
    
    Z = Y0(:, :, j) * T;
    x = X0(:, :, i);
    
    rmsd(k) = sqrt(sum((x(:) - Z(:)) .^ 2) ./ size(X0, dof_idx));
    
    if (~isempty(adjusted))
        adjusted(:, :, j) = Z;
    end
end

rmsd = reshape(rmsd, [size(Y0, frame_idx) size(X0, frame_idx)])';
%rmsd = reshape(rmsd, [size(X0, frame_idx) size(Y0, frame_idx)]);

1;

    function [d, Z, transform] = procrustes(X, Y, varargin)
        if nargin > 2
            [varargin{:}] = convertStringsToChars(varargin{:});
        end
        
        pnames = {   'scaling'  'reflection'};
        dflts =  {       true         'best'};
        [doScaling,doReflection] = internal.stats.parseArgs(pnames, dflts, varargin{:});
        
        if ~isscalar(doScaling) || ~(islogical(doScaling) || isnumeric(doScaling))
            error(message('stats:procrustes:BadScaling'));
        end
        if isequal(doReflection,'best')
            doReflection = [];
        elseif ~isscalar(doReflection) || ~(islogical(doReflection) || isnumeric(doReflection))
            error(message('stats:procrustes:BadReflection'));
        end
        
        [n, m]   = size(X);
        [ny, my] = size(Y);
        
        if ny ~= n
            error(message('stats:procrustes:InputSizeMismatch'));
        elseif my > m
            error(message('stats:procrustes:TooManyColumns'));
        end
        
        % Center at the origin.
        mu = mean(X,1);
        muY = mean(Y,1);
        X0 = X - repmat(mu, n, 1);
        Y0 = Y - repmat(muY, n, 1);
        
        ssqX = sum(X0.^2,1);
        ssqY = sum(Y0.^2,1);
        constX = all(ssqX <= abs(eps(class(X))*n*mu).^2);
        constY = all(ssqY <= abs(eps(class(X))*n*muY).^2);
        ssqX = sum(ssqX);
        ssqY = sum(ssqY);
        
        if ~constX && ~constY
            % The "centered" Frobenius norm.
            normX = sqrt(ssqX); % == sqrt(trace(X0*X0'))
            normY = sqrt(ssqY); % == sqrt(trace(Y0*Y0'))
            
            % Scale to equal (unit) norm.
            X0 = X0 / normX;
            Y0 = Y0 / normY;
            
            % Make sure they're in the same dimension space.
            if my < m
                Y0 = [Y0 zeros(n, m-my)];
            end
            
            % The optimum rotation matrix of Y.
            A = X0' * Y0;
            [L, D, M] = svd(A);
            T = M * L';
            if isempty(doReflection) % 'best'
                % Let the data decide if a reflection is needed.
            else
                haveReflection = (det(T) < 0);
                % If we don't have what was asked for ...
                if (doReflection ~= haveReflection)
                    % ... then either force a reflection, or undo one.
                    M(:,end) = -M(:,end);
                    D(end,end) = -D(end,end);
                    T = M * L';
                end
            end
            
            % The minimized unstandardized distance D(X0,b*Y0*T) is
            % ||X0||^2 + b^2*||Y0||^2 - 2*b*trace(T*X0'*Y0)
            traceTA = sum(diag(D)); % == trace(sqrtm(A'*A)) when doReflection is 'best'
            
            if doScaling
                % The optimum scaling of Y.
                b = traceTA * normX / normY;
                
                % The standardized distance between X and b*Y*T+c.
                d = 1 - traceTA.^2;
                
                if nargout > 1
                    Z = normX*traceTA * Y0 * T + repmat(mu, n, 1);
                end
            else % if ~doScaling
                b = 1;
                
                % The standardized distance between X and Y*T+c.
                d = 1 + ssqY/ssqX - 2*traceTA*normY/normX;
                
                if nargout > 1
                    Z = normY*Y0 * T + repmat(mu, n, 1);
                end
            end
            
            if nargout > 2
                if my < m
                    T = T(1:my,:);
                end
                c = mu - b*muY*T;
                transform = struct('T',T, 'b',b, 'c',repmat(c, n, 1));
            end
            
            % The degenerate cases: X all the same, and Y all the same.
        elseif constX
            d = 0;
            Z = repmat(mu, n, 1);
            T = eye(my,m);
            transform = struct('T',T, 'b',0, 'c',Z);
        else % ~constX & constY
            d = 1;
            Z = repmat(mu, n, 1);
            T = eye(my,m);
            transform = struct('T',T, 'b',0, 'c',Z);
        end
    end
end
