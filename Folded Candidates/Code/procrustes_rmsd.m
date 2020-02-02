function [rmsd, Z, transform] = procrustes_rmsd(X, Y)
% [d, Z, transform] = procrustes_rmsd(X, Y)

[n, m]   = size(X);
[ny, my] = size(Y);

if ny ~= n
    error(message('stats:procrustes:InputSizeMismatch'));
elseif my > m
    error(message('stats:procrustes:TooManyColumns'));
end

% Center at the origin.
muX = mean(X,1);
muY = mean(Y,1);
X0 = X - repmat(muX, n, 1);
Y0 = Y - repmat(muY, n, 1);

ssqX = sum(X0.^2,1);
ssqY = sum(Y0.^2,1);
ssqX = sum(ssqX);
ssqY = sum(ssqY);

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

if (1) % Should not encounter reflections
    haveReflection = (det(T) < 0);
    
    % If we don't have what was asked for ...
    if (haveReflection)
        % ... then either force a reflection, or undo one.
        M(:,end) = -M(:,end);
        D(end,end) = -D(end,end);
        T = M * L';
    end
end

% The minimized unstandardized distance D(X0,b*Y0*T) is
% ||X0||^2 + b^2*||Y0||^2 - 2*b*trace(T*X0'*Y0)
traceTA = sum(diag(D)); % == trace(sqrtm(A'*A)) when doReflection is 'best'

% The standardized distance between X and Y*T+c.
d = 1 + ssqY/ssqX - 2*traceTA*normY/normX;

%if nargout >= 2
    Z = normY*Y0 * T + repmat(muX, n, 1);
%end

rmsd = sqrt(sum((X(:) - Z(:)) .^ 2) ./ size(X, 1));

if nargout >= 3
    if my < m
        T = T(1:my,:);
    end
    
    c = muX - muY*T;
    transform = struct('T', T, 'c', repmat(c, n, 1));
end

end % function
