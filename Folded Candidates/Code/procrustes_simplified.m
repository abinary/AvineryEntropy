function [d, Z, transform] = procrustes(X, Y, varargin)
%PROCRUSTES Procrustes Analysis
%   D = PROCRUSTES(X, Y) determines a linear transformation (translation,
%   reflection, orthogonal rotation, and scaling) of the points in the
%   matrix Y to best conform them to the points in the matrix X.  The
%   "goodness-of-fit" criterion is the sum of squared errors.  PROCRUSTES
%   returns the minimized value of this dissimilarity measure in D.  D is
%   standardized by a measure of the scale of X, given by
%
%      sum(sum((X - repmat(mean(X,1), size(X,1), 1)).^2, 1))
%
%   i.e., the sum of squared elements of a centered version of X.  However,
%   if X comprises repetitions of the same point, the sum of squared errors
%   is not standardized.
%
%   X and Y are assumed to have the same number of points (rows), and
%   PROCRUSTES matches the i'th point in Y to the i'th point in X.  Points
%   in Y can have smaller dimension (number of columns) than those in X.
%   In this case, PROCRUSTES adds columns of zeros to Y as necessary.
%
%   [D, Z] = PROCRUSTES(X, Y) also returns the transformed Y values.
%
%   [D, Z, TRANSFORM] = PROCRUSTES(X, Y) also returns the transformation
%   that maps Y to Z.  TRANSFORM is a structure with fields:
%      c:  the translation component
%      T:  the orthogonal rotation and reflection component
%      b:  the scale component
%   That is, Z = TRANSFORM.b * Y * TRANSFORM.T + TRANSFORM.c.
%
%   [...] = PROCRUSTES(..., 'Scaling',false) computes a procrustes solution
%   that does not include a scale component, that is, TRANSFORM.b == 1.
%   PROCRUSTES(..., 'Scaling',true) computes a procrustes solution that
%   does include a scale component, which is the default.
%
%   [...] = PROCRUSTES(..., 'Reflection',false) computes a procrustes solution
%   that does not include a reflection component, that is, DET(TRANSFORM.T) is
%   1.  PROCRUSTES(..., 'Reflection','best') computes the best fit procrustes
%   solution, which may or may not include a reflection component, 'best' is
%   the default.  PROCRUSTES(..., 'Reflection',true) forces the solution to
%   include a reflection component, that is, DET(TRANSFORM.T) is -1.
%
%   Examples:
%
%      % Create some random points in two dimensions
%      n = 10;
%      X = normrnd(0, 1, [n 2]);
%
%      % Those same points, rotated, scaled, translated, plus some noise
%      S = [0.5 -sqrt(3)/2; sqrt(3)/2 0.5]; % rotate 60 degrees
%      Y = normrnd(0.5*X*S + 2, 0.05, n, 2);
%
%      % Conform Y to X, plot original X and Y, and transformed Y
%      [d, Z, tr] = procrustes(X,Y);
%      plot(X(:,1),X(:,2),'rx', Y(:,1),Y(:,2),'b.', Z(:,1),Z(:,2),'bx');
%
%   See also FACTORAN, CMDSCALE.

%   References:
%     [1] Seber, G.A.F., Multivariate Observations, Wiley, New York, 1984.
%     [2] Gower, J.C. and Dijskterhuis, G.B., Procrustes Problems, Oxford
%         Statistical Science Series, Vol 30. Oxford University Press, 2004.
%     [3] Bulfinch, T., The Age of Fable; or, Stories of Gods and Heroes,
%         Sanborn, Carter, and Bazin, Boston, 1855.

%   Copyright 1993-2009 The MathWorks, Inc.


if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end

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

if (0) % Should not encounter reflections
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

if nargout > 1
    Z = normY*Y0 * T + repmat(muX, n, 1);
end

if nargout > 2
    if my < m
        T = T(1:my,:);
    end
    c = muX - muY*T;
    transform = struct('T',T, 'c',repmat(c, n, 1));
end
