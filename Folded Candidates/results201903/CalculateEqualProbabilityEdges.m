function [edges, x] = CalculateEqualProbabilityEdges(x, nbins, num_spline_breaks)
addpath('splinefit');

if (nargin < 3)
    num_spline_breaks = min([floor(numel(x) / 5) 21]);
end

sz = size(x);
x = x(:);


y = sort(x);
cdf = [1:numel(y)]' * (1/numel(y));

% Limit the process to 10,000 points
if (numel(x) > 1e4)
    [n, edges] = histcounts(double(x), 1e4);
    y = 0.5 * (edges(1:end-1) + edges(2:end));
    cdf = cumsum(n) * (1 ./ numel(x));
end

sp = splinefit(cdf, y, num_spline_breaks);

if (nargin >= 2 || isempty(nbins))
    percentiles = [1/nbins] * [1:nbins-1];
    edges = [-inf fnval(sp, percentiles) inf];
else
    edges = [];
end

if (nargout >= 2)
    sp = splinefit(y, cdf, num_spline_breaks);
    x = reshape(fnval(sp, x), sz);
end

end
