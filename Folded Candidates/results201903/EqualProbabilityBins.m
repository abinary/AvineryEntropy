
addpath('splinefit');

%%

x = randn(1e4, 1);

%%
x = sort(x);
cdf = [1:numel(x)] * (1/numel(x));
figure(1);
plot(x, cdf);

%%
figure(1);
% [n, edges] = histcounts(x, 201);
% centers = 0.5 * (edges(1:end-1) + edges(2:end));
% n = cumsum(n) ./ sum(n);
% bar(centers, n);

breaks = 21;
%sp = splinefit(centers, n, breaks);
sp = splinefit(x, cdf, breaks);
xs = linspace(edges(1), edges(end), 1e3);
ys = fnval(sp, xs);

hold on;
plot(xs, ys, '--g');
hold off;

%sp = splinefit(n, centers, breaks);
sp = splinefit(cdf, x, breaks);
ys = linspace(0, 1, 1e3);
xs = fnval(sp, ys);

hold on;
plot(xs, ys, '--r');
hold off;


x0 = linspace(edges(1), edges(end), 1e4);
y0 = fnval(sp, xs);

k = 11;
percentiles = [1/k] * [1:k-1];

edges = [-inf fnval(sp, percentiles) inf];

[n_equal] = histcounts(x, edges);
figure(2);
bar(n_equal);
