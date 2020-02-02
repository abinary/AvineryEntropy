function [x] = ShortestPeriodicDiff(x, y, period)
x = abs(mod(x - y + period * 0.5, period) - 0.5 * period);
