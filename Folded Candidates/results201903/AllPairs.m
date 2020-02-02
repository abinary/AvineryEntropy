function [i, j] = AllPairs(max_index, order)
% [ij] = AllPairs(max_index)
% [i, j] = AllPairs(max_index)
%

if (nargin < 2)
    order = 0;
end

[i, j] = meshgrid(1:max_index, 1:max_index);
if (order == 0)
    i = triu(i, 1);
    j = triu(j, 1);
else
    i = tril(i, -1);
    j = tril(j, -1);
end

j = j(i ~= 0);
i = i(i ~= 0);

if (nargout < 2)
    i = [i, j];
end

end
