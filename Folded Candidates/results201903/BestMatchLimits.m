function [L, U] = BestMatchLimits(frames, nd)
L = zeros(1, nd);
U = zeros(1, nd);

for d = 1:nd
    x = frames(:, :, d);
    x = x(:);
    l = min(x);
    u = max(x);
    
    l0 = [-pi -1 0];
    [~, ti] = min(abs(l - l0));
    
    if ((l0(ti) - l) > 1e-6)
        warning('The data exceeds the best matching limit.');
    end
    
    L(d) = l0(ti);
    
    u0 = [1 pi 2*pi];
    [~, tj] = min(abs(u - u0));
    
    if ((u - u0(tj)) > 1e-6)
        warning('The data exceeds the best matching limit.');
    end
    
    U(d) = u0(tj);
    
    if ((ti == 2 && tj ~= 1) || (ti == 1 && tj ~= 2) || (tj == 2 && ti ~= 1) || (tj == 3 && ti ~= 3))
        warning('Best matching lower and upper limit pair are weird.');
    end
end
end
