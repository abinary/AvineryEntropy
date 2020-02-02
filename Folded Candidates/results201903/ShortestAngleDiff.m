function [x] = ShortestAngleDiff(x, y)
% x, y - are expected to be between -pi and pi
x = abs(mod(x - y + pi, 2*pi) - pi);

    function [] = debug()
        
        N = 1e6;
        x = (rand([N 1]) - 0.5) * (8 * pi);
        y = (rand([N 1]) - 0.5) * (8 * pi);
        
        t = tic();
        ad1 = abs(atan2(sin(x-y), cos(x-y)));
        toc(t)
        
        %x = mod(x + pi, 2*pi) - pi;
        %y = mod(y + pi, 2*pi) - pi;
        
        %ad2 = [mod(x - y + pi, 2*pi) mod(y - x + pi, 2*pi)];
        
        t = tic();
        ad2 = abs(mod(x - y + pi, 2*pi) - pi);
        toc(t)
        
        figure(1);
        plot(x, y, '.');
        
        figure(2);
        plot(ad1, ad2, '.');
        hold on;
        plot([0 pi], [0 pi], '--g');
        hold off;
    end

end