function [] = SetPlot2to1Rect()
f = gcf();
pos = f.Position;
pos(3) = 2 * pos(4);
f.Position = pos;
end
