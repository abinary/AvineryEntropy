function [] = SetPlotRectRatio(ratio)
f = gcf();
pos = f.Position;
pos(3) = ratio * pos(4);
f.Position = pos;
end
