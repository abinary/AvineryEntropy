function [] = SetMyDefaultFigureSettings()

ax = gca();
ax.LineWidth = 1.5;
ax.FontName = 'Arial';
ax.FontWeight = 'bold';
%ax.FontSize = 14;
ax.FontUnits = 'normalized';
ax.FontSize = 0.06;
ax.FontSmoothing = 'on';
ax.TitleFontSizeMultiplier = 1.5;
ax.TitleFontWeight = 'bold';
ax.LabelFontSizeMultiplier = 1.3;

end
