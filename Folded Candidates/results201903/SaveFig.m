function [] = SaveFig(filename, hfig, resolution)

if (nargin < 3)
    resolution = 600;
end

if (nargin < 2 || isempty(hfig))
    hfig = gcf();
end

[folder, name, ext] = fileparts(filename);

if (isempty(ext))
    ext = '.tif';
    filename = [filename ext];
end

fileType = ext(2:end);

if strcmpi(fileType, 'eps')
    print(hfig, '-depsc2', filename);
    vers = version();
    if ~strcmp(vers(1:3), '8.4')
        fixPSlinestyle(filename);
    end
elseif strcmpi(fileType, 'pdf')
    %print(hfig, '-dpdf', filename);
    print(hfig, '-dpdf', filename, '-fillpage');
elseif strcmpi(fileType, 'jpg') || strcmpi(fileType, 'jpeg')
    print(hfig, '-djpeg', '-opengl', sprintf('-r%d',resolution), filename);
elseif strcmpi(fileType, 'png')
    print(hfig, '-dpng', '-opengl', sprintf('-r%d',resolution), filename);
elseif strcmpi(fileType, 'tiff') || strcmpi(fileType, 'tif')
    print(hfig, '-dtiff', '-opengl', sprintf('-r%d',resolution), filename);
elseif strcmpi(fileType, 'emf')
    print(hfig, '-dmeta', sprintf('-r%d',resolution), filename);
elseif strcmpi(fileType, 'svg')
    print(hfig, '-dsvg', filename);
else
    err = MException('', ...
        '=====> ERROR: File type %s is not supported. ', fileType);
    throw(err);
end
end
