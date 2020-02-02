function [frames, nd] = LoadTrajectoryFileAsDihedralAngles(filepath)

[folder, filename, ext] = fileparts(filepath);

if (strcmpi(ext, '.mat'))
    data = load(filepath);
    
    fn = fieldnames(data);
    if (numel(fn) ~= 1)
        error('Trajectory file contains more than one dataset');
    end
    
    frames = data.(fn{1});
    
    if (size(frames, 1) < size(frames, 2))
        warning('Trajectory''s first dimension is longer than 2nd. Assuming the 1st is the frame dimension.');
    end
    
    nd = size(frames, 3);
    if (nd ~= 1 && nd ~= 2)
        warning('Unexpected number of dimensions in data. Expecting 1 or 2.');
    end
    
%     if (ndims(frames) ~= 3 || size(frames, 3) ~= 2 || size(frames, 1) < size(frames, 2))
%         error('Trajectory should have 3 dimensions: [frames]x[sites]x[2 angles]');
%     end
else
    error('Unrecognized trajectory file type');
end

1;

end
