function [data] = LoadCartesianTrajectoryFile(file_path)

[folder, name, ext] = fileparts(file_path);

switch (ext)
    case '.vtf'
        error('VTF support is not yet implemented');
        1;
        
    case '.mat'
        data = load(file_path);
        
        if (isstruct(data)) % could be a mat binary file
            fn = fieldnames(data);
            
            which_are_trajectories = cellfun(@(f)isnumeric(data.(f)) && (ndims(data.(f)) == 3) && (size(data.(f), 3) == 3), fn);
            
            assert(nnz(which_are_trajectories) == 1, 'Found more than one matrix in the loaded file. The file must contain only one matrix.');
            
            data = data.(fn{which_are_trajectories});
        end
        
    otherwise
        error('Unrecognized file extension: "%s"', ext);
end


end
