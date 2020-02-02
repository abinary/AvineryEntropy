function [data] = LoadMatrixFromFile(file_path)

data = load(file_path);

if (isstruct(data)) % could be a mat binary file
    fn = fieldnames(data);
    
    which_are_matrices = cellfun(@(f)ismatrix(data.(f)), fn);
    
    assert(nnz(which_are_matrices) == 1, 'Found more than one matrix in the loaded file. The file must contain only one matrix.');
    
    data = data.(fn{which_are_matrices});
end

end
