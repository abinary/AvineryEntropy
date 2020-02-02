
addpath('../../simulation/common code/');

fclose('all');
vtf_files = ListFiles('d:\Files\Workspace\Sosnick\201807\1\*.vtf', true);

for i = 1:numel(vtf_files)
    fprintf('Processing: %s \r\n', vtf_files{i});
    convert_vtf_to_mat(vtf_files{i});
end

