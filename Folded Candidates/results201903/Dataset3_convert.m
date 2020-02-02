%%
root = 'd:\Files\Workspace\Sosnick\dataset 3 - 20191222\';

NET.addAssembly('C:\Files\The Lab\Projects\Bridging Physical and Information Entropy\Simulation\PhysiCasino\DcdSharp\bin\x64\Release\DcdSharp.exe');

[pdb_filepaths, filenames] = ListFiles([root 'top\*.pdb']);

%%
for pdb_idx = 1:numel(filenames)
    [~, sim_name] = fileparts(filenames{pdb_idx});
    sim_name
    pdbStruct = pdbread(pdb_filepaths{pdb_idx});
    which_ca = cellfun(@(c)strcmp(c, 'CA'), {pdbStruct.Model.Atom.AtomName});

    for i = 0:27
        i
        dcd = DcdSharp.DCD();
        dcd.Load(sprintf('%sdcd\\%s-%d.dcd', root, sim_name, i));
        numFramesRead = dcd.ReadFrames(dcd.NumOfFramesLeftToRead);
        frames = double(dcd.GetFramesAs3dSinglesArray());
        
        %output_filepath = sprintf('%sdcd\\%s-%d.dcd.mat', root, sim_name, i)
        %save(output_filepath, 'frames');
        
        frames = frames(:, which_ca, :);
        output_filepath = sprintf('%sdcd\\%s-%d.CA.mat', root, sim_name, i)
        save(output_filepath, 'frames');
    end
end

