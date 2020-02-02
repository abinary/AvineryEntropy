
sim = 'a3d';

Filepath = @(i)['d:\Files\Workspace\Sosnick\dataset 2 - 20190308\' sim '.vtf\' sim sprintf('.%d.vtf.mat', i)];
OutFilepath = @(i)['d:\Files\Workspace\Sosnick\dataset 2 - 20190308\' sim '.vtf\' sim sprintf('.%d.CA.mat', i)];

load(['d:\Files\Workspace\Sosnick\Structures\' sim '_mapped_atoms.mat']);
w = arrayfun(@(mapped_atoms)strcmp(mapped_atoms.name, 'CA'), mapped_atoms);
w = find(w);

ref_pdb = pdbread(['d:\Files\Workspace\Sosnick\Structures\' sim '.pdb'])
ref_xyz = [[ref_pdb.Model.Atom.X]' [ref_pdb.Model.Atom.Y]' [ref_pdb.Model.Atom.Z]'];
ref_ca_xyz = ref_xyz([mapped_atoms(w).index_in_pdb], :);

save(['d:\Files\Workspace\Sosnick\Structures\' sim '.ref_ca_xyz.mat'], 'ref_ca_xyz');

%%
for i = 0:27
    i
    frames = load(Filepath(i));
    frames = frames.frames;
    frames = frames(:, w, :);
    
    save(OutFilepath(i), 'frames');
end
