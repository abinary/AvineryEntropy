function [vtfAtoms] = map_vtf_coordinates(filename, pdbFilePath)

%filename = 'nug2_0.vtf';
%pdbFilePath = 'c:\Users\Ram\Dropbox\Shared\Tobin''s trajectories\nug2.pdb';

pdb = pdbread(pdbFilePath);
pdbAtoms = [pdb.Model.Atom];

vtfAtoms = struct();

f = fopen(filename, 'r');

l = fgetl(f);
while (~strcmp(l, 'timestep ordered'))
    
    %%
    if (length(l) > 5 && strcmp(l(1:5), 'atom '))
        matches = regexp(l, '(?<field>[\w\d]+)\s+(?<value>[\w\d]+)\s*', 'names');
        a = struct();
        for i = 1:numel(matches)
            a.(matches(i).field) = matches(i).value;
        end
        
        a.atom = str2double(a.atom);
        a.resid = str2double(a.resid);
        
        vtfAtoms(1 + a.atom).name = a.name;
        vtfAtoms(1 + a.atom).residue = a.resname;
        vtfAtoms(1 + a.atom).residue_id = a.resid;
        vtfAtoms(1 + a.atom).in_pdb = false;
        vtfAtoms(1 + a.atom).index_in_pdb = 0;
        
    end
    
    %%
    l = fgetl(f);
end


fclose(f);

for i = 1:numel(vtfAtoms)

    which = ...
        ([pdbAtoms.resSeq] == vtfAtoms(i).residue_id) & ...
        strcmp({pdbAtoms.AtomName}, vtfAtoms(i).name);
    
    which = find(which);
    if (~isempty(which))
        vtfAtoms(i).in_pdb = true;
        vtfAtoms(i).index_in_pdb = which;
    end
    
    1;
end


%end
