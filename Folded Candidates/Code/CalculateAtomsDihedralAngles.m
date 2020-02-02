function [phi, psi, omega] = CalculateAtomsDihedralAngles(atoms)

atom_name = {atoms.AtomName};
which_atoms_to_keep = strcmp(atom_name, 'C') | strcmp(atom_name, 'CA') | strcmp(atom_name, 'N');
atoms = atoms(which_atoms_to_keep);

atom_name = {atoms.AtomName};
atom_resSeq = [atoms.resSeq];
atom_coord = [[atoms.X]' [atoms.Y]' [atoms.Z]'];

i = 1; % atom index

%=== Preprocessing the Data
% Checking to see if the initial pattern is N-CA-C-N. In this case we
% reject this data point as this would mean that both PHI and PSI are
% not sharing the same CA (PHI would not be defined for the first
% residue, while PSI would be).

if (0)
    if ((isequal(atom_name{i}, 'N')) && (isequal(atom_name{i+1}, 'CA')))
        i = i + 2;
    end
end

%=== get residue numbers
min_atom_resSeq = min(atom_resSeq);
max_atom_resSeq = max(atom_resSeq);
offset = 0;
if min_atom_resSeq <= 0
    offset = -(min_atom_resSeq) + 1;
end

%=== compute torsion angles
%numResidues = max(atom_resSeq) + offset;
numResidues = max_atom_resSeq + offset;
phi = nan(numResidues, 1);
psi = phi;
omega = phi;

%=== angles not defined unless 4 consecutive atoms in the chain
while(i <= numel(atoms) - 3)
    switch(atom_name{i})
        case 'C'
            % does the pattern satisfies for the phi torsion angle?
            if (isequal(atom_name{i+1}, 'N') &&...
                    (isequal(atom_name{i+2}, 'CA')) &&...
                    (isequal(atom_name{i+3}, 'C')))
                % Are atoms are from adjacent amino acids? If not then
                % torsion angle does not exist (is NaN)
                if ((atom_resSeq(i) == atom_resSeq(i+1)-1) && ...
                        (atom_resSeq(i+1) == atom_resSeq(i+2)) &&...
                        (atom_resSeq(i+1)== atom_resSeq(i+3)))
                    phi(atom_resSeq(i+1)+offset) = localCalculateTorsionAngle(atom_coord(i:i+3,:));
                end
                
            end
        case  'N'
            % does the pattern satisfies for the psi torsion angle?
            if (isequal(atom_name{i+1}, 'CA') &&...
                    (isequal(atom_name{i+2}, 'C')) &&...
                    (isequal(atom_name{i+3}, 'N')))
                if ((atom_resSeq(i) == atom_resSeq(i+1)) &&...
                        (atom_resSeq(i) == atom_resSeq(i+2)) &&...
                        (atom_resSeq(i) == atom_resSeq(i+3)-1))
                    psi(atom_resSeq(i)+offset) = localCalculateTorsionAngle(atom_coord(i:i+3,:));
                end
            end
        case 'CA'
            %if (nargout > 0) % compute omega only if user asks for output
            % does the pattern satisfies for the omega torsion angle?
            if (isequal(atom_name{i+1}, 'C') &&...
                    (isequal(atom_name{i+2}, 'N')) &&...
                    (isequal(atom_name{i+3}, 'CA')))
                if ((atom_resSeq(i) == atom_resSeq(i+1)) &&...
                        (atom_resSeq(i) == atom_resSeq(i+2)-1) &&...
                        (atom_resSeq(i) == atom_resSeq(i+3)-1))
                    omega(atom_resSeq(i)+offset) = localCalculateTorsionAngle(atom_coord(i:i+3, :));
                end
            end
            %end
    end
    i = i + 1;
end


if (nargout == 1)
    phi = [phi psi omega];
end

1;

%---------------------------------------------------------------------------

    function [torsionAngle] = localCalculateTorsionAngle(coord)
        % Evaluate the torsion angles.
        
        p1 = coord(1,:);
        p2 = coord(2,:);
        p3 = coord(3,:);
        p4 = coord(4,:);
        
        a = p2 - p1;
        b = p3 - p2;
        c = p4 - p3;
        
        a = a/norm(a,2);
        b = b/norm(b,2);
        c = c/norm(c,2);
        
        b_cross_c = [b(2).*c(3) - b(3).*c(2);
            b(3).*c(1) - b(1).*c(3);
            b(1).*c(2) - b(2).*c(1)];
        
        x = -sum(conj(a).*c) + ((sum(conj(a).*b))*(sum(conj(b).*c)));
        y = sum(conj(a).*b_cross_c');
        torsionAngle = localNewAngle(x,y);
    end

%---------------------------------------------------------------------------

    function ang = localNewAngle(x,y)
        % Calculate the angle represented by (x,y). The angle lies between -pi and
        % pi.
        
        ang = 0; % This is the default value. In case y == 0 and x ~= 0, ang = 0.
        if (x ~= 0) && (y~= 0)
            c = x./sqrt(x.^2 + y.^2);
            ang = sign(y)*acos(c);
        elseif (x == 0)
            if (y > 0)
                ang = pi/2;
            elseif (y < 0)
                ang = -pi/2;
            end
        end
        
    end
end
