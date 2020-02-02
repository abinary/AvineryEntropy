function [rmsd] = CalculateRmsd(xyz1, xyz2)
[~, Z] = procrustes(xyz1, xyz2, 'scaling', false, 'reflection', false);
rmsd = sqrt((sum((xyz1(:) - Z(:)).^2) / size(xyz1, 1)));
end
