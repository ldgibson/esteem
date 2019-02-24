function Vnn = nucnucrepulsion(atoms, xyz_a0)
% Vnn = nucnucrepulsion(atoms, xyz_a0)
%
% Input:
%   atoms   list of element numbers (array with K elements); e.g. [6 8] for CO
%   xyz_a0  K×3 array of Cartesian coordinates of nuclei, in bohr
% Output:
%   Vnn     total nuclear repulsion energy, in hartrees

Vnn = 0;
for i = 1:numel(atoms)-1
    for j = i+1:numel(atoms)
        Vnn = Vnn + atoms(i)*atoms(j)/norm(xyz_a0(i, :) - xyz_a0(j, :));
    end
end
end