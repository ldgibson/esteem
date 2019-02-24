function basis = buildbasis(atoms, xyz_a0, basissetdef)
% basis = buildbasis(atoms, xyz_a0, basissetdef)
% 
% Input:
%     atoms list of element numbers (array with K elements); e.g. [6 8] for CO
%     xyz_a0 K×3 array of Cartesian coordinates of nuclei, in bohr
%     basissetdef basis set definitions, as read in by basisread
% Output:
%     basis M-element structure array, where each element contains:
%         basis(p).atom element number (6 for C, 8 for O, etc.)
%         basis(p).A vector of Cartesian coordinates of the nucleus, in bohr
%         basis(p).a vector of Cartesian exponents ([0 0 0] for s, [1 0 0] for px, etc.)
%         basis(p).alpha array of radial exponents of primitives, in inverse bohr
%         basis(p).d array of contraction coefficients
%         basis(p).N array of normalization constants

basis = {};
for id = 1:numel(atoms)
    shells = unique([basissetdef{atoms(id)}.shelltype]);
    for shell_id = 1:numel(shells)
        basis = [basis getShellBasisFunctions(atoms(id), xyz_a0(id, :), ...
                                              basissetdef, shells(shell_id))];
    end
end
basis = cell2mat(basis);
end


function basis = getShellBasisFunctions(atom, xyz_a0, basissetdef, shelltype)
% Extracts the basis functions for a specified shell type, e.g., S, SP, P,
% D, etc. In the case of shell type SP, it returns in the order (s, px, py,
% pz) for each SP basis function.

reducedbasisdef = getReducedBasisDef(atom, shelltype, basissetdef);
basis_template = struct('atom', [], 'A', zeros(1, 3), ...
                        'a', zeros(1, 3), 'alpha', [], 'd', [], 'N', []);
basis = {};

if strcmp(shelltype, 'S')
    for M = 1:numel(reducedbasisdef)
        basis_function = basis_template;
        basis_function.atom = atom;
        basis_function.A = xyz_a0;
        basis_function.a = [0 0 0];
        basis_function.alpha = reducedbasisdef(M).exponents;
        basis_function.d = reducedbasisdef(M).coeffs;
        basis_function.N = normalization(basis_function.a, basis_function.alpha);
        basis = [basis basis_function];
    end
elseif strcmp(shelltype, 'SP')
    for M = 1:numel(reducedbasisdef)
        % First get S shell basis functions
        basis_function = basis_template; % initialize basis function
        basis_function.a = [0 0 0];
        basis_function.d = reducedbasisdef(M).coeffs(1, :);
        basis_function.atom = atom;
        basis_function.A = xyz_a0;
        basis_function.alpha = reducedbasisdef(M).exponents;
        basis_function.N = normalization(basis_function.a, basis_function.alpha);
        basis = [basis basis_function]; % append
        % Then get P shell basis functions
        for shell_id = 1:3
            basis_function = basis_template; % initialize basis function
            if shell_id == 1
                basis_function.a = [1 0 0];
            elseif shell_id == 2
                basis_function.a = [0 1 0];
            elseif shell_id == 3
                basis_function.a = [0 0 1];
            end
            basis_function.d = reducedbasisdef(M).coeffs(2, :);
            basis_function.atom = atom;
            basis_function.A = xyz_a0;
            basis_function.alpha = reducedbasisdef(M).exponents;
            basis_function.N = normalization(basis_function.a, basis_function.alpha);
            basis = [basis basis_function]; % append
        end
    end
elseif strcmp(shelltype, 'P')
    for M = 1:numel(reducedbasisdef)
        for shell_id = 1:3
            basis_function = basis_template; % initialize basis function
            % add basis functions for both S and P(x,y,z) shells
            if shell_id == 1
                basis_function.a = [1 0 0];
            elseif shell_id == 2
                basis_function.a = [0 1 0];
            elseif shell_id == 3
                basis_function.a = [0 0 1];
            end
            basis_function.atom = atom;
            basis_function.A = xyz_a0;
            basis_function.alpha = reducedbasisdef(M).exponents;
            basis_function.d = reducedbasisdef(M).coeffs;
            basis_function.N = normalization(basis_function.a, basis_function.alpha);
            basis = [basis basis_function]; % append
        end
    end
elseif strcmp(shelltype, 'D')
    for M = 1:numel(reducedbasisdef)
        for shell_id = 1:6
            basis_function = basis_template; % initialize basis function
            % add basis functions for both S and P(x,y,z) shells
            if shell_id == 1
                basis_function.a = [2 0 0];
            elseif shell_id == 2
                basis_function.a = [1 1 0];
            elseif shell_id == 3
                basis_function.a = [1 0 1];
            elseif shell_id == 4
                basis_function.a = [0 2 0];
            elseif shell_id == 5
                basis_function.a = [0 1 1];
            elseif shell_id == 6
                basis_function.a = [0 0 2];
            end
            basis_function.atom = atom;
            basis_function.A = xyz_a0;
            basis_function.alpha = reducedbasisdef(M).exponents;
            basis_function.d = reducedbasisdef(M).coeffs;
            basis_function.N = normalization(basis_function.a, basis_function.alpha);
            basis = [basis basis_function]; % append
        end
    end
end
% end
% basis = cell2mat(basis); % join all basis structures into a single structure
end

function N = normalization(a, alpha)
N = (2/pi)^(3/4) * (2^(a(1)+a(2)+a(3)) * alpha.^((2*(a(1)+a(2)+a(3))+3)/4))...
    /sqrt(fact2(2*a(1)-1) * fact2(2*a(2)-1) * fact2(2*a(3)-1));
end


function reducedbasisdef = getReducedBasisDef(atom, shelltype, basisdef)
reducedbasisdef = {};
for id = 1:numel(basisdef{atom})
    if strcmp(basisdef{atom}(id).shelltype, shelltype)
        reducedbasisdef = [reducedbasisdef basisdef{atom}(id)];
    end
end
reducedbasisdef = cell2mat(reducedbasisdef); % concatenate all the structures into a single structure
end