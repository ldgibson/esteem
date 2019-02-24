function Vne = int_attraction(atoms, xyz_a0, basis)
% Vne = int_attraction(atoms, xyz_a0, basis)
% 
% Input:
%   atoms   list of element numbers (array with K elements); e.g. [6 8] for CO
%   xyz_a0  K×3 array of Cartesian coordinates of nuclei, in bohr
%   basis   basis information, as obtained by buildbasis
% Output:
%   Vne     M×M matrix of attraction energy integrals, in hartrees

M = numel(basis);
Vne = zeros(M);

for i = 1:M
    for j = i:M
        V_ = 0;
        for C_i = 1:numel(atoms)
            % Loop over primitives for each i,j pair
            for k = 1:numel(basis(i).d)
                for l = 1:numel(basis(j).d)
                    % 'k'th primitive basis function in basis function i
                    d_k = basis(i).d(k);
                    N_k = basis(i).N(k);
                    alpha = basis(i).alpha(k);
                    A = basis(i).A;
                    a = basis(i).a;
                    % 'l'th primitive basis function in basis function j
                    d_l = basis(j).d(l);
                    N_l = basis(j).N(l);
                    beta = basis(j).alpha(l);
                    B = basis(j).A;
                    b = basis(j).a;
                    
                    C = xyz_a0(C_i, :);
                    
                    V_ = V_ - atoms(C_i)*d_k*d_l*N_k*N_l*...
                        VRR(0, a, b, alpha, beta, A, B, C);
                end
            end
        end
        Vne(i, j) = V_;
        if i ~= j
            Vne(j, i) = V_;
        end
    end
end

end