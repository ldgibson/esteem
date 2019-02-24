function S = int_overlap(basis)
% S = int_overlap(basis)
%
% Input:
%     basis     basis information, as obtained by buildbasis
% Output:
%     S         MxM matrix of overlap integrals

M = numel(basis);

S = zeros(M);

for i = 1:M
    for j = 1:M
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

                S(i,j) = S(i,j) + d_k*d_l*N_k*N_l* ...
                    overlap_primitive(a, b, alpha, beta, A, B);
            end
        end
    end
end
end