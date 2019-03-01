function [Vxc,Exc,p] = int_xc(basis,P,grid,ExchFunctional,CorrFunctional)
% [Vxc,Exc,rhoInt] = int_xc(basis,P,grid,ExchFunctional,CorrFunctional)
% 
% Input:
%   basis           list of M basis functions
%   P               MxM density matrix
%   grid            molecular integration grid, as obtained by molecular_grid
%   ExchFunctional  name of exchange functional (possible values: 'Slater')
%   CorrFunctional  name of correlation functional (possible values: 'VWN3', 'VWN5')
% Output:
%   Vxc             exchange-correlation energy matrix, MxM, in hartrees
%   Exc             exchange-correlation energy, in hartrees
%   rhoInt          integral of the density over the grid, should be equal to
%                   the number of electrons
M = numel(basis);
% basisprod = zeros(length(grid.xyz), 1);
p = zeros(length(grid.xyz), 1);
rhoInt = 0;
Vxc = zeros(M, M);
Exc = 0;
for i = 1:length(grid.xyz)
    for r = 1:M
        for s = 1:M
%             basisprod(i) = basisprod(i) + eval_bf(basis(r), grid.xyz(i,:))*eval_bf(basis(s),grid.xyz(i,:));
            p(i) = p(i) + eval_bf(basis(r), grid.xyz(i,:))*eval_bf(basis(s),grid.xyz(i,:)) * P(r, s);
        end
    end
%     basisprod(i)
    rhoInt = rhoInt + p(i)
    [ex, Vx] = Slater(p(i));
    [ec, Vc] = VWN(p(i), 5);
    Exc = Exc + (ex + ec) * p(i) * grid.weights(i);
    for mu = 1:numel(basis)
        for nu = 1:numel(basis) 
            Vxc(mu, nu) = Vxc(mu, nu) + eval_bf(basis(mu), grid.xyz(i,:)) * (Vx+Vc) ...
                * eval_bf(basis(nu), grid.xyz(i,:));
        end
    end
end