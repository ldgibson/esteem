function [Vxc,Exc,rhoInt] = int_xc(basis,P,grid,ExchFunctional,CorrFunctional)
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
basisval = zeros(M, length(grid.xyz));
p = zeros(length(grid.xyz), 1);
rhoInt = 0;
Vxc = zeros(M, M);
for i = 1:M
    basisval(i,:) = eval_bf(basis(i), grid.xyz);
end
Exc = 0;
basisprod = zeros(length(grid.xyz), M, M);
for iGrid = 1:length(grid.xyz)
    for r = 1:M
        for s = 1:M
            basisprod(iGrid, r, s) = basisval(r, iGrid)*basisval(s, iGrid);
            p(iGrid) = p(iGrid) + basisprod(iGrid, r, s)*P(r, s);
        end
    end
    if ExchFunctional == 'Slater'
        [ex, Vx] = Slater(p(iGrid));
    end
    if CorrFunctional == 'VWN3'
        [ec, Vc] = VWN(p(iGrid), 3);
    elseif CorrFunctional == 'VWN5'
        [ec, Vc] = VWN(p(iGrid), 5);
    end
    Vxc_ = Vx + Vc;
    for r = 1:M
        for s = 1:M
            iVxc = grid.weights(iGrid)*basisprod(iGrid, r, s)*Vxc_;
            Vxc(r, s) = Vxc(r, s) + iVxc;
        end
    end
    exc = ex + ec;
    Exc = Exc + grid.weights(iGrid) * exc * p(iGrid);
    rhoInt = rhoInt + grid.weights(iGrid)*p(iGrid);
end