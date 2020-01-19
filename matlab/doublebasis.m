function unrestrictedbasis = doublebasis(basis)
M = numel(basis);
unrestrictedbasis = [basis(1) basis(1)];
for iBasis = 2:M
    unrestrictedbasis = [unrestrictedbasis basis(iBasis) basis(iBasis)];
end
for iBasis = 1:2*M
    if mod(iBasis, 2) == 0
        unrestrictedbasis(iBasis).spin = 1/2;
    else
        unrestrictedbasis(iBasis).spin = -1/2;
    end
end
end