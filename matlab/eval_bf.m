function val = eval_bf(basisfun, xyz_a0)
% val = eval_bf(basisfun,xyz_a0)
%
% Input:
%   basisfun    structure with basis function information, from list of basis functions
%   xyz_a0      mx3 list of Cartesian coordinates, one point per row, over which to
%               evaluate the basis function
% Output:
%   val         m-element array, with the values of the basis function 

shape = size(xyz_a0);
m = shape(1);
val = zeros(m, 1);

for i = 1:m
    A = basisfun.A;
    a = basisfun.a;
    for j = 1:numel(basisfun.d)
        d = basisfun.d(j);
        N = basisfun.N(j);
        alpha = basisfun.alpha(j);
        val(i) = val(i) + d*N*exp(-alpha*norm(xyz_a0(i,:)-A)^2);
    end
    if m > 1
        val(i) = (xyz_a0(i,1)-A(1))^a(1) * (xyz_a0(i,2)-A(2))^a(2) *...
            (xyz_a0(i,3)-A(3))^a(3) * val(i);
    elseif m == 1
        val(i) = (xyz_a0(1)-A(1))^a(1) * (xyz_a0(2)-A(2))^a(2) *...
            (xyz_a0(3)-A(3))^a(3) * val(i);
end
end