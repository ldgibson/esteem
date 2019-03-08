function moplot(atoms, xyz_a0, out, iMO, level)

dim = max(abs(xyz_a0), [], 'all') + 2;

[x, y, z] = meshgrid(-dim:0.1:dim);
x0 = reshape(x, 1, numel(x));
y0 = reshape(y, 1, numel(y));
z0 = reshape(z, 1, numel(z));

xyz = zeros(numel(x0), 3);
xyz(:, 1) = x0;
xyz(:, 2) = y0;
xyz(:, 3) = z0;

MO = zeros(size(x0));
for m = 1:numel(out.basis)
    MO(:) = MO(:) + out.C(m, iMO) * eval_bf(out.basis(m), xyz);
end

maxval = max(MO);
MO = reshape(MO, size(x));

isosurface(x, y, z, MO, level*maxval, 'r')
isosurface(x, y, z, MO, -level*maxval, 'blue')
atomscale = 8;  % tuning parameter for atom size in plot
hold on
for at = 1:numel(atoms)
    [a, b, c] = sphere;
    surf(xyz_a0(at, 1)-a/atomscale, ...
        xyz_a0(at, 2)-b/atomscale, ...
        xyz_a0(at, 3)-c/atomscale);
end
axis equal;

end