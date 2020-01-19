%atoms = [6 8 8 6 6 1 1 1 1 8];
% atoms = [6 8 8 8];
% atoms = [8 1 1];
atoms = [6 8];
% xyz = [0 0 0];
xyz = testdata(2).xyz_a0;
% xyz = [0 0 0; 0 0 2.2827892];
% xyz = [0 0 0];

settings.tolDensity = 1e-8;
settings.tolEnergy = 1e-8;
settings.basisset = '6-31G';
settings.method = 'RHF';
settings.multiplicity = 1;
% settings.ExchFunctional = 'Slater';
% settings.CorrFunctional = 'VWN5';
% settings.nRadialPoints = 100;
% settings.nAngularPoints = 302;

charge = 0;

out = mocalc(atoms, xyz, charge, settings)