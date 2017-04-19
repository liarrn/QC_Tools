% set number of iterations
num_iter = 10000;
% set number of atoms
num_atoms = 3;
% set filename
filename = 'initiate.xyz';

for i = 1: num_iter
	% disp(i);
	coor = randLJ(num_atoms, 2.6, 3.4);
	if coor ~= -1
		printxyz(coor, filename, 'a');
	end
end