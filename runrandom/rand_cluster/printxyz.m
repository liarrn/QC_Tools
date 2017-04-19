function printxyz(coor, filename, permission, n)
	if nargin <3
		error('too little args');
	else if nargin == 3
		natoms = size(coor);
		natoms = natoms(1);
	else nargin == 4
		natoms = n;
	end
	fp = fopen(filename, permission);
	fprintf(fp, '%d\nAu%d\n', natoms, natoms);
	for i = 1: natoms
		fprintf(fp, '\tAu\t%.3f\t%.3f\t%.3f\n', coor(i, :));
	end
	fclose(fp);
end	