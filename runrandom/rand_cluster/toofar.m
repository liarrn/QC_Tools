function result = toofar(atom, cluster, criteria, natoms)
	% return 1 if atom is too far away from cluster, return 0 if new atom is not too far away
	result = 1;
	for i = 1: natoms
		dist = sqrt(sum((cluster(i, :) - atom) .^ 2));
		if dist <= criteria
			result = 0;
			break;
		end
	end
end