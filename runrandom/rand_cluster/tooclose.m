function result = tooclose(atom, cluster, criteria, natoms)
	% return 1 if atom is too close to cluster, return 0 if new atom is not too close
	result = 0;
	for i = 1: natoms
		dist = sqrt(sum((cluster(i, :) - atom) .^ 2));
		if dist <= criteria
			result = 1;
			break;
		end
	end
end