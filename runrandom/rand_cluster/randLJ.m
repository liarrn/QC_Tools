function coor = randLJ(natoms, rdmin, rdmax)
    % length = (natoms/6.0)^(1/3);
	% rdmin = 2.6;
	% rdmax = 3.4;
    length = (natoms*4.2)^(1/3);
    
    dis = 0.0;
    coor = zeros(natoms, 3);
    coor(1, :) = ball() * length;
	fail = 0;
    for i = 2: natoms
		if fail == 1
			break;
		end
        satisfied = 0;
		iter = 0;
        while (satisfied == 0)
			% disp(iter)
			iter = iter + 1;
			if iter >10000
				% disp('hello')
				fail = 1;
				break;
			end
            satisfied = 1;
            tmp = ball() * length;
            if tooclose(tmp, coor, rdmin, i-1) || toofar(tmp, coor, rdmax, i-1)
				satisfied = 0;
				% disp('hello')
			end
        end
		coor(i, :) = tmp;
    end
	if fail == 1
		% disp('hello')
		coor = -1;
	end
end
