function coor = ball(n)
	if nargin == 0
		n = 1;
	end
	coor = zeros(n, 3);
	for i = 1: n
		while (1)
			tmpCoor = rand(1, 3)* 2 - 1.0;
			dist = sqrt(tmpCoor*tmpCoor');
			if dist<1.0
				break
			end
		end
		coor(i, :) = tmpCoor;
	end	
end