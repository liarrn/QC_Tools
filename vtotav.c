#include<stdio.h>
#include<stdlib.h>

#define IN_FILE "LOCPOT"
#define NRATOM 146
//number of atoms

int main()
{
	FILE *ifp;
	char s_tmp[500];
	int i, j, k;
	int nx, ny, nz;
	float pot;
	float tmp;
	
	
	ifp = fopen(IN_FILE, "r");
	for(i = 0; i < 8 + NRATOM + 1; i ++) fgets(s_tmp, 499, ifp);
	fscanf(ifp, "%d  %d  %d", &nx, &ny, &nz);
	fgets(s_tmp, 499, ifp);
	for(k = 0; k < nz; k ++)
	{
		pot = 0.0;
		for(i = 0; i < nx; i ++)
		{
			for(j = 0; j < ny; j ++)
			{
				fscanf(ifp, "%f", &tmp);
				pot += tmp;
			}
		}
		pot = pot / nx / ny;
		printf("%d %f\n", k + 1, pot);
	}
	return 0;
}