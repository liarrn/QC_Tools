#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<malloc.h>

//control file: para; geometry file: GEOMETRY.xyz; run dmol script: *sh; input file for dmol: input.input
#define PI 3.1415926535

int Natoms, Ncluster, Nrelax;
int MXCY = 1000;
double rdmin, rdmax;
int TEMP = 3000, NSYS = 1;
int PBC_TOGGLE = 0;
double A_AXIS, B_AXIS, C_AXIS, ALPHA_ANG, BETA_ANG, GAMMA_ANG;
int CNA_TOGGLE = 1;
double EnCut = 0.005, CNAAcpR = 0.5, CUTOFF_DIS = 5.5;
int CNACut = 10;

int initial_param();
int initial_str(char elem[500][3], double x[500], double y[500], double z[500]);
int IODMOL(char elem[500][3], double x[500], double y[500], double z[500], double *ener, int *conv);
int move(char elem[500][3], double x[500], double y[500], double z[500]);
int move_clsuter(char elem[500][3], double x[500], double y[500], double z[500]);
int EXTRACT(double x[500], double y[500], double z[500], double *ener, int *conv);
int CNA_CALC(char elem[500][3], double x[500], double y[500], double z[500], int CNATmp[1000]);
//int CNA_ADD(int m, int CNATmp[1000], int CNACount[MXCY][1000]);
int CNA_CHECK(int CNATmp[1000], int **CNACount, int idenList[100], int *idenNum);
int move_cluster(char elem[500][3], double x[500], double y[500], double z[500]);

int main()
{
	int i, j, k, p, q, m, n;
	//double x[500], y[500], z[500];
	double *x, *y, *z;
	//double xf[MXCY][500], yf[MXCY][500], zf[MXCY][500];
	double **xf, **yf, **zf;
	//double energy, en[MXCY];
	double energy, *en;
	double rms, rms1;
	double ratio;
	double tmpd;
	//char **elem;
	char elem[500][3];
	int conv;
	FILE *fp;
	time_t t;

	int **CNACount, *CNATmp;
	//int CNACount[MXCY][1000], CNATmp[1000];
	int cnaAlike;
	int idenNum, idenList[100];

	initial_param();
	x = (double *)malloc(sizeof(double)* 500);
	y = (double *)malloc(sizeof(double)* 500);
	z = (double *)malloc(sizeof(double)* 500);
	en = (double *)malloc(sizeof(double) * MXCY);
	xf = (double **)malloc(sizeof(double *) * MXCY);
	yf = (double **)malloc(sizeof(double *) * MXCY);
	zf = (double **)malloc(sizeof(double *) * MXCY);
	CNATmp = (int *)malloc(sizeof(int)* 1000);
	CNACount = (int **)malloc(sizeof(int *)* MXCY);
	for (i = 0; i < MXCY; i++)
	{
		xf[i] = (double *)malloc(sizeof(double) * 500);
		yf[i] = (double *)malloc(sizeof(double) * 500);
		zf[i] = (double *)malloc(sizeof(double) * 500);
		CNACount[i] = (int *)malloc(sizeof(int)* 1000);
	}
	// elem = (char **)malloc(sizeof(char *)* 500);
	// for (i = 0; i < 500; i++)
		// elem[i] = (char *)malloc(sizeof(char)* 3);
	
	srand((unsigned) time(&t));
	for (i = 0; i < MXCY; i ++)
		for (j = 0; j < 1000; j ++)
			CNACount[i][j] = 0;
	
	for (i = 0; i < MXCY; i ++)
	{
		for (j = 0; j < 500; j ++)
		{
			xf[i][j]=0.0;
			yf[i][j]=0.0;
			zf[i][j]=0.0;
		}
		en[i] = 0.0;
	}
	energy = 0.0;
	
	fp = fopen("END", "w");
	fprintf(fp, "ENDREAD\n");
	fclose(fp);
	
	initial_str(elem, x, y, z);
	IODMOL(elem, x, y, z, &energy, &conv);
	en[0] = energy;
	for (i = 0; i < Natoms; i ++)
	{
		xf[0][i] = x[i];
		yf[0][i] = y[i];
		zf[0][i] = z[i];
	}
	
	if (CNA_TOGGLE == 1)
	{
		CNA_CALC(elem, x, y, z, CNATmp);
		//CNA_ADD(0, CNATmp, CNACount);
		for (i = 0; i < 1000; i++)
			CNACount[0][i] = CNATmp[i];
		//fp = fopen("cna", "w");
		//for (i = 0; i < 1000; i++)
		//{
		//	if (CNACount[0][i] != 0)
		//		fprintf(fp, "%d	%d\n", i, CNACount[0][i]);
		//}
		//fclose(fp);
	}
	
	n = 0;
	
	for (i = 1; i < MXCY; i ++)
	{
		cnaAlike = 0;
		move(elem, x, y, z);
		IODMOL(elem, x, y, z, &energy, &conv);
		
		m = i;
		en[m] = energy;
		for (j = 0; j < Natoms; j ++)
		{
			xf[m][j] = x[j];
			yf[m][j] = y[j];
			zf[m][j] = z[j];
		}
		
		if (CNA_TOGGLE == 1)
		{
			cnaAlike = 0;
			CNA_CALC(elem, x, y, z, CNATmp);
			CNA_CHECK(CNATmp, CNACount, idenList, &idenNum);
			//CNA_ADD(m, CNATmp, CNACount);
			for (j = 0; j < 1000; j++)
				CNACount[m][j] = CNATmp[j];
			if (idenNum == 0)
				cnaAlike = 0;
			else
			{
				for (p = 0; p < idenNum; p ++)
				{
					if (fabs(en[idenList[p]] - energy) < EnCut)
					{
						if (rand() % 1000 / 1000.0 < CNAAcpR)
						{
							cnaAlike = 1;
							break;
						}
					}
				}
			}
		}
		
		ratio = exp(-(en[m] - en[n]) * 2625499.62 / (8.3144621 * TEMP));
		if (en[m] < en[n])
			n = m;
		else if (rand() % 1000 /1000.0 < ratio && cnaAlike == 0)
			n = m;
		else
		{
			for (j = 0; j < Natoms; j ++)
			{
				x[j] = xf[n][j];
				y[j] = yf[n][j];
				z[j] = zf[n][j];
			}
		}
		
		fp = fopen("out.xyz", "w");
		for (j = 0; j <= i; j ++)
		{
			fprintf(fp, "%d\n", Natoms);
			fprintf(fp, "%dth	%lfHa\n", j, en[j]);
			for (k = 0; k < Natoms; k ++)
				fprintf(fp, "%5s%8.3lf%8.3lf%8.3lf\n", elem[k], xf[j][k], yf[j][k], zf[j][k]);
		}
		fclose(fp);
		
		if (CNA_TOGGLE == 1)
		{
			fp = fopen("cna", "w");
			for (j = 0; j <= i; j ++)
			{
				fprintf(fp, "%dth	%lfHa\n", j, en[j]);
				for (k = 1; k < 1000; k ++)
					if (CNACount[j][k] != 0)
						fprintf(fp, "%d    %d\n", k, CNACount[j][k]);
				fprintf(fp, "\n");
			}
			fclose(fp);
		}
	}
	
	for (i = 1; i < MXCY; i ++)
	{
		for (j = i; j > 0; j --)
		{
			if (en[j] < en[j-1])
			{
				tmpd = en[j-1];
				en[j-1] = en[j];
				en[j] = tmpd;
				for (k = 0; k < Natoms; k ++)
				{
					tmpd = xf[j-1][k];
					xf[j-1][k] = xf[j][k];
					xf[j][k] = tmpd;
					tmpd = yf[j-1][k];
					yf[j-1][k] = yf[j][k];
					yf[j][k] = tmpd;
					tmpd = zf[j-1][k];
					zf[j-1][k] = zf[j][k];
					zf[j][k] = tmpd;
				}
			}
			else break;
		}
	}
	
	fp = fopen("final.xyz", "w");
	for (i = 0; i < MXCY; i ++)
	{
		fprintf(fp, "%d\n", Natoms);
		fprintf(fp, "%dth	%lfHa\n", i, en[i]);
		for (j = 0; j < Natoms; j ++)
			fprintf(fp, "%s    %lf    %lf    %lf\n", elem[j], xf[i][j], yf[i][j], zf[i][j]);
	}
	fclose(fp);
	
	free(x); free(y); free(z);
	free(en);
	free(CNATmp);
	for (i = 0; i < MXCY; i++)
	{
		free(xf[i]);
		free(yf[i]);
		free(zf[i]);
		free(CNACount[i]);
	}
	free(xf); free(yf); free(zf);
	free(CNACount);
	// for (i = 0; i < 500; i++)
		// free(elem[i]);
	// free(elem);
	
	return 0;
}


int initial_str(char elem[500][3], double x[500], double y[500], double z[500])
{
	FILE *geoFile;
	int i, j, k;
	char tmpStr[100];
	
	geoFile = fopen("GEOMETRY.xyz", "r");
	fgets(tmpStr, 100, geoFile);
	fgets(tmpStr, 100, geoFile);
	for(i = 0; i < Natoms; i ++)
	{
		fscanf(geoFile, "%s", elem[i]);
		fscanf(geoFile, "%lf", &x[i]);
		fscanf(geoFile, "%lf", &y[i]);
		fscanf(geoFile, "%lf", &z[i]);
	}
	fclose(geoFile);
	return 0;
}

int IODMOL(char elem[500][3], double x[500], double y[500], double z[500], double *ener, int *conv)
{
	int i, j, k;
	FILE *carFile;
	char *blanks = "                     ";
	
	carFile = fopen("input.car", "w");
	fprintf(carFile, "!BIOSYM archive 3\n");
	if(PBC_TOGGLE == 0)
		fprintf(carFile, "PBC=OFF\n");
	else if (PBC_TOGGLE == 1)
	{
		fprintf(carFile, "PBC=ON\n");
		fprintf(carFile, "%f	%f	%f	%f	%f	%f	(P1)\n", A_AXIS, B_AXIS, C_AXIS, ALPHA_ANG, BETA_ANG, GAMMA_ANG);
	}
	fprintf(carFile, "Materials Studio Generated CAR File\n");
	fprintf(carFile, "!DATE Wed Nov 29 13:31:15 2006\n");
	for (i = 0; i < Ncluster; i ++)
		fprintf(carFile, "%-5s%15.9f%15.9f%15.9f%s%-2s%8.3f\n", elem[i], x[i], y[i], z[i], blanks, elem[i], 0.000);
	for (i = Ncluster; i < Ncluster + Nrelax; i ++)
		fprintf(carFile, "%-5s%15.9f%15.9f%15.9f%s%-2s%8.3f\n", elem[i], x[i], y[i], z[i], blanks, elem[i], 0.000);
	for (i = Ncluster + Nrelax; i < Natoms; i ++)
		fprintf(carFile, "%-5s%15.9f%15.9f%15.9f%s%-2s%8.3f\n", elem[i], x[i], y[i], z[i], blanks, elem[i], 0.000);
	fprintf(carFile, "end\n");
	fprintf(carFile, "end\n");
	fclose(carFile);
	system("rm tmpp");
	system("./*.sh -np 12 input");
	system("cat input.outmol END > tmpp");
	system("cp tmpp input.outmol");
	EXTRACT(x, y, z, ener, conv);
	return 0;
}

int move(char elem[500][3], double x[500], double y[500], double z[500])
{
	if (NSYS == 1)
		move_cluster(elem, x, y, z);
	//else if (NSYS == 2)
	//	exchange_cs(elem, x, y, z);
	//else if (NSYS == 3)
	//	move_cluster_on_sub(elem, x, y, z);
	return 0;
}

int move_cluster(char elem[500][3], double x[500], double y[500], double z[500])
{
	int i, j, k;
	int nfar;
	int flag;
	double xcen_i,ycen_i,zcen_i;
	double xcen_f,ycen_f,zcen_f;
	double r_ave_i,r_ave_f, rij;
	double xx, yy, zz;
	double r[500];
	double psi, fi;
	time_t t;
	
	srand((unsigned) time(&t));
	
	xcen_i = 0.0;
	ycen_i = 0.0;
	zcen_i = 0.0;
	r_ave_i = 0.0;
	r_ave_f = 0.0;
	
	for (i = 0; i < Ncluster; i ++)
	{
		xcen_i += x[i];
		ycen_i += y[i];
		zcen_i += z[i];
	}
	xcen_i=xcen_i/Ncluster;
	ycen_i=ycen_i/Ncluster;
	zcen_i=zcen_i/Ncluster;
	xcen_f = xcen_i;
	ycen_f = ycen_i;
	zcen_f = zcen_i;
	
	for (i = 0; i < Ncluster; i ++)
	{
		r[i] = sqrt(pow(x[i] - xcen_i, 2) + pow(y[i] - ycen_i, 2) + pow(z[i] - zcen_i, 2));
		r_ave_i += r[i];
	}
	r_ave_i /= Ncluster;
	r_ave_f = r_ave_i;
	


	while (fabs(r_ave_f - r_ave_i) <= 0.1)
	{
		r_ave_i = r_ave_f;
		nfar = 0;
		for (i = 0; i < Ncluster; i ++)
			if(r[i] > r[nfar])
				nfar = i;
		
		flag = 1;
		while (flag == 1)
		{
			int toofar;
			
			toofar = 1;
			flag = 0;
			fi = 2 * PI * (rand() % 1000 / 1000.0);
			psi = 2 * PI * (rand() % 1000 / 1000.0);
			
			xx = xcen_f+r[nfar]*(1+pow(-1, (int)(0.5+rand()%1000/1000.0))*0.5*(rand()%1000/1000.0))*cos(fi)*cos(psi);
			yy = ycen_f+r[nfar]*(1+pow(-1, (int)(0.5+rand()%1000/1000.0))*0.5*(rand()%1000/1000.0))*cos(fi)*sin(psi);
			zz = zcen_f+r[nfar]*(1+pow(-1, (int)(0.5+rand()%1000/1000.0))*0.5*(rand()%1000/1000.0))*sin(fi);

			//xx = xcen_f + r[nfar] * (0.5 + rand() % 1000 / 1000.0)*cos(fi)*cos(psi);
			//yy = ycen_f + r[nfar] * (0.5 + rand() % 1000 / 1000.0)*cos(fi)*sin(psi);
			//zz = zcen_f + r[nfar] * (0.5 + rand() % 1000 / 1000.0)*sin(fi);
			
			for (i = 0; i < Ncluster; i ++)
			{
				if (i == nfar) continue;
				rij = sqrt(pow(x[i] - xx, 2) + pow(y[i] - yy, 2) + pow(z[i] - zz, 2));
				if (rij < rdmin)
				{
					flag = 1;
					break;
				}
				if (rij < rdmax) toofar = 0;
			}
			if (toofar == 1)
			{
				flag = 1;
				continue;
			}
		}
		
		x[nfar] = xx;
		y[nfar] = yy;
		z[nfar] = zz;
		xcen_f = 0.0;
		ycen_f = 0.0;
		zcen_f = 0.0;
		
		for (i = 0; i < Ncluster; i ++)
		{
			xcen_f += x[i];
			ycen_f += y[i];
			zcen_f += z[i];
		}
		xcen_f = xcen_f/Ncluster;
		ycen_f = ycen_f/Ncluster;
		zcen_f = zcen_f/Ncluster;
		
		r_ave_f = 0.0;
		for (i = 0; i < Ncluster; i ++)
		{
			r[i] = sqrt(pow(x[i] - xcen_f, 2) + pow(y[i] - ycen_f, 2) + pow(z[i] - zcen_f, 2));
			r_ave_f += r[i];
		}
		r_ave_f /= Ncluster;

	}
	

	
	return 0;
}

int EXTRACT(double x[500], double y[500], double z[500], double *ener, int *conv)
{
	double xf[500], yf[500], zf[500];
	int i, j, k;
	char *k1, *k2, *k3, *k4;
	FILE *outFile;
	char line[500];
	char tmpstr[100];
	char tmpd;
	
	*conv = 0;
	outFile = fopen("input.outmol", "r");
	if (outFile == NULL)
		exit (EXIT_FAILURE);
	while (fgets(line, 500, outFile))
	{
		k1 = strstr(line, "Ef");
		k2 = strstr(line, "Final Coordinates (Angstroms)");
		k3 = strstr(line, "Message: DMol3 job finished successfully");
		k4 = strstr(line, "ENDREAD");
		if (k1 != NULL)
		{
			for(i = 0; i < 13; i ++)
				tmpstr[i] = line[i + 9];
			tmpstr[i] = '\0';
			sscanf(tmpstr, "%lf", ener);
		}
		if (k2 != NULL)
		{
			fgets(line, 500, outFile);
			fgets(line, 500, outFile);
			for (i = 0; i < Natoms; i ++)
			{
				fgets(line, 500, outFile);
				k = 1;
				sscanf(line, "%d%s%lf%lf%lf", &tmpd, &tmpstr, &xf[i], &yf[i], &zf[i]);
			}
			k = 1;
		}
		if (k3 != NULL)
		{
			*conv = 1;
			for (i = 0; i < Natoms; i ++)
			{
				x[i] = xf[i];
				y[i] = yf[i];
				z[i] = zf[i];
			}
			break;
		}
		if ((strstr(line, "Error") != NULL) || (k4 != NULL))
		{
			*conv = 0;
			break;
		}
	}
	fclose(outFile);
	return 0;
}

int CNA_CALC(char elem[500][3], double x[500], double y[500], double z[500], int CNATmp[1000])
{
	int i, j, k;
	int nij, ix1, ix2, ix3;
	int set[1000], Aset[1000][2], setm, setn;
	double rij, rxmn, rymn, rzmn, rmn, rik, rjk;
	int maxlength, m, n, p, q, a, b, flag, c;
	int SubBond[1000], nBond[1000];
	int CNAType;
	
	for (i = 0; i < 1000; i ++)
		CNATmp[i] = 0;
	
	nij = 0;
	for (i = 0; i < Natoms - 1; i ++)
	{
		for (j = i + 1; j < Natoms; j ++)
		{
			rij = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2) + pow(z[i] - z[j], 2));
			if (rij > CUTOFF_DIS)
				continue;
			nij ++;
			ix1 = ix2 = ix3 = 0;
			
			// the first index 
			for (k = 0; k < Natoms; k ++)
			{
				if (k == i || k == j) continue;
				rik = sqrt(pow(x[i] - x[k], 2) + pow(y[i] - y[k], 2) + pow(z[i] - z[k], 2));
				rjk = sqrt(pow(x[j] - x[k], 2) + pow(y[j] - y[k], 2) + pow(z[j] - z[k], 2));
				if (rik < CUTOFF_DIS && rjk < CUTOFF_DIS)
				{
					set[ix1] = k;
					ix1 ++;
				}
			}
			
			// the second index
			for (m = 0; m < ix1 -1; m ++)
			{
				for(n = m + 1; n < ix1; n ++)
				{
					setm = set[m];
					setn = set[n];
					rmn = sqrt(pow(x[setm] - x[setn], 2) + pow(y[setm] - y[setn], 2) + pow(z[setm] - z[setn], 2));
					if (rmn < CUTOFF_DIS)
					{
						Aset[ix2][0] = setm;
						Aset[ix2][1] = setn;
						ix2 ++;
					}
				}
			}
			
			// the third index
			maxlength = 0;
			for (p = 0; p < Natoms; p ++)
			{
				SubBond[p] = p;
				nBond[p] = 0;
			}
			for (p = 0; p < ix2; p ++)
			{
				m = Aset[p][0];
				n = Aset[p][1];
				a = nBond[m];
				b= nBond[n];
				c = SubBond[m];
				if (SubBond[m] == SubBond[n])
				{
					for (q = 0; q < Natoms; q++)
					if (SubBond[q] == c)
						nBond[q] = nBond[q] + 1;
				}
				else
				{
					for (q = 0; q < Natoms; q ++)
					{
						if (SubBond[q] == c)
							SubBond[q] = SubBond[n];
						if(SubBond[q] == SubBond[n])
							nBond[q] = a + b + 1;
					}
				}
			}
			for (q = 0; q < Natoms; q ++)
				if (maxlength < nBond[q])
					maxlength = nBond[q];
			ix3= maxlength;
			
			CNAType = ix1 * 100 + ix2 * 10 + ix3;
			if (CNAType >= 1000)
				CNAType = 999;
			CNATmp[CNAType] += 1;
		}
	}
	return 0;
}

//int CNA_ADD(int m, int CNATmp[1000], int CNACount[MXCY][1000])
//{
//	int i;
//	for (i = 0; i < 1000; i ++)
//		CNACount[m][i] = CNATmp[i];
//	return 0;
//}

int CNA_CHECK(int CNATmp[1000], int **CNACount, int idenList[100], int *idenNum)
{
	int CNADiff;
	int i, j;
	
	*idenNum = 0;
	for (i = 0; i < 100; i ++)
		idenList[i] = 0;
	for (i = 0; i < MXCY; i ++)
	{
		CNADiff = 0;
		for (j = 0; j < 1000; j ++)
			CNADiff += abs(CNACount[i][j] - CNATmp[j]);
		if(CNADiff < CNACut)
		{
			idenList[*idenNum] = i;
			*idenNum += 1;
		}
	}
	return 0;
}

int initial_param()
{
	FILE *fp;
	char line[500];
	char tmpStr[500];
	double tmpd;
	char paraStr[20][500] = {"Natoms", "Ncluster", "Nrelax", "MXCY",\
	"rdmin", "rdmax", "TEMP", "NSYS", "PBC_TOGGLE", "A_AXIS", "B_AXIS", "C_AXIS", \
	"ALPHA_ANG", "BETA_ANG", "GAMMA_ANG", "CNA_TOGGLE", "EnCut", "CNACut", \
	"CNAAcpR", "CUTOFF_DIS"};
	double paraVal[20];
	int i;
	
	fp = fopen("para", "r");
	while (fgets(line, 500, fp))
	{
		if (strcmp(line, "\n") == 0) continue;
		sscanf(line, "%s%lf", tmpStr, &tmpd);
		for (i = 0; i < 20; i ++)
		{
			if (strcmp(tmpStr, paraStr[i]) == 0)
			{
				paraVal[i] = tmpd;
				break;
			}
		}
	}
	Natoms = (int)paraVal[0]; Ncluster = (int)paraVal[1]; Nrelax = (int)paraVal[2];
	MXCY = (int)paraVal[3];
	rdmin = paraVal[4]; rdmax = paraVal[5];
	TEMP = (int)paraVal[6];
	NSYS = (int)paraVal[7];
	PBC_TOGGLE = (int)paraVal[8];
	A_AXIS = paraVal[9]; B_AXIS = paraVal[10]; C_AXIS = paraVal[11];
	ALPHA_ANG = paraVal[12]; BETA_ANG = paraVal[13]; GAMMA_ANG = paraVal[14];
	CNA_TOGGLE = (int)paraVal[15];
	EnCut = paraVal[16]; CNACut = (int)paraVal[17]; CNAAcpR = paraVal[18]; CUTOFF_DIS = paraVal[19];
	fclose(fp);
	return 0;
}




