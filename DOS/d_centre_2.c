#include<stdio.h>
#include<stdlib.h>

#define IN_FILE "DOSCAR"
#define TDOS_FILE "TDOS.dat"
#define PDOS_FILE "PDOS.dat"

int main()
//used only for spin-polarized calculation
{
	char s_tmp[300];
	int i, j, k;
	int nrbands;
	int p_start, p_finish; //define the range of atoms for PDOS print
	float efermi;
	float f_tmp;
	float energy, spin_up, spin_down, total;
	float p_spinup, p_spindown, d_spinup, d_spindown, t_spinup, t_spindown;
	float d_centre_down, d_centre_up, d_dos_down, d_dos_up, d_centre, d_dos;
	struct PDOS_LINE
	{
		float energy;
		float s_spinup;
		float s_spindown;
		float py_spinup;
		float py_spindown;
		float pz_spinup;
		float pz_spindown;
		float px_spinup;
		float px_spindown;
		float dxy_spinup;
		float dxy_spindown;
		float dyz_spinup;
		float dyz_spindown;
		float dz2r2_spinup;
		float dz2r2_spindown;
		float dxz_spinup;
		float dxz_spindown;
		float dx2y2_spinup;
		float dx2y2_spindown;
	};
	struct PDOS_LINE *pdos;
	FILE *ifp, *ofpt, *ofpp;
	
	printf("input the atom number range you want calculate PDOS\n two atom numbers should be separated by a single space\n");
	scanf("%d %d", &p_start, &p_finish);
	
	ifp = fopen(IN_FILE, "r");
	ofpt = fopen(TDOS_FILE, "w");
	ofpp = fopen(PDOS_FILE, "w");
	fgets(s_tmp, 300, ifp);
	fgets(s_tmp, 300, ifp);
	fgets(s_tmp, 300, ifp);
	fgets(s_tmp, 300, ifp);
	fgets(s_tmp, 300, ifp);
	fscanf(ifp, "%f", &f_tmp);
	fscanf(ifp, "%f", &f_tmp);
	fscanf(ifp, "%d", &nrbands);
	fscanf(ifp, "%f", &efermi);
	fscanf(ifp, "%f", &f_tmp);
	fgets(s_tmp, 300, ifp);
	
	//read and print TDOS
	fprintf(ofpt, "energy  spin_up  spin_down  total\n");
	for(i = 0;i < nrbands; i ++)
	{
		fscanf(ifp, "%f", &energy);
		fscanf(ifp, "%f", &spin_up);
		fscanf(ifp, "%f", &spin_down);
		fgets(s_tmp, 300, ifp);
		total = spin_up + spin_down;
		fprintf(ofpt, "%f  %f  %f  %f\n", energy - efermi, spin_up, spin_down, total);
	}
	
	//read and print PDOS
	pdos = (struct PDOS_LINE *)malloc(sizeof(struct PDOS_LINE) * nrbands);
	
	for(i = 0; i < p_start - 1; i ++)
	{
		for(j = 0; j <= nrbands; j ++)
			fgets(s_tmp, 300, ifp);
	}
	//read p_start atom's pdos data
	fgets(s_tmp, 300, ifp);
	for(i = 0; i < nrbands; i++)
	{
		fscanf(ifp, "%f", &pdos[i].energy);
		fscanf(ifp, "%f", &pdos[i].s_spinup);
		fscanf(ifp, "%f", &pdos[i].s_spindown);
		fscanf(ifp, "%f", &pdos[i].py_spinup);
		fscanf(ifp, "%f", &pdos[i].py_spindown);
		fscanf(ifp, "%f", &pdos[i].pz_spinup);
		fscanf(ifp, "%f", &pdos[i].pz_spindown);
		fscanf(ifp, "%f", &pdos[i].px_spinup);
		fscanf(ifp, "%f", &pdos[i].px_spindown);
		fscanf(ifp, "%f", &pdos[i].dxy_spinup);
		fscanf(ifp, "%f", &pdos[i].dxy_spindown);
		fscanf(ifp, "%f", &pdos[i].dyz_spinup);
		fscanf(ifp, "%f", &pdos[i].dyz_spindown);
		fscanf(ifp, "%f", &pdos[i].dz2r2_spinup);
		fscanf(ifp, "%f", &pdos[i].dz2r2_spindown);
		fscanf(ifp, "%f", &pdos[i].dxz_spinup);
		fscanf(ifp, "%f", &pdos[i].dxz_spindown);
		fscanf(ifp, "%f", &pdos[i].dx2y2_spinup);
		fscanf(ifp, "%f", &pdos[i].dx2y2_spindown);
		fgets(s_tmp, 300, ifp);
	}
	
	//read p_start +1 to p_finish atoms' pdos data
	for(i = 0; i < (p_finish - p_start); i ++)
	{
		fgets(s_tmp, 300, ifp);
		for(j = 0; j < nrbands; j ++)
		{
			fscanf(ifp, "%f", &f_tmp);
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].s_spinup += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].s_spindown += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].py_spinup += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].py_spindown += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].pz_spinup += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].pz_spindown += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].px_spinup += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].px_spindown += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dxy_spinup += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dxy_spindown += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dyz_spinup += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dyz_spindown += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dz2r2_spinup += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dz2r2_spindown += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dxz_spinup += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dxz_spindown += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dx2y2_spinup += f_tmp;
			fscanf(ifp, "%f", &f_tmp);
			pdos[j].dx2y2_spindown += f_tmp;
			fgets(s_tmp, 300, ifp);
		}
	}
	d_centre_down = d_centre_up = d_dos_down = d_dos_up = 0.0;
	fprintf(ofpp, "energy  s_spinup  s_spindown  p_spinup  p_spindown  d_spinup  d_spindown  t_spinup  t_spindown s p d total\n");
	for(i = 0; i < nrbands; i++)
	{
		p_spinup = pdos[i].py_spinup + pdos[i].pz_spinup + pdos[i].px_spinup;
		p_spindown = pdos[i].py_spindown + pdos[i].pz_spindown + pdos[i].px_spindown;
		d_spinup = pdos[i].dxy_spinup + pdos[i].dyz_spinup + pdos[i].dz2r2_spinup + pdos[i].dxz_spinup + pdos[i].dx2y2_spinup;
		d_spindown = pdos[i].dxy_spindown + pdos[i].dyz_spindown + pdos[i].dz2r2_spindown + pdos[i].dxz_spindown + pdos[i].dx2y2_spindown;
		t_spinup = pdos[i].s_spinup + p_spinup + d_spinup;
		t_spindown = pdos[i].s_spindown + p_spindown + d_spindown;
		if(i != 0)
		// i == 0 produces wrong dos
		{
			if(pdos[i].energy < efermi)
			{
				d_centre_down = d_centre_down + (d_spinup + d_spindown) * pdos[i].energy;
				d_dos_down = d_dos_down + d_spinup + d_spindown;
			}
			else
			{
				d_centre_up = d_centre_up + (d_spinup + d_spindown) * pdos[i].energy;
				d_dos_up = d_dos_up + d_spinup + d_spindown;
			}
		}
		fprintf(ofpp, "%f  %f  %f  %f  %f  %f  %f  %f  %f %f %f %f %f\n", pdos[i].energy - efermi, pdos[i].s_spinup, pdos[i].s_spindown, p_spinup, p_spindown, d_spinup, d_spindown, t_spinup, t_spindown, pdos[i].s_spinup + pdos[i].s_spindown, p_spinup +  p_spindown, d_spinup + d_spindown, t_spinup + t_spindown);
	}
	d_centre = d_dos = 0.0;
	d_centre = d_centre_up + d_centre_down;
	d_dos = d_dos_up + d_dos_down;
	d_centre /= d_dos;
	d_centre -= efermi;
	d_centre_up /= d_dos_up;
	d_centre_up -= efermi;
	d_centre_down /= d_dos_down;
	d_centre_down -= efermi;
	printf("d_centre with fermi level correction are:\n below fermi level: %f\n upon fermi level: %f\n overall: %f\n", d_centre_down, d_centre_up, d_centre);
	
	free(pdos);
	fclose(ifp);
	fclose(ofpt);
	fclose(ofpp);
	return 0;
}