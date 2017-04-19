// Natoms is the number of total atoms.
// Ncluster is the number of the cluster.
// Nrelax is the number of the relax atoms of the substrate.
// MXCY is the maximum cycle of the basin-hopping MC.
// MXLN is the maximum lines of the DMol3 output.
// rdmin is the minimum bond length of the cluster.
// rdmax is the maximum bond length of the cluster.
// TEMP is the MC temperature.
// NSYS is the option of the MC system:
//   NSYS=1:  pure cluster
//   NSYS=2:  atom exchange between cluster and substrate
//   NSYS=3:  cluster+substrate (MC on cluster only)
//	PBC_TOGGLE=1 enabling the periodic boundary condition
//	CNA_TOGGLE=1 enabling the CNA checking mechanism
//	Encut is the accepable energy criterion
//	CNACut is the accepable CNA difference
//	CNAAcpR is the CNA checking accept ratio
//	CUTOFF_DIS is the cutoff distance for bond
//  This program is written by Yi Gao onh 2014/8/9.

#define Natoms 19
#define Ncluster 19
#define Nrelax 0

#define MXCY 1000
#define MXLN 10000000000

#define rdmin 2.0
#define rdmax 4.0

#define TEMP 3000
#define NSYS 1
#define PI 3.1415926535

#define PBC_TOGGLE 0
#define A_AXIS 10.00
#define B_AXIS 10.00
#define C_AXIS 10.00

#define ALPHA_ANG 90.00
#define  BETA_ANG 90.00
#define GAMMA_ANG 90.00

#define CNA_TOGGLE 1
#define EnCut 0.5
#define CNACut 10
#define CNAAcpR 0.5
#define CUTOFF_DIS 2.7
