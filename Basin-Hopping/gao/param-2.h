! Natoms is the number of total atoms.
! Ncluster is the number of the cluster.
! Nrelax is the number of the relax atoms of the substrate.
! MXCY is the maximum cycle of the basin-hopping MC.
! MXLN is the maximum lines of the DMol3 output.
! rdmin is the minimum bond length of the cluster.
! rdmax is the maximum bond length of the cluster.
! TEMP is the MC temperature.
! NSYS is the option of the MC system:
!   NSYS=1:  pure cluster
!   NSYS=2:  atom exchange between cluster and substrate
!   NSYS=3:  cluster+substrate (MC on cluster only)
!	PBC_TOGGLE=1 enabling the periodic boundary condition
!	CNA_TOGGLE=1 enabling the CNA checking mechanism
!	Encut is the accepable energy criterion
!	CNACut is the accepable CNA difference
!	CNAAcpR is the CNA checking accept ratio
!	CUTOFF_DIS is the cutoff distance for bond
!  This program is written by Yi Gao onh 2014/8/9.
          PARAMETER(Natoms=7,Ncluster=7,Nrelax=0)
          PARAMETER(MXCY=1000,MXLN=10000000,rdmin=2.0,rdmax=3.5)
          PARAMETER(TEMP=3000,pi=3.1415926535,NSYS=1)
          PARAMETER(PBC_TOGGLE=0)
          PARAMETER(A_AXIS = 10.00, B_AXIS = 10.00, C_AXIS = 10.00)
          PARAMETER(ALPHA_ANG = 90.00, BETA_ANG = 90, GAMMA_ANG = 90)
          PARAMETER(CNA_TOGGLE = 0)
          PARAMETER(EnCut = 0.5, CNACut = 10, CNAAcpR = 0.5)
		  PARAMETER(CUTOFF_DIS = 2.7)


