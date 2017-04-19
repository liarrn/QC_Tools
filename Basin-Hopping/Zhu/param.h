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
!  This program is written by Yi Gao onh 2014/8/9.
          PARAMETER(Natoms=7,Ncluster=7,Nrelax=0)
          PARAMETER(MXCY=1000,MXLN=10000000,rdmin=2.0,rdmax=4.0)
          PARAMETER(TEMP=3000,pi=3.1415926535,NSYS=1)

