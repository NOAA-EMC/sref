Module module_data_gocart_dust
  INTEGER, PARAMETER :: ndust=5,ndcls=3,ndsrc=1,maxstypes=100
  INTEGER, PARAMETER :: ngsalt=9
  real, dimension (maxstypes) :: drypoint
  real, dimension (maxstypes) :: porosity

  REAL :: ch_dust(ndust,12)
  REAL, PARAMETER :: dyn_visc = 1.5E-5
  real*8, DIMENSION (5), PARAMETER :: den_dust(5)=(/2500.,2650.,2650.,2650.,2650./)
  real*8, DIMENSION (5), PARAMETER :: reff_dust(5)=(/0.73D-6,1.4D-6,2.4D-6,4.5D-6,8.0D-6/)
  INTEGER, DIMENSION (5), PARAMETER :: ipoint(5)=(/3,2,2,2,2/)
  REAL, DIMENSION (5), PARAMETER :: frac_s(5)=(/0.1,0.25,0.25,0.25,0.25/)

  real*8, DIMENSION (ngsalt), PARAMETER :: reff_salt=(/0.71D-6,1.37D-6,2.63D-6,5.00D-6,9.50D-6,18.1D-6,34.5D-6,65.5D-6,125.D-6/)
  real*8, DIMENSION (ngsalt), PARAMETER :: den_salt=(/2500.,2650.,2650.,2650.,2650.,2650.,2650.,2650.,2650./)
  INTEGER, DIMENSION (ngsalt), PARAMETER :: spoint=(/1,2,2,2,2,2,3,3,3/)
  real*8, DIMENSION (ngsalt), PARAMETER :: frac_salt=(/1.,0.2,0.2,0.2,0.2,0.2,0.333,0.333,0.333/)

  real*8, DIMENSION (ndust), PARAMETER :: distr_dust=(/1.074D-1,1.012D-1,2.078D-1,4.817D-1,1.019D-1/)



END Module module_data_gocart_dust
