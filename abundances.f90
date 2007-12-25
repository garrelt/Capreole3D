module abundances

  ! Module for Capreole (f90) or C2-Ray
  ! Author: Garrelt Mellema
  ! Date: 2007-12-24 (older but no previous date given)
  ! This module is also accepted by the F compiler (checked Dec 24, 2007)

  ! This module contains abundance related parameter definitions

  use precision, only: dp

  real(kind=dp),parameter,public :: abu_he=0.1_dp
  real(kind=dp),parameter,public :: abu_c=3.3e-4_dp
  real(kind=dp),parameter,public :: abu_o=3.0e-4_dp

  ! Mean molecular weight (assuming metallicities of order 10^-4).
  real(kind=dp),parameter,public :: mu=(1.0_dp-abu_he)+4.0_dp*abu_he

  !-----------------------------------------------------------------------
  !     set heavy element abundances
  !     from Osterbrock (1989), Table 5.13
  !-----------------------------------------------------------------------
  !parameter(abu_he=0.1)
  !abu_he=0.074              !cosmic abundance of He (24.2% by mass)
  !parameter(abu_c=7.1e-7)
  !abu_n=1.3e-7
  !abu_o=4.2e-7
  !abu_ne=1.0e-7
  !abu_s=1.0e-8
  
  !abu_c=3.3e-4
  !   abu_n=0.0
  !   abu_o=6.6e-4
  !   abu_ne=8.3e-5
  !   abu_s=1.6e-5
  
end module abundances
