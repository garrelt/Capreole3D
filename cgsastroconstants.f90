module astroconstants

  use precision

  ! A collection of astronomical units and conversion factors
  ! Units: cgs

  real(kind=dp),parameter :: R_SOLAR=6.96e10
  real(kind=dp),parameter :: L_SOLAR=3.86e33

  real(kind=dp),parameter :: YEAR=3.1567e7

  real(kind=dp),parameter :: pc=3.086e18
  real(kind=dp),parameter :: kpc=1e3*pc
  real(kind=dp),parameter :: Mpc=1e6*pc

end module astroconstants
