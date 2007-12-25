module atomic

  ! Module for Capreole (f90)
  ! Author: Garrelt Mellema
  ! Date: 2003-06-01
  ! This module is also accepted by the F compiler (rechecked Dec 24, 2007)

  ! This module contains atomic constants and parameter definitions 

  ! Version: cgs

  use precision, only: dp
  use scaling, only: scvelo
  use cgsconstants, only: kb, m_p

  private

  real(kind=dp),public,parameter :: boltzm = kb/(m_p*scvelo*scvelo)

  ! adiabatic index
  real(kind=dp),public,parameter :: gamma = 5.0_dp/3.0_dp
  real(kind=dp),public,parameter :: gamma1 = gamma - 1.0_dp

end module atomic
