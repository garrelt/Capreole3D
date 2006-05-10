module cgsconstants

  use precision

  ! A collection of physical constants and conversion factors
  ! Units: cgs

  ! proton mass
  real(kind=8), parameter :: m_p=1.672661d-24
  ! speed of light
  real(kind=8), parameter :: c=2.997925d+10
  ! Planck constant
  real(kind=8), parameter :: hplanck=6.6260755d-27
      ! Stefan-Boltzmann constant
  real(kind=8), parameter :: sigmasb=5.670d-5
      ! Boltzmann constant
  real(kind=8), parameter :: kb=1.381d-16

  ! ev2k   - conversion factor between evs and kelvins
  real(kind=8),parameter  :: ev2k=1.0d0/8.617d-05
  ! ev2erg  - conversion factor between evs and ergs
  real(kind=8),parameter  :: ev2erg=1.602d-12
  ! ev2j   - conversion factor between ergs and Joules
  real(kind=8),parameter  :: erg2j=1d-7
  
  ! The following are scaled to frequency scaling

  ! Scaling factor
  ! this scaling parameter is independent of the main program scaling
  ! (see scaling.f90), and is only used in the radiation physics subroutines
  real(kind=8),parameter  :: sclfre=1.0d15
  ! conversion between evs and frequency
  real(kind=8),parameter  :: ev2fr=0.241838d15/sclfre
  ! h/k, Planck/Boltzm
  real(kind=8),parameter  :: hoverk=47979.72484d0
  ! Planck constant scaled
  real(kind=8),parameter  :: hscl=hplanck*sclfre 
  ! mathematical constant
  real(kind=8),parameter :: pi=3.141592654
  ! tpic2  - 2*pi/c^2 times scaling factors needed for the integral cores
  real(kind=8),parameter :: tpic2=2.0*pi/(c*c)
  ! two_pi_c2  - 2*pi/c^2 times scaling factors needed for the integral cores
  real(kind=8),parameter :: two_pi_c2=2.0*pi/(c*c)*sclfre**3

  ! Hydrogen recombination parameters
  real(kind=dp),parameter :: albpow=-0.7_dp
  real(kind=dp),parameter :: bh00=2.59e-13_dp ! OTS value, alpha_B

  ! Hydrogen collisional ionization
  real(kind=8), parameter :: eth0=13.598
  real(kind=dp),parameter :: xih0=1.0
  real(kind=dp),parameter :: fh0=0.83
  real(kind=dp),parameter :: colh0=1.3e-8*fh0*xih0/(eth0*eth0)
  real(kind=dp),parameter :: temph0=eth0*ev2k
  real(kind=dp),parameter :: hionen=eth0*ev2erg

end module cgsconstants


