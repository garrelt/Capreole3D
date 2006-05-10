module ionic

  ! This module contains the interface between the doric module and the
  ! mpi_hydro program. It contains two subroutines:
  ! - init_ionic: initializes all the radiation variables, including
  !         the doric initialization
  ! - update_ionic: updates ionic quantities

  ! Version: dummy routines, no atomic physics applied

  use precision
  use sizes
  use mesh
  use grid
  use hydro

  implicit none

contains

  subroutine init_ionic(r_interface)

    ! This subroutine initializes the atomic/ionic fractions
    ! Garrelt Mellema

    ! Version: dummy

    real(kind=dp),intent(in) :: r_interface

  end subroutine init_ionic

  subroutine update_ionic(dt)

    ! updates the radiation connected quantities one timestep
    ! garrelt mellema
    
    ! Version: dummy for runs without atomic physics

    real(kind=dp),intent(in) :: dt
    
  end subroutine update_ionic
    
end module ionic


