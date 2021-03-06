!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Energy(flow, turb, mult, dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Energy function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  real                          :: dt    ! time step
  integer                       :: ini   ! inner iteration
!==============================================================================!

  end subroutine
