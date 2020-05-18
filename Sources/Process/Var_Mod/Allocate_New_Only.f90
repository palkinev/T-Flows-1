!==============================================================================!
  subroutine Var_Mod_Allocate_New_Only(phi, grid, name_phi)
!------------------------------------------------------------------------------!
!   This is to allocate a simplified uknown, holding only current value,       !
!   such as pressure for example.                                              !
!                                                                              !
!   One could think of storing pointer to the grid as well.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)          :: phi
  type(Grid_Type), target :: grid
  character(len=*)        :: name_phi
  integer                 :: nb, nc, nf
!==============================================================================!

  ! Store grid for which the variable is defined
  phi % pnt_grid => grid
  nb = grid % n_bnd_cells
  nc = grid % n_cells
  nf = grid % n_faces

  ! Store variable name
  phi % name = name_phi

  ! Values in the new (n) time step
  allocate (phi % n(-nb : nc));  phi % n = 0.0

  ! Gradients
  allocate (phi % x(-nb : nc));  phi % x = 0.0
  allocate (phi % y(-nb : nc));  phi % y = 0.0
  allocate (phi % z(-nb : nc));  phi % z = 0.0

  ! Variable's boundary value
  allocate (phi % b(-nb: -1));  phi % b = 0.

  ! Boundary cell type
  if(name_phi == 'P') then
    allocate (phi % bnd_cond_type(-nb : nf));  phi % bnd_cond_type = 0
  end if

  end subroutine
