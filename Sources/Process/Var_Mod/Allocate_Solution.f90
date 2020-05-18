!==============================================================================!
  subroutine Var_Mod_Allocate_Solution(phi, grid, name_phi, name_flux)
!------------------------------------------------------------------------------!
!   This is to allocate a variable for a solution with usual algorithm.        !
!   Variables such as velocities and pressures should be allocated with it.    !
!                                                                              !
!   One could think of storing pointer to the grid as well.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type)          :: phi
  type(Grid_Type), target :: grid
  character(len=*)        :: name_phi
  character(len=*)        :: name_flux
  integer                 :: nc, nb, nf
!==============================================================================!

  ! Store grid for which the variable is defined
  phi % pnt_grid => grid
  nc = grid % n_cells
  nb = grid % n_bnd_cells
  nf = grid % n_faces

  ! Store variable name
  phi % name      = name_phi
  phi % flux_name = name_flux

  ! Values new (n), old (o), and older than old (oo)
  allocate (phi % n (-nb: nc));  phi % n  = 0.
  allocate (phi % o (-nb: nc));  phi % o  = 0.
  allocate (phi % oo(-nb: nc));  phi % oo = 0.

  ! Advection terms
  allocate (phi % a(nc));  phi % a = 0.

  ! Cross diffusion terms
  allocate (phi % c(nc));  phi % c = 0.

  ! Variable's boundary value
  allocate (phi % b(-nb: -1));  phi % b = 0.

  ! Variable's boundary flux
  allocate (phi % q(-nb: -1));  phi % q = 0.

  ! Boundary cell type
  allocate (phi % bnd_cond_type(-nb : nf));  phi % bnd_cond_type = 0

  ! Gradients
  allocate (phi % x(-nb : nc));  phi % x = 0.0
  allocate (phi % y(-nb : nc));  phi % y = 0.0
  allocate (phi % z(-nb : nc));  phi % z = 0.0

  ! Min and max
  allocate (phi % min(-nb : nc));  phi % min = 0.0
  allocate (phi % max(-nb : nc));  phi % max = 0.0

  end subroutine
