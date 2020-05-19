!==============================================================================!
  subroutine Field_Mod_Grad_Pressure_Correction(flow, pp)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure correction.                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  type(Var_Type)   :: pp
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  call Grid_Mod_Exchange_Real(grid, pp % n)

  !---------------------------------!
  !   No correction at boundaries   !
  !---------------------------------!

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if((flow % p % bnd_cond_type(c2) .eq. INFLOW)  .or.  &
       (flow % p % bnd_cond_type(c2) .eq. WALL)    .or.  &
       (flow % p % bnd_cond_type(c2) .eq. WALLFL)  .or.  &
       (flow % p % bnd_cond_type(c2) .eq. OUTFLOW) .or.  &
       (flow % p % bnd_cond_type(c2) .eq. CONVECT) .or.  &
       (flow % p % bnd_cond_type(c2) .eq. SYMMETRY) ) then
      pp % n(c2) = 0.0
    end if  ! .not. PRESSURE
  end do  ! 1, grid % n_faces

  call Field_Mod_Grad_Component(flow, pp % n, 1, pp % x)  ! dp/dx
  call Field_Mod_Grad_Component(flow, pp % n, 2, pp % y)  ! dp/dy
  call Field_Mod_Grad_Component(flow, pp % n, 3, pp % z)  ! dp/dz

  end subroutine
