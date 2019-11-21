!==============================================================================!
  subroutine Vof_Max_Courant_Number(mult, dt, c_d, interf, courant_max, n)
!------------------------------------------------------------------------------!
!   This function aftar VOF is solved                                          !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Field_Type),     pointer :: flow
  type(Grid_Type),      pointer :: grid
  type(Var_Type),       pointer :: vof
  type(Face_Type),      pointer :: v_flux
  real                          :: dt     ! time step
  real                          :: c_d(-mult % pnt_grid % n_bnd_cells:  &
                                        mult % pnt_grid % n_cells) ! Courant n.
!--------------------------------[Locals]--------------------------------------!
  integer                       :: c, c1, c2, s, last_cell, interf, n
  real                          :: b_volume, surface, rise_velocity,           &
                                   circularity, c_position, courant_max
  real                          :: vof_dist
!==============================================================================!
  
  ! Take aliases
  flow     => mult % pnt_flow
  v_flux   => flow % v_flux
  grid     => flow % pnt_grid
  vof      => mult % vof

  courant_max = 0.0

  !Initialize:

  c_d(-mult % pnt_grid % n_bnd_cells:mult % pnt_grid % n_cells) = 0.0

  if (interf == 1) then

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      
      if (c2 > 0) then

        vof_dist = min(max(vof % n(c1), 0.0),1.0)

        vof_dist = (1.0 - vof_dist) * (1.0 - vof_dist)            &
                                    * vof_dist * vof_dist * 16.0    

        c_d(c1) = c_d(c1) + vof_dist * max(- v_flux % n(s)         &
                                     * dt / grid % vol(c1), 0.0)

        vof_dist = min(max(vof % n(c2), 0.0),1.0)

        vof_dist = (1.0 - vof_dist) * (1.0 - vof_dist)            &
                                    * vof_dist * vof_dist * 16.0    

        c_d(c2) = c_d(c2) + vof_dist * max( v_flux % n(s)         &
                                     * dt / grid % vol(c2), 0.0)

      ! Side is on the boundary
      else ! (c2 < 0)
        vof_dist = min(max(vof % n(c1), 0.0),1.0)

        vof_dist = (1.0 - vof_dist) * (1.0 - vof_dist)            &
                                    * vof_dist * vof_dist * 16.0    

        c_d(c1) = c_d(c1) + vof_dist * max(-v_flux % n(s)         &
                                     * dt / grid % vol(c1), 0.0)
      end if

    end do

    do c = 1, grid % n_cells
      courant_max = max(c_d(c), courant_max)    
    end do

  else

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      
      if (c2 > 0) then

        c_d(c1) = c_d(c1) + max(-v_flux % n(s) * dt / grid % vol(c1), 0.0)

        c_d(c2) = c_d(c2) + max( v_flux % n(s) * dt / grid % vol(c2), 0.0)

      ! Side is on the boundary
      else ! (c2 < 0)

        c_d(c1) = c_d(c1) + max(-v_flux % n(s) * dt / grid % vol(c1), 0.0)

      end if

    end do

  end if

  call Comm_Mod_Global_Max_Real(courant_max)
  call Comm_Mod_Exchange_Real(grid, c_d)

  end subroutine
