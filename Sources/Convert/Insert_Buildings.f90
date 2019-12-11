!==============================================================================!
  subroutine Insert_Buildings(grid)
!------------------------------------------------------------------------------!
!   Sorts cells by their height (z coordinate)                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Grid_Mod,  only: Grid_Type
  use Sort_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(grid_type) :: grid
!------------------------------------------------------------------------------!
  include 'Cell_Numbering_Neu.f90'
!-----------------------------------[Locals]-----------------------------------!
  integer              :: run, s, c, cg, nr_cg, cu, n, ni, dir, bc, n_ground_cells
  real                 :: dis2, dis2_min
  integer, allocatable :: ground_cell(:)
  integer, allocatable :: new_c(:)
  integer, allocatable :: old_c(:)
  integer, allocatable :: i_work_1(:)
  integer, allocatable :: i_work_2(:,:)
  integer, allocatable :: i_work_3(:,:)
  integer, allocatable :: criteria(:,:)
  integer              :: n_hor_layers, cnt, ground_n
  character(len=80)    :: bc_name
  real                 :: height
  logical, allocatable :: cell_in_building(:)
  logical, allocatable :: node_on_building(:)
  logical              :: buildings_exist
  integer              :: fn(6,4), n_f_nod, f_nod(4)
!==============================================================================!

  !------------------------------------------------!
  !                                                !
  !   Phase I: Align horizontal cell coordinates   !
  !                                                !
  !------------------------------------------------!

  !---------------------------------------------!
  !   Calculate cell centers because you will   !
  !     be sorting them by height later on      !
  !---------------------------------------------!
  call Calculate_Cell_Centers(grid)

  !----------------------------------------!
  !   Find and store cells on the ground   !
  !----------------------------------------!
  do run = 1, 2
    n_ground_cells = 0
    do c = 1, grid % n_cells
      do dir = 1, 6
        if(grid % cells_bnd_color(dir, c) .ne. 0) then
          bc_name = trim(grid % bnd_cond % name(grid % cells_bnd_color(dir, c)))
          call To_Upper_Case(bc_name)
          if(bc_name(1:8) .eq. 'BUILDING' .or.  &
             bc_name      .eq. 'GROUND') then

            n_ground_cells = n_ground_cells + 1
            if(run .eq. 2) ground_cell(n_ground_cells) = c
          end if
        end if
      end do
    end do
    if(run .eq. 1) then
      allocate(ground_cell(n_ground_cells)); ground_cell(:) = 0
    end if
  end do
  print '(a38,i9)', '# Number of cells on the ground:     ', n_ground_cells
  print '(a38,i9)', '# Number of horizontal layers (1)    ', grid % n_cells  &
                                                           / n_ground_cells

  !---------------------------------------------------------------------!
  !   For each cell, find the nearest on the ground and align with it   !
  !---------------------------------------------------------------------!
  print *, '#=================================================================='
  print *, '# Aligning cell coordinates.  This may take a few minutes          '
  print *, '#------------------------------------------------------------------'
  do c = 1, grid % n_cells
    dis2_min = HUGE

    ! Print progress on the screen
    if(mod(c, (grid % n_cells / 20) ) .eq. 0) then
      print '(a2, f5.0, a14)',   &
            ' #', (100. * c / (1.0*(grid % n_cells))), ' % complete...'
    end if ! each 5%

    do n = 1, n_ground_cells
      cg = ground_cell(n)     ! real ground cell number

      ! Compute planar distance
      dis2 = (  (grid % xc(c) - grid % xc(cg)) ** 2  &
              + (grid % yc(c) - grid % yc(cg)) ** 2)

      ! Store ground cell if nearest
      if(dis2 < dis2_min) then
        dis2_min = dis2
        nr_cg = cg
      end if
    end do

    ! Simply take planar coordinates from the nearest wall cell
    grid % xc(c) = grid % xc(nr_cg)
    grid % yc(c) = grid % yc(nr_cg)
  end do

  !----------------------------------------------------!
  !                                                    !
  !   Phase II: Sort cells in columns in z direction   !
  !                                                    !
  !----------------------------------------------------!

  allocate(i_work_1(   grid % n_cells));     i_work_1(:)   = 0
  allocate(i_work_2(8, grid % n_cells));     i_work_2(:,:) = 0
  allocate(i_work_3(6, grid % n_cells));     i_work_3(:,:) = 0
  allocate(new_c   (   grid % n_cells));     new_c(:)      = 0
  allocate(old_c   (   grid % n_cells));     old_c(:)      = 0
  allocate(criteria(   grid % n_cells, 3));  criteria(:,:) = 0.0

  !------------------------!
  !   Store cells' nodes   !
  !------------------------!
  do c = 1, grid % n_cells
    i_work_1(c)      = grid % cells_n_nodes(c)
    i_work_2(1:8, c) = grid % cells_n(1:8, c)
    i_work_3(1:6, c) = grid % cells_bnd_color(1:6, c)
  end do

  !--------------------------!
  !   Set sorting criteria   !
  !--------------------------!
  do c = 1, grid % n_cells
    criteria(c, 1) = nint(grid % xc(c) * GIGA)
    criteria(c, 2) = nint(grid % yc(c) * GIGA)
    criteria(c, 3) = nint(grid % zc(c) * GIGA)
    old_c(c)       = c
  end do

  !--------------------------------------------------!
  !   Sort new numbers according to three criteria   !
  !--------------------------------------------------!
  call Sort_Mod_3_Int_Carry_Int(criteria(1:grid % n_cells, 1),  &
                                criteria(1:grid % n_cells, 2),  &
                                criteria(1:grid % n_cells, 3),  &
                                old_c   (1:grid % n_cells))

  ! This was a bit of a bluff but it worked
  do c = 1, grid % n_cells
    new_c(old_c(c)) = c
  end do

  !-----------------------------------------------!
  !   Do the sorting of data pertinent to cells   !
  !-----------------------------------------------!
  do c = 1, grid % n_cells
    grid % cells_n_nodes       (new_c(c)) = i_work_1(c)
    grid % cells_n        (1:8, new_c(c)) = i_work_2(1:8, c)
    grid % cells_bnd_color(1:6, new_c(c)) = i_work_3(1:6, c)
  end do
  call Sort_Mod_Real_By_Index(grid % xc   (1), new_c(1), grid % n_cells)
  call Sort_Mod_Real_By_Index(grid % yc   (1), new_c(1), grid % n_cells)
  call Sort_Mod_Real_By_Index(grid % zc   (1), new_c(1), grid % n_cells)

  !----------------------------------------!
  !                                        !
  !   Phase III: Find cells in buildings   !
  !                                        !
  !----------------------------------------!

  ! Estimate number of horizontal layers again (to check if sorting worked)
  do c = 1, grid % n_cells
    if(grid % zc(c+1) < grid % zc(c)) then
      n_hor_layers = c
      goto 1
    end if
  end do
1 continue
  print '(a38,i9)', '# Number of horizontal layers (2):   ', n_hor_layers

  !-----------------------------!
  !   Find cells in buildings   !
  !-----------------------------!
  buildings_exist = .false.
  allocate(cell_in_building(grid % n_cells))
  cell_in_building(:) = .false.
  do c = 1, grid % n_cells
    do dir = 1, 6
      if(grid % cells_bnd_color(dir, c) .ne. 0) then
        bc_name = trim(grid % bnd_cond % name(grid % cells_bnd_color(dir, c)))
        call To_Upper_Case(bc_name)
        if(bc_name(1:8) .eq. 'BUILDING') then
          read(bc_name(10:12), *) height
          height = height / 1000.0
          do cu = c, c + n_hor_layers - 1
            if(grid % zc(cu) - height >= grid % zc(c)) goto 2
            cell_in_building(cu) = .true.        ! mark cell as in building
            buildings_exist = .true.             ! at least one building found
          end do
2         continue
        end if
      end if
    end do
  end do

  ! Count the number of cells in buildings (this is for checking only)
  cnt = 0
  do c = 1, grid % n_cells
    if(cell_in_building(c)) cnt = cnt + 1
  end do
  print '(a38,i9)', '# Number of cells in buildings:      ', cnt

  !-----------------------------!
  !   Mark nodes on buildings   !
  !-----------------------------!
  allocate(node_on_building(grid % n_nodes))
  node_on_building(:) = .false.
  do c = 1, grid % n_cells
    if(cell_in_building(c)) then
      do ni = 1, grid % cells_n_nodes(c)  ! mark node as on building
        n = grid % cells_n(ni, c)
        node_on_building(n) = .true.
      end do
    end if
  end do

  !-----------------------------------------------------------!
  !                                                           !
  !   Phase IV: Manage boundary condition names and numbers   !
  !                                                           !
  !-----------------------------------------------------------!

  ! Find ground b.c. number
  do bc = 1, grid % n_bnd_cond
    bc_name = trim(grid % bnd_cond % name(bc))
    call To_Upper_Case(bc_name)
    if(bc_name .eq. 'GROUND') then
      ground_n = bc
    end if
  end do

  ! Eliminate building_000 b.c.
  do bc = 1, grid % n_bnd_cond
    bc_name = trim(grid % bnd_cond % name(bc))
    call To_Upper_Case(bc_name)
    if(bc_name .eq. 'BUILDING_000') then
      do c = 1, grid % n_cells
        do dir = 1, 6
          if(grid % cells_bnd_color(dir, c) .eq. bc) then
            grid % cells_bnd_color(dir, c) = ground_n
          end if
        end do
      end do
    end if
  end do

  cnt   = 0
  do bc = 1, grid % n_bnd_cond
    bc_name = trim(grid % bnd_cond % name(bc))
    call To_Upper_Case(bc_name)
    if(bc_name(1:8) .ne. 'BUILDING') then
      cnt = cnt + 1
      grid % bnd_cond % name(cnt) = grid % bnd_cond % name(bc)
      do c = 1, grid % n_cells
        do dir = 1, 6
          if(grid % cells_bnd_color(dir, c) .eq. bc) then
            grid % cells_bnd_color(dir, c) = cnt
          end if
        end do
      end do
    end if
  end do
  grid % n_bnd_cond = cnt

  ! Shift all by one up (to be able to insert building walls as first)
  if(buildings_exist) then
    do bc = grid % n_bnd_cond, 1, -1
      grid % bnd_cond % name(bc+1) = grid % bnd_cond % name(bc)
      do c = 1, grid % n_cells
        do dir = 1, 6
          if(grid % cells_bnd_color(dir, c) .eq. bc) then
            grid % cells_bnd_color(dir, c) = bc+1
          end if
        end do
      end do
    end do
    grid % n_bnd_cond = grid % n_bnd_cond + 1

    ! Add first boundary condition for walls
    grid % bnd_cond % name(1) = 'BUILDING_WALLS'
  end if

  !----------------------------------------------------------!
  !                                                          !
  !   Phase V: Store only cells which are not in buildings   !
  !                                                          !
  !----------------------------------------------------------!

  !-----------------------------------------------------!
  !   Store cells' nodes and boundary conditons again   !
  !-----------------------------------------------------!
  do c = 1, grid % n_cells
    i_work_1(c)      = grid % cells_n_nodes(c)
    i_work_2(1:8, c) = grid % cells_n(1:8, c)
    i_work_3(1:6, c) = grid % cells_bnd_color(1:6, c)
  end do

  !------------------------------------------------------!
  !   Renumber the cells without the cells in building   !
  !------------------------------------------------------!
  new_c(:) = 0
  cnt      = 0
  do c = 1, grid % n_cells
    if(.not. cell_in_building(c)) then
      cnt = cnt + 1
      new_c(c) = cnt
    end if
  end do
  print *, '# Old number of cells:                     ', grid % n_cells
  print *, '# Number of cells with excluded buildings: ', cnt
  print *, '# Number of excluded cells: ',                grid % n_cells - cnt

  !--------------------------------------!
  !   Information on cell connectivity   !
  !--------------------------------------!
  do c = 1, grid % n_cells
    if(new_c(c) .ne. 0) then
      grid % cells_n_nodes       (new_c(c)) = i_work_1(c)
      grid % cells_n        (1:8, new_c(c)) = i_work_2(1:8, c)
      grid % cells_bnd_color(1:6, new_c(c)) = i_work_3(1:6, c)
    end if
  end do
  grid % n_cells = cnt

  !-----------------------------------------!
  !                                         !
  !   Phase VI: Insert new boundary faces   !
  !                                         !
  !-----------------------------------------!
  if(buildings_exist) then

    do c = 1, grid % n_cells

      if(grid % cells_n_nodes(c) .eq. 4) fn = neu_tet
      if(grid % cells_n_nodes(c) .eq. 5) fn = neu_pyr
      if(grid % cells_n_nodes(c) .eq. 6) fn = neu_wed
      if(grid % cells_n_nodes(c) .eq. 8) fn = neu_hex

      do dir = 1, 6
        if(grid % cells_bnd_color(dir, c) .eq. 0) then

          n_f_nod    = 0
          f_nod(1:4) = -1

          ! Fill up face nodes
          do ni = 1, 4
            if(fn(dir, ni) > 0) then
              f_nod(ni) = grid % cells_n(fn(dir, ni), c)
              n_f_nod   = n_f_nod + 1
            end if
          end do

          if( n_f_nod > 0 ) then

            ! Qadrilateral face
            if(f_nod(4) > 0) then
              if( node_on_building(f_nod(1)) .and.  &
                  node_on_building(f_nod(2)) .and.  &
                  node_on_building(f_nod(3)) .and.  &
                  node_on_building(f_nod(4)) ) then
                grid % cells_bnd_color(dir, c) = 1
              end if
            else
              if( node_on_building(f_nod(1)) .and.  &
                  node_on_building(f_nod(2)) .and.  &
                  node_on_building(f_nod(3)) ) then
                grid % cells_bnd_color(dir, c) = 1
              end if
            end if
          end if
        end if
      end do
    end do

  end if

  end subroutine