!
! gfortran -cpp claf90/kinds.f90 claf90/cla.f90 claf90/helanal.f95 main.f95 -o test -I./claf90 -llapack
!
!gfortran -g -cpp fortran/claf90/kinds.f90 fortran/claf90/cla.f90 fortran/claf90/helanal.f95 fortran/main.f95 -o test -llapack

program helix_param
  ! *-- modules --*
  use kinds
  use cla
  use helanal
  use utilities
  use io, io_readline=>readline
  use atom_mod
  use group_mod
  use strings

  ! *-- disable implicit type conversion --*
  implicit none

  ! *-- input --*
  real(DP), dimension(:, :), save, allocatable:: points
  real(DP), dimension(:), save, allocatable:: angles

  real(DP), dimension(3) :: reference_axis = (/999.0, 0.0, 0.0/)
  real(DP), dimension(3) :: helix_axis = (/999.0, 0.0, 0.0/)
  !real(DP), dimension(3) :: upper = (/999.0, 0.0, 0.0/)
  !real(DP), dimension(3) :: lower = (/999.0, 0.0, 0.0/)
  !real(DP), dimension(3) :: mem_normal = (/0.0, 0.0, 0.0/)
  !character(len=STRLEN) :: axis_str
  !character(len=STRLEN), dimension(3) :: tokens
  !integer :: ntokens
  integer :: natoms, ncols, i, AllocateStatus, DeAllocateStatus

  !integer :: nin = 50
  !integer :: nout = 6 ! to screen
  integer, parameter :: fout = 99 ! to file
  integer :: natoms_threshold = 7 ! At least 7 a.a. is needed.
  integer :: iostat
  character(LEN=100) :: iomsg

  ! *-- return helix params --*
  real(DP), save, allocatable:: directions(:, :)
  real(DP), save, allocatable:: helix_origins(:, :)
  real(DP) :: tilt
  !real(DP) :: radc, rmsdc, rmsdl, r2 
 
  ! *-- cla variables --*
  character(len=STRLEN) :: input_filename, output_filename
  logical :: flag, verbose
  logical :: quiet

  ! * -- class atom --*
  type(group) :: atomGroup

  ! *-- cli --*
  call cla_init
  !(compact, form)
  call cla_register('-i', '--in', 'The input file name.', cla_char, 'required')
  call cla_register('-o', '--out', 'The output file name.', cla_char, 'helix_out.txt')
!  call cla_register('-r', '--ref', 'x, y, z of a reference axis.', cla_char, '0.0, 0.0, 1.0')
!  call cla_register('-u', '--upper', 'The center of the upper bilayer.', cla_char, '0.0, 0.0, 15.0')
!  call cla_register('-l', '--lower', 'The center of the lower bilayer.', cla_char, '0.0, 0.0, -15.0')
  call cla_register('-q', '--quiet', 'Only output to file.', cla_flag, 'q')
  call cla_register('-v', '--verbose', 'write more to file.', cla_flag, 'v')

  ! *-- processing command line arguments --*
  !call cla_validate
  ! -------- -i input_filename ------------
  flag = cla_key_present('-i')
  if (flag) then
    call cla_get('-i', input_filename)
  else
    call cla_help("biHelix")
    stop
  end if
  ! -------- -r input_filename ------------
!  flag = cla_key_present('-r')
!  if (flag) then
!    call cla_get('-r', axis_str)
!    call parse(axis_str, ',', tokens, ntokens)
!    do i=1, ntokens
!      read(tokens(i), '(F8.3)') reference_axis(i)
!    end do
!  end if
!  ! -------- -u input_filename ------------
!  flag = cla_key_present('-u')
!  if (flag) then
!    call cla_get('-u', axis_str)
!    call parse(axis_str, ',', tokens, ntokens)
!    do i=1, ntokens
!      read(tokens(i), '(F8.3)') upper(i)
!    end do
!  end if
!  ! -------- -l input_filename ------------
!  flag = cla_key_present('-l')
!  if (flag) then
!    call cla_get('-l', axis_str)
!    call parse(axis_str, ',', tokens, ntokens)
!    do i=1, ntokens
!      read(tokens(i), '(F8.3)') lower(i)
!    end do
!  end if
  ! -------- -o output_filename -----------
  flag = cla_key_present('-o')
  if (flag) then
    call cla_get('-o', output_filename)
  else
    output_filename = 'helix_out.txt'
  endif
  ! -------- -s only to screen -----------
  flag = cla_key_present('-q')
  if (flag) then
    quiet = .TRUE.
  else
    quiet = .FALSE.
  endif
  ! -------- -v verbose --------------------
  flag = cla_key_present('-v')
  if (flag) then
    verbose = .TRUE.
  else
    verbose = .FALSE.
  endif

  ! *-- load data --*
  atomGroup = io_readline(input_filename)
  reference_axis = atomGroup%getReferenceAxis()
  !atomGroup%reference_axis = reference_axis
  !atomGroup%mem_info%upcenter = upper
  !atomGroup%mem_info%lowcenter = lower
  !atomGroup%mem_normal = upper - lower
  natoms = atomGroup%getNumAtoms()
  if (natoms < natoms_threshold) then
    print *, input_filename, " - ",natoms, " atoms: # of CA atoms is >=7 a.a."
    stop
  end if
  ncols = 3
  points = atomGroup%getCoords()

  allocate(angles(natoms-6), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for angles'

  allocate(directions(natoms-3, ncols), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for directions'

  allocate(helix_origins(natoms-2, ncols), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for helix_origins'

  if (.not.(quiet)) then
    write(stdout, *) 'Read input file: ', input_filename
    write(stdout, '("natoms:",I6," ncols:",I6)') size(points, 1), size(points, 2)
    write(stdout, '(A15,*(F8.3, 2X))') 'Reference axis', reference_axis(:)
  end if
  ! *-- set distance to the membrane --*
  call atomGroup%findDistToMem() 
  ! *-- fit --*
  call fit(points, angles, directions, helix_origins, helix_axis, &
           tilt, reference_axis=reference_axis, info=quiet)
  ! *-- save result --*
  do i=4, natoms - 3
    call atomGroup%setAtomBendingAngleAt(i, angles(i-3))
  end do
  atomGroup%directions = directions
  atomGroup%helix_origins = helix_origins
  atomGroup%helix_axis = helix_axis
  atomGroup%tilt = tilt

  if (.not.(quiet)) then
    write(*, *) 'Scanning finished.'
  end if

  ! *---- write to file ----
  open(fout, file=output_filename, status='replace', action='write')

  call atomGroup%printf(fout)
  close(fout, iostat=iostat, iomsg=iomsg, status="keep")
  if ( iostat /= 0 ) write(*, *) "Error closing file: ", iomsg
  
  if ( .not.(quiet) ) then
    write(*, *) 'Finish writting.'
  end if

  ! *-- clean memory --*
  if(allocated(points)) then
    deallocate(points, stat=DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop 'Failed to release memory for points'
  end if

  if(allocated(angles)) then
    deallocate(angles, stat=DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop 'Failed to release memory for angles'
  end if

  if(allocated(directions)) then
    deallocate(directions, stat=DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop 'Failed to release memory for directions'
  end if

  if(allocated(helix_origins)) then
    deallocate(helix_origins, stat=DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop 'Failed to release memory for helix_origins'
  end if

end program helix_param
