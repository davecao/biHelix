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
  use PDBLib
  !use memory use directly?

  ! *-- disable implicit type conversion --*
  implicit none

  ! *-- input --*
  real(DP), dimension(:, :), save, allocatable:: points
  real(DP), dimension(:), save, allocatable:: angles
  character(len=3), dimension(:), save, allocatable:: resname
  integer, dimension(:), save, allocatable::resnum
  character(len=1), dimension(:), save, allocatable:: chId

  real(DP), dimension(3) :: reference_axis = (/0.0, 0.0, 1.0/)
  integer :: nrows, ncols, i, AllocateStatus, DeAllocateStatus
  integer :: nin = 50
  integer :: nout = 6 ! to screen
  integer :: fout = 99  ! to file
  integer :: atNum = 0
  character(6) :: recname
  character(4) :: atName
  character :: altLoc, iCode
  integer :: iostat
  character(LEN=100) :: iomsg

  ! *-- return helix params --*
  real(DP), save, allocatable:: directions(:, :)
  real(DP), save, allocatable:: helix_origins(:, :)
  real(DP) :: tilt, radc, rmsdc, rmsdl, r2 

  ! *-- cla variables --*
  character(len=STRLEN) :: input_filename, output_filename
  logical :: flag, verbose
  logical :: quiet

  ! * -- class atom --*
  type(atom), save, allocatable:: ca_atoms(:) 
  ! * -- input format -- *
  ! 1. 1-6(A6): 'ATOM  '
  ! 2. 7-11(I5): Atom serial number
  ! 3.  blank 1x 
  ! 4. 13-16(A4):Atom name
  ! 5. 17(A1): Alternate location indicator
  ! 6. 18-20(A3): Residume name
  ! 7. 21 blank(1x): blank
  ! 8. 22(A1): chain ID
  ! 9. 23-26 (I4): Residue sequence number
  ! 10. 27 (A1): iCode
  !     blank 3x
  ! 11. 31-38(F8.3): x of CA: x coordinate of CA
  ! 12. 39-46(F8.3): y of CA: y coordinate of CA
  ! 13. 47-54(F8.3): z of CA: z coordinate of CA
  ! *--------------------*
10 format(A6,I5,1x,A4,A1,A3,1x,A1,I4,A1,3x,f8.3,f8.3,f8.3)

  ! * -- output format -- *
  ! 1. 1-6(A6): 'ATOM  '
  ! 2. 7-11(I5): Atom serial number 
  ! 3. 13-16(A4):Atom name
  ! 4. 17(A1): Alternate location indicator
  ! 5. 18-20(A3): Residume name
  ! 6. 21 blank(1x): blank
  ! 7. 22(A1): chain ID
  ! 8. 23-26 (I4): Residue sequence number
  ! 9. 27 (A1): iCode
  ! 10. 31-38(F8.3): x of CA: x coordinate of CA
  ! 11. 39-46(F8.3): y of CA: y coordinate of CA
  ! 12. 47-54(F8.3): z of CA: z coordinate of CA
  ! 13. bending angle of the residue: f8.3 angle in degree
  ! 14. distance to upper layer: f8.3 angstroms to the upper layer
  ! 15. distance to lower layer: f8.3 angstroms to the lower layer
20 format(A6,I5,1x,A4,A1,A3,1x,A1,I4,A1,3x,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3)

  ! *-- cli --*
  ! character(len=STRLEN) :: key
  ! character(len=XSTRLEN) :: description
  ! integer(kind=int_kind) :: kkind
  ! character(len=STRLEN) :: default


  call cla_init
  !(compact, form)
  call cla_register('-i', '--in', 'The input file name', cla_char, 'required')
  call cla_register('-o', '--out', 'The output file name', cla_char, 'helix_out.txt')
  call cla_register('-q', '--quiet', 'Only output to file', cla_flag, 'q')
  call cla_register('-v', '--verbose', 'write more to file', cla_flag, 'v')

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
  open(unit=nin, file=input_filename, status='old', action='read', iostat=iostat)
  if (iostat /= 0) then
    print *, "Failed to open ",input_filename, "status:", iostat
    close(unit=nin)
    stop 
  end if
  ! *-- read size of the matrix: first line --*
  read(nin, *) nrows, ncols

  ! *-- allocate memory --*
  allocate(points(nrows, ncols), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for points'

  allocate(resname(nrows), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for resname'

  allocate(resnum(nrows), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for resnum'

  allocate(chId(nrows), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for chId'

  allocate(angles(nrows-6), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for angles'

  allocate(directions(nrows-3, ncols), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for directions'

  allocate(helix_origins(nrows-2, ncols), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for helix_origins'

  allocate(ca_atoms(nrows), stat=AllocateStatus)
  if (AllocateStatus /= 0) stop 'Failed to allocate memory for atoms'

  do i=1, nrows
    ! read line by line
    read(nin, 10, iostat=iostat, iomsg=iomsg) recname, atNum, atName, altLoc, &
                  resname(i), chId(i), resnum(i), iCode, &
                  points(i, 1), points(i, 2), points(i, 3)
    ca_atoms(i) = atom(recname=recname, atNum=atNum, name=atName, &
                       altLoc=altLoc, resname=resname(i), chId=chId(i), &
                       resnum=resnum(i), iCode=iCode, &
                       x=points(i, 1), y=points(i, 2), z=points(i, 3))
  end do

  if (.not.(quiet)) then
    write(nout, *) 'Read input file: ', input_filename
    write(nout, '("nrows:",I6," ncols:",I6)') size(points, 1), size(points, 2)
    write(nout, '(A15,*(F8.3, 2X))') 'Reference axis', reference_axis(:)
  end if
  ! *-- close file --*
  close(nin)

  ! *-- fit --*
  call fit(points, angles, directions, helix_origins, &
           tilt, radc, rmsdc, rmsdl, r2, &
           reference_axis=reference_axis, info=quiet)
  ! *-- save result --*
  do i=4, nrows - 3
      ca_atoms(i)%bending_angle = angles(i-3)
  end do
  if (.not.(quiet)) then
    write(*, *) 'Scanning finished.'
  end if
  ! *---- write to file ----
  open(unit=fout, file=output_filename, status='replace', action='write')
  if ( .not.(quiet) ) then
    write(*, *) 'Prepare to write to ', output_filename
  end if
  ! Write to the file
  call report(angles, fout, msg="Bending angles:")
  call report(directions, fout, msg="Directions:")
  call report(helix_origins, fout, msg="Helix origins:")
  call report(reference_axis, fout, msg="Reference axis:")
  write(fout, '(A, F12.3)') "Tilt angle w.r.t Reference axis:", tilt

  do i=1, nrows
    !write(fout, 20) ca_atoms(i)%recname, ca_atoms(i)%atNum, ca_atoms(i)%name, &
    !                ca_atoms(i)%altLoc, ca_atoms(i)%resname, &
    !                ca_atoms(i)%chId, ca_atoms(i)%resnum, &
    !                ca_atoms(i)%iCode, &
    !                ca_atoms(i)%x, ca_atoms(i)%y, ca_atoms(i)%z, &
    !                ca_atoms(i)%bending_angle, ca_atoms(i)%d2top, &
    !                ca_atoms(i)%d2bot
    !write(fout, NML=ca_atoms(i))
    call ca_atoms(i)%writef(unit=fout, iostat=iostat, iomsg=iomsg)
  end do
  !if (iostat < 0) then
  !  write(6, *) "EOF before end of data "//trim(iomsg)
  !end if
  if ( .not.(quiet) ) then
    write(*, *) 'Finish writting.'
  end if
  close(fout)
  ! *-- deallocate --*
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

  if(allocated(resname)) then
    deallocate(resname, stat=DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop 'Failed to release memory for resname'
  end if

  if(allocated(resnum)) then
    deallocate(resnum, stat=DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop 'Failed to release memory for resnum'
  end if

  if(allocated(chId)) then
    deallocate(chId, stat=DeAllocateStatus)
    if (DeAllocateStatus /= 0) stop 'Failed to release memory for chId'
  end if
end program helix_param
