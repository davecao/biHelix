! * -- ATOM format -- *
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
! 14. 55-62(F8.3): phi of the residue
! 15. 63-70(F8.3): psi of the residue
! 16. 71-78(F8.3): area of the residue
! 17. 79-80(A2): type of the secondary structure(by Stride)
! *--------------------*
! 10 format(A6,I5,1x,A4,A1,A3,1x,A1,I4,A1,3x,f8.3,f8.3,f8.3)
!
! * -- HEADER format -- *
! 1. 1-6(A6): 'HEADER'
! 2. 5x
! 3. 11-50(A40): classification
! 4. 51-59(A9): depDate
! 5. 4x
! 6. 63-66(A4): idCode, pdbid
! ************************************************************************
! 30 format(A6,5x,A40,A9,3x,A4)
! ------------------------------------------------------------------------

module io
  !****************
  ! other modules
  !****************
  use precision 
  use atom_mod
  use group_mod

  implicit none ! use strong type

  private  ! hide the type-bound procedure implementation procedures
  public :: readline

contains
  function readline(filename, kinkOnly) result(group_obj)
    character(*), intent(in) :: filename
    logical, intent(in) :: kinkOnly
    character(len=100) :: buffer, label
    integer, parameter :: fh = 15
    integer :: pos = 6
    integer :: ios = 0
    integer :: line = 0
    integer :: atNum = 0
    integer :: resnum = 0
    real(DP) :: x, y, z
    real(DP) :: phi, psi, area
    character(len=2) :: stype
    character(len=3) :: resname
    character(6) :: recname = 'ATOM'
    character(4) :: atName
    character(len=1) :: chId, altLoc, iCode
    character(len=40) :: cls
    character(len=9) :: depDate
    character(len=4) :: idCode
    logical :: hasUpper = .FALSE.
    logical :: hasLower = .FALSE.

    !class(group), allocatable :: group_obj
    type(group):: group_obj
    type(atom) :: atom_obj
10 format(I5,1x,A4,A1,A3,1x,A1,I4,A1,3x,f8.3,f8.3,f8.3,f8.3,f8.3,f8.3,A2)
20 format(4x,A40,A9,3x,A4)
    ! create a group of atoms: initialize with 10 atoms
    group_obj = group(10)
    ! open files
    open(fh, file=filename, status='old', action='read')
    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    do while (ios == 0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios == 0) then
        line = line + 1
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)
          case ('HEADER')
            read(buffer, 20, iostat=ios) cls, depDate, idCode
            group_obj%classification = cls
            group_obj%depDate = depDate
            group_obj%idCode = idCode
          case ('ATOM  ')
            read(buffer, 10, iostat=ios) atNum, atName, altLoc, resname, &
                                         chId, resnum, iCode, &
                                         x, y, z, phi, psi, area, stype
            atom_obj = atom(recname=recname, atNum=atNum, name=atName, &
                          altLoc=altLoc, resname=resname, chId=chId, &
                          resnum=resnum, iCode=iCode, &
                          x=x, y=y, z=z, phi=phi, psi=psi, area=area, &
                          ss=stype)
            call group_obj%add(atom_obj)
          case ('REFAXS')
            !read(buffer,'(1x,f8.3,f8.3,f8.3)')
          case ('MEMUPP')
            ! membrane upper center
            read(buffer,'(1x,3f8.3)') group_obj%mem_info%upcenter(1:3)
          case ('MEMUNM')
            ! membrane upper normal
            read(buffer,'(1x,3f8.3)') group_obj%mem_info%upNormVec(1:3)
            hasUpper = .TRUE.
          case ('MEMLOW')
            ! membrane lower center
            read(buffer,'(1x,3f8.3)') group_obj%mem_info%lowcenter(1:3)
          case ('MEMLNM')
            ! membrane lower normal
            read(buffer,'(1x,3f8.3)') group_obj%mem_info%lowNormVec(1:3)
            hasLower = .TRUE.
          case default
        end select
      end if
    end do
    if (.not.kinkOnly) then
      call group_obj%setMemLayerNum(hasUpper=hasUpper, hasLower=hasLower)
    end if
    ! close file
    close(fh)
  end function readline

end module io
