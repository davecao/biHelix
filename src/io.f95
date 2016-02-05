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
  function readline(filename) result(group_obj)
    character(*), intent(in) :: filename
    character(len=100) :: buffer, label
    integer, parameter :: fh = 15
    integer :: pos = 6
    integer :: ios = 0
    integer :: line = 0
    integer :: atNum = 0
    integer :: resnum = 0
    real(DP) :: x, y, z
    character(len=3) :: resname
    character(6) :: recname = 'ATOM'
    character(4) :: atName
    character(len=1) :: chId, altLoc, iCode
    !class(group), allocatable :: group_obj
    type(group):: group_obj
    type(atom) :: atom_obj

10 format(I5,1x,A4,A1,A3,1x,A1,I4,A1,3x,f8.3,f8.3,f8.3)
    
    ! create a group of atoms: initialize with 10 atoms
    group_obj = group(10)
    ! open files
    open(fh, file=filename)
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
          case ('ATOM  ')
            read(buffer, 10, iostat=ios) atNum, atName, altLoc, resname, &
                                         chId, resnum, iCode, &
                                         x, y, z
            atom_obj = atom(recname=recname, atNum=atNum, name=atName, &
                          altLoc=altLoc, resname=resname, chId=chId, &
                          resnum=resnum, iCode=iCode, &
                          x=x, y=y, z=z)
            call group_obj%add(atom_obj)
          case ('REFAXS')
            !read(buffer,'(1x,f8.3,f8.3,f8.3)') 
          case default
        end select
      end if
    end do
    ! close file
    close(fh)
  end function readline
end module io
