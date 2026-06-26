!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.04
!
!	The module that contains support functions for input-output used within artemide
!   low-level infrastructure module; changes should be rare and compatibility-preserving
!
!   -- Remove of older and unusing functions
!				A.Vladimirov (19.06.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aTMDe_IO
use aTMDe_numerics
implicit none

private

!!!! colors for coloring of the ansi-output
character(len=1), public, parameter :: c_esc = achar(27)
character(len=2), public, parameter :: c_start = c_esc // '['
character(len=1), public, parameter :: c_end = 'm'
character(len=*), public, parameter :: c_black = '30'
character(len=*), public, parameter :: c_black_bold = '30;1'
character(len=*), public, parameter :: c_red = '31'
character(len=*), public, parameter :: c_red_bold = '31;1'
character(len=*), public, parameter :: c_green = '32'
character(len=*), public, parameter :: c_green_bold = '32;1'
character(len=*), public, parameter :: c_yellow = '33'
character(len=*), public, parameter :: c_yellow_bold = '33;1'
character(len=*), public, parameter :: c_blue = '34'
character(len=*), public, parameter :: c_magenta = '35'
character(len=*), public, parameter :: c_cyan = '36'
character(len=*), public, parameter :: c_white = '37'
character(len=*), public, parameter :: c_clear = c_start // '0' // c_end


public::color, MoveTO, numToStr, int4ToStr
public::WarningString, ErrorString

type, public :: Warning_OBJ
  character(:), allocatable::moduleName
  integer::messageCounter = 0
  integer::messageTrigger = 3
contains
  procedure:: Reset =>Reset_def
  procedure:: WarningRaise => WarningRaise_def
end type


interface numToStr
  module procedure intToStr,spToStr,dpToStr
end interface numToStr

interface ErrorString
  module procedure ErrorString1,ErrorString2
end interface ErrorString
  
contains

!!!!!!!!!!!! procedures that define the class Warning_OBJ
!!!!! resets the counter of warning handler
subroutine Reset_def(this)
class(Warning_OBJ), intent(inout)::this

this%messageCounter=0
end subroutine Reset_def

!!!!! writes the warning message and increase +1 the message counter
!!!!! one the limit is hit, stops show Warnings
subroutine WarningRaise_def(this,str)
class(Warning_OBJ), intent(inout)::this
character(len=*), intent(in) :: str

if(this%messageCounter<this%messageTrigger) then
  write(*,*) WarningString(str,this%moduleName)
  this%messageCounter=this%messageCounter+1

  if(this%messageCounter==this%messageTrigger) then
    write(*,*) WarningString('number of warning messages hits the limit. Further warnings are suppressed',this%moduleName)
  end if
end if

end subroutine WarningRaise_def

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Other functions   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!the function that makes string-ansi-colored
!!!initial code copied from http://fortranwiki.org/fortran/show/ansi_colors
pure function color(str, code) result(out)
  character(len=*), intent(in) :: str
  character(len=*), intent(in) :: code
  character(len=:), allocatable :: out
  out = c_start // code // c_end // str // c_clear
end function color

!!! move the CURRET in stream to the next line that starts from pos (5 char)
!!! this function is universally used in all reading of constants-files
subroutine MoveTO(stream,pos)
  integer,intent(in)::stream
  character(len=5),intent(in)::pos
  character(len=300)::line
  integer::IOstatus
  integer:: i

  do i=1,10000000
      read(stream,'(A)',IOSTAT=IOstatus) line
      if(IOstatus>0) then
          write(*,*) ErrorString("Error in attemt to read the line ("//pos//")", "aTMDe_IO_system")
          error stop
      else if(IOstatus<0) then
          write(*,*) ErrorString("EndOfFile during search of the line ("//pos//")", "aTMDe_IO_system")
          error stop
      else
          if(line(1:5)==pos) exit
      end if
  end do

  if(i>10000000) then
    write(*,*) ErrorString("Line ("//pos//") not found within iteration limit","aTMDe_IO_system")
    error stop
  end if

end subroutine MoveTO

!--------------------convertation
!!! convert a real(dp) number to a string
pure  function dpToStr(num)
real(dp),intent(in)::num
character(len=16)::dpToStr
write(dpToStr,"(G16.10)") num
end function dpToStr

!!! convert a real number to a string
pure function spToStr(num)
real(sp),intent(in)::num
character(len=12)::spToStr
write(spToStr,"(G12.6)") num
end function spToStr

!!! convert an integer number to a string
pure function intToStr(num)
integer,intent(in)::num
character(len=8)::intToStr
write(intToStr,"(I8)") num
end function intToStr
 
!! convert an short integer number to a string
pure function int4ToStr(num)
 integer,intent(in)::num
 character(len=4)::int4ToStr
 write(int4ToStr,"(I4)") num 
end function int4ToStr
 
!!!Common format of Warning line in artemide
pure function WarningString(str, moduleName) result(out)
  character(len=*), intent(in) :: str
  character(len=*), intent(in) :: moduleName
  character(len=:), allocatable :: out
  out = color('WARNING: artemide.'//trim(moduleName)//': '//trim(str),c_red)
end function WarningString

!!!Common format of error line in artemide
pure function ErrorString1(str, moduleName) result(out)
  character(len=*), intent(in) :: str
  character(len=*), intent(in) :: moduleName
  character(len=:), allocatable :: out
  out = color('ERROR: artemide.'//trim(moduleName)//': '//trim(str),c_red_bold)
end function ErrorString1

!!!Common format of error line in artemide
pure function ErrorString2(str, moduleName1, moduleName2) result(out)
  character(len=*), intent(in) :: str
  character(len=*), intent(in) :: moduleName1,moduleName2
  character(len=:), allocatable :: out
  out = color('ERROR: artemide.'//trim(moduleName1)//'.'//trim(moduleName2)//': '//trim(str),c_red_bold)
end function ErrorString2
  
end module aTMDe_IO
