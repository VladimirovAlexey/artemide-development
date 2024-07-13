module LHAPDFreader
use IO_functions
implicit none

public:: ReadInfo

contains

!!! opens the file and search the line that contains lineMark (in the begining of the string)
!!! returns the string (with mark removed)
function ReturnLine(file,lineMark)
character(len=*)::file
character(len=1024)::line !!!! 1024 line size!
OPEN(UNIT=51, FILE=file, ACTION="read", STATUS="old")

    
CLOSE (51, STATUS='KEEP')
end function ReturnLine

subroutine ReadInfo(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path


    if(present(prefix)) then
        path=trim(adjustl(prefix))//trim(adjustr(file))
    else
        path=trim(adjustr(file))
    end if

    

end subroutine ReadInfo

end LHAPDFreader

program example
use LHAPDFreader
implicit none

end program example
