!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           arTeMiDe 3.01
!
! The module that handle the PDF from the LHAPDF-file. 1412.7420
! Reads info-file, interpolate, etc.
! This module is abstract and for several input-PDF one should make several instances of it.
!
!           A.Vladimirov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module LHAPDFreader
character (len=9),parameter :: moduleName="LHAreader"

use IO_functions
implicit none

!Current version of module
character (len=5),parameter :: version="v3.01"
character (len=9),parameter :: moduleName="LHAreader"
integer::outputLevel=2

!!!!! global variables to 
character(len=:), allocatable::ReferenceString !!! contains description of reference
character(len=:), allocatable::PDFname          !!! name of PDF for messages
character(len=:), allocatable::MainPath         !!! the path to the PDF-set without ".000?"
integer:: NumMembers
integer,allocatable,dimension(:):: FlavorArray
real(dp)::XMin,XMax,QMin,QMax

real(dp)::tolerance=10.d-6

!!!!!! Main grid (x,Q,f)
real(dp),dimension(:,:,:),allocatable::MainGrid
real(dp),dimension(:),allocatable::Xnodes,Qnodes,PDF_logQs,PDF_logXs
real(dp),allocatable,dimension(:,:):: PDF_WX,PDF_WQ!! barycentric wieghts
integer,allocatable,dimension(:,:):: PDF_indices!! indices

logical::TablesAreReady=.false. !!! flag that shows that module is ready

public:: ReadInfo

contains

!!! split a line by ":"
!!! attempt to identify the entry, and parse it
subroutine ParseInfoLine(line)
    character(len=*), intent(in)::line
    
    character(len=:), allocatable::linePart1, linePart2, linePart3
    integer::position_i,position_j
    integer::i,j
    
    position_i=SCAN(line,":")
    !!! First part of the line contrain the name
    linePart1=trim(line(:position_i-1))
    !!! second part of the line contains the value
    linePart2=trim(line(position_i+1:))
    
    SELECT CASE(linePart1)
    
    CASE("SetDesc")
        !!!!! Set Description
        if(outputLevel>2) then
            write(*,*) "Description of LHAPDF set:"
            write(*,*) linePart2
        end if
        
    CASE("SetIndex","Authors","FlavorScheme","OrderQCD","ErrorType","DataVersion", &
        "Format", "NumFlavors", "AlphaS_OrderQCD", "MZ","MUp","MDown","MStrange",  &
        "MCharm", "MBottom", "MTop", "AlphaS_MZ", "AlphaS_Type","AlphaS_Qs", "AlphaS_Vals" &
        )
        !!!!! Unused values of LHA-PDF file
        if(outputLevel>3) then
            write(*,*) "Value of "//trim(linePart1)//" in LHAPDF-info :"
            write(*,*) linePart2
        end if
        
    CASE("Reference")
        ReferenceString=linePart2
        if(outputLevel>2) then
            write(*,*) "Reference :",ReferenceString
        end if
    CASE("NumMembers")
        read(linePart2,*) NumMembers
        if(outputLevel>2) then
            write(*,*) "Number of Members of PDF-set :", NumMembers
        end if
    CASE("Flavors")
        !!!! flavors have structure [,,,]
        position_i=SCAN(linePart2,"[")
        position_j=SCAN(linePart2,"]")
        linePart3=linePart2(position_i+1:position_j-1)
        !!! determine lengh of array by counting ","
        j=0
        do i=0,len(linePart3)
            !write(*,*)linePart3(i:i)
            if(trim(linePart3(i:i)) .eq. ",") j=j+1            
        end do
        allocate(FlavorArray(1:j+1))
        read(linePart3,*) FlavorArray
        if(outputLevel>2) then
            write(*,'(A, I4, 999(", ",I4))') "Flavor-enumeration Array :", FlavorArray
        end if
        
    CASE("XMin")
        read(linePart2,*) XMin
        if(outputLevel>2) then
            write(*,'(A,F10.8)') " XMin     :", XMin
        end if
    CASE("XMax")
        read(linePart2,*) XMax
        if(outputLevel>2) then
            write(*,'(A,F10.8)') " XMax     :", XMax
        end if
    CASE("QMin")
        read(linePart2,*) QMin
        if(outputLevel>2) then
            write(*,'(A,F10.2)') " QMin     :", QMin
        end if
    CASE("QMax")
        read(linePart2,*) QMax
        if(outputLevel>2) then
            write(*,'(A,F10.2)') " QMax     :", QMax
        end if
    
    CASE DEFAULT
        if(outputLevel>1) then
        write(*,*) color("Unknown entry in LHAPDF-info file :"//linePart1,c_red)//" ... ignore"
        end if
    
    END SELECT
    
end subroutine ParseInfoLine

!!! opens the info-file and read lines.
subroutine ReadInfo(name,directory)
    character(len=*),intent(in)::name
    character(len=*),intent(in)::directory
    character(len=300)::path
    character(len=1024)::line !!!! 1024 line size!
    integer::ios,i,j
    
    TablesAreReady=.false.

    PDFname=trim(name)

    MainPath=trim(adjustl(directory))//trim(adjustr(name))//"/"//trim(adjustr(name))
    path=trim(MainPath)//".info"

    if(outputLevel>1) write(*,'(A)') color("----- Loading PDF from "//trim(name),c_yellow)

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old", IOSTAT=ios)
    if(ios /= 0) ERROR STOP ErrorString('The info-file is not found at '//trim(path),moduleName)
    
    !!! reading all lines
    do
        read(51,'(A)', iostat=ios) line
        if(ios /= 0) exit !!! end of file
        
        !!!! there could be empty lines.
        if(len(trim(line))>0) call ParseInfoLine(line)
    end do
    
    CLOSE (51, STATUS='KEEP')

    if(outputLevel>1) write(*,*) "Info-file for "//trim(PDFname)//" parsed succesfully"

    call SetReplica(0)

    if(outputLevel>1) write(*,'(A)') color("----- PDF "//trim(PDFname)//" loaded. ",c_yellow)&
        //color("["//trim(ReferenceString)//"]",c_blue)
    
end subroutine ReadInfo

!!! read the table for PDF-replica, and prepare all the Tables
!!! requires only the num. The PDF is searched in the MainPath, which is set by info-file
subroutine SetReplica(num)
    integer,intent(in)::num !!! the number of replica to set
    character(len=500)::path
    character(len=4096)::line !!!! 1024 line size!
    character(len=:),allocatable:: lineToParse
    real(dp),dimension(:),allocatable::listToCheck,PDFentry,Qentry
    integer,dimension(:),allocatable::listOfF
    integer,dimension(-5:5)::FlavorIndex !!! list that contain the enumeration change variable
    logical::XisSet
    integer::Xsize,Qsize
    integer::ios,i,j,k,iX,iQ,iF,indexQ

    TablesAreReady=.false.

    if(allocated(Xnodes)) deallocate(Xnodes)
    if(allocated(Qnodes)) deallocate(Qnodes)
    if(allocated(PDF_logQs)) deallocate(PDF_logQs)
    if(allocated(PDF_logXs)) deallocate(PDF_logXs)
    if(allocated(MainGrid)) deallocate(MainGrid)
    if(allocated(PDF_WQ)) deallocate(PDF_WQ)
    if(allocated(PDF_WX)) deallocate(PDF_WX)
    if(allocated(PDF_indices)) deallocate(PDF_indices)

    XisSet=.false.!!!! trigger fro the first set of X's
    Qsize=0
    allocate(listOfF(0:size(FlavorArray)-1))

    if(num<0 .or. num>NumMembers) &
        ERROR STOP ErrorString('Attempt to call for replica '//trim(intToStr(num))//&
            '. The maximum number of replicas (according to info) is '//trim(intToStr(NumMembers)),moduleName)

    !!!! create the name of replica-file MainPath_000n.dat
    write(path,'(A,"_",I4.4,".dat")') trim(MainPath),num

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old", IOSTAT=ios)
    if(ios /= 0) ERROR STOP ErrorString('The file is not found at '//trim(path),moduleName)

    !!!! the file should be structure as follows [1412.7420]
    !!!PdfType: centra
    !!!Format: lhagrid1
    !!!---
    !!! List of X
    !!! List of Q
    !!! List of flavors
    !!! values per flavor for x1,Q1
    !!! values per flavor for x1,Q2
    !!! ...
    !!! values per flavor for x1,Qn
    !!! values per flavor for x2,Q1
    !!! ...
    !!! ...
    !!! values per flavor for xn,Qn
    !!!---
    !!! List of X
    !!! List of Q
    !!! List of flavors
    !!! values per flavor for x1,Q1
    !!! etc.
    !!!---
    !!!

    !!! To read such structure I do two runs of reading

    !!!! the first run is the check of the consistency
    !!!! 1) all lists of flavors are the same
    !!!! 2) all lists of X are the same
    !!!! 3) Count the number of Q's

    !!! then I can allocate the table for PDF, and other tables
    !!! The second run
    !!! 1) Fill the list of Q's
    !!! 2) Fill the list of PDF's

    !!! reading all lines (first run)
    do
        read(51,'(A)', iostat=ios) line
        if(ios /= 0) exit !!! end of file
        !write(*,*) "---->>>   ",trim(line)

        !!!! search for "---"
        if(trim(line)=="---") then
            read(51,'(A)', iostat=ios) line
            if(ios /= 0) exit !!! end of file

            !!!!!!!!-------------This should be the line of X's
            !!!!!!!! parse and check the correctness
            lineToParse=trim(line)

            !!! the list is space separated (count the number of spaces)
            j=0
            do i=1,len(lineToParse)
                if(lineToParse(i:i) .eq. " ") j=j+1
            end do
            if(XisSet) then
                if(allocated(listToCheck)) deallocate(listToCheck)
                allocate(listToCheck(0:j))
                read(lineToParse,*) listToCheck
                if(size(listToCheck)/=Xsize) then
                    CLOSE (51, STATUS='KEEP')
                    ERROR STOP ErrorString("Different sizes of x-grids in LHAPDF file.",moduleName)
                end if
                if(sum(abs(Xnodes/listToCheck)-1)>tolerance) then
                    CLOSE (51, STATUS='KEEP')
                    ERROR STOP ErrorString("X-grids for subgrids in LHAPDF file do not coincide.",moduleName)
                end if
            else
                XisSet=.true.
                !!! set up the first list of X's
                allocate(Xnodes(0:j))
                Xsize=j+1
                read(lineToParse,*) Xnodes
            end if

            !!!!!!!!-------------This should be the line of Q's
            !!!!!!!! count number of entries
            read(51,'(A)', iostat=ios) line
            lineToParse=trim(line)
            !!! the list is space separated (count the number of spaces)
            j=1 !! (starts from 1 because it number of spaces +1)
            do i=1,len(lineToParse)
                if(lineToParse(i:i) .eq. " ") j=j+1
            end do
            Qsize=Qsize+j

            !!!!!!!!-------------This should be the line of flavors
            !!!!!!!! check that it coincides with the line of flavors from the info-file
            read(51,'(A)', iostat=ios) line
            lineToParse=trim(line)
            read(lineToParse,*) listOfF
            if(sum(listOfF-FlavorArray)/=0) then
                CLOSE (51, STATUS='KEEP')
                ERROR STOP ErrorString("Flavor-grids for subgrids in LHAPDF file do not coincide.",moduleName)
            end if

            !!! then it should be j * Xsize lines of grids
            do i=1, j*Xsize
                read(51,'(A)', iostat=ios) line
                if(ios /= 0) then
                    CLOSE (51, STATUS='KEEP')
                    ERROR STOP ErrorString("Unexpected END-OF-FILE in "//trim(path),moduleName)
                end if

                if(trim(line)=="---") then
                    CLOSE (51, STATUS='KEEP')
                    ERROR STOP ErrorString("Incorrect number of entries subgrid of "//trim(path),moduleName)
                end if
            end do

        end if
    end do

    CLOSE (51, STATUS='KEEP')

    !!!!! now the Qsize, Xsize, and Xnodes are set up

    !!! Allocate varriables
    allocate(Qnodes(0:Qsize-1))
    allocate(MainGrid(0:Xsize-1,0:Qsize-1,-5:5))
    allocate(PDFentry(1:size(FlavorArray)))

    !!!! the list of parcing the flavor to (-5:5)
    k=0
    do i=-5,5
        do j=1,size(FlavorArray)
            if(i==0 .and. FlavorArray(j)==21) then
                FlavorIndex(i)=j
                k=k+1
            else if(i/=0 .and. i>=-5 .and. i<=5 .and. FlavorArray(j)==i) then
                FlavorIndex(i)=j
                k=k+1
            end if
        end do
    end do
    if(k/=11) then
        ERROR STOP ErrorString("The flavor numbering list does not map to (-5:5)",moduleName)
    end if

    !!! now read the second time and write all grids
    !k=1!!!! counter of PDF-entry
    !j=1!!!! counter of Q-entry
    indexQ=0 !!!! index of current position in Qnodes

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old", IOSTAT=ios)
    do
        read(51,'(A)', iostat=ios) line
        if(ios /= 0) exit !!! end of file
        !write(*,*) "---->>>   ",trim(line)

        !!!! search for "---"
        if(trim(line)=="---") then
            read(51,'(A)', iostat=ios) line !!!! this is X-line (skip) or EOF
            !write(*,*) "---->>>",trim(line)
            if(ios /= 0) exit !!! end of file


            !!!!!!!!-------------This should be the line of Q's
            !!!!!!!! re-read entries and add to the Qline
            !!!!!!!! count number of entries
            read(51,'(A)', iostat=ios) line
            lineToParse=trim(line)
            !!! the list is space separated (count the number of spaces)
            j=0
            do i=1,len(lineToParse)
                if(lineToParse(i:i) .eq. " ") j=j+1
            end do
            allocate(Qentry(0:j))
            read(lineToParse,*) Qentry
            !!! add it to the list of Qnodes
            Qnodes(indexQ:indexQ+j)=Qentry(0:j)

            !!!!!!!------------This is line of flavors (skip)
            read(51,'(A)', iostat=ios) line

            !!!!!! reading the PDF-entries
            !!!!!! and insert it into the corresponding Grid-elements
            do i=0,Xsize-1
            do k=0,size(Qentry)-1
                read(51,*) PDFentry
                !!!write(*,*) i,indexQ+k
                do j=-5,5
                    MainGrid(i,indexQ+k,j)=PDFentry(FlavorIndex(j))
                end do
            end do
            end do

            deallocate(Qentry)
            indexQ=indexQ+size(Qentry)  !!!! not forget to update indexQ
        end if

    end do
    CLOSE (51, STATUS='KEEP')

    !!!!!------------ prepare the table of interpolation
    allocate(PDF_logQs(0:size(Qnodes)-1))
    allocate(PDF_logXs(0:size(Xnodes)-1))
    allocate(PDF_indices(1:4,0:size(Qnodes)-2))
    allocate(PDF_WQ(1:4,0:size(Qnodes)-2))
    allocate(PDF_WX(1:4,0:size(Xnodes)-2))

    PDF_LogQs=log(Qnodes)
    PDF_LogXs=log(Xnodes)

    !!! compute the index table for cubic-interpolation
    !!! if the x belongs to (x_{j-1}<x<x_j) the interpolation is done by indices{j}
    do i=0,size(PDF_LogQs)-2
        if(i==0) then
            PDF_indices(1:4,i)=(/0,1,2,3/)
        else if(i==size(PDF_LogQs)-2) then
            PDF_indices(1:4,i)=i-1-(/1,0,-1,-2/)
        else
            PDF_indices(1:4,i)=i-(/1,0,-1,-2/)
            !!! if the value of alpha is repeated these are ends of subgrids.
            !!! I shift the indicing
            if(abs(PDF_LogQs(i+1)-PDF_LogQs(i+2))<10d-4) then
                PDF_indices(1:4,i)=PDF_indices(1:4,i)-1
            else if(abs(PDF_LogQs(i)-PDF_LogQs(i+1))<10d-4) then
                !!! zero-size bin
                !!! this should be never reached. but I assign it to previous one
                PDF_indices(1:4,i)=PDF_indices(1:4,i)-2
            else if(abs(PDF_LogQs(i-1)-PDF_LogQs(i))<10d-4) then
                PDF_indices(1:4,i)=PDF_indices(1:4,i)+1
            end if
        end if
    end do

    !!!! table for X is very simple,  there is no repeating entries.
    !!!! so index is i+(-1,0,1,2), except for i=0 -> (0,1,2,3), and i=S-1 (S-3,S-2,S-1,S)

    !!! precompute barycentric weights for Q
    do i=0,size(PDF_LogQs)-2
        PDF_WQ(1,i)=1.d0/(&
        (PDF_LogQs(PDF_indices(1,i))-PDF_LogQs(PDF_indices(2,i)))*&
        (PDF_LogQs(PDF_indices(1,i))-PDF_LogQs(PDF_indices(3,i)))*&
        (PDF_LogQs(PDF_indices(1,i))-PDF_LogQs(PDF_indices(4,i))))

        PDF_WQ(2,i)=1.d0/(&
        (PDF_LogQs(PDF_indices(2,i))-PDF_LogQs(PDF_indices(1,i)))*&
        (PDF_LogQs(PDF_indices(2,i))-PDF_LogQs(PDF_indices(3,i)))*&
        (PDF_LogQs(PDF_indices(2,i))-PDF_LogQs(PDF_indices(4,i))))

        PDF_WQ(3,i)=1.d0/(&
        (PDF_LogQs(PDF_indices(3,i))-PDF_LogQs(PDF_indices(1,i)))*&
        (PDF_LogQs(PDF_indices(3,i))-PDF_LogQs(PDF_indices(2,i)))*&
        (PDF_LogQs(PDF_indices(3,i))-PDF_LogQs(PDF_indices(4,i))))

        PDF_WQ(4,i)=1.d0/(&
        (PDF_LogQs(PDF_indices(4,i))-PDF_LogQs(PDF_indices(1,i)))*&
        (PDF_LogQs(PDF_indices(4,i))-PDF_LogQs(PDF_indices(2,i)))*&
        (PDF_LogQs(PDF_indices(4,i))-PDF_LogQs(PDF_indices(3,i))))
    end do

    !!! precompute barycentric weights for X
    do i=0,size(PDF_LogXs)-2
        if(i==0) then
            k=1
        else if(i==size(PDF_LogXs)-2) then
            k=i-1
        else
            k=i
        end if

        PDF_WX(1,i)=1.d0/(&
        (PDF_LogXs(k-1)-PDF_LogXs(k))*(PDF_LogXs(k-1)-PDF_LogXs(k+1))*(PDF_LogXs(k-1)-PDF_LogXs(k+2)))

        PDF_WX(2,i)=1.d0/(&
        (PDF_LogXs(k)-PDF_LogXs(k-1))*(PDF_LogXs(k)-PDF_LogXs(k+1))*(PDF_LogXs(k)-PDF_LogXs(k+2)))

        PDF_WX(3,i)=1.d0/(&
        (PDF_LogXs(k+1)-PDF_LogXs(k-1))*(PDF_LogXs(k+1)-PDF_LogXs(k))*(PDF_LogXs(k+1)-PDF_LogXs(k+2)))

        PDF_WX(4,i)=1.d0/(&
        (PDF_LogXs(k+2)-PDF_LogXs(k-1))*(PDF_LogXs(k+2)-PDF_LogXs(k))*(PDF_LogXs(k+2)-PDF_LogXs(k+1)))
    end do

    if(outputLevel>1) write(*,"(A,I4,A)") " "//trim(moduleName)//": replica ",num," of "//trim(PDFname)//" loaded."


end subroutine SetReplica

!!! returns the value of xPDF extrapolated from the grid
function xPDF(x,Q)
real(dp),intent(in)::x,Q
real(dp),dimension(-5:5):: xPDF

real(dp)::logQ,dd,deltaQ(1:4),subQ(1:4,-5:5),logX,deltaX(1:4)
integer::iQ,iX,j
logical::flag

if(Q<Qmin .or. Q>Qmax) ERROR STOP ErrorString('Q is outside of Q-range',moduleName)
if(X<Xmin .or. X>Xmax) ERROR STOP ErrorString('X is outside of X-range',moduleName)

logQ=log(Q)
logX=log(X)

!!! search for the index of Q
do iQ=0,size(PDF_LogQs)-2
     if(logQ<PDF_LogQs(iQ+1)) exit
end do
!!! search for the index of X
if(logX<PDF_LogXs(1)) then
    iX=1
else if(logX>PDF_LogXs(size(PDF_LogXs)-2)) then
    iX=size(PDF_LogXs)-3
else
    do iX=1,size(PDF_LogXs)-3
        if(logX<PDF_LogXs(iX+1)) exit
    end do
end if

!!!! Q-interpolation
do j=1,4
    dd=logQ-PDF_LogQs(PDF_indices(j,iQ))
    !!! I check if the point is close to the node
    if(abs(dd)<tolerance) then
        subQ(1:4,-5:5)=MainGrid(iX-1:iX+2,PDF_indices(j,iQ),-5:5)
        goto 10
    end if
    deltaQ(j)=PDF_WQ(j,iQ)/dd
end do

subQ=(deltaQ(1)*MainGrid(iX-1:iX+2,PDF_indices(1,iQ),-5:5)+&
deltaQ(2)*MainGrid(iX-1:iX+2,PDF_indices(2,iQ),-5:5)+&
deltaQ(3)*MainGrid(iX-1:iX+2,PDF_indices(3,iQ),-5:5)+&
deltaQ(4)*MainGrid(iX-1:iX+2,PDF_indices(4,iQ),-5:5))/sum(deltaQ)

!!! X-interpolation
10 do j=1,4
    dd=logX-PDF_LogXs(iX+j-2)
    !!! I check if the point is close to the node
    if(abs(dd)<tolerance) then
        xPDF(-5:5)=subQ(j,-5:5)
        return
    end if
    deltaX(j)=PDF_WX(j,iX)/dd
end do

xPDF=(deltaX(1)*subQ(1,-5:5)+deltaX(2)*subQ(2,-5:5)+deltaX(3)*subQ(3,-5:5)+deltaX(4)*subQ(4,-5:5))/sum(deltaX)

end function xPDF

end module LHAPDFreader

program example
use LHAPDFreader
implicit none

real*8::aa,QQ,xx
real*8::ff(-5:5)
integer::i

call ReadInfo("MSHT20_REP","Prog/LHAPDFreader/LHAexample/")

call SetReplica(5)

do i=0,120
xx=10**(-i/20.d0)
QQ=i/3.+1.2d0
ff=xPDF(xx,QQ)
write(*,'("{",F12.9,",",F12.8,"},")') xx,ff(1)
end do

!do i=0,60
!QQ=10**(i/30.d0)
!ff=xPDF(0.11d0,QQ)
!write(*,'("{",F10.2,",",F12.8,"},")') QQ,ff(1)
!end do

!ff=xPDF(0.000001d0,1.20000001d0)
!write(*,*) ff

! ff=xPDF(0.0000012d0,1.20000001d0)
! write(*,*) ff
!
! ff=xPDF(0.0000015d0,1.20000001d0)
! write(*,*) ff

end program example
