program example
use aTMDe_control
use uTMDPDF
use uTMDFF
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),mean(:),deviation(:),b(:)
real*8::bMax,step,mu,dd,x,TMD(-5:5)
integer::f,k

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call uTMDPDF_SetLambdaNP(NParray(3:9),.false.,.false.)
call uTMDFF_SetLambdaNP(NParray(10:13),.false.,.false.)



f=1

numB=30
bMax=30d0
mu=4d0

allocate(b(1:numB))
do i=1,numB
  b(i)=i/5d0
end do

do j=0,18

x=1d0-j/20d0

do i=1,numB
    TMD=uTMDPDF_lowScale5(x,b(i),1)
    !TMD=uTMDFF_lowScale5(x,b(i),1)
    if(x>=0.99) then 
        write(*,"('{',F10.7,',',F10.7,',',F14.7,'},')") x,b(i),0d0
    else
        write(*,"('{',F10.7,',',F10.7,',',F14.7,'},')") x,b(i),TMD(f)
    end if
end do

end do

stop
! do j=1,18
! 
! x=0.1d0*(1d0-j/20d0)
! 
! do i=1,numB
!     TMD=uTMDPDF_lowScale5(x,b(i),1)
!     if(x>=0.99) then 
!         write(*,"('{',F10.7,',',F10.7,',',F14.7,'},')") x,b(i),0d0
!     else
!         write(*,"('{',F10.7,',',F10.7,',',F14.7,'},')") x,b(i),TMD(f)
!     end if
! end do
! 
! end do
! stop


x=0.75d0

allocate(central(1:numB))
allocate(mean(1:numB))
allocate(deviation(1:numB))

do i=1,numB
    TMD=uTMDPDF_lowScale5(x,b(i),1)
    !TMD=uTMDFF_lowScale5(x,b(i),1)
    central(i)=TMD(f)
    mean(i)=0d0
    deviation=0d0
end do

do j=1,numR
    call artemide_GetReplicaFromFile(repFILE,j,NParray)
    call uTMDPDF_SetLambdaNP(NParray(3:9),.false.,.false.)
    call uTMDFF_SetLambdaNP(NParray(10:13),.false.,.false.)
    do i=1,numB
        TMD=uTMDPDF_lowScale5(x,b(i),1)
        !TMD=uTMDFF_lowScale5(x,b(i),1)
        mean(i)=mean(i)+TMD(f)
        deviation(i)=deviation(i)+TMD(f)**2
    end do  

end do

do i=1,numB
    mean(i)=mean(i)/numR
    deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
end do


do i=1,numB
write(*,"('{',F10.7,',',F14.7,',',F14.7,',',F14.7,'},')") b(i),mean(i),(mean(i)-deviation(i)),(mean(i)+deviation(i))
end do

end program example
