program example
use aTMDe_control
use TMDX_DY
use TMDs
implicit none

integer::i,j,k

real*8::s
integer,dimension(1:4)::proc

character(len=2)::NAME

integer,parameter::NQ=19
real*8,dimension(1:NQ)::Q,yA,qTA,sA,X
integer,dimension(1:NQ,1:4)::procA
logical,dimension(1:NQ)::iC
real*8,dimension(1:NQ,1:4)::cuts

integer,parameter::Ny=29
real*8,dimension(1:Ny)::y

integer,parameter::NT=41
real*8,dimension(1:NT)::qT

Q=(/80.d0,82.d0,84.d0,86.d0,87.d0,88.d0,89.d0,90.d0,90.5d0,91.d0,91.25d0,91.5d0,92.d0,93.d0,94.d0,95.d0,96.d0,98.d0,100.d0/)

y=(/-3.5d0,-3.25d0,-3.d0,-2.75d0,-2.5d0,-2.25d0,-2.d0,-1.75d0,-1.5d0,-1.25d0,-1.d0,-0.75d0,-0.5d0,-0.25d0,&
0.d0,0.25d0,0.5d0,0.75d0,1.d0,1.25d0,1.5d0,1.75d0,2.d0,2.25d0,2.5d0,2.75d0,3.d0,3.25d0,3.5d0/)

qT=(/0.01d0,0.5d0,1.d0,1.5d0,2.d0,2.5d0,3.d0,3.5d0,4.d0,4.5d0,5.d0,5.5d0,6.d0,6.5d0,7.d0,7.5d0,8.d0,&
8.5d0,9.d0,9.5d0,10.d0,10.5d0,11.d0,110.5d0,12.d0,12.5d0,13.d0,13.5d0,14.d0,14.5d0,15.d0,15.5d0,16.d0,16.5d0,&
17.d0,17.5d0,18.d0,18.5d0,19.d0,19.5d0,20.d0/)


s=8000.d0**2

call artemide_Initialize('KPC2.atmde',prefix='Prog/KPCtest/INI/')

call artemide_SetNPparameters_TMDR((/1.56142d0, 0.0369174d0, 0.0581734d0, 1.0d0/))

call artemide_SetNPparameters_uTMDPDF(&
  (/0.874245d0, 0.913883d0, 0.991563d0, 6.05412d0, 0.353908d0,&
  46.6064d0, 0.115161d0, 1.53235d0, 1.31966d0, 0.434833d0, 0.d0, 0.d0/))


NAME="A4"
proc(1:4)=(/1,1,1,24/)  !!KPC


open (22, FILE='Prog/KPCtest/dat/'//NAME//'.dat', STATUS='REPLACE')
close(22)

do j=1,Ny
do k=1,NT
  do i=1,NQ
  sA(i)=s
  procA(i,1:4)=proc(1:4)
  qTA(i)=qT(k)
  yA(i)=y(j)

  iC(i)=.false.
  cuts(i,1:4)=(/0d0,0d0,-100d0,100d0/)

  end do

  call xSec_DY_List_BINLESS(X,procA,sA,qTA,Q,yA,iC,cuts)

  open (22, FILE='Prog/KPCtest/dat/'//NAME//'.dat', status="old", position="append", action="write")
    do i=1,NQ
      write(22,*) Q(i), yA(i), qTA(i), X(i)
    end do
  close(22)

  write(*,'(A," : ",I3,"/",I3,"  ",I3,"/",I3)') NAME,j,Ny,k,NT

end do
end do


end program example
