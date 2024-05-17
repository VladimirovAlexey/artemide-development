program example
use aTMDe_control
use uTMDPDF
use TMDX_DY
implicit none

integer,parameter::NUM=40

integer::i

real*8,dimension(1:NUM)::X,s,qT,Q,y
integer,dimension(1:NUM,1:4)::proc
logical,dimension(1:NUM)::iC
real*8,dimension(1:NUM,1:4)::cuts
real*8,dimension(-5:5)::kkk
real*8::time1,time2
!$ real*8::omp_get_wtime

  call artemide_Initialize('KPC2.atmde',prefix='Prog/KPCtest/INI/')
  !call artemide_Initialize('LP2.atmde',prefix='Prog/KPCtest/INI/')

  call artemide_SetNPparameters_TMDR((/1.5004d0, 0.05614d0, 0.03862d0, 0.0d0/))

  call artemide_SetNPparameters_uTMDPDF(&
  (/0.565d0, 0.0539d0, 0.5697d0, 6.64d0, 0.565d0, 20.07d0, 0.5697d0, 0.537d0, 1.07d0, 2.39d0, 0.0d0, 0.0d0/))
!!!! set up variables
do i=1,NUM
  Q(i)=4.5d0
  !Q(i)=1.5d0+(i-1)*1
  !s(i)=3841600.d0
  s(i)=1500.d0
  y(i)=0.5d0
  !y(i)=-6.*(real(i)/NUM-0.5)
  qt(i)=0.1d0+(i-1)*0.1
  !qt(i)=10.d0
  proc(i,1:4)=(/2,1,1,20/)  !!KPC
  !proc(i,1:4)=(/1,1,1,3/)   !! LP
  iC(i)=.false.
  cuts(i,1:4)=(/0d0,0d0,-100d0,100d0/)
  !kkk=uTMDPDF_5(0.04d0,qt(i),Q(i),Q(i)*2,1)
  !X(i)=kkk(1)
end do

call cpu_time(time1)
!$ time1=omp_get_wtime()

call xSec_DY_List_BINLESS(X,proc,s,qT,Q,y,iC,cuts)

call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=1,NUM
  !write(*,'("{",F12.8,",",F16.10,"},")') Q(i),4*qT(i)*Q(i)*X(i)
  !write(*,'("{",F12.8,",",F16.10,"},")',advance="no") y(i),4*qT(i)*Q(i)*X(i)
  write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT(i),X(i)!4*qT(i)*Q(i)*X(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME:", time2-time1

write(*,*) "---------------------------------------------------------"
stop
!!!! set up variables
do i=1,NUM
  Q(i)=50.d0
  !Q(i)=1.5d0+(i-1)*1
  s(i)=3841600.d0
  !s(i)=(100.*Q(i))**2
  y(i)=2.5d0
  !y(i)=-6.*(real(i)/NUM-0.5)
  qt(i)=0.1d0+(i-1)*0.5
  !qt(i)=10.d0
  proc(i,1:4)=(/1,1,1,2/)  !!KPC
  !proc(i,1:4)=(/1,1,1,3/)   !! LP
  iC(i)=.false.
  cuts(i,1:4)=(/0d0,0d0,-100d0,100d0/)
  !kkk=uTMDPDF_5(0.04d0,qt(i),Q(i),Q(i)*2,1)
  !X(i)=kkk(1)
end do

call cpu_time(time1)
!$ time1=omp_get_wtime()

call xSec_DY_List_BINLESS(X,proc,s,qT,Q,y,iC,cuts)

call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=1,NUM
  !write(*,'("{",F12.8,",",F16.10,"},")') Q(i),4*qT(i)*Q(i)*X(i)
  !write(*,'("{",F12.8,",",F16.10,"},")',advance="no") y(i),4*qT(i)*Q(i)*X(i)
  write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT(i),4*qT(i)*Q(i)*X(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME:", time2-time1

write(*,*) "---------------------------------------------------------"

!!!! set up variables
do i=1,NUM
  Q(i)=10.d0
  !Q(i)=1.5d0+(i-1)*1
  s(i)=(27.1d0)**2
  !s(i)=(100.*Q(i))**2
  y(i)=0.d0!-6.*(real(i)/NUM-0.5)
  qt(i)=0.1d0+(i-1)*0.05
  proc(i,1:4)=(/1,1,1,2/)  !!KPC
  !proc(i,1:4)=(/1,1,1,3/)   !! LP
  iC(i)=.false.
  cuts(i,1:4)=(/0d0,0d0,-100d0,100d0/)
  !kkk=uTMDPDF_inB(0.04d0,qt(i),Q(i),Q(i)*2,1)
  !(i)=kkk(1)
end do

call cpu_time(time1)
!$ time1=omp_get_wtime()

call xSec_DY_List_BINLESS(X,proc,s,qT,Q,y,iC,cuts)

call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=1,NUM
  !write(*,'("{",F12.8,",",F16.10,"},")') Q(i),4*qT(i)*Q(i)*X(i)
  write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT(i),4*qT(i)*Q(i)*X(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME:", time2-time1

end program example
