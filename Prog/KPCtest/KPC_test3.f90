program example
use aTMDe_control
use TMDX_DY
implicit none

integer,parameter::NUM=4

integer::i

real*8,dimension(1:NUM)::X,s
real*8,dimension(1:NUM,1:2)::qT,Q,y
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
  Q(i,1:2)=(/4.5d0,9.5d0/)
  !Q(i)=1.5d0+(i-1)*1
  s(i)=2*800*0.938d0
  !s(i)=(100.*Q(i))**2
  y(i,1:2)=(/0.d0,0.8d0/)
  !y(i)=-6.*(real(i)/NUM-0.5)
  qt(i,1:2)=(/(i-1)*0.5,i*0.5/)

  proc(i,1:4)=(/2,1,1,2/)  !!KPC
  !proc(i,1:4)=(/1,1,1,3/)   !! LP
  iC(i)=.false.
  cuts(i,1:4)=(/20.d0,15.d0,-1.5d0,1.5d0/)
end do

call cpu_time(time1)
!$ time1=omp_get_wtime()

call xSec_DY_List_APPROXIMATE(X,proc,s,qT,Q,y,iC,Cuts)

call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=1,NUM
  !write(*,'("{",F12.8,",",F16.10,"},")') Q(i),4*qT(i)*Q(i)*X(i)
  !write(*,'("{",F12.8,",",F16.10,"},")',advance="no") y(i),4*qT(i)*Q(i)*X(i)
  write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT(i,1),X(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME:", time2-time1

write(*,*) "---------------------------------------------------------"
stop

! !!!! set up variables
! do i=1,NUM
!   Q(i,1:2)=(/10.d0-2d0,10.d0+2d0/)
!   !Q(i)=1.5d0+(i-1)*1
!   s(i)=(32.5d0)**2
!   !s(i)=(100.*Q(i))**2
!   y(i,1:2)=(/-2.4d0,2.4d0/)
!   !y(i)=-6.*(real(i)/NUM-0.5)
!   qt(i,1:2)=(/0.01d0+(i-1)*0.05,0.1d0+i*0.05/)
!
!   !proc(i,1:4)=(/1,1,1,2/)  !!KPC
!   proc(i,1:4)=(/1,1,1,3/)   !! LP
!   iC(i)=.false.
!   cuts(i,1:4)=(/0d0,0d0,-100d0,100d0/)
! end do
!
! call cpu_time(time1)
! !$ time1=omp_get_wtime()
!
! call xSec_DY_List(X,proc,s,qT,Q,y,iC,Cuts)
!
! call cpu_time(time2)
! !$ time2=omp_get_wtime()
!
! do i=1,NUM
!   !write(*,'("{",F12.8,",",F16.10,"},")') Q(i),4*qT(i)*Q(i)*X(i)
!   !write(*,'("{",F12.8,",",F16.10,"},")',advance="no") y(i),4*qT(i)*Q(i)*X(i)
!   write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT(i,1)+0.25,X(i)
! end do
! write(*,*) " "
! write(*,*) " "
! write(*,*) " COMPUTATION TIME:", time2-time1
!
! write(*,*) "---------------------------------------------------------"

end program example
