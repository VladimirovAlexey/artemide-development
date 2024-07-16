program example
use aTMDe_control
use TMDX_DY
implicit none

integer,parameter::NUM=10

integer::i


real*8,dimension(1:NUM,1:2)::qT,Q,y
integer,dimension(1:NUM,1:4)::proc
logical,dimension(1:NUM)::iC
real*8,dimension(1:NUM,1:4)::cuts

real*8,dimension(1:NUM)::X,s
real*8,dimension(1:NUM)::qT_b,Q_b,y_b

real*8::time1,time2
!$ real*8::omp_get_wtime


!!!!! SETUP ARTEMIDE and NP-parameters

  call artemide_Initialize('xSec.atmde',prefix='Prog/MiniCourse/INI/')

  call artemide_SetNPparameters_TMDR((/1.5004d0, 0.05614d0, 0.03862d0, 0.0d0/))

  call artemide_SetNPparameters_uTMDPDF(&
  (/0.565d0, 0.0539d0, 0.5697d0, 6.64d0, 0.565d0, 20.07d0, 0.5697d0, 0.537d0, 1.07d0, 2.39d0, 0.0d0, 0.0d0/))


!!!!
!!!! set up parameters of cross-section
!!!!
do i=1,NUM

  Q(i,1:2)=(/81.d0,101.d0/)
  Q_b(i)=(Q(i,1)+Q(i,2))/2

  s(i)=8000.d0**2 !! 8TeV

  y(i,1:2)=(/-2.1d0,2.1d0/)
  y_b(i)=(y(i,1)+y(i,2))/2

  qt(i,1:2)=(/(i-1)*1.d0,i*1.d0/)
  qt_b(i)=(qt(i,1)+qt(i,2))/2

  proc(i,1:4)=(/1,1,1,3/)
  iC(i)=.false.
  cuts(i,1:4)=(/20.d0,20.d0,-2.1d0,2.1d0/)
end do

! !!!! Second example
! do i=1,NUM
!
!   Q(i,1:2)=(/6.d0,8.d0/)
!   Q_b(i)=(Q(i,1)+Q(i,2))/2
!
!   s(i)=2*800*0.938d0 !!! E772
!
!   y(i,1:2)=(/0d0,0.8d0/)
!   y_b(i)=(y(i,1)+y(i,2))/2
!
!   qt(i,1:2)=(/(i-1)*0.25d0,i*0.25d0/)
!   qt_b(i)=(qt(i,1)+qt(i,2))/2
!
!   proc(i,1:4)=(/1,1,1,3/)
!   iC(i)=.false.
!   cuts(i,1:4)=(/20.d0,15.d0,-1.5d0,1.5d0/)
! end do

!!!!
!!!! Compute point by point
!!!!

call cpu_time(time1)
!$ time1=omp_get_wtime()

do i=1,NUM
call xSec_DY(X(i),proc(i,:),s(i),qT(i,:),Q(i,:),y(i,:),iC(i),Cuts(i,:))
end do
call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=1,NUM
  write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT_b(i),X(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME WITH Point-by-Point:", time2-time1

write(*,*) "---------------------------------------------------------"

call cpu_time(time1)
!$ time1=omp_get_wtime()

call xSec_DY_List(X,proc,s,qT,Q,y,iC,Cuts)

call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=1,NUM
  write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT_b(i),X(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME WITH LIST:", time2-time1

write(*,*) "---------------------------------------------------------"

call cpu_time(time1)
!$ time1=omp_get_wtime()

call xSec_DY_List_APPROXIMATE(X,proc,s,qT,Q,y,iC,Cuts)

call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=1,NUM
  write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT_b(i),X(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME WITH LIST_APPROXIMATE:", time2-time1

write(*,*) "---------------------------------------------------------"

call cpu_time(time1)
!$ time1=omp_get_wtime()

call xSec_DY_List_BINLESS(X,proc,s,qT_b,Q_b,y_b,iC,Cuts)

call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=1,NUM
  write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT_b(i),(y(i,2)-y(i,1))*(Q(i,2)**2-Q(i,1)**2)*(qT(i,2)**2-qT(i,1)**2)*X(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME WITH LIST_BINLESS:", time2-time1

write(*,*) "---------------------------------------------------------"
end program example
