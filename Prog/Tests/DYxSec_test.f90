!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot for AS-term
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_DY
implicit none

integer,parameter::NUM=100

real*8::Q(0:NUM,1:2),x1(0:NUM,1:2),x2(0:NUM,1:2),qT(0:NUM,1:2),XX(0:NUM),s(0:NUM),y(0:NUM,1:2)
!real*8::Q(0:NUM),x1(0:NUM),x2(0:NUM),qT(0:NUM),XX(0:NUM),s(0:NUM),y(0:NUM)
logical,dimension(0:NUM)::includeCuts
integer::proc(0:NUM,1:4)
real*8,dimension(0:NUM,1:4)::CutParameters
integer::i
real*8::x,mu,time1,time2
!$ real*8::omp_get_wtime

call artemide_Initialize('xSec_DY.atmde',prefix='Prog/Tests/const-files/')
!   Setting the non-perturbative parameters (NP) for TMD evolution
call artemide_SetNPparameters_TMDR((/1.5004d0, 0.05614d0, 0.03862d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse

!   Setting the non-perturbative parameters (NP) for the uTMDPDF
call artemide_SetNPparameters_uTMDPDF(&
(/0.565d0, 0.0539d0, 0.5697d0, 6.64d0, 0.565d0, 20.07d0, 0.5697d0, 0.537d0, 1.07d0, 2.39d0, 0.0d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse

 do i=0,NUM

   Q(i,1:2) =(/75.d0,115.d0/)
   qT(i,1:2) = (/0.d0+i/2.5d0,0.d0+(i+1)/2.5d0/)
   y(i,1:2)=(/-1.d0,1.d0/)

! !
!    Q(i) =90.d0
!    qT(i) = 5.d0!0.d0+i/2.5d0
!    x1(i)=0.001d0+i*0.01d0!0.01d0
!    x2(i)=0.03d0
!    y(i)=0.d0
    s(i) = (1500.d0)**2

   proc(i,1:4)=(/1,1,1,3/)

   includeCuts(i)=.false.
   CutParameters(i,:)=(/5.d0,5.d0,-4.d0,4.d0/)
end do


!  !!!!
!  !!!! Compute point by point
!  !!!!
call cpu_time(time1)
!$ time1=omp_get_wtime()
!
!call xSec_SIDIS_BINLESS_List_forharpy(xx,proc,s,pT,z,x,Q,MM)
!
! !
! do i=0,NUM
!
!   XX(i)=F_toGrid(Q(i)**2,qT(i),x1(i),x2(i),proc(i,2:4))
! end do

!call xSec_DY_List_BINLESS(XX,proc,s,qT,Q,y,includeCuts,CutParameters)
call xSec_DY_List(XX,proc,s,qT,Q,y,includeCuts,CutParameters)

call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=0,NUM
write(*,'("{",F12.6,",",F16.6,"},")') (qT(i,2)+qT(i,1))/2,100000*XX(i)
!write(*,'("{",F12.6,",",F16.6,"},")') x1(i),XX(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME:", time2-time1

write(*,*) "---------------------------------------------------------"

end program example
