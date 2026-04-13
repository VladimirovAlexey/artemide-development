program example
use aTMDe_control
use aTMDe_numerics
use uTMDPDF
use TMDX_DY
implicit none

integer, parameter :: NUM = 25  ! defines the number of points to compute

integer :: i ! to do the loop

real*8, dimension(1:NUM) :: s
real*8, dimension(1:NUM,1:2) :: qT,y,Q
integer, dimension(1:NUM,1:4) :: proc
logical, dimension(1:NUM) :: iC
real*8, dimension(1:NUM,1:4) :: cuts

real*8, dimension(1:NUM) :: xx,xxKPC,xxLP

real*8 :: time1,time2

call artemide_Initialize('ART_KPC.atmde',prefix='Prog/GridCalc/')
!   Setting the non-perturbative parameters (NP) for TMD evolution
call artemide_SetNPparameters_TMDR((/1.5004d0, 0.05614d0, 0.03862d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse

!   Setting the non-perturbative parameters (NP) for the uTMDPDF
call artemide_SetNPparameters_uTMDPDF(&
(/0.565d0, 0.0539d0, 0.5697d0, 6.64d0, 0.565d0, 20.07d0, 0.5697d0, 0.537d0, 1.07d0, 2.39d0, 0.0d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse


do i=1,NUM

   Q(i,1) = 70.d0
   Q(i,2) = 110.d0
   qT(i,1) = 0.d0+0.5d0*(i-1)
   qT(i,2) = 0.d0+0.5d0*i

   s(i) = (1500.d0)**2 !! TeV

   y(i,1) = -1.d0
   y(i,2) = 1.d0

   proc(i,1:4)=(/1,1,1,3/) !KPCs

!    proc(i,1:4)=(/1,1,1,2001/) !LP

   iC(i)=.false.
   cuts(i,1:4)=(/20.d0,20.d0,-2.1d0,2.1d0/)

 end do


!  !!!!
!  !!!! Compute point by point
!  !!!!
call cpu_time(time1)
!
!  call xSec_SIDIS_BINLESS_List_forharpy(xx,proc,s,pT,z,x,Q,MM)
!

call xSec_DY_List(xxKPC,proc,s,qT,Q,y,iC,cuts)
call cpu_time(time2)


do i=1,NUM
!write(*,'("{",F12.6,",",F12.6,",",F12.6,"},")') Q(i),Q(i)**2*pi/2*xxKPC(i),xxLP(i)
write(*,'("{",F12.6,",",F12.6,"},")') (qT(i,1)+qT(i,2))/2,xxKPC(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME:", time2-time1

write(*,*) "---------------------------------------------------------"
end program example
