program example

! Invocation of required modules
use aTMDe_control
use aTMDe_numerics
use uTMDPDF
use TMDX_DY
use TMDF
use TMDF_KPC


implicit none


 integer, parameter :: NUM = 10  ! defines the number of points to compute

 integer :: i ! to do the loop

 real*8, dimension(1:NUM) :: s,qT,x1,x2,Q
 integer, dimension(1:NUM,1:4) :: proc
 logical, dimension(1:NUM) :: iC
 real*8, dimension(1:NUM,1:4) :: cuts

 real*8, dimension(1:NUM) :: xx,xxKPC,xxLP
 real*8, dimension(1:NUM,1:2) :: MM


 real*8 :: time1,time2, t1(1:NUM),t2(1:NUM),t1LP(1:NUM),t2LP(1:NUM)


 !!!!! SETUP ARTEMIDE and NP-parameters
 !
 !   Initialization of arTeMiDe using the constants-file xSec_SIDIS.atmde
   call artemide_Initialize('xSec_DY.atmde',prefix='Prog/KPCtest/INI/')
 !   call artemide_Initialize('KPC_my_test.atmde',prefix='Prog/KPCtest/INI/')
 !
 !   Setting the non-perturbative parameters (NP) for TMD evolution
   call artemide_SetNPparameters_TMDR((/1.5004d0, 0.05614d0, 0.03862d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse

 !   Setting the non-perturbative parameters (NP) for the uTMDPDF
   call artemide_SetNPparameters_uTMDPDF(&
   (/0.565d0, 0.0539d0, 0.5697d0, 6.64d0, 0.565d0, 20.07d0, 0.5697d0, 0.537d0, 1.07d0, 2.39d0, 0.0d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse


!(Q2, qT_in, x1, z1 ,mu, proc1)
  !xxKPC(1)=KPC_SIDISconv(1000.d0,0.001d0,0.5d0,0.5d0,10.d0,(/1,1,999/))
!    xxKPC(1)=KPC_DYconv(52._dp**2,0.1_dp,0.01d0,0.001d0,52._dp,(/1,1,1/))
!     stop
!   write(*,*) ">>>>",xxKPC(1)/0.5d0

  !stop

 do i=1,NUM

   Q(i) = 50.d0+8*(i-1)*2

   s(i) = 1500.d0 !! TeV

   x1(i) = 0.01d0
   x2(i) = 0.001d0

   qT(i) = 10.d0

   proc(i,1:4)=(/1,1,1,2001/) !KPCs

!    proc(i,1:4)=(/1,1,1,2001/) !LP

   iC(i)=.false.
   cuts(i,1:4)=(/20.d0,20.d0,-2.1d0,2.1d0/)

   MM(i,1:2)=(/1d0,1d0/)

 end do


!  !!!!
!  !!!! Compute point by point
!  !!!!
call cpu_time(time1)
!
!  call xSec_SIDIS_BINLESS_List_forharpy(xx,proc,s,pT,z,x,Q,MM)
!

 do i=1,NUM

  call cpu_time(t1(i))
  xxKPC(i)=KPC_DYconv(Q(i)**2,qT(i),x1(i),x2(i),Q(i),(/1,1,1/))
  call cpu_time(t2(i))
  call cpu_time(t1LP(i))
  xxLP(i)=TMDF_F(Q(i)**2,qT(i),x1(i),x2(i),Q(i),Q(i)**2,Q(i)**2,(/1,1,1/))
  call cpu_time(t2LP(i))

end do

call cpu_time(time2)

do i=1,NUM
write(*,'("{",F12.6,",",F12.6,",",F12.6,"},")') Q(i),Q(i)**2*pi/2*xxKPC(i),xxLP(i)
end do
write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME:", time2-time1
do i=1,NUM
write(*,'("{",F12.6,",",F12.6,",",F12.6,"},")') Q(i),t2(i)-t1(i),(t2(i)-t1(i))/(t2LP(i)-t1LP(i))
end do
write(*,*) "---------------------------------------------------------"


end program example
