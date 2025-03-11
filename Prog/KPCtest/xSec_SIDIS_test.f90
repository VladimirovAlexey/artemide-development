program example

! Invocation of required modules
use aTMDe_control
use uTMDPDF
use uTMDFF
use TMDX_SIDIS
use TMDF
use TMDF_KPC


implicit none


 integer, parameter :: NUM = 10  ! defines the number of points to compute

 integer :: i ! to do the loop

 real*8, dimension(1:NUM) :: s,pT,z,x,Q
 integer, dimension(1:NUM,1:4) :: proc
 logical, dimension(1:NUM) :: iC
 real*8, dimension(1:NUM,1:4) :: cuts

 real*8, dimension(1:NUM) :: xx,xxKPC,xxLP
 real*8, dimension(1:NUM,1:2) :: MM


 real*8 :: time1,time2


 !!!!! SETUP ARTEMIDE and NP-parameters
 !
 !   Initialization of arTeMiDe using the constants-file xSec_SIDIS.atmde
   call artemide_Initialize('xSec_SIDIS.atmde',prefix='Prog/KPCtest/INI/')
 !   call artemide_Initialize('KPC_my_test.atmde',prefix='Prog/KPCtest/INI/')
 !
 !   Setting the non-perturbative parameters (NP) for TMD evolution
   call artemide_SetNPparameters_TMDR((/1.5004d0, 0.05614d0, 0.03862d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse

 !   Setting the non-perturbative parameters (NP) for the uTMDPDF
   call artemide_SetNPparameters_uTMDPDF(&
   (/0.565d0, 0.0539d0, 0.5697d0, 6.64d0, 0.565d0, 20.07d0, 0.5697d0, 0.537d0, 1.07d0, 2.39d0, 0.0d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse

 !   como SIDIS involucra un hadrón en el estado final, necesito trabajar con FFs también
 !   Setting the non-perturbative parameters (NP) for the uTMDFF
   call artemide_SetNPparameters_uTMDFF((/0.69769d0, 0.712969d0, -0.133895d0, -0.841651d0, 0.846846d0, &
                       0.774759d0, 1.5565d0, 1.1863d0, 0.692877d0, -0.569062d0, 0.0d0, 0.0d0/)) ! Lo he sacado del test uTMDFF.f90 de artemide

!!(Q2, qT_in, x1, z1 ,mu, proc1)
!   !xxKPC(1)=KPC_SIDISconv(1000.d0,0.001d0,0.5d0,0.5d0,10.d0,(/1,1,999/))
!    xxKPC(1)=KPC_SIDISconv(100.d0,2.d0,0.5d0,0.5d0,10.d0,(/1,1,999/))
!     stop
!   write(*,*) ">>>>",xxKPC(1)/0.5d0
!
!   stop

 do i=1,NUM

   Q(i) = 1.d0+(i-1)*2

   s(i) = 1500.d0 !! TeV

   x(i) = 0.2d0
   z(i) = 0.5d0

   pT(i) = 0.1d0

   proc(i,1:4)=(/1,1,1,2001/) !KPCs

!    proc(i,1:4)=(/1,1,1,2001/) !LP

   iC(i)=.false.
   cuts(i,1:4)=(/20.d0,20.d0,-2.1d0,2.1d0/)

   MM(i,1:2)=(/1d0,1d0/)

 end do


!  !!!!
!  !!!! Compute point by point
!  !!!!
!  call cpu_time(time1)
!
!  call xSec_SIDIS_BINLESS_List_forharpy(xx,proc,s,pT,z,x,Q,MM)
!
!  call cpu_time(time2)

 do i=1,NUM

  xxKPC(i)=KPC_SIDISconv(Q(i)**2,pT(i),x(i),z(i),Q(i),(/1,1,2001/))
  xxLP(i)=TMDF_F(Q(i)**2,pT(i),x(i),z(i),Q(i),Q(i)**2,Q(i)**2,(/1,1,2001/))


  write(*,'("{",F12.6,",",F12.6,",",F12.6,",",F12.6,"},")') Q(i),Q(i)**2*3.14d0/2/z(i)*xxKPC(i),xxLP(i),&
    Q(i)**2*3.14d0/2/z(i)*xxKPC(i)/xxLP(i)

end do

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! Write data to file (KPCs case)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   open(unit=13, file="xSec_KPCs.txt", status="replace", action="write")
!   do i = 1, NUM
!     write(13, *) Q(i),xx(i)
!   end do
!   close(13)
!
!  write(*,*) " "
!  write(*,*) " "
!  write(*,*) " COMPUTATION TIME WITH Point-by-Point:", time2-time1
!
!  write(*,*) "---------------------------------------------------------"
!
!  write(*,*) "Data written to ", "xSec_KPCs.txt"
!  write(*,*) "Use an external plotting program to visualize the data."


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! Write data to file (LP case)
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   open(unit=13, file="xSec_LP.txt", status="replace", action="write")
!   do i = 1, NUM
!     write(13, *) Q(i),xx(i)
!   end do
!   close(13)
!
!  write(*,*) " "
!  write(*,*) " "
!  write(*,*) " COMPUTATION TIME WITH Point-by-Point:", time2-time1
!
!  write(*,*) "---------------------------------------------------------"
!
!  write(*,*) "Data written to ", "xSec_LP.txt"
!  write(*,*) "Use an external plotting program to visualize the data."





end program example
