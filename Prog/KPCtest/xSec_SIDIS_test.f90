program example

! Invocation of required modules
use aTMDe_control
use uTMDPDF
use uTMDFF
use TMDX_SIDIS


implicit none


integer, parameter :: NUM = 4   ! defines the number of points to compute

integer :: i

real*8, dimension(1:NUM,1:2) :: pT,Q,x,z ! (NUM,2) matrix to store qT, Q and y values
integer, dimension(1:NUM,1:4) :: proc ! (NUM,4) matrix to store the process parameters
logical, dimension(1:NUM) :: iC
real*8, dimension(1:NUM,1:4) :: cuts

real*8, dimension(1:NUM) :: xx,s
real*8, dimension(1:NUM) :: pT_b,Q_b,x_b,z_b
real*8, dimension(1:NUM,1:2) :: MM

real*8 :: time1,time2
!$ real*8::omp_get_wtime



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
  call artemide_SetNPparameters_uTMDFF((/0.253434d0, 9.04351d0, 6.999d0, 2.47992d0/)) ! Lo he sacado del test uTMDFF.f90 de artemide


do i=1,NUM

  Q(i,1:2)=(/10.d0,20.d0/)
  Q_b(i)=(Q(i,1)+Q(i,2))/2

  s(i)=1500.d0**2 !! TeV

  x(i,1:2)=(/0.3d0,0.4d0/)
  x_b(i)=(x(i,1)+x(i,2))/2

  z(i,1:2)=(/0.4d0,0.5d0/)
  z_b(i)=(z(i,1)+z(i,2))/2

  pt(i,1:2)=(/(i-1)*0.1d0,i*0.1d0/)
  pt_b(i)=(pt(i,1)+pt(i,2))/2

  proc(i,1:4)=(/1,1,1,2001/)
  iC(i)=.false.
  cuts(i,1:4)=(/20.d0,20.d0,-2.1d0,2.1d0/)

  MM(i,1:2)=(/1d0,1d0/)
end do

!!!!
!!!! Compute point by point
!!!!
call cpu_time(time1)
!$ time1=omp_get_wtime()

!do i=1,NUM
  !call xSec_SIDIS(xx(i),proc(i,:),s(i),pT(i,:),z(i,:),x(i,:),Q(i,:),iC(i),Cuts(i,:))

!end do

call xSec_SIDIS_BINLESS_List_forharpy(xx,proc,s,pT_b,z_b,x_b,Q_b,MM)

call cpu_time(time2)
!$ time2=omp_get_wtime()

do i=1,NUM
  write(*,'("{",F12.8,",",F16.10,"},")',advance="no") pT_b(i),xx(i)
end do

write(*,*) " "
write(*,*) " "
write(*,*) " COMPUTATION TIME WITH Point-by-Point:", time2-time1

write(*,*) "---------------------------------------------------------"



end program example
