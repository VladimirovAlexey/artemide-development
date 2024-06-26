!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Example code for mini-course of artemide (Madrid, june 2024)          !
!             USE ART23 MODEL!                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!
!!! This program makes the plot of CS kernel
!!!

program example
!!! Here we use TMDR module
use TMDR

!!! Variable specification
implicit none

integer::i
real*8::b,R1,R2,zz,CS

!!! Initialization of TMDR-module
!!! It automaticalyl initializes QCDinput
call TMDR_Initialize('CS.atmde',prefix='Prog/MiniCourse/INI/')

!!! Setup parameters of CS kernel (see
call TMDR_SetNPparameters((/1.5d0, 0.0369d0, 0.0582d0, 0.d0/))

do i=-20,10
  !!! value of b
  b=10.d0**(real(i)/10)

  !!! value of evolution factor from optimal line to 4GeV
  R1=TMDR_Rzeta(b,4.d0,16.d0,1)

  !!! value of evolution factor from optimal line to 91GeV
  R2=TMDR_Rzeta(b,91.d0,(91.d0)**2,1)

  !!! value of CS-kernel at 4 GeV
  CS=CS_kernel(4.d0,b,1)

  write(*,'("{",F12.6,",",F12.6,",",F16.6,",",F16.6,"},")') b,R1,R2,CS

end do

!!!! To discuss:
!!!!  to change parameters (TMDR_SetNPparameters)
!!!!  to change PT-order (go to INI-file)
!!!!  to change model (go to model-file and recompile)

end program example
