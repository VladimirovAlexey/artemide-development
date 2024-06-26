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
!!! Here we also use uTMDPDF module
use uTMDPDF

!!! Variable specification
implicit none

integer::i
real*8::b,kT,x
!!!! TMD distributions are vectors in flavor-space (-5:5) = (bBar,cBar,sBar,uBar,dBar,gluon,d,u,s,c,b)
real*8,dimension(-5:5)::tmd1,tmd2,tmd3

!!! Initialization of uTMDPDF-module
!!! It automaticalyl initializes TMDR module and all lower modules
call uTMDPDF_Initialize('uTMDPDF.atmde',prefix='Prog/MiniCourse/INI/')

!!! Setup parameters of CS kernel
call TMDR_SetNPparameters((/1.5d0, 0.0369d0, 0.0582d0, 0.d0/))
!!! Setup parameters for ART23 model
call uTMDPDF_SetLambdaNP((/0.565d0, 0.0539d0, 0.5697d0, 6.64d0, 0.d0, 20.07d0, 0.d0, 0.537d0, 1.07d0, 2.39d0, 0.0d0, 0.0d0/))

!!! value of x
x=0.05d0

write(*,*) '------------ TMDs in b-space ---------------------'
do i=-20,10
  !!! value of b
  b=10.d0**(real(i)/10)

  !!! value of optimal TMDPDF
  tmd1=uTMDPDF_inB(x,b,1)

  !!! value of TMDPDF at 4 GeV
  tmd2=uTMDPDF_inB(x,b,4.d0,(4.d0)**2,1)

  !!! value of TMDPDF at 91 GeV
  tmd3=uTMDPDF_inB(x,b,91.d0,(91.d0)**2,1)

  write(*,'("{",F12.6,",",F12.6,",",F16.6,",",F16.6,",",F16.6,",",F16.6,",",F16.6,"},")') &
        b,tmd1(1),tmd1(2),tmd2(1),tmd2(2),tmd3(1),tmd3(2)

end do


write(*,*) '------------ TMDs in kT-space ---------------------'
do i=-20,10
  !!! value of b
  kT=10.d0**(real(i)/10)

  !!! value of optimal TMDPDF
  tmd1=uTMDPDF_inKT(x,kT,1)

  !!! value of TMDPDF at 4 GeV
  tmd2=uTMDPDF_inKT(x,kT,4.d0,(4.d0)**2,1)

  !!! value of TMDPDF at 91 GeV
  tmd3=uTMDPDF_inKT(x,kT,91.d0,(91.d0)**2,1)

  write(*,'("{",F12.6,",",F12.6,",",F16.6,",",F16.6,",",F16.6,",",F16.6,",",F16.6,"},")') &
        kT,tmd1(1),tmd1(2),tmd2(1),tmd2(2),tmd3(1),tmd3(2)

end do

!!!! To discuss:
!!!!  to change parameters (uTMDPDF_SetLambdaNP)
!!!!  to change PT-order (go to INI-file)
!!!!  to change model (go to model-file and recompile)
!!!!  option of KT-grid, and b-Grid
!!!!  What happens if forget to setup module

end program example
