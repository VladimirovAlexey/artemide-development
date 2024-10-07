!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Example code for mini-course of artemide (Madrid, june 2024)          !
!             USE ART23 MODEL!                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!
!!! This program makes the plot of CS kernel
!!!

program example
use uTMDPDF
use TMDR

!!! Variable specification
implicit none

integer::i,j
real*8::x,b,R1,R2,zz,pdf0(-5:5),pdf1(-5:5),pdf2(-5:5)

!!! Initialization of TMDR-module
!!! It automaticalyl initializes QCDinput
call uTMDPDF_Initialize('TMDPDF.atmde',prefix='Prog/LargeX/INI/')
call TMDR_Initialize('TMDPDF.atmde',prefix='Prog/LargeX/INI/')

!!! Setup parameters of CS kernel (see
call TMDR_SetNPparameters((/1.5d0, 0.0369d0, 0.0582d0, 0.d0/))

call uTMDPDF_SetLambdaNP((/0d0, 0d0, 0d0, 0d0, 0.0d0, 0d0, 0.d0, 0d0, 0d0, 0d0, 0d0, 0.04d0/))

OPEN(UNIT=51, FILE="/data/WorkingFiles/TMD/Fit_Notes/LargeX/Data/uTMDPDF_N2LO_N2LL", ACTION="write", STATUS="replace")


do i=-10,7
b=10**(real(i)/10)
do j=1,8
  !!! value of b
  x=0.1*j
  !!! value of CS-kernel at 4 GeV
  pdf0=uTMDPDF_inB(x,b,1)
  write(51,'(F12.6,", ",F12.6,", ",E16.10,", ",E16.10,", ",E16.10)') x,b,pdf0(1),pdf0(2),pdf0(3)
end do
do j=1,18
  !!! value of b
  x=0.8+0.01*j
  !!! value of CS-kernel at 4 GeV
  pdf0=uTMDPDF_inB(x,b,1)
  write(51,'(F12.6,", ",F12.6,", ",E16.10,", ",E16.10,", ",E16.10)') x,b,pdf0(1),pdf0(2),pdf0(3)
end do
do j=1,19
  !!! value of b
  x=0.98+0.001*j
  !!! value of CS-kernel at 4 GeV
  pdf0=uTMDPDF_inB(x,b,1)
  write(51,'(F12.6,", ",F12.6,", ",E16.10,", ",E16.10,", ",E16.10)') x,b,pdf0(1),pdf0(2),pdf0(3)
end do
end do

CLOSE (51, STATUS='KEEP')


end program example
