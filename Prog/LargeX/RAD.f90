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
use TMD_AD, only : Dpert
use uTMDPDF_model

!!! Variable specification
implicit none

integer::i
real*8::b,R1,R2,zz,CS,CSdown,CSup

!!! Initialization of TMDR-module
!!! It automaticalyl initializes QCDinput
call TMDR_Initialize('CS.atmde',prefix='Prog/LargeX/INI/')

!!! Setup parameters of CS kernel (see
call TMDR_SetNPparameters((/1.5d0, 0.0369d0, 0.0582d0, 0.d0/))

write(*,*) "-------------- b*=b -------------------------"
call ModelUpdate((/0.696635d0, 3.7706d0, 0.477978d0, 0.740207d0, 0.0d0, 35.6754d0, &
  0.0d0, 22.7604d0, 1.4951d0, 0.020935d0, 0.0d0, 0.0d0/))

do i=-20,10
  !!! value of b
  b=10.d0**(real(i)/10)

  !!! value of CS-kernel at 4 GeV
  !CS=Dpert(4.d0,b,1)
  CS=2*Dpert(muOPE(b,1.d0,1.d0,1d0),bSTAR(b,1.d0,1.d0),1)
  CSup=2*Dpert(muOPE(b,1.d0,1.d0,2.d0),bSTAR(b,1.d0,1.d0),1)
  CSdown=2*Dpert(muOPE(b,1.d0,1.d0,0.5d0),bSTAR(b,1.d0,1.d0),1)

  write(*,'("{",F12.6,",",F12.6,",",F12.6,",",F12.6,"},")') b,CS,CSdown,CSup

end do


write(*,*) "-------------- b*=c0/mu -------------------------"
call ModelUpdate((/0.696635d0, 3.7706d0, 0.477978d0, 0.740207d0, 0.0d0, 35.6754d0, &
  0.0d0, 22.7604d0, 1.4951d0, 0.020935d0, 0.0d0, 10000.0d0/))

do i=-20,10
  !!! value of b
  b=10.d0**(real(i)/10)

  !!! value of CS-kernel at 4 GeV
  !CS=Dpert(4.d0,b,1)
  CS=2*Dpert(muOPE(b,1.d0,1.d0,1d0),bSTAR(b,1.d0,1.d0),1)
  CSup=2*Dpert(muOPE(b,1.d0,1.d0,2.d0),bSTAR(b,1.d0,1.d0),1)
  CSdown=2*Dpert(muOPE(b,1.d0,1.d0,0.5d0),bSTAR(b,1.d0,1.d0),1)

  write(*,'("{",F12.6,",",F12.6,",",F12.6,",",F12.6,"},")') b,CS,CSdown,CSup

end do

write(*,*) "-------------- b*=b(0.04) -------------------------"
call ModelUpdate((/0.696635d0, 3.7706d0, 0.477978d0, 0.740207d0, 0.0d0, 35.6754d0, &
  0.0d0, 22.7604d0, 1.4951d0, 0.020935d0, 0.0d0, 0.04d0/))

do i=-20,10
  !!! value of b
  b=10.d0**(real(i)/10)

  !!! value of CS-kernel at 4 GeV
  !CS=Dpert(4.d0,b,1)
  CS=2*Dpert(muOPE(b,1.d0,1.d0,1d0),bSTAR(b,1.d0,1.d0),1)
  CSup=2*Dpert(muOPE(b,1.d0,1.d0,2.d0),bSTAR(b,1.d0,1.d0),1)
  CSdown=2*Dpert(muOPE(b,1.d0,1.d0,0.5d0),bSTAR(b,1.d0,1.d0),1)

  write(*,'("{",F12.6,",",F12.6,",",F12.6,",",F12.6,"},")') b,CS,CSdown,CSup

end do

end program example
