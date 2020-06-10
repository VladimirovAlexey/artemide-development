!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Testing problem of SIDIS
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use uTMDFF
use TMDX_SIDIS
implicit none

  real*8::fff(-5:5)
  real*8::x,b
  integer::i
  real*8::xx

  call artemide_Initialize('const-DY+SIDIS_NNLO-test','/home/vla18041/LinkData2/arTeMiDe_Repository/')
  
  call artemide_SetReplica_TMDR(0)
  call artemide_SetReplica_uTMDPDF(0)
  call artemide_SetNPparameters_uTMDFF((/0.18086d0,0.0003d0,0.0112d0,0.00015d0,0d0/))
  
  call artemide_ShowStatistics()
  
  
!   do i=1,10
!   x=0.3d0
!   b=0.1d0+i*0.5d0
!   fff=uTMDFF_lowScale5(x,b,3)
!   
!   write(*,*) x,b,fff(1),fff(2)
!   
!   end do
!   write(*,*) '------------------------------------------------'
!   do i=1,10
!   x=0.3d0
!   b=0.1d0+i*0.5d0
!   fff=uTMDFF_lowScale5(x,b,4)
!   
!   write(*,*) x,b,fff(1),fff(2)
!   
!   end do

   call xSec_SIDIS(xx,(/1,1,2001/),52.65d0,(/0.15d0,0.25d0/),(/0.375d0,0.475d0/),(/0.35d0,0.6d0/),(/1d0,4.5d0/),&
	.true.,(/0.1d0,0.85d0,10d0/),(/0.938d0,0.139d0/))
	
  write(*,*) '-----',xx
  
end program example