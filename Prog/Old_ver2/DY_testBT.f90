!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calculation of DY cross-section for CMS.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program xSec_DY
use aTMDe_control
use TMDX_DY
implicit none

real*8::s,Q
integer::i,j
integer,parameter::Npoints=40
real*8,parameter::fract=0.95d0 !!!! amount of cross-seciton
real*8,parameter::bMax=5d0
!real*8::qT(1:Npoints),Xexact(1:Npoints)
real*8::qT,Xexact,Xcurrent,bCUT


!!! this program evaluate cross-section for SIDIS



!This is from the artemide ver.1.4
!   call TMDX_DY_Initialize("NNLO")
!   call TMDX_DY_SetNPParameters((/3.3235d0,0.0380d0,0.2204d0, 7.0808d0,351.7950d0, 2.4632d0,-3.8334d0,0.0001d0, 0.0000d0/))

!This is for artemide ver.2.00
  !call artemide_Initialize('constants-file','Models/BSV19.bFIT/')
  call artemide_Initialize('const-DYfit18_NNLO',prefix='/misc/data2/braun/vla18041/arTeMiDe_Repository/')
  call artemide_SetNPparameters_uTMDPDF((/0.2578d0,8.1849d0,300.0770d0,2.4371d0,-4.7716d0,0.0d0,-0.0d0/))
  call artemide_SetNPparameters_TMDR((/3.3021d0, 0.0237d0/))
  
  open(7, file='NP-sense.dat', status="new", action="write")
  
  do j=0,14
  Q=10d0+j*10d0
  
  write(*,*) 'Q=',Q
  s=(100d0*Q)**2
  
  call TMDX_DY_SetProcess(1,1,5)
  call TMDX_DY_SetCuts(.false.,20d0,20d0,-10d0,10d0)
  call TMDX_DY_XSetup(s,Q,0d0)
  
  call artemide_SetNPparameters_uTMDPDF((/0.2578d0,8.1849d0,300.0770d0,2.4371d0,-4.7716d0,0.0d0,-0.0d0/))
    
  do i=1,Npoints
    
    qT=0.25d0*Q*i/Npoints
    
    call artemide_SetNPparameters_uTMDPDF((/0.2578d0,8.1849d0,300.0770d0,2.4371d0,-4.7716d0,0.0d0,-0.0d0/))
    
    call CalcXsec_DY(Xexact,qT)
    !write(*,*) Xexact
    
    call artemide_SetNPparameters_uTMDPDF((/0.0d0,0.0d0,0.0d0,0d0,0d0,0.0d0,100d0/))
    call CalcXsec_DY(Xcurrent,qT)
    
    write(7,*) Q, qT, Xcurrent/Xexact
    
  end do
  end do
  
  CLOSE (7, STATUS='KEEP')

 end program xSec_DY