!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calculation of DY cross-section for CMS.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program xSec_DY
use aTMDe_control
use TMDX_DY
implicit none

real*8::s,Q
integer::i,j
integer,parameter::Npoints=60
real*8,parameter::fract=0.95d0 !!!! amount of cross-seciton
real*8,parameter::bMax=5d0
real*8::qT(1:Npoints),Xexact(1:Npoints),X1(1:Npoints),X2(1:Npoints),X3(1:Npoints),X4(1:Npoints),X5(1:Npoints),X6(1:Npoints)


!!! this program evaluate cross-section for SIDIS



!This is from the artemide ver.1.4
!   call TMDX_DY_Initialize("NNLO")
!   call TMDX_DY_SetNPParameters((/3.3235d0,0.0380d0,0.2204d0, 7.0808d0,351.7950d0, 2.4632d0,-3.8334d0,0.0001d0, 0.0000d0/))

!This is for artemide ver.2.00
  !call artemide_Initialize('constants-file','Models/BSV19.bFIT/')
  call artemide_Initialize('const-DYfit18_NNLO+LP',prefix='/misc/data2/braun/vla18041/arTeMiDe_Repository/')
  !call artemide_SetNPparameters_uTMDPDF((/0.2578d0,8.1849d0,300.0770d0,2.4371d0,-4.7716d0,0.0d0,-0.0d0/))
  call artemide_SetNPparameters_uTMDPDF((/0.228d0,0.306d0/))
  call artemide_SetNPparameters_lpTMDPDF((/0.228d0,0.306d0/))
  call artemide_SetNPparameters_TMDR((/3.3021d0, 0.0237d0/))
  
  Q=125d0
  s=8000d0**2
  
  call TMDX_DY_SetProcess(1,5,20)
  call TMDX_DY_SetCuts(.false.,20d0,20d0,-10d0,10d0)
  call TMDX_DY_XSetup(s,Q,0d0)
  
  call artemide_SetNPparameters_uTMDPDF((/0.2578d0,8.1849d0,300.0770d0,2.4371d0,-4.7716d0,0.0d0,-0.0d0/))
  
  
  
  do i=1,Npoints
    
    qT(i)=0.35d0*Q*i/Npoints
  
  end do
  
    call artemide_SetScaleVariations(1d0,1d0,1d0,1d0)
    call CalcXsec_DY(Xexact,qT)
    
    call artemide_SetScaleVariations(1d0,0.5d0,1d0,1d0)
    call CalcXsec_DY(X1,qT)
    
    call artemide_SetScaleVariations(1d0,2d0,1d0,1d0)
    call CalcXsec_DY(X2,qT)
    
!     call artemide_SetScaleVariations(1d0,1d0,0.5d0,1d0)
!     call CalcXsec_DY(X3,qT)
!     
!     call artemide_SetScaleVariations(1d0,1d0,2d0,1d0)
!     call CalcXsec_DY(X4,qT)
    
    call artemide_SetScaleVariations(1d0,1d0,1d0,0.5d0)
    call CalcXsec_DY(X5,qT)
    
    call artemide_SetScaleVariations(1d0,1d0,1d0,2d0)
    call CalcXsec_DY(X6,qT)

  do i=1,Npoints
  
  write(*,*) '{',qT(i),',',Xexact(i),',',minval((/Xexact(i),X1(i),X2(i),X5(i),X6(i)/)),&
  ',',maxval((/Xexact(i),X1(i),X2(i),X5(i),X6(i)/)),'},'
  end do
 end program xSec_DY