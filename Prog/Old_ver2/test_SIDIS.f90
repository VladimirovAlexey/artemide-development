program example
use TMDX_SIDIS
use TMDR
use uTMDPDF
use uTMDFF
implicit none
  
  integer,parameter::ptN=3
  
  integer::i
  real*8::x1,x2,z1,z2,Q1,Q2,s,xL(1:ptN,1:2),zL(1:ptN,1:2),QL(1:ptN,1:2),sL(1:ptN)
  integer::proc(1:3),procL(1:ptN,1:3)
  real*8::cuts(1:4),cutsL(1:ptN,1:4),massL(1:ptN,1:2)
  
  real*8::pT(1:ptN+1),pTL(1:ptN,1:2)
  real*8::xSec(1:ptN)
  logical::doCutL(1:ptN)
  
  real*8,dimension(1:13)::var,var1,var2,var3
!   
  real*8:: A(1:2),B(1:7),B2(1:4)
  A=(/2.5000d0,0.0376d0/)
  B=(/0.257986d0, 8.2054d0, 301.141d0, 2.43837d0, -4.78378d0,0d0,0d0/)
  B2=(/0.264d0,0.479d0,0.459d0,0.539d0/)

   call TMDX_SIDIS_Initialize('const-DY+SIDIS_LO-DSS',&
      '/misc/data2/braun/vla18041/arTeMiDe_Repository/Constants-files/DY+SIDIS-NNLO/')
   
   
   call TMDR_setNPparameters(A)
   call uTMDPDF_SetLambdaNP(B)
   call uTMDFF_SetLambdaNP(B2)

! 	!!kinematicArray(pT,s,z,x,Q,M2target_in,M2product_in)
!   var=kinematicArray(0.2d0,52d0,0.2d0,0.1d0,2.5d0,0.8d0,0.05d0)
!   
!   var1=kinematicArray(var(1),var(11)+var(12),var(5),var(4),var(2),var(12),var(13))
!   var2=kinematicArray(var1(1),var1(11)+var1(12),var1(5),var1(4),var1(2),var1(12),var1(13))
!   
!   write(*,*) var-var1
!   write(*,*) var-var2
!   
!   stop
  !
  proc=(/1,1,2001/)
  x1=0.2d0
  x2=0.35d0
  z1=0.2d0
  z2=0.25d0
  Q1=1d0
  Q2=3d0
  cuts=(/0.1d0,0.85d0,10d0,10000d0/)
  s=52.65d0
  
  pT=(/0.01d0,0.1d0,0.15d0,0.2d0/)
  
   
   call TMDX_SIDIS_SetProcess(proc)
   
   call TMDX_SIDIS_XSetup(s,z1,x1,Q2,0d0,0d0)   
   call TMDX_SIDIS_SetCuts(.true.,cuts(1),cuts(2),cuts(3),cuts(4))   
   !call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(xSec,pT,z1,z2,x1,x2,Q1,Q2)
   call CalcXsec_SIDIS(xSec,pT(1:ptN))
   write(*,*) xSec
  
    proc=(/1,1,2002/)
  x1=0.2d0
  x2=0.35d0
  z1=0.2d0
  z2=0.25d0
  Q1=1d0
  Q2=3d0
  cuts=(/0.1d0,0.85d0,10d0,1000d0/)
  s=52.65d0
  
  pT=(/0.01d0,0.1d0,0.15d0,0.2d0/)
  
   
   call TMDX_SIDIS_SetProcess(proc)
   
   call TMDX_SIDIS_XSetup(s,z1,x1,Q2,0d0,0d0)   
   call TMDX_SIDIS_SetCuts(.true.,cuts(1),cuts(2),cuts(3),cuts(4))   
   !call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(xSec,pT,z1,z2,x1,x2,Q1,Q2)   
   call CalcXsec_SIDIS(xSec,pT(1:ptN))
   write(*,*) xSec

  
!    call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(xSec(1),pT(1),pT(2),z1,z2,x1,x2,Q1,Q2)
!    write(*,*) xSec(1)
   
   !!
   
   do i=1,ptN
    xL(i,:)=(/x1,x2/)
    zL(i,:)=(/z1,z2/)
    QL(i,:)=(/Q1,Q2/)
    sL(i)=s
    procL(i,:)=proc
    cutsL(i,:)=cuts
    massL(i,:)=(/0.92d0,0.12d0/)
    pTL(i,:)=(/pT(i),pT(i+1)/)
    doCutL(i)=.True.
   end do
   
   call xSec_SIDIS_List(xSec,procL,sL,pTL,zL,xL,QL,doCutL,CutsL,massL)
   write(*,*) xSec
   
end program example