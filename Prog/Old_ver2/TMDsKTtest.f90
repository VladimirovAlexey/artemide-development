program example
use TMDs_inKT
use TMDR
use uTMDPDF
implicit none

real*8:: A(1:2),B(1:7)
real*8::qT,Q2,mu,x1,x2
 A=(/2.5000d0,0.0376d0/)
 B=(/0.257986d0, 8.2054d0, 301.141d0, 2.43837d0, -4.78378d0,0d0,0d0/)
  
 call TMDs_inKT_Initialize('const-DYfit18_NNLO+LP','/home/vla18041/LinkData2/arTeMiDe_Repository/')
 
 qT=0.3d0
 Q2=25d0
 x1=0.1d0
 x2=0.04d0
 mu=Sqrt(Q2)

 call TMDR_setNPparameters(A)
 call uTMDPDF_SetLambdaNP(B)

 write(*,*) uTMDPDF_kT_5(x1,qT,1)
 
 write(*,*) uTMDPDF_kT_5(x1,qT,mu,Q2,1)

end program example