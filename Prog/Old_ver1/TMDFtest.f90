program example
use TMDF
use TMDR
use uTMDPDF
implicit none

real*8:: A(1:2),B(1:7)
real*8::qT,Q2,mu,x1,x2
 A=(/2.5000d0,0.0376d0/)
 B=(/0.257986d0, 8.2054d0, 301.141d0, 2.43837d0, -4.78378d0,0d0,0d0/)
  
 call TMDF_Initialize('const-test')
 
 qT=0.3d0
 Q2=25d0
 x1=0.1d0
 x2=0.04d0
 mu=Sqrt(Q2)

 call TMDR_setNPparameters(A)
 call uTMDPDF_SetLambdaNP(B)

 write(*,*) TMDF_F(Q2,qT,x1,x2,mu,Q2,Q2,2)
 
 mu=Sqrt(Q2)*1.2d0
 write(*,*) TMDF_F(Q2,qT,x1,x2,mu,Q2,Q2,2)

end program example