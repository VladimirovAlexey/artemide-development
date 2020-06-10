program example
use TMDs
use TMDR
implicit none

real*8,dimension(1:2):: A,B
real*8::bT,Q1,Q2
 A=(/2.5000d0,0.0376d0/)
 B=(/2.5000d0,0.0576d0/)
  
  call TMDs_Initialize('const-test')
 
 bT=0.3d0
 Q1=5d0
 Q2=25d0

  call TMDR_setNPparameters(A)

  write(*,*) TMDR_Rzeta(bT,Q1,Q1**2,1),TMDR_Rzeta(bT,Q2,Q2**2,1)
  
  call TMDR_setNPparameters(B)

  write(*,*) TMDR_Rzeta(bT,Q1,Q1**2,1),TMDR_Rzeta(bT,Q2,Q2**2,1)

end program example