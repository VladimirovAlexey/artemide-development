program example
use uTMDPDF
implicit none

real*8,dimension(-5:5):: F1,F2,F3,F4,F5,F6,F7,F8

real*8,dimension(1:9):: A,B,C
real*8::x,bT
 A=(/2.5000d0,0.0376d0,0.2481d0,8.1556d0,276.5180d0,2.5262d0,-4.9586d0, 0.d0,0d0/)
 B=(/2.5000d0,0.0376d0,0.3481d0,4.1556d0,176.5180d0,2.5262d0,-2.9586d0, 0.d0,0d0/)
 C=(/2.5000d0,0.0376d0,0.1481d0,9.1556d0,276.5180d0,2.5262d0,-6.9586d0, 0.d0,0d0/)
  
  call uTMDPDF_Initialize('const-test')
 
 x=0.001d0
 bT=150d0

call uTMDPDF_SetLambdaNP(A,.false.,.false.)

F1=uTMDPDF_lowScale5(x,bT,1)

call uTMDPDF_SetLambdaNP(B,.false.,.false.)

F2=uTMDPDF_lowScale5(x,bT,1)

call uTMDPDF_SetLambdaNP(A,.True.,.false.)

call uTMDPDF_SetLambdaNP(A,.true.,.false.)

F3=uTMDPDF_lowScale5(x,bT,1)

call uTMDPDF_SetLambdaNP(B,.true.,.false.)

F4=uTMDPDF_lowScale5(x,bT,1)

write(*,*) F1(1:5)
write(*,*) F3(1:5)
Write(*,*) (F1(1:5)-F3(1:5))/F1(1:5)
Write(*,*) (F2(1:5)-F4(1:5))/F4(1:5)

end program example