program example
use uTMDFF
implicit none

real*8,dimension(-5:5):: F1,F2,F3,F4,F5,F6,F7,F8

real*8,dimension(1:4):: A,B,C
real*8::x,bT
 A=(/0.5000d0,0.0376d0,0.2481d0,8.1556d0/)
 B=(/0.5000d0,0.0376d0,0.3481d0,4.1556d0/)
 C=(/0.5000d0,0.0376d0,0.1481d0,9.1556d0/)
  
!   A=(/2.5000d0,0.0376d0,0.2481d0/)
!  B=(/2.5000d0,0.0376d0,0.3481d0/)
!  C=(/2.5000d0,0.0376d0,0.1481d0/)
  
  call uTMDFF_Initialize('const-test')
 
 x=0.1d0
 bT=0.3d0

call uTMDFF_SetLambdaNP(A,.false.,.false.)

F1=uTMDFF_lowScale5(x,bT,1)

F2=uTMDFF_lowScale5(x,bT,2)


call uTMDFF_SetLambdaNP(A,.true.,.false.)

F3=uTMDFF_lowScale5(x,bT,1)

F4=uTMDFF_lowScale5(x,bT,2)

write(*,*) F1(1:5)
write(*,*) F3(1:5)
Write(*,*) (F1(1:5)-F3(1:5))/F1(1:5)
write(*,*)
write(*,*) F2(1:5)
write(*,*) F4(1:5)
Write(*,*) (F2(1:5)-F4(1:5))/F2(1:5)

end program example