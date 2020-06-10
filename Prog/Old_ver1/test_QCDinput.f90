program example
use QCDinput
use EWinput
implicit none

  real*8::f1(-5:5),f2(-5:5),f3(-5:5),f4(-5:5)
  
  call EWinput_Initialize('const-test',prefix='/misc/data2/braun/vla18041/arTeMiDe_Repository/artemide/')
  
 call QCDinput_Initialize('const-test',prefix='/misc/data2/braun/vla18041/arTeMiDe_Repository/artemide/')
 
 f1=xPDF(0.1d0,91d0,1)
 f2=xPDF(0.1d0,91d0,2)
 f3=xPDF(0.1d0,91d0,3)
 
 write(*,*) f1(1),f2(1),f3(1)
 
 f1=xPDF(0.1d0,91d0,5)
 f2=xFF(0.1d0,91d0,1)
 f3=xFF(0.1d0,91d0,2)
 
 write(*,*) f1(1),f2(1),f3(1)
 
end program example