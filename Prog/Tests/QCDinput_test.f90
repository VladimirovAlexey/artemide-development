!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
!use aTMDe_control
use QCDinput
implicit none

real*8, dimension(-5:5)::qq

call QCDinput_Initialize('const-1',prefix='Prog/Tests/const-files/')

qq=xPDF(0.1d0,10d0,1)
write(*,*) qq
write(*,*) '---------'
qq=x_hPDF(0.1d0,10d0,5)
write(*,*) qq

end program example
