!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program tests the bin integration for SIDIS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDF
implicit none

real*8::Q1,Q2,x1,x2,b


integer::i,iMax


call artemide_Initialize('const-DY_test',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

Q1=10d0
Q2=15d0
x1=0.2d0
x2=0.03d0

do i=0,20
b=i*0.2+0.01
write(*,'("[",F6.4,",",F8.6,"],")')b, Integrand(Q1**2,b,x1,x2,Q1,Q1**2,Q1**2,5)/Integrand(Q2**2,b,x1,x2,Q2,Q2**2,Q2**2,5)
end do

end program example
