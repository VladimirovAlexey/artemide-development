!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot for AS-term
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use uTMDPDF
implicit none

real*8,dimension(0:200,-5:5)::TT1
integer::i
real*8::x,mu,b

!call artemide_Initialize('uTMDPDF_AS.atmde',prefix='Prog/Tests/const-files/')
call artemide_Initialize('KT.atmde',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_TMDR((/2.d0,0.03d0,0.003d0,1.d0/))
call artemide_SetNPparameters_uTMDPDF((/&
0.874245d0,0.913883d0,0.991563d0,6.05412d0,&
0.353908d0,46.6064d0,0.115161d0,1.53235d0,&
1.31966d0,0.434833d0, 0.0d0, 0.0d0/))

x=0.01d0
mu=100.d0

do i=0,200
b=0.01+0.1*i
! !TT1(i,:)=uTMDPDF_inB(x,b,mu,mu**2,1)
!TT1(i,:)=uTMDPDF_inKT(x,b,mu,mu**2,1)
TT1(i,:)=ExtractFromGrid_inKT(x,b,mu,1)
write(*,'("{",F10.5,",",F16.12,"},")',advance='no') b,TT1(i,-2)
end do
write(*,*) " "



end program example
