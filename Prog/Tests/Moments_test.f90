!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot for AS-term
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use uTMDPDF_OPE
use uTMDPDF
use TMDs_inKT
implicit none

real*8,dimension(0:200,-5:5)::TT1
integer::i
real*8::x,mu

!call artemide_Initialize('uTMDPDF_AS.atmde',prefix='Prog/Tests/const-files/')
call artemide_Initialize('ART23_MSHT_N4LL.atmde',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_uTMDPDF((/&
0.874245d0,0.913883d0,0.991563d0,6.05412d0,&
0.353908d0,46.6064d0,0.115161d0,1.53235d0,&
1.31966d0,0.434833d0, 0.0d0, 0.0d0/))

TT1(0,:)=Moment_G(0,0.1d0,10.d0,uTMDPDF_lowScale5,1)
write(*,*) TT1(0,:)

TT1(0,:)=Moment_G(1,0.1d0,10.d0,uTMDPDF_lowScale5,1)
write(*,*) TT1(0,:)

end program example
