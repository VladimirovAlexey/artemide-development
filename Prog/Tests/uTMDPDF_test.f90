!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use uTMDPDF
use uTMDPDF_model
implicit none

integer::i,f
real*8::x,b
real*8::TT(-5:5)
real*8::lambda(1:12)

call artemide_Initialize('uTMDPDF.atmde',prefix='Prog/Tests/const-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_uTMDPDF((/&
0.874245d0,0.913883d0,0.991563d0,6.05412d0,&
0.353908d0,46.6064d0,0.115161d0,1.53235d0,&
1.31966d0,0.434833d0, 0.0d0, 0.0d0/))


f=1
x=0.3d0

lambda= uTMDPDF_CurrentLambdaNP()

do i=-40,10
  b=10.d0**(i/real(10,8))
  !TT=FNP(x,b,1,lambda)
  TT=uTMDPDF_lowScale5(x,b,1)
    write(*,'("{",F8.4,",",F16.8,"},")') b,TT(f)
end do
!

end program example
