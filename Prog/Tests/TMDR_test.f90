!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDR
implicit none

integer::i
real*8::b,R,zz,CS

call artemide_Initialize('TMDR.atmde',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_TMDR((/1.5d0, 0.049d0, 0.0597d0, 0.d0/))

do i=-10,40
  b=10.d0**(real(i)/20)
  R=TMDR_Rzeta(b,4.d0,16.d0,1)
  CS=CS_kernel(4.d0,b,1)
  zz=zetaNP(4.d0,b,1)
  write(*,'("{",F12.6,",",F12.6,",",F16.6,",",F16.6,"},")') b,R,zz,CS

end do


end program example
