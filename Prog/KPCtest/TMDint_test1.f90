program example
use aTMDe_control
use TMDint_KPC_DY
implicit none

integer::i
real*8::Q,qT,x1,x2,R

  call artemide_Initialize('KPC1.atmde',prefix='Prog/KPCtest/INI/')

  call artemide_SetNPparameters_TMDR((/1.56d0, 0.0639d0, 0.0582d0,0d0/))

  call artemide_SetNPparameters_uTMDPDF(&
  (/0.874245d0, 0.913883d0, 0.991563d0, 6.05412d0, 0.353908d0,&
  46.6064d0, 0.115161d0, 1.53235d0, 1.31966d0, 0.434833d0, 0.d0, 0.d0/))

Q=21.d0
x1=0.2d0
x2=0.3d0


do i=0,9
  qt=0.1d0+i*0.5
  R=KPC_DYconv(Q,qT,x1,x2,Q,(/1,1,1/),2)
  write(*,'("{",F12.8,",",F16.8,"},")') qt,R

end do


end program example
