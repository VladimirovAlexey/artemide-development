program example
use aTMDe_control
use uTMDPDF
use TMDs_inKT
implicit none

integer,parameter::NUM=40

integer::i

real*8,dimension(1:NUM)::X,qT
real*8,dimension(-5:5)::kkk1,kkk2


  call artemide_Initialize('KPC2.atmde',prefix='Prog/KPCtest/INI/')
  !call artemide_Initialize('LP2.atmde',prefix='Prog/KPCtest/INI/')

  call artemide_SetNPparameters_TMDR((/1.56142d0, 0.0369174d0, 0.0581734d0, 1.0d0/))

  call artemide_SetNPparameters_uTMDPDF(&
  (/0.874245d0, 0.913883d0, 0.991563d0, 6.05412d0, 0.353908d0,&
  46.6064d0, 0.115161d0, 1.53235d0, 1.31966d0, 0.434833d0, 0.d0, 0.d0/))

!!!! set up variables
do i=1,NUM
  qT(i)=i*0.2
  kkk1=uTMDPDF_inKT(0.1d0,qT(i),1)
  kkk2=uTMDPDF_kT_5(0.1d0,qT(i),1)

  write(*,*) kkk1(2)-kkk2(2)
end do


end program example
