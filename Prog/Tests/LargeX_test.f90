!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_Numerics
use aTMDe_control
use uTMDPDF_OPE
implicit none

integer::i
real*8::b,x
real*8,dimension(-5:5)::T1,T2

call artemide_Initialize('uTMDPDF.atmde',prefix='Prog/Tests/const-files/')

do i=0,20
!do i=10,10
  x=10.d0**(-real(i)/50)

  T1=CxF_compute(x,3.5d0,1,.false.)
  T2=CxF_largeX_compute(x,3.5d0,1,.false.)
  write(*,'("{",F12.6,",",F12.6,",",F12.6,",",F12.6,"},")') x,T1(1),T2(1),(T1(1)-T2(1))/(abs(T1(1))+1.d-8)

  !T1=uTMDPDF_OPE_convolution(x,0.2d0,1)
  !write(*,'("{",F12.6,",",F12.6,",",F12.6,",",F12.6,"},")') x,T1(1),T1(2),T1(3)

end do
end program example
