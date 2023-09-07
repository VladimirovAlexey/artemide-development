!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use IntegrationRoutines
implicit none

write(*,*) Integrate_SA_2D(F,0.2d0,5.5d0,-2.d0,1.d0,0.0001d0)

contains

function F(x,y)
  real*8::F,x,y

  F=(1-x)*(1+y*Sin(5*y*x+x**2))
end function F

end program example
