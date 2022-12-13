!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use LeptonCutsDY
implicit none

real*8::cc,qT
integer::i

do i=0,200
qT=i*0.1d0
cc=CutFactor4(qT,90d0,-0.1d0,(/0.d0,0.d0,-20.4d0,20.4d0/))

write(*,"('{',F16.10 ,',',F16.10,'},')") qT,cc
!stop
end do



end program example
