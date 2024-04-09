program example
use LeptonCutsDY
implicit none

real*8::CP(1:4),Q,qT,y,RR
integer::i,t

call InitializeLeptonCutDY(0.1d-4,1.d-8)

t=8

CP=(/20.d0,20.d0,-2.1d0,2.1d0/)

Q=91.d0
y=0.1

write(*,*) "------------------ vs KT ---------------------------------"

do i=0,60
qT=i*0.5d0
RR=CutFactor(qT,Q,y,CP,t)

write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT,RR
end do
write(*,*) " "

write(*,*) "------------------ vs y ---------------------------------"

qT=5.d0

do i=0,40
y=-2.d0+0.1*i
RR=CutFactor(qT,Q,y,CP,t)

write(*,'("{",F12.8,",",F16.10,"},")',advance="no") y,RR
end do
write(*,*) " "

end program example
