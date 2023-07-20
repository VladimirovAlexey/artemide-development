!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDPDF_OPE
implicit none

real*8,dimension(0:200,0:200)::TT
real*8,dimension(0:200)::FF,RR
integer::i,j
real::t1,t2,t3,t4

call uTMDPDF_OPE_Initialize("file")


call cpu_time(t1)
TT=Tmatrix(-1,1,1,5d0)
call cpu_time(t2)

! do i=0,10
! write(*,*) real(m1(i,6:10))
! end do

write(*,*) "time1 =", t2-t1

do i=0,300
    FF(i)=TESTF(XatNode(i))
end do

RR=MATMUL(TT,FF)

do i=0,200
    write(*,'("{",F8.6,",",F12.10,"},")')XatNode(i),RR(i)
end do


contains
function TESTF(x)
real*8::x,TESTF
    TESTF=1-Log(x)-x**2
end function TESTF

end program example
