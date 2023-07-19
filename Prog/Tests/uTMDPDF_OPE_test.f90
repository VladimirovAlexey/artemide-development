!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDPDF_OPE
implicit none

real*8,dimension(0:200,0:200)::m1,m2
integer::i,j
real::t1,t2,t3,t4

call cpu_time(t1)
m1=Tmatrix(1,3,1,5d0)
call cpu_time(t2)

! do i=0,5
! write(*,*) m1(i,:)
! end do

write(*,*) "-----------------"

call cpu_time(t3)
m2=Tmatrix2(1,3,1,5d0)
call cpu_time(t4)

! do i=0,5
! write(*,*) m2(i,:)
! end do

write(*,*) "-----------------"
write(*,*) sum(m1-m2)

write(*,*) "time1 =", t2-t1
write(*,*) "time2 =", t4-t3


end program example
