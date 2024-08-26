program example
use aTMDe_control
use QCDinput
implicit none

 real*8 :: time1, time2,Q,R

 integer::i
 !$  real*8::OMP_get_wtime,t1,t2


 
  call artemide_Initialize('test.atmde')
 
 do i=0,100
 Q=i+1.75d0
 write(*,*) Q,(as(Q)-as1(Q))/as(Q)
 end do

 call cpu_time(time1)
 !$  t1=OMP_get_wtime()
 do i=0,10000
 Q=i/50.d0+1.75d0
 R=as(Q)
 end do
 call cpu_time(time2)
 !$  t2=OMP_get_wtime()

 write(*,*) "LHA-time :", time2-time1,t2-t1

 call cpu_time(time1)
 !$  t1=OMP_get_wtime()
 do i=0,10000
 Q=i/50.d0+1.75d0
 R=as1(Q)
 end do
 call cpu_time(time2)
 !$  t2=OMP_get_wtime()

 write(*,*) "MY-time :", time2-time1,t2-t1

end program example
