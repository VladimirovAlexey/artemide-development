program example
use aTMDe_control
use QCDinput
implicit none

 real*8 :: time1, time2,Q,x,rrr,ss
 real*8::R1(-5:5),R2(-5:5)

 integer::i
 !$  real*8::OMP_get_wtime,t1,t2


 
call artemide_Initialize('test.atmde')
!
!  do i=0,100
!  Q=1.7d0*i+1.75d0
!  x=0.0284!10**(-i/1000.d0-1.5d0)
!  R1=xPDF_R(x,Q,1)
!  R2=xPDF(x,Q,1)
!  !write(*,'("{",F12.10,",",F12.8,",",F12.8,"},")') x,R1(2),R2(2)
!  write(*,'("{",F12.4,",",F12.8,",",F12.8,"},")') Q,R1(2),R2(2)
!  end do

ss=0.d0
 do i=0,1000
 CALL RANDOM_NUMBER(rrr)
 Q=1.7d0+rrr*200.d0
 CALL RANDOM_NUMBER(rrr)
 x=10**(-rrr*6)
 R1=xPDF_R(x,Q,1)
 R2=xPDF(x,Q,1)
 ss=ss+sum((R1(-3:3)-R2(-3:3)))
 if(sum((R1(-3:3)-R2(-3:3))/(R1(-3:3)+R2(-3:3)))>0.1) then
 write(*,*) x,Q, R1(-3:3),R2(-3:3)
 end if
 end do
 write(*,*) "--->",ss/7/10001

 call cpu_time(time1)
 !$  t1=OMP_get_wtime()
 do i=0,10000
 CALL RANDOM_NUMBER(rrr)
 Q=i/50.d0+1.75d0
 x=10**(-rrr*6)
 R1=xPDF(x,Q,1)
 end do
 call cpu_time(time2)
 !$  t2=OMP_get_wtime()

 write(*,*) "LHA-time :", time2-time1,t2-t1

 call cpu_time(time1)
 !$  t1=OMP_get_wtime()
 do i=0,10000
 CALL RANDOM_NUMBER(rrr)
 Q=i/50.d0+1.75d0
 x=10**(-rrr*6)
 R1=xPDF_R(x,Q,1)
 end do
 call cpu_time(time2)
 !$  t2=OMP_get_wtime()

 write(*,*) "MY-time :", time2-time1,t2-t1

end program example
