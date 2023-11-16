!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot for AS-term
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDPDF_OPE
implicit none

real*8,dimension(0:200,-5:5)::TT1,TT2
integer::i
real::t1,t2,t3,t4
real*8::x,mu

call uTMDPDF_OPE_Initialize('uTMDPDF_AS.atmde',prefix='Prog/Tests/const-files/')

! x=0.001d0
!
! call cpu_time(t1)
! do i=0,200
!     TT1(i,:)=uTMDPDF_X0_AS(x,(20d0+1d0)**(i/85d0),1)*x
! end do
! call cpu_time(t2)
!
! !write(*,*) TT2
! write(*,*) "------- x= 0.001------------------------"
! do i=0,200
!     write(*,'("{",F12.6,",",F16.10,",",F16.10,",",F16.10,",",F16.10,"},")') &
!         (20d0+1d0)**(i/85d0),TT1(i,-2),TT1(i,-1),TT1(i,1),TT1(i,2)
! end do
! write(*,*) "-------------------------------"
! write(*,*) "t1=",t2-t1,"    t2=",t4-t3
!
!
! x=0.01d0
!
! call cpu_time(t1)
! do i=0,200
!     TT1(i,:)=uTMDPDF_X0_AS(x,(20d0+1d0)**(i/85d0),1)*x
! end do
! call cpu_time(t2)
!
! !write(*,*) TT2
! write(*,*) "------- x= 0.01------------------------"
! do i=0,200
!     write(*,'("{",F12.6,",",F16.10,",",F16.10,",",F16.10,",",F16.10,"},")') &
!         (20d0+1d0)**(i/85d0),TT1(i,-2),TT1(i,-1),TT1(i,1),TT1(i,2)
! end do
! write(*,*) "-------------------------------"
! write(*,*) "t1=",t2-t1,"    t2=",t4-t3
!
!
! x=0.1d0
!
! call cpu_time(t1)
! do i=0,200
!     TT1(i,:)=uTMDPDF_X0_AS(x,(20d0+1d0)**(i/85d0),1)*x
! end do
! call cpu_time(t2)
!
! !write(*,*) TT2
! write(*,*) "------- x= 0.1------------------------"
! do i=0,200
!     write(*,'("{",F12.6,",",F16.10,",",F16.10,",",F16.10,",",F16.10,"},")') &
!         (20d0+1d0)**(i/85d0),TT1(i,-2),TT1(i,-1),TT1(i,1),TT1(i,2)
! end do
! write(*,*) "-------------------------------"
! write(*,*) "t1=",t2-t1,"    t2=",t4-t3

mu=10.d0

call cpu_time(t1)
do i=0,100
    x=(0.0001d0)**(i/100d0)
    TT1(i,:)=uTMDPDF_X0_AS(x,mu,1)*x
end do
call cpu_time(t2)

!write(*,*) TT2
write(*,*) "------- mu= 10------------------------"
do i=0,100
     x=(0.0001d0)**(i/100d0)
    write(*,'("{",F12.6,",",F16.10,",",F16.10,",",F16.10,",",F16.10,"},")') &
        x,TT1(i,-2),TT1(i,-1),TT1(i,1),TT1(i,2)
end do
write(*,*) "-------------------------------"
write(*,*) "t1=",t2-t1,"    t2=",t4-t3

mu=100.d0

call cpu_time(t1)
do i=0,100
    x=(0.0001d0)**(i/100d0)
    TT1(i,:)=uTMDPDF_X0_AS(x,mu,1)*x
end do
call cpu_time(t2)

!write(*,*) TT2
write(*,*) "------- mu= 100------------------------"
do i=0,100
     x=(0.0001d0)**(i/100d0)
    write(*,'("{",F12.6,",",F16.10,",",F16.10,",",F16.10,",",F16.10,"},")') &
        x,TT1(i,-2),TT1(i,-1),TT1(i,1),TT1(i,2)
end do
write(*,*) "-------------------------------"
write(*,*) "t1=",t2-t1,"    t2=",t4-t3

mu=1000.d0

call cpu_time(t1)
do i=0,100
    x=(0.0001d0)**(i/100d0)
    TT1(i,:)=uTMDPDF_X0_AS(x,mu,1)*x
end do
call cpu_time(t2)

!write(*,*) TT2
write(*,*) "------- mu= 1000------------------------"
do i=0,100
     x=(0.0001d0)**(i/100d0)
    write(*,'("{",F12.6,",",F16.10,",",F16.10,",",F16.10,",",F16.10,"},")') &
        x,TT1(i,-2),TT1(i,-1),TT1(i,1),TT1(i,2)
end do
write(*,*) "-------------------------------"
write(*,*) "t1=",t2-t1,"    t2=",t4-t3

end program example
