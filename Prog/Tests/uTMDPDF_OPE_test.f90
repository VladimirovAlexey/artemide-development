!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDPDF_OPE
use uTMDPDF
implicit none

real*8,dimension(0:10,1:5,-5:5)::TT1,TT2
integer::i,j
real::t1,t2,t3,t4

call uTMDPDF_Initialize('const-uTMDPDF',prefix='Prog/Tests/const-files/')
call uTMDPDF_OPE_Initialize("file")


call cpu_time(t1)
do i=0,10
do j=1,5
    TT1(i,j,:)=uTMDPDF_lowScale50(10d0**(-i/5d0)-0.00001d0,j/2d0,1)
end do
end do
call cpu_time(t2)
!write(*,*) TT1

call cpu_time(t3)
do i=0,10
do j=1,5
    TT2(i,j,:)=CxF_compute(10d0**(-i/5d0)-0.00001d0,j/2d0,1,.true.)/(10d0**(-i/5d0)-0.00001d0)
end do
end do

call cpu_time(t4)
!write(*,*) TT2
write(*,*) "-------------------------------"
do i=0,10
do j=1,5
    write(*,*) "x=",10d0**(-i/5d0)-0.00001d0,"  b=",j/2d0,"   delta=",sum(TT2(i,j,:)-TT1(i,j,:))/sum(TT2(i,j,:))

end do
end do
write(*,*) "-------------------------------"
write(*,*) "t1=",t2-t1,"    t2=",t4-t3

end program example
