!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDPDF_OPE
use uTMDPDF
implicit none

real*8,dimension(0:200,-5:5)::TT1,TT2
integer::i
real::t1,t2,t3,t4

call uTMDPDF_Initialize('const-uTMDPDF',prefix='Prog/Tests/const-files/')
call uTMDPDF_OPE_Initialize("file")

call cpu_time(t1)
do i=0,200
    TT1(i,:)=uTMDPDF_lowScale50(0.001d0,(20d0+1d0)**(i/200d0)-1d0,1)*0.001d0
end do
call cpu_time(t2)
!write(*,*) TT1

call cpu_time(t3)
do i=0,200
    TT2(i,:)=CxF_compute(0.001d0,(20d0+1d0)**(i/200d0)-1d0,1,.true.)
end do

call cpu_time(t4)
!write(*,*) TT2
write(*,*) "-------------------------------"
do i=0,200
    write(*,'("{",F10.6,",",F16.10,",",F16.10,"},")') (20d0+1d0)**(i/200d0)-1d0,TT1(i,1),TT2(i,1)
end do
write(*,*) "-------------------------------"
write(*,*) "t1=",t2-t1,"    t2=",t4-t3

end program example
