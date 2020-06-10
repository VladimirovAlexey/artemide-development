program example
use aTMDe_control
use TMDX_DY
implicit none

real*8,parameter::s=1000d0**2
real*8::p1(1:5),p2(1:5),X(1:5)
integer,dimension(1:3)::proc
real*8::y1,y2,Q1,Q2
integer::i

call artemide_Initialize('const-DYfit18_LO-test',prefix='/misc/data2/braun/vla18041/arTeMiDe_Repository/')

call artemide_SetReplica_uTMDPDF(0)
call artemide_SetReplica_TMDR(0)

p1=(/0.0d0,2.2d0,3.4d0,4.6d0,5.8d0/)
p2=(/2.2d0,3.4d0,4.6d0,5.8d0,7.2d0/)

y1=2.0d0
y2=4.5d0
Q1=20d0
Q2=30d0

proc=(/1,2,1/)
do i=1,5
  call xSec_DY(X(i),proc,s,(/p1(i),p2(i)/),(/Q1,Q2/),(/y1,y2/),.false.,(/20d0,20d0,2.0d0,4.5d0/))
end do
write(*,*) 'process=',proc
do i=1,5
write(*,*) p1(i),'--',p2(i),X(i)!/(y2-y1)/(Q2-Q1)/(p2(i)**2-p1(i)**2)
end do

call artemide_ShowStatistics()
! 
proc=(/1,2,101/)
do i=1,5
  call xSec_DY(X(i),proc,s,(/p1(i),p2(i)/),(/Q1,Q2/),(/y1,y2/),.false.,(/20d0,20d0,2.0d0,4.5d0/))
end do
write(*,*) 'process=',proc
do i=1,5
write(*,*) p1(i),'--',p2(i),X(i)!/(y2-y1)/(Q2-Q1)/(p2(i)**2-p1(i)**2)
end do

call artemide_ShowStatistics()

proc=(/1,2,104/)
do i=1,5
  call xSec_DY(X(i),proc,s,(/p1(i),p2(i)/),(/Q1,Q2/),(/y1,y2/),.false.,(/20d0,20d0,2.0d0,4.5d0/))
end do
write(*,*) 'process=',proc
do i=1,5
write(*,*) p1(i),'--',p2(i),X(i)!/(y2-y1)/(Q2-Q1)/(p2(i)**2-p1(i)**2)
end do

call artemide_ShowStatistics()
write(*,*) '---------------------------------------------------------------'
write(*,*) ' '
write(*,*) '---------------------------------------------------------------'
proc=(/1,1,2/)
do i=1,5
  call xSec_DY(X(i),proc,s,(/p1(i),p2(i)/),(/Q1,Q2/),(/y1,y2/),.false.,(/20d0,20d0,2.0d0,4.5d0/))
end do
write(*,*) 'process=',proc
do i=1,5
write(*,*) p1(i),'--',p2(i),X(i)!/(y2-y1)/(Q2-Q1)/(p2(i)**2-p1(i)**2)
end do

call artemide_ShowStatistics()

proc=(/1,1,102/)
do i=1,5
  call xSec_DY(X(i),proc,s,(/p1(i),p2(i)/),(/Q1,Q2/),(/y1,y2/),.false.,(/20d0,20d0,2.0d0,4.5d0/))
end do
write(*,*) 'process=',proc
do i=1,5
write(*,*) p1(i),'--',p2(i),X(i)!/(y2-y1)/(Q2-Q1)/(p2(i)**2-p1(i)**2)
end do

call artemide_ShowStatistics()

proc=(/1,1,103/)
do i=1,5
  call xSec_DY(X(i),proc,s,(/p1(i),p2(i)/),(/Q1,Q2/),(/y1,y2/),.false.,(/20d0,20d0,2.0d0,4.5d0/))
end do
write(*,*) 'process=',proc
do i=1,5
write(*,*) p1(i),'--',p2(i),X(i)!/(y2-y1)/(Q2-Q1)/(p2(i)**2-p1(i)**2)
end do

call artemide_ShowStatistics()

end program example