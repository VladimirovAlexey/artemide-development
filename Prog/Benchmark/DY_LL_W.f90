!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_DY
use TMDF
implicit none

integer,dimension(1:3)::proc
real*8::y,Q,s
real*8::stepqT,maxqT,r
integer::i,maxI
real*8,allocatable::X(:),qT(:),X1(:),X2(:),X3(:),X4(:)

call artemide_Initialize('const-Bench_LL-',prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/DY_Benchmark/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_TMDR((/2.86041d0, 0.229551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.1d0, 0.1d0, 0d0, 0d0, 0d0, 0.d0, 0.d0/))

! call artemide_SetNPparameters_TMDR((/2d0,0.001d0/))
! call artemide_SetNPparameters_uTMDPDF((/0.01d0, 0.01d0, 0d0, 0d0, 0d0, 0.d0, 0.d0/))

proc=(/1,1,5/) !!!! pp DY
s=13000d0**2
Q=91.15348061918276d0
y=0d0

call TMDX_DY_setProcess(proc)
call TMDX_DY_XSetup(s,Q,y)
call TMDX_DY_SetCuts(.false.,0d0,-10d0,10d0)!!! no cuts

!maxqT=10d0
maxqT=100d0
stepqT=0.5d0
maxI=Int(maxqT/stepqT)

allocate(qT(1:maxI))
allocate(X(1:maxI))

do i=1,maxI
  qT(i)=stepqT*i
end do

! call CalcXsec_DY(X,qT)

! do i=1,maxI
!     X(i)=TMDF_F(0d0,qT(i),0d0,0d0,0d0,0d0,0d0,0)
! end do

! do i=1,maxI
!     write(*,'("{",F6.2,",",F14.12,"},")') qT(i),2d0*Q*2d0*qT(i)*X(i)
! !     write(*,'("{",F6.2,",",F14.12,"},")') qT(i),TMDF_F(0d0,qT(i),0d0,0d0,0d0,0d0,0d0,0)
! end do

end program example