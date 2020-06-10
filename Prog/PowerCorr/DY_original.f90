!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_DY
implicit none

integer,dimension(1:3)::proc
real*8::y,Q2,s
integer::i
real*8,allocatable::X(:),qT(:)

call artemide_Initialize('const-TMD-inKT_NNLO',prefix='Prog/PowerCorr/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

! Q2=91d0**2
! s=1800d0**2
s=27.43d0**2
Q2=7d0**2
y=0d0

proc=(/1,-10221191,2/) !!!! pp DY

call TMDX_DY_setProcess(proc)
call TMDX_DY_XSetup(s,sqrt(Q2),y)
call TMDX_DY_SetCuts(.false.,0d0,-10d0,10d0)!!! no cuts


allocate(qT(1:81))
allocate(X(1:81))

do i=0,80
  qT(i+1)=i/80d0*Sqrt(Q2)/2d0+0.0001d0
end do

call CalcXsec_DY(X,qT)

do i=1,81
    write(*,'("{",F8.3,",",F14.12,"},")') qT(i),X(i)
end do

end program example
