!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
!use aTMDe_control
use wglTMDPDF
use QCDinput
implicit none

integer::i,iMax
real*8::bMax,bStep,x
real*8,allocatable::b(:)
real*8::TT(-5:5)

!call artemide_Initialize('wgtTMDPDF.atmde',prefix='Prog/Tests/const-files/')
call wglTMDPDF_Initialize('wglTMDPDF_v3.atmde',prefix='Prog/Tests/const-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call wglTMDPDF_SetLambdaNP((/0.0d0, 0.d0/))
x=0.1

bMax=5d0
bStep=0.1d0
iMax=int(bMax/bStep)

allocate(b(1:iMax))

do i=1,iMax
  b(i)=bStep*i
end do

bMax=0.5d0


write(*,*) '---------------------- D -----------------------'

do i=1,100

     x=i*0.01d0
    TT=wglTMDPDF_inB(x,bMax,1)
    !TT=x_hPDF(i*0.01d0,11.d0,1)
    write(*,'("{",F12.8,",",F14.10,"},")',advance='no') x, x*TT(1)
end do
write(*,*)

write(*,*) '---------------------- U -----------------------'

do i=1,100

     x=i*0.01d0
    TT=wglTMDPDF_inB(x,bMax,1)
    !TT=x_hPDF(i*0.01d0,11.d0,1)
    write(*,'("{",F12.8,",",F14.10,"},")',advance='no') x, x*TT(2)
end do
write(*,*)

! write(*,*) '---------------------- S -----------------------'
!
! do i=1,100
!
!      x=i*0.01d0
!     TT=wglTMDPDF_inB(x,bMax,1)
!     !TT=x_hPDF(i*0.01d0,11.d0,1)
!     write(*,'("{",F12.8,",",F14.10,"},")',advance='no') x, x*TT(3)
! end do
! write(*,*)

end program example
