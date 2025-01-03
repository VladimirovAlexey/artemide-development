!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
!use aTMDe_control
use wgtTMDPDF
use QCDinput
implicit none

integer::i,iMax
real*8::bMax,bStep,x
real*8,allocatable::b(:)
real*8::TT(-5:5)

!call artemide_Initialize('wgtTMDPDF.atmde',prefix='Prog/Tests/const-files/')
call wgtTMDPDF_Initialize('wgtTMDPDF_v3.atmde',prefix='Prog/Tests/const-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call wgtTMDPDF_SetLambdaNP((/0.0d0, 0.d0/))
x=0.1

bMax=5d0
bStep=0.1d0
iMax=int(bMax/bStep)

allocate(b(1:iMax))

do i=1,iMax
  b(i)=bStep*i
end do

! do i=1,iMax
! 
!     TT=wgtTMDPDF_lowScale5(x,b(i),1)
!     write(*,*) "{",b(i),",", TT(1),"},"
! end do

! do i=1,iMax
! 
!     TT=wgtTMDPDF_lowScale5(10**(-b(i)),0.2d0,1)
!     write(*,*) "{",-b(i),",", TT(1),"},"
! end do

write(*,*) '---------------------- D -----------------------'

bMax=0.05d0

! do i=1,100
!     x=i*0.01d0
!     TT=wgtTMDPDF_inB(x,bMax,1)
!     !TT=x_hPDF(i*0.01d0,11.d0,1)
!     write(*,'(A,F12.10,A,F12.10,A)') "{",x,",", x*TT(1),"},"
! end do

write(*,*) '---------------------- U -----------------------'

do i=1,25

     x=i*0.04d0
    TT=wgtTMDPDF_inB(x,bMax,1)
    !TT=x_hPDF(i*0.01d0,11.d0,1)
    write(*,'(A,F12.10,A,F12.10,A)') "{",x,",", x*TT(2),"},"
end do

write(*,*) '---------------------- S -----------------------'

! do i=1,100
!
!      x=i*0.01d0
!     TT=wgtTMDPDF_inB(x,bMax,1)
!     !TT=x_hPDF(i*0.01d0,11.d0,1)
!     write(*,'(A,F12.10,A,F12.10,A)') "{",x,",", x*TT(3),"},"
! end do

end program example
