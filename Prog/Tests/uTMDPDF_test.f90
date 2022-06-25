!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use uTMDPDF
implicit none

integer::i,iMax
real*8::bMax,bStep,x
real*8,allocatable::b(:)
real*8::TT(-5:5)

call artemide_Initialize('const-uTMDPDF',prefix='Prog/Tests/const-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

! call CheckCoefficient(1d0,1d0,3.d0,0.4d0)
! call CheckCoefficient(1d0,1d0,3.d0,0.3d0)
!   call CheckCoefficient(1d0,1d0,3.d0,0.2d0)
!   call CheckCoefficient(1d0,1d0,3.d0,0.1d0)
!   stop

! x=0.1
!
! bMax=5d0
! bStep=0.1d0
! iMax=int(bMax/bStep)
!
! allocate(b(1:iMax))
!
! do i=1,iMax
!   b(i)=bStep*i
! end do
!
! do i=1,iMax
!
!     TT=uTMDPDF_lowScale5(x,b(i),1)
!     write(*,*) "{",b(i),",", TT(1),"},"
! end do

do i=1,100

    TT=uTMDPDF_lowScale50(10d0**(-i/20.),0.1d0,1)
    write(*,*) "{",10d0**(-i/20.),",", TT(-3),",", TT(-2),",", TT(-1),",", TT(0),",", TT(1),",", TT(2),",", TT(3),"},"
    !write(*,*) "{",10d0**(-i/20.),",",TT(0), ",", TT(1),",", TT(2),"},"
end do

end program example
