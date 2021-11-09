!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use wgtTMDPDF
implicit none

integer::i,iMax
real*8::bMax,bStep,x
real*8,allocatable::b(:)
real*8::TT(-5:5)

call artemide_Initialize('const-test-wgt',prefix='Prog/Tests/const-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_wgtTMDPDF((/0.2d0, 1d0,0.01d0/))

x=0.1

bMax=5d0
bStep=0.1d0
iMax=int(bMax/bStep)

allocate(b(1:iMax))

do i=1,iMax
  b(i)=bStep*i
end do

do i=1,iMax

    TT=wgtTMDPDF_lowScale5(x,b(i),1)
    write(*,*) "{",b(i),",", TT(1),"},"
end do

! do i=1,iMax
! 
!     TT=wgtTMDPDF_lowScale5(10**(-b(i)),0.2d0,1)
!     write(*,*) "{",-b(i),",", TT(1),"},"
! end do

end program example
