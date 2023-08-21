!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
!use aTMDe_control
use lpTMDPDF
use QCDinput
implicit none

integer::i,iMax
real*8::bMax,bStep,x
real*8,allocatable::b(:)
real*8::TT(-5:5)

!call artemide_Initialize('wgtTMDPDF.atmde',prefix='Prog/Tests/const-files/')
call lpTMDPDF_Initialize('lpTMDPDF.atmde',prefix='Prog/Tests/const-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call lpTMDPDF_SetLambdaNP((/0.0d0, 1d0/))
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

do i=1,100

    TT=lpTMDPDF_lowScale5(i*0.01d0,0.1d0,1)
    !TT=x_hPDF(i*0.01d0,11.d0,1)
    write(*,'(A,F12.10,A,F12.10,A)') "{",i*0.01d0,",", TT(0),"},"
end do


end program example
