!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDR
use TMDR_model
implicit none

integer::i,iMax,f
real*8::bMax,bStep,TT
real*8,allocatable::b(:)

call artemide_Initialize('const-TMDR',prefix='Prog/Tests/const-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_TMDR((/2d0, 0.04d0/))

bMax=4d0
bStep=0.1d0
iMax=int(bMax/bStep)

allocate(b(1:iMax))

do i=1,iMax
   b(i)=bStep*i
end do

do i=1,iMax

    TT=DNP(2d0,b(i),1)
    !TT=zetaNP(91d0,b(i),1)
    !TT=TMDR_Rzeta(b(i),91d0,91d0**2,1)
    write(*,'("{",F8.4," ,",F12.6,"},")') b(i), TT
end do


!!!!!! simple plot of uTMDPDF at fixed b vs. x
! do i=1,100
!
!     TT=uTMDPDF_lowScale50(10d0**(-i/20.),0.3d0,1)
!     write(*,*) "{",10d0**(-i/20.),",", TT(-3),",", TT(-2),",", TT(-1),",", TT(0),",", TT(1),",", TT(2),",", TT(3),"},"
!     !write(*,*) "{",10d0**(-i/20.),",",TT(0), ",", TT(1),",", TT(2),"},"
! end do




end program example
