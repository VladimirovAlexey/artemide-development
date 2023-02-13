!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDR
use TMDR_model
use TMD_AD
implicit none

integer::i,iMax,f
real*8::bMax,bStep
real*8,allocatable::b(:),TT0(:),TT1(:)

!call artemide_Initialize('const-TMDR2',prefix='Prog/Tests/const-files/')
call artemide_Initialize('ForScale_N4LL',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_TMDR((/1.48571d0, 0.0461374d0, 0.0706658d0,1d0/))

bMax=2.6d0
bStep=0.1d0
iMax=int(bMax/bStep)

allocate(b(-3:iMax))
allocate(TT0(-3:iMax))
allocate(TT1(-3:iMax))

do i=-3,iMax
   b(i)=bStep*i
end do

!call artemide_SetNPparameters_TMDR((/1.5d0, 0.05d0, 0.07d0,1.d0/))
do i=-3,iMax
    TT0(i)=zetaSL(100.d0,b(i),1)
end do

write(*,*) DNP(100.d0,2.d0,1)

do i=-3,iMax
    write(*,'("{",F8.4," ,",F12.8,"},")',advance="no") b(i), TT0(i)
end do

write(*,*) " "
end program example
