program example
use aTMDe_control
use TMDR
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! I set the following model
! ! ! function DNP(mu,b,f)
! ! !  real*8::DNP,mu,b
! ! !  integer::f
! ! !   DNP=Dpert(mu,b,1)+NPparam(2)*b**2/2
! ! !  end function DNP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! with LO+ order of evolution. It seems to agree with their (33) and (32).

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! with the same program one can compute RAD from our 1706.01473
! ! ! DNP=Dpert(mu,b,1)+NPparam(2)*b*b
!! with NNLO order of evolution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/'
character(*),parameter::constFILE='/home/alexey/artemide_Repository/artemide/Prog/SV19_programs/const-Pavia19'
!character(*),parameter::constFILE=prefix//'const-DY_NNLO'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),centralUP(:),centralDOWN(:),b(:)
real*8::bMax,step,mu,dd

call artemide_Initialize(constFILE)

!!!! see table X
call artemide_SetNPparameters_TMDR((/2d0,0.13d0/))

numB=40
bMax=5d0
mu=2d0

allocate(b(1:numB))
do i=1,20
  b(i)=0.005d0+i/20d0
end do
do i=21,numB
  b(i)=bMax+(bMax-1d0)*(i-numB)/(numB-20)
end do

allocate(central(1:numB))
allocate(centralUP(1:numB))
allocate(centralDOWN(1:numB))

!!!! Bacchetta,et al
do i=1,numB
  call artemide_SetNPparameters_TMDR((/2d0,0.13d0/))
  central(i)=DNP(mu, b(i),1)
  call artemide_SetNPparameters_TMDR((/2d0,0.13d0+0.01d0/))
  centralUP(i)=DNP(mu, b(i),1)
  call artemide_SetNPparameters_TMDR((/2d0,0.13d0-0.01d0/))
  centralDOWN(i)=DNP(mu, b(i),1)
end do



!!!! 1706.01473 (model 1)
! do i=1,numB
!   call artemide_SetNPparameters_TMDR((/2d0,0.0073d0/))
!   central(i)=DNP(mu, b(i),1)
!   call artemide_SetNPparameters_TMDR((/2d0,0.0073d0+0.0024d0/))
!   centralUP(i)=DNP(mu, b(i),1)
!   call artemide_SetNPparameters_TMDR((/2d0,0.0073d0-0.0023d0/))
!   centralDOWN(i)=DNP(mu, b(i),1)
! end do

! do i=1,numB
!   call artemide_SetNPparameters_TMDR((/2d0,0.396753d0/))
!   central(i)=DNP(mu, b(i),1)
!   call artemide_SetNPparameters_TMDR((/2d0,0.396753d0+0.003178d0/))
!   centralUP(i)=DNP(mu, b(i),1)
!   call artemide_SetNPparameters_TMDR((/2d0,0.396753d0-0.003178d0/))
!   centralDOWN(i)=DNP(mu, b(i),1)
! end do


do i=1,numB
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") b(i),central(i),centralDOWN(i),centralUP(i)
end do


!  


end program example
