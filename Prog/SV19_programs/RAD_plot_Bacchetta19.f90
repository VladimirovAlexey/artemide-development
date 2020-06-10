program example
use aTMDe_control
use TMDR
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! I set the following model
! ! ! function DNP(mu,b,f)
! ! !  real*8::DNP,mu,b
! ! !  integer::f
! ! !   bSTAR=1.123d0*((1d0-exp(-b**4/(1.123d0)**4))/(1d0-exp(-b**4/(1.123d0/mu)**4)))**(0.25d0)
! ! !   DNP=Dpert(mu,bSTAR,1)+(NPparam(1)+NPparam(2)*b**2)*b**2/4
! ! !  end function DNP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! HERE 4 = mu (??)
!! with NNLO order of evolution. 
!! Since correlation between these parameter is not extreme I just vary each of them

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/'
character(*),parameter::constFILE=prefix//'const-DY_LO'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),centralUP1(:),centralUP2(:),centralDOWN1(:),centralDOWN2(:),b(:)
real*8::bMax,step,mu,dd

call artemide_Initialize(constFILE)

!!!! see table X
call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0/))

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
allocate(centralUP1(1:numB))
allocate(centralDOWN1(1:numB))
allocate(centralUP2(1:numB))
allocate(centralDOWN2(1:numB))

!!!! Bacchetta,et al 1912.07550
do i=1,numB
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0/))
  central(i)=DNP(mu, b(i),1)
  call artemide_SetNPparameters_TMDR((/0.036d0+0.009d0,0.012d0/))
  centralUP1(i)=DNP(mu, b(i),1)
  call artemide_SetNPparameters_TMDR((/0.036d0-0.009d0,0.012d0/))
  centralDOWN1(i)=DNP(mu, b(i),1)
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0+0.003d0/))
  centralUP2(i)=DNP(mu, b(i),1)
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0-0.003d0/))
  centralDOWN2(i)=DNP(mu, b(i),1)
end do

do i=1,numB
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") b(i),central(i),&
	min(centralUP1(i),centralDOWN1(i),centralUP2(i),centralDOWN2(i)),&
	max(centralUP1(i),centralDOWN1(i),centralUP2(i),centralDOWN2(i))
end do


!  


end program example