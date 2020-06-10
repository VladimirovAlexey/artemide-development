program example
use aTMDe_control
use TMDR
use QCDinput
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

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'

integer::numB,i,j
real*8,allocatable::central(:),RADv1(:),RADv2(:),RADv3(:),RADv4(:),b(:),scaleVariationUP(:),scaleVariationDOWN(:)
real*8::bMax,rr

real*8::PV1,PV2,mu0,Pp1,Pp2

call artemide_Initialize(constFILE)

call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0/))


!!!---------------------------------------------------------------
!!!
!!! PV=(p.v)
!!! Pp=p^+
!!!
PV1=1.38375d0 !GeV
PV2=0.922501d0 !GeV
Pp1=2.26819d0 !GeV
Pp2=1.7362d0 !GeV

mu0=(Abs(PV1)+Abs(PV2))/2d0
!!!---------------------------------------------------------------

numB=60
bMax=6d0

allocate(b(1:numB))
do i=1,20
  b(i)=0.005d0+i/20d0
end do
do i=21,numB
  b(i)=bMax+(bMax-1d0)*(i-numB)/(numB-20)
end do

allocate(central(1:numB))
allocate(RADv1(1:numB))
allocate(RADv2(1:numB))
allocate(RADv3(1:numB))
allocate(RADv4(1:numB))
allocate(scaleVariationUP(1:numB))
allocate(scaleVariationDOWN(1:numB))



!!!! Bacchetta,et al 1912.07550
do i=1,numB
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0/))
  central(i)=Ratio1(b(i))
  call artemide_SetNPparameters_TMDR((/0.036d0+0.009d0,0.012d0/))
  RADv1(i)=Ratio1(b(i))
  call artemide_SetNPparameters_TMDR((/0.036d0-0.009d0,0.012d0/))
  RADv2(i)=Ratio1(b(i))
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0+0.003d0/))
  RADv3(i)=Ratio1(b(i))
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0-0.003d0/))
  RADv4(i)=Ratio1(b(i))
end do

!------- Scale variation
call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0/))
!! UP
mu0=2d0*(Abs(PV1)+Abs(PV2))
do i=1,numB
  scaleVariationUP(i)=Ratio1(b(i))
end do
!! DOWN
mu0=(Abs(PV1)+Abs(PV2))/4d0
if(mu0<0.5d0) mu0=0.5d0
do i=1,numB
  scaleVariationDOWN(i)=Ratio1(b(i))
end do

do i=1,numB
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7',',F10.7,',',F10.7,'},')") &
	  b(i),central(i),min(RADv1(i),RADv2(i),RADv3(i),RADv4(i)),max(RADv1(i),RADv2(i),RADv3(i),RADv4(i)),&
	  scaleVariationDOWN(i),scaleVariationUP(i)
end do


contains

!!! Hard coefficient
function CH2(pv,mu)
  real*8::pv,mu,CH2
  real*8::LL
  LL=log((pv)**2/mu**2)
  
  CH2=1d0+As(mu)*4d0/3d0*(-LL**2/2+LL-2d0+0.822467d0)
  
end function CH2

!!! ratio F(pv1,pp1)/F(pv2,pp2) with all the rest being same
!!! =|CH(pv1/m0)|^2(pp1^2/zeta0)^[-D(b,mu0)]/|CH(pv2/m0)|^2(pp2^2/zeta0)^[-D(b,mu0)]
function Ratio1(b_in)
  real*8::Ratio1,b_in
  Ratio1=CH2(PV1,mu0)/CH2(PV2,mu0)*(Abs(pp1/pp2))**(-2d0*DNP(mu0, b_in,1))
end function Ratio1

end program example