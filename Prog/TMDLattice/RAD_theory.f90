program example
use aTMDe_control
use TMDR
use QCDinput
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),mean(:),deviation(:),scaleVariationUP(:),scaleVariationDOWN(:),b(:)
real*8::bMax,rr

real*8::PV1,PV2,mu0,Pp1,Pp2

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters_TMDR(NParray(1:2))


!!!---------------------------------------------------------------
!!!
!!! PV=(p.v)
!!! Pp=p^+
!!!
PV1=1.38375d0 !GeV
PV2=0.461251d0 !GeV
Pp1=2.26819d0 !GeV
Pp2=1.25397d0 !GeV

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
allocate(mean(1:numB))
allocate(deviation(1:numB))
allocate(scaleVariationUP(1:numB))
allocate(scaleVariationDOWN(1:numB))

!-------- central
do i=1,numB
  central(i)=Ratio1(b(i))
  mean(i)=0d0
  deviation=0d0
end do


!------- statistical uncertanty
do j=1,numR
  call artemide_GetReplicaFromFile(repFILE,j,NParray)
  call artemide_SetNPparameters_TMDR(NParray(1:2))
  do i=1,numB
    rr=Ratio1(b(i))
    mean(i)=mean(i)+rr
    deviation(i)=deviation(i)+rr**2
  end do  
  
end do

 do i=1,numB
  mean(i)=mean(i)/numR
  deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
 end do

!------- Scale variation
call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters_TMDR(NParray(1:2))
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
	  b(i),mean(i),mean(i)-deviation(i),mean(i)+deviation(i),scaleVariationDOWN(i),scaleVariationUP(i)
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
