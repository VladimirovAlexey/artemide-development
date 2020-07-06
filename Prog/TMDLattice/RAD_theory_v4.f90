program example
use aTMDe_control
use TMDR
use TMDR_model
use QCDinput
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),mean(:),deviation(:),b(:)
real*8::bMax,rr

real*8::P1,P2,mu0

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters_TMDR(NParray(1:2))


!!!---------------------------------------------------------------
!!!
!!! Momenta
!!! case 1
! P1=2.26818d0 !GeV
! P2=1.25397d0 !GeV
!!! case 2
! P1=2.26818d0 !GeV
! P2=1.7362d0 !GeV
!!! case 3
P1=1.7362d0 !GeV
P2=1.25397d0 !GeV

mu0=sqrt(2*P1*P2)
!!!---------------------------------------------------------------

numB=60
bMax=4.5d0

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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PLOT EVOLUTION DIFFERENCE
!!!! evolution difference = DNP(mu)-DNP(mu0) to mu=2GeV
!call artemide_SetNPparameters_TMDR((/500d0,0d0/))
! do i=1,numB
!   write(*,"('{',F10.7,',',F10.7,'},')") &
! 	  !b(i),DNP(mu0,b(i),1)-DNP(2d0,b(i),1)
! 	  b(i),-2d0*DNP(2d0,b(i),1)
! end do
! stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PLOT RATIO1
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


do i=1,numB
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") &
	  b(i),mean(i),mean(i)-deviation(i),mean(i)+deviation(i)
end do


contains

!!! ratio F(pv1,pp1)/F(pv2,pp2) with all the rest being same
!!! =|CH(pv1/m0)|^2(pp1^2/zeta0)^[-D(b,mu0)]/|CH(pv2/m0)|^2(pp2^2/zeta0)^[-D(b,mu0)]
function Ratio1(b_in)
  real*8::Ratio1,b_in
  Ratio1=-2d0*DNP(2d0,b_in,1)!(P2/P1)**(2d0*DNP(mu0,b_in,1))
end function Ratio1

end program example
