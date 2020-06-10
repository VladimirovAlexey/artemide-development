program example
use aTMDe_control
use TMDR
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'Constants-files/DY-NNLO_spline/const-NNPDF_NNLO'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),mean(:),deviation(:),b(:)
real*8::bMax,step,mu,dd

call artemide_Initialize(constFILE)

!numR=artemide_NumOfReplicasInFile(repFILE)

!call artemide_GetReplicaFromFile(repFILE,0,NParray)
!call artemide_SetNPparameters_TMDR(NParray(1:2))

!call artemide_SetNPparameters_TMDR((/-0.00117d0,0.0136d0,0.0764d0/))
call artemide_SetNPparameters_TMDR((/-0.00102d0,0.01499d0,0.051257d0,0.25796d0/))

numB=40
bMax=5d0
!mu=91d0
mu=2.39d0

allocate(b(1:numB))
do i=1,20
  b(i)=0.005d0+i/20d0
end do
do i=21,numB
  b(i)=bMax+(bMax-1d0)*(i-numB)/(numB-20)
end do

allocate(central(1:numB))

do i=1,numB
  central(i)=DNP(mu, b(i),1)
end do


do i=1,numB
  write(*,"('{',F10.7,',',F10.7,'},')") b(i),central(i)
end do


!  


end program example
