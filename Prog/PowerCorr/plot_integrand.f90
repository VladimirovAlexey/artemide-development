program example
use aTMDe_control
use TMDF
use TMDR
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=&
'/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Prog/PowerCorr/const-DY+SIDIS_NNPDF31+DSS_nnlo'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::b(:)
real*8::bMax,Q2,x1,x2
integer::process

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters(NParray)

numB=150
bMax=10d0





allocate(b(1:numB))
do i=1,100
  b(i)=i/100d0
end do
do i=1,numB-100
  b(i+100)=1.00d0+i/5d0
end do


!!!! CDF 1
Q2=91d0**2
x1=0.05d0
x2=0.05d0
process=6
do i=1,numB
!!!! D vs. b[GeV]
  write(*,"('{',F16.7,',',F16.7,'},')") &
       b(i),Integrand(Q2,b(i),x1,x2,sqrt(Q2),Q2,Q2,process)!TMDR_Rzeta(b(i),sqrt(Q2),Q2,1)
end do

write(*,*) '-------------------------------------------------------------'

!!!! HERMES
Q2=3d0**2
x1=0.1d0
x2=0.4d0
process=2001

do i=1,numB
!!!! D vs. b[GeV]
  write(*,"('{',F16.7,',',F16.7,'},')") &
       b(i),Integrand(Q2,b(i),x1,x2,sqrt(Q2),Q2,Q2,process)!TMDR_Rzeta(b(i),sqrt(Q2),Q2,1)
end do



end program example
