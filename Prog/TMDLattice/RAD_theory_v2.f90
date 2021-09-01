program example
use aTMDe_control
use TMDR
use TMDs
use QCDinput
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),mean(:),deviation(:),b(:)
real*8::bMax,rr,bb,xCUT

real*8::Pp1,Pp2,mu0,vminus,F

integer:: lMoment!!! derivative with respect to ell
!!!! consdered u-d combination

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters(NParray)

!!!!! definotion of l-moment start from 0
lMoment=0
vminus=1d0/sqrt(2d0)
Pp1=2.26819d0 !GeV
!Pp2=1.25397d0 !GeV
Pp2=1.73620d0 !GeV


!Pp1=2.5d0 !GeV
!Pp2=2d0 !GeV

xCUT=1.01d-5

mu0=2d0*Abs(vminus)*Sqrt(Pp1*Pp2)

write(*,*) '--------------M-----------------------------------------------'
do i=1,25
  bb=0.2d0*i
  F=M(bb,mu0,1d-4)
  write(*,*) '{',bb,',', F,'},'
  !write(*,*) '{',bb,',', M(bb,mu0,1d-3),',', M(bb,mu0,1d-4),',', M(bb,mu0,1d-5),'},'
end do

stop
! 
! 
! write(*,*) '----------------RC--------------------------------------------'
! do i=1,25
!    bb=0.2d0*i
!    write(*,*) '{',bb,',', RC(bb,mu0,Pp1/Pp2,1d-3),',', RC(bb,mu0,Pp1/Pp2,1d-4),',',RC(bb,mu0,Pp1/Pp2,1d-5),'},'
! end do
! 
write(*,*) '----------------R1---------------------------------------------'
do i=1,25
   bb=0.2d0*i
   write(*,*) '{',bb,',', R1(bb,mu0,Pp1/Pp2,1d-3),',', R1(bb,mu0,Pp1/Pp2,1d-4),',', R1(bb,mu0,Pp1/Pp2,1d-5),'},'
end do

numB=25
bMax=5d0
! 
allocate(b(1:numB))
do i=1,numB
  b(i)=bMax*i/numB
end do
! 
allocate(central(1:numB))
allocate(mean(1:numB))
allocate(deviation(1:numB))
! 

! xCUT=1.0d-3
! write(*,*) '---------------xcut=10^{-3}----------------------------------------------'
! call artemide_GetReplicaFromFile(repFILE,0,NParray)
! call artemide_SetNPparameters(NParray)
! !-------- central
! do i=1,numB
!   central(i)=R1(b(i),mu0,Pp1/Pp2,xCUT)
!   mean(i)=0d0
!   deviation=0d0
! end do
! ! 
! ! 
! !------- statistical uncertanty
! do j=1,numR
!   call artemide_GetReplicaFromFile(repFILE,j,NParray)
!   call artemide_SetNPparameters_TMDR(NParray(1:2))
!   do i=1,numB
!     rr=R1(b(i),mu0,Pp1/Pp2,xCUT)
!     mean(i)=mean(i)+rr
!     deviation(i)=deviation(i)+rr**2
!   end do  
! !   
! end do
! ! 
!  do i=1,numB
!   mean(i)=mean(i)/numR
!   deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
!  end do
! 
!  do i=1,numB
!   write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") &
! 	  b(i),mean(i),mean(i)-deviation(i),mean(i)+deviation(i)
!  end do


 xCUT=1.00d-4
write(*,*) '------------------xcut=10^{-4}----------------------------------------'
call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters(NParray)
!-------- central
do i=1,numB
  central(i)=R1(b(i),mu0,Pp1/Pp2,xCUT)
  mean(i)=0d0
  deviation=0d0
end do
! 
! 
!------- statistical uncertanty
do j=1,numR
  call artemide_GetReplicaFromFile(repFILE,j,NParray)
  call artemide_SetNPparameters_TMDR(NParray(1:2))
  do i=1,numB
    rr=R1(b(i),mu0,Pp1/Pp2,xCUT)
    mean(i)=mean(i)+rr
    deviation(i)=deviation(i)+rr**2
  end do  
!   
end do
! 
 do i=1,numB
  mean(i)=mean(i)/numR
  deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
 end do

 do i=1,numB
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") &
	  b(i),mean(i),mean(i)-deviation(i),mean(i)+deviation(i)
 end do

 
!  xCUT=1.01d-5
! write(*,*) '--------------------xcut=10^{-5}----------------------------------------'
! call artemide_GetReplicaFromFile(repFILE,0,NParray)
! call artemide_SetNPparameters(NParray)
! !-------- central
! do i=1,numB
!   central(i)=R1(b(i),mu0,Pp1/Pp2,xCUT)
!   mean(i)=0d0
!   deviation=0d0
! end do
! ! 
! ! 
! !------- statistical uncertanty
! do j=1,numR
!   call artemide_GetReplicaFromFile(repFILE,j,NParray)
!   call artemide_SetNPparameters_TMDR(NParray(1:2))
!   do i=1,numB
!     rr=R1(b(i),mu0,Pp1/Pp2,xCUT)
!     mean(i)=mean(i)+rr
!     deviation(i)=deviation(i)+rr**2
!   end do  
! !   
! end do
! ! 
!  do i=1,numB
!   mean(i)=mean(i)/numR
!   deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
!  end do
! 
!  do i=1,numB
!   write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") &
! 	  b(i),mean(i),mean(i)-deviation(i),mean(i)+deviation(i)
!  end do

 
contains

!!!! function is the THe factor RC
function R1(b,mu,PP,xMax)
 real*8::R1,b,mu,PP,xMax,L
  
 R1=PP**(-1-2d0*DNP(mu,b,1)+lMoment)*RC(b,mu,PP,xMax)
end function R1

!!!! function is the THe factor RC
function RC(b,mu,PP,xMax)
 real*8::RC,b,mu,PP,xMax,L
  
 RC=1+16d0/3d0*as(mu)*Log(PP)*(1-2d0*M(b,mu,xMax))
end function RC


!!!! function is the ration of \int dx f(x)|x|^{-2D} TMD /\int dx |x|^{-2D} TMD (for flavour f)
function M(b,mu,xMax)
 real*8::M,b,mu,xMax,L,q1,q2
  
  
 L=-Log(xMax) 
 
 q1=T1(b,mu,L)
 q2=T2(b,mu,L)
 
 M=q2/q1
end function M

!!!! integral of I1 from 0 to L
!!!! simpson with devision by nMax (even)
function T1(b,mu,L)
  real*8::T1,b,mu,L
  integer::i
  real*8:: step
  integer,parameter:: nMAX=200
  
  step=L/nMAX
  
  T1=0d0
  !! i use that x=1 (t=0) the function is null
  T1=T1+4d0*I1(step,b,mu)
  do i=2,nMax-2,2
    T1=T1+2d0*I1(step*i,b,mu)
    T1=T1+4d0*I1(step*(i+1),b,mu)
  end do
  T1=T1+I1(step*nMax,b,mu)
  
  !write(*,*) "T1---->",i,"    ",T2
  
  T1=step/3d0*T1
  
end function T1

!!!! integral of I1 from 0 to L
!!!! simpson with devision by nMax (even)
function T2(b,mu,L)
  real*8::T2,b,mu,L
  integer::i
  real*8:: step
  integer,parameter:: nMAX=200
  
  step=L/nMAX
  
  
  T2=0d0
  !! i use that x=1 (t=0) the function is null
  T2=T2+4d0*I1(step,b,mu)
  do i=2,nMax-2,2
    T2=T2+2d0*I2(step*i,b,mu)
    T2=T2+4d0*I2(step*(i+1),b,mu)
    
    !write(*,*) "T2---->",i,"    ",T2
    
  end do
  T2=T2+I2(step*nMax,b,mu)
  
  T2=step/3d0*T2
  
end function T2

!!! -Exp[-t(1-2D)]TMD[Exp[-t])
!!! correspond to Integrand of log|x| moment with x=exp[-t]
function I2(t,b,mu)
  real*8::I2,t,b,mu
  real*8::TTT(-5:5),ud
  TTT=uTMDPDF_5(exp(-t),b,1)
  
  !! u-d combination
  ud=(TTT(2)-TTT(-2))-(TTT(1)-TTT(-1))
    
  I2=-t*Exp(-t+2d0*t*DNP(mu,b,1))*ud*exp(-t*lMoment)
  
end function I2

!!! -t Exp[-t(1-2D)]TMD[Exp[-t])
!!! correspond to Integrand of 1 moment with x=exp[-t]
function I1(t,b,mu)
  real*8::I1,t,b,mu
  real*8::TTT(-5:5),ud
  TTT=uTMDPDF_5(exp(-t),b,1)
  
  ud=(TTT(2)-TTT(-2))-(TTT(1)-TTT(-1))
  
  I1=Exp(-t+2d0*t*DNP(mu,b,1))*ud*exp(-t*lMoment)
  
end function I1

end program example
