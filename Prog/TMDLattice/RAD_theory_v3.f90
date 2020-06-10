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

real*8::Pp1,Pp2,mu0,vminus

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters(NParray)


vminus=1d0/sqrt(2d0)

Pp1=2.5d0 !GeV
Pp2=2d0 !GeV

xCUT=1.01d-5

mu0=2d0*Abs(vminus)*Sqrt(Pp1*Pp2)

! write(*,*) '-------------------------------------------------------------'
! do i=1,25
!   bb=0.2d0*i
!   write(*,*) '{',bb,',', M(bb,mu0,1,1d-3),',', M(bb,mu0,1,1d-4),',', M(bb,mu0,1,1d-5),'},'
! end do


! write(*,*) '-------------------------------------------------------------'
! do i=1,40
!    bb=0.2d0*i
!    write(*,*) '{',bb,',', RC(bb,mu0,Pp1/Pp2,1,1d-3),',', RC(bb,mu0,Pp1/Pp2,1,1d-4),'},'
! end do

! write(*,*) '-------------------------------------------------------------'
! do i=1,25
!    bb=0.2d0*i
!    write(*,*) '{',bb,',', R1(bb,mu0,Pp1/Pp2,1,1d-3),',', R1(bb,mu0,Pp1/Pp2,1,1d-4),',', R1(bb,mu0,Pp1/Pp2,1,1d-5),'},'
! end do

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

xCUT=1.0d-3
write(*,*) '-------------------------------------------------------------'
call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters(NParray)
!-------- central
do i=1,numB
  central(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
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
    rr=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
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

xCUT=1.00d-4
write(*,*) '-------------------------------------------------------------'
call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters(NParray)
!-------- central
do i=1,numB
  central(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
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
    rr=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
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

 
 xCUT=1.01d-5
write(*,*) '-------------------------------------------------------------'
call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters(NParray)
!-------- central
do i=1,numB
  central(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
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
    rr=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
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

 
contains

!!!! function is the THe factor RC
!!!! here mu=v^-sqrt{P1 P2}
!!!! PP=Sqrt{P1/P2}
function R1(b,mu,PP,f,xMax)
 real*8::R1,b,mu,PP,xMax,L
 integer::f
  
 R1=PP**(-1-2d0*DNP(mu,b,f))*RC(b,mu,PP,f,xMax)
end function R1

!!!! function is the ration of \int dx f(x)|x|^{-2D} TMD /\int dx |x|^{-2D} TMD (for flavour f)
function RC(b,mu,PP,f,xMax)
 real*8::RC,b,mu,xMax,L,PP
 integer::f
  
 L=-Log(xMax)
 RC=T1(b,mu,PP,f,L)/T1(b,mu,1d0/PP,f,L)
end function RC

!!!! integral of I1 from 0 to L
!!!! simpson with devision by nMax (even)
function T1(b,mu,PP,f,L)
  real*8::T1,b,mu,L,PP
  integer::f,i
  real*8:: step
  integer,parameter:: nMAX=200
  
  step=L/nMAX
  
  T1=0d0
  !! i use that x=1 (t=0) the function is null
  T1=T1+4d0*I1(step,b,mu,PP,f)
  do i=2,nMax-2,2
    T1=T1+2d0*I1(step*i,b,mu,PP,f)
    T1=T1+4d0*I1(step*(i+1),b,mu,PP,f)
  end do
  T1=T1+I1(step*nMax,b,mu,PP,f)
  
  T1=step/3d0*T1
  
end function T1

!!! |CH|^2 TMD(x) (x)^{-2D}
!!! correspond to Integrand of 1 moment with x=exp[-t]
function I1(t,b,mu,PP,f)
  real*8::I1,t,b,mu,PP
  integer::f
  real*8::TTT(-5:5),CH2,LL
  TTT=uTMDPDF_5(exp(-t),b,1)
  LL=2d0*(Log(2d0*PP)-t)
  
  CH2=1d0+4d0/3d0*as(mu)*(-LL**2+2d0*LL-4d0+1.6449340668482262d0)
  
  I1=Exp(-t+2d0*t*DNP(mu,b,f))*CH2*(TTT(f)-TTT(-f))
  
end function I1

end program example
