program example
use aTMDe_control
use TMDR
use TMDs
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
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),b(:),c1(:),c2(:),c3(:),c4(:)
real*8::bMax,rr,bb,xCUT

real*8::Pp1,Pp2,mu0,vminus

!!! if true then x->x(1+M^2/x^2 mu^2)
logical,parameter::xMODIFIED=.false.

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters(NParray)
call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0/))


vminus=1d0/sqrt(2d0)
! Pp1=2.26819d0 !GeV
! Pp2=1.25397d0 !GeV

! Pp1=10d0 !GeV
! Pp2=9d0 !GeV

Pp1=2.5d0 !GeV
Pp2=2d0 !GeV

xCUT=1.0d-3
! xCUT=1.0d-4
! xCUT=1.01d-5

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
allocate(c1(1:numB))
allocate(c2(1:numB))
allocate(c3(1:numB))
allocate(c4(1:numB))
! 

xCUT=1.0d-3
write(*,*) '-------------------------------------------------------------'
!-------- central
do i=1,numB
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0/))
  central(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0+0.009d0,0.012d0/))
  c1(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0-0.009d0,0.012d0/))
  c2(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0+0.003d0/))
  c3(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0-0.003d0/))
  c4(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
end do

 do i=1,numB
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") &
	  b(i),central(i),min(central(i),c1(i),c2(i),c3(i),c4(i)),max(central(i),c1(i),c2(i),c3(i),c4(i))
 end do

 xCUT=1.0d-4
write(*,*) '-------------------------------------------------------------'
!-------- central
do i=1,numB
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0/))
  central(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0+0.009d0,0.012d0/))
  c1(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0-0.009d0,0.012d0/))
  c2(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0+0.003d0/))
  c3(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0-0.003d0/))
  c4(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
end do

 do i=1,numB
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") &
	  b(i),central(i),min(central(i),c1(i),c2(i),c3(i),c4(i)),max(central(i),c1(i),c2(i),c3(i),c4(i))
 end do
 
 xCUT=1.01d-5
write(*,*) '-------------------------------------------------------------'
!-------- central
do i=1,numB
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0/))
  central(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0+0.009d0,0.012d0/))
  c1(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0-0.009d0,0.012d0/))
  c2(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0+0.003d0/))
  c3(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
  call artemide_SetNPparameters_TMDR((/0.036d0,0.012d0-0.003d0/))
  c4(i)=R1(b(i),mu0,Pp1/Pp2,1,xCUT)
end do

 do i=1,numB
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") &
	  b(i),central(i),min(central(i),c1(i),c2(i),c3(i),c4(i)),max(central(i),c1(i),c2(i),c3(i),c4(i))
 end do
 
contains

!!!! function is the THe factor RC
function R1(b,mu,PP,f,xMax)
 real*8::R1,b,mu,PP,xMax,L
 integer::f
  
 R1=PP**(-1-2d0*DNP(mu,b,f))*RC(b,mu,PP,f,xMax)
end function R1

!!!! function is the THe factor RC
function RC(b,mu,PP,f,xMax)
 real*8::RC,b,mu,PP,xMax,L
 integer::f
  
 RC=1+16d0/3d0*as(mu)*Log(PP)*(1-2d0*M(b,mu,f,xMax))
end function RC


!!!! function is the ration of \int dx f(x)|x|^{-2D} TMD /\int dx |x|^{-2D} TMD (for flavour f)
function M(b,mu,f,xMax)
 real*8::M,b,mu,xMax,L
 integer::f
  
 L=-Log(xMax)
 M=T2(b,mu,f,L)/T1(b,mu,f,L)
end function M

!!!! integral of I1 from 0 to L
!!!! simpson with devision by nMax (even)
function T1(b,mu,f,L)
  real*8::T1,b,mu,L
  integer::f,i
  real*8:: step
  integer,parameter:: nMAX=200
  
  step=L/nMAX
  
  T1=0d0
  !! i use that x=1 (t=0) the function is null
  T1=T1+4d0*I1(step,b,mu,f)
  do i=2,nMax-2,2
    T1=T1+2d0*I1(step*i,b,mu,f)
    T1=T1+4d0*I1(step*(i+1),b,mu,f)
  end do
  T1=T1+I1(step*nMax,b,mu,f)
  
  T1=step/3d0*T1
  
end function T1

!!!! integral of I1 from 0 to L
!!!! simpson with devision by nMax (even)
function T2(b,mu,f,L)
  real*8::T2,b,mu,L
  integer::f,i
  real*8:: step
  integer,parameter:: nMAX=200
  
  step=L/nMAX
  
  
  T2=0d0
  !! i use that x=1 (t=0) the function is null
  T2=T2+4d0*I1(step,b,mu,f)
  do i=2,nMax-2,2
    T2=T2+2d0*I2(step*i,b,mu,f)
    T2=T2+4d0*I2(step*(i+1),b,mu,f)
  end do
  T2=T2+I2(step*nMax,b,mu,f)
  
  T2=step/3d0*T2
  
end function T2

!!! -Exp[-t(1-2D)]TMD[Exp[-t])
!!! correspond to Integrand of log|x| moment with x=exp[-t]
function I2(t,b,mu,f)
  real*8::I2,t,b,mu
  integer::f
  real*8::TTT(-5:5),t0
  TTT=uTMDPDF_5(exp(-t),b,1)
  
  if(xMODIFIED) then
    t0=t*Log(1+(0.95d0/mu)**2*exp(2d0*t))
  else
    t0=t
  end if
  
  I2=-t*Exp(-t+2d0*t0*DNP(mu,b,f))*(TTT(f)-TTT(-f))
  
end function I2

!!! -t Exp[-t(1-2D)]TMD[Exp[-t])
!!! correspond to Integrand of 1 moment with x=exp[-t]
function I1(t,b,mu,f)
  real*8::I1,t,b,mu
  integer::f
  real*8::TTT(-5:5),t0
  TTT=uTMDPDF_5(exp(-t),b,1)
  
  if(xMODIFIED) then
    t0=t*Log(1+(0.95d0/mu)**2*exp(2d0*t))
  else
    t0=t
  end if
  
  I1=Exp(-t+2d0*t0*DNP(mu,b,f))*(TTT(f)-TTT(-f))
  
end function I1

end program example
