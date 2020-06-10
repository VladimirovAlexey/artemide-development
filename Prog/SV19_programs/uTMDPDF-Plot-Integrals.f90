!!!!
!!! Program coputes integrals of TMD in b-space
!!!!
program example
use aTMDe_control
use uTMDPDF
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'

real*8,allocatable::NParray(:)
integer::numR,i,j
real*8::central0,mean0,deviation0,t0
real*8::central2,mean2,deviation2,t2
real*8::central20,mean20,deviation20
real*8::step,mu,dd,x
integer::f,k
!real*8,parameter::xValues(1:7)=(/.99d0,10d0**(-0.5d0),10d0**(-1d0),10d0**(-1.5d0),10d0**(-2d0),10d0**(-2.5d0),10d0**(-3d0)/)
real*8,parameter::xValues(1:22)=(/&
0.95d0,10d0**(-0.06d0),&
10d0**(-0.125d0),10d0**(-0.25d0),10d0**(-0.375d0),10d0**(-0.5d0),10d0**(-0.625d0),10d0**(-0.75d0),10d0**(-0.875d0),&
10d0**(-1d0),10d0**(-1.125d0),10d0**(-1.25d0),10d0**(-1.375d0),10d0**(-1.5d0),10d0**(-1.625d0),10d0**(-1.75d0),10d0**(-1.875d0),&
10d0**(-2d0),10d0**(-2.25d0),10d0**(-2.5d0),10d0**(-2.75d0),10d0**(-3d0)/)

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call uTMDPDF_SetLambdaNP(NParray(3:9),.false.,.false.)


f=-1

mu=4d0


do k=1,size(xValues)
  x=xValues(k)
!   call Moments(x,mu,f,central0,central2)
  mean0=0d0
  deviation0=0d0
  mean2=0d0
  deviation2=0d0
  mean20=0d0
  deviation20=0d0

  do j=1,numR
    call artemide_GetReplicaFromFile(repFILE,j,NParray)
    call uTMDPDF_SetLambdaNP(NParray(3:9),.false.,.false.)
    call Moments(x,mu,f,t0,t2)
    mean0=mean0+t0
    deviation0=deviation0+(t0)**2
    mean2=mean2+t2
    deviation2=deviation2+(t2)**2
    mean20=mean20+t2/t0
    deviation20=deviation20+(t2/t0)**2
    
  end do

  mean0=mean0/numR
  deviation0=Sqrt(deviation0/numR-mean0**2)
  mean2=mean2/numR
  deviation2=Sqrt(deviation2/numR-mean2**2)
  mean20=mean20/numR
  deviation20=Sqrt(deviation20/numR-mean20**2)

    write(*,"('{'F10.7,',',F10.3,',',F10.3,',',F10.3,'},')") x,mean20,&
	  (mean20-deviation20),(mean20+deviation20)
end do


contains

!!! just and integral int_0^infty TMD(b) db
function NormIntegral(x,mu,f)
  real*8::x,mu,NormIntegral
  integer::f
  real*8,parameter::bMax=10d0
  real*8,parameter::step=0.01d0
  integer::i,imax
  real*8::F1(-5:5),bb,R
  
  imax=int(bMax/step)
  
  bb=step
  F1=uTMDPDF_lowScale5(x,bb,1)
  R=F1(f)
  
  do i=2,imax-1
    bb=step*i
    F1=uTMDPDF_lowScale5(x,bb,1)
    R=R+2d0*F1(f)
  end do
  
  bb=step*imax
  F1=uTMDPDF_lowScale5(x,bb,1)
  R=R+F1(f)
  
  NormIntegral=R*step/2d0
  
end function NormIntegral


!! compute moments of TMD M0=<TMD>, M2=<b^2>
!! <TMD>= 2pi\int_0^\infty  b F db
!! <b^2>= 2pi\int_0^\infty  b^3 F db
subroutine Moments(x,mu,f,M0,M2)
  real*8::x,mu,NormIntegral
  real*8,intent(out)::M0,M2
  integer::f
  real*8::bMax
  real*8,parameter::step=0.05d0
  integer::i,imax
  real*8::F1(-5:5),bb,R0,R2
  
  bMax=10d0+50d0*x
  imax=int(bMax/step)
  
  bb=step
  F1=uTMDPDF_lowScale5(x,bb,1)
  R0=F1(f)*bb
  R2=F1(f)*bb**3
  
  do i=2,imax-1
    bb=step*i
    F1=uTMDPDF_lowScale5(x,bb,1)
    R0=R0+2d0*F1(f)*bb
    R2=R2+2d0*F1(f)*bb**3
  end do
  
  bb=step*imax
  F1=uTMDPDF_lowScale5(x,bb,1)
  R0=R0+F1(f)*bb
  R2=R2+F1(f)*bb**3
  
  M0=R0*step*3.14159d0
  M2=R2*step*3.14159d0
  
end subroutine Moments

end program example