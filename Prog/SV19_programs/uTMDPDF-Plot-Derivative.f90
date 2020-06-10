program example
use aTMDe_control
use uTMDPDF
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'

real*8,allocatable::NParray(:)
integer::numR,i,j
real*8::central,mean,deviation,b,eps
real*8::step,mu,dd,x,TMD(-5:5),TMD1(-5:5),TMD2(-5:5),dTMD(-5:5)
integer::f,k
!real*8,parameter::xValues(1:7)=(/.99d0,10d0**(-0.5d0),10d0**(-1d0),10d0**(-1.5d0),10d0**(-2d0),10d0**(-2.5d0),10d0**(-3d0)/)
real*8,parameter::xValues(1:13)=(/&
0.99d0,10d0**(-0.25d0),10d0**(-0.5d0),10d0**(-0.75d0),&
10d0**(-1d0),10d0**(-1.25d0),10d0**(-1.5d0),10d0**(-1.75d0),&
10d0**(-2d0),10d0**(-2.25d0),10d0**(-2.5d0),10d0**(-2.75d0),10d0**(-3d0)/)

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call uTMDPDF_SetLambdaNP(NParray(3:9),.false.,.false.)



f=3

mu=4d0
b=0.3d0
eps=0.01d0


do k=1,size(xValues)
  x=xValues(k)

  TMD=uTMDPDF_lowScale5(x,b,1)
  TMD1=uTMDPDF_lowScale5(x,b+eps,1)
  TMD2=uTMDPDF_lowScale5(x,b-eps,1)
  dTMD=-1d0*((TMD1-TMD2)/(2*eps*b)+(TMD1-2d0*TMD+TMD2)/eps**2)
  central=dTMD(f)/TMD(f)
  mean=0d0
  deviation=0d0

  do j=1,numR
    call artemide_GetReplicaFromFile(repFILE,j,NParray)
    call uTMDPDF_SetLambdaNP(NParray(3:9),.false.,.false.)
    TMD=uTMDPDF_lowScale5(x,b,1)
    TMD1=uTMDPDF_lowScale5(x,b+eps,1)
    TMD2=uTMDPDF_lowScale5(x,b-eps,1)
    dTMD=-1d0*((TMD1-TMD2)/(2*eps*b)+(TMD1-2d0*TMD+TMD2)/eps**2)
    mean=mean+dTMD(f)/TMD(f)
    deviation=deviation+(dTMD(f)/TMD(f))**2
    
  end do

  mean=mean/numR
  deviation=Sqrt(deviation/numR-mean**2)

    write(*,"('{'F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") x,mean,&
	  (mean-deviation),(mean+deviation)

end do


end program example