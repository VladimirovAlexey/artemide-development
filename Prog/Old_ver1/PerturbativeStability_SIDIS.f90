program PertrubativeStability
use TMDX_SIDIS
implicit none

real*8 :: time1, time2
integer :: dummyINT,j,i
real*8 :: dummyREAL
logical::exist

integer:: process, order,varyMOD,length,Xtype,qqq
real*8:: Q,s,x,z,Qmin,Qmax,zMin,zMax,xMin,xMax
real*8:: pt,eta1,eta2,err1,err2,dpt
real*8,allocatable, dimension(:) ::qtList,qt1,qt2,qt3,xSec,xSec1,xSec2,xSec3

real*8::dS1,dS2,dS3,p1,p2,p3,XXX,minV,maxV,dsigmaMIN,dsigmaMAX

real*8,dimension(1:8) :: errorList

call cpu_time(time1)

process=100!!p+pbar
order=3!!NNLL
s=8000d0**2
x=0.1d0
z=0.1d0
Q=91d0

    SELECT CASE(order)
      CASE (0)
	call TMDX_SIDIS_Initialize("NLL","LO")
      CASE (1)
	call TMDX_SIDIS_Initialize("NLO","NLO")
      CASE (2)
	call TMDX_SIDIS_Initialize("NNLL","NLO")
      CASE (3)
	call TMDX_SIDIS_Initialize("NNLO","NNLO")
      CASE DEFAULT
	write(*,*) 'UNKNOWN ORDER SWITCH TO NLO'
	call TMDX_SIDIS_Initialize("NLO","NLO")
     END SELECT

   call TMDX_SIDIS_XSetup(s,Q,x,z,1,100)
   call InitializeScaleVariations(1d0,1d0,1d0,1d0)  
   
   
   !!!!THESE are parameters for MODEL1
     SELECT CASE(order)
      CASE (1)
	call TMDX_SIDIS_SetNPParameters(1d0,0.0231d0,(/0.1894d0,0d0,0.8500d0,0.1894d0,0d0,0.8500d0/))
      CASE (2)
	call TMDX_SIDIS_SetNPParameters(1d0,0.0127d0,(/0.1747d0,0d0,1.0646d0,0.1747d0,0d0,1.0646d0/))
      CASE (3)
	call TMDX_SIDIS_SetNPParameters(1d0,0.0073d0,(/0.2283d0,0d0,0.6128d0,0.2283d0,0d0,0.6128d0/))
      CASE DEFAULT
	call TMDX_SIDIS_SetNPParameters(1d0,0.0231d0,(/0.1894d0,0d0,0.8500d0,0.1894d0,0d0,0.8500d0/))
     END SELECT
   
   
   
   allocate(qt1(1:2))
   allocate(xSec1(1:2))
   allocate(qt2(1:2))
   allocate(xSec2(1:2))
   allocate(qt3(1:2))
   allocate(xSec3(1:2))
   
   allocate(xSec(1:1))
   
    inquire(file="Logs/PTstable/output.dat", exist=exist)
    if (exist) then
    open(7, file="Logs/PTstable/output.dat", status="replace", action="write")
    else
    open(7, file="Logs/PTstable/output.dat", status="new", action="write")
    end if
   
   
   
   x=0.1d0
   z=0.5d0
   Q=8d0
   Qmin=1.5d0
   Qmax=50d0
   
   zMin=0.1d0
   zMax=0.9d0
   
   xMin=0.01d0
   xMax=0.9d0
   
   write(7,*) "Fixed Q=",q ,"out put as (Q,x,z,pT,sMin,sMax)"
   
   do j=0,50
   
!    Q=Qmin*(Qmax/Qmin)**(j/100d0)
   
   z=(zMax-zMin)*(j/50d0)+zMin
   
   do i=0,50
   
   x=(xMax-xMin)*(i/50d0)+xMin
   
   call InitializeScaleVariations(1d0,1d0,1d0,1d0)
   
   
   call TMDX_SIDIS_XSetup(s,Q,x,z,1,100)
   
   !search for maximum of cross-section
   !initial values
   dpt=Q*0.0000001d0
   p1=0.01d0
   p2=0.2d0*Q
   
   qt1=(/p1,p1+dpt/)   
   call CalcXsec_SIDIS(xSec1,qt1)
   !search for maximum
   qt2=(/p2,p2+dpt/)
   call CalcXsec_SIDIS(xSec2,qt2)
   
   dS1=(qt1(2)*xSec1(2)-qt1(1)*xSec1(1))/dpt !!derivatives
   dS2=(qt2(2)*xSec2(2)-qt2(1)*xSec2(1))/dpt
   
   if(dS2>0) then
   write(*,*) 'Peak is away of TMD factorization region', z*0.2d0*Q
   p3=0.2d0*Q
   dsigmaMIN=-5d0
   dsigmaMAX=+5d0
   goto 100
   end if
   
!    write(*,*) p1,p2,dS1,dS2
   
   do
   
   p3=(p1+p2)/2
   qt3=(/p3,p3+dpt/)   
   call CalcXsec_SIDIS(xSec3,qt3)
   dS3=(qt3(2)*xSec3(2)-qt3(1)*xSec3(1))/dpt
   
   if(dS3>0) then
   dS1=dS3
   p1=p3
   else
   dS2=dS3
   p2=p3
   end if
   
!    write(*,*) p1,p2,dS1,dS2
   
   if((p2-p1)<Q*0.00001d0) exit
   
   end do
   !!! p3 is MAXIMUM!!!
   
!    Evaluate xSec and its errors   
   call CalcXsec_SIDIS(xSec,(/p3/))
   
   
   XXX=xSec(1)
   
   call InitializeScaleVariations(0.5d0,1d0,1d0,1d0)
   call CalcXsec_SIDIS(xSec,(/p3/))
   errorList(1)=xSec(1)
   
   call InitializeScaleVariations(2.d0,1d0,1d0,1d0)
   call CalcXsec_SIDIS(xSec,(/p3/))
   errorList(2)=xSec(1)
     
   call InitializeScaleVariations(1d0,0.5d0,1d0,1d0)
   call CalcXsec_SIDIS(xSec,(/p3/))
   errorList(3)=xSec(1)
   
   call InitializeScaleVariations(1d0,2d0,1d0,1d0)
   call CalcXsec_SIDIS(xSec,(/p3/))
   errorList(4)=xSec(1)
   
   call InitializeScaleVariations(1d0,1d0,0.5d0,1d0)
   call CalcXsec_SIDIS(xSec,(/p3/))
   errorList(5)=xSec(1)
   
   call InitializeScaleVariations(1d0,1d0,2d0,1d0)
   call CalcXsec_SIDIS(xSec,(/p3/))
   errorList(6)=xSec(1)
   
   call InitializeScaleVariations(1d0,1d0,1d0,0.5d0)
   call CalcXsec_SIDIS(xSec,(/p3/))
   errorList(7)=xSec(1)
   
   call InitializeScaleVariations(1d0,1d0,1d0,2d0)
   call CalcXsec_SIDIS(xSec,(/p3/))
   errorList(8)=xSec(1)
   
   minV=MINVAL(errorList)
   maxV=MAXVAL(errorList)
   
   dsigmaMIN=-ABS((minV-XXX)/XXX)
   dsigmaMAX=ABS((maxV-XXX)/XXX)
   
100   write(7,*) Q,x,z,p3,dsigmaMIN,dsigmaMAX
!   write(*,*) Q,x,z,p3,dsigmaMIN,dsigmaMAX
!    do j=1,20
!    p3=(Q*0.2d0)*j/40d0
!    call CalcXsec_DY(xSec,(/p3/))
!    write(*,*) p3*xSec,p3
   end do
   write(*,*) j,"/",50
   end do
  CLOSE (7, STATUS='KEEP')
  write(*,*) 'Evaluation finished. The result saved to output.dat'
call cpu_time(time2)
! 
write(*,*) 'Calculation time=',time2-time1
! 
   
 end program PertrubativeStability  