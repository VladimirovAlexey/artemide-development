program PertrubativeStability
use TMDX_DY
use TMDR
implicit none

real*8 :: time1, time2
integer :: dummyINT,j,jMax,i,iMax
real*8 :: dummyREAL
logical::exist

integer:: process, order,varyMOD,length,Xtype,qqq
real*8:: Q,s,y,Qmin,Qmax,yMin,yMax,xMin,xMax
real*8:: pt,eta1,eta2,err1,err2,dpt,x
real*8,allocatable, dimension(:) ::qtList,qt1,qt2,qt3,xSec,xSec1,xSec2,xSec3

real*8::dS1,dS2,dS3,p1,p2,p3,XXX,minV,maxV

real*8,dimension(1:8) :: errorList
real*8,dimension(1:3) ::lowQ

call cpu_time(time1)

process=7!!p+pbar
order=1
y=0d0
s=8000d0**2
Q=91d0

    SELECT CASE(order)
      CASE (1)
	call TMDX_DY_Initialize("LO")
      CASE (2)
	call TMDX_DY_Initialize("NLO")
      CASE (3)
	call TMDX_DY_Initialize("NNLO")
      CASE DEFAULT
	write(*,*) 'UNKNOWN ORDER SWITCH TO NLO'
	call TMDX_DY_Initialize("NLO")
     END SELECT
     
   call SetCuts(.false.,0d0,-1d0,1d0)
   call TMDX_DY_XSetup(s,Q,y,process)
   call InitializeScaleVariations(1d0,1d0,1d0,1d0)  
   
   
   !!!!THESE are parameters for MODEL1
     SELECT CASE(order)
      CASE (1)
	call TMDX_DY_SetNPParameters((/1.67d0,0.327d0,0.112d0,0.828d0/))
      CASE (2)
	call TMDX_DY_SetNPParameters((/2.43d0,0.106d0,0.183d0,0.328d0/))
      CASE (3)
	call TMDX_DY_SetNPParameters((/2.45d0,0.0087d0,0.246d0,0.307d0/))
      CASE DEFAULT
	call TMDX_DY_SetNPParameters((/2.43d0,0.106d0,0.183d0,0.328d0/))
     END SELECT
   
   
   lowQ=LowestQ()
   write(*,*) 'Lowest Q',lowQ
   
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
   
   
   Qmin=1.5d0
   Qmax=120d0
   xMin=0.01d0
   xMax=0.9d0
   x=0.1d0
   Q=5d0
   s=19.4d0**2!(Q/x)**2   
   ymin=1/2d0*LOG(Q**2/s)+0.000000001d0
   ymax=-1/2d0*LOG(Q**2/s)-0.000000001d0
   y=0d0
   
   write(7,*) "Fixed s=",s, "Fixed x=",x
   write(7,*) 'Lowest Q ={',lowQ,'}'
   
   Qmin=lowQ(3)+0.1d0
   
   jMax=50
!    iMax=50
   do j=0,jmax   
   
   Q=Qmin*(Qmax/Qmin)**(j/Real(jMax))
   
   write(*,*) 'Q=',Q
   
!    do i=0,iMax
!    x=(xmax-xmin)*i/Real(iMax)+xmin
   s=(Q/x)**2  
!    y=(ymax-ymin)*j/Real(jMax)+ymin
   
   call InitializeScaleVariations(1d0,1d0,1d0,1d0)
   
   
   call TMDX_DY_XSetup(s,Q,y,process)
   
   !search for maximum of cross-section
   !initial values
   dpt=Q*0.0000001d0
   p1=0.01d0
   p2=0.1d0*Q
   
   qt1=(/p1,p1+dpt/)   
   call CalcXsec_DY(xSec1,qt1)
   !search for maximum
   qt2=(/p2,p2+dpt/)
   call CalcXsec_DY(xSec2,qt2)
   
   dS1=(qt1(2)*xSec1(2)-qt1(1)*xSec1(1))/dpt !!derivatives
   dS2=(qt2(2)*xSec2(2)-qt2(1)*xSec2(1))/dpt
   
!    write(*,*) p1,p2,dS1,dS2
   
   do
   
   p3=(p1+p2)/2
   qt3=(/p3,p3+dpt/)   
   call CalcXsec_DY(xSec3,qt3)
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
   call CalcXsec_DY(xSec,(/p3/))
   
   XXX=xSec(1)
   
!    call InitializeScaleVariations(0.5d0,1d0,1d0,1d0)
!    call CalcXsec_DY(xSec,(/p3/))
!    errorList(1)=xSec(1)
   
!    call InitializeScaleVariations(2.d0,1d0,1d0,1d0)
!    call CalcXsec_DY(xSec,(/p3/))
!    errorList(2)=xSec(1)
     
   call InitializeScaleVariations(1d0,0.5d0,1d0,1d0)
   call CalcXsec_DY(xSec,(/p3/))
   errorList(3)=xSec(1)
   
   call InitializeScaleVariations(1d0,2d0,1d0,1d0)
   call CalcXsec_DY(xSec,(/p3/))
   errorList(4)=xSec(1)
   
!    call InitializeScaleVariations(1d0,1d0,0.5d0,1d0)
!    call CalcXsec_DY(xSec,(/p3/))
!    errorList(5)=xSec(1)
   
!    call InitializeScaleVariations(1d0,1d0,2d0,1d0)
!    call CalcXsec_DY(xSec,(/p3/))
!    errorList(6)=xSec(1)
   
   call InitializeScaleVariations(1d0,1d0,1d0,0.5d0)
   call CalcXsec_DY(xSec,(/p3/))
   errorList(7)=xSec(1)
   
   call InitializeScaleVariations(1d0,1d0,1d0,2d0)
   call CalcXsec_DY(xSec,(/p3/))
   errorList(8)=xSec(1)
   
   minV=MINVAL((/errorList(3),errorList(4),errorList(7),errorList(8)/))
   maxV=MAXVAL((/errorList(3),errorList(4),errorList(7),errorList(8)/))
!    minV=MINVAL(errorList)
!    maxV=MAXVAL(errorList)
   
!    write(*,*) Sqrt(s),Q,p3,(minV-XXX)/XXX,(maxV-XXX)/XXX!XXX,minV,maxV
!    write(7,*) Sqrt(s),Q,p3,(minV-XXX)/XXX,(maxV-XXX)/XXX!XXX,minV,maxV
!    write(*,*) Sqrt(s),Q,x,(minV-XXX)/XXX,(maxV-XXX)/XXX!XXX,minV,maxV
   write(7,*) Sqrt(s),Q,x,(minV-XXX)/XXX,(maxV-XXX)/XXX!XXX,minV,maxV
!    do j=1,20
!    p3=(Q*0.2d0)*j/40d0
!    call CalcXsec_DY(xSec,(/p3/))
!    write(*,*) p3*xSec,p3
!    end do
!    end do
   write(*,*) j,'/',jMax
   end do
  CLOSE (7, STATUS='KEEP')
  write(*,*) 'Evaluation finished. The result saved to output.dat'
call cpu_time(time2)
! 
write(*,*) 'Calculation time=',time2-time1
! 
   
 end program PertrubativeStability  