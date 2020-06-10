program xSec_AN
use TMDX_DY
implicit none

!!! this program evaluate cross-section for SIDIS

real*8 :: time1, time2
integer :: dummyINT,j
real*8 :: dummyREAL
logical::exist

integer:: process, order
real*8:: Q,s,y
real*8,allocatable, dimension(:) :: ptMin,ptMax, xSec,xSec2


call cpu_time(time1)

process=5!!p+pbar

!!!!! ORDER =0,1,2,3
order=0


    SELECT CASE(order)
      CASE (0)
	call TMDX_DY_Initialize("LO")
      CASE (1)
	call TMDX_DY_Initialize("NLO")
      CASE (2)
	call TMDX_DY_Initialize("NNLO")
      CASE DEFAULT
	write(*,*) 'UNKNOWN ORDER SWITCH TO NLO'
	call TMDX_DY_Initialize("NLO")
     END SELECT

     
s=510d0**2
Q=91d0
y=0d0
     
   call SetCuts(.false.,0d0,-1d0,1d0)
   call TMDX_DY_XSetup(s,Q,y,process)
   call InitializeScaleVariations(1d0,1d0,1d0,1d0) 
   
   
   !!!!THESE are parameters for MODEL1
     SELECT CASE(order)
      CASE (0)
	call TMDX_DY_SetNPParameters((/1.67d0,0.327d0,0.112d0,0.828d0/))
      CASE (1)
	call TMDX_DY_SetNPParameters((/2.43d0,0.106d0,0.183d0,0.328d0/))
      CASE (2)
	call TMDX_DY_SetNPParameters((/2.45d0,0.0087d0,0.246d0,0.307d0/))
      CASE DEFAULT
	call TMDX_DY_SetNPParameters((/2.43d0,0.106d0,0.183d0,0.328d0/))
     END SELECT
   
   
   
   allocate(xSec(10))
   allocate(xSec2(10))
   allocate(ptMin(10))
   allocate(ptMax(10))
   
   ptMin=(/0.5d0,1.d0,1.5d0,2.5d0,3.5d0,4.5d0,5.5d0,6.5d0,7.5d0,8.5d0/)
   ptMax=(/1.d0,1.5d0,2.5d0,3.5d0,4.5d0,5.5d0,6.5d0,7.5d0,8.5d0,10d0/)
   
   do j=0,4
   
   y=-1d0+j/2d0
    
    
   call TMDX_DY_XSetup(s,Q,y,5)
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec,ptMin,ptMax,70d0,110d0,y-0.25d0,y+0.25d0)
   
   call TMDX_DY_XSetup(s,Q,y,100)
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec2,ptMin,ptMax,70d0,110d0,y-0.25d0,y+0.25d0)
   
   write(*,*) 'y=',y,Sum(xSec2),Sum(xSec),Sum(xSec2)/Sum(xSec)
   end do
call cpu_time(time2)
! 
write(*,*) 'Calculation time=',time2-time1
! 
   
 end program xSec_AN 