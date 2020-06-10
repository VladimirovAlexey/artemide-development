program xSec_SIDIS
use TMDX_SIDIS
implicit none

!!! this program evaluate cross-section for SIDIS

real*8 :: time1, time2
integer :: dummyINT,j
real*8 :: dummyREAL
logical::exist

integer:: process, order
real*8:: Q,s,x,z
real*8,allocatable, dimension(:) :: pt, xSec


call cpu_time(time1)

process=100!!p+pbar

!!!!! ORDER =0,1,2,3
order=3!!NNLL


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

     
s=8000d0**2
x=0.1d0
z=0.1d0
Q=91d0
     
   call TMDX_SIDIS_XSetup(s,Q,x,z,1,100)
   call InitializeScaleVariations(1d0,1d0,1d0,1d0)  
   
   
   !!!!THESE are parameters for MODEL1 (PDF=FF)
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
   
   
   
   allocate(pt(1:5))
   allocate(xSec(1:5))
   
   do j=1,5
   
   pt(j)=0.2d0*Q*z*j/5d0
   
   end do
   
   call CalcXsec_SIDIS(xSec,pt)
   
   write(*,*) 'Q=', Q, 'x=',x, 'z=',z
   write(*,*) '  pt 					 Fuu'
   do j=1,5
   write(*,*) pt(j),xSec(j)
   end do
call cpu_time(time2)
! 
write(*,*) 'Calculation time=',time2-time1
! 
   
 end program xSec_SIDIS  