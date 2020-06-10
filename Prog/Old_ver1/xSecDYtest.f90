program xSec_DY
use TMDX_DY
use TMDR
use uTMDPDF
implicit none

!!! this program evaluate cross-section for SIDIS

real*8 :: time1, time2
integer :: dummyINT,j
real*8 :: dummyREAL
logical::exist

integer:: process
real*8:: Qmin,Qmax,s,ymin,ymax
real*8,allocatable, dimension(:) :: pt, xSec1,xSec2

real*8:: A(1:2),B(1:7)
A=(/2.5000d0,0.0376d0/)
B=(/0.257986d0, 8.2054d0, 301.141d0, 2.43837d0, -4.78378d0,0d0,0d0/)
call cpu_time(time1)

process=1!!


   call TMDX_DY_Initialize('const-test')
   
   
   call TMDR_setNPparameters(A)
   call uTMDPDF_SetLambdaNP(B)

     
  s=8000d0**2
  Qmin=66d0
  Qmax=116d0
  ymin=-2.4d0
  ymax=2.4d0
     
   call TMDX_DY_SetProcess(process)
   call TMDX_DY_SetCuts(.true.,20d0,-2.4d0,2.4d0)
   call TMDX_DY_XSetup(s,91d0,0d0)
   
   
   allocate(pt(1:6))
   allocate(xSec1(1:5))
   allocate(xSec2(1:5))
   
   do j=1,6
   
   pt(j)=0.2d0+j
   
   end do
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec1,pt,Qmin,Qmax,ymin,ymax)
   
   write(*,*) Qmin,'<Q<', Qmax,'  ',ymin,'<y<',ymax
   write(*,*) '  pt 					 xSec'
   do j=1,5
   write(*,*) pt(j),'-',pt(j+1),xSec1(j)
   end do

  A=(/1.5000d0,0.0376d0/)
  B=(/0.257986d0, 8.2054d0, 301.141d0, 2.43837d0, -4.78378d0,0d0,0d0/)
  call TMDR_setNPparameters(A)
  call uTMDPDF_SetLambdaNP(B)
  
  call CalcXsec_DY_PTint_Qint_Yint(xSec2,pt,Qmin,Qmax,ymin,ymax)
  
     write(*,*) Qmin,'<Q<', Qmax,'  ',ymin,'<y<',ymax
   write(*,*) '  pt 					 xSec'
   do j=1,5
   write(*,*) pt(j),'-',pt(j+1),xSec1(j),xSec2(j)
   end do
  

call cpu_time(time2)
! 
write(*,*) 'Calculation time=',time2-time1
! 
   
 end program xSec_DY