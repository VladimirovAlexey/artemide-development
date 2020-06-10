program xSec_DY
use TMDX_DY
implicit none

!!! this program evaluate cross-section for SIDIS

real*8 :: time1, time2
integer :: dummyINT,j
real*8 :: dummyREAL
logical::exist

integer:: process
real*8:: Qmin,Qmax,s,ymin,ymax
real*8,allocatable, dimension(:) :: pt, xSec1,xSec2
real*8,dimension(1:7)::A,B,C

A=(/2.8736d0, 0.0543d0, 2.7856d0,  1.4201d0, 31.9944d0,0d0, 0d0/)
B=(/2.8736d0, 0.0543d0, 1.5d0,  1.8d0, 1.9944d0,0d0, 0d0/)
C=(/2.8736d0, 0.0543d0, 3d0,  0.2d0, 50d0,0.4d0, 2.5d0/)

call cpu_time(time1)

process=1!!


   call TMDX_DY_Initialize("LO")
   call TMDX_DY_SetNPParameters(A)

     
  s=8000d0**2
  Qmin=66d0
  Qmax=116d0
  ymin=-2.4d0
  ymax=2.4d0
     
   call TMDX_DY_SetProcess(process)
   call SetCuts(.true.,20d0,-2.4d0,2.4d0)
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

  call TMDX_DY_SetNPParameters(B)
  
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