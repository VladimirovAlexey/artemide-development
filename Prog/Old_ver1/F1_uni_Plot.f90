!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Program that make a plot of unpolarized TMD F1 at universal line
!!  for given x, mu, zeta in the range b in (0,..)
!!  saves into the file as (b, F(-6:6_)
!!
!!			AV: 16.09.2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDPDF
implicit none
 character(*),parameter::name="Fx0001_uni.dat"
 integer :: j
 real*8 :: time1, time2
 
 real*8 :: x,b
 
 x=0.001d0
 
 
 call cpu_time(time1)
 
 call uTMDPDF_Initialize('NNLO')
 call uTMDPDF_SetLambdaNP((/2.7856d0,  1.4201d0, 31.9944d0,1d0, 0d0/))
 
 open(7, file=name, status="replace", action="write")
 
 write(7,*) 'unpolarized TMDPDF at x =', x, '(mu,zeta)= universal'
 
 do j=0,480
  b=real(j)/80d0
  write(7,*) b,uTMDPDF_lowScale5(x,b,1)
 
 end do
 
 
close(7, STATUS='KEEP')
 
 
 call cpu_time(time2)

 write(*,*) 'timing=',time2-time1
 
end program example