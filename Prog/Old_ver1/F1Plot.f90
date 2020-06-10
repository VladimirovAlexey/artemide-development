!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Program that make a plot of unpolarized TMD F1
!!  for given x, mu, zeta in the range b in (0,..)
!!  saves into the file as (b, F(-6:6_)
!!
!!			AV: 16.09.2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use TMDs
implicit none
 character(*),parameter::name="Fx01_2_2.dat"
 integer :: j
 real*8 :: time1, time2
 
 real*8 :: mu, zeta,x,b
 
 mu=2d0
 zeta=mu**2
 x=0.1d0
 
 
 call cpu_time(time1)
 
 call TMDs_Initialize('NNLO')
 call TMDs_SetNPParameters((/2.8736d0, 0.0543d0, 2.7856d0,  1.4201d0, 31.9944d0,1d0, 0d0/))
 
 open(7, file=name, status="replace", action="write")
 
 write(7,*) 'unpolarized TMDPDF at x =', x, '(mu,zeta)=',mu,zeta
 
 do j=0,480
  b=real(j)/80d0
  write(7,*) b,uTMDPDF_5(x,b,mu,zeta,1)
 
 end do
 
 
close(7, STATUS='KEEP')
 
 
 call cpu_time(time2)

 write(*,*) 'timing=',time2-time1
 
end program example