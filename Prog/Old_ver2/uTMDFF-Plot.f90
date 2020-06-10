!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Program that make a plot of unpolarized TMD F1
!!  for given x, mu, zeta in the range b in (0,..)
!!  saves into the file as (b, F(-6:6_)
!!
!!			AV: 16.09.2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDFF
implicit none
 character(*),parameter::name="Fx01_2_2.dat"
 integer :: i,j
 real*8 :: time1, time2
 
 real*8 :: mu, zeta,z,b,TMD(-5:5)
 
 call uTMDFF_Initialize('const-DY+SIDIS_NNLO-DSS',&
      '/misc/data2/braun/vla18041/arTeMiDe_Repository/Constants-files/DY+SIDIS-NNLO/')
 
 
 call cpu_time(time1)
 
 call uTMDFF_Initialize('NNLO')
 call uTMDFF_SetLambdaNP((/0.2631d0,0.4786d0,0.4598d0,0.5400d0/),.false.,.false.)
 
 do i=0,16
  b=real(i)/8d0+0.001d0
 do j=0,25
  z=0.2d0+0.6d0*real(j)/25d0
  TMD=uTMDFF_lowScale5(z,b,1)
  write(*,"('{',F10.7,',',F10.7,',',F10.7,'},')") z,b,TMD(1)
 end do
 end do
 do i=1,6
  b=real(i)/2d0+2.00d0
 do j=0,25
  z=0.2d0+0.6d0*real(j)/25d0
  TMD=uTMDFF_lowScale5(z,b,1)
  write(*,"('{',F10.7,',',F10.7,',',F10.7,'},')") z,b,TMD(1)
 end do
 end do
 
!  open(7, file=name, status="replace", action="write")
!  
!  write(7,*) 'unpolarized TMDPDF at x =', x, '(mu,zeta)=',mu,zeta
!  
!  do j=0,480
!   b=real(j)/80d0
!   write(7,*) b,uTMDPDF_5(x,b,mu,zeta,1)
!  
!  end do
!  
!  
! close(7, STATUS='KEEP')
!  
 
 call cpu_time(time2)

 write(*,*) 'timing=',time2-time1
 
end program example