program PertrubativeStability
use TMDX_SIDIS
implicit none

real*8 :: time1, time2
integer :: dummyINT,j,i
real*8 :: dummyREAL
logical::exist

real*8::x,y,z,Q
real*8,dimension(1:5)::xSec,pT,norm,F2,F2exp
integer::order,process

call cpu_time(time1)



process=101!!SIDIS on deutron
order=2!!NLO

    SELECT CASE(order)
      CASE (0)
	call TMDX_SIDIS_Initialize("LO")
      CASE (1)
	call TMDX_SIDIS_Initialize("NLO")
      CASE (2)
	call TMDX_SIDIS_Initialize("NNLO")
      CASE DEFAULT
	write(*,*) 'UNKNOWN ORDER SWITCH TO NLO'
	call TMDX_SIDIS_Initialize("NLO")
     END SELECT

   !!!!THESE are parameters for MODEL1
     SELECT CASE(order)
      CASE (0)
	call TMDX_SIDIS_SetNPParameters((/1.67d0,0.327d0,0.112d0,0.828d0,0.112d0,0.828d0/))
      CASE (1)
	call TMDX_SIDIS_SetNPParameters((/2.34d0,0.111d0,0.179d0,0.3541d0,0.2791d0,0.3541d0/))
      CASE (2)
	call TMDX_SIDIS_SetNPParameters((/2.52d0,0.0827d0,0.2461d0,0.3111d0,0.2461d0,0.3111d0/))
      CASE DEFAULT
	call TMDX_SIDIS_SetNPParameters((/2.34d0,0.111d0,0.179d0,0.354d0,0.179d0,0.354d0/))
     END SELECT
   
   x=0.157d0
   Q=Sqrt(20d0)
   y=0.439d0
   z=0.2d0
   
   call TMDX_SIDIS_XSetup(Q,x,y,z,process)
   call InitializeScaleVariations(1d0,1d0,1d0,1d0)  
   
   F2exp=(/6.2719d0,6.2175d0,5.0537d0,4.8854d0,4.1757d0/)
   pT=(/0.3d0,0.33166d0,0.36056d0,0.4d0,0.42426d0/)
   norm=0.64329d0*(/1,1,1,1,1/)
   
   call CalcXsec_SIDIS(xSec,pT)
   
   F2=xSec/norm
   write(*,*) '----------------------------------------'
   write(*,*) 'z*pT/Q', z/Q*pT
   write(*,*) 'theor',F2
   write(*,*) 'exp',F2exp
   write(*,*) '----------------------------------------'

   
   x=0.2909d0
   Q=Sqrt(22.1d0)
   y=0.261d0
   z=0.2d0
   
   call TMDX_SIDIS_XSetup(Q,x,y,z,process)
   call InitializeScaleVariations(1d0,1d0,1d0,1d0)  
   
   F2exp=(/6.9848d0,6.2763d0,6.0334d0,5.2609d0,4.781d0/)
   pT=(/0.3d0,0.33166d0,0.36056d0,0.4d0,0.42426d0/)
   norm=0.44934d0*(/1,1,1,1,1/)
   
   call CalcXsec_SIDIS(xSec,pT)
   
   F2=xSec/norm
   write(*,*) '----------------------------------------'
   write(*,*) 'z*pT/Q', z/Q*pT
   write(*,*) 'theor',F2
   write(*,*) 'exp',F2exp
   write(*,*) '----------------------------------------'
   
   x=0.1488d0
   Q=Sqrt(9.8d0)
   y=0.23d0
   z=0.2d0
   
   call TMDX_SIDIS_XSetup(Q,x,y,z,process)
   call InitializeScaleVariations(1d0,1d0,1d0,1d0)  
   
   F2exp=(/6.3571d0,5.9058d0,5.2267d0,4.5806d0,3.9552d0/)
   pT=(/0.3d0,0.33166d0,0.36056d0,0.4d0,0.42426d0/)
   norm=0.66156d0*(/1,1,1,1,1/)
   
   call CalcXsec_SIDIS(xSec,pT)
   
   F2=xSec/norm
   write(*,*) '----------------------------------------'
   write(*,*) 'z*pT/Q', z/Q*pT
   write(*,*) 'theor',F2
   write(*,*) 'exp',F2exp
   write(*,*) '----------------------------------------'
   
call cpu_time(time2)
write(*,*) 'Calculation time=',time2-time1
! 
   
 end program PertrubativeStability  