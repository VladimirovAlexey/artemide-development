

program example
use TMDX_SIDIS
implicit none
  
  integer::i
  real*8::x,z,s,Q,dz,dx,dQ
  real*8::pT(1:6)
  real*8::xSec(1:5)
  integer::proc(1:3)
  real*8::ss(1:5),ppt(1:5,1:2),xx(1:5,1:2),zz(1:5,1:2),QQ(1:5,1:2)
  integer::pproc(1:5,1:3)
  logical::docut,dodocuts(1:5)
  real*8::cuts(1:3),ccuts(1:5,1:3)
  
  real*8::www
!   
  call TMDX_SIDIS_Initialize("LO")
  call TMDX_SIDIS_SetNPParameters((/3.29829d0, 0.0236075d0, 0.257986d0, 8.2054d0, 301.141d0, 2.43837d0, -4.78378d0,0d0,   &
	0.0000d0,0.02d0,0.2d0/))
!   
!   call TMDX_SIDIS_SetProcess((/1,1,2002/))
!   
!   call TMDX_SIDIS_XSetup(52d0,0.2d0,0.25d0,2.5d0)
!   call TMDX_SIDIS_SetCuts(.false.,0.1d0,0.85d0,10d0)
!   
!   call CalcXsec_SIDIS(www,0.15d0)
!   write(*,*) www
!   
!   call xSec_SIDIS(www,(/1,1,2002/),52d0,(/0.0001d0,0.15d0/),(/0.2d0,0.25d0/),&
! 	  (/0.2d0,0.35d0/),(/2.d0,3.d0/),.False.,(/0.1d0,0.85d0,10d0/))
!   write(*,*) 'xSec_SIDIS(false)',www,www/(0.25d0-0.2d0)/(0.15d0**2-0.0001d0**2)
!   
!   call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(www,0.0001d0,0.15d0,0.2d0,0.25d0,0.2d0,0.35d0,2.d0,3.d0)
!   write(*,*) 'CalcXsec_SIDIS_PTint_Zint_Xint_Qint(false)',www,www/(0.25d0-0.2d0)/(0.15d0**2-0.0001d0**2)
! !   
! !   
!   call TMDX_SIDIS_XSetup(52d0,0.3d0,0.2d0,2.5d0)
!   call TMDX_SIDIS_SetCuts(.true.,0.1d0,0.85d0,10d0)
!   
!   call xSec_SIDIS(www,(/1,1,2002/),52d0,(/0.0001d0,0.15d0/),(/0.2d0,0.25d0/),&
! 	  (/0.2d0,0.35d0/),(/2.d0,3.d0/),.True.,(/0.1d0,0.85d0,10d0/))
!   write(*,*) 'xSec_SIDIS(true)',www,www/(0.25d0-0.2d0)/(0.15d0**2-0.0001d0**2)
!   
!   call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(www,0.0001d0,0.15d0,0.2d0,0.25d0,0.2d0,0.35d0,2.d0,3.d0)
!    write(*,*) 'CalcXsec_SIDIS_PTint_Zint_Xint_Qint(true)',www,www/(0.25d0-0.2d0)/(0.15d0**2-0.0001d0**2)
!   stop
  x=0.2d0
  dx=0.0003d0
  z=0.25d0
  dz=0.0001d0
  s=2d0*160d0*0.938d0+0.938d0**2
  Q=4d0
  dQ=0.002d0  
  proc=(/1,1,2002/)
  docut=.true.
  cuts=(/0.01d0,0.9d0,10d0/)
  
  call TMDX_SIDIS_XSetup(s,z,x,Q)
  call TMDX_SIDIS_SetProcess(proc)
  call TMDX_SIDIS_SetCuts(docut,cuts(1),cuts(2),cuts(3))
  
  do i=1,6
    pT(i)=0.01d0+(i-1)*0.1d0
  end do
  
  do i=1,5
  call CalcXsec_SIDIS(xSec(i),pt(i+1))
  write(*,*) 'pt=',pt(i+1),'X=',xSec(i)
  end do
!   
  call CalcXsec_SIDIS(xSec,pt(2:6))
  write(*,*) '---------------------'
  do i=1,5  
  write(*,*) 'pt=',pt(i+1),'X=',xSec(i)
  end do
  write(*,*) '---------------------'
  write(*,*) '---------------------'
  write(*,*) '---------------------'
  
  do i=1,5
  call CalcXsec_SIDIS_Zint_Xint_Qint(xSec(i),pt(i+1),z-dz,z+dz,x-dx,x+dx,Q-dQ,Q+dQ)
  write(*,*) 'pt=',pt(i+1),'X=',xSec(i)/dx/dz/(2d0*Q*dQ)/8d0
  end do
  
  write(*,*) 'q--q'
  
  call CalcXsec_SIDIS_Zint_Xint_Qint(xSec,pt(2:6),z-dz,z+dz,x-dx,x+dx,Q-dQ,Q+dQ)
  write(*,*) '---------------------'
  do i=1,5  
  write(*,*) 'pt=',pt(i+1),'X=',xSec(i)/dx/dz/(2d0*Q*dQ)/8d0
  end do
  
  write(*,*) '---------------------'
  write(*,*) '---------------------'
  write(*,*) '---------------------'
  
  do i=1,5
  call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(xSec(i),pt(i),pt(i+1),z-dz,z+dz,x-dx,x+dx,Q-dQ,Q+dQ)
  write(*,*) 'pt=',pt(i+1),'X=',xSec(i)/dx/dz/(2d0*Q*dQ)/8d0/(pt(i+1)**2-pt(i)**2)
  end do
  
  call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(xSec,pt(1:5),pt(2:6),z-dz,z+dz,x-dx,x+dx,Q-dQ,Q+dQ)
  write(*,*) '---------------------'
  do i=1,5  
  write(*,*) 'pt=',pt(i+1),'X=',xSec(i)/dx/dz/(2d0*Q*dQ)/8d0/(pt(i+1)**2-pt(i)**2)
  end do
  
  call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(xSec,pt,z-dz,z+dz,x-dx,x+dx,Q-dQ,Q+dQ)
  write(*,*) '---------------------'
  do i=1,5  
  write(*,*) 'pt=',pt(i+1),'X=',xSec(i)/dx/dz/(2d0*Q*dQ)/8d0/(pt(i+1)**2-pt(i)**2)
  end do
  
  write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  do i=1,5
  call xSec_SIDIS(xSec(i),proc,s,(/pt(i),pt(i+1)/),(/z-dz,z+dz/),(/x-dx,x+dx/),(/Q-dQ,Q+dQ/),docut,cuts)
  write(*,*) 'pt=',pt(i+1),'X=',xSec(i)/dx/dz/(2d0*Q*dQ)/8d0/(pt(i+1)**2-pt(i)**2)
  end do
  write(*,*) '+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
  do i=1,5
   ss(i)=s
   pproc(i,1:3)=proc
   ppt(i,1:2)=(/pt(i),pt(i+1)/)
   xx(i,1:2)=(/x-dx,x+dx/)
   QQ(i,1:2)=(/Q-dQ,Q+dQ/)
   zz(i,1:2)=(/z-dz,z+dz/)
   dodocuts(i)=docut
   ccuts(i,1:3)=cuts
  end do
  
  call xSec_SIDIS_list(xSec,pproc,ss,ppt,zz,xx,QQ,dodocuts,ccuts)
  do i=1,5
  write(*,*) 'pt=',pt(i+1),'X=',xSec(i)/dx/dz/(2d0*Q*dQ)/8d0/(pt(i+1)**2-pt(i)**2)
  end do
  
  
  call TMDX_SIDIS_ShowStatistic()


end program example