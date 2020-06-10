program example
use TMDX_SIDIS
use TMDF
use QCDinput
implicit none
  
  real*8,parameter::hc2=0.389379338d0*1d9
  real*8::x=0.254d0
  
  real*8::Q,Q2
  real*8::s=52.65
  real*8::y,z,qT
  real*8::fff(-5:5),xDIS,sigma0
  
  real*8::W,U
  
  Q2=5.241270d0
  Q=sqrt(Q2)
  y=Q2/s/x
  
  sigma0=2d0*3.14159d0*(1d0/129d0)**2*(1+(1-y)**2)/Q2**2
  
  call TMDX_SIDIS_Initialize("LO")
  call TMDX_SIDIS_SetNPParameters((/3.29829d0, 0.0236075d0, 0.257986d0, 8.2054d0, 301.141d0, 2.43837d0, -4.78378d0,0d0,0d0,&
	0.02d0,0.01d0,0.144d0,0.3d00/))
  
  fff=xPDF(x,Q,1)/x
  
  xDIS=fff(1)/9.d0+fff(2)*4.d0/9.d0+fff(3)/9.d0+fff(4)*4d0/9.d0+fff(5)/9d0+&
  fff(-1)/9.d0+fff(-2)*4.d0/9.d0+fff(-3)/9.d0+fff(-4)*4d0/9.d0+fff(-5)/9d0
  
  
  
  write(*,*) '--------------'
  
  write(*,*) 'y=',y, 'Nabuo=',x*xDIS
  write(*,*) sigma0, xDIS,sigma0*xDIS,sigma0*xDIS*hc2
  
  write(*,*) '------------'
  
  z=0.23d0
  qT=0.1d0
  
  W=TMDF_F(Q2,qT/z,x,z,Q,Q2,Q2,2001)
  
  write(*,*) W,W*sigma0,W*sigma0*hc2
  write(*,*) W/xDIS,2*qT*W/xDIS
  
  
  write(*,*) '-------------------------'
  
  call TMDX_SIDIS_XSetup(s,z,x,Q,0d0,0d0)
  call TMDX_SIDIS_SetProcess((/1,1,2001/))
  call CalcXsec_SIDIS(U,qT)
  
  write(*,*) U,2*qT*U/(sigma0*xDIS*hc2)

end program example