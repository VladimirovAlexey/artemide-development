program example
use uTMDPDF
use TMDs_inKT
implicit none
  
  real*8,dimension(-5:5)::F1,F2
  integer::i
  real*8::kt,x
  
  call TMDs_inKT_Initialize('const-HERAPDF20_NNLO+pi','/home/vla18041/LinkData2/arTeMiDe_Repository/')
  !call uTMDPDF_Initialize('const-HERAPDF20_NNLO+pi','/home/vla18041/LinkData2/arTeMiDe_Repository/')
  
  call uTMDPDF_SetLambdaNP(0)
  
  x=0.3d0
  
  do i=1,20
  
   kt=i*0.05
   
   F1=uTMDPDF_kT_5(x,kt,1)
   F2=uTMDPDF_kT_5(x,kt,2)
   !F1=MeanTMD(x,kt,2)
   !F2=StdTMD(x,kt,2)
    write(*,*) '{',kt,',',F1(1),',',F2(1),'},'
   
  end do
   do i=1,8
  
   kt=1d0+i*0.5
   
   F1=uTMDPDF_kT_5(x,kt,1)
   F2=uTMDPDF_kT_5(x,kt,2)
   !F1=MeanTMD(x,kt,2)
   !F2=StdTMD(x,kt,2)
    write(*,*) '{',kt,',',F1(1),',',F2(1),'},'
   
  end do
  
  contains
  
  
  !!! evaluetes mean value over replicas
  function MeanTMD(x,kt,h)
  real*8,dimension(-5:5)::MeanTMD
  real*8::x,kt
  integer::i,h
  
  MeanTMD=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  do i=1,100
   call uTMDPDF_SetLambdaNP(i)
   MeanTMD=MeanTMD+uTMDPDF_kT_5(x,kt,h)
  end do
  MeanTMD=MeanTMD/100d0
  
  end function
  
  
  !!! evaluetes standard deviation over replicas
  function StdTMD(x,kt,h)
  real*8,dimension(-5:5)::StdTMD,F,MM
  real*8::x,kt
  integer::i,h
  
  MM=MeanTMD(x,kt,h)
  StdTMD=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  do i=1,100
   call uTMDPDF_SetLambdaNP(i)
   F=uTMDPDF_kT_5(x,kt,h)
   StdTMD=StdTMD+(F-MM)**2
  end do
  StdTMD=sqrt(StdTMD/100d0)
  
  end function
  
end program example