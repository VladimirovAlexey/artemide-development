program example
use uTMDPDF
implicit none
  
  real*8,dimension(-5:5)::F1,F2
  integer::i
  real*8::b,x
  
  call uTMDPDF_Initialize('const-HERAPDF20_NNLO+pi','/home/vla18041/LinkData2/arTeMiDe_Repository/')
  
  call uTMDPDF_SetLambdaNP(0)
  
  x=0.6d0
  
  do i=1,20
  
   b=i*0.05
   
   F1=MeanTMD(x,b,2)
   F2=StdTMD(x,b,2)
    write(*,*) '{',b,',',F1(1),',',F2(1),'},'
   
  end do
   do i=1,18
  
   b=1d0+i*0.5
   
   F1=MeanTMD(x,b,2)
   F2=StdTMD(x,b,2)
    write(*,*) '{',b,',',F1(1),',',F2(1),'},'
   
  end do
  
  contains
  
  
  !!! evaluetes mean value over replicas
  function MeanTMD(x,b,h)
  real*8,dimension(-5:5)::MeanTMD
  real*8::x,b
  integer::i,h
  
  MeanTMD=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  do i=1,100
   call uTMDPDF_SetLambdaNP(i)
   MeanTMD=MeanTMD+uTMDPDF_lowScale5(x,b,h)
  end do
  MeanTMD=MeanTMD/100d0
  
  end function
  
  
  !!! evaluetes standard deviation over replicas
  function StdTMD(x,b,h)
  real*8,dimension(-5:5)::StdTMD,F,MM
  real*8::x,b
  integer::i,h
  
  MM=MeanTMD(x,b,h)
  StdTMD=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  do i=1,100
   call uTMDPDF_SetLambdaNP(i)
   F=uTMDPDF_lowScale5(x,b,h)
   StdTMD=StdTMD+(F-MM)**2
  end do
  StdTMD=sqrt(StdTMD/100d0)
  
  end function
  
end program example