program example
use uTMDPDF
use lpTMDPDF
use TMDs_inKT
implicit none

real*8::W0(0:80),W1(0:80),W2(0:80),F1(-5:5),F2(-5:5)
real*8,dimension(1:2):: A,B,C
real*8::x,bT,qT,lu,llp,v1,v2,v3
real*8::r1,r2,r3
integer::i,j
integer,parameter::NUM=29
real*8::Lpath(0:NUM,1:2)
 !A=(/2.5000d0,0.0376d0,0.2481d0,8.1556d0,276.5180d0,2.5262d0,-4.9586d0, 0.d0,0d0/)
 !B=(/2.5000d0,0.0376d0,0.3481d0,4.1556d0,176.5180d0,2.5262d0,-2.9586d0, 0.d0,0d0/)
 !C=(/2.5000d0,0.0376d0,0.1481d0,9.1556d0,276.5180d0,2.5262d0,-6.9586d0, 0.d0,0d0/)
  
  A=(/0.23,0.11/)
  B=(/0.23,0.31/)
  call uTMDPDF_Initialize('const-DYfit18_NNLO+LP','/home/vla18041/LinkData2/arTeMiDe_Repository/')
  call lpTMDPDF_Initialize('const-DYfit18_NNLO+LP','/home/vla18041/LinkData2/arTeMiDe_Repository/')
  call TMDs_inKT_Initialize('const-DYfit18_NNLO+LP','/home/vla18041/LinkData2/arTeMiDe_Repository/')
 

call lpTMDPDF_SetLambdaNP((/0.228d0,0.306d0/),.false.)
call uTMDPDF_SetLambdaNP((/0.228d0,0.306d0/),.false.,.true.)


call uTMDPDF_SetLambdaNP((/0.5d0,0*0.306d0/),.false.,.true.)
call lpTMDPDF_SetLambdaNP((/0.6d0,0*0.306d0/),.false.)


! do i=1,50
!   x=10**(-i*2d0/25d0)!0.01d0
!   qT=05d0!0.001d0+i*4d0/50
!   F1=lpTMDPDF_kT_50(x,qT,1)
!   F2=uTMDPDF_kT_50(x,qT,1)
!   W0(i)=F1(0)
!   W1(i)=F2(0)
! end do
! 
! 
!  do i=1,50
!  !write(*,*) '{',0.001d0+i*4d0/50,',', Abs(W0(i))/(Abs(W1(i))+0.0001d0),'},'
!  write(*,*) '{',10**(-i*2d0/25d0),',', Abs(W0(i))/(Abs(W1(i))+0.0001d0),'},'
!  end do
  
  
  x=0.001d0
  qT=0.5d0
  
  do i=0,NUM
    
    lu=i*1d0/(NUM+1)+0.025
    
    v1=0.025d0
    v3=1d0-0.025d0
    
    r1= ratioAT(x,qT,lu,v1)
    r3= ratioAT(x,qT,lu,v3)
    do j=1,20
    if(r1>=1d0 .and. r3<1d0) then
      v2=(v1+v3)/2d0
      r2=ratioAT(x,qT,lu,v2)
      
      if(r2>=1d0) then
	v1=v2
	r1=r2
      else
	v3=v2
	r3=r2
      end if
    
    else if(r1<1d0 .and. r3>=1d0) then
      v2=(v1+v3)/2d0
      r2=ratioAT(x,qT,lu,v2)
      
      if(r2<1d0) then
	v1=v2
	r1=r2
      else
	v3=v2
	r3=r2
      end if
    else
      Lpath(i,1:2)=(/-1d0,-1d0/)
      exit
    end if
    end do
    lPath(i,1:2)=(/lu,(v1+v3)/2d0/)
    
    write(*,*) lPath(i,1:2),ratioAT(x,qT,lPath(i,1),lPath(i,2))
  end do

  !write(*,*) '{',lu,',',llp,',', Abs(W0(i))/(Abs(W1(i))+0.0001d0),'},'
  
  do i=0,NUM
    write(*,*) '{',lPath(i,1),',',lPath(i,2),'},'
  end do
  
contains
  !!! function evaluate the ratio
  function ratioAT(x,qT,l1,l2)
    real*8:: F11(-5:5),F22(-5:5)
    real*8::x,qT,l1,l2,ratioAT
    
    call uTMDPDF_SetLambdaNP((/l1,0.0306d0,0d0,0d0,0d0,0d0,0d0/),.false.,.true.)
    call lpTMDPDF_SetLambdaNP((/l2,0*0.306d0/),.false.)
    
    F11=lpTMDPDF_kT_50(x,qT,1)
    F22=uTMDPDF_kT_50(x,qT,1)
    ratioAT=Abs(F11(0))/(Abs(F22(0))+0.000001d0)
  
  end function ratioAT
  
end program example