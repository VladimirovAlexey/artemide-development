program example
use uTMDPDF
use lpTMDPDF
implicit none

real*8,dimension(-5:5):: F1,F2,F3,F4,F5,F6,F7,F8
real*8::W0(0:80),W1(0:80),W2(0:80)
real*8,dimension(1:2):: A,B,C
real*8::x,bT
integer::i
 !A=(/2.5000d0,0.0376d0,0.2481d0,8.1556d0,276.5180d0,2.5262d0,-4.9586d0, 0.d0,0d0/)
 !B=(/2.5000d0,0.0376d0,0.3481d0,4.1556d0,176.5180d0,2.5262d0,-2.9586d0, 0.d0,0d0/)
 !C=(/2.5000d0,0.0376d0,0.1481d0,9.1556d0,276.5180d0,2.5262d0,-6.9586d0, 0.d0,0d0/)
  
  A=(/0.23,0.11/)
  B=(/0.23,0.31/)
  !call uTMDPDF_Initialize('const-DYfit18_NNLO+LP','/home/vla18041/LinkData2/arTeMiDe_Repository/')
  call lpTMDPDF_Initialize('const-DYfit18_NNLO+LP','/home/vla18041/LinkData2/arTeMiDe_Repository/')
 
!  x=0.001d0
!  bT=0.1d0
! 
! call lpTMDPDF_SetLambdaNP(A,.false.)
! 
! F1=lpTMDPDF_lowScale50(x,bT,1)
! write(*,*) F1
! 
! call lpTMDPDF_SetLambdaNP(B,.false.)
! 
! F2=lpTMDPDF_lowScale50(x,bT,1)
! 
! call lpTMDPDF_SetLambdaNP(A,.True.)
! 
! call lpTMDPDF_SetLambdaNP(A,.true.)
! 
! F3=lpTMDPDF_lowScale50(x,bT,1)
! 
! call lpTMDPDF_SetLambdaNP(B,.true.)
! 
! F4=lpTMDPDF_lowScale50(x,bT,1)
! 
! write(*,*) F1(0)
! write(*,*) F3(0)
! Write(*,*) (F1(0)-F3(0))/F1(0)
! Write(*,*) (F2(0)-F4(0))/F4(0)

call lpTMDPDF_SetLambdaNP((/0.228d0,0.306d0/),.false.)
!call uTMDPDF_SetLambdaNP((/0.228d0,0.306d0/),.false.,.true.)
x=0.1d0
do i=0,80
  bT=0.001d0+i*0.05d0
  F1=lpTMDPDF_lowScale50(x,bT,1)
  W0(i)=F1(0)
end do

call lpTMDPDF_SetScaleVariation(0.5d0)

do i=0,80
  bT=0.001d0+i*0.05d0
  F1=lpTMDPDF_lowScale50(x,bT,1)
  W1(i)=F1(0)
end do

call lpTMDPDF_SetScaleVariation(2.0d0)

do i=0,80
  bT=0.001d0+i*0.05d0
  F1=lpTMDPDF_lowScale50(x,bT,1)
  W2(i)=F1(0)
end do
 do i=0,80
 write(*,*) '{',0.001d0+i*0.05d0,',',W0(i),',',W1(i),',',W2(i),'},'
 end do
end program example