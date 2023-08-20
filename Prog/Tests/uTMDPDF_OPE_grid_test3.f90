!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDPDF_OPE
use uTMDPDF
implicit none

real*8,dimension(0:50,0:50,-5:5)::TT1
integer::i,j
real::t1,t2,t3,t4
real*8::x,b

call uTMDPDF_OPE_Initialize('uTMDPDF.atmde',prefix='Prog/Tests/const-files/')

do i=0,50
do j=0,50
    b=i*0.1d0
    x=10**(-j*0.1d0)

    TT1(i,j,-5:5)=uTMDPDF_OPE_convolution(x,b,1)
end do
end do


do i=0,50
do j=0,50
    b=i*0.1d0
    x=10**(-j*0.1d0)

    write(*,'("{",F10.6,",",F16.10,",",F16.10,"},")') x,b,TT1(i,j,2)
end do
end do


end program example
