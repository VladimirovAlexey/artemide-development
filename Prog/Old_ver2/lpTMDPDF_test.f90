program example
use lpTMDPDF
implicit none
 
 real*8::xx(-5:5)
 
 call lpTMDPDF_Initialize('const-DYfit18_NNLO+LP',prefix='/misc/data2/braun/vla18041/arTeMiDe_Repository/')
 
 call lpTMDPDF_SetLambdaNP((/0.228d0,0306d0/))
 
 xx=lpTMDPDF_lowScale50(0.1d0,0.1d0,1)
 write(*,*) xx(0)
 xx=lpTMDPDF_lowScale50(0.01d0,0.1d0,1)
 write(*,*) xx(0)
 xx=lpTMDPDF_lowScale50(0.001d0,0.1d0,1)
 write(*,*) xx(0)
 
end program example