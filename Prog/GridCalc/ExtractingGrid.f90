program example
use aTMDe_control
use aTMDe_numerics
use uTMDPDF_OPE
use TMDX_DY
implicit none

integer:: n,i,f
real*8::x,b

real*8 :: time1,time2

call artemide_Initialize('ART_KPC.atmde',prefix='Prog/GridCalc/')
!   Setting the non-perturbative parameters (NP) for TMD evolution
call artemide_SetNPparameters_TMDR((/1.5004d0, 0.05614d0, 0.03862d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse

!   Setting the non-perturbative parameters (NP) for the uTMDPDF
call artemide_SetNPparameters_uTMDPDF(&
(/0.565d0, 0.0539d0, 0.5697d0, 6.64d0, 0.565d0, 20.07d0, 0.5697d0, 0.537d0, 1.07d0, 2.39d0, 0.0d0, 0.0d0/)) !sacado del ejemplo xSec del MiniCourse

n= uTMDPDF_OPE_GetGridSize()
write(*,*) n
do i=1,120
call uTMDPDF_OPE_GetGridPoint(i,x,b,f)
!write(*,*) i,x,b,f
write(*,'("{",I4,",",F12.6,",",F12.6,",",I3,"},")') i,x,b,f
end do
write(*,*) "---------------------------------------------------------"
end program example
