!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use uTMDPDF_OPE
use uTMDPDF
implicit none

real*8,dimension(-5:5)::f1,f2

call uTMDPDF_Initialize('const-uTMDPDF_MSHT',prefix='Prog/Tests/const-files/')
call uTMDPDF_OPE_Initialize("file")

call MakeGrid()

f1=ExtractFromGrid(0.23d0,2.3d0,1)
f2=CxF_compute(0.23d0,2.3d0,1,.true.)

call TestGrid()

! call uTMDPDF_SetLambdaNP((/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/))
!
! call TestGrid_2()

end program example
