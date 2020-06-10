program example
use aTMDe_control
use TMDs
implicit none

integer::i
real*8::tmd1(-5:5)

call artemide_Initialize('const-test_Sivers',prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/Sivers20/')

call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))
call artemide_SetNPparameters_uTMDFF((/0.264d0,0.479d0,0.459d0,0.539d0/))
call artemide_SetNPparameters_SiversTMDPDF((/1d0,0.3d0,2d0,0.5d0,0.5d0/))

do i=1,10
  tmd1=SiversTMDPDF_5(0.1d0,i/10d0,2d0,2d0,1)
  write(*,'("{",F6.3,",",F12.8,"},")') i/10d0,tmd1(1)
  !write(*,*) qT(i),X(i)
end do

end program example
