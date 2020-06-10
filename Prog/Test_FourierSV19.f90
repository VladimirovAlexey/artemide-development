!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program computes various combinations of structure functions and compare to 
! the values that were get in SV19 fit, and other parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDF
implicit none

integer,parameter::Nt=9
real*8::testString(1:Nt),compareString(1:Nt)
integer::j

call artemide_Initialize('const-TEST',prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Prog/')

!----------------------------------------------------------------------------------------------------------------------------------

call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

write(*,*) '------Q=91 ----- SV19 -----'
call MakeTest(91d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/0.7604178436,     0.2981829881,     0.0869299546,     0.0325154417,     0.0145826214,     0.0072388044,   &
      0.0037788444,     0.0019953303,     0.0010182871/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
  !write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '


write(*,*) '------Q=4 ----- SV19 -----'
call MakeTest(4d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/7.5464385042,     7.3309290002,     6.9526442086,     6.4495777376,     5.8676802015,     5.2490897406,   &
	  4.6296371889,     4.0336928910,     3.4778755852/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!    write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '

write(*,*) '------Q=350 ----- SV19 -----'
call MakeTest(350d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/0.2026195631,     0.0177544877,     0.0040486336,     0.0014530781,     0.0006365438,     0.0003079589,   &
	0.0001553596,     0.0000776231,     0.0000354622/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!    write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '

!----------------------------------------------------------------------------------------------------------------------------------

call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.2d0, 0.2d0, 0d0, 0d0, 0d0, 0.d0, 0.d0/))

write(*,*) '------Q=91 ----- GAUSS 0.2-----'
call MakeTest(91d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/0.9261662664,     0.2974237992,     0.0812610171,     0.0310107579,     0.0141128611,     0.0070584074,   &
      0.0036981681,     0.0019553611,     0.0009971691/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!    write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '


write(*,*) '------Q=4 ----- GAUSS 0.2 -----'
call MakeTest(4d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/10.7167060564,    10.4517486915,     9.9730498036,     9.3116273048,     8.5080435539,     7.6085257676,  &
      6.6604720951,     5.7081474166,     4.7892784814/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!    write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '

write(*,*) '------Q=350 ----- GAUSS 0.2 -----'
call MakeTest(350d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/0.2303226150,     0.0175053940,     0.0040260902,     0.0014483722,     0.0006351366,     0.0003074235,   &
	0.0001551246,     0.0000775036,     0.0000353963/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!    write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '

!----------------------------------------------------------------------------------------------------------------------------------

call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.8d0, 0.8d0, 0d0, 0d0, 0d0, 0.d0, 0.d0/))

write(*,*) '------Q=91 ----- GAUSS 0.8-----'
call MakeTest(91d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/0.6177994285,     0.3068256075,     0.0897330648,     0.0329971612,     0.0147426663,     0.0073055321,   &
      0.0038101246,     0.0020111742,     0.0010267962/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!     write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '


write(*,*) '------Q=4 ----- GAUSS 0.8 -----'
call MakeTest(4d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/3.3610282230,     3.3335465091,     3.2826028861,     3.2092142664,     3.1148156553,     3.0012697608,    &
      2.8707079521,     2.7256001230,     2.5685663897/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!     write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '

write(*,*) '------Q=350 ----- GAUSS 0.8 -----'
call MakeTest(350d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/0.1813385912,     0.0178508197,     0.0040578836,     0.0014550327,     0.0006371416,     0.0003081835,    &
	0.0001554624,     0.0000776721,     0.0000354905/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!     write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '

!----------------------------------------------------------------------------------------------------------------------------------

call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.02d0, 0.02d0, 0d0, 0d0, 0d0, 0.2d0, 0.d0/))

write(*,*) '------Q=91 ----- EXP 0.2-----'
call MakeTest(91d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/1.1938345880,     0.2865425244,     0.0790191888,     0.0304697191,     0.0139354025,     0.0069873815,  &
	0.0036656367,     0.0019390283,     0.0009884717/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!    write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '


write(*,*) '------Q=4 ----- EXP 0.2 -----'
call MakeTest(4d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/46.6582223029,    37.8242988429,    27.0803186592,    18.5254441310,    12.7395146851,     8.9154769233,   &
      6.3187532726,     4.5015261346,     3.1976519474/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!    write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '

write(*,*) '------Q=350 ----- EXP 0.2 -----'
call MakeTest(350d0**2,0.1d0,0.1d0,testString)

 compareString=1d0*(/0.2602117279,     0.0174054372,     0.0040167392,     0.0014463972,     0.0006345353,     0.0003071899,    &
      0.0001550253,     1.0000,     1.0000/)
do j=1,Nt
  write(*,"(F14.6,'   ')",advance='no') testString(j)/compareString(j)
!    write(*,"(F14.10,',   ')",advance='no') testString(j) 
end do
write(*,*) ' '

contains

subroutine MakeTest(Q2,x1,x2,res)
  real*8::Q2,x1,x2
  real*8,intent(out)::res(1:Nt)
  real*8::qT(1:Nt)
  integer::i
  
  do i=1,Nt
    qT(i)=0.1d0+0.35d0*sqrt(Q2)*(i-1)/Nt
    res(i)=TMDF_F(Q2,qT(i),x1,x2,sqrt(Q2),Q2,Q2,2)
  end do
  
end subroutine MakeTest
  

end program example
  