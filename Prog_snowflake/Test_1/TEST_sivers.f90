!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Simple test-code for Snowlake. Checks compilation integrity and simple evaluation        !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SnowflakeTEST
use Snowflake
use Snowflake_Model
use SiversTMDPDF_model
implicit none

real*8::mu0,mu1,x,Q0,ss(-5:5)
integer::i,ff(0:20)

!call  SnowFlake_Initialize("Snowflake.ini")
call  SnowFlake_Initialize("Prog_snowflake/snowflake_forRep.ini")

call SnowFlake_Model_Initialize()

call SetNPparameters(&
    (/5.4000000000000004d0,  1.2000000000000002d0,       -1.1999999999999997d0,       0.78999999999999992d0,&
    4.5000000000000000d0,      -0.12700000000000000d0,        2.7999999999999998d0,      -0.34000000000000002d0,&
    -0.94000000000000006d0,       -5.5190000000000001d0,       -11.000000000000000d0,      -0.80000000000000016d0,&
    -1.0000000000000000d0,        3.0000000000000000d0,       -1.5000000000000000d0,       -1.0000000000000000d0,&
    -1.2000000000000000d0,        0.0000000000000000d0/))

call ModelInitialization(2)

mu0=1.d0
mu1=105.d0
call ComputeEvolution(mu0,mu1,alpha,U1=SplusU,D1=SplusD,S1=SplusS,U2=SminusU,D2=SminusD,S2=SminusS,G1=Tp,G2=Tm,&
    inputQ="C",inputG="T")

x=0.18935864271087760d0
Q0=24.732691d0

write(*,*) "b ",GetPDF(x,0.d0,Q0,5,outputT="T"),GetPDF(-x,0.d0,Q0,5,outputT="T")
write(*,*) "c ",GetPDF(x,0.d0,Q0,4,outputT="T"),GetPDF(-x,0.d0,Q0,4,outputT="T")
write(*,*) "s ",GetPDF(x,0.d0,Q0,3,outputT="T"),GetPDF(-x,0.d0,Q0,3,outputT="T")
write(*,*) "u ",GetPDF(x,0.d0,Q0,2,outputT="T"),GetPDF(-x,0.d0,Q0,2,outputT="T")
write(*,*) "s ",GetPDF(x,0.d0,Q0,1,outputT="T"),GetPDF(-x,0.d0,Q0,1,outputT="T")

ss=FNP(x,0.065007299988276146d0,1,(/0.5d0,0.d0/))

write(*,*) "---", ss


end program SnowflakeTEST
