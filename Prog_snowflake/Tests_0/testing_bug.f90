!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Example code for computation of evolution with Snowflake library                         !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program example
use SnowFlake_Model
use SnowFlake
implicit none

real*8::x1,x2,x3
real*8::mu0,mu1,t
real*8::r1,r2,r3
integer::i

call  SnowFlake_Initialize("snowflake_forRep.ini","Prog_snowflake/Tests_0/")
!call  SnowFlake_Initialize("Snowflake.ini")

call SnowFlake_Model_Initialize()

call SetNPparameters((/3.112d0,   -0.157d0,   &
                           0.128d0,    1.570d0,    0.127d0,    0.794d0,  &
                           -0.180d0,   -1.594d0,   -5.019d0,   -5.872d0,  &
                           0.000d0,    0.000d0,    0.000d0,    0.000d0,    &
                           0.000d0,    0.000d0,    0.000d0,    0.000d0/))

mu0=1.d0

call ComputeEvolution(mu0,25.d0,alpha,U1=SplusU,D1=SplusD,S1=SplusS,U2=SminusU,D2=SminusD,S2=SminusS,G1=Tp,G2=Tm,&
    inputQ="C",inputG="T")

do i=1,50
    t=i*0.1d0
    mu1=exp(t/2)

    !call ComputeEvolution(mu0,mu1,alpha,U1=initialF,U2=initialA,G1=initialG,inputQ="T",inputG="T")

    !!! Here I call for result of Evolution for T(0.2,-0.3,0.1) for u quark
    x1=0.2d0
    x2=-0.d0
    r1=GetPDF(0.011d0,x2,mu1,2,outputT='T')
    !!! Here I call for result of Evolution for S^+(0.2,-0.3,0.1) for u quark
    r2=GetPDF(0.1d0,x2,mu1,2,outputT='T')
    !!! Here I call for result of Evolution for frak{S}^+(0.2,-0.3,0.1) for u quark
    r3=GetPDF(0.18935d0,0.d0,mu1,2,outputT='T')
    !!! The values are
    write(*,'("{", F5.2, ", ", F8.4, ", ", F8.4, ", ", F8.4, "},")') mu1,r1,r2,r3

end do

end program example
