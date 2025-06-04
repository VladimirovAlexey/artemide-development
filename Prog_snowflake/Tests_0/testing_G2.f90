!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Simple test-code for Snowlake. Checks compilation integrity and simple evaluation        !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SnowflakeTEST
use HexGrid
use EvolutionKernels
use SnowFlake

implicit none

integer::i
real*8::t1,t2,mu0,mu1,t,r1,r2,r3

call  SnowFlake_Initialize("TEST.ini","prog/")

mu0=1.d0

call ComputeEvolution(mu0,50.d0,alpha,U1=initialF,U2=initialA,G1=initialG,inputQ="T",inputG="T")

call cpu_time(t1)
!call SaveMatrixG2("kernels/20x25/")

!call ReadMatrixG2("kernels/20x25/")
do i=1,50
    t=i*0.1d0
    mu1=exp(t/2)

    r1=G2(0.1d0,mu1)
    r2=G2(0.3d0,mu1)
    r3=G2(0.6d0,mu1)

    write(*,'("{", F5.2, ", ", F8.4, ", ", F8.4, ", ", F8.4, "},")') mu1,r1,r2,r3
end do

call cpu_time(t2)

write(*,*) ">>>>>>>> Timing: ", t2-t1

contains

    !!!! some symmetric function
    function initialF(x,y)
        real*8::x,y,initialF

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialF=(1-x**2)*(1-y**2)*(1-(x+y)**2)
        else
            initialF=0.d0
        end if
    end function initialF

    !!!! some asymmetric function
    function initialA(x,y)
        real*8::x,y,initialA
        real*8,parameter::pi=3.141592653589793d0

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialA=sin(pi*y)*cos(pi/2*x)*cos(pi/2*(x+y))
        else
            initialA=0.d0
        end if
    end function initialA

    !!!! some symmetric function
    function initialG(x,y)
        real*8::x,y,initialG

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialG=sin(2*x+y)*(1-x**2)*(1-y**2)*(1-(x+y)**2)/sqrt(max(abs(x),abs(y),abs(x+y)))
        else
            initialG=0.d0
        end if
    end function initialG

    !!!!!! Function for alpha_s of QCD.
    !!!!!! Here is a simple model for alpha-s. You can use any other model, or interface it with LHAPDF, or other code
    !!!!!! Preserve the interface.
    pure function alpha(mu)
        real*8,intent(in)::mu
        real*8::alpha
        alpha=12.566370614359172d0/11/(2*log(mu)+3.d0)
    end function alpha

end program SnowflakeTEST
