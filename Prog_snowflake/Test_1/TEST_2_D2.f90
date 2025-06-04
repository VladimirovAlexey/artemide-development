!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Simple test-code for Snowlake. Checks compilation integrity and simple evaluation        !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SnowflakeTEST
use Snowflake
!use EvolutionKernels
!use HexGrid

implicit none

real*8::mu0,mu1,Q,ff(0:24)
integer::i,j
!$ real*8::omp_get_wtime
real*8::t1,t2

!call  SnowFlake_Initialize("Snowflake.ini")
call  SnowFlake_Initialize("prog/TEST.ini")

mu0=1.d0
mu1=105.d0
call ComputeEvolution(mu0,mu1,alpha,U1=Tu,D1=Td,S1=Ts,U2=dTu,D2=dTd,S2=dTs,G1=Tpp,G2=Tmm,inputQ="C",inputG="T")
!call ComputeEvolution(mu0,mu1,alpha,U1=Tu,U2=dTu,inputQ="T",inputG="T")
!
! write(*,*) Tu(0.2d0,0.1d0),Tu(0.2d0,-0.1d0),Tu(0.3d0,0.1d0)
! write(*,*) Td(0.2d0,0.1d0),Td(0.2d0,-0.1d0),Td(0.3d0,0.1d0)
! write(*,*) Ts(0.2d0,0.1d0),Ts(0.2d0,-0.1d0),Ts(0.3d0,0.1d0)
! write(*,*) dTu(0.2d0,0.1d0),dTu(0.2d0,-0.1d0),dTu(0.3d0,0.1d0)
! write(*,*) dTd(0.2d0,0.1d0),dTd(0.2d0,-0.1d0),dTd(0.3d0,0.1d0)
! write(*,*) dTs(0.2d0,0.1d0),dTs(0.2d0,-0.1d0),dTs(0.3d0,0.1d0)
! write(*,*) dTpp(0.2d0,0.1d0),dTpp(0.2d0,-0.1d0),dTpp(0.3d0,0.1d0)
! write(*,*) dTmm(0.2d0,0.1d0),dTmm(0.2d0,-0.1d0),dTmm(0.3d0,0.1d0)

!call cpu_time(t1)
!$ t1=omp_get_wtime()

!OPEN(UNIT=44, FILE=trim("prprpr.dat"), ACTION="write", STATUS="replace")

! do i=0,51
! !
! !     Q=1.d0+i
! !
! !     ff(j)=D2(Q,11)
! !
! !     !write(44,*) x,ff
! !
! ! end do
!
! call cpu_time(t2)
! !$ t2=omp_get_wtime()
! write(*,*) "Time for computation ",t2-t1
!
! write(*,*) "  "
! CLOSE(44, STATUS='KEEP')

do i=0,200
write(*,'("{",F5.2,",",F12.6,"},")') 1.d0+i*0.2,D2(1.d0+i*0.2,11)
end do


contains
        !!!!! functions from the examplei n the paper
        !!!! some symmetric function
    function H(x,y)
        real*8,intent(in)::x,y
        real*8::H

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            H=(1-x**2)*(1-y**2)*(1-(x+y)**2)
        else
            H=0.d0
        end if
    end function H

    !!!! some symmetric function
    function Tu(x,y)
        real*8,intent(in)::x,y
        real*8::Tu
        real*8,parameter::pi=3.141592653589793d0

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            Tu=H(x,y)*cos(4*y)
        else
            Tu=0.d0
        end if
    end function Tu

        !!!! some symmetric function
    function Td(x,y)
        real*8,intent(in)::x,y
        real*8::Td
        real*8,parameter::pi=3.141592653589793d0

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            Td=H(x,y)*(2-cos(3*pi*H(x,y)))
        else
            Td=0.d0
        end if
    end function Td

            !!!! some symmetric function
    function Ts(x,y)
        real*8,intent(in)::x,y
        real*8::Ts

        Ts=-0.3*Td(x,y)
    end function Ts

    !!!! some asymmetric function
    function dTu(x,y)
        real*8,intent(in)::x,y
        real*8::dTu
        real*8,parameter::pi=3.141592653589793d0

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            dTu=H(x,y)*(sin(pi*y)+4*(x**2-(x+y)**2))
        else
            dTu=0.d0
        end if
    end function dTu

        !!!! some asymmetric function
    function dTd(x,y)
        real*8,intent(in)::x,y
        real*8::dTd
        real*8,parameter::pi=3.141592653589793d0

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            dTd=2*sin(pi*y)*(1-cos(H(x,y)))/max(abs(x),abs(y),abs(x+y))
        else
            dTd=0.d0
        end if
    end function dTd

                !!!! some symmetric function
    function dTs(x,y)
        real*8,intent(in)::x,y
        real*8::dTs

        dTs=-0.3*dTd(x,y)
    end function dTs

    !!!! some symmetric function
    function Tpp(x,y)
        real*8,intent(in)::x,y
        real*8::Tpp

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            Tpp=H(x,y)*max(abs(x),abs(y),abs(x+y))*sin(x+x+y)
        else
            Tpp=0.d0
        end if
    end function Tpp

    !!!! some symmetric function
    function Tmm(x,y)
        real*8,intent(in)::x,y
        real*8::Tmm

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            Tmm=H(x,y)*max(abs(x),abs(y),abs(x+y))*cos(x+x+y)
        else
            Tmm=0.d0
        end if
    end function Tmm

    !!!!!! Function for alpha_s of QCD.
    !!!!!! Here is a simple model for alpha-s. You can use any other model, or interface it with LHAPDF, or other code
    !!!!!! Preserve the interface.
    pure function alpha(mu)
        real*8,intent(in)::mu
        real*8::alpha
        alpha=12.566370614359172d0/11/(2*log(mu)+3.d0)
    end function alpha
end program SnowflakeTEST
