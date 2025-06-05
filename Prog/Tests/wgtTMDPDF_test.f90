!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
!use aTMDe_control
use wgtTMDPDF
use QCDinput
use snowflake
implicit none

integer::i,iMax
real*8::bMax,bStep,x
real*8,allocatable::b(:)
real*8::TT(-5:5)

!call artemide_Initialize('wgtTMDPDF.atmde',prefix='Prog/Tests/const-files/')
call wgtTMDPDF_Initialize('wgtTMDPDF_v3.atmde',prefix='Prog/Tests/const-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call wgtTMDPDF_SetLambdaNP((/0.0d0, 0.d0/))


call  SnowFlake_Initialize("Prog_snowflake/TEST_16x8.ini")

call ComputeEvolution(1.d0,105.d0,alpha,U1=Tu,D1=Td,S1=Ts,U2=dTu,D2=dTd,S2=dTs,G1=Tpp,G2=Tmm,inputQ="C",inputG="T")

x=0.1

bMax=5d0
bStep=0.1d0
iMax=int(bMax/bStep)

allocate(b(1:iMax))

do i=1,iMax
  b(i)=bStep*i
end do

bMax=0.5d0

! do i=1,iMax
! 
!     TT=wgtTMDPDF_lowScale5(x,b(i),1)
!     write(*,*) "{",b(i),",", TT(1),"},"
! end do

! do i=1,iMax
! 
!     TT=wgtTMDPDF_lowScale5(10**(-b(i)),0.2d0,1)
!     write(*,*) "{",-b(i),",", TT(1),"},"
! end do

TT=wgtTMDPDF_inB(0.1d0,0.2d0,1)
write(*,*) TT
stop
write(*,*) '---------------------- D -----------------------'

do i=1,100

     x=i*0.01d0
    TT=wgtTMDPDF_inB(x,bMax,1)
    !TT=x_hPDF(i*0.01d0,11.d0,1)
    write(*,'("{",F12.8,",",F14.10,"},")',advance='no') x, x*TT(1)
end do
write(*,*)

write(*,*) '---------------------- U -----------------------'

do i=1,100

     x=i*0.01d0
    TT=wgtTMDPDF_inB(x,bMax,1)
    !TT=x_hPDF(i*0.01d0,11.d0,1)
    write(*,'("{",F12.8,",",F14.10,"},")',advance='no') x, x*TT(2)
end do
write(*,*)

! write(*,*) '---------------------- S -----------------------'
!
! do i=1,100
!
!      x=i*0.01d0
!     TT=wgtTMDPDF_inB(x,bMax,1)
!     !TT=x_hPDF(i*0.01d0,11.d0,1)
!     write(*,'("{",F12.8,",",F14.10,"},")',advance='no') x, x*TT(3)
! end do
! write(*,*)


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
            dTd=2*sin(pi*y)*(1-cos(H(x,y)))/sqrt(1-H(x,y))
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

end program example
