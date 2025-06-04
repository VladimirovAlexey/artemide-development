!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Simple test-code for Snowlake. Checks compilation integrity and simple evaluation        !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use SnowFlake
use EvolutionKernels
use HexGrid
implicit none


real*8::t1,t2
!$ real*8::omp_get_wtime
real*8::x1,x2,x3,ff,ff1,ff2,ff3,mu0,mu1,Q
integer::i,j
real*8,dimension(0:509)::grid,grid2,grid3,grid4

call  SnowFlake_Initialize("TEST.ini","prog/")
!call  Initialize_HexGrid("prog/TEST.ini")

grid=GETgrid(Tu)
grid2=Dgrid_dX2(grid)
!grid3=GETgrid(initialA)
!grid4=Dgrid_dX2(grid3)
!
mu0=1.d0
mu1=25.d0
! ! !call ComputeEvolution(mu0,mu1,alpha,U1=initialF,U2=initialA,G1=initialG,inputQ="T",inputG="T")
call ComputeEvolution(mu0,mu1,alpha,U1=initialF,U2=initialA,D1=initialF2,D2=initialA2,G1=initialG,inputQ="T",inputG="T")

call cpu_time(t1)
!$ t1=omp_get_wtime()

OPEN(UNIT=44, FILE=trim("gr.dat"), ACTION="write", STATUS="replace")
!
! do i=0,size(grid2)-1
! call get_X123_from_1Dindex(i,x1,x2,x3)
!
! !if(abs(x2)<0.00000001d0 .and. x1>0) write(*,'("{",F12.8,",",F16.8,",",F16.8,"},")') x1, grid(i),grid2(i)
!
! write(44,*) x1,x2,grid(i),grid2(i),grid3(i),grid4(i)
!
! end do
!
! stop


do i=-50,50
do j=-50,50

if(i==0 .and. j==0) cycle
x1=i/51.d0
x2=j/51.d0
if(abs(x1+x2)>1.) cycle

ff=GetPDF(x1,x2,10.d0,1,outputT='T')
ff1=GetPDF(x1,x2,10.d0,2,outputT='T')
ff2=GetPDF(x1,x2,25.d0,1,outputT='T')
ff3=GetPDF(x1,x2,25.d0,2,outputT='T')

write(44,*) x1,x2,ff,ff1,ff2,ff3
end do
end do

call cpu_time(t2)
!$ t2=omp_get_wtime()
write(*,*) "Time for computation of evolution",t2-t1

write(*,*) "  "
CLOSE(44, STATUS='KEEP')

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
        real*8::r
        r=max(abs(x),abs(y),abs(x+y))

        if(abs(x)<=1 .and. abs(y)<=1 .and. abs(x+y)<=1) then
            !Tu=Exp(2*(2*x+y)**2)!H(x,y)*cos(4*y)*0+H(x,y)*(2-cos(3*pi*H(x,y)))
            Tu=H(x,y)*(2-cos(3*pi*H(x,y)))
            !Tu=sin(r)+2+r**2
        else
            Tu=0.d0
        end if
    end function Tu

    !!!! some symmetric function
    function initialF(x,y)
        real*8,intent(in)::x,y
        real*8::initialF

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialF=(1-x**2)*(1-y**2)*(1-(x+y)**2)*(1+(x**2+(x+y)**2)/(1+y**2))
        else
            initialF=0.d0
        end if
    end function initialF

    !!!! some symmetric function
    function initialF2(x,y)
        real*8,intent(in)::x,y
        real*8::initialF2

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialF2=(1-x**2)*(1-y**2)*(1-(x+y)**2)*(1+0.2*x*(x+y)+-1.2*x*y**2*(x+y))
        else
            initialF2=0.d0
        end if
    end function initialF2

    !!!! some asymmetric function
    function initialA(x,y)
        real*8,intent(in)::x,y
        real*8::initialA
        real*8,parameter::pi=3.141592653589793d0

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialA=sin(pi*y)*cos(pi/2*x)*cos(pi/2*(x+y))
        else
            initialA=0.d0
        end if
    end function initialA

        !!!! some asymmetric function
    function initialA2(x,y)
        real*8::x,y,initialA2
        real*8,parameter::pi=3.141592653589793d0

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialA2=sin(pi*y)*cos(pi/2*x)*cos(pi/2*(x+y))/(y**2+0.3d0)
        else
            initialA2=0.d0
        end if
    end function initialA2

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
end program example

