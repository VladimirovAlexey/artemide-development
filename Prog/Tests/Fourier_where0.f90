!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot for AS-term
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! module KTtest
! INCLUDE '/data/arTeMiDe_Repository/artemide/src/Code/KTspace/Fourier_Levin.f90'
! end module KTtest


! module KTtest2
! INCLUDE '/data/arTeMiDe_Repository/artemide/src/Code/KTspace/grid_inKT(new).f90'
! end module KTtest2

program example
use aTMDe_control
use uTMDPDF
use Fourier_Levin_uTMDPDF

implicit none

integer::i,j
real*8,dimension(-5:5)::FinQ
real*8,dimension(1:5,0:16,-5:5)::FinQ2
real*8::kT,x,Q,rr,rr2
integer::f

!!!! this code searches for the zero of the Fourier transform

call artemide_Initialize('ART23_MSHT_N4LL.atmde',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_uTMDPDF((/&
0.874245d0,0.913883d0,0.991563d0,6.05412d0,&
0.353908d0,46.6064d0,0.115161d0,1.53235d0,&
1.31966d0,0.434833d0, 0.0d0, 0.0d0/))


x=0.01d0
f=2
Q=3.2d0
!
! do i=1,200
!     Q=1.d0+i*0.5d0
!     rr=root()
!     write(*,'("{",F6.2,",",F8.2,"},")',advance="no") Q,rr
! end do
! write(*,*) " "

do i=1,40

x=10**(-i/10.d0)

Q=5.d0
rr=root()
Q=150.d0
rr2=root()
write(*,'("{",F8.4,",",F8.2,"},")',advance="no") x,(rr2-rr)/(150-5)
end do

contains


!!!! simple division by half
function root()
real*8::root
real*8::up,down,mid,upF,downF,midF
real*8::FF(-5:5)
integer::i

down=0.1d0
FF=uTMDPDF_inKT(x,down,Q,Q**2,1)
downF=FF(f)

up=200d0
FF=uTMDPDF_inKT(x,up,Q,Q**2,1)
upF=FF(f)

if(downF>0 .and. upF<0) then
    do i=1,20
        mid=(up+down)/2
        FF=uTMDPDF_inKT(x,mid,Q,Q**2,1)
        midF=FF(f)

        if(midF<0) then
            up=mid
        else
            down=mid
        end if
    end do
    root=(up+down)/2
else
    root=-1
end if


end function root

end program example
