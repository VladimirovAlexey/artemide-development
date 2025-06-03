program xSec_SIDIS_test

  ! Invocation of required modules
  use aTMDe_control
  use aTMDe_Numerics

  use uTMDPDF
  use uTMDFF

  use TMDX_SIDIS

  use TMDF
  use TMDF_KPC


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                       Declaration of variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Explicit declaration of variables
  implicit none

  ! Defines the number of points to compute
  integer,parameter :: NUM = 80

  ! Loop index
  integer :: i

  real*8, dimension(1:NUM) :: s,pT,z,x,Q,qT
  integer, dimension(1:NUM,1:4) :: proc
  logical, dimension(1:NUM) :: iC
  real*8, dimension(1:NUM,1:4) :: cuts

  real*8, dimension(1:NUM) :: xx, xxLP, xxKPC, FUUT1, FUUT2, FUUL, FUUTLP, FUULLP, FUU0
  real*8, dimension(1:NUM,1:2) :: MM

!   real*8 :: IKPC
  real*8, dimension(1:NUM) :: IKPC

  real*8 :: time1,time2

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  Set-up artemide and NP-parameters
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Initialization of arTeMiDe using the constants-file xSec_SIDIS.atmde
  !call artemide_Initialize('xSec_SIDIS.atmde',prefix='Prog/MySIDISTest/INI/')
  call artemide_Initialize('xSec_SIDIS.atmde',prefix='Prog/KPCtest/INI/')

  ! Setting the non-perturbative parameters (NP) for TMD evolution
  call artemide_SetNPparameters_TMDR((/1.5004d0, 0.05614d0, 0.03862d0, 0.0d0/))

  ! Setting the non-perturbative parameters (NP) for the uTMDPDF
  call artemide_SetNPparameters_uTMDPDF(&
  (/0.565d0, 0.0539d0, 0.5697d0, 6.64d0, 0.565d0, 20.07d0, 0.5697d0, 0.537d0, 1.07d0, 2.39d0, 0.0d0, 0.0d0/))

  ! Setting the non-perturbative parameters (NP) for the uTMDFF
  call artemide_SetNPparameters_uTMDFF((/0.69769d0, 0.712969d0, -0.133895d0, -0.841651d0, 0.846846d0,&
  0.774759d0, 1.5565d0, 1.1863d0, 0.692877d0, -0.569062d0, 0.0d0, 0.0d0/))

!    call artemide_SetNPparameters_uTMDFF((/0.69769d0, 0.d0, 0.d0, 0.d0, 0.d0,&
!    0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.0d0, 0.0d0/))

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!            Program to compute the SIDIS cross-section
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! Loop to choose the values of the variables involved in the xSec computation
  do i = 1, NUM

    ! We want to plot the xSec as a function of Q
!     Q(i) = 1.d0+(i-1)*1.d0
    Q(i) = 20.d0

    s(i) = 1500.d0 ! TeV!

    x(i) = 0.1d0
    z(i) = 0.5d0

    pT(i) = 0.01d0+(i-1)*0.05d0
    !pT(i) = 0.1d0

    ! proc = (p1,p2,p3,p4)
    ! p1 -> selects prefactor2
    ! p2,p3 -> select the hadrons type
    ! p4 -> selects kinematic kernel KK (KERNELpairs_SIDIS)
    proc(i,1:4)=(/1,1,1,2001/)

    iC(i)=.false.
    cuts(i,1:4)=(/20.d0,20.d0,-2.1d0,2.1d0/)

    MM(i,1:2)=(/1d0,1d0/)

  end do

  ! Loop to choose the values of the variables involved in the xSec computation
!     do i = 1, NUM
!
!       ! We want to plot the xSec as a function of qT
!       pT(i) = 0.0001d0+(i-1)*1
!
!       Q(i) = 2.d0
!
!       s(i) = 1500.d0 ! TeV!
!
!       x(i) = 0.1d0
!       z(i) = 0.3d0
!
! !       qT(i) = pT(i)/z(i)
!
!       ! proc = (p1,p2,p3,p4)
!       ! p1 -> selects prefactor2
!       ! p2,p3 -> select the hadrons type
!       ! p4 -> selects kinematic kernel KK (KERNELpairs_SIDIS)
!       proc(i,1:4)=(/1,1,1,2001/)
!
!       iC(i)=.false.
!       cuts(i,1:4)=(/20.d0,20.d0,-2.1d0,2.1d0/)
!
!       MM(i,1:2)=(/1d0,1d0/)
!
!     end do


  ! Loop to choose the values of the variables involved in the xSec computation
!     do i = 1, NUM
!
!       ! We want to plot the xSec as a function of qT
!       pT(i) = 0.0001d0
!
!       Q(i) = 2.d0
!
!       s(i) = 1500.d0 ! TeV!
!
!       x(i) = 0.1d0
!       z(i) = 0.3d0
!
!       !qT(i) = pT(i)/z(i)
!
!       ! proc = (p1,p2,p3,p4)
!       ! p1 -> selects prefactor2
!       ! p2,p3 -> select the hadrons type
!       ! p4 -> selects kinematic kernel KK (KERNELpairs_SIDIS)
!       proc(i,1:4)=(/1,1,1,2001/)
!
!       iC(i)=.false.
!       cuts(i,1:4)=(/20.d0,20.d0,-2.1d0,2.1d0/)
!
!       MM(i,1:2)=(/1d0,1d0/)
!
!     end do


!   IKPC = KPC_SIDISconv(100._dp,0._dp,0.5_dp,0.5_dp,10._dp,(/1,1,2001/))!*100._dp**2*3.14_dp/2/0.5_dp
!
!   write(*,*) "Conv_integral", IKPC

!   call xSec_SIDIS_BINLESS_List_forharpy(xx,(/1,1,1,2001/),1500._dp,0._dp,0.5_dp,0.5_dp,100._dp)


!-------------------------------------------------------------------------------------------
! This part is meant to check the convolution integral (employing some generic gaussians)
!-------------------------------------------------------------------------------------------
!   open(unit = 13, file = "conv_integral_with_gaussian.txt", status = "replace", action = "write")
!     do i = 1, NUM
!
!       IKPC(i) = KPC_SIDISconv(Q(i)**2,pT(i),x(i),z(i),Q(i),(/1,1,999/))
!
!     write(13, *) Q(i), IKPC(i)
!
!     end do
!   close(13)


!------------------------------------------------
! This part computes the SIDIS cross-section
!------------------------------------------------


! OPEN(UNIT=44, FILE=trim("integrand.dat"), ACTION="write", STATUS="replace")
!
! OPEN(UNIT=45, FILE=trim("jacobian.dat"), ACTION="write", STATUS="replace")
!
!
! call cpu_time(time1)
!
!   call xSec_SIDIS_BINLESS_List_forharpy(xx,proc,s,pT,z,x,Q,MM)
!
! call cpu_time(time2)
! !
! write(*,*) " Result---->>:"
!     do i = 1, NUM
!
!       write(*,*) pT(i)/z(i), xx(i)
! !       write(*,*) Q(i), xx(i)
!
!     end do
! !
! !
!   write(*,*) " COMPUTATION TIME:", time2-time1
!
! CLOSE(44, STATUS='KEEP')
!
! CLOSE(45, STATUS='KEEP')


!------------------------------------------------------------------------
!                    qT = 0.1
!------------------------------------------------------------------------
  ! Write data to file (KPCs case)
!   open(unit = 13, file = "xSec_KPCs.txt", status = "replace", action = "write")
!     do i = 1, NUM
!
!       write(13, *) Q(i), xx(i)
!
!     end do
!     close(13)
!
!   write(*,*) "Data written to ", "xSec_KPCs.txt"
!   write(*,*) "Use an external plotting program to visualize the data."
!   write(*,*) " COMPUTATION TIME:", time2-time1


  ! Write data to file (LP case)
!   open(unit = 13, file = "xSec_LP.txt", status = "replace", action = "write")
!     do i = 1, NUM
!       write(13, *) Q(i), xx(i)
!     end do
!   close(13)
!
!   write(*,*) "Data written to ", "xSec_LP.txt"
!   write(*,*) "Use an external plotting program to visualize the data."
!   write(*,*) " COMPUTATION TIME:", time2-time1


!------------------------------------------------------------------------
!                    Q = 2
!------------------------------------------------------------------------
  ! Write data to file (KPCs case)
!   open(unit = 13, file = "xSec_KPCs.txt", status = "replace", action = "write")
!     do i = 1, NUM
!
!       write(13, *) pT(i)/z(i), xx(i)
!
!     end do
!     close(13)
!
!   write(*,*) "Data written to ", "xSec_KPCs.txt"
!   write(*,*) "Use an external plotting program to visualize the data."
!   write(*,*) " COMPUTATION TIME:", time2-time1


    ! Write data to file (LP case)
!   open(unit = 13, file = "xSec_LP.txt", status = "replace", action = "write")
!     do i = 1, NUM
!       write(13, *) qT(i), xx(i)
!     end do
!   close(13)
!
!   write(*,*) "Data written to ", "xSec_LP.txt"
!   write(*,*) "Use an external plotting program to visualize the data."
!   write(*,*) " COMPUTATION TIME:", time2-time1


!------------------------------------------------------------------------
!                    Structure functions ratios
!------------------------------------------------------------------------
do i=1,NUM

    ! To select the FUUT

    ! With the x/2z Q2 factor
    FUUT1(i)=x(i)/2/z(i)*Q(i)**2*KPC_SIDISconv(Q(i)**2,pT(i),x(i),z(i),Q(i),(/1,1,2001/))

    ! Without the x/2z Q2 factor
!     FUUT2(i)=KPC_SIDISconv(Q(i)**2,pT(i),x(i),z(i),Q(i),(/1,1,2001/))

    ! In the LP case the FUUT is simply the LP convolution
    FUUTLP(i)=TMDF_F(Q(i)**2,pT(i),x(i),z(i),Q(i),Q(i)**2,Q(i)**2,(/1,1,2001/))

    FUU0(i)=TMDF_F(Q(i)**2,0.01d0,x(i),z(i),Q(i),Q(i)**2,Q(i)**2,(/1,1,2001/))

    ! With the x/2z Q2 factor
    !write(*,'("{",F12.6,",",F12.6,",",F12.6,",",F12.6,"},")') Q(i), FUUT1(i),FUUTLP(i),FUUT1(i)/FUUTLP(i)/x(i)**2/3.14d0
    write(*,'("{",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,"},")') pT(i), FUUT1(i),FUUTLP(i),FUU0(i),&
      (FUUT1(i)/x(i)**2/3.14d0-FUUTLP(i))/FUU0(i)

    ! Without the x/2z Q2 factor
!     write(*,'("{",F12.6,",",F12.6,",",F12.6,"},")') FUUT2(i),FUUTLP(i),FUUT2(i)/FUUTLP(i)

end do



end program xSec_SIDIS_test

