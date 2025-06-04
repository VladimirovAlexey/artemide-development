!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Simple test-code for Snowlake. Checks compilation integrity and simple evaluation        !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SnowflakeTEST
use HexGrid

implicit none

real*8::x1,x2,x3,v1,v2
integer::n,k,i0,i1,i2,j0,j1,j2
logical:: inter


!public:subgridNUM_from_totalNUM_PHI,subgridNUM_from_totalNUM_RHO,totalNUM_from_subgridNUM_RHO,subgridNUM_from_totalNUM_PHI

call  Initialize_HexGrid("Snowflake.ini")

do n=0,815
    call get_X123_from_1Dindex(n,x1,x2,x3)
end do

! OPEN(UNIT=44, FILE=trim("prprpr.dat"), ACTION="write", STATUS="replace")
!
! do n=0,815
!     call get_X123_from_1Dindex(n,x1,x2,x3)
!     call RP_fromX12(x1,x2,v1,v2)
!     i0=getSubgrid_PHI(v2)
!     write(44,*) n,x1,x2,x3,i0
! end do
! CLOSE(44, STATUS='KEEP')

end program SnowflakeTEST
