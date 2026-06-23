program dump_grid_csv
!! Dumps full grid-vs-exact comparison to CSV for every (x,b) test point.
!!
!! Columns: x, b, ex[-5],gr[-5], ex[-4],gr[-4], ..., ex[5],gr[5]
!! (22 flavor columns = 11 flavors × 2; fl=0 included as zeros when withGluon=F)
!!
!! Identical setup to test_optGrid: same ini-file, NP params, test points.
use aTMDe_IO
use aTMDe_optGrid
use aTMDe_numerics
use aTMDe_control
use uTMDPDF_OPE
implicit none

character(*), parameter :: inifile = 'test.atmde'
character(*), parameter :: prefix  = 'Tests/aTMDe_optGrid/'

integer,  parameter :: Nx = 45, Nb = 100
real(dp), parameter :: xMIN = 0.00001_dp
real(dp), parameter :: bMAX = 15._dp

real(dp) :: xtest(Nx), btest(Nb)
real(dp), dimension(-5:5) :: from_grid, from_exact

integer :: fu, ix, ib, fv
type(optGrid) :: grid

!! ---- build test-point arrays (identical to test_optGrid) ----
do ix = 1, Nx
    xtest(ix) = 10._dp**(log10(xMIN) * (Nx+1-ix) / (Nx+1))
end do
do ib = 1, Nb
    btest(ib) = bMAX * ib / (Nb+1)
end do

!! ---- initialize ----
call artemide_Initialize(inifile, prefix=prefix)
call artemide_SetNPparameters_uTMDPDF([&
    0.874245_dp, 0.913883_dp, 0.991563_dp, 6.05412_dp, &
    0.353908_dp, 46.6064_dp,  0.115161_dp, 1.53235_dp, &
    1.31966_dp,  0.434833_dp, 0.0_dp,      0.0_dp])

!! ---- build grid ----
grid = optGrid(trim(prefix)//trim(inifile), '*4   ', '*E   ', 1, .false., 'dump_grid_csv', 1)
call grid%MakeGrid(F_test)

!! ---- dump CSV ----
open(newunit=fu, file='Tests/aTMDe_optGrid/grid_dump.csv', &
     status='replace', action='write')

write(fu,'(A)') 'x,b,'// &
    'ex[-5],gr[-5],ex[-4],gr[-4],ex[-3],gr[-3],ex[-2],gr[-2],ex[-1],gr[-1],'// &
    'ex[0],gr[0],'// &
    'ex[1],gr[1],ex[2],gr[2],ex[3],gr[3],ex[4],gr[4],ex[5],gr[5]'

do ix = 1, Nx
    do ib = 1, Nb
        from_exact = x_times_OPE(xtest(ix), btest(ib), 1)
        from_grid  = grid%Extract(xtest(ix), btest(ib), 1)
        write(fu, '(ES22.14,",",ES22.14,22(",",ES22.14))') &
            xtest(ix), btest(ib), &
            (from_exact(fv), from_grid(fv), fv=-5, 5)
    end do
end do

close(fu)
write(*, '(A)') 'Written: Tests/aTMDe_optGrid/grid_dump.csv'

contains

function F_test(x, b, h)
    real(dp), intent(in) :: x, b
    integer,  intent(in) :: h
    real(dp), dimension(-5:5) :: F_test
    F_test = x * uTMDPDF_OPE_convolution(x, b, h, .false.)
end function F_test

function x_times_OPE(x, b, h)
    real(dp), intent(in) :: x, b
    integer,  intent(in) :: h
    real(dp), dimension(-5:5) :: x_times_OPE
    x_times_OPE = x * uTMDPDF_OPE_convolution(x, b, h, .false.)
end function x_times_OPE

end program dump_grid_csv
