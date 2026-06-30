!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       arTeMiDe 3.04
!
!   Contains definitions interfaces common for the rest of artemide
!   Used in (almost) each artemide module.
!
!                               A.Vladimirov (13.10.2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aTMDe_interfaces
use aTMDe_numerics
implicit none

public

!!! this is interface for function of 1 variables used in integrations
abstract interface
    function func_1D(x)
        import::dp
        real(dp) :: func_1D
        real(dp), intent(in) ::x
    end function func_1D
end interface

!!! this is interface for function of 2 variables used in integrations
abstract interface
    function func_2D(x,y)
        import::dp
        real(dp) :: func_2D
        real(dp), intent(in) ::x,y
    end function func_2D
end interface

!!! this is interface for function (-5:5) in the integration
abstract interface
    function func_1D_array5(x)
        import::dp
        real(dp),dimension(-5:5) :: func_1D_array5
        real(dp), intent(in) ::x
    end function func_1D_array5
end interface

!!! this is interface for optTMD-like function (-5:5)
abstract interface
    function optTMD_proc(x,b,h)
        import::dp
        real(dp),dimension(-5:5) :: optTMD_proc
        real(dp), intent(in) ::x,b
        integer,intent(in)::h
    end function optTMD_proc
end interface

!!! this is interface for TMDatQ-like function (-5:5)
!!! it can be in b or kT space
abstract interface
    function TMDatQ_proc(x,kT,Q,h)
        import::dp
        real(dp),dimension(-5:5) :: TMDatQ_proc
        real(dp), intent(in) ::x,Q,kT
        integer,intent(in)::h
    end function TMDatQ_proc
end interface

!!! This is a very specific interface, it is used to exchange grids in-between Levin-transform and kt-grids
!!! Basically, it returns values of TMD(-5:5) as a matrix of numKsubgrids, each of kGridSize
abstract interface
    function TMDgrid_inKT(x,Q,h,arraySize1,arraySize2)
        import::dp
        integer,intent(in)::arraySize1,arraySize2
        real(dp),dimension(1:arraySize1,0:arraySize2,-5:5) :: TMDgrid_inKT
        real(dp), intent(in) ::x,Q
        integer,intent(in)::h

    end function TMDgrid_inKT
end interface

!!! this is interface for structure function
!!! with process0 being last 3 numbers of the process-numeration
abstract interface
    function strFUNC_proc(Q2,qT,x1,x2,process0)
        import::dp
        real(dp)::strFUNC_proc
        real(dp),intent(in)::Q2,qT,x1,x2
        integer,dimension(1:3),intent(in)::process0
    end function strFUNC_proc
end interface

end module aTMDe_interfaces
