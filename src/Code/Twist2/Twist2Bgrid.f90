!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD-OPE modules 
!       that operates at twist-2. It is inclucded (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!       The interpolation type is defined by K, 
!       as the segments runs out of grid, the interpolation order is reduced 
!	This part is devoted to the common functions for B-grid
!       BGrid_Initialize()
!       real BatNode(i)      - gives b_i
!       real invB(x)         - gives i such that b_i=b
!       integer NodeForB(x)  - gives i such that b_i<b<b_{i+1}
!
!!!!!---------------------------------
! ADD to the main module
! real(dp) :: BMAX=25._dp !!! maximum B
! real(dp) :: BMIN=1d-6 !!! minimum B
! integer :: NB=250 !!! number of points in grid
! real(dp) :: DeltaB !!! increment of grid
!
!	v.3.00 Created (AV, 18.07.2023)
!
!				A.Vladimirov (18.07.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! sets global variables for the B-griding
subroutine BGrid_Initialize()
    DeltaB=log(BMIN/BMAX)/Nb
end subroutine BGrid_Initialize

!!!! The algorithm is similar to QCD num.
!!!! I.e. it computes the convolution matrices

!!!! Value of B at the node i
pure function BatNode(i)
    real(dp):: BatNode
    integer,intent(in)::i
    
    BatNode=BMIN*exp(-i*DeltaB)
end function BatNode

!!!! Value of B at the node i(real)
pure function BatNode_real(i)
    real(dp):: BatNode_real
    real(dp),intent(in)::i

    BatNode_real=BMIN*exp(-i*DeltaB)
end function BatNode_real

!!!! Inverse B
pure function invB(b)
    real(dp):: invB
    real(dp), intent(in):: b
    
    invB=log(BMIN/b)/DeltaB
end function invB

!!!! Value of low-grid value for given B
pure function NodeForB(b)
    integer:: NodeForB
    real(dp), intent(in):: b
    
    NodeForB=INT(invB(b))
    
end function NodeForB

