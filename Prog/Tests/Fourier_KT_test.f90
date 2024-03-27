!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot for AS-term
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module KTtest
INCLUDE '/data/arTeMiDe_Repository/artemide/src/Code/KTspace/Fourier_overGrid.f90'
end module KTtest


program example
use KTtest
implicit none

integer::i,j


!call PrepareTransformationMatrix(8,(/0.00001d0,0.01d0,0.2d0,2.d0, 8.d0, 25.d0/))
call PrepareTransformationMatrix(16,(/0.00001d0,0.01d0,0.2d0,2.d0, 8.d0, 25.d0/))

end program example
