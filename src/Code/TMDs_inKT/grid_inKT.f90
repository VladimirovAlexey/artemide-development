!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which generates the grid for TMDs in kT-space
!     and restores the value.
!
!	v.3.00 Created (AV, 06.02.2024)
!
!				A.Vladimirov (06.02.2024)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! sets global variables for the X-griding
subroutine Grid_Initialize()
    DeltaX=(acosh(1._dp/xMin)**(1._dp/parX))/Nx
    DeltaQ=Log(Qmax/QMin)/NQ
    DeltaK=Log(parK*Qmax/KMin)/NK
end subroutine Grid_Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The grid in x
!!!! x_i=1/cosh^p(Delta(i-Nx))
!!!! x_0=xMin
!!!! x_Nx=1
!!!!


!!!! Value of X at the node i
pure function XatNode(i)
    real(dp):: XatNode
    integer,intent(in)::i

    XatNode=1._dp/cosh((DeltaX*(Nx-i))**parX)
end function XatNode

!!!! Value of X at the node i(real)
pure function XatNode_real(i)
    real(dp):: XatNode_real
    real(dp),intent(in)::i

    XatNode_real=1._dp/cosh((DeltaX*(Nx-i))**parX)
end function XatNode_real

!!!! Inverse X
pure function invX(x)
    real(dp):: invX
    real(dp), intent(in):: x

    invX=Nx-(acosh(1._dp/x)**(1._dp/parX))/DeltaX
end function invX

!!!! Value of low-grid value for given X
pure function NodeForX(x)
    integer:: NodeForX
    real(dp), intent(in):: x

    NodeForX=INT(invX(x))
end function NodeForX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The grid in Q
!!!! Q_i=Q_min Exp[i DeltaQ]
!!!! Q_0=Q_min
!!!! Q_NQ=Q_max
!!!!


!!!! Value of Q at the node i
pure function QatNode(i)
    real(dp):: QatNode
    integer,intent(in)::i

    QatNode=Qmin*Exp(i*DeltaQ)
end function QatNode

!!!! Value of Q at the node i(real)
pure function QatNode_real(i)
    real(dp):: QatNode_real
    real(dp),intent(in)::i

    QatNode_real=Qmin*Exp(i*DeltaQ)
end function QatNode_real

!!!! Inverse Q
pure function invQ(Q)
    real(dp):: invQ
    real(dp), intent(in):: Q

    invQ=log(Q/Qmin)/DeltaX
end function invQ

!!!! Value of low-grid value for given X
pure function NodeForQ(Q)
    integer:: NodeForQ
    real(dp), intent(in):: Q

    NodeForQ=INT(invQ(Q))
end function NodeForQ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The grid in kT
!!!! k_i=k_min Exp[i DeltaK]
!!!! Q_0=k_min
!!!! k_NK=parK*Q_max
!!!! Meanwhile for each value of Q the kT is limited by k<parK*Q


!!!! Value of kT at the node i
pure function KTatNode(i)
    real(dp):: KTatNode
    integer,intent(in)::i

    KTatNode=Kmin*Exp(i*DeltaK)
end function KTatNode

!!!! Value of kT at the node i(real)
pure function KTatNode_real(i)
    real(dp):: KTatNode_real
    real(dp),intent(in)::i

    KTatNode_real=Kmin*Exp(i*DeltaK)
end function KTatNode_real

!!!! Inverse kT
pure function invKT(KT)
    real(dp):: invKT
    real(dp), intent(in):: KT

    invKT=log(KT/Kmin)/DeltaK
end function invKT

!!!! Value of low-grid value for given X
pure function NodeForKT(KT)
    integer:: NodeForKT
    real(dp), intent(in):: KT

    NodeForKT=INT(invKT(KT))
end function NodeForKT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Make grid function
!!!!
!!!! num is the number of call TMD in the internal
!!!! M is the variable to which the grid is stored
!!!! numH is the number of hadrons to compute
subroutine MakeGrid(M,num,numH)
real(dp),dimension(0:NX,0:NK,0:NQ,-5:5,1:numH)::M
integer,intent(in)::num
integer,intent(in)::numH

integer::xi,ki,Qi,hi
real(dp)::x,k,Q
real*8::time1,time2

!$ real*8::omp_get_wtime

call cpu_time(time1)
!$ time1=omp_get_wtime()

if(outputlevel>2) write(*,*) 'arTeMiDe.',moduleName,' starts to compute grid.'

do hi=1,numH
!$OMP DO
do xi=0,NX
do Qi=0,NQ
do ki=0,NK
    x=XatNode(xi)
    k=KTatNode(ki)
    Q=QatNode(Qi)
    M(xi,ki,Qi,-5:5,hi)=x*(k**2)*Fourier(x,k,Q,Q**2,num,hi)
end do
end do
end do
!$OMP END DO

if(outputLevel>1 .and. numH>1) write(*,'(" ",A,": Grid for hadron ",I3," is done")') moduleName,hi
end do

call cpu_time(time2)
  !$ time2=omp_get_wtime()

if(outputlevel>1) then
if(numH>1) then
    write(*,'(" ",A,": Grids are built  (",I5," x",I5," x",I5,")  calc.time=",F6.2,"s. ")')&
    moduleName, Nx,NK,NQ, time2-time1
else
    write(*,'(" ",A,": Grid is built  (",I5," x",I5," x",I5,")  calc.time=",F6.2,"s. ")')&
    moduleName, NX,NK,NQ, time2-time1
end if
end if

end subroutine MakeGrid
