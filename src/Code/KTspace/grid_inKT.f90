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

!!!!!!!!!!!!!!!!!!!!!!!!!!!! TO BE DEFINED IN THE MAIN CODE
!!! real(dp)::DeltaX_inKT,DeltaK_inKT, DeltaQ_inKT
!!! real(dp)::parX_inKT,parK_inKT
!!! real(dp)::xMin_inKT,KMin_inKT,Qmin_inKT,Qmax_inKT
!!! integer::NX_inKT,NK_inKT,NQ_inKT
!!! real(dp),allocatable::grid_inKT(:,:,:,:,:)

!!! sets global variables for the X-griding
subroutine Grid_Initialize()
    DeltaX_inKT=(acosh(1._dp/xMin_inKT)**(1._dp/parX_inKT))/Nx_inKT
    DeltaQ_inKT=Log(Qmax_inKT/QMin_inKT)/NQ_inKT
    DeltaK_inKT=Log(parK_inKT*Qmax_inKT/KMin_inKT)/NK_inKT
end subroutine Grid_Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The grid in x
!!!! x_i=1/cosh^p(Delta(i-Nx))
!!!! x_0=xMin
!!!! x_Nx=1
!!!!

!!!! Value of X at the node i
pure function XatNode_inKT(i)
    real(dp):: XatNode_inKT
    integer,intent(in)::i

    XatNode_inKT=1._dp/cosh((DeltaX_inKT*(Nx_inKT-i))**parX_inKT)
end function XatNode_inKT

!!!! Value of X at the node i(real)
pure function XatNode_inKT_real(i)
    real(dp):: XatNode_inKT_real
    real(dp),intent(in)::i

    XatNode_inKT_real=1._dp/cosh((DeltaX_inKT*(Nx_inKT-i))**parX_inKT)
end function XatNode_inKT_real

!!!! Inverse X
pure function invX_inKT(x)
    real(dp):: invX_inKT
    real(dp), intent(in):: x

    invX_inKT=Nx_inKT-(acosh(1._dp/x)**(1._dp/parX_inKT))/DeltaX_inKT
end function invX_inKT

!!!! Value of low-grid value for given X
pure function NodeForX_inKT(x)
    integer:: NodeForX_inKT
    real(dp), intent(in):: x

    NodeForX_inKT=INT(invX_inKT(x))
end function NodeForX_inKT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The grid in Q
!!!! Q_i=Q_min Exp[i DeltaQ]
!!!! Q_0=Q_min
!!!! Q_NQ=Q_max
!!!!

!!!! Value of Q at the node i
pure function QatNode_inKT(i)
    real(dp):: QatNode_inKT
    integer,intent(in)::i

    QatNode_inKT=Qmin_inKT*Exp(i*DeltaQ_inKT)
end function QatNode_inKT

!!!! Value of Q at the node i(real)
pure function QatNode_inKT_real(i)
    real(dp):: QatNode_inKT_real
    real(dp),intent(in)::i

    QatNode_inKT_real=Qmin_inKT*Exp(i*DeltaQ_inKT)
end function QatNode_inKT_real

!!!! Inverse Q
pure function invQ_inKT(Q)
    real(dp):: invQ_inKT
    real(dp), intent(in):: Q

    invQ_inKT=log(Q/Qmin_inKT)/DeltaQ_inKT
end function invQ_inKT

!!!! Value of low-grid value for given Q
pure function NodeForQ_inKT(Q)
    integer:: NodeForQ_inKT
    real(dp), intent(in):: Q

    NodeForQ_inKT=INT(invQ_inKT(Q))
end function NodeForQ_inKT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The grid in kT
!!!! k_i=k_min Exp[i DeltaK]
!!!! Q_0=k_min
!!!! k_NK=parK*Q_max
!!!! Meanwhile for each value of Q the kT is limited by k<parK*Q

!!!! Value of kT at the node i
pure function KTatNode_inKT(i)
    real(dp):: KTatNode_inKT
    integer,intent(in)::i

    KTatNode_inKT=Kmin_inKT*Exp(i*DeltaK_inKT)
end function KTatNode_inKT

!!!! Value of kT at the node i(real)
pure function KTatNode_inKT_real(i)
    real(dp):: KTatNode_inKT_real
    real(dp),intent(in)::i

    KTatNode_inKT_real=Kmin_inKT*Exp(i*DeltaK_inKT)
end function KTatNode_inKT_real

!!!! Inverse kT
pure function invKT_inKT(KT)
    real(dp):: invKT_inKT
    real(dp), intent(in):: KT

    invKT_inKT=log(KT/Kmin_inKT)/DeltaK_inKT
end function invKT_inKT

!!!! Value of low-grid value for given X
pure function NodeForKT_inKT(KT)
    integer:: NodeForKT_inKT
    real(dp), intent(in):: KT

    NodeForKT_inKT=INT(invKT_inKT(KT))
end function NodeForKT_inKT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Make-grid function in KT space
!!!!
!!!! The stored values, are multiplied by (x k^2)-factor
subroutine ComputeGrid_inKT()

integer::xi,ki,Qi,hi
real(dp)::x,k,Q,Q_add
real*8::time1,time2

!$ real*8::omp_get_wtime

call cpu_time(time1)
!$ time1=omp_get_wtime()

if(outputlevel>1) write(*,*) 'arTeMiDe.',moduleName,' starts to compute grid in KT space.'

do hi=1,numOfHadrons
!$OMP DO
do xi=0,NX_inKT
do Qi=0,NQ_inKT
do ki=0,NK_inKT
    x=XatNode_inKT(xi)
    k=KTatNode_inKT(ki)
    Q=QatNode_inKT(Qi)
    Q_add=QatNode_inKT(Qi+1)!!!! used for the check

    if(k<parK_inKT*Q_add) then!!!! store only kT<p*Q
        grid_inKT(xi,ki,Qi,-5:5,hi)=x*(k**2)*Fourier_ev(x,k,Q,Q**2,hi)
    else
        grid_inKT(xi,ki,Qi,-5:5,hi)=0._dp
    end if
end do
end do
end do
!$OMP END DO

if(outputLevel>1 .and. numOfHadrons>1) write(*,'(" ",A,": Grid for hadron ",I3," is done")') moduleName,hi
end do

call cpu_time(time2)
  !$ time2=omp_get_wtime()

if(outputlevel>1) then
if(numOfHadrons>1) then
    write(*,'(" ",A,": Grids are built  (",I5," x",I5," x",I5,")  calc.time=",F6.2,"s. ")')&
    moduleName, Nx_inKT,NK_inKT,NQ_inKT, time2-time1
else
    write(*,'(" ",A,": Grid is built  (",I5," x",I5," x",I5,")  calc.time=",F6.2,"s. ")')&
    moduleName, Nx_inKT,NK_inKT,NQ_inKT, time2-time1
end if
end if

!!!!! TEST GRID!!!!
  if(makeTest_inKT) then
    call TestGrid_inKT()
  end if
end subroutine ComputeGrid_inKT

function ExtractFromGrid_inKT(x,kT,Q,h)
real(dp),intent(in)::x,kT,Q
integer,intent(in)::h
real(dp),dimension(-5:5)::ExtractFromGrid_inKT

real(dp)::indexX,indexK,indexQ,fX,fK,fQ
integer::iX,iK,IQ
real(dp)::var1,var2,var3,var4 !!dummyvariables
real(dp)::interI(0:3,0:3,-5:5),interI2(0:3,-5:5),inter3(-5:5)
integer::dummyI

!!! checking exeptions
  if(h==0 .or. h>numOfHadrons) then
    write(*,*) ErrorString('the hadron '//numToStr(h)//' is not found in the KT-grid',moduleName)
    write(*,*) 'arTeMiDe: evaluation STOP'
    stop
  end if

  if(x<xMin_inKT) then
   write(*,*) ErrorString('The TMD with x ='//numToStr(x)//'is called. Current KT-grid size is up to '//&
   numToStr(xMin_inKT)//'. Enlarge boundaries.',moduleName)
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if
  if(x>1d0) then
   write(*,*) ErrorString('The TMD with x >1 ('//numToStr(x)//') is called from KT-grid.',moduleName)
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if
  if(kT<0d0) then
   write(*,*) ErrorString('The TMD with kT <0 ('//numToStr(kT)//') is called from KT-grid.',moduleName)
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if
  !!! If kT is larger than pQ
  if(kT>parK_inKT*Q) then
    write(*,*) ErrorString('The TMD with kT > pQ ('//numToStr(parK_inKT*Q)//') is called from KT-grid.',moduleName)
    write(*,*) 'arTeMiDe: evaluation STOP'
    stop
  end if
  if(QMin_inKT<0d0) then
   write(*,*) ErrorString('The TMD with Q <Qmin ('//numToStr(Qmin_inKT)//') is called from KT-grid.',moduleName)
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if
  if(QMax_inKT<0d0) then
   write(*,*) ErrorString('The TMD with Q >QMax ('//numToStr(Qmax_inKT)//') is called from KT-grid.',moduleName)
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if

    !!!!!!!!!Finding X index
    indexX=invX_inKT(x)
    iX=INT(indexX)
    if(iX>0) then
        if(iX<NX_inKT-1) then !!normal restoration X
            iX=INT(indexX)
        else !!X in the last interval
            iX=INT(indexX)-2
        end if
    else !! X in the first interval
        iX=INT(indexX)+1
    end if
    fX=indexX-iX !!! fraction part (automatically shifted by +- 1 if needed)

    !!!!!!!!!Finding Q index
    indexQ=invQ_inKT(Q)
    iQ=INT(indexQ)
    if(iQ>0) then
        if(iQ<NQ_inKT-1) then !!normal restoration Q
            iQ=INT(indexQ)
        else !!Q in the last interval
            iQ=INT(indexQ)-2
        end if
    else !! Q in the first interval
        iQ=INT(indexQ)+1
    end if
    fQ=indexQ-iQ !!! fraction part (automatically shifted by +- 1 if needed)



    !!! 1) Interpolate in KT
    if(kT<KMin_inKT) then
    !!!! no need to interpolate.
        interI(0,0,-5:5)=grid_inKT(iX-1,0,iQ-1,-5:5,h)
        interI(1,0,-5:5)=grid_inKT(iX+0,0,iQ-1,-5:5,h)
        interI(2,0,-5:5)=grid_inKT(iX+1,0,iQ-1,-5:5,h)
        interI(3,0,-5:5)=grid_inKT(iX+2,0,iQ-1,-5:5,h)
        interI(0,1,-5:5)=grid_inKT(iX-1,0,iQ+0,-5:5,h)
        interI(1,1,-5:5)=grid_inKT(iX+0,0,iQ+0,-5:5,h)
        interI(2,1,-5:5)=grid_inKT(iX+1,0,iQ+0,-5:5,h)
        interI(3,1,-5:5)=grid_inKT(iX+2,0,iQ+0,-5:5,h)
        interI(0,2,-5:5)=grid_inKT(iX-1,0,iQ+1,-5:5,h)
        interI(1,2,-5:5)=grid_inKT(iX+0,0,iQ+1,-5:5,h)
        interI(2,2,-5:5)=grid_inKT(iX+1,0,iQ+1,-5:5,h)
        interI(3,2,-5:5)=grid_inKT(iX+2,0,iQ+1,-5:5,h)
        interI(0,3,-5:5)=grid_inKT(iX-1,0,iQ+2,-5:5,h)
        interI(1,3,-5:5)=grid_inKT(iX+0,0,iQ+2,-5:5,h)
        interI(2,3,-5:5)=grid_inKT(iX+1,0,iQ+2,-5:5,h)
        interI(3,3,-5:5)=grid_inKT(iX+2,0,iQ+2,-5:5,h)
    else
        indexK=invKT_inKT(kT)
        iK=INT(indexK)
        if(iK>0) then
            if(iK<NK_inKT-1) then !!normal restoration K
                iK=INT(indexK)
            else !!K in the last interval
                iK=INT(indexK)-2
            end if
        else !! K in the first inteval
            iK=INT(indexK)+1
        end if
        fK=indexK-iK !!! fraction part (automatically shifted by +- 1 if needed)

        ! 1) intepolation procedure over B (at grids of X)
        var1=-fK*(fK-1d0)*(fK-2d0)
        var2=3d0*(fK+1d0)*(fK-1d0)*(fK-2d0)
        var3=-3d0*(fK+1d0)*fK*(fK-2d0)
        var4=(fK+1d0)*fK*(fK-1d0)

        do dummyI=0,4
            interI(0,dummyI,-5:5)= &
                (var1*grid_inKT(iX-1,iK-1,iQ-1+dummyI,-5:5,h)&
                +var2*grid_inKT(iX-1,iK  ,iQ-1+dummyI,-5:5,h)&
                +var3*grid_inKT(iX-1,iK+1,iQ-1+dummyI,-5:5,h)&
                +var4*grid_inKT(iX-1,iK+2,iQ-1+dummyI,-5:5,h))/6d0

            interI(1,dummyI,-5:5)= &
                (var1*grid_inKT(iX,iK-1,iQ-1+dummyI,-5:5,h)&
                +var2*grid_inKT(iX,iK  ,iQ-1+dummyI,-5:5,h)&
                +var3*grid_inKT(iX,iK+1,iQ-1+dummyI,-5:5,h)&
                +var4*grid_inKT(iX,iK+2,iQ-1+dummyI,-5:5,h))/6d0

            interI(2,dummyI,-5:5)=&
                (var1*grid_inKT(iX+1,iK-1,iQ-1+dummyI,-5:5,h)&
                +var2*grid_inKT(iX+1,iK  ,iQ-1+dummyI,-5:5,h)&
                +var3*grid_inKT(iX+1,iK+1,iQ-1+dummyI,-5:5,h)&
                +var4*grid_inKT(iX+1,iK+2,iQ-1+dummyI,-5:5,h))/6d0
            if(iX+2==Nx_inKT) then
                interI(3,dummyI,-5:5)=(/0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp/)
            else
                interI(3,dummyI,-5:5)=&
                (var1*grid_inKT(iX+2,iK-1,iQ-1+dummyI,-5:5,h)&
                +var2*grid_inKT(iX+2,iK  ,iQ-1+dummyI,-5:5,h)&
                +var3*grid_inKT(iX+2,iK+1,iQ-1+dummyI,-5:5,h)&
                +var4*grid_inKT(iX+2,iK+2,iQ-1+dummyI,-5:5,h))/6d0
            end if
        end do
    end if

    !! 2) intepolation procedure over Q (from values computed at KT)
    var1=-fQ*(fQ-1d0)*(fQ-2d0)
    var2=3d0*(fQ+1d0)*(fQ-1d0)*(fQ-2d0)
    var3=-3d0*(fQ+1d0)*fQ*(fQ-2d0)
    var4=(fQ+1d0)*fQ*(fQ-1d0)

    interI2(0,-5:5)=(var1*interI(0,0,-5:5)+var2*interI(0,1,-5:5)+var3*interI(0,2,-5:5)+var4*interI(0,3,-5:5))/6
    interI2(1,-5:5)=(var1*interI(1,0,-5:5)+var2*interI(1,1,-5:5)+var3*interI(1,2,-5:5)+var4*interI(1,3,-5:5))/6
    interI2(2,-5:5)=(var1*interI(2,0,-5:5)+var2*interI(2,1,-5:5)+var3*interI(2,2,-5:5)+var4*interI(2,3,-5:5))/6
    interI2(3,-5:5)=(var1*interI(3,0,-5:5)+var2*interI(3,1,-5:5)+var3*interI(3,2,-5:5)+var4*interI(3,3,-5:5))/6

    !! 3) intepolation procedure over X (from values computed at KT+Q)
    var1=-fX*(fX-1d0)*(fX-2d0)
    var2=3d0*(fX+1d0)*(fX-1d0)*(fX-2d0)
    var3=-3d0*(fX+1d0)*fX*(fX-2d0)
    var4=(fX+1d0)*fX*(fX-1d0)

    inter3(-5:5)=(var1*interI2(0,-5:5)+var2*interI2(1,-5:5)+var3*interI2(2,-5:5)+var4*interI2(3,-5:5))/6

    !!! finally we devide by factor x k^2
    ExtractFromGrid_inKT=inter3/x/kT**2

    !write(*,*) ">>>>>>",indexX,indexK,indexQ

do dummyI=-5,5
    if(ISNAN(ExtractFromGrid_inKT(dummyI))) then

        write(*,*) ErrorString('KT-grid extraction produced NaN. EVALUATION STOP',moduleName)
        write(*,*) '----- information on last call -----'
        write(*,*) 'x=',x,'kT=',kT,'Q=',Q,' flavor=',dummyI
        write(*,*) 'fx=',fx,'fk=',fk,'fQ=',fQ
        write(*,*) 'ix=',ix,'ik=',ik,'iQ=',iQ
        write(*,*) 'extraction=',ExtractFromGrid_inKT(dummyI)
        write(*,*) 'interI2=',interI2(0:3,dummyI)

        stop
    end if
end do

end function ExtractFromGrid_inKT

!!!! This subroutine test grid by computing the values of function at middle points and compare to exact values
subroutine TestGrid_inKT()
    integer::h,i,j,k
    real(dp),dimension(-5:5):: test_local

    character*5,dimension(0:5),parameter::xSTR=(/"    0","10^-4"," 0.01","  0.1","  0.8","    1"/)
    character*5,dimension(0:4),parameter::kSTR=(/"    0","  0.1","   1.","   5."," rest"/)
    character*5,dimension(0:5),parameter::QSTR=(/"    1","    2","    5","   15","   25","  100"/)
    integer,dimension(0:5)::xVal
    integer,dimension(0:4)::kVal
    integer,dimension(0:5)::QVal

    xVal=(/0,NodeForX_inKT(0.0001d0),NodeForX_inKT(0.01d0),NodeForX_inKT(0.1d0),NodeForX_inKT(0.8d0),NX_inKT-1/)
    kVal=(/0,NodeForKT_inKT(0.1d0),NodeForKT_inKT(1.d0),NodeForKT_inKT(5.d0),NK_inKT/)
    QVal=(/NodeForQ_inKT(1.d0),NodeForQ_inKT(2.d0),NodeForQ_inKT(5.d0),&
        NodeForQ_inKT(15.d0),NodeForQ_inKT(25.d0),NodeForQ_inKT(100.d0)/)

   write(*,*) color("  ---  TEST OF GRIDS ---" , c_yellow_bold)
   write(*,*) "* Avarage values of [interpolation/exact-1] in the middle points of grid (max deviation)"
   do h=1,numOfHadrons
    write(*,*) color("  -----  Hadron "//trim(int4ToStr(h)) , c_yellow_bold)

    write(*,&
    '(" |              |     sBar    |     uBar    |     dBar    |   gluon     |      d      |      u      |      s      |")')
    write(*,*)&
       "|----------------------------------------------------------------------------------------------------------------|"
      do i=1,size(xVal)-1
      if(xMin_inKT<xVal(i)) then
        write(*,*)"|",xSTR(i-1),"<x<",xSTR(i),"                                                      "
        do k=1,size(QVal)-1
        if(QMin_inKT<QVal(k) .and. QVal(k)<QMax_inKT) then
        write(*,*)"|",QSTR(i-1),"<Q<",QSTR(i),"                                                      "
        write(*,*)&
    "|----------------------------------------------------------------------------------------------------------------|"
            do j=1,size(kVal)-1
      test_local=TestPartOfGrid_inKT(xVal(i-1),xVal(i),kVal(j-1),kVal(j),QVal(k-1),QVal(k),h)
      write(*,'(" |",A5,"<kT<",A5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |")') &
      kSTR(j-1),kSTR(j),test_local(-3:3)
            end do
    write(*,*)&
    "|----------------------------------------------------------------------------------------------------------------|"
        end if
        end do
        end if
    end do
     write(*,*) color("  -----  test for hadron "//trim(int4ToStr(h))//" complete." , c_yellow_bold)
  end do
end subroutine TestGrid_inKT

!!! actual test of the grid. Runs around the square, compute and perform avarage.
function TestPartOfGrid_inKT(Ix_Low,Ix_high,Ik_low,Ik_high,IQ_low,IQ_high,h)
integer,intent(in)::Ix_Low,Ix_high,Ik_low,Ik_high,IQ_low,IQ_high,h
real(dp),dimension(-5:5)::TestPartOfGrid_inKT
real(dp),dimension(1:(Ix_high-Ix_Low),1:(Ik_high-Ik_low),1:(IQ_high-IQ_low),-5:5)::f1,f2
integer::ix,ik,iQ,i
real(dp)::x_local,k_local,Q_local
real(dp),dimension(-5:5)::comulant
integer::numOfNodes

comulant=0._dp
numOfNodes=0

do ix=1,Ix_high-Ix_Low
    x_local=XatNode_inKT_real(Ix_Low+ix-0.5_dp)
    do ik=1,Ik_high-Ik_Low
        k_local=KTatNode_inKT_real(Ik_Low+ik-0.5_dp)
        do iQ=1,IQ_high-IQ_Low
            Q_local=QatNode_inKT_real(IQ_Low+iQ-0.5_dp)
            if(k_local<parK_inKT*Q_local) then

                f1(ix,ik,iQ,-5:5)=ExtractFromGrid_inKT(x_local,k_local,Q_local,h)
                f2(ix,ik,iQ,-5:5)=Fourier_ev(x_local,k_local,Q_local,Q_local**2,h)
                if(.not.includeGluon) then
                    f1(ix,ik,iQ,0)=1._dp
                    f2(ix,ik,iQ,0)=1._dp
                end if
                numOfNodes=numOfNodes+1
                comulant=comulant+abs(f1(ix,ik,iQ,-5:5)/f2(ix,ik,iQ,-5:5)-1._dp)

!                 do i=-3,3
!                 if(abs(f1(ix,ik,iQ,i)/f2(ix,ik,iQ,i)-1._dp)>1._dp) then
!                     write(*,*) "EXCEPTIONALLY BAD POINT -->",i,x_local,k_local,Q_local,f1(ix,ik,iQ,i),f2(ix,ik,iQ,i)
!                 end if
!                 end do
            end if
        end do
    end do
end do

if(numOfNodes>0) then
    TestPartOfGrid_inKT=comulant/numOfNodes
else
    write(*,*) ">>>>>>> =0"
    TestPartOfGrid_inKT=0._dp
end if

end function TestPartOfGrid_inKT
