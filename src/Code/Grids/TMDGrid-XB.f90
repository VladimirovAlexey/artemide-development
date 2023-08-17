!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD-evaluation modules
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!	Be AWARE of possible clash of variable names.
!
!   Significant part of this code migrated from the artemide 2.01
!	This part is devoted to the Grid evaluation
!
!				A.Vladimirov (16.08.2023)
!---------------------------------------------------------------------------------------
!   To the code of main module ADD
! ! !!!------------------------- SPECIAL VARIABLES FOR GRID (used by TMDGrid-XB)------------------
! ! real(dp), dimension(:,:,:,:), allocatable :: gridMain !!!! THIS IS HUGE(!) matrix for the grid
! ! real(dp), dimension(:,:,:,:), allocatable :: interpolationParameters !!!! for b>bGrid_Max we
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Griding functions  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!! prepare variables
!!!! hList is the list of hadrons to prepare grid
subroutine TMDGrid_XB_Initialize(hList)
  integer,dimension(:),intent(in)::hList
  if(useGrid) then
    allocate(gridMain(0:NX,0:NB,-5:5,1:numberOfHadrons))
    !!! these are paramerer bE and aE
    allocate(interpolationParameters(1:2,0:NX,-5:5,1:numberOfHadrons))
  end if

  !!! save the list of hadrons
  numberOfHadrons=size(hList)
  allocate(hadronsInGRID(1:numberOfHadrons))
  hadronsInGRID=hList

end subroutine TMDGrid_XB_Initialize

!!! this subroutine create grid.
!!! it stores the value of CxF (which is OPE multiplied by x)
subroutine MakeGrid()
  real(dp):: x_local,b_local
  integer:: iX,iB,j,h,h_local
  
  real(dp),dimension(0:Nx-1,-5:5)::f1,f2,aE!! for interpolation computation
  real(dp)::b1,b2,time1,time2

  !$ real*8::omp_get_wtime
  
  call cpu_time(time1)
  !$ time1=omp_get_wtime()
  if(outputlevel>2) write(*,*) 'arTeMiDe.',moduleName,' starts to compute grid.'

  do h=1,numberOfHadrons
   h_local=hadronsInGRID(h)
   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(x_local,b_local)
   do iX=0,Nx
    do iB=0,Nb
     x_local=XatNode(iX)
     b_local=BatNode(iB)
     
     if(iX==Nx) then
      !!! Nx corresponds to x=1, all PDF MUST be 0.
      gridMain(iX,iB,-5:5,h)=0d0
     else
       gridMain(iX,iB,-5:5,h)=CxF_compute(x_local,b_local,h_local,withGluon)
     end if
    end do    
    end do
    !$OMP END PARALLEL DO
    if(outputLevel>1 .and. numberOfHadrons>1) write(*,'(" ",A,": Grid for hadron ",I3," is done")') moduleName,h_local
   end do

  !!!!!!!!!!!!!**************************************************************************************************************
  !!!!!!!!!!!!!*****************************TALE PART************************************************************************
  !!!!!!!!!!!!!**************************************************************************************************************

  !!! the large-b is interpolated by the power tale.
  !!! for it we compute 2 parameters aE=interpolationParameters(1), and bE=interpolationParameters(2)
  !!! the interpolation is done by f=(b/bE)**aE
  do h=1,numberOfHadrons
   h_local=hadronsInGRID(h)
   b1=BatNode(Nb-1)
   f1=gridMain(0:Nx-1,Nb-1,-5:5,h_local)
   b2=BatNode(Nb)
   f2=gridMain(0:Nx-1,Nb,-5:5,h_local)
   aE=log(f2/f1)/log(b2/b1)
   interpolationParameters(1,0:Nx-1,-5:5,h)=aE
   interpolationParameters(2,0:Nx-1,-5:5,h)=b1*(f1)**(-1_dp/aE)
  end do
  !!! for i=Nx, ae=1000, and b1=10d-6 (so the result is 0)
  interpolationParameters(1,Nx,-5:5,h)=1000_dp
  interpolationParameters(2,Nx,-5:5,h)=0.000000001_dp
  call cpu_time(time2)
  !$ time2=omp_get_wtime()

  if(outputlevel>1) then
    if(numberOfHadrons>1) then
      write(*,'(" ",A,": Grids are built  (",I5," x",I5,")  calc.time=",F6.2,"s. ")')&
        moduleName, Nx,NB, time2-time1
    else
      write(*,'(" ",A,": Grid is built  (",I5," x",I5,")  calc.time=",F6.2,"s. ")')&
        moduleName, NX,NB, time2-time1
    end if
  end if
    
end subroutine MakeGrid


!!!! this code is largerly migrated from artemide 2.01
function ExtractFromGrid(x,bT,hadron)
  real(dp),intent(in)::x,bT
  integer,intent(in)::hadron
  real(dp),dimension(-5:5)::ExtractFromGrid

  real(dp),dimension(0:3,-5:5):: interI
  real(dp)::indexX,indexB,fX,fB
  integer::i,iX,iB,h
  real(dp)::var1,var2,var3,var4 !!dummyvariables

  !!!searching for hadron
  h=0
  do i=1,numberOfHadrons
    if(hadronsInGRID(i)==hadron) then
      h=i
      exit
    end if
  end do


  if(h==0) then
    write(*,*) ErrorString('the hadron '//numToStr(hadron)//' is not found in the grid',moduleName)
    write(*,*) 'arTeMiDe: evaluation STOP'
    stop
  end if

  if(x<XMin) then
   write(*,*) ErrorString('The TMD with x ='//numToStr(x)//'is called. Current grid size is up to '//&
   numToStr(XMin)//'. Enlarge boundaries.',moduleName)
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if
  if(x>1d0) then
   write(*,*) ErrorString('The TMD with x >1 ('//numToStr(x)//') is called.',moduleName)
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if
  if(bT<0d0) then
   write(*,*) ErrorString('The TMD with bT <0 ('//numToStr(bT)//') is called.',moduleName)
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if

  !!!!!!!!!Finding X index
  indexX=invX(x)
  iX=INT(indexX)
    if(iX>0) then
     if(iX<NX-1) then !!normal restoration X
      iX=INT(indexX)
     else !!X in the last interval
      iX=INT(indexX)-2
      end if
    else !! X in the first interval
      iX=INT(indexX)+1
    end if
    fX=indexX-iX !!! fraction part (automatically shifted by +- 1 if needed)

    !!!! finding extrapolation over b, for nodes in X
   if(bT>bMax) then

      !! 1) extrapolate
      interI(0,-5:5)=(bT/interpolationParameters(2,iX-1,-5:5,h))**(interpolationParameters(1,iX-1,-5:5,h))
      interI(1,-5:5)=(bT/interpolationParameters(2,iX,-5:5,h))**(interpolationParameters(1,iX,-5:5,h))
      interI(2,-5:5)=(bT/interpolationParameters(2,iX+1,-5:5,h))**(interpolationParameters(1,iX+1,-5:5,h))
      if(iX+2==Nx) then
      interI(3,-5:5)=(/0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp/)
      else
      interI(3,-5:5)=(bT/interpolationParameters(2,iX+2,-5:5,h))**(interpolationParameters(1,iX+2,-5:5,h))
      end if
   else if(bT<bMIN) then
     !! 1) freezed value (at IB=0)
      interI(0,-5:5)=gridMain(iX-1,0,-5:5,h)
      interI(1,-5:5)=gridMain(iX,0,-5:5,h)
      interI(2,-5:5)=gridMain(iX+1,0,-5:5,h)
      if(iX+2==Nx) then
      interI(3,-5:5)=(/0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp/)
      else
      interI(3,-5:5)=gridMain(iX+2,0,-5:5,h)
      end if
   else!!! b inside the main region

      indexB=invB(bT)
      iB=INT(indexB)
      if(iB>0) then
        if(iB<NB-1) then !!normal restoration B
        iB=INT(indexB)
        else !!B in the last interval
        iB=INT(indexB)-2
        end if
      else !! B in the first inteval
        iB=INT(indexB)+1
      end if
      fB=indexB-iB !!! fraction part (automatically shifted by +- 1 if needed)

      !! 1) intepolation procedure over B (at grids of X)
      var1=-fB*(fB-1d0)*(fB-2d0)
      var2=3d0*(fB+1d0)*(fB-1d0)*(fB-2d0)
      var3=-3d0*(fB+1d0)*fB*(fB-2d0)
      var4=(fB+1d0)*fB*(fB-1d0)

      interI(0,-5:5)=(var1*gridMain(iX-1,iB-1,-5:5,h)+var2*gridMain(iX-1,iB,-5:5,h)&
            +var3*gridMain(iX-1,iB+1,-5:5,h)+var4*gridMain(iX-1,iB+2,-5:5,h))/6d0
      interI(1,-5:5)=(var1*gridMain(iX,iB-1,-5:5,h)+var2*gridMain(iX,iB,-5:5,h)&
            +var3*gridMain(iX,iB+1,-5:5,h)+var4*gridMain(iX,iB+2,-5:5,h))/6d0
      interI(2,-5:5)=(var1*gridMain(iX+1,iB-1,-5:5,h)+var2*gridMain(iX+1,iB,-5:5,h)&
            +var3*gridMain(iX+1,iB+1,-5:5,h)+var4*gridMain(iX+1,iB+2,-5:5,h))/6d0
      if(iX+2==Nx) then
      interI(3,-5:5)=(/0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp,0_dp/)
      else
      interI(3,-5:5)=(var1*gridMain(iX+2,iB-1,-5:5,h)+var2*gridMain(iX+2,iB,-5:5,h)&
            +var3*gridMain(iX+2,iB+1,-5:5,h)+var4*gridMain(iX+2,iB+2,-5:5,h))/6d0
      end if
    !
    !      !!! linear interpolation procedure
    ! ! 	interI(0,-5:5)=fB*gridMain(iX-1,iB,-5:5,h)+(1d0-fB)*gridMain(iX-1,iB+1,-5:5,h)
    ! ! 	interI(1,-5:5)=fB*gridMain(iX,iB,-5:5,h)+(1d0-fB)*gridMain(iX,iB+1,-5:5,h)
    ! ! 	interI(2,-5:5)=fB*gridMain(iX+1,iB,-5:5,h)+(1d0-fB)*gridMain(iX+1,iB+1,-5:5,h)
    ! ! 	interI(3,-5:5)=fB*gridMain(iX+2,iB,-5:5,h)+(1d0-fB)*gridMain(iX+2,iB+1,-5:5,h)
    end if

    !! 2) intepolation procedure over X (from values computed at B)
	var1=-fX*(fX-1d0)*(fX-2d0)
	var2=3d0*(fX+1d0)*(fX-1d0)*(fX-2d0)
	var3=-3d0*(fX+1d0)*fX*(fX-2d0)
	var4=(fX+1d0)*fX*(fX-1d0)
	ExtractFromGrid=(var1*interI(0,-5:5)+var2*interI(1,-5:5)&
			  +var3*interI(2,-5:5)+var4*interI(3,-5:5))/6d0

  do i=-5,5
   if(ISNAN(ExtractFromGrid(i))) then

    write(*,*) ErrorString('grid extraction produced NaN. EVALUATION STOP',moduleName)
    write(*,*) '----- information on last call -----'
    write(*,*) 'bT=',bT,' i=',i, ' extraction=',ExtractFromGrid(i)
    write(*,*) 'interI=',interI(0:3,i)

    stop
   end if
  end do

end function ExtractFromGrid


!!!! This subroutine test grid by computing the values of function at middle points and compare to exact values
subroutine TestGrid()
    integer::h,h_local,i,j
    real(dp),dimension(-5:5):: test_local

    character*5,dimension(0:5),parameter::xSTR=(/"    0","10^-4"," 0.01","  0.1","  0.8","    1"/)
    character*5,dimension(0:4),parameter::bSTR=(/"    0","  0.1","   1.","   5."," BMAX"/)
    integer,dimension(0:5)::xVal
    integer,dimension(0:4)::bVal

    xVal=(/0,NodeForX(0.0001d0),NodeForX(0.01d0),NodeForX(0.1d0),NodeForX(0.8d0),NX-1/)
    bVal=(/0,NodeForB(0.1d0),NodeForB(1.d0),NodeForB(5.d0),NB/)

   write(*,*) color("  ---  TEST OF GRIDS ---" , c_yellow_bold)
   write(*,*) "* Avarage values of [interpolation/exact-1] in the middle points of grid (max deviation)"
   do h=1,numberOfHadrons
    h_local=hadronsInGRID(h)
    write(*,*) color("  -----  Hadron "//trim(int4ToStr(h)) , c_yellow_bold)

    write(*,&
    '(" |              |     sBar    |     uBar    |     dBar    |   gluon     |      d      |      u      |      s      |")')
    write(*,*)&
       "|----------------------------------------------------------------------------------------------------------------|"
      do i=1,size(xVal)-1
      if(xMin<xVal(i)) then
        write(*,*)"|",xSTR(i-1),"<x<",xSTR(i),"                                                      "
        write(*,*)&
    "|----------------------------------------------------------------------------------------------------------------|"
      do j=1,size(bVal)-1
      test_local=TestPartOfGrid(xVal(i-1),xVal(i),bVal(j-1),bVal(j),h_local)
      write(*,'(" |",A5,"<b<",A5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |")') &
      bSTR(j-1),bSTR(j),test_local(-3:3)
    end do
    write(*,*)&
    "|----------------------------------------------------------------------------------------------------------------|"
    end if
    end do
     write(*,*) color("  -----  test for hadron "//trim(int4ToStr(h))//" complete." , c_yellow_bold)
  end do
end subroutine TestGrid

function TestPartOfGrid(Ix_Low,Ix_high,Ib_low,Ib_high,h)
  integer,intent(in)::Ix_Low,Ix_high,Ib_low,Ib_high,h
  real(dp),dimension(-5:5)::TestPartOfGrid
  real(dp),dimension(1:(Ix_high-Ix_Low),1:(Ib_high-Ib_low),-5:5)::f1,f2
  integer::ix,ib
  real(dp)::x_local,b_local
  real(dp),dimension(-5:5)::comulant

  comulant=0._dp

  do ix=1,Ix_high-Ix_Low
    x_local=XatNode_real(Ix_Low+ix-0.5_dp)
  do ib=1,Ib_high-Ib_Low
    b_local=BatNode_real(Ib_Low+ib-0.5_dp)

    f1(ix,ib,-5:5)=ExtractFromGrid(x_local,b_local,h)
    f2(ix,ib,-5:5)=CxF_compute(x_local,b_local,h,withGluon)
    if(.not.withGluon) then
      f1(ix,ib,0)=1._dp
      f2(ix,ib,0)=1._dp
    end if

    comulant=comulant+abs(f1(ix,ib,-5:5)/f2(ix,ib,-5:5)-1._dp)
  end do
  end do

  TestPartOfGrid=comulant/(Ix_high-Ix_Low)/(Ib_high-Ib_low)

end function TestPartOfGrid
