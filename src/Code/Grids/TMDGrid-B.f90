!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	This file contains the part of the code, which is common for all TMD-evaluation modules that operates at twist-2
!	It shares the variables with the module where it is inlucded (as a text)
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part is devoted to the Grid evaluation
!
!	v.2.00 Large-b evaluation changed AV (25.03.2019)
!	v.2.01 Cosmetic changes		  AV (25.04.2019)
!
!				A.Vladimirov (08.10.2018)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------
!Used global variables:
! outputlevel, moduleName
! numberOfHadrons, hadronsInGRID
! lambdaNP_grid, lambdaNP, FNP
! aTMDe_Numerics, IO_functions functions
! + variables defined in TMDGrid-B-VAR.f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Griding functions  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

!!! this subroutine create grid. It is called only from the uTMDPDF_SetLambdaNP function.
subroutine MakeGrid()
    real(dp):: x_local,b_local

    !!!size of grid is power of 2, to make it testable and reusable. order is the current power.
    integer:: iX,iB,j,h

    real(dp),dimension(-5:5)::checkValue,dummyfNP
    real(dp)::maxValue,eps,epsB,epsX,B_BAD,x_BAD,eps_av,time1,time2
    integer::grid_size
    real(dp),dimension(0:GridSizeX,1:2,-5:5,1:numberOfHadrons)::smallGrid
    real(dp)::v1_BAD,v2_BAD

    if(numberOfHadrons==1 .and. hadronsInGRID(1)==0) return

    !!!!!!!!!!!!!**************************************************************************************************************
    !!initial message
    call cpu_time(time1)
    if(outputlevel>1) then
        if(numberOfHadrons>1) then
            write(*,*) 'arTeMiDe.',moduleName,' starts to build grids for',numToStr(numberOfHadrons), ' hadrons. At ',time1
        else
            write(*,*) 'arTeMiDe.',moduleName,' starts to build grid. At ',time1            
        end if
    end if
    
    !!!!!!!!!!!!!**************************************************************************************************************
    !!main part
    
    !!!save the parameters on which the grid is build
    lambdaNP_grid=lambdaNP

    bFactor=bGrid_Max/(EXP(slope)-1d0) !!! often used combination
    slopeFactor=slope/REAL(GridSizeB)
    !  goto 111 
    !!!We use the following grid
    !! b=(e^{a i/imax}-1)/(e^a-1)*b_Max 
    !! where a=slope. If slope higher more points concentrated at b->0
    !! For x we use the grid
    !! x=(xMin)^(i/imax)
    do h=1,numberOfHadrons
        do iX=0,GridSizeX
            do iB=0,GridSizeB
                x_local=xGrid_Min**(REAL(iX)/REAL(GridSizeX))
                b_local=(EXP(slopeFactor*REAL(iB))-1d0)*bFactor

                if(x_local==1d0) then
                    gridMain(iX,iB,-5:5,h)=0d0
                else
                    if(withGluon) then
                        gridMain(iX,iB,-5:5,h)=x_local*Common_lowScale50(x_local,b_local,hadronsInGRID(h))
                    else
                        gridMain(iX,iB,-5:5,h)=x_local*Common_lowScale5(x_local,b_local,hadronsInGRID(h))
                    end if
                end if
            end do    
        end do
        if(outputLevel>1 .and. numberOfHadrons>1) write(*,'(" ",A,": Grid for hadron ",I3," is done")') moduleName,hadronsInGRID(h)
    end do
    !!!!!!!!!!!!!**************************************************************************************************************
    !!!!!!!!!!!!!**************************************************************************************************************
    !!!!!!!!!!!!!*****************************TALE PART************************************************************************
    !!!!!!!!!!!!!**************************************************************************************************************

    !!!The high B (b>bGrid_Max) tales are approximated as value(bGrid_Max)/fNP(x,x,bGrid_Max) * fNP(x,x,b)
    !!! this approximation is exact if coeff-function freeze at large b.
111 do h=1,numberOfHadrons
        do iX=1,GridSizeX
            x_local=xGrid_Min**(REAL(iX)/REAL(GridSizeX))
            b_local=bGrid_Max

            if(withGluon) then
                boundaryValues(iX,-5:5,h)=x_local*Common_lowScale50(x_local,b_local,hadronsInGRID(h))
            else
                boundaryValues(iX,-5:5,h)=x_local*Common_lowScale5(x_local,b_local,hadronsInGRID(h))
                boundaryValues(iX,0,h)=0d0
            end if

            dummyfNP=FNP(x_local,x_local,b_local,hadronsInGRID(h),lambdaNP_grid)

            !!!! if fNP=0, here, we check is C x f fNP=0, if so, we save 1, otherwise exception
            do j=-5,5
                if(dummyfNP(j)==0d0) then
                if(boundaryValues(iX,j,h)/=0d0 .and. outputLevel>0) then  
                    call Warning_Raise('error in evaluation of boundary for grid',messageCounter,messageTrigger,moduleName)
                    if(messageCounter<messageTrigger) then
                        write(*,*) '----- information on last call -----'
                        write(*,*) 'fNP=',dummyfNP(j),'boundary=',boundaryValues(iX,j,h), 'for (h,f)=(',hadronsInGRID(h),j,')'
                    write(*,*) 'Continue with value=1'
                    end if
                end if
                boundaryValues(iX,j,h)=1d0

                else
                    !!!normal case
                    boundaryValues(iX,j,h)=boundaryValues(iX,j,h)/dummyfNP(j)
                end if
            end do

        end do
    end do
    
    !!!!!!!!!!!!!**************************************************************************************************************
    !!check the grid
    !! we pass though grid restoring a point from near-by points
    !! in this way we estimate error
    if(outputLevel>2) then
        do h=1,numberOfHadrons
            epsB=0d0
            epsX=0d0
            eps_av=0d0
            grid_size=0
            !check b-grid
            do iX=0,GridSizeX
                maxValue=MAXVAL(ABS(gridMain(iX,2:GridSizeB,1:3,h)))
                if(maxValue>0d0) then
                    do iB=2,GridSizeB-4
                        checkValue=((-gridMain(iX,iB,-5:5,h)+4d0*gridMain(iX,iB+1,-5:5,h)+4d0*gridMain(iX,iB+3,-5:5,h)&
                        -gridMain(iX,iB+4,-5:5,h))/6d0-gridMain(iX,iB+2,-5:5,h))

                        if(.not.withGluon) checkValue(0)=0

                        eps=MAXVAL(ABS(checkValue))/maxValue
                        eps_av=eps_av+checkValue(MAXLOC(ABS(checkValue),1)-6)/maxValue      
                        grid_size=grid_size+1

                        if(eps>epsB) then
                            epsB=eps
                            B_BAD=(EXP(slopeFactor*REAL(iB+2))-1d0)*bFactor
                            x_BAD=xGrid_Min**(REAL(iX)/REAL(GridSizeX))
                        end if
                    end do    
                end if
            end do
            !     write(*,*) 'The worst B reconstraction with ',epsB, 'at (x,b)=',x_BAD,B_BAD

            !   check x-grid
            do iB=0,GridSizeB
                if(withGluon) then 
                    maxValue=MAXVAL(ABS(gridMain(0:GridSizeX,iB,0:3,h)))
                else
                    maxValue=MAXVAL(ABS(gridMain(0:GridSizeX,iB,1:3,h)))
                end if
                if(maxValue>0d0) then
                    do iX=0,GridSizeX-4
                        checkValue=((-gridMain(iX,iB,-5:5,h)+4d0*gridMain(iX+1,iB,-5:5,h)+4d0*gridMain(iX+3,iB,-5:5,h)&
                        -gridMain(iX+4,iB,-5:5,h))/6d0-gridMain(iX+2,iB,-5:5,h))

                        if(.not.withGluon) checkValue(0)=0

                        eps=MAXVAL(ABS(checkValue))/maxValue
                        eps_av=eps_av+eps
                        grid_size=grid_size+1
                        if(eps>epsX) then
                            epsX=eps
                            B_BAD=(EXP(slopeFactor*REAL(iB))-1d0)*bFactor
                            x_BAD=xGrid_Min**(REAL(iX+2)/REAL(GridSizeX))
                        end if
                    end do    
                end if
            end do
            !     write(*,*) 'The worst X reconstraction with ',epsX, 'at (x,b)=',x_BAD,B_BAD

            write(*,'(" ",A,": Grid (for hadron ",I3,") av.badness =",ES12.3)')&
            moduleName,hadronsInGRID(h),eps_av/REAL(grid_size)

        end do
    end if
    
    !!!!!!!!!!!!!**************************************************************************************************************
    !!final message
    call cpu_time(time2)
    if(outputlevel>1) then
        if(numberOfHadrons>1) then
            write(*,'(" ",A,": Grids are built  (",I5," x",I5,")  calc.time=",F6.2,"s. ")')&
            moduleName, GridSizeX,GridSizeB, time2-time1
        else
            write(*,'(" ",A,": Grid is built  (",I5," x",I5,")  calc.time=",F6.2,"s. ")')&
            moduleName, GridSizeX,GridSizeB, time2-time1
        end if
    end if

end subroutine MakeGrid
  
function ExtractFromGrid(x,bT,hadron)
  real(dp),intent(in)::x,bT
  integer,intent(in)::hadron
  real(dp),dimension(-5:5)::ExtractFromGrid
  real(dp),dimension(0:3,-5:5):: interI
  real(dp),dimension(-5:5)::dummyfNP
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
  
  if(x<xGrid_Min) then
   write(*,*) ErrorString('The TMD with x ='//numToStr(x)//'is called. Current grid size is up to '//&
   numToStr(xGrid_Min)//'. Enlarge boundaries.',moduleName)
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
  indexX=GridSizeX*LOG(x)/LOG(xGrid_Min)
  iX=INT(indexX)
    if(iX>0) then
     if(iX<GridSizeX-1) then !!normal restoration X
      iX=INT(indexX)
     else !!X in the last interval
      iX=INT(indexX)-2      
      end if
    else !! X in the first inteval
      iX=INT(indexX)+1
    end if
    fX=indexX-iX !!! fraction part (automatically shifted by +- 1 if needed)
  
   if(bT>bGrid_Max) then
    
   !!!!here will code for large b
   
   !!! we interpolate over x the boundary value
    var1=-fX*(fX-1d0)*(fX-2d0)
    var2=3d0*(fX+1d0)*(fX-1d0)*(fX-2d0)
    var3=-3d0*(fX+1d0)*fX*(fX-2d0)
    var4=(fX+1d0)*fX*(fX-1d0)
    ExtractFromGrid=FNP(x,x,bT,hadronsInGRID(h),lambdaNP_grid)*&
		    (var1*boundaryValues(iX-1,-5:5,h)+var2*boundaryValues(iX,-5:5,h)&
			  +var3*boundaryValues(iX+1,-5:5,h)+var4*boundaryValues(iX+2,-5:5,h))/6d0/x
   
   else!!! b inside the main region
    
    indexB=LOG(1d0+bT/bFactor)/slopeFactor
    iB=INT(indexB)
    if(iB>0) then
     if(iB<GridSizeB-1) then !!normal restoration B
      iB=INT(indexB)
     else !!B in the last interval
      iB=INT(indexB)-2      
      end if
    else !! B in the first inteval
      iB=INT(indexB)+1
    end if
    fB=indexB-iB !!! fraction part (automatically shifted by +- 1 if needed)
     
     
     
     !! intepolation procedure
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
	interI(3,-5:5)=(var1*gridMain(iX+2,iB-1,-5:5,h)+var2*gridMain(iX+2,iB,-5:5,h)&
	      +var3*gridMain(iX+2,iB+1,-5:5,h)+var4*gridMain(iX+2,iB+2,-5:5,h))/6d0
	
     !!! linear interpolation procedure
! 	interI(0,-5:5)=fB*gridMain(iX-1,iB,-5:5,h)+(1d0-fB)*gridMain(iX-1,iB+1,-5:5,h)
! 	interI(1,-5:5)=fB*gridMain(iX,iB,-5:5,h)+(1d0-fB)*gridMain(iX,iB+1,-5:5,h)
! 	interI(2,-5:5)=fB*gridMain(iX+1,iB,-5:5,h)+(1d0-fB)*gridMain(iX+1,iB+1,-5:5,h)
! 	interI(3,-5:5)=fB*gridMain(iX+2,iB,-5:5,h)+(1d0-fB)*gridMain(iX+2,iB+1,-5:5,h)
	
	
	var1=-fX*(fX-1d0)*(fX-2d0)
	var2=3d0*(fX+1d0)*(fX-1d0)*(fX-2d0)
	var3=-3d0*(fX+1d0)*fX*(fX-2d0)
	var4=(fX+1d0)*fX*(fX-1d0)
	ExtractFromGrid=(var1*interI(0,-5:5)+var2*interI(1,-5:5)&
			  +var3*interI(2,-5:5)+var4*interI(3,-5:5))/6d0/x
	
     end if
	
    
  do i=-5,5
   if(ISNAN(ExtractFromGrid(i))) then
    
    write(*,*) ErrorString('grid extraction produced NaN. EVALUATION STOP',moduleName)
    write(*,*) '----- information on last call -----'
    write(*,*) 'bT=',bT,' i=',i, ' extraction=',ExtractFromGrid(i)
    write(*,*) 'interI=',interI(0:3,i)
    
    dummyfNP=FNP(x,x,bT,hadronsInGRID(h),lambdaNP_grid)
    write(*,*) 'fNP=', dummyfNP(i)
    
    stop
   end if
  end do
			  
end function ExtractFromGrid 


!!!! This subroutine test grid by computing the values of function at middle points and compare to exact values
subroutine TestGrid_2()
    integer::h,h_local,i,j
    real(dp),dimension(-5:5):: test_local

    character*5,dimension(0:5),parameter::xSTR=(/"    0","10^-4"," 0.01","  0.1","  0.8","    1"/)
    character*5,dimension(0:4),parameter::bSTR=(/"    0","  0.1","   1.","   5."," BMAX"/)
    integer,dimension(0:5)::xVal
    integer,dimension(0:4)::bVal

    xVal=(/GridSizeX-1,INT(GridSizeX*LOG(0.0001d0)/LOG(xGrid_Min)),INT(GridSizeX*LOG(0.01d0)/LOG(xGrid_Min)),&
    INT(GridSizeX*LOG(0.1d0)/LOG(xGrid_Min)),INT(GridSizeX*LOG(0.8d0)/LOG(xGrid_Min)),0/)
    bVal=(/0,INT(LOG(1d0+0.1d0/bFactor)/slopeFactor),INT(LOG(1d0+1.d0/bFactor)/slopeFactor),&
    INT(LOG(1d0+5.d0/bFactor)/slopeFactor),GridSizeB/)

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
      if(xGrid_Min<xVal(i)) then
        write(*,*)"|",xSTR(i-1),"<x<",xSTR(i),"                                                      "
        write(*,*)&
    "|----------------------------------------------------------------------------------------------------------------|"
      do j=1,size(bVal)-1
      test_local=TestPartOfGrid(xVal(i),xVal(i-1),bVal(j-1),bVal(j),h_local)
      write(*,'(" |",A5,"<b<",A5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |",E12.5," |")') &
      bSTR(j-1),bSTR(j),test_local(-3:3)
    end do
    write(*,*)&
    "|----------------------------------------------------------------------------------------------------------------|"
    end if
    end do
     write(*,*) color("  -----  test for hadron "//trim(int4ToStr(h))//" complete." , c_yellow_bold)
  end do
end subroutine TestGrid_2

function TestPartOfGrid(Ix_Low,Ix_high,Ib_low,Ib_high,h)
  integer,intent(in)::Ix_Low,Ix_high,Ib_low,Ib_high,h
  real(dp),dimension(-5:5)::TestPartOfGrid
  real(dp),dimension(1:(Ix_high-Ix_Low),1:(Ib_high-Ib_low),-5:5)::f1,f2
  integer::ix,ib
  real(dp)::x_local,b_local
  real(dp),dimension(-5:5)::comulant

  comulant=0._dp

  do ix=1,Ix_high-Ix_Low

  x_local=xGrid_Min**((Ix_high-ix+0.5_dp)/REAL(GridSizeX))

  do ib=1,Ib_high-Ib_Low
  b_local=(EXP(slopeFactor*(Ib_Low+ib+0.5_dp))-1d0)*bFactor
  !if(ix==1) write(*,*) "---->",Ib_Low,Ib_high,b_local

    f1(ix,ib,-5:5)=ExtractFromGrid(x_local,b_local,h)
    f2(ix,ib,-5:5)=Common_lowScale5(x_local,b_local,h)
    if(.not.withGluon) then
      f1(ix,ib,0)=1._dp
      f2(ix,ib,0)=1._dp
    end if

    comulant=comulant+abs(f1(ix,ib,-5:5)/f2(ix,ib,-5:5)-1._dp)
  end do
  end do

  TestPartOfGrid=comulant/(Ix_high-Ix_Low)/(Ib_high-Ib_low)

end function TestPartOfGrid
