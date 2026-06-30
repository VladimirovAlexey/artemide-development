!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.05
!
!    Evaluation of the TMD cross-section for DY-like cross-sections
!    This module provide general interface and parralel structure.
!
!    if you use this module please, quote 1706.01473
!
!    ver 3.0: created from v.2.06 (AV, 05.09.2023)
!    ver 3.04: Polishing (AV, 26.06.2026)
!    ver 3.05: The computational part of module is split into TMDX_DY_point and TMDX_DY_bin. (AV, 30.06.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_DY
use aTMDe_Numerics
use aTMDe_IO
use TMDX_DY_point
use TMDX_DY_bin

implicit none
private

!Current version of module
character (len=7),parameter :: moduleName="TMDX-DY"
character (len=5),parameter :: version="v3.05"
!Last appropriate version of constants-file
integer,parameter::inputver=38

integer::outputlevel=2
type(Warning_OBJ)::Warning_Handler

!!!!! numeric precisions
real(dp) :: toleranceGEN=0.0000001d0

!!! number of sections for PT-integral by default
integer::NumPTdefault=6
!!! Partition the qT bins for consecutive evaluation
logical::doPartitioning_byDefault=.false.
real(dp)::MaxQT_range_toPartite=15._dp  !!!!! if the range is bigger than this, it is cut

logical::started=.false.

public:: TMDX_DY_Initialize,TMDX_DY_IsInitialized
public:: TMDX_DY_ResetCounters, TMDX_DY_SetScaleVariation
public:: xSec_DY,xSec_DY_List


contains
  
function TMDX_DY_IsInitialized()
  logical::TMDX_DY_IsInitialized
  TMDX_DY_IsInitialized=started
end function TMDX_DY_IsInitialized

!! Initialization of the package
subroutine TMDX_DY_Initialize(file,prefix)
  character(len=*),intent(in)::file
  character(len=*),intent(in),optional::prefix
  character(len=:),allocatable::path
  logical::initRequired
  integer::i,FILEver,messageTrigger
  !$ integer:: omp_get_thread_num

  if(started) return

  if(present(prefix)) then
    path=trim(adjustl(prefix))//trim(adjustl(file))
  else
    path=trim(adjustl(file))
  end if

  OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

  !!! Check the file version
  call MoveTO(51,'*0   ')
  call MoveTO(51,'*A   ')
  call MoveTO(51,'*p1  ')
  read(51,*) FILEver
  if(FILEver<inputver) then
    write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
    write(*,*) '             Update the const-file with artemide.setup'
    write(*,*) '  '
    CLOSE (51, STATUS='KEEP')
    error stop
  end if

  !!! Fill the message system
  call MoveTO(51,'*p2  ')
  read(51,*) outputLevel
  if(outputLevel>2) write(*,*) '--------------------------------------------- '
  if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
  call MoveTO(51,'*p3  ')
  read(51,*) messageTrigger

  !$    if(outputLevel>1) write(*,*) '    artemide.TMDX_DY: parallel evaluation of cross-sections is to be used'
  !$    call MoveTO(51,'*C   ')
  !$    call MoveTO(51,'*p1  ')
  !$    read(51,*) i
  !$    call OMP_set_num_threads(i)
  !$    if(outputLevel>1) write(*,*) '    artemide.TMDX_DY: number of threads for parallel evaluation is set to ', i

  !!! go to section TMD-DY
  call MoveTO(51,'*9   ')
  call MoveTO(51,'*p1  ')
  read(51,*) initRequired
  if(.not.initRequired) then
    if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
    started=.false.
    CLOSE (51, STATUS='KEEP')
    return
  end if

  !!!------ parameters of numerics
  call MoveTO(51,'*B   ')
  call MoveTO(51,'*p1  ')
  read(51,*) toleranceGEN
  call MoveTO(51,'*p3  ')
  read(51,*) NumPTdefault
  call MoveTO(51,'*p6  ')
  read(51,*) doPartitioning_byDefault
  call MoveTO(51,'*p8  ')
  read(51,*) MaxQT_range_toPartite

  CLOSE (51, STATUS='KEEP')

  Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

  !$     if(outputLevel>2) write(*,*) '------TEST OF PARALLEL PROCESSING ----------'
  !$OMP PARALLEL
  !$     if(outputLevel>2) write(*,*) '   artemide.TMDX_DY:thread num ',  omp_get_thread_num(), ' ready.'
  !$OMP END PARALLEL

  if(.not.TMDX_DY_bin_IsInitialized()) then
    if(outputLevel>1) write(*,*) '.. initializing TMDX_DY_bin (from ',moduleName,')'
    if(present(prefix)) then
      call TMDX_DY_bin_Initialize(file,prefix)
    else
      call TMDX_DY_bin_Initialize(file)
    end if
  end if

  started=.true.

  write(*,*)  color('----- arTeMiDe.TMD_DY '//trim(version)//': .... initialized',c_green)
end subroutine TMDX_DY_Initialize


!!!!! Reset counters here and in other TMDX sub modules
subroutine TMDX_DY_ResetCounters()
    call TMDX_DY_1pt_ResetCounters()
    call TMDX_DY_bin_ResetCounters()
    call Warning_Handler%Reset()
end subroutine TMDX_DY_ResetCounters

!!!!Call this after TMD initialization but before NP, and X parameters
subroutine TMDX_DY_SetScaleVariation(c2_in)
  real(dp),intent(in)::c2_in

  call TMDX_DY_1pt_SetScaleVariation(c2_in)

end subroutine TMDX_DY_SetScaleVariation

!!! function determines the best value of PT-sections from PT-bin size, and Q
!!! it is determined by formula Q/PT< val/ (2 k) => def+2K
function NumPT_auto(dPT,Q)
  real,parameter::val=40.
  real(dp)::dPT,Q,rat
  integer::i,NumPT_auto
  rat=Q/dPT

  if(rat>40.) then
      NumPT_auto=NumPTdefault
      return
  else
      do i=1,25
          if(rat>(40d0/2d0/i)) then
              NumPT_auto=NumPTdefault+2*i
              return
          end if
      end do
  end if
  if(outputlevel>1) then
    write(*,*) WarningString('Fail to automatically determine number of Pt-section for a bin.',moduleName)
    write(*,*) '  ... Possibly Pt-bin is too large', dPT
  end if
  NumPT_auto=NumPTdefault+12
end function NumPT_auto

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!THE MAIN INTERFACE TO CROSS-SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! interface for array,s,array,array,array,logical,optional, optional
subroutine xSec_DY(X,process,s,qT,Q,y,includeCuts,CutParameters,Num)
!   function xSec_DY(process,s,qT,Q,y,includeCuts,CutParameters,Num)
  integer,intent(in),dimension(1:4)::process            !the number of process
  real(dp),intent(in)::s                          !Mandelshtam s
  real(dp),intent(in),dimension(1:2)::qT           !(qtMin,qtMax)
  real(dp),intent(in),dimension(1:2)::Q            !(Qmin,Qmax)
  real(dp),intent(in),dimension(1:2)::y             !(ymin,ymax)
  logical,intent(in)::includeCuts                !include cuts
  real(dp),intent(in),dimension(1:4),optional::CutParameters    !(p1,p2,eta1,eta2)
  integer,intent(in),optional::Num                !number of sections

  real(dp)::X

  integer::nn
  real(dp),dimension(1:4)::CutParam

  !! determine number of sections
  if(present(Num)) then
    nn=Num
  else
    nn=NumPT_auto(qT(2)-qT(1),(Q(2)+Q(1))/2d0)
  end if

  !!! determine cut parameters
  if(includeCuts) then
    if(present(CutParameters)) then
      CutParam=CutParameters
    else
      error stop ErrorString('called includeCuts=true, while CutParameters are undefined',moduleName)
    end if
  else
    CutParam=(/0d0,0d0,0d0,0d0/)
  end if

  !!!! evaluation
  X=Xsec_PTint_Qint_Yint(process,includeCuts,CutParam,&
                  s,qT(1),qT(2),Q(1),Q(2),y(1),y(2),nn)

end subroutine xSec_DY

!!!!! this is the most used subroutine
!!!!! It takes lists of bin-sizes and sends them in parallel mode to computation
!!!!! If the partitioning is .true. it attempts to find consecutive lists of pt-bins (all same except the pt-bin)
subroutine xSec_DY_List(X,process,s,qT,Q,y,includeCuts,CutParameters,Num,doPartitioning)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                !Mandelshtam s
  real(dp),intent(in),dimension(:,:)::qT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(:,:)::Q                !(Qmin,Qmax)
  real(dp),intent(in),dimension(:,:)::y                !(ymin,ymax)
  logical,intent(in),dimension(:)::includeCuts        !include cuts
  real(dp),intent(in),dimension(:,:)::CutParameters            !(p1,p2,eta1,eta2)
  integer,intent(in),dimension(:),optional::Num        !number of sections
  logical,intent(in),optional::doPartitioning        !if .true. intents to split the pt-integrate to sebsections
  real(dp),dimension(:),intent(out)::X
  integer :: i,length
  integer,allocatable::nn(:)

  logical::doP
  integer::k,j,numberOfP,listOfParts(1:size(s)),n_private
  integer,allocatable:: partI1(:),partSize(:)

if(.not.started) error stop ErrorString('The module is not initialized. Check INI-file.',moduleName)

  length=size(s)

!!! cheking sizes
if(size(X)/=length) then
  error stop ErrorString('xSec_DY_List: sizes of xSec and s lists are not equal.',moduleName)
end if
if(size(process,1)/=length) then
  error stop ErrorString('xSec_DY_List: sizes of process and s lists are not equal.',moduleName)
end if
if(size(qT,1)/=length) then
  error stop ErrorString('xSec_DY_List: sizes of qT and s lists are not equal.',moduleName)
end if
if(size(y,1)/=length) then
  error stop ErrorString('xSec_DY_List: sizes of y and s lists are not equal.',moduleName)
end if
if(size(Q,1)/=length) then
  error stop ErrorString('xSec_DY_List: sizes of Q and s lists are not equal.',moduleName)
end if
if(size(includeCuts)/=length) then
  error stop ErrorString('xSec_DY_List: sizes of includeCuts and s lists are not equal.',moduleName)
end if
if(size(CutParameters,1)/=length) then
  error stop ErrorString('xSec_DY_List: sizes of CutParameters and s lists are not equal.',moduleName)
end if
if(size(process,2)/=4) then
  error stop ErrorString('xSec_DY_List: process list must be (:,1:4).',moduleName)
end if
if(size(qT,2)/=2) then
  error stop ErrorString('xSec_DY_List: qt list must be (:,1:2).',moduleName)
end if
if(size(y,2)/=2) then
  error stop ErrorString('xSec_DY_List: y list must be (:,1:2).',moduleName)
end if
if(size(Q,2)/=2) then
  error stop ErrorString('xSec_DY_List: Q list must be (:,1:2).',moduleName)
end if

  if(present(doPartitioning)) then
    doP=doPartitioning
  else
    doP=doPartitioning_byDefault
  end if

  if(doP) then
  !!!! attempt to make partitioning into qT sectors
  !!!! basically, I run though the whole list and compare the terms if all (except qT) are the same, they are marked by the same number
  !!!! Only consecutive ranges are marked, and there is also upper cut
    k=1
    listOfParts(1)=k
    do i=2,length
      if(&
      qT(i,2)>MaxQT_range_toPartite &  !!! check the max size
      .or.(abs(y(i,1)-y(i-1,1))>toleranceGEN) .or. (abs(y(i,2)-y(i-1,2))>toleranceGEN) &
      .or.(abs(Q(i,1)-Q(i-1,1))>toleranceGEN) .or. (abs(Q(i,2)-Q(i-1,2))>toleranceGEN) &
      .or.(process(i,1)/=process(i-1,1)) .or. (process(i,2)/=process(i-1,2)) .or. (process(i,3)/=process(i-1,3)) &
      .or.(process(i,4)/=process(i-1,4)) .or. (abs(s(i-1)-s(i))>toleranceGEN) .or. (includeCuts(i-1).neqv.includeCuts(i)) &
      .or.(abs(CutParameters(i,1)-CutParameters(i-1,1))>toleranceGEN) &
      .or.(abs(CutParameters(i,2)-CutParameters(i-1,2))>toleranceGEN) &
      .or.(abs(CutParameters(i,3)-CutParameters(i-1,3))>toleranceGEN) &
      .or.(abs(CutParameters(i,4)-CutParameters(i-1,4))>toleranceGEN) &
      .or.(qT(i-1,2)>qT(i,1)+toleranceGEN) & !!!!! the bins are successive
      ) k=k+1

      listOfParts(i)=k
    end do
    numberOfP=k
    !!!! partitions is a list that contain initial numbers of each partition
    allocate(partI1(1:numberOfP))
    partI1(1)=1
    k=2
    do i=2,length
      if(listOfParts(i)/=listOfParts(i-1)) then
        partI1(k)=i
        k=k+1
      end if
    end do

    !!!! partSize is a list that contain sizes of partI1
    allocate(partSize(1:numberOfP))

    do i=1,numberOfP-1
      partSize(i)=partI1(i+1)-partI1(i)
    end do
    partSize(numberOfP)=length-partI1(numberOfP)+1


    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(n_private)
    do i=1,numberOfP
      !!!!! if bin is not partitioned, it is computed usually
      if(partSize(i)==1) then
        n_private=NumPT_auto(qT(partI1(i),2)-qT(partI1(i),1),(Q(partI1(i),2)+Q(partI1(i),1))/2d0)
        X(partI1(i))=Xsec_PTint_Qint_Yint(&
              process(partI1(i),1:4),includeCuts(partI1(i)),CutParameters(partI1(i),1:4),&
              s(partI1(i)),qT(partI1(i),1),qT(partI1(i),2),&
              Q(partI1(i),1),Q(partI1(i),2),y(partI1(i),1),y(partI1(i),2),n_private)
      else
        !!! actual computation
        X(partI1(i):partI1(i)+partSize(i)-1)=Xsec_PTspectrum_Qint_Yint(&
            process(partI1(i),1:4),includeCuts(partI1(i)),CutParameters(partI1(i),1:4),&
            s(partI1(i)),&
            qT(partI1(i):partI1(i)+partSize(i)-1,1),qT(partI1(i):partI1(i)+partSize(i)-1,2),&
            Q(partI1(i),1),Q(partI1(i),2),y(partI1(i),1),y(partI1(i),2))

      end if
    end do
    !$OMP END PARALLEL DO
  else
  !!!!!--------------- compute in the usual way
    allocate(nn(1:length))
    if(present(Num)) then
        if(size(Num,1)/=length) then
      error stop ErrorString('xSec_DY_List: sizes of Num and s lists are not equal.',moduleName)
    end if
    nn=Num
    else
        do i=1,length
            nn(i)=NumPT_auto(qT(i,2)-qT(i,1),(Q(i,2)+Q(i,1))/2d0)
        end do
    end if
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED)

      do i=1,length
        X(i)=Xsec_PTint_Qint_Yint(process(i,1:4),includeCuts(i),CutParameters(i,1:4),&
                s(i),qT(i,1),qT(i,2),Q(i,1),Q(i,2),y(i,1),y(i,2),nn(i))
      end do
    !$OMP END PARALLEL DO
    deallocate(nn)
  end if
end subroutine xSec_DY_List
  
end module TMDX_DY
