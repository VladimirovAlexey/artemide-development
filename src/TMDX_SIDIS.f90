!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.05
!
!    Evaluation of the TMD cross-section for SIDIS-like cross-sections
!
!    if you use this module please, quote 1912.06532
!
!    ver 1.2: release (AV, 15.12.2017)
!    ver 1.32: part of functions migrated to TMDF, rest updated (AV, 16.08.2018)
!    ver 2.02:                            (AV,16.08.2019)
!    ver 3.03:                            (AV,14.02.2026)
!    ver 3.04:                            (AV,26.06.2026)
!    ver 3.05: Parts split into _point and _bin (AV,02.07.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_SIDIS
use aTMDe_Numerics
use aTMDe_IO
use TMDX_SIDIS_point
use TMDX_SIDIS_bin

implicit none
private

!Current version of module
character (len=10),parameter :: moduleName="TMDX-SIDIS"
character (len=5),parameter :: version="v3.05"
!Last appropriate version of constants-file
integer,parameter::inputver=39

real(dp) :: toleranceGEN=0.0000001d0

integer::outputlevel=2
type(Warning_OBJ)::Warning_Handler

logical::started=.false.


!!! number of sections for PT-integral by default
integer::NumPTdefault=6

logical::doPartitioning_byDefault=.false.
real(dp)::MaxQT_range_toPartite=18._dp  !!!!! if the range is bigger than this, it is cut


public::TMDX_SIDIS_Initialize,TMDX_SIDIS_IsInitialized,TMDX_SIDIS_ResetCounters,TMDX_SIDIS_SetScaleVariation

public::xSec_SIDIS,xSec_SIDIS_List



contains

function TMDX_SIDIS_IsInitialized()
    logical::TMDX_SIDIS_IsInitialized
    TMDX_SIDIS_IsInitialized=started
end function TMDX_SIDIS_IsInitialized

!! Initialization of the package
subroutine TMDX_SIDIS_Initialize(file,prefix)
  character(len=*),intent(in)::file
  character(len=*),intent(in),optional::prefix
  character(len=:),allocatable::path
  logical::initRequired
  character(len=8)::orderMain
  integer::i,j,FILEver,messageTrigger
  !$ integer:: omp_get_thread_num

  if(started) return

  if(present(prefix)) then
      path=trim(adjustl(prefix))//trim(adjustl(file))
  else
      path=trim(adjustl(file))
  end if

  OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

  call MoveTO(51,'*0   ')
  call MoveTO(51,'*A   ')
  call MoveTO(51,'*p1  ')
  read(51,*) FILEver
  if(FILEver<inputver) then
      CLOSE (51, STATUS='KEEP')
      error stop ErrorString('const-file version is too old. Update with artemide.setup.',moduleName)
  end if
  call MoveTO(51,'*p2  ')
  read(51,*) outputLevel
  if(outputLevel>2) write(*,*) '--------------------------------------------- '
  if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
  call MoveTO(51,'*p3  ')
  read(51,*) messageTrigger

  !$    if(outputLevel>1) write(*,*) '    artemide.TMDX_SIDIS: parallel evaluation of cross-sections is to be used'
  !$    call MoveTO(51,'*C   ')
  !$    call MoveTO(51,'*p1  ')
  !$    read(51,*) i
  !$    call OMP_set_num_threads(i)
  !$    if(outputLevel>1) write(*,*) '    artemide.TMDX_SIDIS: number of threads for parallel evaluation is set to ', i

  call MoveTO(51,'*10  ')
  call MoveTO(51,'*p1  ')
  read(51,*) initRequired
  if(.not.initRequired) then
      if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
      started=.false.
      CLOSE (51, STATUS='KEEP')
      return
  end if

  call MoveTO(51,'*B   ')
  call MoveTO(51,'*p1  ')
  read(51,*) toleranceGEN
  call MoveTO(51,'*p3  ')
  read(51,*) NumPTdefault
  call MoveTO(51,'*p5  ')
  read(51,*) doPartitioning_byDefault
  call MoveTO(51,'*p7  ')
  read(51,*) MaxQT_range_toPartite

  CLOSE (51, STATUS='KEEP')

  Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

  !$     if(outputLevel>2) write(*,*) '------TEST OF PARALLEL PROCESSING ----------'
  !$OMP PARALLEL
  !$     if(outputLevel>2) write(*,*) '   artemide.TMDX_SIDIS:thread num ',  omp_get_thread_num(), ' ready.'
  !$OMP END PARALLEL

  if(.not.TMDX_SIDIS_bin_IsInitialized()) then
    if(outputLevel>1) write(*,*) '.. initializing TMDX_SIDIS_bin (from ',moduleName,')'
    if(present(prefix)) then
      call TMDX_SIDIS_bin_Initialize(file,prefix)
    else
      call TMDX_SIDIS_bin_Initialize(file)
    end if
  end if

  started=.true.
  write(*,*)  color('----- arTeMiDe.TMD_SIDIS '//trim(version)//': .... initialized',c_green)
end subroutine TMDX_SIDIS_Initialize

!!!!! Reset counters here and in other TMDX sub modules
subroutine TMDX_SIDIS_ResetCounters()
    call TMDX_SIDIS_1pt_ResetCounters()
    call TMDX_SIDIS_bin_ResetCounters()
    call Warning_Handler%Reset()
end subroutine TMDX_SIDIS_ResetCounters

!!!!Call this after TMD initialization but before NP, and X parameters
subroutine TMDX_SIDIS_SetScaleVariation(c2_in)
  real(dp),intent(in)::c2_in

  call TMDX_SIDIS_1pt_SetScaleVariation(c2_in)

end subroutine TMDX_SIDIS_SetScaleVariation

!!! function determines the best value of PT-sections from PT-bin size, and Q
!!! it is determined by formula Q/PT< val/ (2 k) => def+2K
function NumPT_auto(dPT,Q)
real(dp),parameter::val=40d0
real(dp),intent(in)::dPT,Q
integer::NumPT_auto

real(dp)::rat
integer::i

rat=Q/dPT

if(rat>40d0) then
  NumPT_auto=NumPTdefault
  return
else
  do i=1,5
    if(rat>(40d0/2d0/i)) then
      NumPT_auto=NumPTdefault+2*i
      return
    end if
  end do
end if
if(outputlevel>1) then
write(*,*) WarningString('Fail to automatically determine number of Pt-section for a bin.',moduleName)
write(*,*) '>>  Possibly Pt-bin is too large', dPT
end if
NumPT_auto=NumPTdefault+12

end function NumPT_auto


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  MAIN INTERFACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! single value interface
subroutine xSec_SIDIS(xx,process,s,pT,z,x,Q,doCut,Cuts,masses)
  integer,intent(in),dimension(1:4)::process            !the number of process
  real(dp),intent(in)::s                                !Mandelshtam s
  real(dp),intent(in),dimension(1:2)::pT                !(qtMin,qtMax)
  real(dp),intent(in),dimension(1:2)::z                 !(zmin,zmax)
  real(dp),intent(in),dimension(1:2)::x                 !(xmin,xmax)
  real(dp),intent(in),dimension(1:2)::Q                 !(Qmin,Qmax)
  logical,intent(in)::doCut                             !triger cuts
  real(dp),intent(in),dimension(1:4)::Cuts              !(ymin,ymax,W2min,W2max)
  real(dp),intent(in),dimension(1:2),optional::masses   !(mass_target,mass_product)GeV
  real(dp),intent(out)::xx

  integer::Num

if(.not.started) error stop ErrorString('The module is not initialized. Check INI-file.',moduleName)

  Num=NumPT_auto(pt(2)-pt(1),(Q(1)+Q(2))/2d0)

  if(PRESENT(masses)) then
    xx=Xsec_PTint_Qint_Xint_Zint(process,s,z(1),z(2),x(1),x(2),Q(1),Q(2),pT(1),pT(2),doCut,Cuts,Num,masses(1),masses(2))
  else
    xx=Xsec_PTint_Qint_Xint_Zint(process,s,z(1),z(2),x(1),x(2),Q(1),Q(2),pT(1),pT(2),doCut,Cuts,Num,0._dp,0._dp)
  end if

end subroutine xSec_SIDIS

subroutine xSec_SIDIS_List(xx,process,s,pT,z,x,Q,doCut,Cuts,masses,doPartitioning)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                !Mandelshtam s
  real(dp),intent(in),dimension(:,:)::pT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(:,:)::z                !(zmin,zmax)
  real(dp),intent(in),dimension(:,:)::x                !(xmin,xmax)
  real(dp),intent(in),dimension(:,:)::Q                !(Qmin,Qmax)
  logical,intent(in),dimension(:)::doCut            !triger cuts
  real(dp),intent(in),dimension(:,:)::Cuts            !(ymin,yMax,W2)
  real(dp),intent(in),dimension(:,:),optional::masses        !(mass_target,mass-product)GeV
  logical, intent(in), optional :: doPartitioning
  real(dp),dimension(:),intent(out)::xx

  real(dp),dimension(1:size(s),1:2)::m2_list
  integer :: i,length

  logical::doP
  integer::k,numberOfP,listOfParts(1:size(s)),n_private
  integer,allocatable:: partI1(:),partSize(:)

if(.not.started) error stop ErrorString('The module is not initialized. Check INI-file.',moduleName)

  length=size(s)

  !!! checking sizes
  if(size(xx)/=length) then
    error stop ErrorString('xSec_SIDIS_List: sizes of xSec and s lists are not equal.',moduleName)
  end if
  if(size(process,1)/=length) then
    error stop ErrorString('xSec_SIDIS_List: sizes of process and s lists are not equal.',moduleName)
  end if
  if(size(pT,1)/=length) then
    error stop ErrorString('xSec_SIDIS_List: sizes of pT and s lists are not equal.',moduleName)
  end if
  if(size(x,1)/=length) then
    error stop ErrorString('xSec_SIDIS_List: sizes of x and s lists are not equal.',moduleName)
  end if
  if(size(Q,1)/=length) then
    error stop ErrorString('xSec_SIDIS_List: sizes of Q and s lists are not equal.',moduleName)
  end if
  if(size(z,1)/=length) then
    error stop ErrorString('xSec_SIDIS_List: sizes of z and s lists are not equal.',moduleName)
  end if
  if(size(doCut)/=length) then
    error stop ErrorString('xSec_SIDIS_List: sizes of doCut and s lists are not equal.',moduleName)
  end if
  if(size(Cuts,1)/=length) then
    error stop ErrorString('xSec_SIDIS_List: sizes of Cuts and s lists are not equal.',moduleName)
  end if
  if(size(process,2)/=4) then
    error stop ErrorString('xSec_SIDIS_List: process list must be (:,1:4).',moduleName)
  end if
  if(size(pT,2)/=2) then
    error stop ErrorString('xSec_SIDIS_List: pt list must be (:,1:2).',moduleName)
  end if
  if(size(x,2)/=2) then
    error stop ErrorString('xSec_SIDIS_List: x list must be (:,1:2).',moduleName)
  end if
  if(size(Q,2)/=2) then
    error stop ErrorString('xSec_SIDIS_List: Q list must be (:,1:2).',moduleName)
  end if
  if(size(z,2)/=2) then
    error stop ErrorString('xSec_SIDIS_List: z list must be (:,1:2).',moduleName)
  end if
  if(size(Cuts,2)/=4) then
    error stop ErrorString('xSec_SIDIS_List: cuts list must be (:,1:4).',moduleName)
  end if

  if(PRESENT(masses)) then
    if(size(masses,1)/=length) then
      error stop ErrorString('xSec_SIDIS_List: sizes of masses and s lists are not equal.',moduleName)
    end if
    if(size(masses,2)/=2) then
      error stop ErrorString('xSec_SIDIS_List: mass list must be (:,1:2).',moduleName)
    end if
    m2_list(1:length,1:2)=masses(1:length,1:2)
  else
    m2_list(1:length,1:2)=0._dp
  end if

  if(present(doPartitioning)) then
    doP=doPartitioning
  else
    doP=doPartitioning_byDefault
  end if

  if(doP) then
  !!!! attempt to make partitioning into pt-sectors
  !!!! basically, I run though the whole list and compare the terms if all (except pT) are the same, they are marked by the same number
  !!!! Only consecutive ranges are marked, and there is also an upper cut
    k=1
    listOfParts(1)=k
    do i=2,length
      if(&
      pT(i,2)*2>(z(i,1)+z(i,2))*MaxQT_range_toPartite &  !!! check the max size by formula pt/z=qT>maxPT
      .or.(abs(x(i,1)-x(i-1,1))>toleranceGEN) .or. (abs(x(i,2)-x(i-1,2))>toleranceGEN) &
      .or.(abs(z(i,1)-z(i-1,1))>toleranceGEN) .or. (abs(z(i,2)-z(i-1,2))>toleranceGEN) &
      .or.(abs(Q(i,1)-Q(i-1,1))>toleranceGEN) .or. (abs(Q(i,2)-Q(i-1,2))>toleranceGEN) &
      .or.(process(i,1)/=process(i-1,1)) .or. (process(i,2)/=process(i-1,2)) .or. (process(i,3)/=process(i-1,3)) &
      .or.(process(i,4)/=process(i-1,4)) .or. (abs(s(i-1)-s(i))>toleranceGEN) .or. (doCut(i-1).neqv.doCut(i)) &
      .or.(abs(Cuts(i,1)-Cuts(i-1,1))>toleranceGEN) &
      .or.(abs(Cuts(i,2)-Cuts(i-1,2))>toleranceGEN) &
      .or.(abs(Cuts(i,3)-Cuts(i-1,3))>toleranceGEN) &
      .or.(abs(Cuts(i,4)-Cuts(i-1,4))>toleranceGEN) &
      .or.(abs(m2_list(i,1)-m2_list(i-1,1))>toleranceGEN) .or. (abs(m2_list(i,2)-m2_list(i-1,2))>toleranceGEN) &
      .or.(pT(i-1,2)>pT(i,1)+toleranceGEN) & !!!!! the bins are succesive
      ) k=k+1

      listOfParts(i)=k
    end do
    numberOfP=k
    !!!! partI1 is a list that contains the initial index of each partition
    allocate(partI1(1:numberOfP))
    partI1(1)=1
    k=2
    do i=2,length
      if(listOfParts(i)/=listOfParts(i-1)) then
        partI1(k)=i
        k=k+1
      end if
    end do

    !!!! partSize is a list that contains the size of each partition
    allocate(partSize(1:numberOfP))

    do i=1,numberOfP-1
      partSize(i)=partI1(i+1)-partI1(i)
    end do
    partSize(numberOfP)=length-partI1(numberOfP)+1

   !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(n_private)
    do i=1,numberOfP
      !!!!! if bin is not partitioned, it is computed usually
      if(partSize(i)==1) then
        n_private=NumPT_auto(pT(partI1(i),2)-pT(partI1(i),1),(Q(partI1(i),2)+Q(partI1(i),1))/2d0)

        XX(partI1(i))=Xsec_PTint_Qint_Xint_Zint(process(partI1(i),1:4),s(partI1(i)),&
              z(partI1(i),1),z(partI1(i),2),x(partI1(i),1),x(partI1(i),2),&
              Q(partI1(i),1),Q(partI1(i),2),pT(partI1(i),1),pT(partI1(i),2),&
              doCut(partI1(i)),Cuts(partI1(i),1:4),n_private,&
              m2_list(partI1(i),1),m2_list(partI1(i),2))

      else
        !!! actual computation
        XX(partI1(i):partI1(i)+partSize(i)-1)=Xsec_pTspectrum_Qint_Xint_Zint(process(partI1(i),1:4),s(partI1(i)),&
              z(partI1(i),1),z(partI1(i),2),x(partI1(i),1),x(partI1(i),2),&
              Q(partI1(i),1),Q(partI1(i),2),pT(partI1(i):partI1(i)+partSize(i)-1,1),pT(partI1(i):partI1(i)+partSize(i)-1,2),&
              doCut(partI1(i)),Cuts(partI1(i),1:4),&
              m2_list(partI1(i),1),m2_list(partI1(i),2))

      end if
    end do
    !$OMP END PARALLEL DO
  else
  !!!!!--------------- compute in the usual way
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) DEFAULT(SHARED) PRIVATE(n_private)
    do i=1,length

      n_private=NumPT_auto(pT(i,2)-pT(i,1),(Q(i,2)+Q(i,1))/2d0)

      XX(i)=Xsec_PTint_Qint_Xint_Zint(process(i,1:4),s(i),&
          z(i,1),z(i,2),x(i,1),x(i,2),Q(i,1),Q(i,2),pT(i,1),pT(i,2),&
          doCut(i),Cuts(i,1:4),n_private,m2_list(i,1),m2_list(i,2))

    end do
    !$OMP END PARALLEL DO
  end if

end subroutine xSec_SIDIS_List

end module TMDX_SIDIS
