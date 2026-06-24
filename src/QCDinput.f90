!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.04
!
! Interface module to the user defined alpha-s, PDF, FF, etc.
!
!    ver.2.0: 28.03.2019 AV
!    ver.2.01: 14.06.2019 AV (+lpPDF)
!    ver.2.07: 09.11.2021 AV (+g1T)
!    ver.3.01: 22.08.2024 AV (introduction of LHA)
!    ver.3.03: 16.10.2025 AV (LHA_PDF replaced by a class)
!    ver.3.04: 24.06.2026 AV (cleaning of the code)
!                A.Vladimirov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module QCDinput
use aTMDe_Numerics
use aTMDe_IO
use LHA_alpha, only : ReadInfo_alpha => ReadInfo, AlphaS_fromLHA => AlphaS
use LHA_PDF
!!
implicit none

private

public::QCDinput_Initialize,QCDinput_IsInitialized
!!!! IMPORTANT: As retuns \alpha_s(Q)/(4pi)
public::As,activeNf
!!!! Interface of all PDFs is (x,Q,h) where
!!!! x = Bjorken x; Q = scale;
!!!! h = is the "hadron" number, i.e. the id number for this PDF, it actually defines which hadron is used
public::xPDF,xFF,x_lp_PDF,x_gPDF,x_hPDF
public:: QCDinput_SetPDFreplica, QCDinput_SetFFreplica, QCDinput_SetlpPDFreplica,QCDinput_SetgPDFreplica,QCDinput_SethPDFreplica

character (len=8),parameter :: moduleName="QCDinput"
!Current version of module
character (len=5),parameter :: version="v3.04"
!Last appropriate version of constants-file
integer,parameter::inputver=35
!--- general
logical:: started=.false.
integer::outputLevel

real(dp),public::mCHARM,mBOTTOM,mTOP

!---uPDFs
integer::num_of_uPDFs
integer,allocatable::current_replica_uPDFs(:)
type(LHAPDFgridReader),allocatable::uLHAPDF(:)

!---uFFs
integer::num_of_uFFs
integer,allocatable::current_replica_uFFs(:)
type(LHAPDFgridReader),allocatable::uLHAFF(:)

!---lpPDFs
integer::num_of_lpPDFs
integer,allocatable::current_replica_lpPDFs(:)
type(LHAPDFgridReader),allocatable::lpLHAPDF(:)

!---gPDFs
integer::num_of_gPDFs
integer,allocatable::current_replica_gPDFs(:)
type(LHAPDFgridReader),allocatable::gLHAPDF(:)

!---hPDFs
integer::num_of_hPDFs
integer,allocatable::current_replica_hPDFs(:)
type(LHAPDFgridReader),allocatable::hLHAPDF(:)

contains
 
!!!!! Function indicates if the module is initialized or not.
!!!!! Used by other higher modules to check -- if the initialization is needed
function QCDinput_IsInitialized()
logical::QCDinput_IsInitialized
QCDinput_IsInitialized=started
end function QCDinput_IsInitialized


!!!! initialization of the module.
!!!! It also initialize all LHA_alpha and LHA_PDF modules and objects, and set all replicas to 0.
subroutine QCDinput_Initialize(file,prefix)
character(len=*),intent(in)::file
character(len=*),optional,intent(in)::prefix
character(len=:),allocatable::path
character(len=600)::lineToRead

character(len=:),allocatable::pathToLHA
character(len=:),allocatable::alphaNAME
character(len=64),allocatable::names_uPDF(:),names_uFF(:),names_lpPDF(:),names_gPDF(:),names_hPDF(:)
integer::i,FILEver,j

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
    write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
    write(*,*) '             Update the const-file with artemide.setup'
    ERROR STOP '  '
  end if
  call MoveTO(51,'*p2  ')
  read(51,*) outputLevel
  if(outputLevel>1) write(*,*) '--------------------------------------------- '
  if(outputLevel>1) write(*,*) 'artemide.QCDinput: initialization started ... '

  !!! Search for QCDinput initialization options
  call MoveTO(51,'*1   ')

  !!!! reading the path to LHAPDF
  call MoveTO(51,'*p1  ')
  read(51,"(A)") lineToRead
  pathToLHA=trim(adjustl(lineToRead))//"/"
  call MoveTO(51,'*p2  ')
  read(51,"(A)") lineToRead
  alphaNAME=trim(adjustl(lineToRead))

  !!! Search for parameters
  call MoveTO(51,'*A   ')
  call MoveTO(51,'*p1  ')
  read(51,*) mCHARM
  call MoveTO(51,'*p2  ')
  read(51,*) mBOTTOM
  call MoveTO(51,'*p3  ')
  read(51,*) mTOP

  if(outputLevel>2) write(*,*) '    mass of charm : ',mCHARM
  if(outputLevel>2) write(*,*) '    mass of bottom: ',mBOTTOM
  if(outputLevel>2) write(*,*) '    mass of top   : ',mTOP


  !!!---------------------------------------- Search for uPDF initialization options
  call MoveTO(51,'*B   ')
  call MoveTO(51,'*p1  ')
  read(51,*) num_of_uPDFs

  if(num_of_uPDFs>0) then

    allocate(names_uPDF(1:num_of_uPDFs))
    call MoveTO(51,'*p2  ')
    do i=1,num_of_uPDFs
      read(51,"(A)") lineToRead
      names_uPDF(i)=trim(adjustl(lineToRead))
    end do

  else
    !!! initialization is not needed
    if(outputLevel>2) write(*,*)'    no uPDFs to initialize...'
  end if

  !!! -------------------------------------Search for uFF initialization options
  call MoveTO(51,'*C   ')
  call MoveTO(51,'*p1  ')
  read(51,*) num_of_uFFs

  if(num_of_uFFs>0) then

    allocate(names_uFF(1:num_of_uFFs))
    call MoveTO(51,'*p2  ')
    do i=1,num_of_uFFs
      read(51,"(A)") lineToRead
      names_uFF(i)=trim(adjustl(lineToRead))
    end do

  else
    !!! initialization is not needed
    if(outputLevel>2) write(*,*)'    no uFFs to initialize...'
  end if

    !!!---------------------------------------- Search for lpPDF initialization options
  call MoveTO(51,'*D   ')
  call MoveTO(51,'*p1  ')
  read(51,*) num_of_lpPDFs
  if(num_of_lpPDFs>0) then

    allocate(names_lpPDF(1:num_of_lpPDFs))
    call MoveTO(51,'*p2  ')
    do i=1,num_of_lpPDFs
      read(51,"(A)") lineToRead
      names_lpPDF(i)=trim(adjustl(lineToRead))
    end do

  else
    !!! initialization is not needed
    if(outputLevel>2) write(*,*)'    no lpPDFs to initialize...'
  end if


  !!!---------------------------------------- Search for gPDF initialization options
  call MoveTO(51,'*E   ')
  call MoveTO(51,'*p1  ')
  read(51,*) num_of_gPDFs
  if(num_of_gPDFs>0) then
    allocate(names_gPDF(1:num_of_gPDFs))
    call MoveTO(51,'*p2  ')
    do i=1,num_of_gPDFs
      read(51,"(A)") lineToRead
      names_gPDF(i)=trim(adjustl(lineToRead))
    end do

  else
    !!! initialization is not needed
    if(outputLevel>2) write(*,*)'    no gPDFs to initialize...'
  end if

  !!!---------------------------------------- Search for hPDF initialization options
  call MoveTO(51,'*F   ')
  call MoveTO(51,'*p1  ')
  read(51,*) num_of_hPDFs
  if(num_of_hPDFs>0) then
    allocate(names_hPDF(1:num_of_hPDFs))
    call MoveTO(51,'*p2  ')
    do i=1,num_of_hPDFs
      read(51,"(A)") lineToRead
      names_hPDF(i)=trim(adjustl(lineToRead))
    end do

  else
    !!! initialization is not needed
    if(outputLevel>2) write(*,*)'    no hPDFs to initialize...'
  end if

CLOSE (51, STATUS='KEEP')

!   !-----------------------------------------------------
!   !=====LHAPDF======
!   call InitPDFsetByNameM(1,alphaNAME)
!   call InitPDFM(1,0)
!   !=================
!   !-----------------------------------------------------

!!!! initialize the alphaS-LHA table
call ReadInfo_alpha(alphaNAME,pathToLHA,outputLevel)

!!!! initialization of uPDFs (by default they are set to replica 0)
allocate(current_replica_uPDFs(1:num_of_uPDFs))
allocate(uLHAPDF(1:num_of_uPDFs))
do j=1,num_of_uPDFs
  uLHAPDF(j)=LHAPDFgridReader(trim(names_uPDF(j)),pathToLHA,outputLevel)
  current_replica_uPDFs(j)=0
end do

!!!! initialization of uFFs (by default they are set to replica 0)
allocate(current_replica_uFFs(1:num_of_uFFs))
allocate(uLHAFF(1:num_of_uFFs))
do j=1,num_of_uFFs
  uLHAFF(j)=LHAPDFgridReader(trim(names_uFF(j)),pathToLHA,outputLevel)
  current_replica_uFFs(j)=0
end do

!!!! initialization of lpPDFs (by default they are set to replica 0)
allocate(current_replica_lpPDFs(1:num_of_lpPDFs))
allocate(lpLHAPDF(1:num_of_lpPDFs))
do j=1,num_of_lpPDFs
  lpLHAPDF(j)=LHAPDFgridReader(trim(names_lpPDF(j)),pathToLHA,outputLevel)
  current_replica_lpPDFs(j)=0
end do

!!!! initialization of gPDFs (by default they are set to replica 0)
allocate(current_replica_gPDFs(1:num_of_gPDFs))
allocate(gLHAPDF(1:num_of_gPDFs))
do j=1,num_of_gPDFs
  gLHAPDF(j)=LHAPDFgridReader(trim(names_gPDF(j)),pathToLHA,outputLevel)
  current_replica_gPDFs(j)=0
end do

!!!! initialization of hPDFs (by default they are set to replica 0)
allocate(current_replica_hPDFs(1:num_of_hPDFs))
allocate(hLHAPDF(1:num_of_hPDFs))
do j=1,num_of_hPDFs
  hLHAPDF(j)=LHAPDFgridReader(trim(names_hPDF(j)),pathToLHA,outputLevel)
  current_replica_hPDFs(j)=0
end do

  if(outputLevel>0) write(*,*) color('----- arTeMiDe.QCDinput '//trim(version)//': .... initialized',c_green)
  if(outputLevel>1) write(*,*) ' '
started=.true.
end subroutine QCDinput_Initialize

!!! set a different replica of uPDF (pointing to hadron)
!!! check whether it is the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SetPDFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_uPDFs) &
      ERROR STOP ErrorString('SetPDFreplica. Called with non-existent hadron. h= '//trim(numToStr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_uPDFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_uPDFs(hadron)=rep
      call uLHAPDF(hadron)%SetReplica(rep)

      newPDF=.true.
    end if
end subroutine QCDinput_SetPDFreplica

!!! set a different replica of uFF (pointing to hadron)
!!! check whether it is the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SetFFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_uFFs) &
      ERROR STOP ErrorString('SetFFreplica. Called with non-existent hadron. h= '//trim(numToStr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_uFFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_uFFs(hadron)=rep
      call uLHAFF(hadron)%SetReplica(rep)

      newPDF=.true.
    end if
end subroutine QCDinput_SetFFreplica

!!! set a different replica of lp_PDF (pointing to hadron)
!!! check whether it is the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SetlpPDFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_lpPDFs) &
      ERROR STOP ErrorString('SetlpPDFreplica. Called with non-existent hadron. h= '//trim(numToStr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_lpPDFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_lpPDFs(hadron)=rep
      call lpLHAPDF(hadron)%SetReplica(rep)

      newPDF=.true.
    end if
end subroutine QCDinput_SetlpPDFreplica

!!! set a different replica of gPDF (pointing to hadron)
!!! check whether it is the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SetgPDFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_gPDFs) &
      ERROR STOP ErrorString('SetgPDFreplica. Called with non-existent hadron. h= '//trim(numToStr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_gPDFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_gPDFs(hadron)=rep
      call gLHAPDF(hadron)%SetReplica(rep)

      newPDF=.true.
    end if
end subroutine QCDinput_SetgPDFreplica

!!! set a different replica of hPDF (pointing to hadron)
!!! check whether it is the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SethPDFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_hPDFs) &
      ERROR STOP ErrorString('SethPDFreplica. Called with non-existent hadron. h= '//trim(numToStr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_hPDFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_hPDFs(hadron)=rep
      call hLHAPDF(hadron)%SetReplica(rep)
      newPDF=.true.
    end if
end subroutine QCDinput_SethPDFreplica
 
!!number of active flavors at a given mu
elemental function activeNf(mu)
  real(dp),intent(in)::mu
  integer::activeNf

  if(mu>mBOTTOM) then
     activeNf=5
    else if(mu>mCHARM) then
     activeNf=4
    else
     activeNf=3
  end if
end function activeNf
 
!!!!alphas(Q)/4pi
!!! DO NOT FORGET 4pi !!!
pure function As(Q)
real(dp),intent(in)::Q
real(dp)::As

As=AlphaS_fromLHA(Q)/pix4
!  as=1d0/(2d0*23d0/3d0*Log(Q/0.08782708014552364d0))  !!! Nf=5 LO solution


! !-------------------------
! !======LHAPDF========
! real(dp)::alphasPDF
! As=alphasPDF(Q)/pix4
! !======LHAPDF========
! !-------------------------


end function As

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!! unpolarized PDF used in uTMDPDF
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
pure function xPDF(x,Q,hadron)
    real(dp),intent(in) :: x,Q
    integer,intent(in):: hadron
    real(dp), dimension(-5:5):: xPDF

!     !-------------------------
!     !======LHAPDF========
!     real(dp), dimension(-6:6):: inputPDF
!     call evolvePDFM(1,x,Q,inputPDF)
!     xPDF=inputPDF(-5:5)
!     !======LHAPDF========
!     !-------------------------
    xPDF=uLHAPDF(hadron)%xPDF(x,Q)

end function xPDF
  
!!!! return x*F(x,Q)
!!!! enumeration of flavors
!!!!  f = -5,-4, -3,  -2,  -1,0,1,2,3,4,5
!!!!    = bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b
!!! unpolarized FF used in uTMDFF
!!!! enumeration of hadrons
pure function xFF(x,Q,hadron)
    integer,intent(in) :: hadron
    real(dp),intent(in) :: x,Q
    real(dp),dimension(-5:5):: xFF

    xFF=uLHAFF(hadron)%xPDF(x,Q)

end function xFF
 
!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!! unpolarized PDF used in lpTMDPDF
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
pure function x_lp_PDF(x,Q,hadron)
    real(dp),intent(in) :: x,Q
    integer,intent(in):: hadron
    real(dp), dimension(-5:5):: x_lp_PDF

    x_lp_PDF=lpLHAPDF(hadron)%xPDF(x,Q)

end function x_lp_PDF

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!! helicity PDF
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
pure function x_gPDF(x,Q,hadron)
    real(dp),intent(in) :: x,Q
    integer,intent(in):: hadron
    real(dp), dimension(-5:5):: x_gPDF

    x_gPDF=gLHAPDF(hadron)%xPDF(x,Q)

end function x_gPDF

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!! transversity PDF
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
pure function x_hPDF(x,Q,hadron)
    real(dp),intent(in) :: x,Q
    integer,intent(in):: hadron
    real(dp), dimension(-5:5):: x_hPDF

    x_hPDF=hLHAPDF(hadron)%xPDF(x,Q)

end function x_hPDF
 
end module QCDinput
