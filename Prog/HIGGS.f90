!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at given s,Q,x,z, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_DY
implicit none

integer,dimension(1:3)::proc
real*8::x,z,Q,s,y
real*8::stepqT,maxqT
integer::i,maxI
real*8,allocatable::Xsec(:),qT(:),XsecQ(:),Qvalues(:)
real*8,dimension(1:2)::pTvec,zvec,xvec,Qvec
real*8::Xs,ptMin,ptMax,zMin,zMax,xMin,xMax,Qmin,Qmax


write(*,*) "Initialize artemide at LO. It is fast ~1 sek. NNLO could be much longer ~5 min."
call artemide_Initialize('Prog/const-file_higgs')
!call artemide_Initialize('const-file_higgs',prefix='/Users/dd46mi/Documents/work/artemide-public-v2p03/constant-files/')

call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

proc=(/1,5,22/) !!!! pp -> Higgs (unpol.part)
s=7000d0**2
Q=125d0
y=0d0

call TMDX_DY_setProcess(proc)
call TMDX_DY_XSetup(s,Q,y)
call TMDX_DY_SetCuts(.false.,0d0,-10d0,10d0)!!! no cuts

maxqT=0.4d0*Q
stepqT=maxqT/25d0
maxI=Int(maxqT/stepqT)

allocate(qT(1:maxI))
allocate(Xsec(1:maxI))

do i=1,maxI
  qT(i)=stepqT*i
end do

call CalcXsec_DY(Xsec,qT)
write(*,*) 'Hello', Xsec(1)

do i=1,maxI
  write(*,'("{",F6.3,",",F12.8,"},")') qT(i),Xsec(i)/1000d0
end do



!zMin=z
!zMax=zMin+0.001d0
!xMin=x
!xMax=xMin+0.001d0
!ptMin=0.64d0
!ptMax=ptMin+0.001d0
!Qmin=Q
!Qmax=Qmin+0.001d0
!!call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(XsecQ,ptMin,ptMax,zMin,zMax,xMin,xMax,Qmin,Qmax)
!call CalcXsec_SIDIS_PTint_Zint_Xint_Qint(Xs,ptMin,ptMax,zMin,zMax,xMin,xMax,Qmin,Qmax)
!write(*,*) Qmin,Xs

end program example
