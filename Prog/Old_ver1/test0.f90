program example
use aTMDe_setup
implicit none
  
 call artemide_Setup_Default('LO')
 
 !call artemide_include('TMDX_SIDIS','uTMDFF','TMDs_inKT')
 
 call Set_uTMDPDF(1,'MHHT14nlo68cl')
 
 call Set_uTMDPDF_order('LO')
 call Set_uTMDPDF_gridEvaluation(.true.,.false.)
 call Set_uTMDPDF_lengthNParray(7)
 
!  call Set_uTMDFF(1,'dsspipNLO')
!  call Set_uTMDFF(2,'dsspimNLO')
!  call Set_uTMDFF(3,'dssKpNLO')
!  call Set_uTMDFF(4,'dssKmNLO')
 
!  call Set_uTMDFF_order('LO')
!  call Set_uTMDFF_gridEvaluation(.true.,.false.)
!  call Set_uTMDFF_lengthNParray(4)
 
 call CreateConstantsFile('const-test')


! !------------------------------------------------
!  call artemide_Setup_Default('LO')
!  
!  call Set_uTMDPDF(1,'NNPDF31_nnlo_as_0118')
!  
!  call Set_uTMDPDF_gridEvaluation(.true.,.false.)
!  call Set_uTMDPDF_lengthNParray(7)
!  
!  call Set_TMDR_lengthNParray(2)
!  
!  call CreateConstantsFile('const-DYfit18_LO',prefix='/misc/data2/braun/vla18041/arTeMiDe_Repository/artemide/')
 
end program example