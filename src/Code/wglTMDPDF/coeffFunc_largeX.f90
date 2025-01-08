!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of wgtTMDPDF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the mathing coefficient--------------------------------------
!!!--------------The order is accumulative pertrubative order of coefficient =1,2,3 (LO,NLO,NNLO)---------
!!!-------------------------------------------------------------------------------------------------------

!!!! parametrizationString=(/1d0,z,Log(z),Log(1d0-z)/)

!!!!!coefficient function q<-q regular-part
pure function Coeff_q_q_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_q_q_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the Leading order is 1, it is WW-part of worm-gear function
    !! This 1 is accounded in the large-x convolution
    Coeff_q_q_reg_largeX=(/0d0,0d0,0d0,0d0/) !
    !!! WGL absorbes all finite parts into the resummation formula
end function Coeff_q_q_reg_largeX

!!!!!coefficient function q<-g regular-part  
pure function Coeff_q_g_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_q_g_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the gluon-mixing is always zero, for transversily polarized quarks
    Coeff_q_g_reg_largeX=(/0d0,0d0,0d0,0d0/)
end function Coeff_q_g_reg_largeX

!!!!!coefficient function g<-q regular-part  
pure function Coeff_g_q_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_g_q_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! there is no gluon, for transversily polarized quarks
    Coeff_g_q_reg_largeX=(/0d0,0d0,0d0,0d0/)
end function Coeff_g_q_reg_largeX

    !!!!!coefficient function g<-g regular-part  
function Coeff_g_g_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_g_g_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! there is no gluon, for transversily polarized quarks
    Coeff_g_g_reg_largeX=(/0d0,0d0,0d0,0d0/)

end function Coeff_g_g_reg_largeX

!!!!!coefficient function q<-qb regular-part  
pure function Coeff_q_qb_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_q_qb_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the Leading order is always zero, therefore calculation should be done only for order >=1
    Coeff_q_qb_reg_largeX=(/0d0,0d0,0d0,0d0/)!

end function Coeff_q_qb_reg_largeX

    !!!!!coefficient function q<-qp regular-part  
pure function Coeff_q_qp_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_q_qp_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the Leading order is always zero, therefore calculation should be done only for order >=1
    Coeff_q_qp_reg_largeX=(/0d0,0d0,0d0,0d0/)
end function Coeff_q_qp_reg_largeX
