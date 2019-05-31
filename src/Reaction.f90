MODULE Reaction_MOD
  USE Type_MOD
  USE Config_MOD
  USE Util_MOD
  USE DEF_MOD
  USE NEUTRON_MOD
  IMPLICIT NONE
  

  
  CONTAINS
  
  SUBROUTINE GetDTC(Neutron)
  IMPLICIT NONE
  INTEGER :: ir,ig,imix,icelltype
  REAL :: Sigtr, DTC, xi
  TYPE(Neutron_TYPE) :: Neutron
  
  CALL RANDOM_NUMBER(xi)
  
  ig = Neutron%g
  
  icelltype = Neutron%mycell%icelltype
  ir = Neutron%ring
  
  imix = CellInfo(icelltype)%mat(ir)
  
  Sigtr = XS(imix)%Sigtr(ig)
  Neutron%DTC = -log(xi)/Sigtr
  
  END SUBROUTINE
  
  SUBROUTINE SelectReactionType(Neutron, reactiontype)
  IMPLICIT NONE
  INTEGER :: ReactionType, i, imix, g
  REAL :: xi
  TYPE(Neutron_TYPE) :: Neutron

  g = Neutron%g
  
  CALL RANDOM_NUMBER(xi)
  
  imix = CellInfo(Neutron%mycell%icelltype)%mat(Neutron%ring)
  DO i = 1, 3
    IF (xi .LT. CDF(imix)%reactiontype(i,g)) EXIT
  END DO
  reactiontype = i
  END SUBROUTINE
  
  SUBROUTINE FissionReaction(Neutron)
  IMPLICIT NONE
  INTEGER :: i, nu 
  TYPE(Neutron_TYPE) :: Neutron

  CALL Samplenu(Neutron,nu)
  
  DO i = 1, nu
    ! direction
    CALL SampleAngle(Neutron)
    ! Energy
    CALL SampleGroup_Fis(Neutron)
    CALL Bank_Fis%Push(Neutron)
  END DO
  END SUBROUTINE
  
  SUBROUTINE ScatteringReaction(Neutron)
  IMPLICIT NONE
  INTEGER :: i
  TYPE(Neutron_TYPE) :: Neutron
  
  CALL SampleAngle(Neutron)
  CALL SampleGroup_scat(Neutron)

  END SUBROUTINE
  
  SUBROUTINE SampleGroup_Fis(Neutron)
  IMPLICIT NONE
  TYPE(Neutron_TYPE) :: Neutron
  REAL :: xi
  INTEGER :: ig, imix
  
  imix = CellInfo(Neutron%mycell%icelltype)%mat(Neutron%ring)  
  CALL RANDOM_NUMBER(xi)
  DO ig = 1, Ng
    IF (xi .LT. CDF(imix)%chi(ig)) EXIT
  END DO
  
  Neutron%g = ig
  
  END SUBROUTINE
  
  SUBROUTINE SampleGroup_Scat(Neutron)
  IMPLICIT NONE
  INTEGER :: jg, ig, imix
  REAL :: xi
  TYPE(Neutron_TYPE) :: Neutron
  imix = CellInfo(Neutron%mycell%icelltype)%mat(Neutron%ring)  
  ig = Neutron%g
  
  CALL RANDOM_NUMBER(xi)
  
  DO jg = 1, Ng
    IF (xi .LT. CDF(imix)%scat(ig,jg)) EXIT
  END DO
  Neutron%g = jg
  END SUBROUTINE
  
  SUBROUTINE Samplenu(neutron,nu)
  IMPLICIT NONE
  REAL :: yield
  INTEGER :: nu, g
  INTEGER :: imix
  REAL :: w, xi
  TYPE(Neutron_TYPE) :: Neutron
  
  imix = CellInfo(Neutron%mycell%icelltype)%mat(Neutron%ring) 
  g = neutron%g
  w = neutron%w
  yield = XS(imix)%nu(g)
  CALL RANDOM_NUMBER(xi)
  nu = INT(w*yield/keff+xi)
  
  END SUBROUTINE
  
  SUBROUTINE SampleAngle(neutron)
  IMPLICIT NONE
  REAL :: theta, xi
  TYPE(Neutron_TYPE) :: Neutron
  CALL RANDOM_NUMBER(xi)
  theta = 2*pi*xi
  Neutron%theta = theta
  neutron%azic = cos(theta); Neutron%azis = sin(theta)
  
!  CALL RANDOM_NUMBER(xi)
!  Neutron%polars = - 1 + 2*xi
!  neutron%polarc = sqrt(1. - Neutron%polars**2)
  
  END SUBROUTINE
   
  SUBROUTINE SetCDF
  IMPLICIT NONE
  INTEGER :: ig, jg, i, imix
  REAL :: sum_reaction
  REAL, ALLOCATABLE :: scatsum_ig(:), chisum_ig(:)
  ALLOCATE(scatsum_ig(Nmixture), chisum_ig(Nmixture))
  DO imix = 1,Nmixture 
    DO ig = 1, Ng
      scatsum_ig(imix) = SUM(XS(imix)%Sigs(ig,:))
      sum_reaction = XS(imix)%Sigr(ig) + XS(imix)%Sigf(ig) + scatsum_ig(imix)
      CDF(imix)%reactiontype(cap,ig) = XS(imix)%Sigr(ig)/sum_reaction
      CDF(imix)%reactiontype(fis,ig) = CDF(imix)%reactiontype(cap,ig) + XS(imix)%Sigf(ig)/sum_reaction
      CDF(imix)%reactiontype(scat,ig) = 1
    END DO
  END DO
  
  DO imix = 1, Nmixture 
    chisum_ig(imix) = SUM(XS(imix)%chi)
    DO ig = 1, Ng
      CDF(imix)%chi(ig) = SUM(XS(imix)%chi(1:ig))/chisum_ig(imix)
      scatsum_ig(imix) = SUM(XS(imix)%Sigs(ig,:))
      DO jg = 1, Ng
        CDF(imix)%scat(ig,jg) = SUM(XS(imix)%Sigs(ig,1:jg))/scatsum_ig(imix)
      END DO        
    END DO
  END DO
  END SUBROUTINE
  
  SUBROUTINE UpdtWeight(Neutron)
  IMPLICIT NONE
  INTEGER :: imix, ig
  REAL :: Sigtr, Siga, Sfac
  TYPE(Neutron_TYPE) :: Neutron
  
  imix = CellInfo(Neutron%mycell%icelltype)%mat(Neutron%ring)
  ig = Neutron%g
  
  Sigtr = XS(imix)%Sigtr(ig)
  Siga = XS(imix)%Siga(ig)
  
  Sfac = 1 - Siga/Sigtr
  Neutron%w = Neutron%w * Sfac
  END SUBROUTINE
  
  SUBROUTINE CountOutgoing(Neutron,Status,incoming_idx)
  IMPLICIT NONE
  TYPE(Neutron_TYPE) :: Neutron
  TYPE(Status_TYPE) :: Status
  
  INTEGER :: ig_out, iazi_out, dir, incoming_idx
  REAL :: theta, theta_out, normvec(3), x, y, mu, hpitch, qpitch
  
  theta = Neutron%theta
  
  normvec = Status%nout
  x = normvec(1); y = normvec(2)
  
  qpitch = 0.25*PITCH
  
  theta_out = theta + (pi/2.*x - INT(2./(3.*pi)*theta)*2*pi*abs(x)) + (pi/2.*y - pi/2.*abs(y))
  iazi_out = INT(16. * theta_out / pi) + 1.
  ig_out = Neutron%g
  mu = x*Neutron%azic + y*Neutron%azis
  dir = 2.*x + 3.*abs(x) - 2.*y + 11.*abs(y) + INT((Neutron%x*abs(y) + Neutron%y*abs(x))/qpitch)
  Outgoing(ig_out, iazi_out, incoming_idx, dir)%count = Outgoing(ig_out, iazi_out, incoming_idx, dir)%count + Neutron%w
  END SUBROUTINE
  END MODULE
  
  