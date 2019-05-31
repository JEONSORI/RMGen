MODULE SOLVER_MOD
  USE Type_MOD
  USE Config_MOD
  USE Util_MOD
  USE DEF_MOD
  USE NEUTRON_MOD
  USE REACTION_MOD
  IMPLICIT NONE
  
  CONTAINS
  
  SUBROUTINE Alloc
  IMPLICIT NONE
  
  IF (problem .EQ. 1) THEN
    ALLOCATE(Incoming(Ng,Nazi,4), Outgoing(Ng,Nazi,4*Ng*Nazi,16), phi(Ng,4,4*Ng*Nazi))
    ALLOCATE(Response(Ng*Nazi*16,Ng*Nazi*4)); Response = 0.
  ELSE IF (problem .EQ. 2) THEN
    ALLOCATE(Incoming(Ng,Nazi,4), Outgoing(Ng,Nazi,4*Ng*Nazi,16), phi(Ng,4,4*Ng*Nazi))
    ALLOCATE(Response(Nfuel*Ng,Ng*Nazi*4)); Response = 0.
  ELSE IF (problem .EQ. 3) THEN
    ALLOCATE(Incoming(Ng,1,Nfuel), Outgoing(Ng,Nazi,Nfuel*Ng,16), phi(Ng,4,Nfuel*Ng))
    ALLOCATE(Response(Ng*Nazi*16,Nfuel*Ng), Response2(3*Ng,3*Ng)); Response = 0.; Response2 = 0.
  END IF
  
  
  END SUBROUTINE
  
  SUBROUTINE Solver
  IMPLICIT NONE
  TYPE(Neutron_TYPE), POINTER :: Neutron
  TYPE(Status_TYPE) :: Status
  INTEGER :: i, reactiontype, icyc, ig, imix, iazi, res1 = 14, res2 = 15, dir = 0, idir = 0, irvol, irphi
  INTEGER :: incoming_idx, outgoing_idx, phi_idx, outgoing_idx0, phi_idx0
  REAL :: k_tl, w, length, nusigf, sumkeff = 0., sumkeffsqr = 0., mean = 0., std = 0., mean_keffsqr = 0., meansqr = 0., mean_temp, mu
  REAL :: hpitch, qpitch
  REAL :: std_temp = 0., mean_keffsqr_temp = 0., meansqr_temp = 0.
  TYPE(Flux_TYPE), POINTER :: phi0(:,:,:), incoming0(:,:,:), outgoing0(:,:,:,:)
  REAL, POINTER :: response0(:,:)
  phi0 => phi; incoming0 => incoming; outgoing0 => outgoing; response0 => response
  
  CALL SetCDF
!  Bank_Head%next => Bank_Fis%next
!  Bank_Head%index = Bank_Fis%index
!  NULLIFY(Bank_Fis%next); Bank_Fis%index = 0
  hpitch = 0.5*PITCH; qpitch = 0.25*pitch
  
  IF (problem .EQ. 1) THEN ! Surface to Surface
    Incoming(:,:,:)%count = 0.; outgoing(:,:,:,:)%count = 0.; phi(:,:,:)%count = 0.
    DO i = 1, Nparticle
      Neutron => NeutronType(i)
      iazi = Neutron%iazi; ig = Neutron%g; mu = Neutron%azis
      idir = INT(Neutron%x/qpitch) + 1
      Incoming(ig,iazi,idir)%count = Incoming(ig,iazi,idir)%count + 1.
      incoming_idx = (idir - 1)*Ng*Nazi + (iazi - 1)*Ng + ig
      Neutron%kill = .FALSE.
      DO WHILE(.NOT. Neutron%kill)
        Status = GetIntersection(Neutron)  ! --- Get DTS
      !  Status = Neutron%Cast()   
        CALL GetDTC(Neutron)
        IF (Status%distance .LT. Neutron%DTC) THEN ! DTS < DTC
          CALL MovePosition(Neutron, Status)
      !    CALL Neutron%Move(Status)
          IF (status%boundary) CALL CountOutgoing(Neutron,Status,incoming_idx)
          IF (Neutron%kill) EXIT
        ELSE  ! DTS > DTC
          Neutron%x = Neutron%x + Neutron%DTC * Neutron%azic  ! relocation  
          Neutron%y = Neutron%y + Neutron%DTC * Neutron%azis  ! relocation
          CALL UpdtWeight(Neutron)
          CALL ScatteringReaction(Neutron)
        END IF
      END DO
    END DO
    
    DO idir = 1, 4
      DO iazi = 1, Nazi
        DO ig = 1, Ng
          incoming_idx = (idir - 1)*Ng*Nazi + (iazi - 1) * Ng + ig
          Outgoing(:,:,incoming_idx,:)%count = Outgoing(:,:,incoming_idx,:)%count/Incoming(ig,iazi,idir)%count
        END DO
      END DO
    END DO
    
    DO dir = 1, 16
      DO iazi = 1, Nazi
        DO ig = 1, Ng
            Outgoing_idx = (dir - 1) * Ng * Nazi + (iazi - 1) * Ng + ig
            Response(outgoing_idx,:) = Response(outgoing_idx,:) + Outgoing(ig,iazi,:,dir)%count
        END DO
      END DO
    END DO

    
  ELSE IF (problem .EQ. 2) THEN ! Surface to Flux
    Incoming(:,:,:)%count = 0.; outgoing(:,:,:,:)%count = 0.; phi(:,:,:)%count = 0.
    DO i = 1, Nparticle
      Neutron => NeutronType(i)
      iazi = Neutron%iazi; ig = Neutron%g; mu = Neutron%azis
      idir = INT(Neutron%x/qpitch) + 1
      Incoming(ig,iazi,idir)%count = Incoming(ig,iazi,idir)%count + 1.
      incoming_idx = (idir - 1)*Ng*Nazi + (iazi - 1)*Ng + ig
      Neutron%kill = .FALSE.
      DO WHILE(.NOT. Neutron%kill)
        Status = GetIntersection(Neutron)  ! --- Get DTS
        CALL GetDTC(Neutron)
        length = min(Status%distance, Neutron%DTC)
        phi(Neutron%g,Neutron%ring,incoming_idx)%count = phi(Neutron%g,Neutron%ring,incoming_idx)%count + Neutron%w * length
        IF (Status%distance .LT. Neutron%DTC) THEN ! DTS < DTC
          CALL MovePosition(Neutron, Status)
          IF (Neutron%kill) EXIT
        ELSE  ! DTS > DTC
          Neutron%x = Neutron%x + Neutron%DTC * Neutron%azic  ! relocation  
          Neutron%y = Neutron%y + Neutron%DTC * Neutron%azis  ! relocation
          CALL UpdtWeight(Neutron)
          CALL ScatteringReaction(Neutron)
        END IF
      END DO
    END DO
    
    DO idir = 1, 4
      DO iazi = 1, Nazi
        DO ig = 1, Ng
          incoming_idx = (idir - 1)*Ng*Nazi + (iazi - 1) * Ng + ig
          phi(:,:,incoming_idx)%count = phi(:,:,incoming_idx)%count/Incoming(ig,iazi,idir)%count
        END DO
      END DO
    END DO
    
    DO irphi = 1, 3
      DO ig = 1, Ng
        phi_idx = (irphi - 1)*Ng + ig
        Response(phi_idx,:) = Response(phi_idx, :) + phi(ig,irphi,:)%count
      END DO
    END DO
    

  ELSE IF (problem .EQ. 3) THEN ! Volume to Surface
    Incoming(:,:,:)%count = 0.; outgoing(:,:,:,:)%count = 0.; phi(:,:,:)%count = 0.
    DO i = 1, Nparticle
      Neutron => NeutronType(i)
      ig = Neutron%g
      irvol = Neutron%ring
      Incoming(ig,1,irvol)%count = Incoming(ig,1,irvol)%count + 1
      incoming_idx = (irvol - 1)*Ng + ig
      Neutron%kill = .FALSE.
      DO WHILE(.NOT. Neutron%kill)
        Status = GetIntersection(Neutron)  ! --- Get DTS
        CALL GetDTC(Neutron)
        length = min(Status%distance, Neutron%DTC)
        phi(Neutron%g,Neutron%ring,incoming_idx)%count = phi(Neutron%g,Neutron%ring,incoming_idx)%count + Neutron%w * length
        IF (Status%distance .LT. Neutron%DTC) THEN ! DTS < DTC
          CALL MovePosition(Neutron, Status)
          IF (status%boundary) CALL CountOutgoing(Neutron,Status,incoming_idx)
          IF (Neutron%kill) EXIT
        ELSE  ! DTS > DTC
          Neutron%x = Neutron%x + Neutron%DTC * Neutron%azic  ! relocation  
          Neutron%y = Neutron%y + Neutron%DTC * Neutron%azis  ! relocation
          CALL UpdtWeight(Neutron)
          CALL ScatteringReaction(Neutron)
        END IF
      END DO
    END DO

    DO irvol = 1, 3
      DO ig = 1, Ng
        incoming_idx = (irvol - 1)*Ng + ig
        Outgoing(:,:,incoming_idx,:)%count = Outgoing(:,:,incoming_idx,:)%count / Incoming(ig,1,irvol)%count
        phi(:,:,incoming_idx)%count = phi(:,:,incoming_idx)%count / Incoming(ig,1,irvol)%count
      END DO
    END DO
    
    DO dir = 1, 16 
      DO iazi = 1, Nazi
        DO ig = 1, Ng
          outgoing_idx = (dir - 1) * Nazi * Ng + (iazi - 1) * Ng + ig
          Response(outgoing_idx, :) = Response(outgoing_idx, :) + Outgoing(ig,iazi,:,dir)%count
        END DO
      END DO
    END DO
    
    DO irphi = 1, 3
      DO ig = 1, Ng
        phi_idx = (irphi - 1)*Ng + ig
        Response2(phi_idx, :) = Response2(phi_idx, :) + phi(ig,irphi,:)%count
      END DO
    END DO
    
  ELSE
    PRINT *, "Set the type of problem (1, 2, 3)"
  END IF
  
  END SUBROUTINE
  END MODULE
  