SUBROUTINE SetGeom
  USE Type_MOD
  USE Config_MOD
  USE Util_MOD
  USE DEF_MOD
  IMPLICIT NONE
  
  INTEGER :: i, j, ipin, ixa, iya, iasy, ir, ic, ix, iy, iz, icel, ipintype, icelltype, nfxr, nfsr
  CHARACTER*300 :: ONELINE
  
  
  Core%nx = Core%nxa*npins
  Core%ny = Core%nya*npins
  Core%nxy = Core%nxya*npins**2
  Core%ntotcell = Core%nxy*Nz
  
  ALLOCATE(Core%GlobalPinMap(0:Core%ny+1, 0:Core%nx+1))
  Core%GlobalPinMap = voidcell
  ALLOCATE(Core%GlobalPinType(Core%ny,Core%nx))
  Core%GlobalPinType = 0
  
  IF (albedo(1) .EQ. 0.) Core%GlobalPinMap(:,0) = refcell
  IF (albedo(2) .EQ. 0.) Core%GlobalPinMap(:,Core%nx+1) = refcell
  IF (albedo(3) .EQ. 0.) Core%GlobalPinMap(0,:) = refcell
  IF (albedo(4) .EQ. 0.) Core%GlobalPinMap(Core%ny+1,:) = refcell
  
  ipin = 0
  
  DO iya = 1, Core%nya
    DO ixa = 1, Core%nxa
      iasy = Core%CoreIdx(iya,ixa) !!!!!
      IF (iasy .EQ. 0) CYCLE
      DO ir = 1, npins
        DO ic = 1, npins
          ix = ic + (ixa-1)*npins
          iy = ir + (iya-1)*npins
          ipin = ipin + 1
          Core%GlobalPinMap(iy,ix) = ipin
          Core%GlobalPinType(iy,ix) = AsyInfo(iasy)%AsyConfig(ir,ic)
        END DO
      END DO
    END DO
  END DO
  
  ALLOCATE(Core%GlobalPlaneMap(Core%ntotcell))
  ALLOCATE(Core%GlobalCell(Core%ntotcell))
  
  icel = 0
  
  DO iz = 1, Nz
    DO iy = 1, Core%ny
      DO ix = 1, Core%nx
        ipin = Core%GlobalPinMap(iy,ix)
        ipintype = Core%GlobalPinType(iy,ix)
        IF (ipin .EQ. voidcell) CYCLE
        icel = ipin + (iz-1)*Core%nxy
        icelltype = PinInfo(ipintype)%cell(iz)
        Core%GlobalPlaneMap(icel) = iz
        Core%GlobalCell(icel)%icelltype = icelltype
        Core%GlobalCell(icel)%ix = ix
        Core%GlobalCell(icel)%iy = iy
!        Core%GlobalCell(icel)%iz = iz
        Core%GlobalCell(icel)%ipin = ipin
        Core%GlobalCell(icel)%xbd(1) = (ix-1)*PITCH
        Core%GlobalCell(icel)%xbd(2) = ix*PITCH
        Core%GlobalCell(icel)%ybd(1) = (Core%ny-iy+1)*PITCH
        Core%GlobalCell(icel)%ybd(2) = (Core%ny-iy)*PITCH
!        Core%GlobalCell(icel)%zbd(1) = Plane(iz-1)
!        Core%GlobalCell(icel)%zbd(2) = Plane(iz)
        Core%GlobalCell(icel)%center(1) = Core%GlobalCell(icel)%xbd(1) + 0.5 * PITCH
        Core%GlobalCell(icel)%center(2) = Core%GlobalCell(icel)%ybd(1) - 0.5 * PITCH
        nfxr = CellInfo(icelltype)%nfxr
        Core%GlobalCell(icel)%nfsr = CellInfo(icelltype)%nfsr
        Core%GlobalCell(icel)%nring = CellInfo(icelltype)%nring
        ALLOCATE(Core%GlobalCell(icel)%rr(Nfxr-1))
        Core%GlobalCell(icel)%rr = CellInfo(icelltype)%rr
        Core%GlobalCell(icel)%lrec = CellInfo(icelltype)%lrec
        
        Core%GlobalCell(icel)%neighbor(1) = Core%GlobalPinMap(iy,ix-1)
        IF (Core%GlobalPinMap(iy,ix-1) .GT. 0) THEN
          Core%GlobalCell(icel)%neighbor(1) = Core%GlobalPinMap(iy,ix-1) + (iz-1)*Core%nxy
        END IF
        
        
        Core%GlobalCell(icel)%neighbor(2) = Core%GlobalPinMap(iy,ix+1)
        IF (Core%GlobalPinMap(iy,ix+1) .GT. 0) THEN
          Core%GlobalCell(icel)%neighbor(2) = Core%GlobalPinMap(iy,ix+1) + (iz-1)*Core%nxy
        END IF
        
        Core%GlobalCell(icel)%neighbor(3) = Core%GlobalPinMap(iy-1,ix)
        IF (Core%GlobalPinMap(iy-1,ix) .GT. 0) THEN
          Core%GlobalCell(icel)%neighbor(3) = Core%GlobalPinMap(iy-1,ix) + (iz-1)*Core%nxy
        END IF
        
        Core%GlobalCell(icel)%neighbor(4) = Core%GlobalPinMap(iy+1,ix)
        IF (Core%GlobalPinMap(iy+1,ix) .GT. 0) THEN
          Core%GlobalCell(icel)%neighbor(4) = Core%GlobalPinMap(iy+1,ix) + (iz-1)*Core%nxy
        END IF
        
!        Core%GlobalCell(icel)%neighbor(5) = icel - Core%nxy
!        Core%GlobalCell(icel)%neighbor(6) = icel + Core%nxy
        
!        IF (iz .EQ. 1) THEN
!          IF (albedo(5) .EQ. 0.) THEN
!            Core%GlobalCell(icel)%neighbor(5) = refcell
!          ELSE IF (albedo(5) .EQ. 0.5) THEN
!            Core%GlobalCell(icel)%neighbor(5) = voidcell
!          END IF
!        END IF
!        IF (iz .EQ. Nz) THEN
!          IF (albedo(6) .EQ. 0.) THEN
!            Core%GlobalCell(icel)%neighbor(6) = refcell
!          ELSE IF (albedo(6) .EQ. 0.5) THEN
!            Core%GlobalCell(icel)%neighbor(6) = voidcell
!          END IF
!        END IF
      END DO
    END DO
  END DO
  
  END SUBROUTINE
  
          
        
        
        