MODULE Type_MOD
  IMPLICIT NONE
  
  
  TYPE CoreInfo_TYPE
    INTEGER :: nxa = 0, nya = 0, nxya = 0
    INTEGER :: nx = 0, ny = 0, nxy = 0, ntotcell = 0
    INTEGER, POINTER :: CoreMap(:), CoreIdx(:,:)
    INTEGER, POINTER :: GlobalPinMap(:,:), GlobalPinType(:,:), GlobalPlaneMap(:)
    TYPE(GlobalCell_TYPE), POINTER :: GlobalCell(:)
  END TYPE
  
  TYPE GlobalCell_TYPE
    INTEGER :: icelltype, nfxr, nfsr, nring
    INTEGER :: ix, iy, iz, ipin, neighbor(6) = 0
    REAL :: xbd(2), ybd(2), zbd(2), center(2)
    REAL, POINTER :: rr(:)
    LOGICAL :: lrec
  END TYPE

  TYPE CDF_TYPE
    REAL, POINTER :: reactiontype(:,:), chi(:), scat(:,:)
  END TYPE
  
  TYPE XS_TYPE
    REAL, POINTER :: Sigt(:), Sigtr(:),Sigr(:), Siga(:), Sigf(:), nu(:), chi(:)
    REAL, POINTER :: Sigs(:,:)
  END TYPE
  
  TYPE Status_TYPE
    LOGICAL :: boundary = .FALSE.
    LOGICAL :: reflection = .FALSE.
    LOGICAL :: void = .FALSE.
    INTEGER :: ring, neighbor, dir
    REAL :: distance
    REAL :: x, y, z
    REAL :: n(1:3), nout(1:3)
  END TYPE
  
  TYPE ConfigInfo_TYPE
    TYPE(CoreInfo_TYPE) :: Core
    TYPE(CellInfo_TYPE), POINTER :: CellInfo(:)
    TYPE(PinInfo_TYPE), POINTER :: PinInfo(:)
    TYPE(AsyInfo_TYPE), POINTER :: AsyInfo(:)
  END TYPE
  
  TYPE CellInfo_TYPE
    INTEGER :: icel, nfxr, nfsr, nring
    INTEGER, POINTER :: mat(:) 
    REAL, POINTER :: rr(:)
    LOGICAL :: lrec = .FALSE.
  END TYPE
  
  TYPE PinInfo_TYPE
    INTEGER :: ipin
    INTEGER, POINTER :: Cell(:)
  END TYPE
  
  TYPE AsyInfo_TYPE
    INTEGER :: iasy, nx, ny, nxy
    INTEGER, POINTER :: Pin(:), Pin2Idx(:,:), AsyConfig(:,:)
  END TYPE
  
  TYPE Neutron_TYPE
    LOGICAL :: kill = .FALSE.
    INTEGER :: ring, g, iazi, myFsr = 0, FsrIdxInFxr = 0
    REAL :: x, y, z, polarc, polars, azic, azis, w = 1., DTC, theta
    TYPE(GlobalCell_TYPE), POINTER :: myCell
  END TYPE
  
      
  TYPE Flux_TYPE
    REAL :: count = 0., mu
  END TYPE
  
  END MODULE
  
  