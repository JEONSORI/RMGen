MODULE Config_MOD
  USE TYPE_MOD
  INTEGER :: Ntotcyc = 0, Nactcyc = 0, Ninactcyc = 0, Nparticle = 0, Nazi = 16
  INTEGER :: Nmixture, Ng, NFuel
  INTEGER :: nCellType, nPinType, nAsyType, Npins, Nz = 1
  REAL :: ALBEDO(6)
  REAL :: PITCH
  REAL, POINTER :: hz(:), Plane(:)
  TYPE(XS_TYPE), POINTER :: XS(:)
  END MODULE
  
  