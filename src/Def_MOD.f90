MODULE DEF_MOD
  USE Type_MOD
  
  INTEGER, PARAMETER :: infid = 28, outfid = 29, outfid_tempk = 30
  INTEGER, PARAMETER :: voidcell = -1, refcell = 0, cap = 1, fis = 2, scat = 3
  INTEGER, PARAMETER :: west = 1, east = 2, north = 3, south = 4, DOWN = 5, UP = 6
  INTEGER :: problem, Nrun = 1
  INTEGER :: outfid_SS = 10, outfid_Sphi = 11, outfid_VS = 12, outfid_Vphi = 13
  REAL, PARAMETER :: pi = 4*atan(1.)
  REAL :: normalVector(3, 6), normalVector_out(3, 6), keff = 1.
  REAL, POINTER :: Response(:,:), Response2(:,:)
  DATA normalVector /  1,  0,  0,  &
                      -1,  0,  0,  &
                       0, -1,  0,  &
                       0,  1,  0,  &
                       0,  0,  1,  &
                       0,  0, -1   /
  DATA normalVector_out /  1,  0,  0,  &
                          -1,  0,  0,  &
                           0, -1,  0,  &
                           0,  1,  0,  &
                           0,  0,  1,  &
                           0,  0, -1   /
  CHARACTER*100 :: filename
  
  TYPE(PinInfo_TYPE), POINTER :: PinInfo(:)
  TYPE(CellInfo_TYPE), POINTER :: CellInfo(:)
  TYPE(AsyInfo_TYPE), POINTER :: AsyInfo(:)
  TYPE(CoreInfo_TYPE) :: Core
  TYPE(ConfigInfo_TYPE) :: ConfigInfo
  TYPE(Neutron_TYPE), POINTER :: NeutronType(:)
  TYPE(Flux_TYPE), POINTER :: Incoming(:,:,:), Outgoing(:,:,:,:), phi(:,:,:)
  TYPE(CDF_TYPE), POINTER :: CDF(:)
  
  END MODULE
