SUBROUTINE Scaninput
  USE Util_MOD
  USE Config_MOD
  USE Type_MOD
  USE DEF_MOD
  IMPLICIT NONE
  
  INTEGER :: nLineField, Nword_, ipin, iasy
  CHARACTER*300 :: ONELINE
  CHARACTER*100 :: cardname
  OPEN(infid, file = "input.txt")
  
  DO WHILE(.TRUE.)
    READ(infid, '(A)'), ONELINE
    IF(ONELINE(1:1) .EQ. "!") CYCLE; IF(ONELINE .EQ. " ") CYCLE
    IF(ONELINE(1:1) .EQ. ".") EXIT
    READ(ONELINE,*), cardname; CALL toupper(cardname)
    nLineField = nfields(ONELINE) - 1
    SELECT CASE(cardname)
    CASE ("CELL")
      nCellType = nCellType + 1
    CASE ("PIN")
      nPinType = nPinType + 1
    CASE ("ASSEMBLY")
      nAsyType = nAsyType + 1
    CASE ("AX_MESH")
      Nz = nLineField
      ALLOCATE(Hz(0:Nz))
      Hz(0) = 0.
    CASE ("RAD_CONF")
      DO WHILE(.TRUE.)
        READ(infid, '(A)'), ONELINE
        IF (.NOT. IFnumeric(ONELINE)) EXIT
        Nword_ = nfields(ONELINE)
        Core%nxa = max(Core%nxa, Nword_)
        Core%nya = Core%nya + 1
        Core%nxya = Core%nxya + Nword_
      END DO
    CASE ("PARTICLE")
      READ(ONELINE,*), cardname, Nparticle
    CASE ("ACTIVE")
      READ(ONELINE,*), cardname, Nactcyc
    CASE ("INACTIVE")
      READ(ONELINE,*), cardname, Ninactcyc
    CASE ("NRUN")
      READ(ONELINE,*), cardname, nrun
    CASE ("RM_TYPE")
      READ(ONELINE,*), cardname, problem
    END SELECT
  END DO
  
  Ntotcyc = Nactcyc + Ninactcyc
  ALLOCATE(CellInfo(nCellType), PinInfo(nPinType), AsyInfo(nAsyType))
  ALLOCATE(Core%CoreIdx(Core%nya,Core%nxa), Core%CoreMap(Core%nxya))
  DO ipin = 1, nPinType
    ALLOCATE(PinInfo(ipin)%Cell(Nz))
    PinInfo(ipin)%Cell = 0.
  END DO
  
  END SUBROUTINE
  
    