MODULE ReadInput_MOD
  USE Type_MOD
  USE Config_MOD
  USE Util_MOD
  USE DEF_MOD
  IMPLICIT NONE
  
  INTEGER :: nLineField, Nword_
  CHARACTER*300 :: ONELINE
  CHARACTER*100 :: blockname, cardname
  CONTAINS
  
  SUBROUTINE Readinput
  IMPLICIT NONE
  
  
  REWIND(infid)
  
  DO WHILE(.TRUE.)
    READ(infid, '(A)'), ONELINE
    IF(ONELINE(1:1) .EQ. "!") CYCLE; IF(ONELINE .EQ. " ") CYCLE
    IF(ONELINE(1:1) .EQ. " ") CYCLE
    IF(ONELINE(1:1) .EQ. ".") EXIT
    READ(ONELINE,*), blockname; CALL toupper(blockname)
    SELECT CASE(blockname)
    CASE ("CASEID")
      READ(ONELINE,*), blockname, filename
    CASE ("GEOM")
      CALL ReadGeomBlock
    CASE ("XS")
      CALL ReadXSBlock
    END SELECT
  END DO
  
  END SUBROUTINE
  
  SUBROUTINE ReadGeomBlock
  IMPLICIT NONE
  INTEGER :: i, icel, nfxr, nfsr, nr, nspt, ipos(100), ipin, ix, iy, iasy, k, nya, nxa, nxya, jfr,jto, ireg(100) = 0, ifxr, ifr, ito, mixture = 0
  CHARACTER*300 :: front, back
  TYPE(CellInfo_Type), POINTER :: CELL(:)
  TYPE(PinInfo_TYPE), POINTER :: PIN(:)
  TYPE(AsyInfo_TYPE), POINTER :: ASY(:)
  
  CELL => CellInfo
  PIN => PinInfo
  ASY => AsyInfo
  
  DO WHILE(.TRUE.)
    READ(infid, '(A)'), ONELINE
    IF(ONELINE(1:1) .EQ. "!") CYCLE; IF(ONELINE .EQ. " ") CYCLE
    IF(ONELINE(1:1) .NE. " ") GO TO 20
    IF(ONELINE(1:1) .EQ. ".") EXIT
    READ(ONELINE,*), cardname; CALL toupper(cardname)
    nLineField = nfields(ONELINE) - 1
    SELECT CASE(cardname)
    CASE ("NPINS")
      READ(ONELINE,*), cardname, npins
    CASE ("PITCH")
      READ(ONELINE,*), cardname, PITCH
    CASE ("AX_MESH")
      READ(ONELINE,*), cardname, (Hz(i), i=1,Nz)
      ALLOCATE(Plane(0:Nz))
      Plane(0) = 0.
      DO i = 1, Nz
        Plane(i) = SUM(Hz(1:i))
      END DO
    CASE ("CELL")
      READ(ONELINE,*), cardname, icel
      CellInfo(icel)%icel = icel
      CALL fndchara(ONELINE,ipos,nspt,"/")
      IF (nspt .EQ. 1) THEN
        CALL Ctrim("/",ONELINE,front,back)
        nfxr = nfields(back)
        CellInfo(icel)%nfxr = nfxr
        CellInfo(icel)%nfsr = 4*nfxr
        nr = nfxr-1
        CellInfo(icel)%nring = nr
        Nfuel = nr
        ALLOCATE(CellInfo(icel)%rr(nr), CellInfo(icel)%mat(CellInfo(icel)%nfsr))
        READ(front,*), cardname, icel, (CellInfo(icel)%rr(i), i = 1,nr)
        READ(back,*), ireg(1:nfxr)
        DO ifxr = 1, nfxr
          ifr = (ifxr-1)*4 + 1; ito = ifxr*4
          mixture = ireg(ifxr)
          CellInfo(icel)%mat(ifr:ito) = mixture
        END DO
      ELSE IF (nspt .EQ. 2) THEN
        CellInfo(icel)%lrec = .TRUE.
        nfxr = 1; nfsr = 1
        CellInfo(icel)%nfxr = nfxr; CellInfo(icel)%nfsr = nfsr
        ALLOCATE(CellInfo(icel)%mat(nfsr))
        READ(ONELINE(ipos(2)+1:300),*), CellInfo(icel)%mat
      ELSE
        PRINT '(A)', "Input error - Cell card"
      END IF
    CASE ("PIN")
      READ(ONELINE,*), cardname, ipin
      REAd(ONELINE,*), cardname, PinInfo(ipin)%ipin, PinInfo(ipin)%Cell
    CASE ("ASSEMBLY")
      READ(ONELINE,*), cardname, iasy
      AsyInfo(iasy)%iasy = iasy; AsyInfo(iasy)%nx = Npins; AsyInfo(iasy)%ny = Npins; AsyInfo(iasy)%nxy = Npins**2
      ALLOCATE(AsyInfo(iasy)%Pin(Npins**2), AsyInfo(iasy)%Pin2Idx(Npins,Npins), AsyInfo(iasy)%AsyConfig(Npins,Npins))
      AsyInfo(iasy)%Pin = 0; AsyInfo(iasy)%Pin2Idx = 0; AsyInfo(iasy)%AsyConfig = 0
      DO iy = 1, Npins
        READ(infid,'(A)'), ONELINE
        READ(ONELINE,*), (AsyInfo(iasy)%AsyConfig(iy,ix), ix = 1, Npins)
      END DO
      k = 0
      DO iy = 1, Npins
        DO ix = 1, Npins
          k = k + 1
          AsyInfo(iasy)%Pin(k) = AsyInfo(iasy)%AsyConfig(iy,ix)
          AsyInfo(iasy)%Pin2Idx(iy,ix) = k
        END DO
      END DO
    CASE ("RAD_CONF")
      jfr = 0; jto = 0
      DO iy = 1, Core%nya
        READ(infid,'(A)'), ONELINE
        nxa = nfields(ONELINE)
        READ(ONELINE,*), (Core%CoreIdx(iy,ix), ix = 1, nxa)
        jfr = jto + 1; jto = jfr + nxa - 1
        Core%CoreMap(jfr:jto) = Core%CoreIdx(iy,:)
      END DO
    CASE ("ALBEDO")
      READ(ONELINE,*), cardname, (albedo(i), i = 1, 6)    ! W E N S B T
    END SELECT
  END DO
  
20 CONTINUE
   BACKSPACE(infid)
  
  END SUBROUTINE
  
  SUBROUTINE ReadXSBlock
  IMPLICIT NONE
  INTEGER :: infid_xs = 1, imix, g, ig
  
  OPEN(infid_xs, file = "XS_C5G7.lib")
  DO WHILE(.TRUE.)
    READ(infid_xs, '(A)', END = 21), ONELINE
    IF(ONELINE(1:1) .EQ. "!") CYCLE; IF(ONELINE .EQ. " ") CYCLE
    READ(ONELINE,*), cardname; CALL toupper(cardname)
    SELECT CASE(cardname)
    CASE("NMIXTURE")
      READ(ONELINE,*), cardname, Nmixture
      ALLOCATE(XS(Nmixture), CDF(Nmixture))
    CASE("NG")
      READ(ONELINE,*), cardname, Ng
      DO imix = 1, Nmixture
        ALLOCATE(CDF(imix)%reactiontype(3,Ng), CDF(imix)%chi(Ng), CDF(imix)%scat(Ng,Ng))
        ALLOCATE(XS(imix)%Sigt(Ng), Xs(imix)%Sigtr(Ng), Xs(imix)%Siga(Ng), Xs(imix)%Sigr(Ng), Xs(imix)%Sigf(Ng), Xs(imix)%nu(Ng), Xs(imix)%chi(Ng), XS(imix)%Sigs(Ng,Ng))
        XS(imix)%Sigt = 0.; Xs(imix)%Sigtr = 0.; Xs(imix)%Sigr = 0.; Xs(imix)%Siga = 0.; Xs(imix)%Sigf = 0.; Xs(imix)%nu = 0.; Xs(imix)%chi = 0.;  XS(imix)%Sigs = 0.
      END DO
    CASE("MIXTURE")
      READ(ONELINE,*), cardname, imix
      DO ig = 1, Ng
        READ(infid_xs,*), g, XS(imix)%Sigt(ig), Xs(imix)%Sigtr(ig), Xs(imix)%Siga(ig), Xs(imix)%Sigr(ig), Xs(imix)%Sigf(ig), Xs(imix)%nu(ig), Xs(imix)%chi(ig)
      END DO
      DO ig = 1, Ng
        READ(infid_xs,*), g, XS(imix)%Sigs(ig,:)
      END DO
    END SELECT
  END DO
21 CONTINUE  
  END SUBROUTINE
  
  END MODULE
  
  
  
    