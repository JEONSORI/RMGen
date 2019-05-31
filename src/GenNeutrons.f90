MODULE NEUTRON_MOD
  USE Type_MOD
  USE Config_MOD
  USE Util_MOD
  USE DEF_MOD
  IMPLICIT NONE
    
 !   TYPE Neutron_TYPE
 !     LOGICAL :: kill = .FALSE.
 !     INTEGER :: ring, g, iazi, myFsr = 0, FsrIdxInFxr = 0
 !     REAL :: x, y, z, polarc, polars, azic, azis, w = 1., DTC, theta
 !     TYPE(GlobalCell_TYPE), POINTER :: myCell
 !   CONTAINS
 !     PROCEDURE :: Cast => GetIntersection
 !     PROCEDURE :: Move => MovePosition
 !     PROCEDURE :: Copy => CopyNeutron
 !     GENERIC :: ASSIGNMENT(=) => Copy
 !   END TYPE
    
    TYPE Bank_TYPE
      INTEGER :: index = 0
      TYPE(Neutron_TYPE) :: neutron
      TYPE(Bank_TYPE), POINTER :: next
    CONTAINS
      PROCEDURE :: Push
      PROCEDURE :: Pop
    END TYPE
    
    TYPE(Bank_TYPE) :: Bank_Head
    TYPE(BanK_TYPE) :: Bank_Fis
    
    CONTAINS
    
    SUBROUTINE GenNeutrons    ! 190214 (Gen "Incoming" Neutrons)
    
    IMPLICIT NONE
    
    INTEGER :: i, icel, iz, ir, ig, iazi
    REAL :: xi, r, theta, cdf
    REAL :: rad, radsqr
!    TYPE(Neutron_TYPE) :: Neutron
    
    CALL RANDOM_SEED()
    
    ALLOCATE(NeutronType(Nparticle))
    
    IF (problem .NE.  3) THEN   ! starting from a surface
      DO i = 1,  Nparticle
        NeutronType(i)%kill = .FALSE.
        CALL RANDOM_NUMBER(xi)
        icel = int(xi*Core%ntotcell) + 1
  !      iz = Core%GlobalPlaneMap(icel)
        NeutronType(i)%myCell => Core%GlobalCell(icel)
        CALL RANDOM_NUMBER(xi); NeutronType(i)%x = Core%GlobalCell(icel)%xbd(1) + PITCH*xi
        CALL RANDOM_NUMBER(xi); NeutronType(i)%y = 0   ! Core%GlobalCell(icel)%ybd(1) - PITCH*xi !!!!! 190214
  !      CALL RANDOM_NUMBER(xi); Neutron%z = Core%GlobalCell(icel)%zbd(1) + Hz(iz)*xi
        CALL RANDOM_NUMBER(xi); theta = pi*xi;   ! --- 190214 2pi*xi -> pi*xi
        NeutronType(i)%azic = cos(theta); NeutronType(i)%azis = sin(theta)
        NeutronType(i)%theta = theta
        iazi = INT(xi*Nazi) + 1   ! --- 190214 
        NeutronType(i)%iazi = iazi       ! --- //
        IF (NeutronType(i)%kill) CYCLE
        
        ! --- location
        r = sqrt((NeutronType(i)%x - Core%GlobalCell(icel)%center(1))**2 + (NeutronType(i)%y - Core%GlobalCell(icel)%center(2))**2) !!!
        DO ir = 1, Core%GlobalCell(icel)%nring !!!!!
          IF (r .LT. Core%GlobalCell(icel)%rr(ir)) EXIT
        END DO
        
        NeutronType(i)%ring = ir
        
        ! --- energy
        CALL RANDOM_NUMBER(xi)
        DO ig = 1, Ng
          cdf = dble(ig)/7
          IF (xi .LT. cdf) EXIT
        END DO
        
        NeutronType(i)%g = ig
        
!        CALL Bank_Fis%Push(Neutron)
      ENDDO
    ELSE
      DO i = 1, Nparticle
        NeutronType(i)%kill = .FALSE.
        CALL RANDOM_NUMBER(xi)
        icel = int(xi*Core%ntotcell) + 1
  !      iz = Core%GlobalPlaneMap(icel)
        NeutronType(i)%myCell => Core%GlobalCell(icel) 
        CALL RANDOM_NUMBER(xi); radsqr = xi * Core%GlobalCell(icel)%rr(3)**2; rad = sqrt(radsqr)
        CALL RANDOM_NUMBER(xi); theta = 2*pi*xi
        NeutronType(i)%x = 0.5 * PITCH + rad * cos(theta)
        NeutronType(i)%y = 0.5 * PITCH + rad * sin(theta)
        IF (7./4.*pi .GT. theta .OR. theta .LT. pi/4.) THEN
          NeutronType(i)%FsrIdxInFxr = 1
        ELSE IF (theta .GT. pi/4. .AND. theta .LT. 3./4.*pi) THEN
          NeutronType(i)%FsrIdxInFxr = 2
        ELSE IF (theta .GT. 3./4.*pi .AND. theta .LT. 5./4.*pi) THEN
          NeutronType(i)%FsrIdxInFxr = 3
        ELSE IF (theta .GT. 5./4.*pi .AND. theta .LT. 7./4.*pi) THEN
          NeutronType(i)%FsrIdxInFxr = 4
        ELSE
          PRINT *, "Neutron generated at the boundary !"
        END IF
      
        CALL RANDOM_NUMBER(xi); theta = 2*pi*xi;
        NeutronType(i)%azic = cos(theta); NeutronType(i)%azis = sin(theta)
        NeutronType(i)%theta = theta
  
        IF (NeutronType(i)%kill) CYCLE
        
        ! --- location
        r = sqrt((NeutronType(i)%x - Core%GlobalCell(icel)%center(1))**2 + (NeutronType(i)%y - Core%GlobalCell(icel)%center(2))**2) !!!
        DO ir = 1, Core%GlobalCell(icel)%nring !!!!!
          IF (r .LT. Core%GlobalCell(icel)%rr(ir)) EXIT
        END DO
        
        NeutronType(i)%ring = ir
        NeutronType(i)%myFsr = (ir-1)*4 + NeutronType(i)%FsrIdxInFxr
        
        ! --- energy
        CALL RANDOM_NUMBER(xi)
        DO ig = 1, Ng
          cdf = dble(ig)/7
          IF (xi .LT. cdf) EXIT
        END DO
        
        NeutronType(i)%g = ig
      ENDDO
    END IF
    
    
    END SUBROUTINE
    

    FUNCTION GetIntersection(this) RESULT(status)
    USE DEF_MOD
    IMPLICIT NONE
    
    CLASS(Neutron_TYPE) :: this
    TYPE(Status_TYPE) :: status
    TYPE(GlobalCell_TYPE), POINTER :: Cell
    INTEGER :: ir, dir
    REAL :: d(2), r(2), P(2), H(2), I(2)
    REAL :: dx, dy, dz, dl, dist, DTS
    
    !--- d : Direction Vector
    !--- r : Position Vector to Center
    !--- P : Neutron Point
    !--- H : Perpendicular Point
    !--- I : Intersection Point
    
    IF (this%kill) RETURN
    
    Cell => this%myCell
    DTS = huge(DTS)
    
    IF (Cell%nring .EQ. 0) GOTO 1
    
    P = (/ this.x, this.y /)
    d = (/ this.azic, this.azis /)
    dist = abs(d(2) * Cell%center(1) - d(1) * Cell%center(2) - d(2) * P(1) + d(1) * P(2));  
    
    r = Cell%center - P
    H = P + dot_product(Cell%center - P, d) * d
    
    ir = this%ring - 1
    IF (ir .GT. 0) THEN
      IF (dot_product(d, r) .GT. 0) THEN
        IF (dist .LT. Cell%rr(ir)) THEN
          I = H - sqrt(Cell%rr(ir) ** 2 - dist ** 2) * d
          dl = norm2(P - I) ! / this%polarc
          IF (dl .LT. DTS) THEN
            DTS = dl
            status%ring = ir
            status%boundary = .FALSE.
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    
    ir = this%ring
    IF (ir .LE. Cell%nring) THEN
      IF (dist .LT. Cell%rr(ir)) THEN
        I = H + sqrt(Cell%rr(ir) ** 2 - dist ** 2) * d
        dl = norm2(P - I) ! / this%polarc
        IF (dl .LT. DTS) THEN
            DTS = dl
            status%ring = ir + 1
            status%boundary = .FALSE.
          ENDIF
      ENDIF
    ENDIF
    
1   CONTINUE

    IF (this%azic .LT. 0.0) THEN
      dx = Cell%xbd(1) - this%x
      dir = WEST
    ELSE
      dx = Cell%xbd(2) - this%x
      dir = EAST
    ENDIF
    dl = abs(dx / this%azic) ! / this%polarc !!!!!!! gma
    IF (dl .LT. DTS) THEN
      DTS = dl
      status%void = .FALSE.
      status%reflection = .FALSE.
      status%dir = dir
      status%boundary = .TRUE.
      status%neighbor = Cell%neighbor(dir)
      status%nout(1:3) = - normalVector_out(:, dir)
      IF (Cell%neighbor(dir) .EQ. REFCELL) THEN
        status%reflection = .TRUE.
        status%n = normalVector(:, dir)
      ELSEIF (Cell%neighbor(dir) .EQ. VOIDCELL) THEN
        status%void = .TRUE.
      ENDIF
    ENDIF
          
    IF (this%azis .LT. 0.0) THEN
      dy = Cell%ybd(2) - this%y
      dir = SOUTH
    ELSE
      dy = Cell%ybd(1) - this%y
      dir = NORTH
    ENDIF
    dl = abs(dy / this%azis) ! / this%polarc !!!!!!! gma
    IF (dl .LT. DTS) THEN
      DTS = dl
      status%void = .FALSE.
      status%reflection = .FALSE.
      status%dir = dir
      status%boundary = .TRUE.
      status%neighbor = Cell%neighbor(dir)
      status%nout = - normalVector_out(:, dir)
      IF (Cell%neighbor(dir) .EQ. REFCELL) THEN
        status%reflection = .TRUE.
        status%n = normalVector(:, dir)
      ELSEIF (Cell%neighbor(dir) .EQ. VOIDCELL) THEN
        status%void = .TRUE.
      ENDIF
    ENDIF
    
    status%distance = DTS
    status%x = this%x + DTS * this%azic           !!!!!!! gma
    status%y = this%y + DTS * this%azis           !!!!!!! gma
!    status%z = this%z + DTS * this%polars
    
    END FUNCTION

    SUBROUTINE SetLocation(this)
    
    IMPLICIT NONE
    
    CLASS(Neutron_TYPE) :: this
    TYPE(GlobalCell_TYPE), POINTER :: Cell
    INTEGER :: ir
    REAL :: dist
    
    IF (this%kill) RETURN
    
    Cell => this%myCell
    
    DO ir = 1, Cell%nring
      dist = sqrt((this%x - Cell%center(1)) ** 2 + (this%y - Cell%center(2)) ** 2)
      IF (dist .LT. Cell%rr(ir)) EXIT
    ENDDO
    
    this%ring = ir
    
    END SUBROUTINE
    SUBROUTINE MovePosition(this, status)
    
    IMPLICIT NONE
    
    CLASS(Neutron_TYPE) :: this
    TYPE(Status_TYPE) :: status
    REAL :: d(2)
    
    IF (this%kill) RETURN
    
    this%x = status%x
    this%y = status%y
!    this%z = status%z
    
    IF (Status%boundary) THEN
      this%kill = .TRUE.
      RETURN
    END IF
    
    
!    IF (status%void) THEN
!      this%kill = .TRUE.
!      RETURN
!    ENDIF
!    
!    IF (status%reflection) THEN
!!      IF (status%n(3) .EQ. 0.0) THEN       ! --- xy
!        d = (/ this%azic, this%azis /)
!        d = d - 2.0 * dot_product(d, status%n(1 : 2)) * status%n(1 : 2)
!        this%azic = d(1); this%azis = d(2)
!!      ELSE
!!        this%polars = - this%polars        ! --- z
!!      ENDIF
!      RETURN
!    ENDIF
!    
!    IF (status%boundary) THEN
!      this%myCell => Core%GlobalCell(status%neighbor)
!      CALL SetLocation(this)
!      RETURN
!    ENDIF
    
    this%ring = status%ring
    
    END SUBROUTINE
    
    SUBROUTINE Push(this, Neutron)
    
    IMPLICIT NONE
    
    CLASS(Bank_TYPE) :: this
    TYPE(Neutron_TYPE) :: Neutron
    TYPE(Bank_TYPE), POINTER :: ptr
    
    this%index = this%index + 1
    
    ptr => this%next
    ALLOCATE(this%next)
    this%next%next => ptr
    
    this%next%neutron = neutron
    this%next%index = this%index
    
    END SUBROUTINE
    
    FUNCTION Pop(this) RESULT(Neutron)
    
    IMPLICIT NONE
    
    CLASS(Bank_TYPE) :: this
    TYPE(Neutron_TYPE) :: Neutron
    TYPE(Bank_TYPE), POINTER :: ptr
    
    IF (this%index .EQ. 0) RETURN
    
    this%index = this%index - 1
    
    neutron = this%next%neutron
    
    ptr => this%next%next
    DEALLOCATE(this%next)
    this%next => ptr
    
    END FUNCTION
    
    SUBROUTINE CopyNeutron(this, n)
    
    IMPLICIT NONE
    
    CLASS(Neutron_TYPE), INTENT(OUT) :: this
    TYPE(Neutron_TYPE), INTENT(IN) :: n
    
    this%kill = n%kill
    this%ring = n%ring
    this%x = n%x
    this%y = n%y
    this%z = n%z
    this%azic = n%azic 
    this%azis = n%azis
    this%polarc = n%polarc
    this%polars = n%polars
    this%myCell => n%myCell
    this%g = n%g
    this%w = n%w
    this%iazi = n%iazi
    this%theta = n%theta
    END SUBROUTINE
    
  END MODULE
  