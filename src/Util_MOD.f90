MODULE Util_MOD
  IMPLICIT NONE
  
  CONTAINS
  FUNCTION Nword(oneline)

  CHARACTER*300 :: ONELINE
  INTEGER :: Nword, na, nb, i
  
  na = 0; nb = 0; Nword = 0
  DO i = 1, 300
    IF(oneline(i:i) .EQ. "!") EXIT
    IF((oneline(i:i) .NE. " ") .AND. (oneline(i:i) .NE. "/")) nb = 1
    IF(((oneline(i:i) .NE. " ") .AND. (oneline(i:i) .NE. "/")) .AND. (na .NE. nb)) Nword = Nword+1
    na = nb; nb = 0
  END DO
  END FUNCTION
  
  SUBROUTINE Ctrim(tarr, inarr, front, back) !'   ab / cde  ' -> 'ab' + 'cde  '
     IMPLICIT NONE
     CHARACTER :: tarr
     CHARACTER*300 :: inarr, front, back
     INTEGER :: i, j
     i=2
     front = " "; back = " "
     DO WHILE( inarr(i:i) .NE. tarr )
         i=i+1
         IF( i .GT. 100 )THEN
             front='!'
             back=inarr
             RETURN
         ENDIF
     ENDDO
     j=i+1
     i=i-1
     DO WHILE( inarr(j:j) .EQ. ' ' )
         j=j+1
     ENDDO
     front=' '
     front(1:i)=inarr(1:i)
     back=inarr(j:300)
  END SUBROUTINE
  
  SUBROUTINE Ccount(tarr, inarr, count)
      IMPLICIT NONE
      CHARACTER*300:: inarr
      CHARACTER :: tarr
      INTEGER :: count
      INTEGER :: i
      count=0
      DO i = 1, 300        
          IF( inarr(i:i) .EQ. tarr )THEN
              count=count+1
          ENDIF
      ENDDO
  END SUBROUTINE
  
  SUBROUTINE DIV_arr(tarr, ONELINE, Ntarr, itarr, front, mid, back) ! 'ab/cde/fg' -> 'ab/' + 'cde' + '/fg'
  IMPLICIT NONE
  
  INTEGER :: Ntarr, itarr, a0, a, i
  CHARACTER(1) :: tarr
  CHARACTER*300 :: ONELINE, front, mid, back, ONELINE_
  
  a0 = 0
  ONELINE_ = ONELINE
  DO i = 1, Ntarr
    a = INDEX(ONELINE_, tarr)
    IF (i .NE. itarr) THEN
      front(a0+1 : a0+a) = ONELINE_(1:a)
      back = ONELINE_(a+1:300)
      ONELINE_ = back
      a0 = a0+a
    ELSE
      mid = ONELINE_(1:a-1)
      back = ONELINE_(a:300)
      EXIT
    END IF
  END DO
  END SUBROUTINE
  
  FUNCTION IFnumeric(aline)
  logical :: IFnumeric
  CHARACTER*300 :: aline, ONELINE
  INTEGER :: i, mxncol, iascii
  
  mxncol = LEN_TRIM(aline)
  
  IFnumeric=.FALSE.
  oneline=aline
  DO i=1,mxncol
     iascii=ichar(oneline(i:i))
     IF(oneline(i:i).ne." " .and. iascii.ne.9 ) then  !determine IF the first character is numeric
        IF((iascii-48)*(iascii-57) .le. 0) IFnumeric=.TRUE.
        IF(oneline(i:i) .EQ. '+') IFnumeric=.TRUE.
        IF(oneline(i:i) .EQ. '-') IFnumeric=.TRUE.
        return
     ENDIF
  ENDDO
  
  RETURN
  END FUNCTION
  
  FUNCTION IFnumeric1(a)
  logical :: IFnumeric1
  CHARACTER*1 :: a
  INTEGER :: i, mxncol, iascii
  
  
  IFnumeric1=.FALSE.
  iascii=ichar(a)
  IF(a.ne." " .and. iascii.ne.9 ) then  !determine IF the first character is numeric
     IF((iascii-48)*(iascii-57) .le. 0) IFnumeric1=.TRUE.
     IF(a .EQ. '+') IFnumeric1=.TRUE.
     IF(a .EQ. '-') IFnumeric1=.TRUE.
     return
  ENDIF  
  END FUNCTION
  
FUNCTION nfields(aline)
character*300  aline
logical nonblankd,nonblank,multidata
INTEGER :: nfields, NCOL, N ,NMULT, MULTIDCOL, i
!
nonblankd=.FALSE.
multidata=.FALSE.
!oneline=aline
ncol=len_trim(aline)
n=0
DO i=ncol,1,-1
   IF(aline(i:i).eq.' ' .or. ichar(aline(i:i)).eq.9) then !ichar(tab)=9
      nonblank=.false.
   else
      IF(aline(i:i).eq.'*') then
         multidata=.true.
         multidcol=i
      ENDIF
      nonblank=.true.
      IF((aline(i:i).eq.'!') .OR. (aline(i:i) .eq. "/")) then
         n=-1
         nonblankd=.true.
      ENDIF
   ENDIF
   IF((.not.nonblankd.and.nonblank) .or. (nonblankd.and..not.nonblank)) then
      n=n+1
      IF(multidata.and.(nonblankd.and. .not.nonblank)) then
         read(aline(i+1:multidcol-1),*) nmult
         n=n+(nmult-1)*2
         multidata=.false.
      ENDIF
   ENDIF   
   nonblankd=nonblank
ENDDO
IF(mod(n,2).ne.0) then
  nfields=n/2+1
else
  nfields=n/2
ENDIF
!
return
END FUNCTION
  
  
  SUBROUTINE toupper(aa)
  INTEGER :: lenaa, i, ia
  INTEGER, PARAMETER :: INDXA=97,IDNXZ=122
  CHARACTER  aa*(*)

  lenaa=len_trim(aa)
  i=1
  DO while (aa(i:i).ne.' ' .and. i.le.lenaa)
     ia=ichar(aa(i:i))
     IF(ia.ge.INDXA) aa(i:i)=char(ia-32)
     i=i+1
     IF(i.gt.lenaa) return
  ENDDO
  END SUBROUTINE
  
FUNCTION icolfield(aline,ithfield)
  character(300)  aline
  INTEGER :: n, ncol, ifield, i, ithfield, icolfield
  logical nonblankd,nonblank,multidata
  nonblankd=.false.
  multidata=.false.
  ncol=len_trim(aline)
  n=0
  IF(aline(1:1).ne.' ' .and. ichar(aline(1:1)).ne.9) then !ichar(tab)=9
    IField=1
  else
    IField=0
  ENDIF
  DO i=2,256
     IF((aline(i-1:i-1).eq.' ' .or. ichar(aline(i-1:i-1)).eq.9).and. &
       (aline(i:i).ne.' ' .and. ichar(aline(i:i)).ne.9)) then
        IField=IField+1
     ENDIF
     IF(IField.eq.ithfield) exit
  ENDDO
  icolfield=min(i,256)
  return
END FUNCTION

SUBROUTINE fndchara(aline,ipos,nchar,onec)
INTEGER :: nchar, ipos(100), i, nstr
character(300) aline
character(1) onec
nchar=0
nstr=len_trim(aline)
DO i=1,nstr
  IF(aline(i:i).eq.'!') exit
  IF(aline(i:i).eq.onec) then
    nchar=nchar+1
    ipos(nchar)=i
  ENDIF
ENDDO
ipos(nchar+1)=nstr+1
return
END SUBROUTINE

END MODULE