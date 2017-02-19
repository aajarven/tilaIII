  LOGICAL FUNCTION LSAME (CA, CB)
!
!! LSAME tests two characters to determine if they are the same ...
!            letter, except for case.
!
!***LIBRARY   SLATEC
!***CATEGORY  R, N3
!***TYPE      LOGICAL (LSAME-L)
!***KEYWORDS  CHARACTER COMPARISON, LEVEL 2 BLAS, LEVEL 3 BLAS
!***AUTHOR  Hanson, R., (SNLA)
!           Du Croz, J., (NAG)
!***DESCRIPTION
!
!  LSAME  tests if CA is the same letter as CB regardless of case.
!  CB is assumed to be an upper case letter. LSAME returns .TRUE. if
!  CA is either the same as CB or the equivalent lower case letter.
!
!  N.B. This version of the code is correct for both ASCII and EBCDIC
!       systems.  Installers must modify the routine for other
!       character-codes.
!
!       For CDC systems using 6-12 bit representations, the system-
!       specific code in comments must be activated.
!
!  Parameters
!  ==========
!
!  CA     - CHARACTER*1
!  CB     - CHARACTER*1
!           On entry, CA and CB specify characters to be compared.
!           Unchanged on exit.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   860720  DATE WRITTEN
!   910606  Modified to meet SLATEC prologue standards.  Only comment
!           lines were modified.  (BKS)
!   910607  Modified to handle ASCII and EBCDIC codes.  (WRB)
!   930201  Tests for equality and equivalence combined.  (RWC and WRB)
!***END PROLOGUE  LSAME
!     .. Scalar Arguments ..
  CHARACTER CA*1, CB*1
!     .. Local Scalars ..
  INTEGER IOFF
  LOGICAL FIRST
!     .. Intrinsic Functions ..
  INTRINSIC ICHAR
!     .. Save statement ..
  SAVE FIRST, IOFF
!     .. Data statements ..
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  LSAME
  if (FIRST) IOFF = ICHAR('a') - ICHAR('A')
!
  FIRST = .FALSE.
!
!     Test if the characters are equal or equivalent.
!
  LSAME = (CA == CB) .OR. (ICHAR(CA)-IOFF == ICHAR(CB))
!
  return
!
!  The following comments contain code for CDC systems using 6-12 bit
!  representations.
!
!     .. Parameters ..
!     INTEGER                ICIRFX
!     PARAMETER            ( ICIRFX=62 )
!     .. Scalar Arguments ..
!     CHARACTER*1            CB
!     .. Array Arguments ..
!     CHARACTER*1            CA(*)
!     .. Local Scalars ..
!     INTEGER                IVAL
!     .. Intrinsic Functions ..
!     INTRINSIC              ICHAR, CHAR
!     .. Executable Statements ..
!     INTRINSIC              ICHAR, CHAR
!
!     See if the first character in string CA equals string CB.
!
!     LSAME = CA(1)  ==  CB .AND. CA(1)  /=  CHAR(ICIRFX)
!
!     if (LSAME) RETURN
!
!     The characters are not identical. Now check them for equivalence.
!     Look for the 'escape' character, circumflex, followed by the
!     letter.
!
!     IVAL = ICHAR(CA(2))
!     if (IVAL >= ICHAR('A') .AND. IVAL <= ICHAR('Z')) THEN
!        LSAME = CA(1)  ==  CHAR(ICIRFX) .AND. CA(2)  ==  CB
!     ENDIF
!
!     return
!
!     End of LSAME.
!
end
