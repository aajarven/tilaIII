subroutine XERBLA (SRNAME, INFO)
!
!! XERBLA is the error handler for the Level 2 and Level 3 BLAS Routines.
!
!***LIBRARY   SLATEC
!***CATEGORY  R3
!***TYPE      ALL (XERBLA-A)
!***KEYWORDS  ERROR MESSAGE
!***AUTHOR  Dongarra, J. J., (ANL)
!***DESCRIPTION
!
!  Purpose
!  =======
!
!  It is called by Level 2 and 3 BLAS routines if an input parameter
!  is invalid.
!
!  Parameters
!  ==========
!
!  SRNAME - CHARACTER*6.
!           On entry, SRNAME specifies the name of the routine which
!           called XERBLA.
!
!  INFO   - INTEGER.
!           On entry, INFO specifies the position of the invalid
!           parameter in the parameter-list of the calling routine.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   860720  DATE WRITTEN
!   910610  Routine rewritten to serve as an interface between the
!           Level 2 and Level 3 BLAS routines and the SLATEC error
!           handler XERMSG.  (BKS)
!***END PROLOGUE  XERBLA
!
!     ..    Scalar Arguments ..
  INTEGER            INFO
  CHARACTER*6        SRNAME
  CHARACTER*2        XERN1
!
!***FIRST EXECUTABLE STATEMENT  XERBLA
!
  WRITE (XERN1, '(I2)') INFO
  call XERMSG ('SLATEC', SRNAME, 'On entry to '//SRNAME// &
               ' parameter number '//XERN1//' had an illegal value', &
               INFO,1)
!
  return
!
!     End of XERBLA.
!
end
