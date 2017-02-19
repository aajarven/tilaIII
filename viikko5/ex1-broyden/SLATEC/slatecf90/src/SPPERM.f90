subroutine SPPERM (X, N, IPERM, IER)
!
!! SPPERM rearranges a given array according to a prescribed permutation vector.
!
!***LIBRARY   SLATEC
!***CATEGORY  N8
!***TYPE      SINGLE PRECISION (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H)
!***KEYWORDS  APPLICATION OF PERMUTATION TO DATA VECTOR
!***AUTHOR  McClain, M. A., (NIST)
!           Rhoads, G. S., (NBS)
!***DESCRIPTION
!
!         SPPERM rearranges the data vector X according to the
!         permutation IPERM: X(I) <--- X(IPERM(I)).  IPERM could come
!         from one of the sorting routines IPSORT, SPSORT, DPSORT or
!         HPSORT.
!
!     Description of Parameters
!         X - input/output -- real array of values to be rearranged.
!         N - input -- number of values in real array X.
!         IPERM - input -- permutation vector.
!         IER - output -- error indicator:
!             =  0  if no error,
!             =  1  if N is zero or negative,
!             =  2  if IPERM is not a valid permutation.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   901004  DATE WRITTEN
!   920507  Modified by M. McClain to revise prologue text.
!***END PROLOGUE  SPPERM
  INTEGER N, IPERM(*), I, IER, INDX, INDX0, ISTRT
  REAL X(*), TEMP
!***FIRST EXECUTABLE STATEMENT  SPPERM
  IER=0
  if ( N < 1)THEN
     IER=1
     call XERMSG ('SLATEC', 'SPPERM', &
      'The number of values to be rearranged, N, is not positive.', &
      IER, 1)
     return
  end if
!
!     CHECK WHETHER IPERM IS A VALID PERMUTATION
!
  DO 100 I=1,N
     INDX=ABS(IPERM(I))
     if ( (INDX >= 1).AND.(INDX <= N))THEN
        if ( IPERM(INDX) > 0)THEN
           IPERM(INDX)=-IPERM(INDX)
           GOTO 100
        ENDIF
     ENDIF
     IER=2
     call XERMSG ('SLATEC', 'SPPERM', &
      'The permutation vector, IPERM, is not valid.', IER, 1)
     return
  100 CONTINUE
!
!     REARRANGE THE VALUES OF X
!
!     USE THE IPERM VECTOR AS A FLAG.
!     if IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION
!
  DO 330 ISTRT = 1 , N
     if (IPERM(ISTRT)  >  0) GOTO 330
     INDX = ISTRT
     INDX0 = INDX
     TEMP = X(ISTRT)
  320    CONTINUE
     if (IPERM(INDX)  >=  0) GOTO 325
        X(INDX) = X(-IPERM(INDX))
        INDX0 = INDX
        IPERM(INDX) = -IPERM(INDX)
        INDX = IPERM(INDX)
        GOTO 320
  325    CONTINUE
     X(INDX0) = TEMP
  330 CONTINUE
!
  return
end
