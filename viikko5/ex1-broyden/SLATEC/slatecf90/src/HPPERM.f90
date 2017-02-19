subroutine HPPERM (HX, N, IPERM, WORK, IER)
!
!! HPPERM rearranges an array according to a permutation vector.
!
!***LIBRARY   SLATEC
!***CATEGORY  N8
!***TYPE      CHARACTER (SPPERM-S, DPPERM-D, IPPERM-I, HPPERM-H)
!***KEYWORDS  APPLICATION OF PERMUTATION TO DATA VECTOR
!***AUTHOR  McClain, M. A., (NIST)
!           Rhoads, G. S., (NBS)
!***DESCRIPTION
!
!         HPPERM rearranges the data vector HX according to the
!         permutation IPERM: HX(I) <--- HX(IPERM(I)).  IPERM could come
!         from one of the sorting routines IPSORT, SPSORT, DPSORT or
!         HPSORT.
!
!     Description of Parameters
!         HX - input/output -- character array of values to be
!                 rearranged.
!         N - input -- number of values in character array HX.
!         IPERM - input -- permutation vector.
!         WORK - character variable which must have a length
!                   specification at least as great as that of HX.
!         IER - output -- error indicator:
!             =  0  if no error,
!             =  1  if N is zero or negative,
!             =  2  if work array is not long enough,
!             =  3  if IPERM is not a valid permutation.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   901004  DATE WRITTEN
!   920507  Modified by M. McClain to revise prologue text and to add
!           check for length of work array.
!***END PROLOGUE  HPPERM
  INTEGER N, IPERM(*), I, IER, INDX, INDX0, ISTRT
  CHARACTER*(*) HX(*), WORK
!***FIRST EXECUTABLE STATEMENT  HPPERM
  IER=0
  if ( N < 1)THEN
     IER=1
     call XERMSG ('SLATEC', 'HPPERM', &
      'The number of values to be rearranged, N, is not positive.', &
      IER, 1)
     return
  end if
  if ( LEN(WORK) < LEN(HX(1)))THEN
     IER=2
     call XERMSG ('SLATEC', 'HPPERM', &
      'The length of the work variable, WORK, is too short.',IER,1)
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
     IER=3
     call XERMSG ('SLATEC', 'HPPERM', &
      'The permutation vector, IPERM, is not valid.', IER, 1)
     return
  100 CONTINUE
!
!     REARRANGE THE VALUES OF HX
!
!     USE THE IPERM VECTOR AS A FLAG.
!     if IPERM(I) > 0, THEN THE I-TH VALUE IS IN CORRECT LOCATION
!
  DO 330 ISTRT = 1 , N
     if (IPERM(ISTRT)  >  0) GOTO 330
     INDX = ISTRT
     INDX0 = INDX
     WORK = HX(ISTRT)
  320    CONTINUE
     if (IPERM(INDX)  >=  0) GOTO 325
        HX(INDX) = HX(-IPERM(INDX))
        INDX0 = INDX
        IPERM(INDX) = -IPERM(INDX)
        INDX = IPERM(INDX)
        GOTO 320
  325    CONTINUE
     HX(INDX0) = WORK
  330 CONTINUE
!
  return
end
