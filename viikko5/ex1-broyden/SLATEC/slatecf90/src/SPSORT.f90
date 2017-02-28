subroutine SPSORT (X, N, IPERM, KFLAG, IER)
!
!! SPSORT returns the permutation vector generated by sorting a given array ...
!  and, optionally, rearrange the elements of the array.
!            The array may be sorted in increasing or decreasing order.
!            A slightly modified quicksort algorithm is used.
!
!***LIBRARY   SLATEC
!***CATEGORY  N6A1B, N6A2B
!***TYPE      SINGLE PRECISION (SPSORT-S, DPSORT-D, IPSORT-I, HPSORT-H)
!***KEYWORDS  NUMBER SORTING, PASSIVE SORTING, SINGLETON QUICKSORT, SORT
!***AUTHOR  Jones, R. E., (SNLA)
!           Rhoads, G. S., (NBS)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   SPSORT returns the permutation vector IPERM generated by sorting
!   the array X and, optionally, rearranges the values in X.  X may
!   be sorted in increasing or decreasing order.  A slightly modified
!   quicksort algorithm is used.
!
!   IPERM is such that X(IPERM(I)) is the Ith value in the rearrangement
!   of X.  IPERM may be applied to another array by calling IPPERM,
!   SPPERM, DPPERM or HPPERM.
!
!   The main difference between SPSORT and its active sorting equivalent
!   SSORT is that the data are referenced indirectly rather than
!   directly.  Therefore, SPSORT should require approximately twice as
!   long to execute as SSORT.  However, SPSORT is more general.
!
!   Description of Parameters
!      X - input/output -- real array of values to be sorted.
!          If ABS(KFLAG) = 2, then the values in X will be
!          rearranged on output; otherwise, they are unchanged.
!      N - input -- number of values in array X to be sorted.
!      IPERM - output -- permutation array such that IPERM(I) is the
!              index of the value in the original order of the
!              X array that is in the Ith location in the sorted
!              order.
!      KFLAG - input -- control parameter:
!            =  2  means return the permutation vector resulting from
!                  sorting X in increasing order and sort X also.
!            =  1  means return the permutation vector resulting from
!                  sorting X in increasing order and do not sort X.
!            = -1  means return the permutation vector resulting from
!                  sorting X in decreasing order and do not sort X.
!            = -2  means return the permutation vector resulting from
!                  sorting X in decreasing order and sort X also.
!      IER - output -- error indicator:
!          =  0  if no error,
!          =  1  if N is zero or negative,
!          =  2  if KFLAG is not 2, 1, -1, or -2.
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified by John A. Wisniewski to use the Singleton
!           quicksort algorithm.
!   870423  Modified by Gregory S. Rhoads for passive sorting with the
!           option for the rearrangement of the original data.
!   890620  Algorithm for rearranging the data vector corrected by R.
!           Boisvert.
!   890622  Prologue upgraded to Version 4.0 style by D. Lozier.
!   891128  Error when KFLAG < 0 and N=1 corrected by R. Boisvert.
!   920507  Modified by M. McClain to revise prologue text.
!   920818  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (SMR, WRB)
!***END PROLOGUE  SPSORT
!     .. Scalar Arguments ..
  INTEGER IER, KFLAG, N
!     .. Array Arguments ..
  REAL X(*)
  INTEGER IPERM(*)
!     .. Local Scalars ..
  REAL R, TEMP
  INTEGER I, IJ, INDX, INDX0, ISTRT, J, K, KK, L, LM, LMT, M, NN
!     .. Local Arrays ..
  INTEGER IL(21), IU(21)
!     .. External Subroutines ..
  EXTERNAL XERMSG
!     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  SPSORT
  IER = 0
  NN = N
  if (NN  <  1) THEN
     IER = 1
     call XERMSG ('SLATEC', 'SPSORT', &
      'The number of values to be sorted, N, is not positive.', &
      IER, 1)
     return
  end if
  KK = ABS(KFLAG)
  if (KK /= 1 .AND. KK /= 2) THEN
     IER = 2
     call XERMSG ('SLATEC', 'SPSORT', &
      'The sort control parameter, KFLAG, is not 2, 1, -1, or -2.', &
      IER, 1)
     return
  end if
!
!     Initialize permutation vector
!
  DO 10 I=1,NN
     IPERM(I) = I
   10 CONTINUE
!
!     Return if only one value is to be sorted
!
  if (NN  ==  1) RETURN
!
!     Alter array X to get decreasing order if needed
!
  if (KFLAG  <=  -1) THEN
     DO 20 I=1,NN
        X(I) = -X(I)
   20    CONTINUE
  end if
!
!     Sort X only
!
  M = 1
  I = 1
  J = NN
  R = .375E0
!
   30 if (I  ==  J) go to 80
  if (R  <=  0.5898437E0) THEN
     R = R+3.90625E-2
  ELSE
     R = R-0.21875E0
  end if
!
   40 K = I
!
!     Select a central element of the array and save it in location L
!
  IJ = I + INT((J-I)*R)
  LM = IPERM(IJ)
!
!     If first element of array is greater than LM, interchange with LM
!
  if (X(IPERM(I))  >  X(LM)) THEN
     IPERM(IJ) = IPERM(I)
     IPERM(I) = LM
     LM = IPERM(IJ)
  end if
  L = J
!
!     If last element of array is less than LM, interchange with LM
!
  if (X(IPERM(J))  <  X(LM)) THEN
     IPERM(IJ) = IPERM(J)
     IPERM(J) = LM
     LM = IPERM(IJ)
!
!        If first element of array is greater than LM, interchange
!        with LM
!
     if (X(IPERM(I))  >  X(LM)) THEN
        IPERM(IJ) = IPERM(I)
        IPERM(I) = LM
        LM = IPERM(IJ)
     ENDIF
  end if
  go to 60
   50 LMT = IPERM(L)
  IPERM(L) = IPERM(K)
  IPERM(K) = LMT
!
!     Find an element in the second half of the array which is smaller
!     than LM
!
   60 L = L-1
  if (X(IPERM(L))  >  X(LM)) go to 60
!
!     Find an element in the first half of the array which is greater
!     than LM
!
   70 K = K+1
  if (X(IPERM(K))  <  X(LM)) go to 70
!
!     Interchange these elements
!
  if (K  <=  L) go to 50
!
!     Save upper and lower subscripts of the array yet to be sorted
!
  if (L-I  >  J-K) THEN
     IL(M) = I
     IU(M) = L
     I = K
     M = M+1
  ELSE
     IL(M) = K
     IU(M) = J
     J = L
     M = M+1
  end if
  go to 90
!
!     Begin again on another portion of the unsorted array
!
   80 M = M-1
  if (M  ==  0) go to 120
  I = IL(M)
  J = IU(M)
!
   90 if (J-I  >=  1) go to 40
  if (I  ==  1) go to 30
  I = I-1
!
  100 I = I+1
  if (I  ==  J) go to 80
  LM = IPERM(I+1)
  if (X(IPERM(I))  <=  X(LM)) go to 100
  K = I
!
  110 IPERM(K+1) = IPERM(K)
  K = K-1
!
  if (X(LM)  <  X(IPERM(K))) go to 110
  IPERM(K+1) = LM
  go to 100
!
!     Clean up
!
  120 if (KFLAG  <=  -1) THEN
     DO 130 I=1,NN
        X(I) = -X(I)
  130    CONTINUE
  end if
!
!     Rearrange the values of X if desired
!
  if (KK  ==  2) THEN
!
!        Use the IPERM vector as a flag.
!        If IPERM(I) < 0, then the I-th value is in correct location
!
     DO 150 ISTRT=1,NN
        if (IPERM(ISTRT)  >=  0) THEN
           INDX = ISTRT
           INDX0 = INDX
           TEMP = X(ISTRT)
  140          if (IPERM(INDX)  >  0) THEN
              X(INDX) = X(IPERM(INDX))
              INDX0 = INDX
              IPERM(INDX) = -IPERM(INDX)
              INDX = ABS(IPERM(INDX))
              go to 140
           ENDIF
           X(INDX0) = TEMP
        ENDIF
  150    CONTINUE
!
!        Revert the signs of the IPERM values
!
     DO 160 I=1,NN
        IPERM(I) = -IPERM(I)
  160    CONTINUE
!
  end if
!
  return
end