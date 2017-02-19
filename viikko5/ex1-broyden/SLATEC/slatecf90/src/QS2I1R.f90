subroutine QS2I1R (IA, JA, A, N, KFLAG)
!
!! QS2I1R sorts an integer array, and adjusts two companion arrays.
!
!***SUBSIDIARY
!***PURPOSE  Sort an integer array, moving an integer and real array.
!            This routine sorts the integer array IA and makes the same
!            interchanges in the integer array JA and the real array A.
!            The array IA may be sorted in increasing order or decreas-
!            ing order.  A slightly modified QUICKSORT algorithm is
!            used.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  N6A2A
!***TYPE      SINGLE PRECISION (QS2I1R-S, QS2I1D-D)
!***KEYWORDS  SINGLETON QUICKSORT, SLAP, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Kahaner, D. K., (NBS)
!           Seager, M. K., (LLNL) seager@llnl.gov
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!     Written by Rondall E Jones
!     Modified by John A. Wisniewski to use the Singleton QUICKSORT
!     algorithm. date 18 November 1976.
!
!     Further modified by David K. Kahaner
!     National Bureau of Standards
!     August, 1981
!
!     Even further modification made to bring the code up to the
!     Fortran 77 level and make it more readable and to carry
!     along one integer array and one real array during the sort by
!     Mark K. Seager
!     Lawrence Livermore National Laboratory
!     November, 1987
!     This routine was adapted from the ISORT routine.
!
!     ABSTRACT
!         This routine sorts an integer array IA and makes the same
!         interchanges in the integer array JA and the real array A.
!         The array IA may be sorted in increasing order or decreasing
!         order.  A slightly modified quicksort algorithm is used.
!
!     DESCRIPTION OF PARAMETERS
!        IA - Integer array of values to be sorted.
!        JA - Integer array to be carried along.
!         A - Real array to be carried along.
!         N - Number of values in integer array IA to be sorted.
!     KFLAG - Control parameter
!           = 1 means sort IA in INCREASING order.
!           =-1 means sort IA in DECREASING order.
!
!***SEE ALSO  SS2Y
!***REFERENCES  R. C. Singleton, Algorithm 347, An Efficient Algorithm
!                 for Sorting With Minimal Storage, Communications ACM
!                 12:3 (1969), pp.185-7.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761118  DATE WRITTEN
!   890125  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   900805  Changed XERROR calls to calls to XERMSG.  (RWC)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910506  Made subsidiary to SS2Y and corrected reference.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!   920929  Corrected format of reference.  (FNF)
!   921012  Added E0's to f.p. constants.  (FNF)
!***END PROLOGUE  QS2I1R
!VD$R NOVECTOR
!VD$R NOCONCUR
!     .. Scalar Arguments ..
  INTEGER KFLAG, N
!     .. Array Arguments ..
  REAL A(N)
  INTEGER IA(N), JA(N)
!     .. Local Scalars ..
  REAL R, TA, TTA
  INTEGER I, IIT, IJ, IT, J, JJT, JT, K, KK, L, M, NN
!     .. Local Arrays ..
  INTEGER IL(21), IU(21)
!     .. External Subroutines ..
  EXTERNAL XERMSG
!     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  QS2I1R
  NN = N
  if (NN < 1) THEN
     call XERMSG ('SLATEC', 'QS2I1R', &
        'The number of values to be sorted was not positive.', 1, 1)
     return
  end if
  if (  N == 1 ) RETURN
  KK = ABS(KFLAG)
  if ( KK /= 1 ) THEN
     call XERMSG ('SLATEC', 'QS2I1R', &
        'The sort control parameter, K, was not 1 or -1.', 2, 1)
     return
  end if
!
!     Alter array IA to get decreasing order if needed.
!
  if (  KFLAG < 1 ) THEN
     DO 20 I=1,NN
        IA(I) = -IA(I)
 20      CONTINUE
  end if
!
!     Sort IA and carry JA and A along.
!     And now...Just a little black magic...
  M = 1
  I = 1
  J = NN
  R = .375E0
 210  if (  R <= 0.5898437E0 ) THEN
     R = R + 3.90625E-2
  ELSE
     R = R-.21875E0
  end if
 225  K = I
!
!     Select a central element of the array and save it in location
!     it, jt, at.
!
  IJ = I + INT ((J-I)*R)
  IT = IA(IJ)
  JT = JA(IJ)
  TA = A(IJ)
!
!     If first element of array is greater than it, interchange with it.
!
  if (  IA(I) > IT ) THEN
     IA(IJ) = IA(I)
     IA(I)  = IT
     IT     = IA(IJ)
     JA(IJ) = JA(I)
     JA(I)  = JT
     JT     = JA(IJ)
     A(IJ)  = A(I)
     A(I)   = TA
     TA     = A(IJ)
  end if
  L=J
!
!     If last element of array is less than it, swap with it.
!
  if (  IA(J) < IT ) THEN
     IA(IJ) = IA(J)
     IA(J)  = IT
     IT     = IA(IJ)
     JA(IJ) = JA(J)
     JA(J)  = JT
     JT     = JA(IJ)
     A(IJ)  = A(J)
     A(J)   = TA
     TA     = A(IJ)
!
!     If first element of array is greater than it, swap with it.
!
     if ( IA(I) > IT ) THEN
        IA(IJ) = IA(I)
        IA(I)  = IT
        IT     = IA(IJ)
        JA(IJ) = JA(I)
        JA(I)  = JT
        JT     = JA(IJ)
        A(IJ)  = A(I)
        A(I)   = TA
        TA     = A(IJ)
     ENDIF
  end if
!
!     Find an element in the second half of the array which is
!     smaller than it.
!
  240 L=L-1
  if (  IA(L) > IT ) go to 240
!
!     Find an element in the first half of the array which is
!     greater than it.
!
  245 K=K+1
  if (  IA(K) < IT ) go to 245
!
!     Interchange these elements.
!
  if (  K <= L ) THEN
     IIT   = IA(L)
     IA(L) = IA(K)
     IA(K) = IIT
     JJT   = JA(L)
     JA(L) = JA(K)
     JA(K) = JJT
     TTA   = A(L)
     A(L)  = A(K)
     A(K)  = TTA
     GOTO 240
  end if
!
!     Save upper and lower subscripts of the array yet to be sorted.
!
  if (  L-I > J-K ) THEN
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
  go to 260
!
!     Begin again on another portion of the unsorted array.
!
  255 M = M-1
  if (  M == 0 ) go to 300
  I = IL(M)
  J = IU(M)
  260 if (  J-I >= 1 ) go to 225
  if (  I == J ) go to 255
  if (  I == 1 ) go to 210
  I = I-1
  265 I = I+1
  if (  I == J ) go to 255
  IT = IA(I+1)
  JT = JA(I+1)
  TA =  A(I+1)
  if (  IA(I) <= IT ) go to 265
  K=I
  270 IA(K+1) = IA(K)
  JA(K+1) = JA(K)
  A(K+1)  =  A(K)
  K = K-1
  if (  IT < IA(K) ) go to 270
  IA(K+1) = IT
  JA(K+1) = JT
  A(K+1)  = TA
  go to 265
!
!     Clean up, if necessary.
!
  300 if (  KFLAG < 1 ) THEN
     DO 310 I=1,NN
        IA(I) = -IA(I)
 310     CONTINUE
  end if
  return
end
