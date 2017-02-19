subroutine SSORT (X, Y, N, KFLAG)
!
!! SSORT sorts an array and optionally make the same interchanges in ...
!            an auxiliary array.  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      Y - array to be (optionally) carried along
!      N - number of values in array X to be sorted
!      KFLAG - control parameter
!            =  2  means sort X in increasing order and carry Y along.
!            =  1  means sort X in increasing order (ignoring Y)
!            = -1  means sort X in decreasing order (ignoring Y)
!            = -2  means sort X in decreasing order and carry Y along.
!
!***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
!                 for sorting with minimal storage, Communications of
!                 the ACM, 12, 3 (1969), pp. 185-187.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761101  DATE WRITTEN
!   761118  Modified to use the Singleton quicksort algorithm.  (JAW)
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement labels.  (WRB)
!   891024  Changed category.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901012  Declared all variables; changed X,Y to SX,SY. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  SSORT
!     .. Scalar Arguments ..
  INTEGER KFLAG, N
!     .. Array Arguments ..
  REAL X(*), Y(*)
!     .. Local Scalars ..
  REAL R, T, TT, TTY, TY
  INTEGER I, IJ, J, K, KK, L, M, NN
!     .. Local Arrays ..
  INTEGER IL(21), IU(21)
!     .. External Subroutines ..
  EXTERNAL XERMSG
!     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  SSORT
  NN = N
  if (NN  <  1) THEN
     call XERMSG ('SLATEC', 'SSORT', &
        'The number of values to be sorted is not positive.', 1, 1)
     return
  end if
!
  KK = ABS(KFLAG)
  if (KK /= 1 .AND. KK /= 2) THEN
     call XERMSG ('SLATEC', 'SSORT', &
        'The sort control parameter, K, is not 2, 1, -1, or -2.', 2, &
        1)
     return
  end if
!
!     Alter array X to get decreasing order if needed
!
  if (KFLAG  <=  -1) THEN
     DO 10 I=1,NN
        X(I) = -X(I)
   10    CONTINUE
  end if
!
  if (KK  ==  2) go to 100
!
!     Sort X only
!
  M = 1
  I = 1
  J = NN
  R = 0.375E0
!
   20 if (I  ==  J) go to 60
  if (R  <=  0.5898437E0) THEN
     R = R+3.90625E-2
  ELSE
     R = R-0.21875E0
  end if
!
   30 K = I
!
!     Select a central element of the array and save it in location T
!
  IJ = I + INT((J-I)*R)
  T = X(IJ)
!
!     If first element of array is greater than T, interchange with T
!
  if (X(I)  >  T) THEN
     X(IJ) = X(I)
     X(I) = T
     T = X(IJ)
  end if
  L = J
!
!     If last element of array is less than than T, interchange with T
!
  if (X(J)  <  T) THEN
     X(IJ) = X(J)
     X(J) = T
     T = X(IJ)
!
!        If first element of array is greater than T, interchange with T
!
     if (X(I)  >  T) THEN
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
     ENDIF
  end if
!
!     Find an element in the second half of the array which is smaller
!     than T
!
   40 L = L-1
  if (X(L)  >  T) go to 40
!
!     Find an element in the first half of the array which is greater
!     than T
!
   50 K = K+1
  if (X(K)  <  T) go to 50
!
!     Interchange these elements
!
  if (K  <=  L) THEN
     TT = X(L)
     X(L) = X(K)
     X(K) = TT
     go to 40
  end if
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
  go to 70
!
!     Begin again on another portion of the unsorted array
!
   60 M = M-1
  if (M  ==  0) go to 190
  I = IL(M)
  J = IU(M)
!
   70 if (J-I  >=  1) go to 30
  if (I  ==  1) go to 20
  I = I-1
!
   80 I = I+1
  if (I  ==  J) go to 60
  T = X(I+1)
  if (X(I)  <=  T) go to 80
  K = I
!
   90 X(K+1) = X(K)
  K = K-1
  if (T  <  X(K)) go to 90
  X(K+1) = T
  go to 80
!
!     Sort X and carry Y along
!
  100 M = 1
  I = 1
  J = NN
  R = 0.375E0
!
  110 if (I  ==  J) go to 150
  if (R  <=  0.5898437E0) THEN
     R = R+3.90625E-2
  ELSE
     R = R-0.21875E0
  end if
!
  120 K = I
!
!     Select a central element of the array and save it in location T
!
  IJ = I + INT((J-I)*R)
  T = X(IJ)
  TY = Y(IJ)
!
!     If first element of array is greater than T, interchange with T
!
  if (X(I)  >  T) THEN
     X(IJ) = X(I)
     X(I) = T
     T = X(IJ)
     Y(IJ) = Y(I)
     Y(I) = TY
     TY = Y(IJ)
  end if
  L = J
!
!     If last element of array is less than T, interchange with T
!
  if (X(J)  <  T) THEN
     X(IJ) = X(J)
     X(J) = T
     T = X(IJ)
     Y(IJ) = Y(J)
     Y(J) = TY
     TY = Y(IJ)
!
!        If first element of array is greater than T, interchange with T
!
     if (X(I)  >  T) THEN
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
        Y(IJ) = Y(I)
        Y(I) = TY
        TY = Y(IJ)
     ENDIF
  end if
!
!     Find an element in the second half of the array which is smaller
!     than T
!
  130 L = L-1
  if (X(L)  >  T) go to 130
!
!     Find an element in the first half of the array which is greater
!     than T
!
  140 K = K+1
  if (X(K)  <  T) go to 140
!
!     Interchange these elements
!
  if (K  <=  L) THEN
     TT = X(L)
     X(L) = X(K)
     X(K) = TT
     TTY = Y(L)
     Y(L) = Y(K)
     Y(K) = TTY
     go to 130
  end if
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
  go to 160
!
!     Begin again on another portion of the unsorted array
!
  150 M = M-1
  if (M  ==  0) go to 190
  I = IL(M)
  J = IU(M)
!
  160 if (J-I  >=  1) go to 120
  if (I  ==  1) go to 110
  I = I-1
!
  170 I = I+1
  if (I  ==  J) go to 150
  T = X(I+1)
  TY = Y(I+1)
  if (X(I)  <=  T) go to 170
  K = I
!
  180 X(K+1) = X(K)
  Y(K+1) = Y(K)
  K = K-1
  if (T  <  X(K)) go to 180
  X(K+1) = T
  Y(K+1) = TY
  go to 170
!
!     Clean up
!
  190 if (KFLAG  <=  -1) THEN
     DO 200 I=1,NN
        X(I) = -X(I)
  200    CONTINUE
  end if
  return
end
