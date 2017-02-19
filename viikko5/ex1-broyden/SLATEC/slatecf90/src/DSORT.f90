subroutine DSORT (DX, DY, N, KFLAG)
!
!! DSORT sorts an array and optionally make the same interchanges in ...
!  an auxiliary array.
!
!  The array may be sorted in increasing
!            or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!
!***LIBRARY   SLATEC
!***CATEGORY  N6A2B
!***TYPE      DOUBLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
!***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!
!   DSORT sorts array DX and optionally makes the same interchanges in
!   array DY.  The array DX may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      DX - array of values to be sorted   (usually abscissas)
!      DY - array to be (optionally) carried along
!      N  - number of values in array DX to be sorted
!      KFLAG - control parameter
!            =  2  means sort DX in increasing order and carry DY along.
!            =  1  means sort DX in increasing order (ignoring DY)
!            = -1  means sort DX in decreasing order (ignoring DY)
!            = -2  means sort DX in decreasing order and carry DY along.
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
!   901012  Declared all variables; changed X,Y to DX,DY; changed
!           code to parallel SSORT. (M. McClain)
!   920501  Reformatted the REFERENCES section.  (DWL, WRB)
!   920519  Clarified error messages.  (DWL)
!   920801  Declarations section rebuilt and code restructured to use
!           IF-THEN-ELSE-ENDIF.  (RWC, WRB)
!***END PROLOGUE  DSORT
!     .. Scalar Arguments ..
  INTEGER KFLAG, N
!     .. Array Arguments ..
  DOUBLE PRECISION DX(*), DY(*)
!     .. Local Scalars ..
  DOUBLE PRECISION R, T, TT, TTY, TY
  INTEGER I, IJ, J, K, KK, L, M, NN
!     .. Local Arrays ..
  INTEGER IL(21), IU(21)
!     .. External Subroutines ..
  EXTERNAL XERMSG
!     .. Intrinsic Functions ..
  INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  DSORT
  NN = N
  if (NN  <  1) THEN
     call XERMSG ('SLATEC', 'DSORT', &
        'The number of values to be sorted is not positive.', 1, 1)
     return
  end if
!
  KK = ABS(KFLAG)
  if (KK /= 1 .AND. KK /= 2) THEN
     call XERMSG ('SLATEC', 'DSORT', &
        'The sort control parameter, K, is not 2, 1, -1, or -2.', 2, &
        1)
     return
  end if
!
!     Alter array DX to get decreasing order if needed
!
  if (KFLAG  <=  -1) THEN
     DO 10 I=1,NN
        DX(I) = -DX(I)
   10    CONTINUE
  end if
!
  if (KK  ==  2) go to 100
!
!     Sort DX only
!
  M = 1
  I = 1
  J = NN
  R = 0.375D0
!
   20 if (I  ==  J) go to 60
  if (R  <=  0.5898437D0) THEN
     R = R+3.90625D-2
  ELSE
     R = R-0.21875D0
  end if
!
   30 K = I
!
!     Select a central element of the array and save it in location T
!
  IJ = I + INT((J-I)*R)
  T = DX(IJ)
!
!     If first element of array is greater than T, interchange with T
!
  if (DX(I)  >  T) THEN
     DX(IJ) = DX(I)
     DX(I) = T
     T = DX(IJ)
  end if
  L = J
!
!     If last element of array is less than than T, interchange with T
!
  if (DX(J)  <  T) THEN
     DX(IJ) = DX(J)
     DX(J) = T
     T = DX(IJ)
!
!        If first element of array is greater than T, interchange with T
!
     if (DX(I)  >  T) THEN
        DX(IJ) = DX(I)
        DX(I) = T
        T = DX(IJ)
     ENDIF
  end if
!
!     Find an element in the second half of the array which is smaller
!     than T
!
   40 L = L-1
  if (DX(L)  >  T) go to 40
!
!     Find an element in the first half of the array which is greater
!     than T
!
   50 K = K+1
  if (DX(K)  <  T) go to 50
!
!     Interchange these elements
!
  if (K  <=  L) THEN
     TT = DX(L)
     DX(L) = DX(K)
     DX(K) = TT
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
  T = DX(I+1)
  if (DX(I)  <=  T) go to 80
  K = I
!
   90 DX(K+1) = DX(K)
  K = K-1
  if (T  <  DX(K)) go to 90
  DX(K+1) = T
  go to 80
!
!     Sort DX and carry DY along
!
  100 M = 1
  I = 1
  J = NN
  R = 0.375D0
!
  110 if (I  ==  J) go to 150
  if (R  <=  0.5898437D0) THEN
     R = R+3.90625D-2
  ELSE
     R = R-0.21875D0
  end if
!
  120 K = I
!
!     Select a central element of the array and save it in location T
!
  IJ = I + INT((J-I)*R)
  T = DX(IJ)
  TY = DY(IJ)
!
!     If first element of array is greater than T, interchange with T
!
  if (DX(I)  >  T) THEN
     DX(IJ) = DX(I)
     DX(I) = T
     T = DX(IJ)
     DY(IJ) = DY(I)
     DY(I) = TY
     TY = DY(IJ)
  end if
  L = J
!
!     If last element of array is less than T, interchange with T
!
  if (DX(J)  <  T) THEN
     DX(IJ) = DX(J)
     DX(J) = T
     T = DX(IJ)
     DY(IJ) = DY(J)
     DY(J) = TY
     TY = DY(IJ)
!
!        If first element of array is greater than T, interchange with T
!
     if (DX(I)  >  T) THEN
        DX(IJ) = DX(I)
        DX(I) = T
        T = DX(IJ)
        DY(IJ) = DY(I)
        DY(I) = TY
        TY = DY(IJ)
     ENDIF
  end if
!
!     Find an element in the second half of the array which is smaller
!     than T
!
  130 L = L-1
  if (DX(L)  >  T) go to 130
!
!     Find an element in the first half of the array which is greater
!     than T
!
  140 K = K+1
  if (DX(K)  <  T) go to 140
!
!     Interchange these elements
!
  if (K  <=  L) THEN
     TT = DX(L)
     DX(L) = DX(K)
     DX(K) = TT
     TTY = DY(L)
     DY(L) = DY(K)
     DY(K) = TTY
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
  T = DX(I+1)
  TY = DY(I+1)
  if (DX(I)  <=  T) go to 170
  K = I
!
  180 DX(K+1) = DX(K)
  DY(K+1) = DY(K)
  K = K-1
  if (T  <  DX(K)) go to 180
  DX(K+1) = T
  DY(K+1) = TY
  go to 170
!
!     Clean up
!
  190 if (KFLAG  <=  -1) THEN
     DO 200 I=1,NN
        DX(I) = -DX(I)
  200    CONTINUE
  end if
  return
end
