subroutine DNBFS (ABE, LDA, N, ML, MU, V, ITASK, IND, WORK, IWORK)
!
!! DNBFS solves a general nonsymmetric banded system of linear equations.
!
!***LIBRARY   SLATEC
!***CATEGORY  D2A2
!***TYPE      DOUBLE PRECISION (SNBFS-S, DNBFS-D, CNBFS-C)
!***KEYWORDS  BANDED, LINEAR EQUATIONS, NONSYMMETRIC
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!    Subroutine DNBFS solves a general nonsymmetric banded NxN
!    system of double precision real linear equations using
!    SLATEC subroutines DNBCO and DNBSL.  These are adaptations
!    of the LINPACK subroutines DGBCO and DGBSL which require
!    a different format for storing the matrix elements.  If
!    A  is an NxN double precision matrix and if  X  and  B  are
!    double precision N-vectors, then DNBFS solves the equation
!
!                          A*X=B.
!
!    A band matrix is a matrix whose nonzero elements are all
!    fairly near the main diagonal, specifically  A(I,J) = 0
!    if  I-J is greater than  ML  or  J-I  is greater than
!    MU .  The integers ML and MU are called the lower and upper
!    band widths and  M = ML+MU+1  is the total band width.
!    DNBFS uses less time and storage than the corresponding
!    program for general matrices (DGEFS) if 2*ML+MU  <   N .
!
!    The matrix A is first factored into upper and lower tri-
!    angular matrices U and L using partial pivoting.  These
!    factors and the pivoting information are used to find the
!    solution vector X.  An approximate condition number is
!    calculated to provide a rough estimate of the number of
!    digits of accuracy in the computed solution.
!
!    If the equation A*X=B is to be solved for more than one vector
!    B, the factoring of A does not need to be performed again and
!    the option to only solve (ITASK  >  1) will be faster for
!    the succeeding solutions.  In this case, the contents of A,
!    LDA, N and IWORK must not have been altered by the user follow-
!    ing factorization (ITASK=1).  IND will not be changed by DNBFS
!    in this case.
!
!
!    Band Storage
!
!          If  A  is a band matrix, the following program segment
!          will set up the input.
!
!                  ML = (band width below the diagonal)
!                  MU = (band width above the diagonal)
!                  DO 20 I = 1, N
!                     J1 = MAX(1, I-ML)
!                     J2 = MIN(N, I+MU)
!                     DO 10 J = J1, J2
!                        K = J - I + ML + 1
!                        ABE(I,K) = A(I,J)
!               10    CONTINUE
!               20 CONTINUE
!
!          This uses columns  1  through  ML+MU+1  of ABE .
!          Furthermore,  ML  additional columns are needed in
!          ABE  starting with column  ML+MU+2  for elements
!          generated during the triangularization.  The total
!          number of columns needed in  ABE  is  2*ML+MU+1 .
!
!    Example:  If the original matrix is
!
!          111213  0  0  0
!          21222324  0  0
!           032333435  0
!           0  043444546
!           0  0  0545556
!           0  0  0  06566
!
!     then  N = 6, ML = 1, MU = 2, LDA  >=  5  and ABE should contain
!
!           * 111213  +     , * = not used
!          21222324  +     , + = used for pivoting
!          32333435  +
!          43444546  +
!          545556  *  +
!          6566  *  *  +
!
!
!  Argument Description ***
!
!    ABE    DOUBLE PRECISION(LDA,NC)
!             on entry, contains the matrix in band storage as
!               described above.  NC  must not be less than
!               2*ML+MU+1 .  The user is cautioned to specify  NC
!               with care since it is not an argument and cannot
!               be checked by DNBFS.  The rows of the original
!               matrix are stored in the rows of  ABE  and the
!               diagonals of the original matrix are stored in
!               columns  1  through  ML+MU+1  of  ABE .
!             on return, contains an upper triangular matrix U and
!               the multipliers necessary to construct a matrix L
!               so that A=L*U.
!    LDA    INTEGER
!             the leading dimension of array ABE.  LDA must be great-
!             er than or equal to N.  (terminal error message IND=-1)
!    N      INTEGER
!             the order of the matrix A.  N must be greater
!             than or equal to 1 .  (terminal error message IND=-2)
!    ML     INTEGER
!             the number of diagonals below the main diagonal.
!             ML  must not be less than zero nor greater than or
!             equal to  N .  (terminal error message IND=-5)
!    MU     INTEGER
!             the number of diagonals above the main diagonal.
!             MU  must not be less than zero nor greater than or
!             equal to  N .  (terminal error message IND=-6)
!    V      DOUBLE PRECISION(N)
!             on entry, the singly subscripted array(vector) of di-
!               mension N which contains the right hand side B of a
!               system of simultaneous linear equations A*X=B.
!             on return, V contains the solution vector, X .
!    ITASK  INTEGER
!             If ITASK=1, the matrix A is factored and then the
!               linear equation is solved.
!             If ITASK  >  1, the equation is solved using the existing
!               factored matrix A and IWORK.
!             If ITASK  <  1, then terminal error message IND=-3 is
!               printed.
!    IND    INTEGER
!             GT. 0  IND is a rough estimate of the number of digits
!                     of accuracy in the solution, X.
!             LT. 0  See error message corresponding to IND below.
!    WORK   DOUBLE PRECISION(N)
!             a singly subscripted array of dimension at least N.
!    IWORK  INTEGER(N)
!             a singly subscripted array of dimension at least N.
!
!  Error Messages Printed ***
!
!    IND=-1  terminal   N is greater than LDA.
!    IND=-2  terminal   N is less than 1.
!    IND=-3  terminal   ITASK is less than 1.
!    IND=-4  terminal   The matrix A is computationally singular.
!                         A solution has not been computed.
!    IND=-5  terminal   ML is less than zero or is greater than
!                         or equal to N .
!    IND=-6  terminal   MU is less than zero or is greater than
!                         or equal to N .
!    IND=-10 warning    The solution has no apparent significance.
!                         The solution may be inaccurate or the matrix
!                         A may be poorly scaled.
!
!               Note-  The above terminal(*fatal*) error messages are
!                      designed to be handled by XERMSG in which
!                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
!                      for warning error messages from XERMSG.  Unless
!                      the user provides otherwise, an error message
!                      will be printed followed by an abort.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  D1MACH, DNBCO, DNBSL, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800812  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Convert XERRWV calls to XERMSG calls, changed GOTOs to
!           IF-THEN-ELSEs.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DNBFS
!
  INTEGER LDA,N,ITASK,IND,IWORK(*),ML,MU
  DOUBLE PRECISION ABE(LDA,*),V(*),WORK(*),D1MACH
  DOUBLE PRECISION RCOND
  CHARACTER*8 XERN1, XERN2
!***FIRST EXECUTABLE STATEMENT  DNBFS
  if (LDA < N) THEN
     IND = -1
     WRITE (XERN1, '(I8)') LDA
     WRITE (XERN2, '(I8)') N
     call XERMSG ('SLATEC', 'DNBFS', 'LDA = ' // XERN1 // &
        ' IS LESS THAN N = ' // XERN2, -1, 1)
     return
  end if
!
  if (N <= 0) THEN
     IND = -2
     WRITE (XERN1, '(I8)') N
     call XERMSG ('SLATEC', 'DNBFS', 'N = ' // XERN1 // &
        ' IS LESS THAN 1', -2, 1)
     return
  end if
!
  if (ITASK < 1) THEN
     IND = -3
     WRITE (XERN1, '(I8)') ITASK
     call XERMSG ('SLATEC', 'DNBFS', 'ITASK = ' // XERN1 // &
        ' IS LESS THAN 1', -3, 1)
     return
  end if
!
  if (ML < 0 .OR. ML >= N) THEN
     IND = -5
     WRITE (XERN1, '(I8)') ML
     call XERMSG ('SLATEC', 'DNBFS', &
        'ML = ' // XERN1 // ' IS OUT OF RANGE', -5, 1)
     return
  end if
!
  if (MU < 0 .OR. MU >= N) THEN
     IND = -6
     WRITE (XERN1, '(I8)') MU
     call XERMSG ('SLATEC', 'DNBFS', &
        'MU = ' // XERN1 // ' IS OUT OF RANGE', -6, 1)
     return
  end if
!
  if (ITASK == 1) THEN
!
!        FACTOR MATRIX A INTO LU
!
     call DNBCO(ABE,LDA,N,ML,MU,IWORK,RCOND,WORK)
!
!        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
!
     if (RCOND == 0.0D0) THEN
        IND = -4
        call XERMSG ('SLATEC', 'DNBFS', &
           'SINGULAR MATRIX A - NO SOLUTION', -4, 1)
        return
     ENDIF
!
!        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
!        AND CHECK FOR IND GREATER THAN ZERO
!
     IND = -LOG10(D1MACH(4)/RCOND)
     if (IND <= 0) THEN
        IND = -10
        call XERMSG ('SLATEC', 'DNBFS', &
           'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
     ENDIF
  end if
!
!     SOLVE AFTER FACTORING
!
  call DNBSL(ABE,LDA,N,ML,MU,IWORK,V,0)
  return
end
