subroutine CPOIR (A, LDA, N, V, ITASK, IND, WORK)
!
!! CPOIR solves a positive definite Hermitian system of linear equations. ...
! Iterative refinement is used to obtain an error estimate.
!
!***LIBRARY   SLATEC
!***CATEGORY  D2D1B
!***TYPE      COMPLEX (SPOIR-S, CPOIR-C)
!***KEYWORDS  HERMITIAN, LINEAR EQUATIONS, POSITIVE DEFINITE, SYMMETRIC
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!    Subroutine CPOIR solves a complex positive definite Hermitian
!    NxN system of single precision linear equations using LINPACK
!    subroutines CPOFA and CPOSL.  One pass of iterative refine-
!    ment is used only to obtain an estimate of the accuracy.  That
!    is, if A is an NxN complex positive definite Hermitian matrix
!    and if X and B are complex N-vectors, then CPOIR solves the
!    equation
!
!                          A*X=B.
!
!    Care should be taken not to use CPOIR with a non-Hermitian
!    matrix.
!
!    The matrix A is first factored into upper and lower
!    triangular matrices R and R-TRANSPOSE.  These
!    factors are used to calculate the solution, X.
!    Then the residual vector is found and used
!    to calculate an estimate of the relative error, IND.
!    IND estimates the accuracy of the solution only when the
!    input matrix and the right hand side are represented
!    exactly in the computer and does not take into account
!    any errors in the input data.
!
!    If the equation A*X=B is to be solved for more than one vector
!    B, the factoring of A does not need to be performed again and
!    the option to only solve (ITASK  >  1) will be faster for
!    the succeeding solutions.  In this case, the contents of A,
!    LDA, N, and WORK must not have been altered by the user
!    following factorization (ITASK=1).  IND will not be changed
!    by CPOIR in this case.
!
!  Argument Description ***
!    A       COMPLEX(LDA,N)
!             the doubly subscripted array with dimension (LDA,N)
!             which contains the coefficient matrix.  Only the
!             upper triangle, including the diagonal, of the
!             coefficient matrix need be entered.  A is not
!             altered by the routine.
!    LDA    INTEGER
!             the leading dimension of the array A.  LDA must be great-
!             er than or equal to N.  (terminal error message IND=-1)
!    N      INTEGER
!             the order of the matrix A.  N must be greater than
!             or equal to one.  (terminal error message IND=-2)
!    V      COMPLEX(N)
!             on entry, the singly subscripted array(vector) of di-
!               mension N which contains the right hand side B of a
!               system of simultaneous linear equations A*X=B.
!             on return, V contains the solution vector, X .
!    ITASK  INTEGER
!             if ITASK = 1, the matrix A is factored and then the
!               linear equation is solved.
!             if ITASK  >  1, the equation is solved using the existing
!               factored matrix A (stored in WORK).
!             if ITASK  <  1, then terminal terminal error IND=-3 is
!               printed.
!    IND    INTEGER
!             GT. 0  IND is a rough estimate of the number of digits
!                     of accuracy in the solution, X.  IND=75 means
!                     that the solution vector X is zero.
!             LT. 0  see error message corresponding to IND below.
!    WORK   COMPLEX(N*(N+1))
!             a singly subscripted array of dimension at least N*(N+1).
!
!  Error Messages Printed ***
!
!    IND=-1  terminal   N is greater than LDA.
!    IND=-2  terminal   N is less than one.
!    IND=-3  terminal   ITASK is less than one.
!    IND=-4  terminal   The matrix A is computationally singular
!                         or is not positive definite.
!                         A solution has not been computed.
!    IND=-10 warning    The solution has no apparent significance.
!                         the solution may be inaccurate or the matrix
!                         a may be poorly scaled.
!
!               NOTE-  the above terminal(*fatal*) error messages are
!                      designed to be handled by XERMSG in which
!                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
!                      for warning error messages from XERMSG.  Unless
!                      the user provides otherwise, an error message
!                      will be printed followed by an abort.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CCOPY, CPOFA, CPOSL, DCDOT, R1MACH, SCASUM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800530  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Convert XERRWV calls to XERMSG calls, cvt GOTO's to
!           IF-THEN-ELSE.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPOIR
!
  INTEGER LDA,N,ITASK,IND,INFO,J
  COMPLEX A(LDA,*),V(*),WORK(N,*)
  REAL SCASUM,XNORM,DNORM,R1MACH
  DOUBLE PRECISION DR1,DI1,DR2,DI2
  CHARACTER*8 XERN1, XERN2
!***FIRST EXECUTABLE STATEMENT  CPOIR
  if (LDA < N) THEN
     IND = -1
     WRITE (XERN1, '(I8)') LDA
     WRITE (XERN2, '(I8)') N
     call XERMSG ('SLATEC', 'CPOIR', 'LDA = ' // XERN1 // &
        ' IS LESS THAN N = ' // XERN2, -1, 1)
     return
  end if
!
  if (N <= 0) THEN
     IND = -2
     WRITE (XERN1, '(I8)') N
     call XERMSG ('SLATEC', 'CPOIR', 'N = ' // XERN1 // &
        ' IS LESS THAN 1', -2, 1)
     return
  end if
!
  if (ITASK < 1) THEN
     IND = -3
     WRITE (XERN1, '(I8)') ITASK
     call XERMSG ('SLATEC', 'CPOIR', 'ITASK = ' // XERN1 // &
        ' IS LESS THAN 1', -3, 1)
     return
  end if
!
  if (ITASK == 1) THEN
!
!        MOVE MATRIX A TO WORK
!
     DO 10 J=1,N
        call CCOPY(N,A(1,J),1,WORK(1,J),1)
   10    CONTINUE
!
!        FACTOR MATRIX A INTO R
!
     call CPOFA(WORK,N,N,INFO)
!
!        CHECK FOR  SINGULAR OR NOT POS.DEF. MATRIX
!
     if (INFO /= 0) THEN
        IND = -4
        call XERMSG ('SLATEC', 'CPOIR', &
           'SINGULAR OR NOT POSITIVE DEFINITE - NO SOLUTION', -4, 1)
        return
     ENDIF
  end if
!
!     SOLVE AFTER FACTORING
!     MOVE VECTOR B TO WORK
!
  call CCOPY(N,V(1),1,WORK(1,N+1),1)
  call CPOSL(WORK,N,N,V)
!
!     FORM NORM OF X0
!
  XNORM = SCASUM(N,V(1),1)
  if (XNORM == 0.0) THEN
     IND = 75
     return
  end if
!
!     COMPUTE  RESIDUAL
!
  DO 40 J=1,N
     call DCDOT(J-1,-1.D0,A(1,J),1,V(1),1,DR1,DI1)
     call DCDOT(N-J+1,1.D0,A(J,J),LDA,V(J),1,DR2,DI2)
     DR1 = DR1+DR2-DBLE(REAL(WORK(J,N+1)))
     DI1 = DI1+DI2-DBLE(AIMAG(WORK(J,N+1)))
     WORK(J,N+1) = CMPLX(REAL(DR1),REAL(DI1))
   40 CONTINUE
!
!     SOLVE A*DELTA=R
!
  call CPOSL(WORK,N,N,WORK(1,N+1))
!
!     FORM NORM OF DELTA
!
  DNORM = SCASUM(N,WORK(1,N+1),1)
!
!     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
!     AND CHECK FOR IND GREATER THAN ZERO
!
  IND = -LOG10(MAX(R1MACH(4),DNORM/XNORM))
  if (IND <= 0) THEN
     IND = -10
     call XERMSG ('SLATEC', 'CPOIR', &
        'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
  end if
  return
end
