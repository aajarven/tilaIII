subroutine CPOFS (A, LDA, N, V, ITASK, IND, WORK)
!
!! CPOFS solves a positive definite symmetric complex linear system.
!
!***LIBRARY   SLATEC
!***CATEGORY  D2D1B
!***TYPE      COMPLEX (SPOFS-S, DPOFS-D, CPOFS-C)
!***KEYWORDS  HERMITIAN, LINEAR EQUATIONS, POSITIVE DEFINITE, SYMMETRIC
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!    Subroutine CPOFS solves a  positive definite symmetric
!    NxN system of complex linear equations using LINPACK
!    subroutines CPOCO and CPOSL.  That is, if A is an NxN
!    complex positive definite symmetric matrix and if X and B
!    are complex N-vectors, then CPOFS solves the equation
!
!                          A*X=B.
!
!    Care should be taken not to use CPOFS with a non-Hermitian
!    matrix.
!
!    The matrix A is first factored into upper and lower tri-
!    angular matrices R and R-TRANSPOSE.  These factors are used to
!    find the solution vector X.  An approximate condition number is
!    calculated to provide a rough estimate of the number of
!    digits of accuracy in the computed solution.
!
!    If the equation A*X=B is to be solved for more than one vector
!    B, the factoring of a does not need to be performed again and
!    the option to only solve (ITASK  >  1) will be faster for
!    the succeeding solutions.  In this case, the contents of A,
!    LDA, and N must not have been altered by the user following
!    factorization (ITASK=1).  IND will not be changed by CPOFS
!    in this case.
!
!  Argument Description ***
!
!    A      COMPLEX(LDA,N)
!             on entry, the doubly subscripted array with dimension
!               (LDA,N) which contains the coefficient matrix.  Only
!               the upper triangle, including the diagonal, of the
!               coefficient matrix need be entered and will subse-
!               quently be referenced and changed by the routine.
!             on return, contains in its upper triangle an upper
!               triangular matrix R such that  A = (R-TRANSPOSE) * R .
!    LDA    INTEGER
!             the leading dimension of the array A.  LDA must be great-
!             er than or equal to N.  (terminal error message IND=-1)
!    N      INTEGER
!             the order of the matrix A.  N must be greater
!             than or equal to 1.  (terminal error message IND=-2)
!    V      COMPLEX(N)
!             on entry the singly subscripted array(vector) of di-
!               mension N which contains the right hand side B of a
!               system of simultaneous linear equations A*X=B.
!             on return, V contains the solution vector, X .
!    ITASK  INTEGER
!             if ITASK = 1, the matrix A is factored and then the
!               linear equation is solved.
!             if ITASK  >  1, the equation is solved using the existing
!               factored matrix A.
!             if ITASK  <  1, then terminal error message IND=-3 is
!               printed.
!    IND    INTEGER
!             GT. 0  IND is a rough estimate of the number of digits
!                     of accuracy in the solution, X.
!             LT. 0  see error message corresponding to IND below.
!    WORK   COMPLEX(N)
!             a singly subscripted array of dimension at least N.
!
!  Error Messages Printed ***
!
!    IND=-1  terminal   N is greater than LDA.
!    IND=-2  terminal   N is less than 1.
!    IND=-3  terminal   ITASK is less than 1.
!    IND=-4  terminal   The matrix A is computationally singular or
!                         is not positive definite.  A solution
!                         has not been computed.
!    IND=-10 warning    The solution has no apparent significance.
!                         The solution may be inaccurate or the
!                         matrix A may be poorly scaled.
!
!               NOTE-  The above terminal(*fatal*) error messages are
!                      designed to be handled by XERMSG in which
!                      LEVEL=1 (recoverable) and IFLAG=2 .  LEVEL=0
!                      for warning error messages from XERMSG.  Unless
!                      the user provides otherwise, an error message
!                      will be printed followed by an abort.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CPOCO, CPOSL, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800516  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Convert XERRWV calls to XERMSG calls, cvt GOTO's to
!           IF-THEN-ELSE.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPOFS
!
  INTEGER LDA,N,ITASK,IND,INFO
  COMPLEX A(LDA,*),V(*),WORK(*)
  REAL R1MACH
  REAL RCOND
  CHARACTER*8 XERN1, XERN2
!***FIRST EXECUTABLE STATEMENT  CPOFS
  if (LDA < N) THEN
     IND = -1
     WRITE (XERN1, '(I8)') LDA
     WRITE (XERN2, '(I8)') N
     call XERMSG ('SLATEC', 'CPOFS', 'LDA = ' // XERN1 // &
        ' IS LESS THAN N = ' // XERN2, -1, 1)
     return
  end if
!
  if (N <= 0) THEN
     IND = -2
     WRITE (XERN1, '(I8)') N
     call XERMSG ('SLATEC', 'CPOFS', 'N = ' // XERN1 // &
        ' IS LESS THAN 1', -2, 1)
     return
  end if
!
  if (ITASK < 1) THEN
     IND = -3
     WRITE (XERN1, '(I8)') ITASK
     call XERMSG ('SLATEC', 'CPOFS', 'ITASK = ' // XERN1 // &
        ' IS LESS THAN 1', -3, 1)
     return
  end if
!
  if (ITASK == 1) THEN
!
!        FACTOR MATRIX A INTO R
!
     call CPOCO(A,LDA,N,RCOND,WORK,INFO)
!
!        CHECK FOR POSITIVE DEFINITE MATRIX
!
     if (INFO /= 0) THEN
        IND = -4
        call XERMSG ('SLATEC', 'CPOFS', &
           'SINGULAR OR NOT POSITIVE DEFINITE - NO SOLUTION', -4, 1)
        return
     ENDIF
!
!        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
!        AND CHECK FOR IND GREATER THAN ZERO
!
     IND = -LOG10(R1MACH(4)/RCOND)
     if (IND <= 0) THEN
        IND = -10
        call XERMSG ('SLATEC', 'CPOFS', &
           'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
     ENDIF
  end if
!
!     SOLVE AFTER FACTORING
!
  call CPOSL(A,LDA,N,V)
  return
end
