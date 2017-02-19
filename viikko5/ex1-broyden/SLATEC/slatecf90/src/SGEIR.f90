subroutine SGEIR (A, LDA, N, V, ITASK, IND, WORK, IWORK)
!
!! SGEIR solves a general system of linear equations.  Iterative refinement ...
!  is used to obtain an error estimate.
!
!***LIBRARY   SLATEC
!***CATEGORY  D2A1
!***TYPE      SINGLE PRECISION (SGEIR-S, CGEIR-C)
!***KEYWORDS  COMPLEX LINEAR EQUATIONS, GENERAL MATRIX,
!             GENERAL SYSTEM OF LINEAR EQUATIONS
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!    Subroutine SGEIR solves a general NxN system of single
!    precision linear equations using LINPACK subroutines SGEFA and
!    SGESL.  One pass of iterative refinement is used only to obtain
!    an estimate of the accuracy.  That is, if A is an NxN real
!    matrix and if X and B are real N-vectors, then SGEIR solves
!    the equation
!
!                          A*X=B.
!
!    The matrix A is first factored into upper and lower tri-
!    angular matrices U and L using partial pivoting.  These
!    factors and the pivoting information are used to calculate
!    the solution, X.  Then the residual vector is found and
!    used to calculate an estimate of the relative error, IND.
!    IND estimates the accuracy of the solution only when the
!    input matrix and the right hand side are represented
!    exactly in the computer and does not take into account
!    any errors in the input data.
!
!    If the equation A*X=B is to be solved for more than one vector
!    B, the factoring of A does not need to be performed again and
!    the option to solve only (ITASK  >  1) will be faster for
!    the succeeding solutions.  In this case, the contents of A,
!    LDA, N, WORK, and IWORK must not have been altered by the
!    user following factorization (ITASK=1).  IND will not be
!    changed by SGEIR in this case.
!
!  Argument Description ***
!
!    A      REAL(LDA,N)
!             the doubly subscripted array with dimension (LDA,N)
!             which contains the coefficient matrix.  A is not
!             altered by the routine.
!    LDA    INTEGER
!             the leading dimension of the array A.  LDA must be great-
!             er than or equal to N.  (terminal error message IND=-1)
!    N      INTEGER
!             the order of the matrix A.  The first N elements of
!             the array A are the elements of the first column of
!             matrix A.  N must be greater than or equal to 1.
!             (terminal error message IND=-2)
!    V      REAL(N)
!             on entry, the singly subscripted array(vector) of di-
!               mension N which contains the right hand side B of a
!               system of simultaneous linear equations A*X=B.
!             on return, V contains the solution vector, X .
!    ITASK  INTEGER
!             If ITASK=1, the matrix A is factored and then the
!               linear equation is solved.
!             If ITASK  >  1, the equation is solved using the existing
!               factored matrix A (stored in WORK).
!             If ITASK  <  1, then terminal error message IND=-3 is
!               printed.
!    IND    INTEGER
!             GT. 0  IND is a rough estimate of the number of digits
!                     of accuracy in the solution, X.  IND=75 means
!                     that the solution vector X is zero.
!             LT. 0  see error message corresponding to IND below.
!    WORK   REAL(N*(N+1))
!             a singly subscripted array of dimension at least N*(N+1).
!    IWORK  INTEGER(N)
!             a singly subscripted array of dimension at least N.
!
!  Error Messages Printed ***
!
!    IND=-1  terminal   N is greater than LDA.
!    IND=-2  terminal   N is less than one.
!    IND=-3  terminal   ITASK is less than one.
!    IND=-4  terminal   The matrix A is computationally singular.
!                         A solution has not been computed.
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
!***ROUTINES CALLED  R1MACH, SASUM, SCOPY, SDSDOT, SGEFA, SGESL, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800430  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SGEIR
!
  INTEGER LDA,N,ITASK,IND,IWORK(*),INFO,J
  REAL A(LDA,*),V(*),WORK(N,*),XNORM,DNORM,SDSDOT,SASUM,R1MACH
  CHARACTER*8 XERN1, XERN2
!***FIRST EXECUTABLE STATEMENT  SGEIR
  if (LDA < N) THEN
     IND = -1
     WRITE (XERN1, '(I8)') LDA
     WRITE (XERN2, '(I8)') N
     call XERMSG ('SLATEC', 'SGEIR', 'LDA = ' // XERN1 // &
        ' IS LESS THAN N = ' // XERN2, -1, 1)
     return
  end if
!
  if (N <= 0) THEN
     IND = -2
     WRITE (XERN1, '(I8)') N
     call XERMSG ('SLATEC', 'SGEIR', 'N = ' // XERN1 // &
        ' IS LESS THAN 1', -2, 1)
     return
  end if
!
  if (ITASK < 1) THEN
     IND = -3
     WRITE (XERN1, '(I8)') ITASK
     call XERMSG ('SLATEC', 'SGEIR', 'ITASK = ' // XERN1 // &
        ' IS LESS THAN 1', -3, 1)
     return
  end if
!
  if (ITASK == 1) THEN
!
!        MOVE MATRIX A TO WORK
!
     DO 10 J=1,N
        call SCOPY(N,A(1,J),1,WORK(1,J),1)
   10    CONTINUE
!
!        FACTOR MATRIX A INTO LU
!
     call SGEFA(WORK,N,N,IWORK,INFO)
!
!        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
!
     if (INFO /= 0) THEN
        IND = -4
        call XERMSG ('SLATEC', 'SGEIR', &
           'SINGULAR MATRIX A - NO SOLUTION', -4, 1)
        return
     ENDIF
  end if
!
!     SOLVE WHEN FACTORING COMPLETE
!     MOVE VECTOR B TO WORK
!
  call SCOPY(N,V(1),1,WORK(1,N+1),1)
  call SGESL(WORK,N,N,IWORK,V,0)
!
!     FORM NORM OF X0
!
  XNORM=SASUM(N,V(1),1)
  if (XNORM == 0.0) THEN
     IND = 75
     return
  end if
!
!     COMPUTE  RESIDUAL
!
  DO 40 J=1,N
     WORK(J,N+1) = SDSDOT(N,-WORK(J,N+1),A(J,1),LDA,V,1)
   40 CONTINUE
!
!     SOLVE A*DELTA=R
!
  call SGESL(WORK,N,N,IWORK,WORK(1,N+1),0)
!
!     FORM NORM OF DELTA
!
  DNORM = SASUM(N,WORK(1,N+1),1)
!
!     COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
!     AND CHECK FOR IND GREATER THAN ZERO
!
  IND = -LOG10(MAX(R1MACH(4),DNORM/XNORM))
  if (IND <= 0) THEN
     IND = -10
     call XERMSG ('SLATEC', 'SGEIR', &
        'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
  end if
  return
end
