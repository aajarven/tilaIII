subroutine CGEFS (A, LDA, N, V, ITASK, IND, WORK, IWORK)
!
!! CGEFS solves a general system of linear equations.
!
!***LIBRARY   SLATEC
!***CATEGORY  D2C1
!***TYPE      COMPLEX (SGEFS-S, DGEFS-D, CGEFS-C)
!***KEYWORDS  COMPLEX LINEAR EQUATIONS, GENERAL MATRIX,
!             GENERAL SYSTEM OF LINEAR EQUATIONS
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!    Subroutine CGEFS solves A general NxN system of complex
!    linear equations using LINPACK subroutines CGECO
!    and CGESL.  That is, if A is an NxN complex matrix
!    and if X  and B are complex  N-vectors, then CGEFS
!    solves the equation
!
!                          A*X=B.
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
!    ing factorization (ITASK=1).  IND will not be changed by CGEFS
!    in this case.
!
!  Argument Description ***
!
!    A      COMPLEX(LDA,N)
!             on entry, the doubly subscripted array with dimension
!               (LDA,N) which contains the coefficient matrix.
!             on return, an upper triangular matrix U and the
!               multipliers necessary to construct a matrix L
!               so that A=L*U.
!    LDA    INTEGER
!             the leading dimension of the array A.  LDA must be great-
!             er than or equal to N.  (Terminal error message IND=-1)
!    N      INTEGER
!             the order of the matrix A.  The first N elements of
!             the array A are the elements of the first column of
!             the matrix A.  N must be greater than or equal to 1.
!             (Terminal error message IND=-2)
!    V      COMPLEX(N)
!             on entry, the singly subscripted array(vector) of di-
!               mension N which contains the right hand side B of a
!               system of simultaneous linear equations A*X=B.
!             on return, V contains the solution vector, X .
!    ITASK  INTEGER
!             if ITASK=1, the matrix A is factored and then the
!               linear equation is solved.
!             if ITASK  >  1, the equation is solved using the existing
!               factored matrix A and IWORK.
!             if ITASK  <  1, then terminal error message IND=-3 is
!               printed.
!    IND    INTEGER
!             GT.0  IND is a rough estimate of the number of digits
!                     of accuracy in the solution, X.
!             LT.0  see error message corresponding to IND below.
!    WORK   COMPLEX(N)
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
!    IND=-10 warning    The solution has no apparent significance.
!                         The solution may be inaccurate or the matrix
!                         A may be poorly scaled.
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
!***ROUTINES CALLED  CGECO, CGESL, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800328  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Convert XERRWV calls to XERMSG calls, cvt GOTO's to
!           IF-THEN-ELSE.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CGEFS
!
  INTEGER LDA,N,ITASK,IND,IWORK(*)
  COMPLEX A(LDA,*),V(*),WORK(*)
  REAL R1MACH
  REAL RCOND
  CHARACTER*8 XERN1, XERN2
!***FIRST EXECUTABLE STATEMENT  CGEFS
  if (LDA < N) THEN
     IND = -1
     WRITE (XERN1, '(I8)') LDA
     WRITE (XERN2, '(I8)') N
     call XERMSG ('SLATEC', 'CGEFS', 'LDA = ' // XERN1 // &
        ' IS LESS THAN N = ' // XERN2, -1, 1)
     return
  end if
!
  if (N <= 0) THEN
     IND = -2
     WRITE (XERN1, '(I8)') N
     call XERMSG ('SLATEC', 'CGEFS', 'N = ' // XERN1 // &
        ' IS LESS THAN 1', -2, 1)
     return
  end if
!
  if (ITASK < 1) THEN
     IND = -3
     WRITE (XERN1, '(I8)') ITASK
     call XERMSG ('SLATEC', 'CGEFS', 'ITASK = ' // XERN1 // &
        ' IS LESS THAN 1', -3, 1)
     return
  end if
!
!     FACTOR MATRIX A INTO LU
!
  if (ITASK == 1) THEN
     call CGECO(A,LDA,N,IWORK,RCOND,WORK)
!
!        CHECK FOR COMPUTATIONALLY SINGULAR MATRIX
!
     if (RCOND == 0.0) THEN
        IND = -4
        call XERMSG ('SLATEC', 'CGEFS', &
           'SINGULAR MATRIX A - NO SOLUTION', -4, 1)
        return
     ENDIF
!
!        COMPUTE IND (ESTIMATE OF NO. OF SIGNIFICANT DIGITS)
!
     IND = -LOG10(R1MACH(4)/RCOND)
!
!        CHECK FOR IND GREATER THAN ZERO
!
     if (IND <= 0) THEN
        IND = -10
        call XERMSG ('SLATEC', 'CGEFS', &
           'SOLUTION MAY HAVE NO SIGNIFICANCE', -10, 0)
     ENDIF
  end if
!
!     SOLVE AFTER FACTORING
!
  call CGESL(A,LDA,N,IWORK,V,0)
  return
end
