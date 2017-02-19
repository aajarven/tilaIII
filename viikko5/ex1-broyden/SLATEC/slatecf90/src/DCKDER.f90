subroutine DCKDER (M, N, X, FVEC, FJAC, LDFJAC, XP, FVECP, MODE, ERR)
!
!! DCKDER checks the gradients of M nonlinear functions in N variables...
!  evaluated at a point X, for consistency with the functions themselves.
!
!***LIBRARY   SLATEC
!***CATEGORY  F3, G4C
!***TYPE      DOUBLE PRECISION (CHKDER-S, DCKDER-D)
!***KEYWORDS  GRADIENTS, JACOBIAN, MINPACK, NONLINEAR
!***AUTHOR  Hiebert, K. L. (SNLA)
!***DESCRIPTION
!
!   This subroutine is a companion routine to DNSQ and DNSQE. It may
!   be used to check the coding of the Jacobian calculation.
!
!     SUBROUTINE DCKDER
!
!     This subroutine checks the gradients of M nonlinear functions
!     in N variables, evaluated at a point X, for consistency with
!     the functions themselves. The user must call DCKDER twice,
!     first with MODE = 1 and then with MODE = 2.
!
!     MODE = 1. On input, X must contain the point of evaluation.
!               On output, XP is set to a neighboring point.
!
!     MODE = 2. On input, FVEC must contain the functions and the
!                         rows of FJAC must contain the gradients
!                         of the respective functions each evaluated
!                         at X, and FVECP must contain the functions
!                         evaluated at XP.
!               On output, ERR contains measures of correctness of
!                          the respective gradients.
!
!     The subroutine does not perform reliably if cancellation or
!     rounding errors cause a severe loss of significance in the
!     evaluation of a function. Therefore, none of the components
!     of X should be unusually small (in particular, zero) or any
!     other value which may cause loss of significance.
!
!     The SUBROUTINE statement is
!
!       SUBROUTINE DCKDER(M,N,X,FVEC,FJAC,LDFJAC,XP,FVECP,MODE,ERR)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of functions.
!
!       N is a positive integer input variable set to the number
!         of variables.
!
!       X is an input array of length N.
!
!       FVEC is an array of length M. On input when MODE = 2,
!         FVEC must contain the functions evaluated at X.
!
!       FJAC is an M by N array. On input when MODE = 2,
!         the rows of FJAC must contain the gradients of
!         the respective functions evaluated at X.
!
!       LDFJAC is a positive integer input parameter not less than M
!         which specifies the leading dimension of the array FJAC.
!
!       XP is an array of length N. On output when MODE = 1,
!         XP is set to a neighboring point of X.
!
!       FVECP is an array of length M. On input when MODE = 2,
!         FVECP must contain the functions evaluated at XP.
!
!       MODE is an integer input variable set to 1 on the first call
!         and 2 on the second. Other values of MODE are equivalent
!         to MODE = 1.
!
!       ERR is an array of length M. On output when MODE = 2,
!         ERR contains measures of correctness of the respective
!         gradients. If there is no severe loss of significance,
!         then if ERR(I) is 1.0 the I-th gradient is correct,
!         while if ERR(I) is 0.0 the I-th gradient is incorrect.
!         For values of ERR between 0.0 and 1.0, the categorization
!         is less certain. In general, a value of ERR(I) greater
!         than 0.5 indicates that the I-th gradient is probably
!         correct, while a value of ERR(I) less than 0.5 indicates
!         that the I-th gradient is probably incorrect.
!
!***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
!                 tions. In Numerical Methods for Nonlinear Algebraic
!                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
!                 1988.
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DCKDER
  INTEGER I, J, LDFJAC, M, MODE, N
  DOUBLE PRECISION D1MACH, EPS, EPSF, EPSLOG, EPSMCH, ERR(*), &
       FACTOR, FJAC(LDFJAC,*), FVEC(*), FVECP(*), ONE, TEMP, X(*), &
       XP(*), ZERO
  SAVE FACTOR, ONE, ZERO
  DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
!***FIRST EXECUTABLE STATEMENT  DCKDER
  EPSMCH = D1MACH(4)
!
  EPS = SQRT(EPSMCH)
!
  if (MODE  ==  2) go to 20
!
!        MODE = 1.
!
     DO 10 J = 1, N
        TEMP = EPS*ABS(X(J))
        if (TEMP  ==  ZERO) TEMP = EPS
        XP(J) = X(J) + TEMP
   10       CONTINUE
     go to 70
   20 CONTINUE
!
!        MODE = 2.
!
     EPSF = FACTOR*EPSMCH
     EPSLOG = LOG10(EPS)
     DO 30 I = 1, M
        ERR(I) = ZERO
   30       CONTINUE
     DO 50 J = 1, N
        TEMP = ABS(X(J))
        if (TEMP  ==  ZERO) TEMP = ONE
        DO 40 I = 1, M
           ERR(I) = ERR(I) + TEMP*FJAC(I,J)
   40          CONTINUE
   50       CONTINUE
     DO 60 I = 1, M
        TEMP = ONE
        if (FVEC(I)  /=  ZERO .AND. FVECP(I)  /=  ZERO &
            .AND. ABS(FVECP(I)-FVEC(I))  >=  EPSF*ABS(FVEC(I))) &
           TEMP = EPS*ABS((FVECP(I)-FVEC(I))/EPS-ERR(I)) &
                  /(ABS(FVEC(I)) + ABS(FVECP(I)))
        ERR(I) = ONE
        if (TEMP  >  EPSMCH .AND. TEMP  <  EPS) &
           ERR(I) = (LOG10(TEMP) - EPSLOG)/EPSLOG
        if (TEMP  >=  EPS) ERR(I) = ZERO
   60       CONTINUE
   70 CONTINUE
!
  return
end
