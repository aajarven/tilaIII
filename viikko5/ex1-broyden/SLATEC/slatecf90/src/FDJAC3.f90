subroutine FDJAC3 (FCN, M, N, X, FVEC, FJAC, LDFJAC, IFLAG, &
     EPSFCN, WA)
!
!! FDJAC3 is subsidiary to SNLS1 and SNLS1E.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (FDJAC3-S, DFDJC3-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine computes a forward-difference approximation
!     to the M by N Jacobian matrix associated with a specified
!     problem of M functions in N variables.
!
!     The subroutine statement is
!
!       SUBROUTINE FDJAC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA)
!
!     where
!
!       FCN is the name of the user-supplied subroutine which
!         calculates the functions. FCN must be declared
!         in an external statement in the user calling
!         program, and should be written as follows.
!
!         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!         INTEGER LDFJAC,M,N,IFLAG
!         REAL X(N),FVEC(M),FJAC(LDFJAC,N)
!         ----------
!         When IFLAG == 1 calculate the functions at X and
!         return this vector in FVEC.
!         ----------
!         return
!         END
!
!         The value of IFLAG should not be changed by FCN unless
!         the user wants to terminate execution of FDJAC3.
!         In this case set IFLAG to a negative integer.
!
!       M is a positive integer input variable set to the number
!         of functions.
!
!       N is a positive integer input variable set to the number
!         of variables. N must not exceed M.
!
!       X is an input array of length N.
!
!       FVEC is an input array of length M which must contain the
!         functions evaluated at X.
!
!       FJAC is an output M by N array which contains the
!         approximation to the Jacobian matrix evaluated at X.
!
!       LDFJAC is a positive integer input variable not less than M
!         which specifies the leading dimension of the array FJAC.
!
!       IFLAG is an integer variable which can be used to terminate
!         THE EXECUTION OF FDJAC3. See description of FCN.
!
!       EPSFCN is an input variable used in determining a suitable
!         step length for the forward-difference approximation. This
!         approximation assumes that the relative errors in the
!         functions are of the order of EPSFCN. If EPSFCN is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       WA is a work array of length M.
!
!***SEE ALSO  SNLS1, SNLS1E
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  FDJAC3
  INTEGER M,N,LDFJAC,IFLAG
  REAL EPSFCN
  REAL X(*),FVEC(*),FJAC(LDFJAC,*),WA(*)
  INTEGER I,J
  REAL EPS,EPSMCH,H,TEMP,ZERO
  REAL R1MACH
  SAVE ZERO
  DATA ZERO /0.0E0/
!***FIRST EXECUTABLE STATEMENT  FDJAC3
  EPSMCH = R1MACH(4)
!
  EPS = SQRT(MAX(EPSFCN,EPSMCH))
!      SET IFLAG=1 TO INDICATE THAT FUNCTION VALUES
!           ARE TO BE RETURNED BY FCN.
  IFLAG = 1
  DO 20 J = 1, N
     TEMP = X(J)
     H = EPS*ABS(TEMP)
     if (H  ==  ZERO) H = EPS
     X(J) = TEMP + H
     call FCN(IFLAG,M,N,X,WA,FJAC,LDFJAC)
     if (IFLAG  <  0) go to 30
     X(J) = TEMP
     DO 10 I = 1, M
        FJAC(I,J) = (WA(I) - FVEC(I))/H
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
  return
!
!     LAST CARD OF SUBROUTINE FDJAC3.
!
end
