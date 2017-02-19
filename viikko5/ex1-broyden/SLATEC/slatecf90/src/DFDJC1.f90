subroutine DFDJC1 (FCN, N, X, FVEC, FJAC, LDFJAC, IFLAG, ML, MU, &
     EPSFCN, WA1, WA2)
!
!! DFDJC1 computes a forward difference approximation to an N by N Jacobian.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (FDJAC1-S, DFDJC1-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine computes a forward-difference approximation
!     to the N by N Jacobian matrix associated with a specified
!     problem of N functions in N variables. If the Jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     The subroutine statement is
!
!       SUBROUTINE DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
!                         WA1,WA2)
!
!     where
!
!       FCN is the name of the user-supplied subroutine which
!         calculates the functions. FCN must be declared
!         in an EXTERNAL statement in the user calling
!         program, and should be written as follows.
!
!         SUBROUTINE FCN(N,X,FVEC,IFLAG)
!         INTEGER N,IFLAG
!         DOUBLE PRECISION X(N),FVEC(N)
!         ----------
!         Calculate the functions at X and
!         return this vector in FVEC.
!         ----------
!         return
!
!         The value of IFLAG should not be changed by FCN unless
!         the user wants to terminate execution of DFDJC1.
!         In this case set IFLAG to a negative integer.
!
!       N is a positive integer input variable set to the number
!         of functions and variables.
!
!       X is an input array of length N.
!
!       FVEC is an input array of length N which must contain the
!         functions evaluated at X.
!
!       FJAC is an output N by N array which contains the
!         approximation to the Jacobian matrix evaluated at X.
!
!       LDFJAC is a positive integer input variable not less than N
!         which specifies the leading dimension of the array FJAC.
!
!       IFLAG is an integer variable which can be used to terminate
!         the execution of DFDJC1. See description of FCN.
!
!       ML is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         Jacobian matrix. If the Jacobian is not banded, set
!         ML to at least N - 1.
!
!       EPSFCN is an input variable used in determining a suitable
!         step length for the forward-difference approximation. This
!         approximation assumes that the relative errors in the
!         functions are of the order of EPSFCN. If EPSFCN is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       MU is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         Jacobian matrix. If the Jacobian is not banded, set
!         MU to at least N - 1.
!
!       WA1 and WA2 are work arrays of length N. If ML + MU + 1 is at
!         least N, then the Jacobian is considered dense, and WA2 is
!         not referenced.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DFDJC1
  DOUBLE PRECISION D1MACH
  INTEGER I, IFLAG, J, K, LDFJAC, ML, MSUM, MU, N
  DOUBLE PRECISION EPS, EPSFCN, EPSMCH, FJAC(LDFJAC,*), &
       FVEC(*), H, TEMP, WA1(*), WA2(*), X(*), ZERO
  SAVE ZERO
  DATA ZERO /0.0D0/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
!***FIRST EXECUTABLE STATEMENT  DFDJC1
  EPSMCH = D1MACH(4)
!
  EPS = SQRT(MAX(EPSFCN,EPSMCH))
  MSUM = ML + MU + 1
  if (MSUM  <  N) go to 40
!
!        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
!
     DO 20 J = 1, N
        TEMP = X(J)
        H = EPS*ABS(TEMP)
        if (H  ==  ZERO) H = EPS
        X(J) = TEMP + H
        call FCN(N,X,WA1,IFLAG)
        if (IFLAG  <  0) go to 30
        X(J) = TEMP
        DO 10 I = 1, N
           FJAC(I,J) = (WA1(I) - FVEC(I))/H
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
     go to 110
   40 CONTINUE
!
!        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
!
     DO 90 K = 1, MSUM
        DO 60 J = K, N, MSUM
           WA2(J) = X(J)
           H = EPS*ABS(WA2(J))
           if (H  ==  ZERO) H = EPS
           X(J) = WA2(J) + H
   60          CONTINUE
        call FCN(N,X,WA1,IFLAG)
        if (IFLAG  <  0) go to 100
        DO 80 J = K, N, MSUM
           X(J) = WA2(J)
           H = EPS*ABS(WA2(J))
           if (H  ==  ZERO) H = EPS
           DO 70 I = 1, N
              FJAC(I,J) = ZERO
              if (I  >=  J - MU .AND. I  <=  J + ML) &
                 FJAC(I,J) = (WA1(I) - FVEC(I))/H
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
  return
!
!     LAST CARD OF SUBROUTINE DFDJC1.
!
end
