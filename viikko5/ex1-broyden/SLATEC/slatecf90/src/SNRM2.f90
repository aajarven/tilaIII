FUNCTION SNRM2 (N, SX, INCX)
!
!! SNRM2 computes the Euclidean length (L2 norm) of a vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A3B
!***TYPE      SINGLE PRECISION (SNRM2-S, DNRM2-D, SCNRM2-C)
!***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
!             LINEAR ALGEBRA, UNITARY, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!
!     --Output--
!    SNRM2  single precision result (zero if N  <=  0)
!
!     Euclidean norm of the N-vector stored in SX with storage
!     increment INCX .
!     If N  <=  0, return with result = 0.
!     If N  >=  1, then INCX must be  >=  1
!
!     Four Phase Method using two built-in constants that are
!     hopefully applicable to all machines.
!         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
!         CUTHI = minimum of  SQRT(V)      over all known machines.
!     where
!         EPS = smallest no. such that EPS + 1.  >  1.
!         U   = smallest positive no.   (underflow limit)
!         V   = largest  no.            (overflow  limit)
!
!     Brief Outline of Algorithm.
!
!     Phase 1 scans zero components.
!     Move to phase 2 when a component is nonzero and  <=  CUTLO
!     Move to phase 3 when a component is  >  CUTLO
!     Move to phase 4 when a component is  >=  CUTHI/M
!     where M = N for X() real and M = 2*N for complex.
!
!     Values for CUTLO and CUTHI.
!     From the environmental parameters listed in the IMSL converter
!     document the limiting values are as follows:
!     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
!                   Univac and DEC at 2**(-103)
!                   Thus CUTLO = 2**(-51) = 4.44089E-16
!     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
!                   Thus CUTHI = 2**(63.5) = 1.30438E19
!     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
!                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
!     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
!     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
!     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SNRM2
  real SNRM2
  INTEGER NEXT
  REAL SX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO, ONE
  SAVE CUTLO, CUTHI, ZERO, ONE
  DATA ZERO, ONE /0.0E0, 1.0E0/
!
  DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
!***FIRST EXECUTABLE STATEMENT  SNRM2
  if (N  >  0) go to 10
     SNRM2  = ZERO
     go to 300
!
   10 ASSIGN 30 TO NEXT
  SUM = ZERO
  NN = N * INCX
!
!                                                 BEGIN MAIN LOOP
!
  I = 1
   20    go to NEXT,(30, 50, 70, 110)
   30 if (ABS(SX(I))  >  CUTLO) go to 85
  ASSIGN 50 TO NEXT
  XMAX = ZERO
!
!                        PHASE 1.  SUM IS ZERO
!
   50 if (SX(I)  ==  ZERO) go to 200
  if (ABS(SX(I))  >  CUTLO) go to 85
!
!                                PREPARE FOR PHASE 2.
!
  ASSIGN 70 TO NEXT
  go to 105
!
!                                PREPARE FOR PHASE 4.
!
  100 I = J
  ASSIGN 110 TO NEXT
  SUM = (SUM / SX(I)) / SX(I)
  105 XMAX = ABS(SX(I))
  go to 115
!
!                   PHASE 2.  SUM IS SMALL.
!                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
!
   70 if (ABS(SX(I))  >  CUTLO) go to 75
!
!                     COMMON CODE FOR PHASES 2 AND 4.
!                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
!
  110 if (ABS(SX(I))  <=  XMAX) go to 115
     SUM = ONE + SUM * (XMAX / SX(I))**2
     XMAX = ABS(SX(I))
     go to 200
!
  115 SUM = SUM + (SX(I)/XMAX)**2
  go to 200
!
!                  PREPARE FOR PHASE 3.
!
   75 SUM = (SUM * XMAX) * XMAX
!
!     FOR REAL OR D.P. SET HITEST = CUTHI/N
!     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
!
   85 HITEST = CUTHI / N
!
!                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
!
  DO 95 J = I,NN,INCX
  if (ABS(SX(J))  >=  HITEST) go to 100
   95    SUM = SUM + SX(J)**2
  SNRM2 = SQRT( SUM )
  go to 300
!
  200 CONTINUE
  I = I + INCX
  if (I  <=  NN) go to 20
!
!              END OF MAIN LOOP.
!
!              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
!
  SNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
  return
end
