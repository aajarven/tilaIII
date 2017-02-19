subroutine DPCHBS (N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF, &
     NDIM, KORD, IERR)
!
!! DPCHBS is a piecewise Cubic Hermite to B-Spline converter.
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3
!***TYPE      DOUBLE PRECISION (PCHBS-S, DPCHBS-D)
!***KEYWORDS  B-SPLINES, CONVERSION, CUBIC HERMITE INTERPOLATION,
!             PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Computing and Mathematics Research Division
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
! *Usage:
!
!        INTEGER  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
!        PARAMETER  (INCFD = ...)
!        DOUBLE PRECISION  X(nmax), F(INCFD,nmax), D(INCFD,nmax),
!       *      T(2*nmax+4), BCOEF(2*nmax)
!
!        call DPCHBS (N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF,
!       *             NDIM, KORD, IERR)
!
! *Arguments:
!
!     N:IN  is the number of data points, N.ge.2 .  (not checked)
!
!     X:IN  is the real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1)  <  X(I),  I = 2(1)N.   (not checked)
!           nmax, the dimension of X, must be .ge.N.
!
!     F:IN  is the real array of dependent variable values.
!           F(1+(I-1)*INCFD) is the value corresponding to X(I).
!           nmax, the second dimension of F, must be .ge.N.
!
!     D:IN  is the real array of derivative values at the data points.
!           D(1+(I-1)*INCFD) is the value corresponding to X(I).
!           nmax, the second dimension of D, must be .ge.N.
!
!     INCFD:IN  is the increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           It may have the value 1 for one-dimensional applications,
!           in which case F and D may be singly-subscripted arrays.
!
!     KNOTYP:IN  is a flag to control the knot sequence.
!           The knot sequence T is normally computed from X by putting
!           a double knot at each X and setting the end knot pairs
!           according to the value of KNOTYP:
!              KNOTYP = 0:  Quadruple knots at X(1) and X(N).  (default)
!              KNOTYP = 1:  Replicate lengths of extreme subintervals:
!                           T( 1 ) = T( 2 ) = X(1) - (X(2)-X(1))  ;
!                           T(M+4) = T(M+3) = X(N) + (X(N)-X(N-1)).
!              KNOTYP = 2:  Periodic placement of boundary knots:
!                           T( 1 ) = T( 2 ) = X(1) - (X(N)-X(N-1));
!                           T(M+4) = T(M+3) = X(N) + (X(2)-X(1))  .
!              Here M=NDIM=2*N.
!           If the input value of KNOTYP is negative, however, it is
!           assumed that NKNOTS and T were set in a previous call.
!           This option is provided for improved efficiency when used
!           in a parametric setting.
!
!     NKNOTS:INOUT  is the number of knots.
!           If KNOTYP >= 0, then NKNOTS will be set to NDIM+4.
!           If KNOTYP < 0, then NKNOTS is an input variable, and an
!              error return will be taken if it is not equal to NDIM+4.
!
!     T:INOUT  is the array of 2*N+4 knots for the B-representation.
!           If KNOTYP >= 0, T will be returned by DPCHBS with the
!              interior double knots equal to the X-values and the
!              boundary knots set as indicated above.
!           If KNOTYP < 0, it is assumed that T was set by a
!              previous call to DPCHBS.  (This routine does **not**
!              verify that T forms a legitimate knot sequence.)
!
!     BCOEF:OUT  is the array of 2*N B-spline coefficients.
!
!     NDIM:OUT  is the dimension of the B-spline space.  (Set to 2*N.)
!
!     KORD:OUT  is the order of the B-spline.  (Set to 4.)
!
!     IERR:OUT  is an error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -4  if KNOTYP > 2 .
!              IERR = -5  if KNOTYP < 0 and NKNOTS /= (2*N+4).
!
! *Description:
!     DPCHBS computes the B-spline representation of the PCH function
!     determined by N,X,F,D.  To be compatible with the rest of PCHIP,
!     DPCHBS includes INCFD, the increment between successive values of
!     the F and D arrays.
!
!     The output is the B-representation for the function:  NKNOTS, T,
!     BCOEF, NDIM, KORD.
!
! *Caution:
!     Since it is assumed that the input PCH function has been
!     computed by one of the other routines in the package PCHIP,
!     input arguments N, X, INCFD are **not** checked for validity.
!
! *Restrictions/assumptions:
!     1. N >= 2 .  (not checked)
!     2. X(i) < X(i+1), i=1,...,N .  (not checked)
!     3. INCFD > 0 .  (not checked)
!     4. KNOTYP <= 2 .  (error return if not)
!    *5. NKNOTS = NDIM+4 = 2*N+4 .  (error return if not)
!    *6. T(2*k+1) = T(2*k) = X(k), k=1,...,N .  (not checked)
!
!       * Indicates this applies only if KNOTYP < 0 .
!
! *Portability:
!     Argument INCFD is used only to cause the compiler to generate
!     efficient code for the subscript expressions (1+(I-1)*INCFD) .
!     The normal usage, in which DPCHBS is called with one-dimensional
!     arrays F and D, is probably non-Fortran 77, in the strict sense,
!     but it works on all systems on which DPCHBS has been tested.
!
! *See Also:
!     PCHIC, PCHIM, or PCHSP can be used to determine an interpolating
!        PCH function from a set of data.
!     The B-spline routine DBVALU can be used to evaluate the
!        B-representation that is output by DPCHBS.
!        (See BSPDOC for more information.)
!
!***REFERENCES  F. N. Fritsch, "Representations for parametric cubic
!                 splines," Computer Aided Geometric Design 6 (1989),
!                 pp.79-82.
!***ROUTINES CALLED  DPCHKT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   870701  DATE WRITTEN
!   900405  Converted Fortran to upper case.
!   900405  Removed requirement that X be dimensioned N+1.
!   900406  Modified to make PCHKT a subsidiary routine to simplify
!           usage.  In the process, added argument INCFD to be com-
!           patible with the rest of PCHIP.
!   900410  Converted prologue to SLATEC 4.0 format.
!   900410  Added calls to XERMSG and changed constant 3. to 3 to
!           reduce single/double differences.
!   900411  Added reference.
!   900430  Produced double precision version.
!   900501  Corrected declarations.
!   930317  Minor cosmetic changes.  (FNF)
!   930514  Corrected problems with dimensioning of arguments and
!           clarified DESCRIPTION.  (FNF)
!   930604  Removed  NKNOTS from DPCHKT call list.  (FNF)
!***END PROLOGUE  DPCHBS
!
!*Internal Notes:
!
!**End
!
!  Declare arguments.
!
  INTEGER  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
  DOUBLE PRECISION  X(*), F(INCFD,*), D(INCFD,*), T(*), BCOEF(*)
!
!  Declare local variables.
!
  INTEGER  K, KK
  DOUBLE PRECISION  DOV3, HNEW, HOLD
  CHARACTER*8  LIBNAM, SUBNAM
!***FIRST EXECUTABLE STATEMENT  DPCHBS
!
!  Initialize.
!
  NDIM = 2*N
  KORD = 4
  IERR = 0
  LIBNAM = 'SLATEC'
  SUBNAM = 'DPCHBS'
!
!  Check argument validity.  Set up knot sequence if OK.
!
  if ( KNOTYP > 2 )  THEN
     IERR = -1
     call XERMSG (LIBNAM, SUBNAM, 'KNOTYP GREATER THAN 2', IERR, 1)
     return
  end if
  if ( KNOTYP < 0 )  THEN
     if ( NKNOTS /= NDIM+4 )  THEN
        IERR = -2
        call XERMSG (LIBNAM, SUBNAM, &
                      'KNOTYP < 0 AND NKNOTS /= (2*N+4)', IERR, 1)
        return
     ENDIF
  ELSE
!          Set up knot sequence.
     NKNOTS = NDIM + 4
     call DPCHKT (N, X, KNOTYP, T)
  end if
!
!  Compute B-spline coefficients.
!
  HNEW = T(3) - T(1)
  DO 40  K = 1, N
     KK = 2*K
     HOLD = HNEW
!          The following requires mixed mode arithmetic.
     DOV3 = D(1,K)/3
     BCOEF(KK-1) = F(1,K) - HOLD*DOV3
!          The following assumes T(2*K+1) = X(K).
     HNEW = T(KK+3) - T(KK+1)
     BCOEF(KK) = F(1,K) + HNEW*DOV3
   40 CONTINUE
!
!  Terminate.
!
  return
!------------- LAST LINE OF DPCHBS FOLLOWS -----------------------------
end
