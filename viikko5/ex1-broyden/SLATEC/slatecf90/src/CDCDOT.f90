FUNCTION CDCDOT (N, CB, CX, INCX, CY, INCY)
!
!! CDCDOT computes the inner product of two vectors with extended ...
!            precision accumulation.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      COMPLEX (SDSDOT-S, CDCDOT-C)
!***KEYWORDS  BLAS, DOT PRODUCT, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
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
!       CB  complex scalar to be added to inner product
!       CX  complex vector with N elements
!     INCX  storage spacing between elements of CX
!       CY  complex vector with N elements
!     INCY  storage spacing between elements of CY
!
!     --Output--
!   CDCDOT  complex dot product (CB if N  <=  0)
!
!     Returns complex result with dot product accumulated in D.P.
!     CDCDOT = CB + sum for I = 0 to N-1 of CX(LX+I*INCY)*CY(LY+I*INCY)
!     where LX = 1 if INCX  >=  0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CDCDOT
  COMPLEX CDCDOT
  INTEGER N, INCX, INCY, I, KX, KY
  COMPLEX CX(*), CY(*), CB
  DOUBLE PRECISION DSDOTR, DSDOTI, DT1, DT2, DT3, DT4
!***FIRST EXECUTABLE STATEMENT  CDCDOT
  DSDOTR = DBLE(REAL(CB))
  DSDOTI = DBLE(AIMAG(CB))
  if (N  <=  0) go to 10
  KX = 1
  KY = 1
  if ( INCX < 0) KX = 1+(1-N)*INCX
  if ( INCY < 0) KY = 1+(1-N)*INCY
  DO 5 I = 1,N
    DT1 = DBLE(REAL(CX(KX)))
    DT2 = DBLE(REAL(CY(KY)))
    DT3 = DBLE(AIMAG(CX(KX)))
    DT4 = DBLE(AIMAG(CY(KY)))
    DSDOTR = DSDOTR+(DT1*DT2)-(DT3*DT4)
    DSDOTI = DSDOTI+(DT1*DT4)+(DT3*DT2)
    KX = KX+INCX
    KY = KY+INCY
    5 CONTINUE
   10 CDCDOT = CMPLX(REAL(DSDOTR),REAL(DSDOTI))
  return
end
