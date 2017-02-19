  DOUBLE PRECISION FUNCTION DQDOTI (N, DB, QC, DX, INCX, DY, INCY)
!
!! DQDOTI computes the inner product of two vectors with extended precision.
!
!***LIBRARY   SLATEC
!***CATEGORY  D1A4
!***TYPE      DOUBLE PRECISION (DQDOTI-D)
!***KEYWORDS  DOT PRODUCT, INNER PRODUCT
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       DB  double precision scalar to be added to inner product
!       QC  extended precision scalar to be added
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!   DQDOTI  double precision result
!       QC  extended precision result
!
!     D.P. dot product with extended precision accumulation (and result)
!     QC and DQDOTI are set = DB + sum for I = 0 to N-1 of
!       DX(LX+I*INCX) * DY(LY+I*INCY),  where QC is an extended
!       precision result which can be used as input to DQDOTA,
!       and LX = 1 if INCX  >=  0, else LX = (-INCX)*N, and LY is
!       defined in a similar way using INCY.  The MP package by
!       Richard P. Brent is used for the extended precision arithmetic.
!
!     Fred T. Krogh,  JPL,  1977,  June 1
!
!     The common block for the MP package is named MPCOM.  If local
!     variable I1 is zero, DQDOTI calls MPBLAS to initialize the MP
!     package and reset I1 to 1.
!
!    The argument QC(*), and the local variables QX and QY are INTEGER
!    arrays of size 30.  See the comments in the routine MPBLAS for the
!    reason for this choice.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  MPADD, MPBLAS, MPCDM, MPCMD, MPMUL
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   930124  Increased Array sizes for SUN -r8.  (RWC)
!***END PROLOGUE  DQDOTI
  DOUBLE PRECISION DX(*), DY(*), DB
  INTEGER  QC(30), QX(30), QY(30)
  COMMON /MPCOM/  MPB, MPT, MPM, MPLUN, MPMXR, MPR(30)
  SAVE I1
  DATA  I1 / 0 /
!***FIRST EXECUTABLE STATEMENT  DQDOTI
  if (I1  ==  0) call MPBLAS(I1)
  QC(1) = 0
  if (DB  ==  0.D0) go to 60
  call MPCDM(DB, QX)
  call MPADD(QC, QX, QC)
   60 if (N  ==  0) go to 80
  IX = 1
  IY = 1
  if (INCX  <  0) IX = (-N + 1) * INCX + 1
  if (INCY  <  0) IY = (-N + 1) * INCY + 1
  DO  70  I = 1,N
     call MPCDM(DX(IX), QX)
     call MPCDM(DY(IY), QY)
     call MPMUL(QX, QY, QX)
     call MPADD(QC, QX, QC)
     IX = IX + INCX
     IY = IY + INCY
   70 CONTINUE
   80 call MPCMD(QC, DQDOTI)
  return
end
