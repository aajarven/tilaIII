subroutine DROTM (N, DX, INCX, DY, INCY, DPARAM)
!
!! DROTM applies a modified Givens rotation.
!
!***PURPOSE  Apply a modified Givens transformation.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A8
!***TYPE      DOUBLE PRECISION (SROTM-S, DROTM-D)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR
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
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!   DPARAM  5-element D.P. vector.  DPARAM(1) is DFLAG described below.
!           Locations 2-5 of SPARAM contain elements of the
!           transformation matrix H described below.
!
!     --Output--
!       DX  rotated vector (unchanged if N  <=  0)
!       DY  rotated vector (unchanged if N  <=  0)
!
!     Apply the modified Givens transformation, H, to the 2 by N matrix
!     (DX**T)
!     (DY**T) , where **T indicates transpose.  The elements of DX are
!     in DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX  >=  0, else
!     LX = 1+(1-N)*INCX, and similarly for DY using LY and INCY.
!
!     With DPARAM(1)=DFLAG, H has one of the following forms:
!
!     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
!
!       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
!     H=(          )    (          )    (          )    (          )
!       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
!
!     See DROTMG for a description of data storage in DPARAM.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DROTM
  DOUBLE PRECISION DFLAG, DH12, DH22, DX, TWO, Z, DH11, DH21, &
                   DPARAM, DY, W, ZERO
  DIMENSION DX(*), DY(*), DPARAM(5)
  SAVE ZERO, TWO
  DATA ZERO, TWO /0.0D0, 2.0D0/
!***FIRST EXECUTABLE STATEMENT  DROTM
  DFLAG=DPARAM(1)
  if (N <= 0 .OR. (DFLAG+TWO == ZERO)) go to 140
      if (.NOT.(INCX == INCY.AND. INCX  > 0)) go to 70
!
           NSTEPS=N*INCX
           if (DFLAG) 50,10,30
   10          CONTINUE
           DH12=DPARAM(4)
           DH21=DPARAM(3)
                DO 20 I = 1,NSTEPS,INCX
                W=DX(I)
                Z=DY(I)
                DX(I)=W+Z*DH12
                DY(I)=W*DH21+Z
   20               CONTINUE
           go to 140
   30          CONTINUE
           DH11=DPARAM(2)
           DH22=DPARAM(5)
                DO 40 I = 1,NSTEPS,INCX
                W=DX(I)
                Z=DY(I)
                DX(I)=W*DH11+Z
                DY(I)=-W+DH22*Z
   40               CONTINUE
           go to 140
   50          CONTINUE
           DH11=DPARAM(2)
           DH12=DPARAM(4)
           DH21=DPARAM(3)
           DH22=DPARAM(5)
                DO 60 I = 1,NSTEPS,INCX
                W=DX(I)
                Z=DY(I)
                DX(I)=W*DH11+Z*DH12
                DY(I)=W*DH21+Z*DH22
   60               CONTINUE
           go to 140
   70     CONTINUE
      KX=1
      KY=1
      if (INCX  <  0) KX = 1+(1-N)*INCX
      if (INCY  <  0) KY = 1+(1-N)*INCY
!
      if (DFLAG) 120,80,100
   80     CONTINUE
      DH12=DPARAM(4)
      DH21=DPARAM(3)
           DO 90 I = 1,N
           W=DX(KX)
           Z=DY(KY)
           DX(KX)=W+Z*DH12
           DY(KY)=W*DH21+Z
           KX=KX+INCX
           KY=KY+INCY
   90          CONTINUE
      go to 140
  100     CONTINUE
      DH11=DPARAM(2)
      DH22=DPARAM(5)
           DO 110 I = 1,N
           W=DX(KX)
           Z=DY(KY)
           DX(KX)=W*DH11+Z
           DY(KY)=-W+DH22*Z
           KX=KX+INCX
           KY=KY+INCY
  110          CONTINUE
      go to 140
  120     CONTINUE
      DH11=DPARAM(2)
      DH12=DPARAM(4)
      DH21=DPARAM(3)
      DH22=DPARAM(5)
           DO 130 I = 1,N
           W=DX(KX)
           Z=DY(KY)
           DX(KX)=W*DH11+Z*DH12
           DY(KY)=W*DH21+Z*DH22
           KX=KX+INCX
           KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
      return
end
