subroutine DROTMG (DD1, DD2, DX1, DY1, DPARAM)
!
!! DROTMG constructs a modified Givens rotation.
!
!***PURPOSE  Construct a modified Givens transformation.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B10
!***TYPE      DOUBLE PRECISION (SROTMG-S, DROTMG-D)
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
!      DD1  double precision scalar
!      DD2  double precision scalar
!      DX1  double precision scalar
!      DX2  double precision scalar
!   DPARAM  D.P. 5-vector. DPARAM(1)=DFLAG defined below.
!           Locations 2-5 contain the rotation matrix.
!
!     --Output--
!      DD1  changed to represent the effect of the transformation
!      DD2  changed to represent the effect of the transformation
!      DX1  changed to represent the effect of the transformation
!      DX2  unchanged
!
!     Construct the modified Givens transformation matrix H which zeros
!     the second component of the 2-vector  (SQRT(DD1)*DX1,SQRT(DD2)*
!     DY2)**T.
!     With DPARAM(1)=DFLAG, H has one of the following forms:
!
!     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
!
!       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
!     H=(          )    (          )    (          )    (          )
!       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
!
!     Locations 2-5 of DPARAM contain DH11, DH21, DH12, and DH22,
!     respectively.  (Values of 1.D0, -1.D0, or 0.D0 implied by the
!     value of DPARAM(1) are not stored in DPARAM.)
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   780301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920316  Prologue corrected.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DROTMG
  DOUBLE PRECISION GAM, ONE, RGAMSQ, DD1, DD2, DH11, DH12, DH21, &
                   DH22, DPARAM, DP1, DP2, DQ1, DQ2, DU, DY1, ZERO, &
                   GAMSQ, DFLAG, DTEMP, DX1, TWO
  DIMENSION DPARAM(5)
  SAVE ZERO, ONE, TWO, GAM, GAMSQ, RGAMSQ
  DATA ZERO, ONE, TWO /0.0D0, 1.0D0, 2.0D0/
  DATA GAM, GAMSQ, RGAMSQ /4096.0D0, 16777216.D0, 5.9604645D-8/
!***FIRST EXECUTABLE STATEMENT  DROTMG
  if (.NOT. DD1  <  ZERO) go to 10
!       GO ZERO-H-D-AND-DX1..
      go to 60
   10 CONTINUE
!     CASE-DD1-NONNEGATIVE
  DP2=DD2*DY1
  if (.NOT. DP2  ==  ZERO) go to 20
      DFLAG=-TWO
      go to 260
!     REGULAR-CASE..
   20 CONTINUE
  DP1=DD1*DX1
  DQ2=DP2*DY1
  DQ1=DP1*DX1
!
  if (.NOT. ABS(DQ1)  >  ABS(DQ2)) go to 40
      DH21=-DY1/DX1
      DH12=DP2/DP1
!
      DU=ONE-DH12*DH21
!
      if (.NOT. DU  <=  ZERO) go to 30
!         GO ZERO-H-D-AND-DX1..
           go to 60
   30     CONTINUE
           DFLAG=ZERO
           DD1=DD1/DU
           DD2=DD2/DU
           DX1=DX1*DU
!         GO SCALE-CHECK..
           go to 100
   40 CONTINUE
      if (.NOT. DQ2  <  ZERO) go to 50
!         GO ZERO-H-D-AND-DX1..
           go to 60
   50     CONTINUE
           DFLAG=ONE
           DH11=DP1/DP2
           DH22=DX1/DY1
           DU=ONE+DH11*DH22
           DTEMP=DD2/DU
           DD2=DD1/DU
           DD1=DTEMP
           DX1=DY1*DU
!         GO SCALE-CHECK
           go to 100
!     PROCEDURE..ZERO-H-D-AND-DX1..
   60 CONTINUE
      DFLAG=-ONE
      DH11=ZERO
      DH12=ZERO
      DH21=ZERO
      DH22=ZERO
!
      DD1=ZERO
      DD2=ZERO
      DX1=ZERO
!         return..
      go to 220
!     PROCEDURE..FIX-H..
   70 CONTINUE
  if (.NOT. DFLAG  >=  ZERO) go to 90
!
      if (.NOT. DFLAG  ==  ZERO) go to 80
      DH11=ONE
      DH22=ONE
      DFLAG=-ONE
      go to 90
   80     CONTINUE
      DH21=-ONE
      DH12=ONE
      DFLAG=-ONE
   90 CONTINUE
  go to IGO,(120,150,180,210)
!     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
      if (.NOT. DD1  <=  RGAMSQ) go to 130
           if (DD1  ==  ZERO) go to 160
           ASSIGN 120 TO IGO
!              FIX-H..
           go to 70
  120          CONTINUE
           DD1=DD1*GAM**2
           DX1=DX1/GAM
           DH11=DH11/GAM
           DH12=DH12/GAM
      go to 110
  130 CONTINUE
  140     CONTINUE
      if (.NOT. DD1  >=  GAMSQ) go to 160
           ASSIGN 150 TO IGO
!              FIX-H..
           go to 70
  150          CONTINUE
           DD1=DD1/GAM**2
           DX1=DX1*GAM
           DH11=DH11*GAM
           DH12=DH12*GAM
      go to 140
  160 CONTINUE
  170     CONTINUE
      if (.NOT. ABS(DD2)  <=  RGAMSQ) go to 190
           if (DD2  ==  ZERO) go to 220
           ASSIGN 180 TO IGO
!              FIX-H..
           go to 70
  180          CONTINUE
           DD2=DD2*GAM**2
           DH21=DH21/GAM
           DH22=DH22/GAM
      go to 170
  190 CONTINUE
  200     CONTINUE
      if (.NOT. ABS(DD2)  >=  GAMSQ) go to 220
           ASSIGN 210 TO IGO
!              FIX-H..
           go to 70
  210          CONTINUE
           DD2=DD2/GAM**2
           DH21=DH21*GAM
           DH22=DH22*GAM
      go to 200
  220 CONTINUE
      if (DFLAG) 250,230,240
  230     CONTINUE
           DPARAM(3)=DH21
           DPARAM(4)=DH12
           go to 260
  240     CONTINUE
           DPARAM(2)=DH11
           DPARAM(5)=DH22
           go to 260
  250     CONTINUE
           DPARAM(2)=DH11
           DPARAM(3)=DH21
           DPARAM(4)=DH12
           DPARAM(5)=DH22
  260 CONTINUE
      DPARAM(1)=DFLAG
      return
end
