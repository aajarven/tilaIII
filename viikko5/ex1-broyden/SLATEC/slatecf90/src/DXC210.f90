subroutine DXC210 (K, Z, J, IERROR)
!
!! DXC210 provides double-precision floating-point arithmetic ...
!            with an extended exponent range.
!
!***LIBRARY   SLATEC
!***CATEGORY  A3D
!***TYPE      DOUBLE PRECISION (XC210-S, DXC210-D)
!***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
!***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
!           Smith, John M., (NBS and George Mason University)
!***DESCRIPTION
!     INTEGER K, J
!     DOUBLE PRECISION Z
!
!                  GIVEN K THIS SUBROUTINE COMPUTES J AND Z
!                  SUCH THAT  RADIX**K = Z*10**J, WHERE Z IS IN
!                  THE RANGE 1/10  <=  Z  <  1.
!                  THE VALUE OF Z WILL BE ACCURATE TO FULL
!                  DOUBLE-PRECISION PROVIDED THE NUMBER
!                  OF DECIMAL PLACES IN THE LARGEST
!                  INTEGER PLUS THE NUMBER OF DECIMAL
!                  PLACES CARRIED IN DOUBLE-PRECISION DOES NOT
!                  EXCEED 60. DXC210 IS CALLED BY SUBROUTINE
!                  DXCON WHEN NECESSARY. THE USER SHOULD
!                  NEVER NEED TO call DXC210 DIRECTLY.
!
!***SEE ALSO  DXSET
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***COMMON BLOCKS    DXBLK3
!***REVISION HISTORY  (YYMMDD)
!   820712  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXC210
  DOUBLE PRECISION Z
  INTEGER K, J
  INTEGER NLG102, MLG102, LG102
  COMMON /DXBLK3/ NLG102, MLG102, LG102(21)
  SAVE /DXBLK3/
!
!   THE CONDITIONS IMPOSED ON NLG102, MLG102, AND LG102 BY
! THIS SUBROUTINE ARE
!
!     (1) NLG102  >=  2
!
!     (2) MLG102  >=  1
!
!     (3) 2*MLG102*(MLG102 - 1)  <=  2**NBITS - 1
!
! THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
! IN SUBROUTINE DXSET.
!
!***FIRST EXECUTABLE STATEMENT  DXC210
  IERROR=0
  if (K == 0) go to 70
  M = MLG102
  KA = ABS(K)
  KA1 = KA/M
  KA2 = MOD(KA,M)
  if (KA1 >= M) go to 60
  NM1 = NLG102 - 1
  NP1 = NLG102 + 1
  IT = KA2*LG102(NP1)
  IC = IT/M
  ID = MOD(IT,M)
  Z = ID
  if (KA1 > 0) go to 20
  DO 10 II=1,NM1
    I = NP1 - II
    IT = KA2*LG102(I) + IC
    IC = IT/M
    ID = MOD(IT,M)
    Z = Z/M + ID
   10 CONTINUE
  JA = KA*LG102(1) + IC
  go to 40
   20 CONTINUE
  DO 30 II=1,NM1
    I = NP1 - II
    IT = KA2*LG102(I) + KA1*LG102(I+1) + IC
    IC = IT/M
    ID = MOD(IT,M)
    Z = Z/M + ID
   30 CONTINUE
  JA = KA*LG102(1) + KA1*LG102(2) + IC
   40 CONTINUE
  Z = Z/M
  if (K > 0) go to 50
  J = -JA
  Z = 10.0D0**(-Z)
  go to 80
   50 CONTINUE
  J = JA + 1
  Z = 10.0D0**(Z-1.0D0)
  go to 80
   60 CONTINUE
!   THIS ERROR OCCURS if K EXCEEDS  MLG102**2 - 1  IN MAGNITUDE.
!
  call XERMSG ('SLATEC', 'DXC210', 'K too large', 208, 1)
  IERROR=208
  return
   70 CONTINUE
  J = 0
  Z = 1.0D0
   80 RETURN
end
