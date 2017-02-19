subroutine DPLPDM (MRELAS, NVARS, LMX, LBM, NREDC, INFO, IOPT, &
     IBASIS, IMAT, IBRC, IPR, IWR, IND, IBB, ANORM, EPS, UU, GG, &
     AMAT, BASMAT, CSC, WR, SINGLR, REDBAS)
!
!! DPLPDM is subsidiary to DSPLP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SPLPDM-S, DPLPDM-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     THIS SUBPROGRAM IS FROM THE DSPLP( ) PACKAGE.  IT PERFORMS THE
!     TASK OF DEFINING THE ENTRIES OF THE BASIS MATRIX AND
!     DECOMPOSING IT USING THE LA05 PACKAGE.
!     IT IS THE MAIN PART OF THE PROCEDURE (DECOMPOSE BASIS MATRIX).
!
!***SEE ALSO  DSPLP
!***ROUTINES CALLED  DASUM, DPNNZR, LA05AD, XERMSG
!***COMMON BLOCKS    LA05DD
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890605  Added DASUM to list of DOUBLE PRECISION variables.
!   890605  Removed unreferenced labels.  (WRB)
!   891009  Removed unreferenced variable.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls, convert do-it-yourself
!           DO loops to DO loops.  (RWC)
!***END PROLOGUE  DPLPDM
  INTEGER IBASIS(*),IMAT(*),IBRC(LBM,2),IPR(*),IWR(*),IND(*),IBB(*)
  DOUBLE PRECISION AIJ,AMAT(*),BASMAT(*),CSC(*),WR(*),ANORM,DASUM, &
   EPS,GG,ONE,SMALL,UU,ZERO
  LOGICAL SINGLR,REDBAS
  CHARACTER*16 XERN3
!
!     COMMON BLOCK USED BY LA05 () PACKAGE..
  COMMON /LA05DD/ SMALL,LP,LENL,LENU,NCP,LROW,LCOL
!
!***FIRST EXECUTABLE STATEMENT  DPLPDM
  ZERO = 0.D0
  ONE = 1.D0
!
!     DEFINE BASIS MATRIX BY COLUMNS FOR SPARSE MATRIX EQUATION SOLVER.
!     THE LA05AD() SUBPROGRAM REQUIRES THE NONZERO ENTRIES OF THE MATRIX
!     TOGETHER WITH THE ROW AND COLUMN INDICES.
!
  NZBM = 0
!
!     DEFINE DEPENDENT VARIABLE COLUMNS. THESE ARE
!     COLS. OF THE IDENTITY MATRIX AND IMPLICITLY GENERATED.
!
  DO 20 K = 1,MRELAS
     J = IBASIS(K)
     if (J > NVARS) THEN
        NZBM = NZBM+1
        if (IND(J) == 2) THEN
           BASMAT(NZBM) = ONE
        ELSE
           BASMAT(NZBM) = -ONE
        ENDIF
        IBRC(NZBM,1) = J-NVARS
        IBRC(NZBM,2) = K
     ELSE
!
!           DEFINE THE INDEP. VARIABLE COLS.  THIS REQUIRES RETRIEVING
!           THE COLS. FROM THE SPARSE MATRIX DATA STRUCTURE.
!
        I = 0
   10       call DPNNZR(I,AIJ,IPLACE,AMAT,IMAT,J)
        if (I > 0) THEN
           NZBM = NZBM+1
           BASMAT(NZBM) = AIJ*CSC(J)
           IBRC(NZBM,1) = I
           IBRC(NZBM,2) = K
           go to 10
        ENDIF
     ENDIF
   20 CONTINUE
!
  SINGLR = .FALSE.
!
!     RECOMPUTE MATRIX NORM USING CRUDE NORM  =  SUM OF MAGNITUDES.
!
  ANORM = DASUM(NZBM,BASMAT,1)
  SMALL = EPS*ANORM
!
!     GET AN L-U FACTORIZATION OF THE BASIS MATRIX.
!
  NREDC = NREDC+1
  REDBAS = .TRUE.
  call LA05AD(BASMAT,IBRC,NZBM,LBM,MRELAS,IPR,IWR,WR,GG,UU)
!
!     CHECK RETURN VALUE OF ERROR FLAG, GG.
!
  if (GG >= ZERO) RETURN
  if (GG == (-7.)) THEN
     call XERMSG ('SLATEC', 'DPLPDM', &
        'IN DSPLP, SHORT ON STORAGE FOR LA05AD.  ' // &
        'USE PRGOPT(*) TO GIVE MORE.', 28, IOPT)
     INFO = -28
  ELSEIF (GG == (-5.)) THEN
     SINGLR = .TRUE.
  ELSE
     WRITE (XERN3, '(1PE15.6)') GG
     call XERMSG ('SLATEC', 'DPLPDM', &
        'IN DSPLP, LA05AD RETURNED ERROR FLAG = ' // XERN3, &
        27, IOPT)
     INFO = -27
  end if
  return
end
