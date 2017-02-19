subroutine WNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM, &
     IDOPE, DOPE, DONE)
!
!! WNLIT is subsidiary to WNNLS.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (WNLIT-S, DWNLIT-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     This is a companion subprogram to WNNLS( ).
!     The documentation for WNNLS( ) has complete usage instructions.
!
!     Note  The M by (N+1) matrix W( , ) contains the rt. hand side
!           B as the (N+1)st col.
!
!     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with
!     col interchanges.
!
!***SEE ALSO  WNNLS
!***ROUTINES CALLED  H12, ISAMAX, SCOPY, SROTM, SROTMG, SSCAL, SSWAP,
!                    WNLT1, WNLT2, WNLT3
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and revised.  (WRB & RWC)
!   890620  Revised to make WNLT1, WNLT2, and WNLT3 subroutines.  (RWC)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  WNLIT
  INTEGER IDOPE(*), IPIVOT(*), ITYPE(*), L, M, MDW, N
  REAL             DOPE(*), H(*), RNORM, SCALE(*), W(MDW,*)
  LOGICAL DONE
!
  EXTERNAL H12, ISAMAX, SCOPY, SROTM, SROTMG, SSCAL, SSWAP, WNLT1, &
     WNLT2, WNLT3
  INTEGER ISAMAX
  LOGICAL WNLT2
!
  REAL             ALSQ, AMAX, EANORM, FACTOR, HBAR, RN, SPARAM(5), &
     T, TAU
  INTEGER I, I1, IMAX, IR, J, J1, JJ, JP, KRANK, L1, LB, LEND, ME, &
     MEND, NIV, NSOLN
  LOGICAL INDEP, RECALC
!
!***FIRST EXECUTABLE STATEMENT  WNLIT
  ME    = IDOPE(1)
  NSOLN = IDOPE(2)
  L1    = IDOPE(3)
!
  ALSQ   = DOPE(1)
  EANORM = DOPE(2)
  TAU    = DOPE(3)
!
  LB     = MIN(M-1,L)
  RECALC = .TRUE.
  RNORM  = 0.E0
  KRANK  = 0
!
!     We set FACTOR=1.0 so that the heavy weight ALAMDA will be
!     included in the test for column independence.
!
  FACTOR = 1.E0
  LEND = L
  DO 180 I=1,LB
!
!        Set IR to point to the I-th row.
!
     IR = I
     MEND = M
     call WNLT1 (I, LEND, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE, &
                  W)
!
!        Update column SS and find pivot column.
!
     call WNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
!
!        Perform column interchange.
!        Test independence of incoming column.
!
  130    if (WNLT2(ME, MEND, IR, FACTOR, TAU, SCALE, W(1,I))) THEN
!
!           Eliminate I-th column below diagonal using modified Givens
!           transformations applied to (A B).
!
!           When operating near the ME line, use the largest element
!           above it as the pivot.
!
        DO 160 J=M,I+1,-1
           JP = J-1
           if (J == ME+1) THEN
              IMAX = ME
              AMAX = SCALE(ME)*W(ME,I)**2
              DO 150 JP=J-1,I,-1
                 T = SCALE(JP)*W(JP,I)**2
                 if (T > AMAX) THEN
                    IMAX = JP
                    AMAX = T
                 ENDIF
  150             CONTINUE
              JP = IMAX
           ENDIF
!
           if (W(J,I) /= 0.E0) THEN
              call SROTMG (SCALE(JP), SCALE(J), W(JP,I), W(J,I), &
                           SPARAM)
              W(J,I) = 0.E0
              call SROTM (N+1-I, W(JP,I+1), MDW, W(J,I+1), MDW, &
                          SPARAM)
           ENDIF
  160       CONTINUE
     ELSE if (LEND > I) THEN
!
!           Column I is dependent.  Swap with column LEND.
!           Perform column interchange,
!           and find column in remaining set with largest SS.
!
        call WNLT3 (I, LEND, M, MDW, IPIVOT, H, W)
        LEND = LEND - 1
        IMAX = ISAMAX(LEND-I+1, H(I), 1) + I - 1
        HBAR = H(IMAX)
        go to 130
     ELSE
        KRANK = I - 1
        go to 190
     ENDIF
  180 CONTINUE
  KRANK = L1
!
  190 if (KRANK < ME) THEN
     FACTOR = ALSQ
     DO 200 I=KRANK+1,ME
        call sinit ( L, 0.E0, W(I,1), MDW)
  200    CONTINUE
!
!        Determine the rank of the remaining equality constraint
!        equations by eliminating within the block of constrained
!        variables.  Remove any redundant constraints.
!
     RECALC = .TRUE.
     LB = MIN(L+ME-KRANK, N)
     DO 270 I=L+1,LB
        IR = KRANK + I - L
        LEND = N
        MEND = ME
        call WNLT1 (I, LEND, ME, IR, MDW, RECALC, IMAX, HBAR, H, &
                     SCALE, W)
!
!           Update col ss and find pivot col
!
        call WNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
!
!           Perform column interchange
!           Eliminate elements in the I-th col.
!
        DO 240 J=ME,IR+1,-1
           if (W(J,I) /= 0.E0) THEN
             call SROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I), &
                          SPARAM)
              W(J,I) = 0.E0
              call SROTM (N+1-I, W(J-1,I+1), MDW,W(J,I+1), MDW, &
                          SPARAM)
           ENDIF
  240       CONTINUE
!
!           I=column being eliminated.
!           Test independence of incoming column.
!           Remove any redundant or dependent equality constraints.
!
        if (.NOT.WNLT2(ME, MEND, IR, FACTOR,TAU,SCALE,W(1,I))) THEN
           JJ = IR
           DO 260 IR=JJ,ME
              call sinit ( N, 0.E0, W(IR,1), MDW)
              RNORM = RNORM + (SCALE(IR)*W(IR,N+1)/ALSQ)*W(IR,N+1)
              W(IR,N+1) = 0.E0
              SCALE(IR) = 1.E0
!
!                 Reclassify the zeroed row as a least squares equation.
!
              ITYPE(IR) = 1
  260          CONTINUE
!
!              Reduce ME to reflect any discovered dependent equality
!              constraints.
!
           ME = JJ - 1
           go to 280
        ENDIF
  270    CONTINUE
  end if
!
!     Try to determine the variables KRANK+1 through L1 from the
!     least squares equations.  Continue the triangularization with
!     pivot element W(ME+1,I).
!
  280 if (KRANK < L1) THEN
     RECALC = .TRUE.
!
!        Set FACTOR=ALSQ to remove effect of heavy weight from
!        test for column independence.
!
     FACTOR = ALSQ
     DO 350 I=KRANK+1,L1
!
!           Set IR to point to the ME+1-st row.
!
        IR = ME+1
        LEND = L
        MEND = M
        call WNLT1 (I, L, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE, &
                     W)
!
!           Update column SS and find pivot column.
!
        call WNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
!
!           Perform column interchange.
!           Eliminate I-th column below the IR-th element.
!
        DO 320 J=M,IR+1,-1
           if (W(J,I) /= 0.E0) THEN
             call SROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I), &
                          SPARAM)
              W(J,I) = 0.E0
              call SROTM (N+1-I, W(J-1,I+1),  MDW, W(J,I+1), MDW, &
                          SPARAM)
           ENDIF
  320       CONTINUE
!
!           Test if new pivot element is near zero.
!           If so, the column is dependent.
!           Then check row norm test to be classified as independent.
!
        T = SCALE(IR)*W(IR,I)**2
        INDEP = T  >  (TAU*EANORM)**2
        if (INDEP) THEN
           RN = 0.E0
           DO 340 I1=IR,M
              DO 330 J1=I+1,N
                 RN = MAX(RN, SCALE(I1)*W(I1,J1)**2)
  330             CONTINUE
  340          CONTINUE
           INDEP = T  >  RN*TAU**2
        ENDIF
!
!           If independent, swap the IR-th and KRANK+1-th rows to
!           maintain the triangular form.  Update the rank indicator
!           KRANK and the equality constraint pointer ME.
!
        if (.NOT.INDEP) go to 360
        call SSWAP(N+1, W(KRANK+1,1), MDW, W(IR,1), MDW)
        call SSWAP(1, SCALE(KRANK+1), 1, SCALE(IR), 1)
!
!           Reclassify the least square equation as an equality
!           constraint and rescale it.
!
        ITYPE(IR) = 0
        T = SQRT(SCALE(KRANK+1))
        call SSCAL(N+1, T, W(KRANK+1,1), MDW)
        SCALE(KRANK+1) = ALSQ
        ME = ME+1
        KRANK = KRANK+1
  350    CONTINUE
  end if
!
!     If pseudorank is less than L, apply Householder transformation.
!     from right.
!
  360 if (KRANK < L) THEN
     DO 370 J=KRANK,1,-1
        call H12 (1, J, KRANK+1, L, W(J,1), MDW, H(J), W, MDW, 1, &
                  J-1)
  370    CONTINUE
  end if
!
  NIV = KRANK + NSOLN - L
  if (L == N) DONE = .TRUE.
!
!     End of initial triangularization.
!
  IDOPE(1) = ME
  IDOPE(2) = KRANK
  IDOPE(3) = NIV
  return
end
