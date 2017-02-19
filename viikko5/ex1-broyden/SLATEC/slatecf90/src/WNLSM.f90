subroutine WNLSM (W, MDW, MME, MA, N, L, PRGOPT, X, RNORM, MODE, &
     IPIVOT, ITYPE, WD, H, SCALE, Z, TEMP, D)
!
!! WNLSM is subsidiary to WNNLS
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (WNLSM-S, DWNLSM-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     This is a companion subprogram to WNNLS.
!     The documentation for WNNLS has complete usage instructions.
!
!     In addition to the parameters discussed in the prologue to
!     subroutine WNNLS, the following work arrays are used in
!     subroutine WNLSM  (they are passed through the calling
!     sequence from WNNLS for purposes of variable dimensioning).
!     Their contents will in general be of no interest to the user.
!
!         IPIVOT(*)
!            An array of length N.  Upon completion it contains the
!         pivoting information for the cols of W(*,*).
!
!         ITYPE(*)
!            An array of length M which is used to keep track
!         of the classification of the equations.  ITYPE(I)=0
!         denotes equation I as an equality constraint.
!         ITYPE(I)=1 denotes equation I as a least squares
!         equation.
!
!         WD(*)
!            An array of length N.  Upon completion it contains the
!         dual solution vector.
!
!         H(*)
!            An array of length N.  Upon completion it contains the
!         pivot scalars of the Householder transformations performed
!         in the case KRANK < L.
!
!         SCALE(*)
!            An array of length M which is used by the subroutine
!         to store the diagonal matrix of weights.
!         These are used to apply the modified Givens
!         transformations.
!
!         Z(*),TEMP(*)
!            Working arrays of length N.
!
!         D(*)
!            An array of length N that contains the
!         column scaling for the matrix (E).
!                                       (A)
!
!***SEE ALSO  WNNLS
!***ROUTINES CALLED  H12, ISAMAX, R1MACH, SASUM, SAXPY, SCOPY, SNRM2,
!                    SROTM, SROTMG, SSCAL, SSWAP, WNLIT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and revised.  (WRB & RWC)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   900510  Fixed an error message.  (RWC)
!***END PROLOGUE  WNLSM
  INTEGER IPIVOT(*), ITYPE(*), L, MA, MDW, MME, MODE, N
  REAL             D(*), H(*), PRGOPT(*), RNORM, SCALE(*), TEMP(*), &
     W(MDW,*), WD(*), X(*), Z(*)
!
  EXTERNAL H12, ISAMAX, R1MACH, SASUM, SAXPY, SCOPY, SNRM2, SROTM, &
     SROTMG, SSCAL, SSWAP, WNLIT, XERMSG
  REAL             R1MACH, SASUM, SNRM2
  INTEGER ISAMAX
!
  REAL             ALAMDA, ALPHA, ALSQ, AMAX, BLOWUP, BNORM, &
     DOPE(3), EANORM, FAC, SM, SPARAM(5), SRELPR, T, TAU, WMAX, Z2, &
     ZZ
  INTEGER I, IDOPE(3), IMAX, ISOL, ITEMP, ITER, ITMAX, IWMAX, J, &
     JCON, JP, KEY, KRANK, L1, LAST, LINK, M, ME, NEXT, NIV, NLINK, &
     NOPT, NSOLN, NTIMES
  LOGICAL DONE, FEASBL, FIRST, HITCON, POS
!
  SAVE SRELPR, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  WNLSM
!
!     Initialize variables.
!     SRELPR is the precision for the particular machine
!     being used.  This logic avoids resetting it every entry.
!
  if (FIRST) SRELPR = R1MACH(4)
  FIRST = .FALSE.
!
!     Set the nominal tolerance used in the code.
!
  TAU = SQRT(SRELPR)
!
  M = MA + MME
  ME = MME
  MODE = 2
!
!     To process option vector
!
  FAC = 1.E-4
!
!     Set the nominal blow up factor used in the code.
!
  BLOWUP = TAU
!
!     The nominal column scaling used in the code is
!     the identity scaling.
!
  call sinit ( N, 1.E0, D, 1)
!
!     Define bound for number of options to change.
!
  NOPT = 1000
!
!     Define bound for positive value of LINK.
!
  NLINK = 100000
  NTIMES = 0
  LAST = 1
  LINK = PRGOPT(1)
  if (LINK <= 0 .OR. LINK > NLINK) THEN
     call XERMSG ('SLATEC', 'WNLSM', &
        'WNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
     return
  end if
!
  100 if (LINK > 1) THEN
     NTIMES = NTIMES + 1
     if (NTIMES > NOPT) THEN
     call XERMSG ('SLATEC', 'WNLSM', &
        'WNNLS, THE LINKS IN THE OPTION VECTOR ARE CYCLING.', 3, 1)
        return
     ENDIF
!
     KEY = PRGOPT(LAST+1)
     if (KEY == 6 .AND. PRGOPT(LAST+2) /= 0.E0) THEN
        DO 110 J = 1,N
           T = SNRM2(M,W(1,J),1)
           if (T /= 0.E0) T = 1.E0/T
           D(J) = T
  110       CONTINUE
     ENDIF
!
     if (KEY == 7) call SCOPY (N, PRGOPT(LAST+2), 1, D, 1)
     if (KEY == 8) TAU = MAX(SRELPR,PRGOPT(LAST+2))
     if (KEY == 9) BLOWUP = MAX(SRELPR,PRGOPT(LAST+2))
!
     NEXT = PRGOPT(LINK)
     if (NEXT <= 0 .OR. NEXT > NLINK) THEN
        call XERMSG ('SLATEC', 'WNLSM', &
           'WNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
        return
     ENDIF
!
     LAST = LINK
     LINK = NEXT
     go to 100
  end if
!
  DO 120 J = 1,N
     call SSCAL (M, D(J), W(1,J), 1)
  120 CONTINUE
!
!     Process option vector
!
  DONE = .FALSE.
  ITER = 0
  ITMAX = 3*(N-L)
  MODE = 0
  NSOLN = L
  L1 = MIN(M,L)
!
!     Compute scale factor to apply to equality constraint equations.
!
  DO 130 J = 1,N
     WD(J) = SASUM(M,W(1,J),1)
  130 CONTINUE
!
  IMAX = ISAMAX(N,WD,1)
  EANORM = WD(IMAX)
  BNORM = SASUM(M,W(1,N+1),1)
  ALAMDA = EANORM/(SRELPR*FAC)
!
!     Define scaling diagonal matrix for modified Givens usage and
!     classify equation types.
!
  ALSQ = ALAMDA**2
  DO 140 I = 1,M
!
!        When equation I is heavily weighted ITYPE(I)=0,
!        else ITYPE(I)=1.
!
     if (I <= ME) THEN
        T = ALSQ
        ITEMP = 0
     ELSE
        T = 1.E0
        ITEMP = 1
     ENDIF
     SCALE(I) = T
     ITYPE(I) = ITEMP
  140 CONTINUE
!
!     Set the solution vector X(*) to zero and the column interchange
!     matrix to the identity.
!
  call sinit ( N, 0.E0, X, 1)
  DO 150 I = 1,N
     IPIVOT(I) = I
  150 CONTINUE
!
!     Perform initial triangularization in the submatrix
!     corresponding to the unconstrained variables.
!     Set first L components of dual vector to zero because
!     these correspond to the unconstrained variables.
!
  call sinit ( L, 0.E0, WD, 1)
!
!     The arrays IDOPE(*) and DOPE(*) are used to pass
!     information to WNLIT().  This was done to avoid
!     a long calling sequence or the use of COMMON.
!
  IDOPE(1) = ME
  IDOPE(2) = NSOLN
  IDOPE(3) = L1
!
  DOPE(1) = ALSQ
  DOPE(2) = EANORM
  DOPE(3) = TAU
  call WNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM, &
              IDOPE, DOPE, DONE)
  ME    = IDOPE(1)
  KRANK = IDOPE(2)
  NIV   = IDOPE(3)
!
!     Perform WNNLS algorithm using the following steps.
!
!     Until(DONE)
!        compute search direction and feasible point
!        when (HITCON) add constraints
!        else perform multiplier test and drop a constraint
!        fin
!     Compute-Final-Solution
!
!     To compute search direction and feasible point,
!     solve the triangular system of currently non-active
!     variables and store the solution in Z(*).
!
!     To solve system
!     Copy right hand side into TEMP vector to use overwriting method.
!
  160 if (DONE) go to 330
  ISOL = L + 1
  if (NSOLN >= ISOL) THEN
     call SCOPY (NIV, W(1,N+1), 1, TEMP, 1)
     DO 170 J = NSOLN,ISOL,-1
        if (J > KRANK) THEN
           I = NIV - NSOLN + J
        ELSE
           I = J
        ENDIF
!
        if (J > KRANK .AND. J <= L) THEN
           Z(J) = 0.E0
        ELSE
           Z(J) = TEMP(I)/W(I,J)
           call SAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
        ENDIF
  170    CONTINUE
  end if
!
!     Increment iteration counter and check against maximum number
!     of iterations.
!
  ITER = ITER + 1
  if (ITER > ITMAX) THEN
     MODE = 1
     DONE = .TRUE.
  end if
!
!     Check to see if any constraints have become active.
!     If so, calculate an interpolation factor so that all
!     active constraints are removed from the basis.
!
  ALPHA = 2.E0
  HITCON = .FALSE.
  DO 180 J = L+1,NSOLN
     ZZ = Z(J)
     if (ZZ <= 0.E0) THEN
        T = X(J)/(X(J)-ZZ)
        if (T < ALPHA) THEN
           ALPHA = T
           JCON = J
        ENDIF
        HITCON = .TRUE.
     ENDIF
  180 CONTINUE
!
!     Compute search direction and feasible point
!
  if (HITCON) THEN
!
!        To add constraints, use computed ALPHA to interpolate between
!        last feasible solution X(*) and current unconstrained (and
!        infeasible) solution Z(*).
!
     DO 190 J = L+1,NSOLN
        X(J) = X(J) + ALPHA*(Z(J)-X(J))
  190    CONTINUE
     FEASBL = .FALSE.
!
!        Remove column JCON and shift columns JCON+1 through N to the
!        left.  Swap column JCON into the N th position.  This achieves
!        upper Hessenberg form for the nonactive constraints and
!        leaves an upper Hessenberg matrix to retriangularize.
!
  200    DO 210 I = 1,M
        T = W(I,JCON)
        call SCOPY (N-JCON, W(I, JCON+1), MDW, W(I, JCON), MDW)
        W(I,N) = T
  210    CONTINUE
!
!        Update permuted index vector to reflect this shift and swap.
!
     ITEMP = IPIVOT(JCON)
     DO 220 I = JCON,N - 1
        IPIVOT(I) = IPIVOT(I+1)
  220    CONTINUE
     IPIVOT(N) = ITEMP
!
!        Similarly permute X(*) vector.
!
     call SCOPY (N-JCON, X(JCON+1), 1, X(JCON), 1)
     X(N) = 0.E0
     NSOLN = NSOLN - 1
     NIV = NIV - 1
!
!        Retriangularize upper Hessenberg matrix after adding
!        constraints.
!
     I = KRANK + JCON - L
     DO 230 J = JCON,NSOLN
        if (ITYPE(I) == 0 .AND. ITYPE(I+1) == 0) THEN
!
!              Zero IP1 to I in column J
!
           if (W(I+1,J) /= 0.E0) THEN
              call SROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J), &
                           SPARAM)
              W(I+1,J) = 0.E0
              call SROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW, &
                          SPARAM)
           ENDIF
        ELSEIF (ITYPE(I) == 1 .AND. ITYPE(I+1) == 1) THEN
!
!              Zero IP1 to I in column J
!
           if (W(I+1,J) /= 0.E0) THEN
              call SROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J), &
                           SPARAM)
              W(I+1,J) = 0.E0
              call SROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW, &
                          SPARAM)
           ENDIF
        ELSEIF (ITYPE(I) == 1 .AND. ITYPE(I+1) == 0) THEN
           call SSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
           call SSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
           ITEMP = ITYPE(I+1)
           ITYPE(I+1) = ITYPE(I)
           ITYPE(I) = ITEMP
!
!              Swapped row was formerly a pivot element, so it will
!              be large enough to perform elimination.
!              Zero IP1 to I in column J.
!
           if (W(I+1,J) /= 0.E0) THEN
              call SROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J), &
                           SPARAM)
              W(I+1,J) = 0.E0
              call SROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW, &
                          SPARAM)
           ENDIF
        ELSEIF (ITYPE(I) == 0 .AND. ITYPE(I+1) == 1) THEN
           if (SCALE(I)*W(I,J)**2/ALSQ > (TAU*EANORM)**2) THEN
!
!                 Zero IP1 to I in column J
!
              if (W(I+1,J) /= 0.E0) THEN
                 call SROTMG (SCALE(I), SCALE(I+1), W(I,J), &
                              W(I+1,J), SPARAM)
                 W(I+1,J) = 0.E0
                 call SROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW, &
                             SPARAM)
              ENDIF
           ELSE
              call SSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
              call SSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
              ITEMP = ITYPE(I+1)
              ITYPE(I+1) = ITYPE(I)
              ITYPE(I) = ITEMP
              W(I+1,J) = 0.E0
           ENDIF
        ENDIF
        I = I + 1
  230    CONTINUE
!
!        See if the remaining coefficients in the solution set are
!        feasible.  They should be because of the way ALPHA was
!        determined.  If any are infeasible, it is due to roundoff
!        error.  Any that are non-positive will be set to zero and
!        removed from the solution set.
!
     DO 240 JCON = L+1,NSOLN
        if (X(JCON) <= 0.E0) go to 250
  240    CONTINUE
     FEASBL = .TRUE.
  250    if (.NOT.FEASBL) go to 200
  ELSE
!
!        To perform multiplier test and drop a constraint.
!
     call SCOPY (NSOLN, Z, 1, X, 1)
     if (NSOLN < N) call sinit ( N-NSOLN, 0.E0, X(NSOLN+1), 1)
!
!        Reclassify least squares equations as equalities as necessary.
!
     I = NIV + 1
  260    if (I <= ME) THEN
        if (ITYPE(I) == 0) THEN
           I = I + 1
        ELSE
           call SSWAP (N+1, W(I,1), MDW, W(ME,1), MDW)
           call SSWAP (1, SCALE(I), 1, SCALE(ME), 1)
           ITEMP = ITYPE(I)
           ITYPE(I) = ITYPE(ME)
           ITYPE(ME) = ITEMP
           ME = ME - 1
        ENDIF
        go to 260
     ENDIF
!
!        Form inner product vector WD(*) of dual coefficients.
!
     DO 280 J = NSOLN+1,N
        SM = 0.E0
        DO 270 I = NSOLN+1,M
           SM = SM + SCALE(I)*W(I,J)*W(I,N+1)
  270       CONTINUE
        WD(J) = SM
  280    CONTINUE
!
!        Find J such that WD(J)=WMAX is maximum.  This determines
!        that the incoming column J will reduce the residual vector
!        and be positive.
!
  290    WMAX = 0.E0
     IWMAX = NSOLN + 1
     DO 300 J = NSOLN+1,N
        if (WD(J) > WMAX) THEN
           WMAX = WD(J)
           IWMAX = J
        ENDIF
  300    CONTINUE
     if (WMAX <= 0.E0) go to 330
!
!        Set dual coefficients to zero for incoming column.
!
     WD(IWMAX) = 0.E0
!
!        WMAX  >  0.E0, so okay to move column IWMAX to solution set.
!        Perform transformation to retriangularize, and test for near
!        linear dependence.
!
!        Swap column IWMAX into NSOLN-th position to maintain upper
!        Hessenberg form of adjacent columns, and add new column to
!        triangular decomposition.
!
     NSOLN = NSOLN + 1
     NIV = NIV + 1
     if (NSOLN /= IWMAX) THEN
        call SSWAP (M, W(1,NSOLN), 1, W(1,IWMAX), 1)
        WD(IWMAX) = WD(NSOLN)
        WD(NSOLN) = 0.E0
        ITEMP = IPIVOT(NSOLN)
        IPIVOT(NSOLN) = IPIVOT(IWMAX)
        IPIVOT(IWMAX) = ITEMP
     ENDIF
!
!        Reduce column NSOLN so that the matrix of nonactive constraints
!        variables is triangular.
!
     DO 320 J = M,NIV+1,-1
        JP = J - 1
!
!           When operating near the ME line, test to see if the pivot
!           element is near zero.  If so, use the largest element above
!           it as the pivot.  This is to maintain the sharp interface
!           between weighted and non-weighted rows in all cases.
!
        if (J == ME+1) THEN
           IMAX = ME
           AMAX = SCALE(ME)*W(ME,NSOLN)**2
           DO 310 JP = J - 1,NIV,-1
              T = SCALE(JP)*W(JP,NSOLN)**2
              if (T > AMAX) THEN
                 IMAX = JP
                 AMAX = T
              ENDIF
  310          CONTINUE
           JP = IMAX
        ENDIF
!
        if (W(J,NSOLN) /= 0.E0) THEN
           call SROTMG (SCALE(JP), SCALE(J), W(JP,NSOLN), &
                        W(J,NSOLN), SPARAM)
           W(J,NSOLN) = 0.E0
           call SROTM (N+1-NSOLN, W(JP,NSOLN+1), MDW, W(J,NSOLN+1), &
                       MDW, SPARAM)
        ENDIF
  320    CONTINUE
!
!        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
!        this is nonpositive or too large.  If this was true or if the
!        pivot term was zero, reject the column as dependent.
!
     if (W(NIV,NSOLN) /= 0.E0) THEN
        ISOL = NIV
        Z2 = W(ISOL,N+1)/W(ISOL,NSOLN)
        Z(NSOLN) = Z2
        POS = Z2  >  0.E0
        if (Z2*EANORM >= BNORM .AND. POS) THEN
           POS = .NOT. (BLOWUP*Z2*EANORM >= BNORM)
        ENDIF
!
!           Try to add row ME+1 as an additional equality constraint.
!           Check size of proposed new solution component.
!           Reject it if it is too large.
!
     ELSEIF (NIV <= ME .AND. W(ME+1,NSOLN) /= 0.E0) THEN
        ISOL = ME + 1
        if (POS) THEN
!
!              Swap rows ME+1 and NIV, and scale factors for these rows.
!
           call SSWAP (N+1, W(ME+1,1), MDW, W(NIV,1), MDW)
           call SSWAP (1, SCALE(ME+1), 1, SCALE(NIV), 1)
           ITEMP = ITYPE(ME+1)
           ITYPE(ME+1) = ITYPE(NIV)
           ITYPE(NIV) = ITEMP
           ME = ME + 1
        ENDIF
     ELSE
        POS = .FALSE.
     ENDIF
!
     if (.NOT.POS) THEN
        NSOLN = NSOLN - 1
        NIV = NIV - 1
     ENDIF
     if (.NOT.(POS.OR.DONE)) go to 290
  end if
  go to 160
!
!     Else perform multiplier test and drop a constraint.  To compute
!     final solution.  Solve system, store results in X(*).
!
!     Copy right hand side into TEMP vector to use overwriting method.
!
  330 ISOL = 1
  if (NSOLN >= ISOL) THEN
     call SCOPY (NIV, W(1,N+1), 1, TEMP, 1)
     DO 340 J = NSOLN,ISOL,-1
        if (J > KRANK) THEN
           I = NIV - NSOLN + J
        ELSE
           I = J
        ENDIF
!
        if (J > KRANK .AND. J <= L) THEN
           Z(J) = 0.E0
        ELSE
           Z(J) = TEMP(I)/W(I,J)
           call SAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
        ENDIF
  340    CONTINUE
  end if
!
!     Solve system.
!
  call SCOPY (NSOLN, Z, 1, X, 1)
!
!     Apply Householder transformations to X(*) if KRANK < L
!
  if (KRANK < L) THEN
     DO 350 I = 1,KRANK
        call H12 (2, I, KRANK+1, L, W(I,1), MDW, H(I), X, 1, 1, 1)
  350    CONTINUE
  end if
!
!     Fill in trailing zeroes for constrained variables not in solution.
!
  if (NSOLN < N) call sinit ( N-NSOLN, 0.E0, X(NSOLN+1), 1)
!
!     Permute solution vector to natural order.
!
  DO 380 I = 1,N
     J = I
  360    if (IPIVOT(J) == I) go to 370
     J = J + 1
     go to 360
!
  370    IPIVOT(J) = IPIVOT(I)
     IPIVOT(I) = J
     call SSWAP (1, X(J), 1, X(I), 1)
  380 CONTINUE
!
!     Rescale the solution using the column scaling.
!
  DO 390 J = 1,N
     X(J) = X(J)*D(J)
  390 CONTINUE
!
  DO 400 I = NSOLN+1,M
     T = W(I,N+1)
     if (I <= ME) T = T/ALAMDA
     T = (SCALE(I)*T)*T
     RNORM = RNORM + T
  400 CONTINUE
!
  RNORM = SQRT(RNORM)
  return
end
