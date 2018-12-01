    MODULE SetEqn
    USE IntrType
    USE Data
    USE Eos
    USE ScalarDat
    CONTAINS

    SUBROUTINE SetCoeffs
!
!   Set coefficients in the finite volume conduction
!
!                 Programmed by John Mahaffy   1/01
!
    IMPLICIT NONE

    REAL(sdk) :: cols(5), coeffs(5), cxm, cxp, cyp, cym
    INTEGER(sik) :: i, j, k, ieq, ieqxp, ieqxm, ieqym, ieqyp, ii, ic, ieqt
    INTEGER(sik) :: nCoef,iDiag
!
    IF (numStepsTried .GT. 0) THEN
       d = cnd*delt/(cpm*rhom)  !  coefficient of diffusion per time step
       pm = qm*delt/(cpm*rhom)    !  temperature change per step due to power source
    ELSE
       d = cnd
       pm = qm
    END IF
!
!   Load the sparse version of the matrix
!
!                   Set the location of the Diagonal Coefficient
    as%iDiag = 3
    as(1:nx)%iDiag = as(1:nx)%iDiag - 1   ! one less lower index coefficient
    as(1:n:nx)%iDiag = as(1:n:nx)%iDiag - 1
!
!   Count the number of non-zero coefficients in each row
!
    as(1:n)%nCoef = 5
    as(1:nx)%nCoef = as(1:nx)%nCoef -1            !  Bottom  Boundary
    as(n-nx+1:n)%nCoef = as(n-nx+1:n)%nCoef -1    !  Top  Boundary
    as(1:n:nx)%nCoef = as(1:n:nx)%nCoef -1        !  Left Boundary
    as(nx:n:nx)%nCoef = as(nx:n:nx)%nCoef - 1     !  Right Boundary
!
    DO i = 1,n
      ii = as(i)%nCoef
      ALLOCATE(as(i)%a(ii), as(i)%ax(ii), as(i)%iCol(ii), as(i)%iDim(ii))
    ENDDO
!                    Initialize information on the diagonal
    DO ieq = 1,n
      iDiag = as(ieq)%iDiag
      as(ieq)%a(iDiag) = 0.0_sdk
      as(ieq)%iCol(iDiag) = ieq
      as(ieq)%iDim(iDiag) = 0
    ENDDO
!                    Set information on for coupling in the lower x direction
    DO j = 1,ny
      DO i = 2,nx
        ieq = (j-1)*nx + i
        iDiag = as(ieq)%iDiag
        as(ieq)%a(iDiag-1) = -d/dxm**2
        as(ieq)%a(iDiag) = as(ieq)%a(iDiag) + d/dxm**2
        as(ieq)%iDim(iDiag-1) = 1
        as(ieq)%iCol(iDiag-1) = ieq - 1
      ENDDO
    ENDDO

!                    Set information on for coupling in the upper x direction
    DO j = 1,ny
      DO i = 1,nx-1
        ieq = (j-1)*nx + i
        iDiag = as(ieq)%iDiag
        as(ieq)%a(iDiag+1) = -d/dxm**2
        as(ieq)%a(iDiag) = as(ieq)%a(iDiag) + d/dxm**2
        as(ieq)%iDim(iDiag+1) = 1
        as(ieq)%iCol(iDiag+1) = ieq + 1
      ENDDO
    ENDDO
!                    Set information on for coupling in the lower y direction
    DO j = 2,ny
      DO i = 1,nx
        ieq = (j-1)*nx + i
        r_o = r_i + (j-0.5_sdk)*dy
        iDiag = as(ieq)%iDiag
        as(ieq)%a(1) = -d*(r_o - dy*0.5_sdk)/(r_o*dy**2)
        as(ieq)%a(iDiag) = as(ieq)%a(iDiag) + d*(r_o-dy*0.5_sdk)/(r_o*dy**2) 
        as(ieq)%iDim(1) = 2
        as(ieq)%iCol(1) = ieq - nx
      ENDDO
    ENDDO
!                    Set information on for coupling in the upper y direction
    DO j = 1,ny-1
      DO i = 1,nx
        ieq = (j-1)*nx + i
        r_o = r_i + (j-0.5_sdk)*dy
        iDiag = as(ieq)%iDiag
        ii = as(ieq)%nCoef
        as(ieq)%a(ii) = -d*(r_o + dy*0.5_sdk)/(r_o*dy**2)
        as(ieq)%a(iDiag) = as(ieq)%a(iDiag) + d*(r_o + dy*0.5_sdk)/(r_o*dy**2)
        as(ieq)%iDim(ii) = 2
        as(ieq)%iCol(ii) = ieq + nx
      ENDDO
    ENDDO
!
! Adjustments for the y face boundary conditions
!
    DO ieq=1,nx
      iDiag = as(ieq)%iDiag
      r_o = r_i + (0.5_sdk)*dy
      as(ieq)%a(iDiag) = as(ieq)%a(iDiag)+9*d*hy*(r_o-dy*0.5_sdk)/(r_o*dy*(3*hy*dy+8*cnd)) !Bottom
      ii = as(ieq)%nCoef
      as(ieq)%a(ii) = as(ieq)%a(ii)-d*hy*(r_o-dy*0.5_sdk)/(r_o*dy*(3*hy*dy+8*cnd)) !Bottom
      ieqt = ieq + n - nx
      iDiag = as(ieqt)%iDiag
      as(ieqt)%a(iDiag) = as(ieqt)%a(iDiag)!+9*d*hy/(dy*(3*hy*dy+8*cnd)) !Top
      as(ieqt)%a(1) = as(ieqt)%a(1)!-d*hy/(dy*(3*hy*dy+8*cnd)) !Top
    ENDDO
!
! Adjustments for the x face boundary conditions
!
    DO ieq=1,n,nx
      iDiag = as(ieq)%iDiag
      ii = iDiag + 1
      as(ieq)%a(iDiag) = as(ieq)%a(iDiag)!+9*d*hx/(dx*(3*hx*dx+8*cnd)) !Left
      as(ieq)%a(ii) = as(ieq)%a(ii)!-d*hx/(dx*(3*hx*dx+8*cnd))    !Left
      ieqt = ieq + nx - 1
      iDiag = as(ieqt)%iDiag
      ii = iDiag - 1
      as(ieqt)%a(iDiag) = as(ieqt)%a(iDiag)!+9*d*hx/(dx*(3*hx*dx+8*cnd))  !Right
      as(ieqt)%a(ii) = as(ieqt)%a(ii)!-d*hx/(dx*(3*hx*dx+8*cnd))    !Right
    ENDDO
!
!          Load the Band Matrix for use by dgbfa and dgbsl
!
    abd = 0.0_sdk          !   Initialize to zero
                           !   Then set nonzero terms
    DO ieq=1,n
      DO i = 1, as(ieq)%nCoef
        k = ieq - as(ieq)%iCol(i) + nx + nx + 1
        abd(k,as(ieq)%iCol(i)) = as(ieq)%a(i)
      ENDDO
    ENDDO
!
!   Set Explicit Method coefficients
!
    DO ieq=1,n
        as(ieq)%ax = (1-f)*as(ieq)%a
    ENDDO

!
!   Set Implicit Method coefficients
!
    ai = 0.0_sdk
    DO ieq=1,n
      DO i = 1, as(ieq)%nCoef
        k = ieq - as(ieq)%iCol(i) + nx + nx + 1
        ai(k,as(ieq)%iCol(i)) = f*as(ieq)%a(i)
        IF(i.NE.as(ieq)%iDiag) CYCLE
        ai(k,as(ieq)%iCol(i)) = ai(k,as(ieq)%iCol(i)) + 1.0_sdk
      ENDDO
    ENDDO
!
    IF (fullMatrix) THEN
!
!      Load full matrix for comparison with student results
!
       af=0
       DO ieq=1,n
         DO i = 1, as(ieq)%nCoef
            af(ieq,as(ieq)%iCol(i)) = as(ieq)%a(i)
         ENDDO
       ENDDO
    ENDIF

    END SUBROUTINE SetCoeffs

    SUBROUTINE SetRHS

!   Set the right hand side of the equation.  Note that when dton is 1.0
!   you get the right hand side of a transient equation.  When dton is 0.0
!   you get the right hand side of a steady state equation.
!
!                     Programmed by John Mahaffy 1/01
    IMPLICIT NONE
    INTEGER(sik) :: j, i, k, kp
!    pm = qm*delt/(cpm*rhom)                         ! This is required in the test case
! 
!     Contributions common to all volumes 
!
    DO j=1,ny
      DO i=1,nx
        k = (j-1)*nx + i
        bm(k) = pm + dton*ts(i,j) !  (dton allows use of coding for transients
      ENDDO
    ENDDO
!
! Adjustments for the y face boundary conditions
!
    DO k=1,nx
       r_o = r_i + (0.5_sdk)*dy
       bm(k) = bm(k)+8*d*hy*txka(k)*(r_o-dy*0.5_sdk)/(r_o*dy*(3*hy*dy+8*cnd))
    !  kp = k + n - nx
    !  bm(kp) = bm(kp)+8*d*hy*tf/(dy*(3*hy*dy+8*cnd))
    ENDDO

! Adjustments for the x face boundary conditions

    DO k=1,n,nx
      bm(k) = bm(k)!+8*d*hx*tf/(dx*(3*hx*dx+8*cnd))
      kp = k + nx - 1
      bm(kp) = bm(kp)!+8*d*hx*tf/(dx*(3*hx*dx+8*cnd))
    ENDDO
!
    IF(dton.EQ.0.0_sdk) RETURN      !  Leave now for the steady state case
!                                      Otherwise add transient terms
    DO i = 1,n
      DO j = 1,as(i)%nCoef
        k = as(i)%iCol(j)
        bm(i) = bm(i) - as(i)%ax(j)*tt(k)
      ENDDO
    ENDDO

    RETURN
    END SUBROUTINE SetRHS

    END MODULE
