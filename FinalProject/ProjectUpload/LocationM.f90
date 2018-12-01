    MODULE Location
!
!   This contains several types of locaters.  The first gives locations and
!   associated weighting used to perform averages or differentials of
!   quantities.  The second provides pointers to move information from
!   the data structure associated with the full system of equations and
!   unknowns, to the data structure associated with physical state variables
!   and spatial locations.  The third is a set of translation indices to
!   provide the full system variable index for each independent state
!   variable. The fourth lists independent variables associated with

    USE Intrtype

    IMPLICIT NONE
!=======================================================================
!                            Information for averages (and differences)
    TYPE averageT
      INTEGER(sik) :: n                !   number of points used 
      INTEGER(sik), POINTER :: i(:)    !   cell or face index
      REAL(sdk), POINTER :: wf(:)      !   weight factor
    END TYPE
!
    TYPE (averageT), ALLOCATABLE, TARGET:: fluxAv(:), dPdx(:), dV(:)
!    fluxAv  -  information used to construct averages in cell edge
!               mass and energy fluxes
!    dPdx    -  information used to evaluate the pressure gradient
!               in the momentum equation
!    dV      -  information used to evaluate the velocity derivative
!               in momentum transport term
!
!========================================================================
!                           Translation from system to physical variables
    TYPE varPtrT
      REAL(sdk), POINTER :: var
      INTEGER(sik) :: i, itype
    END TYPE
!       %var    -   pointer to the variable storage in the fluid arrays
!       %i      -   index in the fluid array
!       %itype  -   variable type
!                   1   -  pressure
!                   2   -  temperature
!                   3   -  velocity
!
    TYPE (varPtrT), ALLOCATABLE :: loc(:)
!
!=======================================================================
!                Translation for physical to system variables
!
    INTEGER(sik), ALLOCATABLE ::  ivarP(:), ivarT(:), ivarV(:)
!
!        ivarP  -  index in the linear solution array "b" for the change
!                  in the corresponding element of pN
!        ivarT  -  index in the linear solution array "b" for the change
!                  in the corresponding element of TN
!        ivarV  -  index in the linear solution array "b" for the change
!                  in the corresponding element of vN
!
!========================================================================
!                 List of active varibles in each equation
!                 Also can be considered as a list of column numbers
!                 in each matrix row with potential for non-zero elements
    TYPE varListT
      INTEGER(sik) :: n, iloc
      INTEGER(sik), POINTER :: i(:)
      INTEGER(sik) :: eqnType
    END TYPE
!              n   -   number of columns with potential for non-zero elements
!              iloc  - index of cell or face location for the equation
!              i   -   list of column indices
!              eqnType -  Equation type associated with this row
!                         1  -  mass
!                         2  -  energy
!                         3  -  momentum
!
    TYPE (varListT), ALLOCATABLE :: varList(:)

    CONTAINS

    SUBROUTINE SetLocAr
!
!   Allocate and load indices in the locator arrays
!
!   Note that I only record potential nonzero contributions to matricies
!   Depending on direction of velocity, some coefficients marked as potentially
!   non-zero in structures below may actually be zero.
!
    USE ScalarDat
    USE FluidArrays

    IMPLICIT NONE

    INTEGER(sik) :: i, icol, irow, ivar, iLB, iUB, j, jj, nv
!
!                     Lower and Upper Bounds of velocity (cell edge) arrays
    ivLB = 1
    ivUB = ncell+1
                              !  Calculate the total number of unknowns
    nvar = 3*ncell + 1        !  p & T in each volume, velocity at each edge
    IF (lBC.EQ.'v') THEN
      ivLB = 2                !  velocity fixed at left edge, do not evaluate
      nvar = nvar - 1         !  momentum equation at face 1
    ENDIF

    IF (uBC.EQ.'v') THEN
      ivUB = ncell            !  velocity fixed at right edge, do not evaluate
      nvar = nvar -1          !  momentum equation at face ncell+1
    ENDIF

    ALLOCATE(loc(nvar),varList(nvar))
    ALLOCATE(fluxAv(ncell+1), dPdx(ivLB:ivUB), dV(ivLB:ivUB))
!
!        Load tables used for weight factors
!
    DO i = 1,ncell+1
!                 Number of volumes that might contribute to an edge average
      fluxAv(i)%n = 2*nbnd
      ALLOCATE (fluxAv(i)%i(fluxAv(i)%n))
      ALLOCATE (fluxAv(i)%wf(fluxAv(i)%n))
!           Load indices of volumes that might contribute to an edge average
      DO j = -nbnd, nbnd-1
        fluxAv(i)%i(j+nbnd+1) = i + j
      ENDDO
    ENDDO

    DO i = ivLB,ivUB
!                  Number of volumes contributing to the pressure derivative
      dPdx(i)%n = 2           
      ALLOCATE (dPdx(i)%i(dPdx(i)%n))
      ALLOCATE (dPdx(i)%wf(dPdx(i)%n))
      dPdx(i)%i(1) = i - 1               ! Indicies of volumes contributing
      dPdx(i)%i(2) = i                   ! to the derivative

!                   Weighting factors to generate the derivative.
!                   Note that calculation here requires a fixed mesh.
      dPdx(i)%wf(2) = 2.0_sdk/(dx(i)+dx(i-1))
      dPdx(i)%wf(1) = - dPdx(i)%wf(2)
    ENDDO

    DO i = ivLB, ivUB
!                  Number of edges contributing to the velocity derivative
      dV(i)%n = 1+ 2*nbnd
      ALLOCATE (dV(i)%i(dV(i)%n))
      ALLOCATE (dV(i)%wf(dV(i)%n))
!           Load indices of edges that might contribute 
      DO j = -nbnd, nbnd
        dV(i)%i(j+nbnd+1) = i + j
      ENDDO
    ENDDO
!                             !  Translation of state to system variables
    iLB = LBOUND(p,DIM=1)
    iUB = UBOUND(p,DIM=1)

    ALLOCATE ( ivarP(iLB:iUB), ivarT(iLB:iUB))

    iLB = LBOUND(v,DIM=1)
    iUB = UBOUND(v,DIM=1)

    ALLOCATE ( ivarV(iLB:iUB))

    ivarP = 0
    ivarT = 0
    ivarV = 0

    ivar = 0              !  variable index in the full system of equations

    IF (ivLB.EQ.1) THEN  
      ivar = ivar+1
      loc(ivar)%var => vN(1)
      loc(ivar)%i = 1
      loc(ivar)%itype = 3
      ivarV(1) = ivar
    ENDIF

    DO i  = 1,ncell
      ivar = ivar+1               !  Pressure
      loc(ivar)%var => pN(i)
      loc(ivar)%i = i
      loc(ivar)%itype = 1
      ivarP(i) = ivar

      ivar = ivar+1               !  Temperature
      loc(ivar)%var => TN(i)
      loc(ivar)%i = i
      loc(ivar)%itype = 2
      ivarT(i) = ivar
                                  !  No system variable at right edge
      IF(i+1.GT.ivUB) EXIT        !  if it has a constant velocity BC

      ivar = ivar+1               !  Velocity
      loc(ivar)%var => vN(i+1)
      loc(ivar)%i = i+1
      loc(ivar)%itype = 3
      ivarV(i+1) = ivar

    ENDDO
!                     Calculate and store indices of all variables used
!                     in each equation.  This is a list of columns in
!                     each row of the matrix with potential for non-zero
!                     elements.
!
    irow = 0
    IF(ivarV(1).NE.0) THEN
      irow = irow + 1
      varList(irow)%n = 1      ! Velocity at this edge appears in the equation
      varList(irow)%iloc = 1
      varList(irow)%eqnType = 1
!                              ! Velocity at edge 2 may appear in the equation
      IF( ivarV(i+1).NE.0) varList(irow)%n = varList(irow)%n + 1
      varList(irow)%n = varList(irow)%n + 2*nbnd ! Variables in Dp & 1/rho terms

      ALLOCATE(varList(irow)%i(varList(irow)%n))
      j = 1

      varList(irow)%i(j) = ivarV(1)            !   Local Velocity
      varList(irow)%i(j+1) = ivarT(1)          !   Temp in 1/rho
      varList(irow)%i(j+2) = ivarP(1)          !   press in 1/rho and Grad P
      IF( ivarV(2).NE.0) varList(irow)%i(j+3) = ivarV(2)
    ENDIF
!
    DO i = 1,ncell
      irow = irow + 1                     !  Mass Equation
      varList(irow)%iloc = i
      varList(irow)%eqnType = 1
      varList(irow)%n = 2*( 2*nbnd+1)     !  Temperatures and Pressures
      IF( ivarV(i).NE.0) varList(irow)%n = varList(irow)%n + 1
      IF( ivarV(i+1).NE.0) varList(irow)%n = varList(irow)%n + 1
      ALLOCATE(varList(irow)%i(varList(irow)%n))
      nv = 0
      IF(ivarV(i).NE.0) THEN               !  Velocity for left edge flux
        nv = nv+1
        varList(irow)%i(nv) = ivarV(i)
      ENDIF
!                         Pressure and Temperatures in time derivative and
!                         Flux terms
      DO j = -nbnd,nbnd
        nv = nv+1
        varList(irow)%i(nv) = ivarP(i+j)
        nv = nv+1
        varList(irow)%i(nv) = ivarT(i+j)
      ENDDO

      IF(ivarV(i+1).NE.0) THEN             !  Velocity for right edge
        nv = nv+1
        varList(irow)%i(nv) = ivarV(i+1)
      ENDIF
      irow = irow + 1                 !  Energy Equation
                                      !  Same footprint as Mass Equation
      ALLOCATE(varList(irow)%i(varList(irow-1)%n))
      varList(irow) = varList(irow-1)
      varList(irow)%eqnType = 2
!
!                                Momentum Equation at the Right cell face
      IF(ivarV(i+1).EQ.0) CYCLE
      irow = irow + 1
      varList(irow)%iloc = i+1
      varList(irow)%eqnType = 3
      varList(irow)%n = 1
      IF( ivarV(i+2).NE.0) varList(irow)%n = varList(irow)%n + 1
      IF( ivarV(i).NE.0) varList(irow)%n = varList(irow)%n + 1
      varList(irow)%n = varList(irow)%n +  2*nbnd  ! Variables in 1/rho term

      IF(ivarP(i+1).NE.0) THEN            ! Dp and 1/rho contributions to right
         varList(irow)%n = varList(irow)%n + 2*nbnd
      ENDIF
      ALLOCATE(varList(irow)%i(varList(irow)%n))
      varList(irow)%i = 0
      nv = 0
      IF( ivarV(i).NE.0) THEN
        nv = nv +1
        varList(irow)%i(nv) =  ivarV(i)
      ENDIF
      nv =  nv + 1
      varList(irow)%i(nv) = ivarV(i+1)        !   Local Velocity

      DO j = -nbnd+1,nbnd                    !  1/rho and Grad P terms
        IF(ivarP(i+j).EQ.0) CYCLE
        nv = nv + 1
        varList(irow)%i(nv) = ivarT(i+j)
        nv = nv + 1
        varList(irow)%i(nv) = ivarP(i+j)
      ENDDO
      IF( ivarV(i+2).NE.0) varList(irow)%i(nv+1) = ivarV(i+2)
    ENDDO


    END SUBROUTINE SetLocAr

    END MODULE Location
