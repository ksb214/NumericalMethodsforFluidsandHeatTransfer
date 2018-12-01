      MODULE Data
      USE IntrType
      USE Eos
      USE ScalarDat

      REAL(sdk), ALLOCATABLE :: abd(:,:)  !  Steady State Matrix in band storage form
      REAL(sdk), ALLOCATABLE :: ai(:,:)   !  Implicit portion of the Transient matrix 
      REAL(sdk), ALLOCATABLE :: bm(:)      !  Right hand side (source and boundary terms)
      REAL(sdk), ALLOCATABLE :: af(:,:)   !  Full storage of the Steady State Matrix
!
      INTEGER(sik), ALLOCATABLE :: mipvt(:) !  Pivoting information for matrix solvers

      REAL(sdk), ALLOCATABLE :: r(:),re(:),ts(:,:),tt(:),txs(:,:),tys(:,:),p_r(:), err_r(:), exact_r(:)
      REAL(sdk), ALLOCATABLE :: tx(:), ty(:), resid(:), tss(:,:), txka(:)
!
!     Data Structure to start dealing with Matrix Sparcity
!
      TYPE SparseMatrixT
        INTEGER(sik) :: nCoef                   !  number of nonzero coefficients in the row
        REAL(sdk), POINTER :: a(:)=>NULL()      !  nonzero coefficients
        REAL(sdk), POINTER :: ai(:)=>NULL()     !  implicit nonzero coefficients
        REAL(sdk), POINTER :: ax(:)=>NULL()     !  explicit nonzero coefficients
        INTEGER(sik),POINTER :: iCol(:)=>NULL() !  column index for the nonzero coefficient
        INTEGER(sik) :: iDiag                   !  index of the diagonal coefficient
        INTEGER(sik), POINTER:: iDim(:)=>NULL() !   index of the 3-D dimension associated
                                                !  with the coefficient (e.g 1=x, 2=y, 3=z)
                                                !  (used only for ADI)
      END TYPE SparseMatrixT
!
!     The following derived type array will contain the coefficients of the steady 
!     state matrix.
!
      TYPE(SparseMatrixT), ALLOCATABLE :: as(:)
!
!     Data For running with multiple time steps and analyzing relative errors
!     Set numStepsTried = 0 to only do steady state calculations.
!
      INTEGER(sik), PARAMETER :: numStepsTried = 1
!
!      REAL(sdk) :: timeSteps(3) =(/0.01_sdk, 0.02_sdk, 0.04_sdk/)           ! For implicit method
!                   0.05_sdk, 0.1_sdk, 0.2_sdk, 0.4_sdk/)
      REAL(sdk) :: timeSteps(3) =(/0.1_sdk, 0.2_sdk, 0.4_sdk/)        

!      REAL(sdk) :: timeSteps(3) =(/0.001_sdk, 0.002_sdk, 0.004_sdk/) !, & ! For explicit method
!                   0.0005_sdk, 0.0008_sdk, 0.001_sdk, 0.002_sdk/)
!
!     Storage for temperatures at the last time step for each stepsize in
!     timeSteps
!
      REAL(sdk), ALLOCATABLE :: tLast(:,:)
!
!     Data for testing mesh refinement and analyzing relative errors
!
      INTEGER(sik), PARAMETER :: numMeshTried = 1 
      INTEGER(sik) :: imesh
!
!     Number of cells in both directions (needs to be multiples of 3
!     for easy analysis.
!
      INTEGER(sik) :: nCells(5) =(/3 , 3 , 27, 81, 243/)
!
!     Storage for temperatures used in Richardson Error analysis
!
      REAL(sdk) :: tCenter(numMeshTried), txwc(numMeshTried), tywc(numMeshTried)
!
!     Various Scalar variables
!

!
!     Set the following flag to true to compare results with student homeworks.
!
      LOGICAL ::  fullMatrix = .TRUE.
!
      REAL(sdk) :: ysize , cnd, tf, dxm, dy, qm, t11ss, r_i, r_o, delt, Tik, r_ot
      REAL(sdk) :: cpm, rhom, pm, d, dton=0.0_sdk, delt_exp
      REAL(sdk) :: hx
      INTEGER(sik):: nx, ny, n
      REAL(sdk) :: f                   ! method selector
      REAL(sdk) :: rTime               ! Time step refinement ratio
      REAL(sdk) :: rSpace              ! Mesh refinement ratio
      Real(sdk):: Tana                 ! Analytically calculated temperature for the case of zero heat transfer coefficientss 
      INTEGER(sik) :: iLUlevel, resetGMRES, methodLB, methodUB
      LOGICAL :: doDirect = .TRUE.
      CONTAINS

      SUBROUTINE SetData

      IMPLICIT NONE
      REAL(sdk) :: rLast
      INTEGER(sik) :: i
!
      r_i = 0.002_sdk
      ysize = 0.001_sdk      !  height of the bar
      r_ot = r_i + ysize
      cnd = 80.0_sdk           !  conductivity
      tf = 300.0_sdk           !  fluid temperature
      cpm = 500.0_sdk         !  specific heat (not needed for this problem)
      rhom = 8000.0_sdk       !  density       (not needed for this problem)
      qm = 6366197.0_sdk        !  power source per unit volume

      nx = ncell ! nCells(imesh)   !  number of volumes in the x direction
      ny = 30    ! nCells(imesh)   !  number of volumes in the y direction
  
      dxm = xsize/nx
      dy = ysize/ny
      n=nx*ny                !  Total number of volumes

      f = 1.0_sdk            ! numerical method selector
                             ! = 1     Fully Implicit
                             ! = 0     Fully Explicit
                             ! = 0.5   Crank-Nicholson
      delt = dt      
      hx =  0.0_sdk        !  heat transfer coefficient on faces at fixed x
   
      delt_exp = 0.5_sdk*dxm**2*dy**2*rhom*cpm/(cnd*(dxm**2+dy**2))   !  time step limitation for explicit method

    ! Limiting the time step of the explicit method to the limiting time step when time step set is greater than the stability limit
  
      IF ((f.EQ.0.0_sdk).AND.(delt.GT.delt_exp))THEN
         delt = delt_exp
      END IF

      IF (numStepsTried .EQ. 0) THEN
         delt = 1.0_sdk
         rhom =  1.0_sdk
         cpm  =  1.0_sdk
      ELSE IF (numStepsTried .EQ. 1) THEN
         timeSteps(1) = delt
         rTime = 1.0_sdk
      ELSE

!     Set and check Time Refinement Ratio

         rTime = timeSteps(2)/timeSteps(1)
         DO i = 3, numStepsTried
            rLast = rTime
            rTime = timeSteps(i)/timeSteps(i-1)
            IF (ABS(rTime-rLast)/ABS(rTime) .LT. 1.e-6_sdk) CYCLE
            WRITE (*,*) 'Selected Time Steps do not have a constant refinement ratio.'
            STOP
         END DO

!     Set and check Spatial Mesh Refinement Ratio

      END IF
      IF (numMeshTried .EQ. 1) THEN
         rSpace = 1.0_sdk
      ELSE

         rSpace = REAL(nCells(2),sdk)/REAL(nCells(1),sdk)
         DO i = 3, numStepsTried
            rLast = rSpace
            rSpace = REAL(nCells(i),sdk)/REAL(nCells(i-1),sdk)
            IF (ABS(rSpace-rLast)/ABS(rSpace) .LT. 1.e-6_sdk) CYCLE
            WRITE (*,*) 'Selected Meshes do not have a constant refinement ratio.'
            STOP
         END DO
      END IF

      iLUlevel = 2           !  level of incomplete factorization
      resetGMRES = 10        !  Number iterations before GMRES is restarted
      doDirect = .TRUE.      !  Include a direct matrix solution?


      IF ( n .GT. 100) THEN
         fullMatrix = .FALSE.
      END IF
!
      END SUBROUTINE SetData

      SUBROUTINE AllocArrays
!
!     Allocate almost all arrays needed for this problem, and initialize
!     two of them.  Note that components of the derived type array "as"
!     must be allocated later, when more information has been generated
!     about the spatial mesh.
!     
!     Programmed by John Mahaffy 1/6/05
!
      IMPLICIT NONE
      INTEGER(sik) :: n2
!
      ALLOCATE (ai(3*nx+1,n), mipvt(n),abd(3*nx+1,n))
      ALLOCATE ( bm(n), tss(nx,ny))
      ALLOCATE (ts(nx,ny),tt(n),txs(ny,2), tys(nx,2), tx(nx), ty(ny), p_r(n), err_r(n), exact_r(n),txka(nx) )

      IF (fullMatrix) THEN
         ALLOCATE ( af(n,n) )
      ENDIF
!
!     Last set of temperatures for each transient
!
      n2 = MAX(1,numStepsTried)
      ALLOCATE (tLast(n,n2))
!
       bm = 0.0                !  Clear the right hand side array
       ts = tf                  !  Initial temperature condition
!
!     Sparse Matrix data structure
!
      ALLOCATE (as(n))

!
      END SUBROUTINE AllocArrays
!
      SUBROUTINE DeallocArrays
!
!     Deallocate all arrays use for this problem.
!     
!     Programmed by John Mahaffy 1/6/05
!
      IMPLICIT NONE
      INTEGER(sik) :: i
      DEALLOCATE (ai, mipvt,abd)
      DEALLOCATE ( bm, tss)
      DEALLOCATE (ts,tt,txs, tys, tx, ty,txka)

      IF (fullMatrix) THEN
         DEALLOCATE ( af )
      ENDIF
!
!     Deallocate components of the sparse matrix derived type first to
!     prevent a memory leak
!
      DO i = 1,n
         DEALLOCATE(as(i)%a, as(i)%ax, as(i)%iCol, as(i)%iDim)
      ENDDO
!
!     Now the sparse matrix derived type array can be safely deallocated
!
      DEALLOCATE (as)
!
      END SUBROUTINE DeallocArrays

      SUBROUTINE SaveTemperatures
!
!     Saves temperatures at the mesh center and halfway along 
!     one x and one y boundary.
!     This only works if nx and ny are odd integers
!
!     Programmed by John Mahaffy  1/6/05
!
      IMPLICIT NONE
!
!      tCenter(imesh) = tLast(n/2+1,1)
!
!     Surface temperatures from Simple First Order Surface Flux Match
!
!      txwc(imesh) = (hx*tf+2*cnd/dx*tLast(1+ (ny/2)*nx,1))/(hx+2*cnd/dx)
!      tywc(imesh) = (hy*tf+2*cnd/dy*tLast(nx/2+1,1))/(hy+2*cnd/dy)
!
!     Surface temperatures from Second Order Surface Flux Match
!

      txwc(imesh) = (3*hx*tf+9*cnd/dxm*tLast(1+(ny/2)*nx,1)-cnd/dxm*tLast(2+(ny/2)*nx,1))/(3*hx+8*cnd/dxm)
      tywc(imesh) = (3*hy*tf+9*cnd/dy*tLast(nx/2+1,1)-cnd/dy*tLast(nx/2+1+nx,1))/(3*hy+8*cnd/dy)
!
      END SUBROUTINE SaveTemperatures
!
      END MODULE Data


