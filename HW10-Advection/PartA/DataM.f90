    MODULE Data

    USE IntrType
    REAL(sdk) :: v, dx , c, dt, time, length
    INTEGER(sik) :: n, np, iBC, lda
    INTEGER(sik), ALLOCATABLE :: ipvt(:)
    REAL(sdk), ALLOCATABLE :: rho(:), rho1(:), rho2(:), rhon(:), rhom(:)
    REAL(sdk), ALLOCATABLE :: rhoDonorExp(:), rhoDonorImp(:), rhoLeith(:), rhoQuick(:)
    REAL(sdk), ALLOCATABLE :: rhoQuickExp(:), rhoQuickImp(:)
    REAL(sdk), ALLOCATABLE :: rhoQuickest(:), aImp(:,:)

    CHARACTER*1 :: tab=ACHAR(9)


    CONTAINS

    SUBROUTINE Setup
!
!   Set initial and boundary conditions and initialize array space
!
    IMPLICIT NONE

    length = 5.0_sdk
    n = 100           !  number of volumes evaluating the advection Equation
    v = 1.0_sdk      !  Velocity
    dx = length/n    !  Mesh spacing
    c = 0.5_sdk
              !  Courant Number
    dt = c*dx/v      !  Time Step

!    dt = 0.001_sdk
!   c = dt*v/dx

    np= n + 1        !  number of volume edges
!
!     Exit Boundary condition flag
!     iBC = 0    Constant value
!           1    zero derivative
!           2    zero curvature
!
    iBC = 2

!
!
!   Allocate storage to permit two boundary points at and below x=0 and one
!   boundary point above x=5.0
!
    ALLOCATE(rho(-1:np),rhon(-1:np),rho1(0:np), rho2(0:np), rhom(1:np))
    ALLOCATE(rhoLeith(0:n), rhoQuickExp(0:n), rhoQuickImp(0:n))
    ALLOCATE(rhoDonorExp(0:n), rhoDonorImp(0:n), rhoQuickest(0:n))
!
!   ALLOCATE aImp with enough space to handle an implicit QUICK method
!
    lda = 6
    ALLOCATE (aImp(lda,n))
    ALLOCATE (ipvt(n))
!
    time = 0.0             !  Initial Time

    END SUBROUTINE Setup

    END MODULE Data
