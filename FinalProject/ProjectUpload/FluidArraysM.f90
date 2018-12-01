    MODULE FluidArrays
    USE IntrType
    USE ScalarDat
!
!   Arrays containing state information about the fluid in 1-D flow
!   All use SI units
!
!   Arrays names ending with "N" are evaluated at the new (n+1) time  level
!   others are values at the start of the time step.
!
    REAL(sdk), ALLOCATABLE :: p(:)     ! Start of step pressure
    REAL(sdk), ALLOCATABLE :: pN(:)    ! End of step pressure
    REAL(sdk), ALLOCATABLE :: T(:)     ! Start of step temperature
    REAL(sdk), ALLOCATABLE :: TN(:)    ! End of step temperature
    REAL(sdk), ALLOCATABLE :: e(:)     ! Start of step specific internal energy
    REAL(sdk), ALLOCATABLE :: eN(:)    ! End of step specific internal energy
    REAL(sdk), ALLOCATABLE :: rho(:)   ! Start of step density
    REAL(sdk), ALLOCATABLE :: rhoN(:)  ! End of step density
    REAL(sdk), ALLOCATABLE :: rhoe(:)  ! Start of step internal energy per vol
    REAL(sdk), ALLOCATABLE :: rhoeN(:) ! End of step internal energy per vol
    REAL(sdk), ALLOCATABLE :: v(:)     ! Start of step cell edge velocity
    REAL(sdk), ALLOCATABLE :: vN(:)    ! End of step velocity
    REAL(sdk), ALLOCATABLE :: q(:)     ! power source per volume
    REAL(sdk), ALLOCATABLE :: fric(:)  ! wall friction coefficient per volume
    REAL(sdk), ALLOCATABLE :: Tm(:)    ! Metal temperature per unit volume
!
!   Geometry arrays
!
    REAL(sdk), ALLOCATABLE :: dx(:)    ! length of each computational volume
    REAL(sdk), ALLOCATABLE :: vol(:)   ! fluid volume of each finite volume
    REAL(sdk), ALLOCATABLE :: fa(:)    ! fluid flow area at each volume edge
    REAL(sdk), ALLOCATABLE :: hd(:)    ! hydraulic diameter at each volume edge
!
!   Pressure, temperature and velocity are independent variables.  The
!   TARGET attribute makes it easier to add perturbations from the solution
!   arrays at the end of each iteration.
!
    TARGET :: pN, TN, vN

    CONTAINS

    SUBROUTINE StartStep
    IMPLICIT NONE
!
!   Copy end-of-step variables from the last step into start-of-step
!   variables for this time step
!
    rho = rhoN
    rhoe = rhoeN
    p = pN
    v = vN
    e = eN
    T = TN

    END SUBROUTINE StartStep

    SUBROUTINE InitFluid
    USE Eos
    IMPLICIT NONE
    INTEGER(sik) :: i
!
!   Initialize all dependent fluid state variables
!
    DO i = LBOUND(p,DIM=1), UBOUND(p,DIM=1)
      rhoN(i) = rhoLiq (TN(i), pN(i))
      eN(i) = SpEnergy(TN(i), pN(i))
      rhoeN(i) = rhoN(i)*eN(i)
    ENDDO
!
    END SUBROUTINE InitFluid

    SUBROUTINE InitFluid2
    USE Eos
    IMPLICIT NONE
    INTEGER(sik) :: i
!
!   Initialize all dependent fluid state variables
!
    DO i = -nbnd+1, 0
      rhoN(i) = rhoLiq (TN(i), pN(i))
      eN(i) = SpEnergy(TN(i), pN(i))
      rhoeN(i) = rhoN(i)*eN(i)
    ENDDO
!
    END SUBROUTINE InitFluid2

    SUBROUTINE AllocFluidAr
    USE ScalarDat

!    Allocate state variable arrays, space is reserved at each end for 
!    storage of boundary information.

    ALLOCATE(pN(1-nbnd:ncell+nbnd))
    ALLOCATE(p(1-nbnd:ncell+nbnd))
    ALLOCATE(TN(1-nbnd:ncell+nbnd))
    ALLOCATE(T(1-nbnd:ncell+nbnd))
    ALLOCATE(rho(1-nbnd:ncell+nbnd))
    ALLOCATE(rhoN(1-nbnd:ncell+nbnd))
    ALLOCATE(e(1-nbnd:ncell+nbnd))
    ALLOCATE(eN(1-nbnd:ncell+nbnd))
    ALLOCATE(rhoe(1-nbnd:ncell+nbnd))
    ALLOCATE(rhoeN(1-nbnd:ncell+nbnd))
    ALLOCATE(v(1-nbnd:ncell+nbnd+1))
    ALLOCATE(vN(1-nbnd:ncell+nbnd+1))
    ALLOCATE(q(1:ncell))
    ALLOCATE(Tm(1:ncell))
    ALLOCATE(fric(1:ncell+1))
    ALLOCATE(dx(1-nbnd:ncell+nbnd))
    ALLOCATE(fa(1-nbnd:ncell+nbnd+1))
    ALLOCATE(hd(1-nbnd:ncell+nbnd+1))
    ALLOCATE(vol(1-nbnd:ncell+nbnd))
    

    END SUBROUTINE AllocFluidAr

    END MODULE FluidArrays
