    MODULE Matrix
    USE IntrType
!
!   Array storage for the linear system used in the Newton iteration
!
    REAL(sdk), ALLOCATABLE :: a(:,:), x(:), b(:)
    INTEGER(sik), ALLOCATABLE :: ipvt(:)

    CONTAINS

    SUBROUTINE AllocMat
!
!   Allocate space for the linear system
!
    USE ScalarDat

    ALLOCATE ( a(nvar,nvar), x(nvar), b(nvar), ipvt(nvar) )

    END SUBROUTINE AllocMat

    SUBROUTINE Solve
!                     Apply a linear solver appropriate to the matrix
    USE ScalarDat
    USE LUsolve
!
    IMPLICIT NONE
    INTEGER(sik) :: info

    CALL sgefat(a,nvar,nvar,ipvt,info)
    CALL sgeslt(a,nvar,nvar,ipvt,b,0)

    END SUBROUTINE Solve

    END MODULE Matrix
