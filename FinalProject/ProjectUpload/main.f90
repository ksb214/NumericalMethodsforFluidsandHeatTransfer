    PROGRAM testFluid
!
!   Program to test numerical methods
!
    USE Input
    USE FluidArrays
    USE Trans
    USE Matrix
    USE Location
    IMPLICIT NONE
!
    CALL SetIC            !  Set initial conditions
    CALL SetLocAr         !  Set arrays to cross reference data structures
    CALL AllocMat         !  Allocate Matrices
    CALL InitFluid        !  Initialize dependent variables
    CALL Transient        !  Run a transient
!
    STOP
!
    END
