    PROGRAM steady
!                 Steady State 2-D Cartesian conduction
!
!                 Programmed by John Mahaffy   1/01
!
    USE Data
    USE LUsolve
    USE SetEqn

    IMPLICIT NONE

    INTEGER(sik) :: i,info, nstep
!                             
    CALL SetData    ! Initialize data describing the problem
!                    
    CALL SetCoeffs  ! Load the coefficients for the linear system
    CALL SetRHS     ! Load the right hand side of the equation system

!                                Scale terms to match a student's formulation
!    a = -a/cnd*dx*dy
!    b = -b/cnd*dx*dy
    aa =  a                      !  Save a copy of the original equations for
    tt =  b                      !  later comparisons during debugging
!
!   Factor the Matrix
    CALL dgefa(aa,n,n,ipvt,info)
!   Solve the factored system
    CALL dgesl(aa,n,n,ipvt,tt,0)
!                                  Load a 2-D temperature array from the 
!                                  1-D solution vector
    t = RESHAPE(tt,(/nx,ny/))
    CALL edit
    STOP                        !  Output of key variables

    END

