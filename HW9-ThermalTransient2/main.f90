    PROGRAM Richardson
!                 Steady  2-D Cartesian conduction
!                 with a spatial error analysis
!                 Programmed by John Mahaffy   2/05
!
    USE Data
    USE LinAlg
    USE SetEqn
    USE Output, ONLY : SpatialRichardson, edit

    IMPLICIT NONE

    INTEGER(sik) :: info, nstep
    REAL(sdk) :: tstart, tend
!
!   Open a file to store error information
!
    OPEN (22, FILE='errors.txt')
!
!   Loop over all meshes being tested
!
    DO imesh = 1, numMeshTried
!                             
       CALL SetData     ! Initialize data describing the problem
       CALL AllocArrays ! Allocate arrays for state and equation information                  
       CALL SetCoeffs   ! Load the coefficients for the linear system
       CALL SetRHS      ! Load the right hand side of the equation system
!
       tt = b                
!
       CALL CPU_TIME(tstart)

!      Factor the Matrix
       CALL dgbfa(abd,3*nx+1,n,nx,nx,ipvt,info)
!      Solve the factored system
       CALL dgbsl(abd,3*nx+1,n,nx,nx,ipvt,tt,0)

       IF (fullMatrix) THEN
!          af = -af/cnd
 !         b  = -b/cnd
       END IF

       CALL CPU_TIME(tend)
       WRITE(*,*)   'Direct Solver CPU = ',tend-tstart,' seconds for ', n, 'volumes'
!
!                                  Load a 2-D temperature array from the 
!                                  1-D solution vector
       tss = RESHAPE(tt,(/nx,ny/))
!
       IF (numStepsTried .NE. 0) THEN
          CALL Transient      !  Run the transient
                              !  Look at trans.txt for output
       ELSE
          tLast(:,1) = tt     !  Load tLast if the transient is not run
          t = tss
          CALL edit
       END IF
!
       CALL SaveTemperatures  !  Save some temperatures for Richardson analysis    
!
       CALL DeallocArrays
    END DO

    CALL SpatialRichardson
!
    STOP

    END

    SUBROUTINE Transient
!                 Transient  2-D Cartesian conduction
!                 Explicit (Forward) Time Difference
!                 Programmed by John Mahaffy   2/01
!
    USE Data
    USE SetEqn
    USE Output
    USE LinAlg

    IMPLICIT NONE

    INTEGER(sik) :: info, istep, nstep

!
    dton = 1.0_sdk       !    Activate Transient terms in equations
!
!   The outer loop allows several runs with different time steps for
!   sensitivity studies
!
    DO istep = 1,numStepsTried
       delt = timeSteps(istep)
       CALL SetCoeffs          ! Load the coefficients for the linear system
!
!   Factor the Implicit Matrix
!   This is done first because no elements of the implicit matrix
!   change with time
!
       CALL dgbfa(ai,3*nx+1,n,nx,nx,ipvt,info)

       time = 0.0_sdk
!
       tt = tf                     !  Set Initial Conditions
       t  = tf
!
       DO nstep = 1,100000   !    Transient Loop

         CALL SetRHS
!
         tt = b
!        Solve the factored system
         CALL dgbsl(ai,3*nx+1,n,nx,nx,ipvt,tt,0)
         t = RESHAPE(tt,(/nx,ny/))
         time = time + delt
         CALL edit
         IF(time.GE.endTime-0.00001*delt) EXIT
       ENDDO
       
       tLast(:,istep) = tt(:)

    END DO
!
!   Process results and output several error estimates for all time steps used.
!
    CALL TimeErrorEstimate
!
    END SUBROUTINE
