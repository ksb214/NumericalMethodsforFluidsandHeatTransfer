    MODULE Trans
    USE IntrType
    USE Eos
!    USE TransCond
    !
!   Subprograms used in the flow transient
!
    CONTAINS

    SUBROUTINE Transient
!
!   Driver for the transient
!
    USE ScalarDat
    USE FluidArrays
    USE Matrix
    USE Outputf

    USE Data
    USE LinAlg
    USE SetEqn
    USE Moutput, ONLY : medit
    
    IMPLICIT NONE
    INTEGER(sik) :: info    
    CALL SetData     ! Initialize conduction data describing the problem
    CALL AllocArrays ! Allocate arrays for state and equation information                  
    CALL SetCoeffs   ! Load the coefficients for the linear system

!    CALL SetRHS      ! Load the right hand side of the equation system
!    tt = bm                
!    Factor the Matrix
       CALL dgbfa(abd,3*nx+1,n,nx,nx,mipvt,info)

       tt = tf            ! Initialise the solution 
       ts = tf          

    time = 0.0_sdk
    DO nstep = 1,nsmax 
       
      CALL StartStep
      converged = .FALSE.
      CALL SetWf
      CALL Residuals(0)
      DO it = 1,itmax
 
        CALL SetRHS
         tt=bm
         CALL dgbsl(abd,3*nx+1,n,nx,nx,mipvt,tt,0)
         ts = RESHAPE(tt,(/nx,ny/)) 
       !  IF(time.GT.100.0_sdk)then 
       !     call medit                 !Activate this for recording the metal temperature at every iteration during the test case or every time step and in the test case
       !  end if
            CALL SetWf
         CALL Jacobian
         CALL Solve
         CALL EvalVar
         CALL Residuals(1)
         IF(converged) EXIT
      ENDDO
      IF(.NOT.converged) THEN
        WRITE (*, '(a,f7.2)') 'Iteration failed at time ',time
        STOP
      ENDIF
      time = time + dt
      
      ! Block for the test case
      !    IF (time.GT.100.0_sdk)THEN     ! Test case 
      !       qm=0.0_sdk
      !    END IF

      ! Applying zero curvature boundary condition for velocity and other properties are set same as the outlet to the ghost cells
     
      DO z = ncell+1,ncell+nbnd 
          
         rhoN(z) = rhoN(ncell)              
         rhoeN(z)= rhoeN(ncell)             
                                            
         vN(z) = 2*vN(z-1)-vN(z-2)          ! Zero curvature boundary to velocity only
         eN(z) = eN(ncell)                  
         TN(z) = TN(ncell)                  
      ENDDO
      vN(ncell+nbnd+1) = 2*vN(ncell+nbnd)-vN(ncell+nbnd-1)      
      
         
     ! Updating the transient input condition
         
      IF (time .GT. 2.0000001_sdk) THEN
            TN(nbnd-2:nbnd-1) = 400.0_sdk   
         ELSE
            TN(0) = 300.0_sdk + 75.0_sdk*time**2-25.0_sdk*time**3                             
            TN(-1)= 300.0_sdk + 75.0_sdk*time**2-25.0_sdk*time**3                             
         ENDIF

         !  Calculating metal surface temperature for heat transfer to fluid in fluid solver
         DO zi = 1, ncell
            Tm(zi) = (3.0_sdk*dy/(3.0_sdk*hy*dy+8.0_sdk*cnd))*(hy*TN(zi)+ (3.0_sdk*cnd)/dy*ts(zi,1)-(cnd/(3.0_sdk*dy))*ts(zi,2))
         end do
                  
         DO zi = 1, ncell
            q(zi)= hy*( Tm(zi) - TN(zi))*pi*dx(zi)*df/vol(zi)  
            txka(zi) = TN(zi)                                 ! Copying fluid temperature for calculating heat transfer from metal to fluid in metal solver
         end do
           
         call InitFluid2
         
         IF( time.GT. endtime-0.01*dt)then
            CALL EDITF         ! Writing the output file of fluid
            call medit         ! Writting output file of the metal after transient iteration complets
            EXIT
         end if
         
      ENDDO


    END SUBROUTINE Transient

    SUBROUTINE EvalVar
    USE FluidArrays
    USE Location
    USE Matrix
    USE ScalarDat
    USE Eos

    IMPLICIT NONE

    INTEGER(sik) :: i
!
!   Update the state variable arrays
!
    DO i = 1,nvar
      loc(i)%var = loc(i)%var + b(i)
    ENDDO

    DO i = 1,ncell
      rhoN(i) = rhoLiq (TN(i), pN(i))
      eN(i) = SpEnergy(TN(i), pN(i))
      rhoeN(i) = rhoN(i)*eN(i)
    ENDDO

    END SUBROUTINE EvalVar

    SUBROUTINE Residuals(mode)
!
!    Evaluate Mass, Energy, and Momentum equation residuals
!
!    Input
!
!    mode    -  0  only evaluate residuals
!               1  find maximum scaled residual
!
    USE Matrix
    USE FluidArrays
    USE Location
    USE FlowEqn
    USE ScalarDat
!
    IMPLICIT NONE
    INTEGER(sik), INTENT(IN) :: mode
    INTEGER(sik) :: i, ivar
!
    DO i = 1,nvar
      SELECT CASE (varList(i)%eqnType)
        CASE(1)                              !  Mass Equation
          b(i) =  MassEqn(varList(i)%iloc)
        CASE(2)                              !  Energy Equation
          b(i) =  EnergyEqn(varList(i)%iloc)
        CASE(3)                              !  Momentum Equation
          b(i) =  MomenEqn(varList(i)%iloc)
      END SELECT
    ENDDO
 
!
    IF(mode.EQ.0) RETURN
!
    ivar = 0
    residMax = 0.0_sdk
    IF(ivarV(1).NE.0) THEN
      ivar = ivar+1
      residMax = MAX(residMax,b(ivar)/MAX(0.01, v(1), vN(1))*dt)
    ENDIF
    DO i = 1,ncell
      ivar = ivar+1
      residMax = MAX(residMax,b(ivar)/(vol(i)*MAX(0.01, rho(i), rhoN(i)))*dt)
      ivar = ivar+1
      residMax = MAX(residMax,b(ivar)/(vol(i)*MAX(0.01, rhoe(i), rhoeN(i)))*dt)
      IF(ivarV(i+1).EQ.0) CYCLE
      ivar = ivar+1
      residMax = MAX(residMax,b(ivar)/MAX(0.01, v(i+1), vN(i+1))*dt)
    ENDDO
!
    IF(residMax.LT.eps) converged = .TRUE.
!
    END SUBROUTINE Residuals

    SUBROUTINE Jacobian
!
!    Evaluate the Jacobian Matrix using finite approximations to each
!    partial derivative
!
    USE Matrix
    USE FluidArrays
    USE ScalarDat
    USE FlowEqn
    USE Eos
    USE Location
!
    IMPLICIT NONE
    INTEGER(sik) :: i, ii, ivar, j
    REAL(sdk) :: fbase,varBase, rhoBase, rhoeBase, fracDvar=0.001_sdk, dvar
    REAL(sdk) :: fpert
!                                First clear the Jacobian Matrix
    a  = 0.0_sdk
!                                Now the set non-zero Jacobian elements
!                                Loop over all rows in  the matrix
    DO i = 1,nvar
      DO j = 1,varList(i)%n     !   Loop over all non-zero columns in the row
        ivar = varList(i)%i(j)
        IF(ivar.EQ.0) CYCLE
        varBase = loc(ivar)%var   !   hold the old variable value
!   
!       Perturb the variable
!
        IF(loc(ivar)%itype.EQ.3) THEN   !  velocity
          IF(loc(ivar)%var.GE.0) THEN
            dvar = MAX(fracDvar, fracDvar*loc(ivar)%var)
            loc(ivar)%var = loc(ivar)%var + dvar
          ELSE
            dvar = MIN(-fracDvar, fracDvar*loc(ivar)%var)
            loc(ivar)%var = loc(ivar)%var + dvar
          ENDIF
        ELSE                            !  pressure or temperature
            dvar = fracDvar*loc(ivar)%var
            loc(ivar)%var = loc(ivar)%var + dvar
            ii = loc(ivar)%i
            rhoBase = rhoN(ii)
            rhoeBase = rhoeN(ii)
            rhoN(ii) = RhoLiq(TN(ii),pN(ii))
            rhoeN(ii) = rhoN(ii)*SpEnergy(TN(ii),pN(ii))
        ENDIF
!                                 Reevaluate the function
        SELECT CASE (varList(i)%eqnType)
          CASE(1)                              !  Mass Equation
            fpert =  MassEqn(varList(i)%iloc)
          CASE(2)                              !  Energy Equation
            fpert =  EnergyEqn(varList(i)%iloc)
          CASE(3)                              !  Momentum Equation
            fpert =  MomenEqn(varList(i)%iloc)
        END SELECT
!                                                 Jacobian
        a(i,ivar) = (b(i) - fpert)/dvar
!
        loc(ivar)%var = varBase                !  Restore variables
        IF(loc(ivar)%itype.NE.3) THEN
          ii = loc(ivar)%i
          rhoN(ii) = rhoBase
          rhoeN(ii) = rhoeBase
        ENDIF
      ENDDO
    ENDDO

    END SUBROUTINE Jacobian

    SUBROUTINE SetWf
!
!   Set weighting factors for edge averages and differences
!   Contents of this subroutine are specific to the 1st order upwind
!   method used as a sample
!
    USE FluidArrays
    USE ScalarDat
    USE Location
!
    IMPLICIT NONE
    INTEGER(sik) :: i

    DO i = 1,ncell+1
      IF(vN(i).LT.0.0_sdk) THEN
        fluxAv(i)%wf(1) = 0.0_sdk
        fluxAv(i)%wf(2) = -0.125_sdk
        fluxAv(i)%wf(3) = 0.75_sdk
        fluxAv(i)%wf(4) = 0.375_sdk

        IF(i.LT.ivLB.OR.i.GT.ivUB) CYCLE
        dV(i)%wf(1) = 0.125_sdk
        dV(i)%wf(2) = 0.875_sdk
        dV(i)%wf(3) = 0.375_sdk
        dV(i)%wf(4) = 0.375_sdk
        dV(i)%wf(5) = 0.0_sdk
      ELSE
         fluxAv(i)%wf(1) = 0.375_sdk
         fluxAv(i)%wf(2) = 0.75_sdk
         fluxAv(i)%wf(3) = -0.125_sdk
         fluxAv(i)%wf(4) = 0.0_sdk
         IF(i.LT.ivLB.OR.i.GT.ivUB) CYCLE
         dV(i)%wf(1) = 0.0_sdk
         dV(i)%wf(2) = 0.375_sdk
         dV(i)%wf(3) = 0.75_sdk
         dV(i)%wf(4) = -0.125_sdk
         dV(i)%wf(5) = 0.0_sdk
      ENDIF
    ENDDO

    END SUBROUTINE SetWf


    END MODULE Trans
