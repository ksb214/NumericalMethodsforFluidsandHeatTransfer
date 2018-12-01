    MODULE FlowEqn
!
!   This contains functions needed to evaluate flow equation residuals
!
    USE IntrType
    USE Location
    USE FluidArrays
    USE ScalarDat
    Use Eos

    CONTAINS

    FUNCTION MassEqn(i)
!
!   Residual of  the mass equation at cell i
!
    IMPLICIT NONE
!
!   Input
!
!   i   -   cell index at which the equation is evaluated
!
    REAL(sdk) :: MassEqn
    INTEGER(sik), INTENT(IN) :: i

    MassEqn = (rhoN(i) - rho(i))*vol(i)/dt + Flux(rhoN,i+1) -  &
                                         Flux(rhoN,i)

    END FUNCTION MassEqn
!
    FUNCTION EnergyEqn(i)
!
!   Residual of  the energy equation at cell i
!
    IMPLICIT NONE
!
!   Input
!
!   i   -   cell index at which the equation is evaluated
!
    REAL(sdk) :: EnergyEqn
    INTEGER(sik), INTENT(IN) :: i

    EnergyEqn = (rhoeN(i)- rhoe(i))*vol(i)/dt +   &
                Flux(rhoeN,i+1) - Flux(rhoeN,i)   &
                + pN(i)*(fa(i+1)*vN(i+1)-fa(i)*vN(i)) -q(i)*vol(i)

    END FUNCTION EnergyEqn
 
    FUNCTION MomenEqn(i)
!
!   Residual of  the momentum equation at face i (low edge of cell i)
!
    IMPLICIT NONE
!
!   Input
!
!   i   -   edge index at which the equation is evaluated
!
    REAL(sdk) :: MomenEqn
    INTEGER(sik), INTENT(IN) :: i

    MomenEqn = (vN(i)- v(i))/dt + Vdv(i) + GradP(i)/EdgeAv(rhoN,i) +  &
                fric(i)*vN(i)*ABS(vN(i))

    END FUNCTION MomenEqn
!
    FUNCTION EdgeAv(x,j)
!
!    Calculate the value of scalar field x at cell face j based on
!    weighting coefficients in fluxAv
!
    IMPLICIT NONE
    REAL(sdk) :: EdgeAv
!                            quantity being fluxed
    REAL(sdk), INTENT(IN) :: x(LBOUND(rho,1):UBOUND(rho,1))          
    INTEGER(sik), INTENT(IN) :: j          !  edge index
    TYPE (averageT), POINTER :: a
    INTEGER(sik) :: i

    a => fluxAv(j)
    EdgeAv = 0
    DO i = 1,a%n
      EdgeAv = EdgeAv + a%wf(i)*x(a%i(i))
    ENDDO

    END FUNCTION EdgeAv
!
    FUNCTION Flux(x,j)
!
!    Calculate the flux of scalar field x at cell face j based on
!    weighting coefficients in fluxAv
!
    IMPLICIT NONE
    REAL(sdk) :: Flux
!                         Quantity being fluxed           
    REAL(sdk), INTENT(IN) :: x(LBOUND(rho,1):UBOUND(rho,1))
    INTEGER(sik), INTENT(IN) :: j                !  edge index
    TYPE (averageT), POINTER :: a
    INTEGER(sik) :: i

    a => fluxAv(j)
    Flux = 0                             ! First calculate edge value of x
    DO i = 1,a%n
      Flux = Flux + a%wf(i)*x(a%i(i))
    ENDDO

    Flux = Flux*vN(j)*fa(j)    !  Multiply edge value by edge area and velocity

    END FUNCTION Flux
!
    FUNCTION Vdv(j)
!
!    Calculate the momentum flux term at face j based upon weighting 
!    information in derived type array dV
!    
    IMPLICIT NONE

    REAL(sdk) :: Vdv
    INTEGER(sik), INTENT(IN) :: j  ! Face where momentum Equation is evaluated
    INTEGER(sik) :: i
    TYPE (averageT), POINTER :: a

    a => dV(j)
    Vdv = 0
    DO i = 1,a%n
      Vdv = Vdv + a%wf(i)*vN(a%i(i))
    ENDDO
      Vdv = vN(j)*Vdv

    END FUNCTION Vdv

    FUNCTION GradP(j)
!
!    Calculate an approximation to the pressure gradient at face j
!    based on information in the dPdx array
!
    IMPLICIT NONE
    REAL(sdk) :: GradP
    TYPE (averageT), POINTER :: a
    INTEGER(sik), INTENT(IN) :: j  ! Face where momentum Equation is evaluated
    INTEGER(sik) :: i

    a => dPdx(j)
    GradP = 0
    DO i = 1,a%n
      GradP = GradP + a%wf(i)*pN(a%i(i))
    ENDDO

    END FUNCTION GradP

    END MODULE FlowEqn
