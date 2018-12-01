    MODULE Eos
    USE IntrType
!
!   Equation of state information and other physical properties
!
    REAL(sdk) :: mu = 0.00086_sdk    ! water viscosity  (Pa*s)
    REAL(sdk) :: cond = 0.62_sdk     ! water conductivity (w/m/K)
    REAL(sdk) :: cv = 3720.76_sdk     ! water specific heat (J/kg/K)
    REAL(sdk) :: dRhoDt =-0.35809_sdk ! derivative of density with
                                     ! respect to temperature (kg/m**3/K)
    CONTAINS

    FUNCTION RhoLiq(T,p)
!
!     Water density
!
!   Input
!      T   - Temperature (K)
!      p   - pressure (Pa)
!   Output
!      RhoLiq   -  density  (kg/m**3)
!
    IMPLICIT NONE

    REAL(sdk), INTENT(IN) :: T, p
    REAL(sdk) :: RhoLiq
!
    RhoLiq = 999.7838 + dRhoDt*T
!
    END FUNCTION RhoLiq

    FUNCTION SpEnergy(T,p)
!
!     Water Specific Internal Energy
!
!   Input
!      T   - Temperature (K)
!      p   - pressure (Pa)
!   Output
!      SpEnergy  -  Specific Internal Energy (J/kg)
!
    IMPLICIT NONE

    REAL(sdk), INTENT(IN) :: T, p
    REAL(sdk) :: SpEnergy
!
    SpEnergy = 41990.0_sdk + cv*T
!
    END FUNCTION SpEnergy
 
    END MODULE Eos
