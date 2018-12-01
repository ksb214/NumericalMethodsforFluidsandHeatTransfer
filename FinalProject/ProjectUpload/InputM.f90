    MODULE Input
    USE Intrtype
    USE ScalarDat
    USE FluidArrays
    USE Eos
!
!    This module sets initial conditions directly that would normally be
!    obtained from an input file
!
    CONTAINS

    SUBROUTINE SetIC
   ! dt = 0.1_sdk
    CALL AllocFluidAr

    DO i = LBOUND(v,DIM=1), UBOUND(v,DIM=1)
      vN(i) = 0.1_sdk
      fa(i) = (pi/4.0_sdk)*df**2
    ENDDO

    DO i = LBOUND(v,DIM=1), UBOUND(p,DIM=1)
      pN(i) = 20.0e5_sdk
      TN(i) = 300.0_sdk
      dx(i) = xsize/ncell
      vol(i) = dx(i)*fa(i)
    ENDDO
    
  
    DO zi = 1, ncell
       xloc= (zi-0.5)*dx(zi)
       Tm(zi)=  TN(1) + 4.0_sdk*7957.7_sdk*xloc/(df*uav*(rhoLiq(T(1),pN(1)))*cv) + 11.76_sdk
       q(zi)= hy*( Tm(zi) - TN(zi))*pi*dx(zi)*df/vol(zi)   
    end do

!    q(1:ncell)= hy*(  T(1)+ 4.0_sdk*7957.7_sdk/(2*r_i*uav)*(1/(rhoLiq (T(1), pN(i))*cv))*pi*dx*df/vol(1:ncell))
!    q(1:ncell)= hy*(Tm(1:ncell) - T(1:ncell))*pi*dx(1:ncell)*df/vol(1:ncell) 
    

    fric = 17.3_sdk               ! Constant value of the friction coefficient

    lBC = 'v'
    uBC = 'p'

    END SUBROUTINE SetIC

    END MODULE Input
