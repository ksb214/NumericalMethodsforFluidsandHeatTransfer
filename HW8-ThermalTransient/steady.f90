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
    


CALL CPU_TIME(tstart)
! Transient Solution                              
    itmax = (tend-tstart)/delttra
    aatra =  atra                      !  Save a copy of the original equations for
    tttra =  btra                      !  later comparisons during debugging
    time = 0                           !  Setting the initial time zero 
    
    do i= 1,itmax                                  ! Transient loop         
       inte = matmul(aatra,temp)
       temp_n = inte + tttra
       temp = temp_n
       t2 = RESHAPE(temp,(/nx,ny/))
       Tanal = (q/(rho*cp))*time + tf           ! Analytical solution    
       CALL edit    
       time= time + delttra  
    end do
    CALL CPU_TIME(tend)
    WRITE(*,*)   'Direct Solver CPU = ',tend-tstart,' seconds for ', n, 'volumes'
    STOP                        !  Output of key variables

    END

