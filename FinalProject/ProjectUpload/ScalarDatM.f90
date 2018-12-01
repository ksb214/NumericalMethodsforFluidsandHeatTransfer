    MODULE ScalarDat
    USE IntrType
!
!   Global scalar data needed by the program
!
    REAL(sdk) :: dt = 0.1_sdk       !   time step
    REAL(sdk) :: time =  0.0_sdk    !   time
    REAL(sdk) :: endTime = 10.0_sdk !   end time for the transient
    REAL(sdk) :: eps = 1.e-8_sdk    !   convergence criterion
    REAL(sdk) :: residMax           !   maximum factional residual
    INTEGER(sik) :: it              !   current iteration number
    INTEGER(sik) :: itmax = 10      !   number of iterations
    INTEGER(sik) :: nsmax = 100000  !   maximum number of time steps
    INTEGER(sik) :: nstep = 0       !   current time step number
    INTEGER(sik) :: ncell = 300      !   number of computational volumes
    INTEGER(sik) :: nvar            !   total number of independent variables
    INTEGER(sik) :: nbnd = 2        !   number of boundary cells
    INTEGER(sik) :: ivLB            !   lower bound of velocity indices
    INTEGER(sik) :: ivUB            !   upper bound of velocity indices
    CHARACTER(LEN=1) :: lBC, uBC    !   Boundary condition flags
                                    !   'p' for pressure boundary condition
                                    !   'v' for velocity boundary condition
    LOGICAL :: converged            !   test on iteration convergence
    INTEGER(sik)::z, zi 
    REAL(sdk) ::hy = 675.8_sdk      ! Heat transfer coeficient for laminar flow, from Nu = hD/k = 4.36
    REAL(sdk) :: xloc               
    REAL(sdk) :: xsize=2.0_sdk, uav = 0.1_sdk       ! Length of the tube  
    REAL(sdk) :: df= 0.004_sdk                      ! ID of the tube  
    REAL(sdk) :: pi=3.14159265358979_sdk
    END MODULE ScalarDat
