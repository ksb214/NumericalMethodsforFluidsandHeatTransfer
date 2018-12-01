    PROGRAM AdvectStep
!
!  Advect a step function using 1st Order Upwind, Leith, Quick, and Quickest
!
    USE Data
!
    CALL Setup                       !   Set initial conditions

    CALL Transient                   !   Run the transient

    CALL Output                      !   Write results to a file for plotting

    STOP
    END
