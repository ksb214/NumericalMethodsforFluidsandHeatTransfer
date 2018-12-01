    MODULE Moutput
    USE IntrType
    USE Data
    USE FluidArrays

    CHARACTER*1 :: tab = ACHAR(9)   ! easy way use tabs in output
    CONTAINS

    SUBROUTINE medit
!                Programmed by John Mahaffy  1/01
!
    IMPLICIT NONE
    INTEGER :: i,j
    REAL(sdk) :: qsurf, qvol, tmax, tmin, tmana, tfana
    LOGICAL :: isOpen

! Surface temperatures and heat loss at the surface
    INQUIRE (UNIT=25,OPENED=isOpen)
    IF(.NOT.isopen) OPEN(25,FILE='HB_steadyM.txt')
    txs(1:ny,1) = (3*hx*tf+9*cnd/dxm*ts(1,1:ny)-cnd/dxm*ts(2,1:ny))/(3*hx+8*cnd/dxm)
    txs(1:ny,2) = (3*hx*tf+9*cnd/dxm*ts(nx,1:ny)-cnd/dxm*ts(nx-1,1:ny))/(3*hx+8*cnd/dxm)
    tys(1:nx,1) = (3*hy*tf+9*cnd/dy*ts(1:nx,1)-cnd/dy*ts(1:nx,2))/(3*hy+8*cnd/dy)
    tys(1:nx,2) = (3*hy*tf+9*cnd/dy*ts(1:nx,ny)-cnd/dy*ts(1:nx,ny-1))/(3*hy+8*cnd/dy)
    
    qsurf = SUM(hy*(ts(1:nx,1)-txka(1:nx))*dxm)*2*pi*r_i
    qvol =  qm*pi*((r_i+ysize)**2-r_i**2)*xsize                   !  Total heat generation
    
!       Check for a match of total predicted surface heat flux with
!       total heat produced in the region.
!
    IF (numStepsTried .EQ. 1) THEN
       WRITE(25,2000) qsurf, qvol
2000   FORMAT(/,5x, 'Surface Energy Loss', 3x,'Energy Generated', /,   &
              1x, 1P, 2E22.12)
    END IF
    close(25)

! Block for writting the analytical solution in x direction to the file for the comparison
    INQUIRE (UNIT=100,OPENED=isOpen)
    IF(.NOT.isopen) OPEN(100,FILE='Analytical.txt')
    Do j=1,nx
       tfana = 400 + 25.29_sdk*(j-0.5_sdk)*dxm
       tmana= tfana + 11.76_sdk
       WRITE(100,101)(j-0.5_sdk)*dxm, tab,tmana, tab, tfana
    END DO
101 FORMAT(1x, F7.3, 2(A, E18.12))
    close(100)



! Block for writting the solution in r direction to the file for comparison

    INQUIRE (UNIT=100,OPENED=isOpen)
    IF(.NOT.isopen) OPEN(100,FILE='Analyticalr.txt')
     Do j=1,ny
        r_o = r_i+(j-0.5_sdk)*dy
        tmana = (qm/(4.0_sdk*cnd))*(r_i**2-r_o**2) + (qm/(2.0_sdk*cnd))*(r_ot**2)*(log(r_o/r_i)) + (qm/(2.0_sdk*hy))*((r_ot**2-r_i**2)/r_i) + TN(1)    
       WRITE(100,102) r_o, tab, tmana
    END DO
102 FORMAT(1x, E18.10, (A, E18.12))
    close(100)



!   OPEN the output file
!    INQUIRE (UNIT=11,OPENED=isOpen)
!    IF(.NOT.isopen)OPEN(11,FILE='wall_variation.txt')    !,'POSITION='APPEND')  ! 'transient30.txt'      ! Use this file for writting the outlet wall temperature
!    Do j = 1,nx
!       WRITE(11,20)(j-0.5_sdk)*dxm, tab, Ts(j,1)                                                   !    Block for writting the wall temperature variation along x
!    END DO




! this is for writting @outlet cell temperature at every transient iteration for checking steady state

!    WRITE(11,20)time,tab,(9.0_sdk*ts(1,ny)-ts(1,ny-1))/8.0_sdk,tab,(9.0_sdk*ts(nx,ny)-ts(nx,ny-1))/8.0_sdk,tab,(9.0_sdk*ts(nx/2,ny)-ts(nx/2,ny-1))/8.0_sdk !t(nx,ny),tab(nx/2,ny)       

! Block for writting the temperature variation in r direction


    IF(.NOT.isopen)OPEN(11,FILE='r_variation.txt')
    Do j = 1,ny
        r_o = r_i+(j-0.5_sdk)*dy
       WRITE(11,20) r_o, tab, ts(1,j)                 !    Block for writting the wall temperature variation along x
    END DO

20  FORMAT(1x, E18.10, 3(A, E18.12))!, 2(A, 1P, E22.12)) 
    close (11)
       


!    Tabs are added to data to ease pickup by spreadsheets
!

!
    RETURN
  END SUBROUTINE medit

    SUBROUTINE TimeErrorEstimate
    IMPLICIT NONE
    INTEGER(sik) :: i, iCenter
    REAL(sdk) :: avgErr, cenErr, maxErr, rmsErr
    LOGICAL :: isOpen
!     
    INQUIRE (UNIT=22,OPENED=isOpen)
    IF(.NOT.isopen) OPEN(22,FILE='Richardson_exp.txt')
    IF (numStepsTried .LE. 1) RETURN
    WRITE(22,'(/,a,i3)') 'Time Error Estimates for nx=ny=',nx
    iCenter = n/2 + 1
    WRITE(22,*) 'dt        ',tab,'cenErr    ',tab,'maxErr    ',tab,'rmsErr    ',tab,'avgErr'
!
!   Compare results of larger timesteps to those of the smallest
!
    DO i = 1, numStepsTried
       maxErr = MAXVAL(ABS(tLast(:,i)-tLast(:,1)))
       rmsErr = SQRT(SUM((tLast(:,i)-tLast(:,1))**2)/n)
       avgErr = ABS(SUM(tLast(:,i)-tLast(:,1))/n)
       cenErr = ABS(tLast(iCenter,i)-tLast(iCenter,1))
       WRITE(22,'(f8.5,1P,4(a,e10.4))')  &
            timeSteps(i),tab,cenErr,tab,maxErr,tab,rmsErr,tab,avgErr
    END DO

    END SUBROUTINE TimeErrorEstimate

  
  END MODULE Moutput



