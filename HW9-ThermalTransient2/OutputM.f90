    MODULE Output
    USE IntrType
    USE Data

    CHARACTER*1 :: tab = ACHAR(9)   ! easy way use tabs in output
!
!     Subroutines associated with output
!
    CONTAINS
!
    SUBROUTINE edit
!
!     Edit key data 
!     Calculate surface conditions, and energy balance
!     Output is to UNIT 11
!
!                Programmed by John Mahaffy  1/01
!
    IMPLICIT NONE
    INTEGER :: i,j
    REAL(sdk) :: qsurf, qvol, tmax, tmin
    LOGICAL :: isOpen
!
!                                     OPEN the output file
    INQUIRE (UNIT=11,OPENED=isOpen)
    IF(.NOT.isopen) OPEN(11,FILE='trans_CN_0.4.txt')
!
! Surface temperatures and heat loss at the surface
!
!    txs(1:ny,1) = (3*hx*tf+9*cnd/dx*t(1,1:ny)-cnd/dx*t(2,1:ny))/(3*hx+8*cnd/dx)
!    txs(1:ny,2) = (3*hx*tf+9*cnd/dx*t(nx,1:ny)-cnd/dx*t(nx-1,1:ny))/(3*hx+8*cnd/dx)

 !   tys(1:nx,1) = (3*hy*tf+9*cnd/dy*t(1:nx,1)-cnd/dy*t(1:nx,2))/(3*hy+8*cnd/dy)
 !   tys(1:nx,2) = (3*hy*tf+9*cnd/dy*t(1:nx,ny)-cnd/dy*t(1:nx,ny-1))/(3*hy+8*cnd/dy)

 !   qsurf = SUM(hx*(txs(1:ny,1)-tf)*dy) +  &
 !          SUM(hx*(txs(1:ny,2)-tf)*dy) +  &              
 !           SUM(hy*(tys(1:nx,1)-tf)*dx) +  &
 !           SUM(hy*(tys(1:nx,2)-tf)*dx) 
 !   qvol =  q*xsize*ysize                   !  Total heat generation
!
!       Check for a match of total predicted surface heat flux with
!       total heat produced in the region.
!
!    IF (numStepsTried .EQ. 0) THEN
!       WRITE(11,2000) qsurf, qvol
!2000   FORMAT(/,5x, 'Surface Energy Loss', 3x,'Energy Generated', /,   &
!              1x, 1P, 2E22.12)
!       RETURN
!    END IF

!    WRITE(11,20) time, tab, MAXVAL(abs(tss-t)), tab, tss(ny/2+1,nx/2+1)-t(ny/2+1,nx/2+1)!, tab,MINVAL(t), &
                 !tab, qsurf, tab, qvol
    WRITE(11,20) time, tab,t(ny/2+1,nx/2+1)
20  FORMAT(1x, F7.3, 2(A, F18.12))!, 2(A, 1P, E22.12)) 

!  OPEN(12,FILE='transient2.out',POSITION='APPEND')
!    WRITE(12,'(a,/)') 'time'
!    WRITE(12,22)time
!    WRITE(12,22)
!    WRITE(12,'(a,/)') 'Temperature is'
!    DO j = 1,n
!    WRITE(12,*)tt(j)
!    ENDDO
!22  FORMAT(10F8.4)
 
!    Tana = (q/(rho*cp))*time + tf  
!    !Uncomment this for comparison with the analytical solution with zero heat flux   
!    OPEN(14,FILE='analytical.out',POSITION='APPEND')
!    WRITE(14,19)time, Tana, t(ny/2+1,nx/2+1)
!19  FORMAT(1X,f7.3,2(4X,E22.12)) 


!
!    Tabs are added to data to ease pickup by spreadsheets
!

!
    RETURN
    END SUBROUTINE edit

    SUBROUTINE TimeErrorEstimate
    IMPLICIT NONE
    INTEGER(sik) :: i, iCenter
    REAL(sdk) :: avgErr, cenErr, maxErr, rmsErr
!     
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

    SUBROUTINE SpatialRichardson
!
!   Calculate estimate of spatial order of accuracy and error estimates
!   for each mesh.
!   This requires evaluation on at least 3 spatial meshes, and monotonic
!   error change.  If error is not monotonic the value returned for p is 0
!
!   Programmed by John Mahaffy PSU 1/05
!
    IMPLICIT NONE
    INTEGER(sik) :: i
    REAL(sdk) :: pC, pWx, pWy, errC, errWx, errWy, exactC, exactWx, exactWy
    REAL(sdk) :: error
!
     DO i = 3, numMeshTried
!
!      Center point
!
       pC = (tCenter(i-2)-tCenter(i-1))/(tCenter(i-1)-tCenter(i))
       IF (pC .GT. 0.0_sdk) THEN
          pC = log(pC)/log(rSpace)
          exactC = tCenter(i)-(tCenter(i-1)-tCenter(i))/(rSpace**pC-1)
          errC = tCenter(i)-exactC
       ELSE
          pC = 0.0_sdk
          exactC = 999
          errC = 999
       ENDIF
!
!     Center of a surface facing in the negative x direction     
!
       pWx = (txwc(i-2)-txwc(i-1))/(txwc(i-1)-txwc(i))
       IF (pWx .GT. 0.0_sdk) THEN
          pWx = log(pWx)/log(rSpace)
          exactWx = txwc(i)-(txwc(i-1)-txwc(i))/(rSpace**pWx-1)
          errWx = txwc(i)-exactWx
       ELSE
          pWx = 0.0_sdk
          exactWx = 999
          errWx = 999
       ENDIF
!
!     Center of a surface facing in the negative y direction     
!
       pWy = (tywc(i-2)-tywc(i-1))/(tywc(i-1)-tywc(i))
       IF (pWy .GT. 0.0_sdk) THEN
          pWy = log(pWy)/log(rSpace)
          exactWy = tywc(i)-(tywc(i-1)-tywc(i))/(rSpace**pWy-1)
          errWy = tywc(i)-exactWy
       ELSE
          pWy = 0.0_sdk
          exactWy = 999
          errWy = 999
       ENDIF
       IF (i .EQ. 3) THEN
          WRITE (*,*)  'Richardson Based Analysis for Mesh'
          WRITE (*,'(/,a,i3,a)')'   For ',nCells(1),' mesh cells in each direction'
          error = tCenter(1)-exactC
          WRITE (*,'(1p,3(a,e15.7))')'   Center temperature ', tCenter(1), ', error ', &
                error 
          WRITE (*,'(1p,3(a,e15.7))')'      Order ', pC, ', "exact =', exactC 
          error = txwc(1)-exactWx
          WRITE (*,'(1p,3(a,e15.7))')'X Surface temperature ', txwc(1), ', error ', &
                error 
          WRITE (*,'(1p,3(a,e15.7))')'      Order ', pWx, ', "exact =', exactWx 
          error = tywc(1)-exactWy
          WRITE (*,'(1p,3(a,e15.7))')'Y Surface temperature ', tywc(1), ', error ', &
             error 
          WRITE (*,'(1p,3(a,e15.7))')'      Order ', pWy, ', "exact =', exactWy 
!
          WRITE (*,'(/,a,i3,a)')'   For ',nCells(2),' mesh cells in each direction'
          error = tCenter(2)-exactC
          WRITE (*,'(1p,3(a,e15.7))')'   Center temperature ', tCenter(2), ', error ', &
             error 
          WRITE (*,'(1p,3(a,e15.7))')'      Order ', pC, ', "exact =', exactC 
          error = txwc(2)-exactWx
          WRITE (*,'(1p,3(a,e15.7))')'X Surface temperature ', txwc(2), ', error ', &
             error 
          WRITE (*,'(1p,3(a,e15.7))')'      Order ', pWx, ', "exact =', exactWx 
          error = tywc(2)-exactWy
          WRITE (*,'(1p,3(a,e15.7))')'Y Surface temperature ', tywc(2), ', error ', &
                error 
          WRITE (*,'(1p,3(a,e15.7))')'      Order ', pWy, ', "exact =', exactWy 

          WRITE (22,*)  'Richardson Based Analysis for Mesh'
          WRITE (22,'(/,a,i3,a)')'   For ',nCells(1),' mesh cells in each direction'
          error = tCenter(1)-exactC
          WRITE (22,'(1p,3(a,e15.7))')'   Center temperature ', tCenter(1), ', error ', &
                error 
          WRITE (22,'(1p,3(a,e15.7))')'      Order ', pC, ', "exact =', exactC 
          error = txwc(1)-exactWx
          WRITE (22,'(1p,3(a,e15.7))')'X Surface temperature ', txwc(1), ', error ', &
                error 
          WRITE (22,'(1p,3(a,e15.7))')'      Order ', pWx, ', "exact =', exactWx 
          error = tywc(1)-exactWy
          WRITE (22,'(1p,3(a,e15.7))')'Y Surface temperature ', tywc(1), ', error ', &
             error 
          WRITE (22,'(1p,3(a,e15.7))')'      Order ', pWy, ', "exact =', exactWy 
!
          WRITE (22,'(/,a,i3,a)')'   For ',nCells(2),' mesh cells in each direction'
          error = tCenter(2)-exactC
          WRITE (22,'(1p,3(a,e15.7))')'   Center temperature ', tCenter(2), ', error ', &
             error 
          WRITE (22,'(1p,3(a,e15.7))')'      Order ', pC, ', "exact =', exactC 
          error = txwc(2)-exactWx
          WRITE (22,'(1p,3(a,e15.7))')'X Surface temperature ', txwc(2), ', error ', &
             error 
          WRITE (22,'(1p,3(a,e15.7))')'      Order ', pWx, ', "exact =', exactWx 
          error = tywc(2)-exactWy
          WRITE (22,'(1p,3(a,e15.7))')'Y Surface temperature ', tywc(2), ', error ', &
                error 
          WRITE (22,'(1p,3(a,e15.7))')'      Order ', pWy, ', "exact =', exactWy 
       END IF
!
       WRITE (*,*)  'Richardson Based Analysis for Mesh'
       WRITE (*,'(/,a,i3,a)')'   For ',nCells(i),' mesh cells in each direction'
       error = tCenter(i)-exactC
       WRITE (*,'(1p,3(a,e15.7))')'   Center temperature ', tCenter(i), ', error ', &
             error 
       WRITE (*,'(1p,3(a,e15.7))')'      Order ', pC, ', "exact =', exactC 
       error = txwc(i)-exactWx
       WRITE (*,'(1p,3(a,e15.7))')'X Surface temperature ', txwc(i), ', error ', &
             error 
       WRITE (*,'(1p,3(a,e15.7))')'      Order ', pWx, ', "exact =', exactWx 
       error = tywc(i)-exactWy
       WRITE (*,'(1p,3(a,e15.7))')'Y Surface temperature ', tywc(i), ', error ', &
          error 
       WRITE (*,'(1p,3(a,e15.7))')'      Order ', pWy, ', "exact =', exactWy 

       WRITE (22,*)  'Richardson Based Analysis for Mesh'
       WRITE (22,'(/,a,i3,a)')'   For ',nCells(i),' mesh cells in each direction'
       error = tCenter(i)-exactC
       WRITE (22,'(1p,3(a,e15.7))')'   Center temperature ', tCenter(i), ', error ', &
             error 
       WRITE (22,'(1p,3(a,e15.7))')'      Order ', pC, ', "exact =', exactC 
       error = txwc(i)-exactWx
       WRITE (22,'(1p,3(a,e15.7))')'X Surface temperature ', txwc(i), ', error ', &
             error 
       WRITE (22,'(1p,3(a,e15.7))')'      Order ', pWx, ', "exact =', exactWx 
       error = tywc(i)-exactWy
       WRITE (22,'(1p,3(a,e15.7))')'Y Surface temperature ', tywc(i), ', error ', &
          error 
       WRITE (22,'(1p,3(a,e15.7))')'      Order ', pWy, ', "exact =', exactWy 
!

    END DO
    END SUBROUTINE SpatialRichardson

    END MODULE Output



