      MODULE Data
      USE IntrType
! Variable declaration part
      REAL(sdk), ALLOCATABLE :: a(:,:),aa(:,:),b(:), h(:)
      REAL(sdk), ALLOCATABLE :: r(:),re(:),t(:,:),tt(:),txs(:,:),tys(:,:)
      REAL(sdk), ALLOCATABLE :: tx(:), ty(:)

      INTEGER(sik), ALLOCATABLE :: ipvt(:)
      REAL(sdk) :: pi=3.14159265358979_sdk

      REAL(sdk) :: xsize, ysize , cnd, tf, dx, dy, q, t11ss, time
      REAL(sdk) :: cp, rho, delt, p, d, s, dton=0.0_sdk
      REAL(sdk) :: hx,hy
      INTEGER(sik):: nx, ny, n

      CONTAINS

      SUBROUTINE SetData
! 
!     Basic description of the problem
!     This is derived from a Transient equation solver, so has some extras
!
      IMPLICIT NONE

      ysize =0.03_sdk        !  height of the bar
      xsize = 0.004_sdk      !  width of the bar
      cnd=85.0_sdk           !  conductivity
      tf=300.0_sdk           !  fluid temperature
      cp = 500.0_sdk         !  specific heat (not needed for this problem)
      rho = 8000.0_sdk       !  density       (not needed for this problem)
      q=3.0e+09_sdk          !  power source per unit volume

      delt = 1.0_sdk         !  time step

      s = delt/(cp*rho)      !  scale factor for transient temperature change 
      s = 1.0_sdk            !  Steady state

      d = cnd*s              !  coefficient of diffusion per time step

      p = q*s                !  temperature change per step due to power source

      nx = 27                 !  number of volumes in the x direction
      ny = 27              !  number of volumes in the y direction

      hx = 2.5e4_sdk         !  heat transfer coefficient on faces at fixed x
      hy = 2.5e4_sdk          !  heat transfer coefficient on faces at fixed y

      n=nx*ny                !  Total number of volumes
!
      ALLOCATE (a(n,n),h(4),aa(n,n),b(n),ipvt(n))
      ALLOCATE (t(nx,ny),tt(n),txs(ny,2), tys(nx,2), tx(nx), ty(ny))

      a = 0.0                !  Clear the entire matrix
      b = 0.0                !  Clear the right hand side array

      dx=xsize/nx
      dy=ysize/ny

      END SUBROUTINE SetData

    SUBROUTINE edit
!
!     Edit key data 
!     Calculate surface conditions, and energy balance
!     Output is to UNIT 11
!
!                Programmed by John Mahaffy  1/01
!
    IMPLICIT NONE
    CHARACTER*1 :: tab = ACHAR(9)   ! easy way use tabs in output
    INTEGER :: i,j
    REAL(sdk) :: qsurf, qvol
!
!                                     OPEN the output file
    OPEN(11,FILE='steady.out')
!
! Surface temperatures and heat loss at the surface

    txs(1:ny,1) = (9*cnd*t(1,1:ny)+3*dx*hx*tf-cnd*t(2,1:ny))/(8*cnd+3*hx*dx) !Left side
    txs(1:ny,2) = (9*cnd*t(nx,1:ny)+3*dx*hx*tf-cnd*t(nx-1,1:ny))/(8*cnd+3*hx*dx)  !Right side
    tys(1:nx,1) = (9*cnd*t(1:nx,1)+3*dy*hy*tf-cnd*t(1:nx,2))/(8*cnd+3*hy*dy)      !Top  
    tys(1:nx,2) = (9*cnd*t(1:nx,ny)+3*dy*hy*tf-cnd*t(1:nx,ny-1))/(8*cnd+3*hy*dy)  !Bottom
    !   txs(1:ny,1) = (hx*tf+2*cnd/dx*t(nx,1:ny))/(hx+2*cnd/dx)
    !   txs(1:ny,2) = (hx*tf+2*cnd/dx*t(nx,1:ny))/(hx+2*cnd/dx)
    !   tys(1:nx,1) = (hy*tf+2*cnd/dy*t(1:nx,1))/(hy+2*cnd/dy)
    !   tys(1:nx,2) = (hy*tf+2*cnd/dy*t(1:nx,ny))/(hy+2*cnd/dy)


    qsurf = SUM(hx*(txs(1:ny,1)-tf)*dy) +  &
            SUM(hx*(txs(1:ny,2)-tf)*dy) +  &              
            SUM(hy*(tys(1:nx,1)-tf)*dx) +  &
            SUM(hy*(tys(1:nx,2)-tf)*dx) 
    qvol =  q*xsize*ysize                   !  Total heat generation
!
!       Check for a match of total predicted surface heat flux with
!       total heat produced in the region.
!
    WRITE(11,20) qsurf, qvol
20  FORMAT(/,5x, 'Surface Energy Loss', 3x,'Energy Generated', /,   &
           1x, 1P, 2E22.12)
!
!   Special Edit for solution with no y dependence
!
    IF(hx.NE.0.and.hy.EQ.0) THEN
      WRITE(11,'(//)')
      tx = tf + q*xsize/(2.0_sdk*hx) + q/(2*cnd)*(xsize/2.0_sdk)**2 -   &
             q/(2*cnd)*(( (/(i,i=1,nx)/)-0.5_sdk)*dx - xsize/2.0_sdk)**2
      WRITE(11,*) ' x analytic solution'
      WRITE(11,'("  x(m)  ", "    T(K)",/,(f8.5,f9.3))')    &
             ( (i-0.5_sdk)*dx, tx(i), i=1,nx)  
!      WRITE(11,21) tx(1:nx)
      WRITE(11,'(//)')
    ENDIF
!
!   Special Edit for solution with no x dependence
!
    IF(hy.NE.0.and.hx.EQ.0) THEN
      WRITE(11,'(//)')
      ty = tf + q*ysize/(2*hy) + q/(2*cnd)*(ysize/2.0_sdk)**2 -   &
           q/(2*cnd)*(((/(i,i=1,ny)/)-0.5_sdk)*dy - ysize/2.0_sdk)**2
       WRITE(11,*) ' y analytic solution'
      WRITE(11,21) ty(1:ny)
      WRITE(11,'(//)')
    ENDIF
!
!   Print full temperature field
!
    WRITE(11,'(a,/)') 'Numerical Solution'

    DO j = 1,ny
      WRITE(11,21) t(1:nx,j)
    ENDDO
21  FORMAT(10F8.2)
    
    WRITE(11,'(a,/)') 'Temperature at the rod center  '
     WRITE(11,21)t((ny+1)/2,(nx+1)/2)
  ! WRITE(11,21)t(ny/2+1,nx/2+1)
    
    WRITE(11,'(a,/)') 'Temperature at the mid of the bottom face   '
    WRITE(11,21)tys((nx+1)/2,1)
    
    WRITE(11,'(a,/)') 'Temperature at the mid of RHS side   '
    WRITE(11,21)txs((ny+1)/2,2)

    WRITE(11,'(a,/)') 'dx is    '
    WRITE(11,*)dx

     WRITE(11,'(a,/)') 'dy is    '
     WRITE(11,*)dy

    RETURN
    END SUBROUTINE edit


      END MODULE Data



