      MODULE Data
      USE IntrType
! Variable declaration part
      REAL(sdk), ALLOCATABLE :: a(:,:),aa(:,:),b(:), h(:),atra(:,:),btra(:),aatra(:,:),tttra(:),inte(:),temp(:),temp_n(:),t2(:,:)
      REAL(sdk), ALLOCATABLE :: r(:),re(:),t(:,:),tt(:),txs(:,:),tys(:,:)
      REAL(sdk), ALLOCATABLE :: tx(:), ty(:)

      INTEGER(sik), ALLOCATABLE :: ipvt(:)
      REAL(sdk) :: pi=3.14159265358979_sdk

      REAL(sdk) :: xsize, ysize , cnd, tf, dx, dy, q, t11ss, time
      REAL(sdk) :: cp, rho, delt, p, d, s, dton=0.0_sdk
      REAL(sdk) :: hx,hy, Tanal, tstart, tend, dtra, delttra, ptra, qtra
      INTEGER(sik):: nx, ny, n, itmax

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

      nx = 99                 !  number of volumes in the x direction
      ny = 99               !  number of volumes in the y direction

      hx = 2.5e4_sdk         !  heat transfer coefficient on faces at fixed x
      hy = 2.5e4_sdk          !  heat transfer coefficient on faces at fixed y

      n=nx*ny                !  Total number of volumes
!
      ALLOCATE (a(n,n),h(4),aa(n,n),b(n),ipvt(n))
      ALLOCATE (t(nx,ny),tt(n),txs(ny,2), tys(nx,2), tx(nx), ty(ny),aatra(n,n),tttra(n),atra(n,n),btra(n),temp(n),inte(n),temp_n(n),t2(nx,ny))

      a = 0.0                !  Clear the entire matrix
      b = 0.0                !  Clear the right hand side array

      dx=xsize/nx
      dy=ysize/ny
      temp=tf
      
! Setting up transient data
      
      qtra= 3.0e+09_sdk 
      delttra=(0.5_sdk*dx**2*dy**2*rho*cp)/(cnd*(dx**2 + dy**2))    
      tstart = 0.0_sdk
      tend = 5.0_sdk

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

!

    OPEN(11,FILE='transient.out',POSITION='APPEND')
    WRITE(11,'(a,/)') 'time'
    WRITE(11,22)time
    WRITE(11,22)
    WRITE(11,'(a,/)') 'Temperature is'
    DO j = 1,n
    WRITE(11,*)temp(j)
    ENDDO
22  FORMAT(10F8.4)
    close (11)

  ! Uncomment this for comparison with the analytical solution with zero heat flux   
!       OPEN(14,FILE='analytical.out',POSITION='APPEND')
!       WRITE(14,19)time, Tanal, t2(ny/2+1,nx/2+1)
!19     FORMAT(1X,f7.3,2(4X,E22.12)) 

!    Writting down the abs absolute maximum difference between the transient solution and steady state with time 
!    OPEN(12,FILE='diff_9by9_1.05.out',POSITION='APPEND')
!    WRITE(12,10)time,maxval(abs(tt-temp))
!10  FORMAT(1X,f7.3,4X,E22.12)
!   close (12)
                                                                             
    ! Writting down the center temperature and time   
!    OPEN(13,FILE='center_9by9_1.05.out',POSITION='APPEND')
!    WRITE(13,18)time,t(ny/2+1,nx/2+1)-t2(ny/2+1,nx/2+1)
!18  FORMAT(1X,f7.3,4X,E22.12)
!    close (13)

    RETURN
    END SUBROUTINE edit

 
  END MODULE Data



