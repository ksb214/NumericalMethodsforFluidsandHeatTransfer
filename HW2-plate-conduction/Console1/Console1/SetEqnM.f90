    MODULE SetEqn
    USE IntrType
    USE Data

    CONTAINS

    SUBROUTINE SetCoeffs
!
!   Set coefficients in the finite volume conduction
!
!                 Programmed by Kamalesh S. Bhambare
!
    IMPLICIT NONE
!
    INTEGER(sik) :: i,j, k
!
!   Zero the coefficient matrix so that only non-zero coefficients need
!   to be set in the remainder of the subroutine
!
    a = 0.0_sdk
!
!   First calculate equation coefficients for all interior volumes
!   Loop over the index for the y location of the interior volumes (j) and 
!   Loop over the index for the x location of the interior volumes (i)
!   These indices are used for the convenience of isolating volumes that
!   don't have surface convective flux at any edge.
!
    DO j = 2,ny-1
      DO i = 2,nx-1
!
!       Set the unique volume index (row index in the matrix)
!
         k = (j-1)*nx + i
!         
         a(k,k) = -(2*cnd/dx**2 + 2*cnd/dy**2)
         a(k,k+1) = cnd/dx**2
         a(k,k-1) = cnd/dx**2
         a(k,k-nx) = cnd/dy**2
         a(k,k+nx) = cnd/dy**2
      ENDDO
    ENDDO
!
! Add code for volumes at the edges of the rod
!

    i=1
    j=1
    k = (j-1)*nx + i
    a(k,k) = -(3*cnd/dx**2 + 3*cnd/dy**2)+ (2*cnd/dx)**2*(1/(2*cnd + hx*dx)) +  (2*cnd/dy)**2*(1/(2*cnd + hy*dy))
    a(k,k+1) = cnd/dx**2
    a(k,k+nx) = cnd/dy**2
    

    i=nx
    j=1
    k = (j-1)*nx + i
    a(k,k) = -(3*cnd/dx**2 + 3*cnd/dy**2)+ (2*cnd/dx)**2*(1/(2*cnd + hx*dx)) +  (2*cnd/dy)**2*(1/(2*cnd + hy*dy))
    a(k,k-1) = cnd/dx**2
    a(k,k+nx) = cnd/dy**2

    i=1
    j=ny
    k = (j-1)*nx + i
    a(k,k) = -(3*cnd/dx**2 + 3*cnd/dy**2)+ (2*cnd/dx)**2*(1/(2*cnd + hx*dx)) +  (2*cnd/dy)**2*(1/(2*cnd + hy*dy))
    a(k,k+1) = cnd/dx**2
    a(k,k-nx) = cnd/dy**2

    i=nx
    j=ny
    k = (j-1)*nx + i
    a(k,k) = -(3*cnd/dx**2 + 3*cnd/dy**2)+ (2*cnd/dx)**2*(1/(2*cnd + hx*dx)) +  (2*cnd/dy)**2*(1/(2*cnd + hy*dy))
    a(k,k-1) = cnd/dx**2
    a(k,k-nx) = cnd/dy**2

! Code for the volumes adjacent to the wall

!Top wall
      j= ny
      DO i=2,nx-1
      k = (j-1)*nx + i
      a(k,k) = -(2*cnd/dx**2 + 3*cnd/dy**2)+ (2*cnd/dy)**2*(1/(2*cnd + hy*dy))
      a(k,k-1) = cnd/dx**2
      a(k,k+1) = cnd/dx**2
      a(k,k-nx) = cnd/dy**2   
      END DO

! Bottom wall
       j= 1
      DO i=2,nx-1
      k = (j-1)*nx + i
      a(k,k) = -(2*cnd/dx**2 + 3*cnd/dy**2)+ (2*cnd/dy)**2*(1/(2*cnd + hy*dy))
      a(k,k-1) = cnd/dx**2
      a(k,k+1) = cnd/dx**2
      a(k,k+nx) = cnd/dy**2   
      END DO

! RHS wall
      i= 1
      DO j=2,ny-1
      k = (j-1)*nx + i
      a(k,k) = -(3*cnd/dx**2 + 2*cnd/dy**2)+ (2*cnd/dx)**2*(1/(2*cnd + hx*dx))
      a(k,k+1) = cnd/dx**2
      a(k,k-nx) = cnd/dy**2   
      a(k,k+nx) = cnd/dy**2   
      END DO

! LHS wall
      i= nx
      DO j=2,ny-1
      k = (j-1)*nx + i
      a(k,k) = -(3*cnd/dx**2 + 2*cnd/dy**2)+ (2*cnd/dx)**2*(1/(2*cnd + hx*dx))
      a(k,k-1) = cnd/dx**2
      a(k,k-nx) = cnd/dy**2   
      a(k,k+nx) = cnd/dy**2   
      END DO

!
    END SUBROUTINE SetCoeffs

    SUBROUTINE SetRHS

!   Set the right hand side of the equation.  Note that when dton is 1.0
!   you get the right hand side of a transient equation.  When dton is 0.0
!   you get the right hand side of a steady state equation.
!
!                     Programmed by John Mahaffy 1/01
    IMPLICIT NONE
    INTEGER(sik) :: j, i, k, kp
!
!     Contributions common to all volumes
!
    DO j=1,ny
      DO i=1,nx
        k = (j-1)*nx + i
        b(k) = -p + dton*t(i,j) !  (dton allows use of coding for transients
      ENDDO
    ENDDO
!
! Adjustments for the y face boundary conditions
!
    i=1
    DO j= 2,ny-1
        k = (j-1)*nx + i
        b(k) = b(k) - 2*cnd*hx*tf/(2*cnd*dx + hx*dx**2)
     ENDDO

      i=nx
      DO j= 2,ny-1
        k = (j-1)*nx + i
        b(k) = b(k) - 2*cnd*hx*tf/(2*cnd*dx + hx*dx**2)
     ENDDO

! Adjustments for the x face boundary conditions

     j=1
    DO i= 2,nx-1
        k = (j-1)*nx + i
        b(k) = b(k) - 2*cnd*hy*tf/(2*cnd*dy + hy*dy**2)
     ENDDO

     j=ny
      DO i= 2,nx-1
        k = (j-1)*nx + i
        b(k) = b(k) - 2*cnd*hy*tf/(2*cnd*dy + hy*dy**2)
     ENDDO

! Adustments for the corners

     i=1
     j=1
     k = (j-1)*nx + i
     b(k) = b(k) - 2*cnd*hy*tf/(2*cnd*dy+ hy*dy**2) - 2*cnd*hx*tf/(2*cnd*dx + hx*dx**2)    

    i=nx
    j=1
    k = (j-1)*nx + i
    b(k) = b(k) - 2*cnd*hy*tf/(2*cnd*dy+ hy*dy**2) - 2*cnd*hx*tf/(2*cnd*dx + hx*dx**2)    

    i=1
    j=ny
    k = (j-1)*nx + i
    b(k) = b(k) - 2*cnd*hy*tf/(2*cnd*dy+ hy*dy**2) - 2*cnd*hx*tf/(2*cnd*dx + hx*dx**2)    

    i=nx
    j=ny
    k = (j-1)*nx + i
    b(k) = b(k) - 2*cnd*hy*tf/(2*cnd*dy+ hy*dy**2) - 2*cnd*hx*tf/(2*cnd*dx + hx*dx**2)    

 
    RETURN
    END SUBROUTINE SetRHS

    END MODULE
