    MODULE Iterate3
    USE IntrType
    USE Data

    CONTAINS
!Subroutine for Jacobi method
    SUBROUTINE Jsolve3(a,b,n1,iter1,Tn1,t)

    IMPLICIT NONE

    INTEGER(sik) :: i ,j ,k
    INTEGER(sik),INTENT(IN) :: n1, iter1
    REAL(sdk),INTENT(IN OUT) :: a(:,:),b(:),Tn1(:), t(:)
    REAL(sdk):: sum, eps
    REAL(sdk), allocatable:: Tm1(:)
    ALLOCATE (Tm1(n)) !A dummy vector
       
    Tm1 = 0
    sum=0
    Tn=0
    do k = 1,iter1
      
       do i = 1, n1
          sum=0
          do j = 1, n1
             if(i.NE.j)then 
                sum= sum + a(i,j)*Tm1(j)
                
             end if
          end do
         
         Tn1(i) = (b(i)-sum)/a(i,i)
      end do
      eps=ERR1(t,Tn1,n1) !Call to a function to calculate the maximum error in the Tn1
      OPEN (7, FILE = 'Rplot_Jac_1450.dat', ACCESS = 'APPEND') 
      WRITE(7,*)k,LOG(eps)
      CLOSE(7)
  
      Tm1 = Tn1
      
   end do

  RETURN

    END SUBROUTINE Jsolve3
 
 !Subroutine Gauss Sidel   
    SUBROUTINE GS(a,b,n1,iter1,Tn2,t)
      IMPLICIT NONE

    INTEGER(sik) :: i ,j ,k
    INTEGER(sik),INTENT(IN) :: n1, iter1
    REAL(sdk),INTENT(IN OUT) :: a(:,:),b(:),Tn2(:),t(:)
    REAL(sdk):: sum,eps
    REAL(sdk), allocatable:: Tm2(:)
    ALLOCATE (Tm2(n))
       
    Tm2 = 0
    sum=0
    Tn2=0
    do k = 1,iter1
      
       do i = 1, n1
          sum=0
          do j = 1, i-1
             sum= sum + a(i,j)*Tn2(j)
          end do     
          
          do j = i+1,n1
          sum= sum + a(i,j)*Tm2(j)
         end do
          
          Tn2(i) = (b(i)-sum)/a(i,i)
       end do

       !Error calculation
     eps=ERR1(t,Tn2,n1)
     OPEN (7, FILE = 'Rplot_GS_1450.dat', ACCESS = 'APPEND')
     WRITE(7,*)k,LOG(eps)
     CLOSE(7)
    
     Tm2 = Tn2
      
  end do
         
   RETURN

    
    END SUBROUTINE GS
 


    !Gauss sidel other way   
    ! do j = 1, n1
    !         if(j.LT.i)then 
    !           sum= sum + a(i,j)*Tn2(j)
    !           end if
    !       if(j.GT.i)then 
    !         sum= sum + a(i,j)*Tm2(j)
    !     end if
    ! end do
    

!Subroutine SOR   

    SUBROUTINE SOR(a,b,n3,iter3,Tn3,t)
      IMPLICIT NONE

    INTEGER(sik) :: i ,j ,k 
    INTEGER(sik),INTENT(IN) :: n3, iter3
    REAL(sdk),INTENT(IN OUT) :: a(:,:),b(:),Tn3(:),t(:)
    REAL(sdk):: sum,w=1.528, eps
    REAL(sdk), allocatable:: Tm3(:)
    ALLOCATE (Tm3(n))
       
    Tm3 = 0
    Tn3 = 0
    do k = 1,iter3
      
       do i = 1, n3
        sum=0
          do j = 1, i-1
             sum= sum + a(i,j)*Tn3(j)
          end do     
          
          do j = i+1,n3
          sum= sum + a(i,j)*Tm3(j)
          end do
          
        sum= (b(i)-sum)/a(i,i)
        Tn3(i)= Tm3(i) + w*(sum-Tm3(i))
     end do
     
     !Error calculation
     eps=ERR1(t,Tn3,n3)
     OPEN (7, FILE = 'Rplot_SOR_1450.dat', ACCESS = 'APPEND')
     WRITE(7,*)k,LOG(eps)
     CLOSE(7)
  
     Tm3 = Tn3
      
  end do
   !Writting the output
   ! do i= 1,n
   !   WRITE(*,*)Tn3(i)
   ! end do
   
   RETURN
    
    END SUBROUTINE SOR
 
   PURE FUNCTION ERR1(te,it,n) result (Maxerr)
   IMPLICIT NONE
   INTEGER(sik),INTENT(IN) :: n
   REAL(sdk),INTENT(IN):: te(:),it(:)
   REAL(sdk), allocatable:: err(:)
   REAL(sdk):: Maxerr
   ALLOCATE (err(n))  !Vector to store the errors 
       
    err=te-it
    Maxerr = MAXVAL(err)
    Return

END FUNCTION ERR1




   END MODULE Iterate3
