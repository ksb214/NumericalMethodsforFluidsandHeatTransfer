    SUBROUTINE Output
    USE IntrType
    USE Data

    IMPLICIT NONE
    INTEGER(sik) :: i

    OPEN(10, FILE='rho_0.5_100.txt')
    WRITE(10,'(1x,a8,6(a,a9))') '#    x   ',tab,'Donor-Exp',tab,'Donor-Imp'
    WRITE(10,'(1x,F8.3,6(A,F9.5))') 0.0_sdk,tab,rhoDonorExp(0), & 
            tab, rhoDonorImp(0)
    DO i = 1,n
      WRITE(10,'(1x,F8.3,6(A,F9.5))') (i-0.5_sdk)*dx,tab,rhoDonorExp(i), & 
            tab, rhoDonorImp(i)
    ENDDO


    STOP
    END


