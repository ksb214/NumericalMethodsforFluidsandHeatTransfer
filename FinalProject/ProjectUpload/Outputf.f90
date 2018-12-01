       MODULE OUTPUTF
         USE IntrType
         USE Eos
         USE Data
         
         !   Subprograms used in the flow transient for printing the data to file
         !
       CONTAINS
         
         SUBROUTINE EDITF
           
           USE ScalarDat
           USE FluidArrays
           
           IMPLICIT NONE

           INTEGER(sik)::kc,j
           LOGICAL :: isOpen
           CHARACTER*1 :: tab = ACHAR(9)
           REAL(sdk)::min, mout, qsurf, q_gain 

           min = rho(1)*fa(1)*vN(1)
           mout = rho(ncell)*fa(ncell)*vN(ncell)
         


             INQUIRE (UNIT=91,OPENED=isOpen)
           IF(.NOT.isopen) OPEN(91,FILE='HB_steadyF.txt')
           Write(91,*)'Inlet mass flow is',min
           Write(91,*)'Outlet mass flow is',mout
           
           q_gain = min*cv*(TN(ncell)-TN(1))
           Write(91,*)'Total heat gained by the fluid',q_gain         

           qsurf = sum((hy*(Tm(1:ncell) - TN(1:ncell))*pi*dx(1:ncell)*df))                    
           Write(91,*)'Total heat transfer to the fluid',qsurf
           Write(91,*)'Percentage imbalance',(qsurf-q_gain)/qsurf*100
910        FORMAT(1x, E18.12)
           

! Writting the flow variables to the file
 
           INQUIRE (UNIT=99,OPENED=isOpen)
           IF(.NOT.isopen) OPEN(99,FILE='ansR0.03.txt')
           DO kc = 1, ncell
              Write(99,10)(kc-0.5)*dx(kc),tab,T(kc)
10            FORMAT(1x, F17.4, (A, E18.12))
           END DO

! Block for writting steady state metal temperature
           OPEN the output file
           INQUIRE (UNIT=11,OPENED=isOpen)
           IF(.NOT.isopen)OPEN(11,FILE='wall_variation.txt')                                
           Do j = 1,nx
              WRITE(11,30)(j-0.5_sdk)*dxm, tab, Tm(j)                 
           END DO
30         FORMAT(1x, E18.10, 3(A, E18.12))

           RETURN
         END SUBROUTINE EDITF
       END MODULE OUTPUTF
