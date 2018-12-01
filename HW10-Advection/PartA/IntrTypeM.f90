      MODULE IntrType
        IMPLICIT NONE
!
!       These same definitions are repeated in all FUNCTION
!       declarations, for the NagWare F90 compiler
!
        INTEGER, PARAMETER :: sdk = selected_real_kind(13,307)
        INTEGER, PARAMETER :: sik = kind(10000000)
      END MODULE IntrType
