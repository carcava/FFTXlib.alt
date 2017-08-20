MODULE fft_helper_subroutines

  IMPLICIT NONE
  SAVE

  INTERFACE tg_reduce_rho
    MODULE PROCEDURE tg_reduce_rho_1,tg_reduce_rho_2,tg_reduce_rho_3
  END INTERFACE

CONTAINS


  SUBROUTINE tg_reduce_rho_1( rhos, tg_rho_nc, tg_rho, ispin, noncolin, domag, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(IN) :: ispin
     LOGICAL, INTENT(IN) :: noncolin, domag
     REAL(DP), INTENT(INOUT)  :: tg_rho(:)
     REAL(DP), INTENT(INOUT)  :: tg_rho_nc(:,:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, ioff, idx, ir, ipol

     IF ( desc%nogrp > 1 ) THEN
#ifdef __MPI
        IF( noncolin) THEN
           CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_rho_nc, SIZE(tg_rho_nc), MPI_DOUBLE_PRECISION, MPI_SUM, desc%ogrp_comm, ierr )
        ELSE
           CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_rho, SIZE(tg_rho), MPI_DOUBLE_PRECISION, MPI_SUM, desc%ogrp_comm, ierr )
        END IF
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     !If the current processor is not the "first" processor in its
     !orbital group then does a local copy (reshuffling) of its data
     !
     ioff = 0
     DO idx = 1, desc%nogrp
        IF( desc%mype == desc%nolist( idx ) ) EXIT
        ioff = ioff + desc%nr1x * desc%nr2x * desc%npp( desc%nolist( idx ) + 1 )
     END DO
     !
     ! copy the charge back to the processor location
     !
     IF (noncolin) THEN
!$omp parallel do
        DO ir = 1, desc%nnr
           rhos(ir,1) = rhos(ir,1) + tg_rho_nc(ir+ioff,1)
        END DO
!$omp end parallel do
        IF (domag) THEN
!$omp parallel do
           DO ipol=2,4
              DO ir = 1, desc%nnr
                 rhos(ir,ipol) = rhos(ir,ipol) + tg_rho_nc(ir+ioff,ipol)
              END DO
           END DO
!$omp end parallel do
        ENDIF
     ELSE
!$omp parallel do
        DO ir = 1, desc%nnr
           rhos(ir,ispin) = rhos(ir,ispin) + tg_rho(ir+ioff)
        END DO
!$omp end parallel do
     END IF

  END SUBROUTINE



  SUBROUTINE tg_reduce_rho_2( rhos, tmp_rhos, ispin, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(IN) :: ispin
     REAL(DP), INTENT(INOUT)  :: tmp_rhos(:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, ioff, idx, ir

     IF ( desc%nogrp > 1 ) THEN
#ifdef __MPI
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%ogrp_comm, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     !If the current processor is not the "first" processor in its
     !orbital group then does a local copy (reshuffling) of its data
     !
     ioff = 0
     DO idx = 1, desc%nogrp
        IF( desc%mype == desc%nolist( idx ) ) EXIT
        ioff = ioff + desc%nr1x * desc%nr2x * desc%npp( desc%nolist( idx ) + 1 )
     END DO
     !
     ! copy the charge back to the processor location
     !
     DO ir = 1, desc%nnr
        rhos(ir,ispin) = rhos(ir,ispin) + tmp_rhos(ir+ioff)
     END DO

  END SUBROUTINE


  SUBROUTINE tg_reduce_rho_3( rhos, tmp_rhos, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     REAL(DP), INTENT(INOUT)  :: tmp_rhos(:,:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, from, ii, ir

     IF ( desc%nogrp > 1 ) THEN
#ifdef __MPI
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%ogrp_comm, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     !If the current processor is not the "first" processor in its
     !orbital group then does a local copy (reshuffling) of its data
     !

     from = 1
     DO ii = 1, desc%nogrp
        IF ( desc%nolist( ii ) == desc%mype ) EXIT !Exit the loop
        from = from +  desc%nr1x*desc%nr2x*desc%npp( desc%nolist( ii ) + 1 )! From where to copy initially
     ENDDO
     !
     DO ir = 1, SIZE(rhos,2)
         CALL dcopy( desc%nr1x*desc%nr2x*desc%npp(desc%mype+1), tmp_rhos(from,ir), 1, rhos(1,ir), 1)
     ENDDO

  END SUBROUTINE


  SUBROUTINE tg_get_nnr( desc, right_nnr )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: right_nnr
     right_nnr = desc%tg_nnr
  END SUBROUTINE


  SUBROUTINE tg_get_local_nr3( desc, val )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%npl
  END SUBROUTINE

  SUBROUTINE tg_get_group_nr3( desc, val )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     IF( desc%nogrp > 1 ) THEN
        val = desc%tg_npp( desc%mype + 1 )
     ELSE
        val = desc%my_nr3p
     END IF
  END SUBROUTINE

  SUBROUTINE tg_get_recip_inc( desc, val )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%nr3x * desc%nsw(desc%mype+1)
  END SUBROUTINE


END MODULE fft_helper_subroutines
