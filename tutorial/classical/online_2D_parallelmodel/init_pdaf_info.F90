!$Id: init_pdaf_info.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: init_pdaf_info - Screen output on assimilation configuration
!
! !INTERFACE:
SUBROUTINE init_pdaf_info()

! !DESCRIPTION:
! This routine performs a model-sided screen output about
! the coniguration of the data assimilation system.
! Using this output is optional. Most of the information
! is also displayed by PDAF itself when it is initialized
! in PDAF_init. Not displayed by PDAF is the assimilation
! interval (delt_obs), which is unknown to PDAF.
!
! !REVISION HISTORY:
! 2011-05 - Lars Nerger - Initial code extracted from init_pdaf
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation, & ! Variables for assimilation
       ONLY: filtertype, subtype, dim_ens, delt_obs, model_error, &
       model_err_amp, forget, rank_ana_enkf, int_rediag

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: init_pdaf
!EOP


! *****************************
! *** Initial Screen output ***
! *****************************

  IF (filtertype == 0) THEN
     WRITE (*, '(/21x, a)') 'Filter: SEEK'
     IF (subtype == 2) THEN
        WRITE (*, '(6x, a)') '-- fixed basis filter with update of matrix U'
        WRITE (*, '(6x, a)') '-- no re-diagonalization of VUV^T'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(6x, a)') '-- fixed basis filter & no update of matrix U'
        WRITE (*, '(6x, a)') '-- no re-diagonalization of VUV^T'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(13x, a, i5)') 'number of EOFs:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (subtype /= 5) THEN
        IF ((int_rediag > 0) .AND. ((subtype /= 2) .OR. (subtype /= 3))) &
             WRITE (*, '(10x, a, i4, a)') &
             'Re-diag each ', int_rediag, '-th analysis step'
     ELSE
        IF (int_rediag == 1) THEN
           WRITE (*, '(10x, a)') 'Perform re-diagonalization'
        ELSE
           WRITE (*, '(10x, a)') 'No re-diagonalization'
        END IF
     END IF
  ELSE IF (filtertype == 1) THEN
     WRITE (*, '(21x, a)') 'Filter: SEIK'
     IF (subtype == 2) THEN
        WRITE (*, '(6x, a)') '-- fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(6x, a)') '-- use ensemble transformation'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 2) THEN
     WRITE (*, '(21x, a)') 'Filter: EnKF'
     IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
     IF (rank_ana_enkf > 0) THEN
        WRITE (*, '(6x, a, i5)') &
             'analysis with pseudo-inverse of HPH, rank:', rank_ana_enkf
     END IF
  ELSE IF (filtertype == 3) THEN
     WRITE (*, '(21x, a)') 'Filter: LSEIK'
     IF (subtype == 2) THEN
        WRITE (*, '(6x, a)') '-- fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(6x, a)') '-- fixed state covariance matrix'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(6x, a)') '-- use ensemble transformation'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 4) THEN
     WRITE (*, '(21x, a)') 'Filter: ETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 5) THEN
     WRITE (*, '(21x, a)') 'Filter: LETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(6x, a)') '-- Variant following Hunt et al. (2007)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 6) THEN
     WRITE (*, '(21x, a)') 'Filter: ESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 7) THEN
     WRITE (*, '(21x, a)') 'Filter: LESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(6x, a)') '-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(6x, a)') '-- Offline mode'
     END IF
     WRITE (*, '(14x, a, i5)') 'ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(6x, a, i5)') 'Assimilation interval:', delt_obs
     WRITE (*, '(10x, a, f5.2)') 'forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(6x, a, f5.2)') 'model error amplitude:', model_err_amp
     END IF
  END IF     

END SUBROUTINE init_pdaf_info
