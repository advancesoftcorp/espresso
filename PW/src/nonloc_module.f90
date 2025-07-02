!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE nonloc_module
  !--------------------------------------------------------------------------
  !
  ! ... the module for Non-Local Energy Derivative
  !
  USE cell_base,     ONLY : at, alat
  USE constants,     ONLY : e2
  USE control_flags, ONLY : isolve, rmm_conv
  USE ener,          ONLY : ef
  USE fft_base,      ONLY : dffts, dfftp
  USE force_mod,     ONLY : lstres
  USE io_files,      ONLY : tmp_dir, prefix
  USE io_global,     ONLY : ionode, stdout
  USE ions_base,     ONLY : nat, ityp, atm, tau
  USE lsda_mod,      ONLY : nspin
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum, mp_barrier
  USE mp_images,     ONLY : intra_image_comm
  USE paw_variables, ONLY : okpaw
  USE sannp_module,  ONLY : do_sannp
  USE scatter_mod,   ONLY : gather_grid
  USE scf,           ONLY : vrs
  USE uspp_param,    ONLY : nhm
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL  :: do_nonloc      = .FALSE.
  REAL(DP) :: nonloc_kin_tf  = 1.0_DP
  REAL(DP) :: nonloc_kin_vw  = 0.2_DP
  LOGICAL  :: nonloc_kin_py  = .FALSE.
  REAL(DP) :: nonloc_perturb = 0.0_DP
  INTEGER  :: nonloc_nprint  = 0
  INTEGER  :: iunnonloc
  !
  PUBLIC :: do_nonloc
  PUBLIC :: nonloc_kin_tf
  PUBLIC :: nonloc_kin_vw
  PUBLIC :: nonloc_kin_py
  PUBLIC :: nonloc_perturb
  PUBLIC :: nonloc_nprint
  PUBLIC :: nonloc_print
  PUBLIC :: nonloc_add_perturb
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE nonloc_open(idx)
    !----------------------------------------------------------------------------
    !
    ! ... open file of Non-Local Energy Derivative
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: idx
    !
    INTEGER            :: ios
    CHARACTER(LEN=8)   :: str
    CHARACTER(LEN=256) :: filename
    !
    INTEGER, EXTERNAL  :: find_free_unit
    !
    IF (.NOT. do_nonloc) THEN
      RETURN
    END IF
    !
    iunnonloc = find_free_unit()
    !
    IF (idx >= 0) THEN
      !
      WRITE(str, '(I8)') idx
      filename = TRIM(tmp_dir) // TRIM(prefix) // '.nld.' // TRIM(ADJUSTL(str))
      !
    ELSE
      !
      filename = TRIM(tmp_dir) // TRIM(prefix) // '.nld'
      !
    END IF
    !
    IF (ionode) THEN
      !
      OPEN(unit=iunnonloc, file=TRIM(filename), &
         & status='unknown', form='formatted', action='write', iostat=ios)
      !
      ios = ABS(ios)
      !
    ELSE
      !
      ios = 0
      !
    END IF
    !
    CALL mp_sum(ios, intra_image_comm)
    !
    CALL errore('nonloc_open', 'cannot open file: ' // TRIM(filename), ios)
    !
  END SUBROUTINE nonloc_open
  !
  !----------------------------------------------------------------------------
  SUBROUTINE nonloc_close()
    !----------------------------------------------------------------------------
    !
    ! ... close file of Non-Local Energy Derivative
    !
    IMPLICIT NONE
    !
    LOGICAL :: opnd
    !
    IF (.NOT. do_nonloc) THEN
      RETURN
    END IF
    !
    IF (ionode) THEN
      !
      INQUIRE(unit=iunnonloc, opened=opnd)
      !
      IF (opnd) CLOSE(unit=iunnonloc)
      !
    END IF
    !
  END SUBROUTINE nonloc_close
  !
  !----------------------------------------------------------------------------
  SUBROUTINE nonloc_print(idx)
    !----------------------------------------------------------------------------
    !
    ! ... print data for Non-Local Energy Derivative
    !
    IMPLICIT NONE
    !
    INTEGER, OPTIONAL, INTENT(IN) :: idx
    !
    INTEGER  :: idx_
    INTEGER  :: ia, it
    INTEGER  :: ir1, ir2, ir3
    INTEGER  :: nr1, nr2, nr3
    INTEGER  :: nr1x, nr2x, nr3x
    INTEGER  :: nfft
    !
    REAL(DP), ALLOCATABLE :: becsum(:,:,:) ! \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
    REAL(DP), ALLOCATABLE :: eneNL (:)
    REAL(DP), ALLOCATABLE :: dvdr  (:)
    REAL(DP), ALLOCATABLE :: dvdr_g(:)
    !
    IF (PRESENT(idx)) THEN
      idx_ = idx
    ELSE
      idx_ = -1
    END IF
    !
    IF (.NOT. do_nonloc) THEN
      RETURN
    END IF
    !
    IF (ionode .AND. idx_ < 0) THEN
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Writing Non-Local Energy and its Derivative to file.")')
      !
    END IF
    !
    IF (nspin /= 1) THEN
      !
      CALL errore('nonloc_print', 'Non-Local Energy Derivative supports only nspin = 1', 1)
      !
    END IF
    !
    IF (lstres) THEN
      !
      CALL errore('nonloc_print', 'Non-Local Energy Derivative does not support Stress', 1)
      !
    END IF
    !
    IF (okpaw) THEN
      !
      CALL errore('nonloc_print', 'Non-Local Energy Derivative does not support PAW', 1)
      !
    END IF
    !
    IF (do_sannp) THEN
      !
      CALL errore('nonloc_print', 'Non-Local Energy Derivative does not support SANNP', 1)
      !
    END IF
    !
    IF (isolve == 4 .AND. .NOT. rmm_conv) THEN
      !
      CALL errore('nonloc_print', &
      'diago_rmm_conv must be .TRUE., when calculate Non-Local Energy Derivative', 1)
      !
    END IF
    !
    ! ... allocate memory
    !
    nr1  = dffts%nr1
    nr2  = dffts%nr2
    nr3  = dffts%nr3
    nr1x = dffts%nr1x
    nr2x = dffts%nr2x
    nr3x = dffts%nr3x
    nfft = nr1x * nr2x * nr3x
    !
    ALLOCATE(becsum(nhm, nhm, nat)) ! w/o spin
    ALLOCATE(eneNL (nat))
    ALLOCATE(dvdr  (dffts%nnr))
    ALLOCATE(dvdr_g(nfft))
    !
    ! ... calculate Non-Local Energy and its Derivative
    !
    CALL occmat_sum_band(becsum)
    !
    CALL nonloc_energy(eneNL, becsum)
    !
    IF (ionode .AND. idx_ < 0) THEN
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,A)') '  # atom     Enl (Ry)'
      !
      DO ia = 1, nat
        !
        WRITE(stdout, '(5X,I3,2X,A4,0PE14.5)') &
        & ia, ADJUSTL(atm(ityp(ia))) // '    ', eneNL(ia)
        !
      END DO
      !
    END IF
    !
    dvdr(1:dffts%nnr) = ef - vrs(1:dffts%nnr, 1)
    !
#if defined(__MPI)
    dvdr_g = 0.0_DP
    CALL gather_grid(dffts, dvdr, dvdr_g)
#else
    dvdr_g = dvdr
#endif
    !
    ! ... print data, in Hartree unit
    !
    CALL nonloc_open(idx_)
    !
    IF (ionode) THEN
      !
      WRITE(iunnonloc, '("#Mesh")')
      WRITE(iunnonloc, "(3I8)") nr1, nr2, nr3
      !
      WRITE(iunnonloc, '("#Lattice")')
      WRITE(iunnonloc, '(3E25.16)') alat * at(1, 1), alat * at(2, 1), alat * at(3, 1)
      WRITE(iunnonloc, '(3E25.16)') alat * at(1, 2), alat * at(2, 2), alat * at(3, 2)
      WRITE(iunnonloc, '(3E25.16)') alat * at(1, 3), alat * at(2, 3), alat * at(3, 3)
      !
      WRITE(iunnonloc, '("#Number of Atoms")')
      WRITE(iunnonloc, '(I5)') nat
      WRITE(iunnonloc, '("#Atoms")')
      !
      DO ia = 1, nat
        !
        it = ityp(ia)
        !
        WRITE(iunnonloc, '(I5,A6,4E25.16)') ia, atm(it), &
        alat * tau(1, ia), alat * tau(2, ia), alat * tau(3, ia), eneNL(ia) / e2
        !
      END DO
      !
      WRITE(iunnonloc, '("#Non-Local Energy Derivative")')
      !
      DO ir1 = 1, nr1
        !
        DO ir2 = 1, nr2
          !
          WRITE(iunnonloc,'(6E25.16)') &
          & (dvdr_g(ir1 + (ir2 - 1) * nr1x + (ir3 - 1) * nr1x * nr2x) / e2, ir3 = 1, nr3)
          !
        END DO
        !
      END DO
      !
    END IF
    !
    CALL nonloc_close()
    !
    CALL mp_barrier(intra_image_comm)
    !
    ! ... deallocate memory
    !
    DEALLOCATE(becsum)
    DEALLOCATE(eneNL)
    DEALLOCATE(dvdr)
    DEALLOCATE(dvdr_g)
    !
  END SUBROUTINE nonloc_print
  !
  !----------------------------------------------------------------------------
  SUBROUTINE nonloc_add_perturb(aux)
    !----------------------------------------------------------------------------
    !
    ! ... add perturbation potential to local potential, in G-space
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: aux(dfftp%nnr)
    !
    IF (.NOT. do_nonloc) THEN
      RETURN
    END IF
    !
    IF (nonloc_perturb > 0.0_DP) THEN
      !
      CALL add_local_perturb(nonloc_perturb, aux)
      !
    END IF
    !
  END SUBROUTINE nonloc_add_perturb
  !
END MODULE nonloc_module
