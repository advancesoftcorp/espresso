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
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL  :: do_nonloc      = .FALSE.
  REAL(DP) :: nonloc_kin_tf  = 1.0_DP
  REAL(DP) :: nonloc_kin_vw  = 0.2_DP
  REAL(DP) :: nonloc_perturb = 0.0_DP
  INTEGER  :: nonloc_nprint  = 0
  INTEGER  :: iunnonloc
  !
  PUBLIC :: do_nonloc
  PUBLIC :: nonloc_coef_tf
  PUBLIC :: nonloc_coef_vw
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
    INTEGER  :: ir
    INTEGER  :: nr1x, nr2x, nr3x
    INTEGER  :: nfft
    REAL(DP) :: fac
    REAL(DP) :: eneTauG
    REAL(DP) :: eneTauL
    REAL(DP) :: eneDtdr
    !
    REAL(DP), ALLOCATABLE :: rhor(:)
    REAL(DP), ALLOCATABLE :: tauG(:)
    REAL(DP), ALLOCATABLE :: tauL(:)
    REAL(DP), ALLOCATABLE :: dtdr(:)
    REAL(DP), ALLOCATABLE :: rho1_g(:)
    REAL(DP), ALLOCATABLE :: rho2_g(:)
    REAL(DP), ALLOCATABLE :: tauG_g(:)
    REAL(DP), ALLOCATABLE :: tauL_g(:)
    REAL(DP), ALLOCATABLE :: dtdr_g(:)
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
      WRITE(stdout, '(5X,"Writing Non-Local Energy Derivative to file.")')
      !
    END IF
    !
    IF (nspin /= 1) THEN
      !
      CALL errore('nonloc_print', 'Non-Local Energy Derivative supports only nspin = 1', 1)
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
    nr1x = dffts%nr1x
    nr2x = dffts%nr2x
    nr3x = dffts%nr3x
    nfft = nr1x * nr2x * nr3x
    !
    ALLOCATE(rhor(dffts%nnr))
    ALLOCATE(tauG(dffts%nnr))
    ALLOCATE(tauL(dffts%nnr))
    ALLOCATE(dtdr(dffts%nnr))
    ALLOCATE(rho1_g(nfft))
    ALLOCATE(rho2_g(nfft))
    ALLOCATE(tauG_g(nfft))
    ALLOCATE(tauL_g(nfft))
    ALLOCATE(dtdr_g(nfft))
    !
    ! ... calculate Kinetic Energy Density
    !
    CALL kinetic_sum_band(rhor, tauG, tauL, dtdr, .TRUE.)
    !
    eneTauG = 0.0_DP
    eneTauL = 0.0_DP
    eneDtdr = 0.0_DP
    !
    fac = omega / DBLE(dffts%nr1 * dffts%nr2 * dffts%nr3)
    !
    DO ir = 1, dffts%nnr
      !
      eneTauG = eneTauG + fac * tauG(ir)
      eneTauL = eneTauL + fac * tauL(ir)
      eneDtdr = eneDtdr + fac * dtdr(ir) * rhor(ir)
      !
    END DO
    !
    CALL mp_sum(eneTauG, intra_bgrp_comm)
    CALL mp_sum(eneTauL, intra_bgrp_comm)
    CALL mp_sum(eneDtdr, intra_bgrp_comm)
    !
    IF (ionode .AND. idx_ < 0) THEN
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Kinetic energy (by Gradient)  =",F17.8," Ry")') eneTauG
      WRITE(stdout, '(5X,"Kinetic energy (by Laplacian) =",F17.8," Ry")') eneTauL
      WRITE(stdout, '(5X,"Integral [ dT/drho * rho ]    =",F17.8," Ry")') eneDtdr
      !
    END IF
    !
#if defined(__MPI)
    rho1_g = 0.0_DP
    rho2_g = 0.0_DP
    tauG_g = 0.0_DP
    tauL_g = 0.0_DP
    dtdr_g = 0.0_DP
    CALL gather_grid(dffts, rhor,           rho1_g)
    CALL gather_grid(dffts, rho%of_r(:, 1), rho2_g)
    CALL gather_grid(dffts, tauG,           tauG_g)
    CALL gather_grid(dffts, tauL,           tauL_g)
    CALL gather_grid(dffts, dtdr,           dtdr_g)
#else
    rho1_g = rhor
    rho2_g = rho%of_r(:, 1)
    tauG_g = tauG
    tauL_g = tauL
    dtdr_g = dtdr
#endif
    !
    ! ... print data, in Hartree unit
    !
    CALL kinetic_open(idx_)
    !
    IF (ionode) THEN
      !
      WRITE(iunkinetic, '("#Mesh")')
      WRITE(iunkinetic, "(3I8)") dffts%nr1, dffts%nr2, dffts%nr3
      !
      WRITE(iunkinetic, '("#Lattice")')
      WRITE(iunkinetic, '(3E25.16)') alat * at(1, 1), alat * at(2, 1), alat * at(3, 1)
      WRITE(iunkinetic, '(3E25.16)') alat * at(1, 2), alat * at(2, 2), alat * at(3, 2)
      WRITE(iunkinetic, '(3E25.16)') alat * at(1, 3), alat * at(2, 3), alat * at(3, 3)
      !
      WRITE(iunkinetic, '("#Including Kinetic Energy Derivative")')
      !
      IF (with_dtdr) THEN
        WRITE(iunkinetic, "(I8)") 1
      ELSE
        WRITE(iunkinetic, "(I8)") 0
      END IF
      !
      WRITE(iunkinetic, '("#Charge (IN)")')
      CALL density_print(iunkinetic, nr1x, nr2x, nr3x, 1.0_DP, rho2_g)
      !
      WRITE(iunkinetic, '("#Charge (OUT)")')
      CALL density_print(iunkinetic, nr1x, nr2x, nr3x, 1.0_DP, rho1_g)
      !
      WRITE(iunkinetic, '("#Kinetic Energy Density (by Gradient)")')
      CALL density_print(iunkinetic, nr1x, nr2x, nr3x, 1.0_DP / e2, tauG_g)
      !
      WRITE(iunkinetic, '("#Kinetic Energy Density (by Laplacian)")')
      CALL density_print(iunkinetic, nr1x, nr2x, nr3x, 1.0_DP / e2, tauL_g)
      !
      IF (with_dtdr) THEN
        !
        WRITE(iunkinetic, '("#Kinetic Energy Derivative")')
        CALL density_print(iunkinetic, nr1x, nr2x, nr3x, 1.0_DP / e2, dtdr_g)
        !
      END IF
      !
    END IF
    !
    CALL kinetic_close()
    !
    CALL mp_barrier(intra_image_comm)
    !
    ! ... deallocate memory
    !
    DEALLOCATE(rhor)
    DEALLOCATE(tauG)
    DEALLOCATE(tauL)
    DEALLOCATE(dtdr)
    DEALLOCATE(rho1_g)
    DEALLOCATE(rho2_g)
    DEALLOCATE(tauG_g)
    DEALLOCATE(tauL_g)
    DEALLOCATE(dtdr_g)
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

