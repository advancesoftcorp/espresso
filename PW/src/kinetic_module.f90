!
! Copyright (C) 2024 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE kinetic_module
  !--------------------------------------------------------------------------
  !
  ! ... the module for Kinetic Energy Density
  !
  USE cell_base,   ONLY : at, alat, omega
  USE constants,   ONLY : e2
  USE fft_base,    ONLY : dffts, dfftp
  USE kinds,       ONLY : DP
  USE io_files,    ONLY : tmp_dir, prefix
  USE io_global,   ONLY : ionode, stdout
  USE mp,          ONLY : mp_sum, mp_barrier
  USE mp_bands,    ONLY : intra_bgrp_comm
  USE mp_images,   ONLY : intra_image_comm
  USE scatter_mod, ONLY : gather_grid
  USE scf,         ONLY : rho
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL :: do_kinetic = .FALSE.
  INTEGER :: iunkinetic
  !
  PUBLIC :: do_kinetic
  PUBLIC :: kinetic_print
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE kinetic_open()
    !----------------------------------------------------------------------------
    !
    ! ... open file of Kinetic Energy Density
    !
    IMPLICIT NONE
    !
    INTEGER            :: ios
    CHARACTER(LEN=256) :: filename
    !
    INTEGER, EXTERNAL  :: find_free_unit
    !
    IF (.NOT. do_kinetic) THEN
      RETURN
    END IF
    !
    iunkinetic = find_free_unit()
    !
    filename = TRIM(tmp_dir) // TRIM(prefix) // '.ked'
    !
    IF (ionode) THEN
      !
      OPEN(unit=iunkinetic, file=TRIM(filename), &
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
    CALL errore('kinetic_open', 'cannot open file: ' // TRIM(filename), ios)
    !
  END SUBROUTINE kinetic_open
  !
  !----------------------------------------------------------------------------
  SUBROUTINE kinetic_close()
    !----------------------------------------------------------------------------
    !
    ! ... close file of Kinetic Energy Density
    !
    IMPLICIT NONE
    !
    LOGICAL :: opnd
    !
    IF (.NOT. do_kinetic) THEN
      RETURN
    END IF
    !
    IF (ionode) THEN
      !
      INQUIRE(unit=iunkinetic, opened=opnd)
      !
      IF (opnd) CLOSE(unit=iunkinetic)
      !
    END IF
    !
  END SUBROUTINE kinetic_close
  !
  !----------------------------------------------------------------------------
  SUBROUTINE kinetic_print()
    !----------------------------------------------------------------------------
    !
    ! ... print data for Kinetic Energy Density
    !
    IMPLICIT NONE
    !
    INTEGER  :: ir
    INTEGER  :: nx1, nx2, nx3
    INTEGER  :: nfft
    REAL(DP) :: fac
    REAL(DP) :: eneTauG
    REAL(DP) :: eneTauL
    REAL(DP) :: eneDtdr
    !
    REAL(DP), ALLOCATABLE :: tauG(:)
    REAL(DP), ALLOCATABLE :: tauL(:)
    REAL(DP), ALLOCATABLE :: dtdr(:)
    REAL(DP), ALLOCATABLE :: rhor_g(:)
    REAL(DP), ALLOCATABLE :: tauG_g(:)
    REAL(DP), ALLOCATABLE :: tauL_g(:)
    REAL(DP), ALLOCATABLE :: dtdr_g(:)
    !
    IF (.NOT. do_kinetic) THEN
      RETURN
    END IF
    !
    IF (ionode) THEN
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Writing Kinetic Energy Density to file.")')
      !
    END IF
    !
    IF (dfftp%nr1  /= dffts%nr1  .OR. dfftp%nr2  /= dffts%nr2  .OR. dfftp%nr3  /= dffts%nr3 .OR. \
        dfftp%nr1x /= dffts%nr1x .OR. dfftp%nr2x /= dffts%nr2x .OR. dfftp%nr3x /= dffts%nr3x) THEN
      !
      CALL errore('kinetic_print', 'Kinetic Energy Density does not support dual FFT-mesh', 1)
      !
    END IF
    !
    ! ... allocate memory
    !
    nx1  = dffts%nr1x
    nx2  = dffts%nr2x
    nx3  = dffts%nr3x
    nfft = nx1 * nx2 * nx3
    !
    ALLOCATE(tauG(dffts%nnr))
    ALLOCATE(tauL(dffts%nnr))
    ALLOCATE(dtdr(dffts%nnr))
    ALLOCATE(rhor_g(nfft))
    ALLOCATE(tauG_g(nfft))
    ALLOCATE(tauL_g(nfft))
    ALLOCATE(dtdr_g(nfft))
    !
    ! ... calculate Kinetic Energy Density
    !
    CALL kinetic_sum_band(tauG, tauL, dtdr)
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
      eneDtdr = eneDtdr + fac * dtdr(ir) * rho%of_r(ir, 1)
      !
    END DO
    !
    CALL mp_sum(eneTauG, intra_bgrp_comm)
    CALL mp_sum(eneTauL, intra_bgrp_comm)
    CALL mp_sum(eneDtdr, intra_bgrp_comm)
    !
    WRITE(stdout, '()')
    WRITE(stdout, '(5X,"Kinetic energy (by Gradient)  =",F17.8," Ry")') eneTauG
    WRITE(stdout, '(5X,"Kinetic energy (by Laplacian) =",F17.8," Ry")') eneTauL
    WRITE(stdout, '(5X,"Integral [ dT/drho * rho ]    =",F17.8," Ry")') eneDtdr
    !
#if defined(__MPI)
    rhor_g = 0.0_DP
    tauG_g = 0.0_DP
    tauL_g = 0.0_DP
    dtdr_g = 0.0_DP
    CALL gather_grid(dffts, rho%of_r(:, 1), rhor_g)
    CALL gather_grid(dffts, tauG, tauG_g)
    CALL gather_grid(dffts, tauL, tauL_g)
    CALL gather_grid(dffts, dtdr, dtdr_g)
#else
    rhor_g = rho%of_r(:, 1)
    tauG_g = tauG
    tauL_g = tauL
    dtdr_g = dtdr
#endif
    !
    ! ... print data, in Hartree unit
    !
    CALL kinetic_open()
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
      WRITE(iunkinetic, '("#Charge")')
      CALL density_print(iunkinetic, dffts%nr1x, dffts%nr2x, dffts%nr3x, 1.0_DP, rhor_g)
      !
      WRITE(iunkinetic, '("#Kinetic Energy Density (by Gradient)")')
      CALL density_print(iunkinetic, dffts%nr1x, dffts%nr2x, dffts%nr3x, 1.0_DP / e2, tauG_g)
      !
      WRITE(iunkinetic, '("#Kinetic Energy Density (by Laplacian)")')
      CALL density_print(iunkinetic, dffts%nr1x, dffts%nr2x, dffts%nr3x, 1.0_DP / e2, tauL_g)
      !
      WRITE(iunkinetic, '("#Kinetic Energy Derivative")')
      CALL density_print(iunkinetic, dffts%nr1x, dffts%nr2x, dffts%nr3x, 1.0_DP / e2, dtdr_g)
      !
    END IF
    !
    CALL kinetic_close()
    !
    CALL mp_barrier(intra_image_comm)
    !
    ! ... deallocate memory
    !
    DEALLOCATE(tauG)
    DEALLOCATE(tauL)
    DEALLOCATE(dtdr)
    DEALLOCATE(rhor_g)
    DEALLOCATE(tauG_g)
    DEALLOCATE(tauL_g)
    DEALLOCATE(dtdr_g)
    !
  END SUBROUTINE kinetic_print
  !
END MODULE kinetic_module
!
!----------------------------------------------------------------------------
SUBROUTINE density_print(iun, nr1x, nr2x, nr3x, fac, rhor)
  !----------------------------------------------------------------------------
  !
  USE fft_base, ONLY : dffts
  USE kinds,    ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: iun
  INTEGER,  INTENT(IN) :: nr1x, nr2x, nr3x
  REAL(DP), INTENT(IN) :: fac
  REAL(DP), INTENT(IN) :: rhor(nr1x, nr2x, nr3x)
  !
  INTEGER :: i1,  i2,  i3
  INTEGER :: nr1, nr2, nr3
  !
  nr1 = dffts%nr1
  nr2 = dffts%nr2
  nr3 = dffts%nr3
  !
  DO i1 = 1, nr1
    !
    DO i2 = 1, nr2
      !
      WRITE(iun,'(6E25.16)') (fac * rhor(i1, i2, i3), i3 = 1, nr3)
      !
    END DO
    !
  END DO
  !
END SUBROUTINE density_print

