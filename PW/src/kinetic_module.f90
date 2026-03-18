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
  USE cell_base,      ONLY : at, alat, omega
  USE constants,      ONLY : e2
  USE control_flags,  ONLY : isolve, rmm_conv
  USE fft_base,       ONLY : dffts, dfftp
  USE fft_types,      ONLY : fft_index_to_3d
  USE gvect,          ONLY : ngm, igtongl, mill, eigts1, eigts2, eigts3
  USE kinds,          ONLY : DP
  USE io_files,       ONLY : tmp_dir, prefix
  USE io_global,      ONLY : ionode, stdout
  USE ions_base,      ONLY : nat, ityp, zv, atm
  USE lsda_mod,       ONLY : nspin
  USE mp,             ONLY : mp_sum, mp_barrier
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE mp_images,      ONLY : intra_image_comm
  USE scatter_mod,    ONLY : gather_grid
  USE scf,            ONLY : rho
  USE vlocal,         ONLY : vloc
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL  :: do_kinetic      = .FALSE.
  LOGICAL  :: kinetic_dtdr    = .TRUE.
  REAL(DP) :: kinetic_perturb = 0.0_DP
  INTEGER  :: kinetic_nprint  = 0
  INTEGER  :: iunkinetic
  !
  PUBLIC :: do_kinetic
  PUBLIC :: kinetic_dtdr
  PUBLIC :: kinetic_perturb
  PUBLIC :: kinetic_nprint
  PUBLIC :: kinetic_print
  PUBLIC :: kinetic_add_perturb
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE kinetic_open(idx)
    !----------------------------------------------------------------------------
    !
    ! ... open file of Kinetic Energy Density
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
    IF (.NOT. do_kinetic) THEN
      RETURN
    END IF
    !
    iunkinetic = find_free_unit()
    !
    IF (idx >= 0) THEN
      !
      WRITE(str, '(I8)') idx
      filename = TRIM(tmp_dir) // TRIM(prefix) // '.ked.' // TRIM(ADJUSTL(str))
      !
    ELSE
      !
      filename = TRIM(tmp_dir) // TRIM(prefix) // '.ked'
      !
    END IF
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
  SUBROUTINE kinetic_print(with_dtdr, idx)
    !----------------------------------------------------------------------------
    !
    ! ... print data for Kinetic Energy Density
    !
    IMPLICIT NONE
    !
    LOGICAL,           INTENT(IN) :: with_dtdr
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
    IF (.NOT. do_kinetic) THEN
      RETURN
    END IF
    !
    IF (ionode .AND. idx_ < 0) THEN
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Writing Kinetic Energy Density to file.")')
      !
    END IF
    !
    IF (dfftp%nr1  /= dffts%nr1  .OR. dfftp%nr2  /= dffts%nr2  .OR. dfftp%nr3  /= dffts%nr3 .OR. &
        dfftp%nr1x /= dffts%nr1x .OR. dfftp%nr2x /= dffts%nr2x .OR. dfftp%nr3x /= dffts%nr3x) THEN
      !
      CALL errore('kinetic_print', 'Kinetic Energy Density does not support dual FFT-mesh', 1)
      !
    END IF
    !
    IF (nspin /= 1) THEN
      !
      CALL errore('kinetic_print', 'Kinetic Energy Density supports only nspin = 1', 1)
      !
    END IF
    !
    IF (isolve == 4 .AND. .NOT. rmm_conv) THEN
      !
      CALL errore('kinetic_print', &
      'diago_rmm_conv must be .TRUE., when calculate Kinetic Energy Density', 1)
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
      IF (kinetic_dtdr .AND. with_dtdr) THEN
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
      IF (kinetic_dtdr .AND. with_dtdr) THEN
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
  END SUBROUTINE kinetic_print
  !
  !----------------------------------------------------------------------------
  SUBROUTINE kinetic_add_perturb(aux)
    !----------------------------------------------------------------------------
    !
    ! ... add perturbation potential to local potential, in G-space
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: aux(dfftp%nnr)
    !
    INTEGER  :: ia, it
    INTEGER  :: ig
    REAL(DP) :: fac
    REAL(DP) :: ztot
    REAL(DP) :: qtot
    !
    REAL(DP), ALLOCATABLE :: za(:)
    REAL(DP), ALLOCATABLE :: qa(:)
    !
    INTEGER  :: seed(2)
    INTEGER  :: clock
    REAL(DP) :: rand_value
    LOGICAL, SAVE :: rand_init = .FALSE.
    !
    IF (.NOT. do_kinetic) THEN
      RETURN
    END IF
    !
    IF (kinetic_perturb <= 0.0_DP) THEN
      RETURN
    END IF
    !
    ALLOCATE(za(nat))
    ALLOCATE(qa(nat))
    za = 0.0_DP
    qa = 0.0_DP
    !
    IF (ionode) THEN
      !
      IF (.NOT. rand_init) THEN
        rand_init = .TRUE.
        CALL system_clock(count=clock)
        seed(1) = clock
        seed(2) = clock / 2
        CALL random_seed(put=seed)
      END IF
      !
      DO ia = 1, nat
        !
        it = ityp(ia)
        !
        IF (zv(it) > 0.0_DP) THEN
          CALL random_number(rand_value)
          fac = 2.0_DP * (rand_value - 0.5_DP) * kinetic_perturb / zv(it)
          fac = MAX(fac, -1.0_DP)
        ELSE
          fac = 0.0_DP
        END IF
        !
        qa(ia) = fac * zv(it)
        za(ia) = zv(it) + qa(ia)
        !
      END DO
      !
    END IF
    !
    CALL mp_sum(za, intra_image_comm)
    CALL mp_sum(qa, intra_image_comm)
    !
    IF (ionode) THEN
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Perturbed Atomic Charge:")')
    END IF
    !
    ztot = 0.0_DP
    qtot = 0.0_DP
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      !
      ztot = ztot + za(ia)
      qtot = qtot + qa(ia)
      !
      IF (ionode) &
      WRITE(stdout, '(5X,I3,2X,A4,2F10.4)') ia, ADJUSTL(atm(it)) // '    ', za(ia), qa(ia)
      !
      IF (zv(it) > 0.0_DP) THEN
        fac = qa(ia) / zv(it)
      ELSE
        fac = 0.0_DP
      END IF
      !
      DO ig = 1, ngm
        !
        aux(dfftp%nl(ig)) = aux(dfftp%nl(ig)) &
                          + fac * vloc(igtongl(ig), it) &
                          * eigts1(mill(1, ig), ia) &
                          * eigts2(mill(2, ig), ia) &
                          * eigts3(mill(3, ig), ia)
        !
      END DO
      !
    END DO
    !
    IF (ionode) &
    WRITE(stdout, '(10X,A4,2F10.4)') "Sum ", ztot, qtot
    !
    DEALLOCATE(za)
    DEALLOCATE(qa)
    !
  END SUBROUTINE kinetic_add_perturb
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

