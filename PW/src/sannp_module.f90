!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE sannp_module
  !--------------------------------------------------------------------------
  !
  ! ... the module for Single Atom Neural Network Potential (SANNP)
  !
  USE bp,                ONLY : lelfield
  USE cell_base,         ONLY : at, alat
  USE constants,         ONLY : eps16
  USE control_flags,     ONLY : textfor, llondon, ldftd3, lxdm, ts_vdw
  USE Coul_cut_2D,       ONLY : do_cutoff_2D
  USE ener,              ONLY : etot, demet
  USE esm,               ONLY : do_comp_esm
  USE extfield,          ONLY : tefield, dipfield, gate
  USE fcp_module,        ONLY : lfcp
  USE fft_base,          ONLY : dfftp
  USE force_mod,         ONLY : lforce
  USE funct,             ONLY : dft_is_hybrid, dft_is_meta, dft_is_nonlocc
  USE gcscf_module,      ONLY : lgcscf
  USE io_files,          ONLY : tmp_dir, prefix
  USE io_global,         ONLY : ionode, stdout
  USE ions_base,         ONLY : nat, ityp, atm, tau
  USE kinds,             ONLY : DP
  USE ldaU,              ONLY : lda_plus_u
  USE lsda_mod,          ONLY : nspin
  USE mp,                ONLY : mp_sum, mp_barrier
  USE mp_bands,          ONLY : inter_bgrp_comm
  USE mp_images,         ONLY : intra_image_comm
  USE mp_pools,          ONLY : inter_pool_comm
  USE martyna_tuckerman, ONLY : do_comp_mt
  USE noncollin_module,  ONLY : noncolin
  USE paw_variables,     ONLY : okpaw
  USE qmmm,              ONLY : qmmm_mode
  USE rism_module,       ONLY : lrism
  USE spin_orb,          ONLY : lspinorb
  USE uspp_param,        ONLY : nhm
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL :: do_sannp = .FALSE.
  INTEGER :: iunsannp
  !
  PUBLIC :: do_sannp
  PUBLIC :: sannp_check
  PUBLIC :: sannp_open
  PUBLIC :: sannp_close
  PUBLIC :: sannp_print
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE sannp_check()
    !----------------------------------------------------------------------------
    !
    ! ... check condition for SANNP
    !
    IMPLICIT NONE
    !
    IF (.NOT. do_sannp) THEN
      RETURN
    END IF
    !
    IF (.NOT. lforce) THEN
      CALL errore('sannp_check', 'you have to calculate force for SANNP', 1)
    END IF
    !
    IF (qmmm_mode > -1) THEN
      CALL errore('sannp_check', 'SANNP does not support QM/MM', 1)
    END IF
    !
    ! NB: the next version will support it.
    IF (noncolin) THEN
      CALL errore('sannp_check', 'SANNP does not support non-colinear spin', 1)
    END IF
    !
    ! NB: the next version will support it.
    IF (lspinorb) THEN
      CALL errore('sannp_check', 'SANNP does not support Spin-Orbit spin', 1)
    END IF
    !
  END SUBROUTINE sannp_check
  !
  !----------------------------------------------------------------------------
  LOGICAL FUNCTION sannp_alternative()
    !----------------------------------------------------------------------------
    !
    ! ... check alternative mode for HDNNP
    !
    IMPLICIT NONE
    !
    sannp_alternative = .FALSE.
    !
    IF (dft_is_hybrid()) THEN
      CALL infomsg('sannp_print', 'SANNP does not support Hybrid-GGA')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (dft_is_meta()) THEN
      CALL infomsg('sannp_print', 'SANNP does not support meta-GGA')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (dft_is_nonlocc()) THEN
      CALL infomsg('sannp_print', 'SANNP does not support vdW-DF/rVV10')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (lelfield) THEN
      CALL infomsg('sannp_print', 'SANNP does not support electric field (lelfield)')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (tefield) THEN
      CALL infomsg('sannp_print', 'SANNP does not support electric field (tefield)')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (dipfield) THEN
      CALL infomsg('sannp_print', 'SANNP does not support dipole correction')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (gate) THEN
      CALL infomsg('sannp_print', 'SANNP does not support gate potential')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (textfor) THEN
      CALL infomsg('sannp_print', 'SANNP does not support external force')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (llondon .OR. ldftd3) THEN
      CALL infomsg('sannp_print', 'SANNP does not support DFT-D')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (lxdm) THEN
      CALL infomsg('sannp_print', 'SANNP does not support XDM')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (ts_vdw) THEN
      CALL infomsg('sannp_print', 'SANNP does not support TS-vdW')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (do_comp_mt) THEN
      CALL infomsg('sannp_print', 'SANNP does not support Martyna-Tuckerman')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (do_comp_esm) THEN
      CALL infomsg('sannp_print', 'SANNP does not support ESM method')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (do_cutoff_2D) THEN
      CALL infomsg('sannp_print', 'SANNP does not support 2D-framework')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (lgcscf) THEN
      CALL infomsg('sannp_print', 'SANNP does not support GC-SCF')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (lfcp) THEN
      CALL infomsg('sannp_print', 'SANNP does not support FCP')
      sannp_alternative = .TRUE.
    END IF
    !
    IF (lrism) THEN
      CALL infomsg('sannp_print', 'SANNP does not support 3D-RISM')
      sannp_alternative = .TRUE.
    END IF
    !
    ! NB: the next version will support it.
    IF (okpaw) THEN
      CALL infomsg('sannp_print', 'SANNP does not support PAW')
      sannp_alternative = .TRUE.
    END IF
    !
    ! NB: the next version will support it.
    IF (lda_plus_u) THEN
      CALL infomsg('sannp_print', 'SANNP does not support LDA+U')
      sannp_alternative = .TRUE.
    END IF
    !
  END FUNCTION sannp_alternative
  !
  !----------------------------------------------------------------------------
  SUBROUTINE sannp_open()
    !----------------------------------------------------------------------------
    !
    ! ... open file of SANNP
    !
    IMPLICIT NONE
    !
    INTEGER            :: ios
    CHARACTER(LEN=256) :: filename
    !
    INTEGER, EXTERNAL  :: find_free_unit
    !
    IF (.NOT. do_sannp) THEN
      RETURN
    END IF
    !
    iunsannp = find_free_unit()
    !
    filename = TRIM(tmp_dir) // TRIM(prefix) // '.sannp'
    !
    IF (ionode) THEN
      !
      OPEN(unit=iunsannp, file=TRIM(filename), &
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
    CALL errore('sannp_open', 'cannot open file: ' // TRIM(filename), ios)
    !
  END SUBROUTINE sannp_open
  !
  !----------------------------------------------------------------------------
  SUBROUTINE sannp_close()
    !----------------------------------------------------------------------------
    !
    ! ... close file of SANNP
    !
    IMPLICIT NONE
    !
    LOGICAL :: opnd
    !
    IF (.NOT. do_sannp) THEN
      RETURN
    END IF
    !
    IF (ionode) THEN
      !
      INQUIRE(unit=iunsannp, opened=opnd)
      !
      IF (opnd) CLOSE(unit=iunsannp)
      !
    END IF
    !
  END SUBROUTINE sannp_close
  !
  !----------------------------------------------------------------------------
  SUBROUTINE sannp_print(force)
    !----------------------------------------------------------------------------
    !
    ! ... print data for SANNP
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: force(3, nat)
    !
    INTEGER  :: i
    INTEGER  :: ia, ia_start, ia_end
    INTEGER  :: ja, ja_start, ja_end
    LOGICAL  :: alter
    REAL(DP) :: x, y, z
    REAL(DP) :: fx,  fy,  fz
    REAL(DP) :: fx0, fy0, fz0
    REAL(DP) :: e, etot0, e0
    REAL(DP) :: q, qsum
    !
    REAL(DP), ALLOCATABLE :: energy(:)
    REAL(DP), ALLOCATABLE :: charge(:)
    REAL(DP), ALLOCATABLE :: energy_q(:)
    REAL(DP), ALLOCATABLE :: force_q(:, :)
    REAL(DP), ALLOCATABLE :: rhokin(:)
    REAL(DP), ALLOCATABLE :: rhobec(:,:,:)
    !
    IF (.NOT. do_sannp) THEN
      RETURN
    END IF
    !
    IF (ionode) THEN
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Atomic energies & charges (w/ Hirshfeld method)")')
      !
    END IF
    !
    alter = sannp_alternative()
    !
    IF (alter) THEN
      !
      CALL infomsg('sannp_print', &
      & 'WARNING: atomic energies are NOT correct, you cannot perform SANNP but HDNNP.')
      !
    END IF
    !
    CALL divide(inter_pool_comm, nat, ia_start, ia_end)
    !
    ja = ia_end - ia_start + 1
    !
    IF (ja > 0) THEN
      !
      CALL divide(inter_bgrp_comm, ja, ja_start, ja_end)
      !
      ia_start = ia_start + ja_start - 1
      ia_end   = ia_start + ja_end   - 1
      !
    END IF
    !
    ! ... allocate memory
    !
    ALLOCATE(energy(nat))
    ALLOCATE(charge(nat))
    ALLOCATE(energy_q(nat))
    ALLOCATE(force_q(3, nat))
    ALLOCATE(rhokin(dfftp%nnr))
    ALLOCATE(rhobec(nhm * (nhm + 1) / 2, nat, nspin))
    !
    ! ... calculate atomic energies and charges
    !
    DO ia = 1, nat
      !
      energy(ia) = 0.0_DP
      !
      charge(ia) = 0.0_DP
      !
    END DO
    !
    CALL sannp_sum_band(rhokin, rhobec)
    !
    CALL sannp_energy_lc(energy, charge, rhokin, ia_start, ia_end)
    !
    CALL sannp_energy_us(energy, rhobec, ia_start, ia_end)
    !
    CALL sannp_energy_ew(energy, ia_start, ia_end)
    !
    IF (okpaw)      CALL sannp_energy_paw(energy, ia_start, ia_end)
    !
    IF (lda_plus_u) CALL sannp_energy_hub(energy, ia_start, ia_end)
    !
    DO ia = ia_start, ia_end
      !
      energy(ia) = energy(ia) + demet / DBLE(nat)
      !
    END DO
    !
    qsum = 0.0_DP
    !
    DO ia = ia_start, ia_end
      !
      qsum = qsum + charge(ia)
      !
    END DO
    !
    CALL mp_sum(qsum, inter_bgrp_comm)
    CALL mp_sum(qsum, inter_pool_comm)
    !
    IF (ABS(qsum) >= eps16) THEN
      !
      DO ia = ia_start, ia_end
        !
        charge(ia) = charge(ia) - qsum / DBLE(nat)
        !
      END DO
      !
    END IF
    !
    CALL mp_sum(energy, inter_bgrp_comm)
    CALL mp_sum(energy, inter_pool_comm)
    !
    CALL mp_sum(charge, inter_bgrp_comm)
    CALL mp_sum(charge, inter_pool_comm)
    !
    IF (ionode) THEN
      !
      etot0 = 0.0_DP
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,A)') '  # atom     E (Ry)         Q (e)'
      !
      DO ia = 1, nat
        !
        etot0 = etot0 + energy(ia)
        !
        WRITE(stdout, '(5X,I3,2X,A4,0PE14.5,F10.4)') ia, &
        & ADJUSTL(atm(ityp(ia))) // '    ', energy(ia), charge(ia)
        !
      END DO
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"total energy of SCF       =",F17.8," Ry")') etot
      WRITE(stdout, '(5X,"sum of atomic energies    =",F17.8," Ry")') etot0
      !
    END IF
    !
    ! ... Ewald energy and force from atomic charges
    !
    DO ia = 1, nat
      !
      energy_q(ia) = 0.0_DP
      !
      force_q(1, ia) = 0.0_DP
      force_q(2, ia) = 0.0_DP
      force_q(3, ia) = 0.0_DP
      !
    END DO
    !
    CALL sannp_charge_ew(energy_q, force_q, charge, ia_start, ia_end)
    !
    CALL mp_sum(energy_q, inter_bgrp_comm)
    CALL mp_sum(energy_q, inter_pool_comm)
    !
    CALL mp_sum(force_q, inter_bgrp_comm)
    CALL mp_sum(force_q, inter_pool_comm)
    !
    ! ... print data
    !
    IF (ionode) THEN
      !
      IF (alter) THEN
        i = 0
      ELSE
        i = 1
      END IF
      !
      WRITE(iunsannp, "(2I8,4X,E20.12)") nat, i, etot
      !
      DO i = 1, 3
        !
        x = alat * at(1, i)
        y = alat * at(2, i)
        z = alat * at(3, i)
        !
        WRITE(iunsannp, "(3E20.12)") x, y, z
        !
      END DO
      !
      DO ia = 1, nat
        !
        x = alat * tau(1, ia)
        y = alat * tau(2, ia)
        z = alat * tau(3, ia)
        !
        fx = force(1, ia)
        fy = force(2, ia)
        fz = force(3, ia)
        !
        e = energy(ia)
        q = charge(ia)
        !
        e0 = energy_q(ia)
        !
        fx0 = force_q(1, ia)
        fy0 = force_q(2, ia)
        fz0 = force_q(3, ia)
        !
        WRITE(iunsannp, "(A5,12E20.12)") &
        & ADJUSTL(atm(ityp(ia))), x, y, z, e, fx, fy, fz, q, e0, fx0, fy0, fz0
        !
      END DO
      !
      FLUSH(iunsannp)
      !
    END IF
    !
    CALL mp_barrier(intra_image_comm)
    !
    ! ... deallocate memory
    !
    DEALLOCATE(energy)
    DEALLOCATE(charge)
    DEALLOCATE(energy_q)
    DEALLOCATE(force_q)
    DEALLOCATE(rhokin)
    DEALLOCATE(rhobec)
    !
  END SUBROUTINE sannp_print
  !
END MODULE sannp_module
