!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE atomnl_module
  !--------------------------------------------------------------------------
  !
  ! ... the module for Atomic Non-Local Energy
  !
  USE atom,             ONLY : rgrid
  USE cell_base,        ONLY : at, alat
  USE constants,        ONLY : e2, eps16
  USE control_flags,    ONLY : isolve, rmm_conv
  USE fft_base,         ONLY : dffts, dfftp
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : ionode, stdout
  USE ions_base,        ONLY : nat, ityp, atm, tau
  USE lsda_mod,         ONLY : nspin
  USE mp,               ONLY : mp_sum, mp_barrier
  USE mp_images,        ONLY : intra_image_comm
  USE noncollin_module, ONLY : noncolin
  USE scatter_mod,      ONLY : gather_grid
  USE spin_orb,         ONLY : lspinorb
  USE uspp_param,       ONLY : upf, nhm
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL  :: do_atomnl     = .FALSE.
  INTEGER  :: atomnl_nprint = 0
  INTEGER  :: iunatomnl
  !
  PUBLIC :: do_atomnl
  PUBLIC :: atomnl_nprint
  PUBLIC :: atomnl_print
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE atomnl_open(idx)
    !----------------------------------------------------------------------------
    !
    ! ... open file of Atomic Non-Local Energy
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
    IF (.NOT. do_atomnl) THEN
      RETURN
    END IF
    !
    iunatomnl = find_free_unit()
    !
    IF (idx >= 0) THEN
      !
      WRITE(str, '(I8)') idx
      filename = TRIM(tmp_dir) // TRIM(prefix) // '.anl.' // TRIM(ADJUSTL(str))
      !
    ELSE
      !
      filename = TRIM(tmp_dir) // TRIM(prefix) // '.anl'
      !
    END IF
    !
    IF (ionode) THEN
      !
      OPEN(unit=iunatomnl, file=TRIM(filename), &
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
    CALL errore('atomnl_open', 'cannot open file: ' // TRIM(filename), ios)
    !
  END SUBROUTINE atomnl_open
  !
  !----------------------------------------------------------------------------
  SUBROUTINE atomnl_close()
    !----------------------------------------------------------------------------
    !
    ! ... close file of Atomic Non-Local Energy
    !
    IMPLICIT NONE
    !
    LOGICAL :: opnd
    !
    IF (.NOT. do_atomnl) THEN
      RETURN
    END IF
    !
    IF (ionode) THEN
      !
      INQUIRE(unit=iunatomnl, opened=opnd)
      !
      IF (opnd) CLOSE(unit=iunatomnl)
      !
    END IF
    !
  END SUBROUTINE atomnl_close
  !
  !----------------------------------------------------------------------------
  SUBROUTINE atomnl_check()
    !----------------------------------------------------------------------------
    !
    ! ... check condition for Atomic Non-Local Energy
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER :: eps = eps16
    !
    IF (.NOT. do_atomnl) THEN
      RETURN
    END IF
    !
    IF (noncolin) THEN
      CALL errore('atomnl_print', 'Atom-NL does not support non-colinear spin', 1)
    END IF
    !
    IF (lspinorb) THEN
      CALL errore('atomnl_print', 'Atom-NL does not support Spin-Orbit spin', 1)
    END IF
    !
    IF (ABS(at(1, 1) - 1.0_DP) > eps .OR. ABS(at(1, 2)) > eps .OR. ABS(at(1, 3)) > eps .OR. &
      & ABS(at(2, 2) - 1.0_DP) > eps .OR. ABS(at(2, 1)) > eps .OR. ABS(at(2, 3)) > eps .OR. &
      & ABS(at(3, 3) - 1.0_DP) > eps .OR. ABS(at(3, 1)) > eps .OR. ABS(at(3, 2)) > eps) THEN
      CALL errore('atomnl_print', 'Atom-NL works only for cubic cell', 1)
    END IF
    !
    IF (dffts%nr1 /= dffts%nr2 .OR. dffts%nr1 /= dffts%nr3 .OR. dffts%nr2 /= dffts%nr3) THEN
      CALL errore('atomnl_print', 'Atom-NL works only for isotropic FFT-mesh', 1)
    END IF
    !
    IF (nat /= 1) THEN
      CALL errore('atomnl_print', 'Atom-NL works only for single atomic system', 1)
    END IF
    !
    IF (ABS(tau(1, 1)) > eps .OR. ABS(tau(2, 1)) > eps .OR. ABS(tau(3, 1)) > eps) THEN
      CALL errore('atomnl_print', 'You have to put the atom on the position (0,0,0) for Atom-NL', 1)
    END IF
    !
    IF (isolve == 4 .AND. .NOT. rmm_conv) THEN
      CALL errore('atomnl_print', 'diago_rmm_conv must be .TRUE., when calculate Atom-NL', 1)
    END IF
    !
  END SUBROUTINE atomnl_check
  !
  !----------------------------------------------------------------------------
  SUBROUTINE atomnl_print(idx)
    !----------------------------------------------------------------------------
    !
    ! ... print data for Atomic Non-Local Energy
    !
    IMPLICIT NONE
    !
    INTEGER, OPTIONAL, INTENT(IN) :: idx
    !
    INTEGER  :: idx_
    INTEGER  :: it
    INTEGER  :: ir1, ir2, ir3
    INTEGER  :: nr1, nr2, nr3
    INTEGER  :: nr1x, nr2x, nr3x
    INTEGER  :: nfft
    !
    REAL(DP), ALLOCATABLE :: energy(:)
    REAL(DP), ALLOCATABLE :: rhokin(:)  ! this is just a dummy
    REAL(DP), ALLOCATABLE :: rhobec(:,:,:)
    !
    REAL(DP), ALLOCATABLE :: rhor(:)
    REAL(DP), ALLOCATABLE :: tauG(:)    ! this is just a dummy
    REAL(DP), ALLOCATABLE :: tauL(:)    ! this is just a dummy
    REAL(DP), ALLOCATABLE :: dtdr(:)    ! this is just a dummy
    REAL(DP), ALLOCATABLE :: rhor_g(:)
    !
    IF (PRESENT(idx)) THEN
      idx_ = idx
    ELSE
      idx_ = -1
    END IF
    !
    IF (.NOT. do_atomnl) THEN
      RETURN
    END IF
    !
    IF (ionode .AND. idx_ < 0) THEN
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Writing Atomic Non-Local Energy to file.")')
      !
    END IF
    !
    CALL atomnl_check()
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
    ALLOCATE(energy(nat))
    ALLOCATE(rhokin(dfftp%nnr))
    ALLOCATE(rhobec(nhm * (nhm + 1) / 2, nat, nspin))
    !
    ALLOCATE(rhor(dffts%nnr))
    ALLOCATE(tauG(dffts%nnr))
    ALLOCATE(tauL(dffts%nnr))
    ALLOCATE(dtdr(dffts%nnr))
    ALLOCATE(rhor_g(nfft))
    !
    ! ... calculate Non-Local Energy
    !
    energy = 0.0_DP
    !
    CALL sannp_sum_band(rhokin, rhobec)
    !
    CALL sannp_energy_us(energy, rhobec, 1, 1)
    !
    ! ... calculate Charge density from Wave function
    !
    CALL kinetic_sum_band(rhor, tauG, tauL, dtdr, .TRUE.)
    !
#if defined(__MPI)
    rhor_g = 0.0_DP
    CALL gather_grid(dffts, rhor, rhor_g)
#else
    rhor_g = rhor
#endif
    !
    ! ... print data, in Hartree/Bohr unit
    !
    CALL atomnl_open(idx_)
    !
    IF (ionode) THEN
      !
      WRITE(iunatomnl, '("#Mesh")')
      WRITE(iunatomnl, "(I8)") nr1
      !
      WRITE(iunatomnl, '("#Lattice")')
      WRITE(iunatomnl, '(E25.16)') alat
      !
      it = ityp(1) ! using only the first atom
      !
      WRITE(iunatomnl, '("#Element")')
      WRITE(iunatomnl, '(A5)') ADJUSTL(atm(it))
      !
      WRITE(iunatomnl, '("#Cutoff Radius")')
      WRITE(iunatomnl, '(E25.16)') rgrid(it)%r(upf(it)%kkbeta)
      !
      WRITE(iunatomnl, '("#Non-Local Energy")')
      WRITE(iunatomnl, '(E25.16)') energy(1) / e2
      !
      WRITE(iunatomnl, '("#Charge Density")')
      !
      DO ir1 = 1, nr1
        !
        DO ir2 = 1, nr2
          !
          WRITE(iunatomnl,'(6E25.16)') &
          & (rhor_g(ir1 + (ir2 - 1) * nr1x + (ir3 - 1) * nr1x * nr2x), ir3 = 1, nr3)
          !
        END DO
        !
      END DO
      !
    END IF
    !
    CALL atomnl_close()
    !
    CALL mp_barrier(intra_image_comm)
    !
    ! ... deallocate memory
    !
    DEALLOCATE(energy)
    DEALLOCATE(rhokin)
    DEALLOCATE(rhobec)
    !
    DEALLOCATE(rhor)
    DEALLOCATE(tauG)
    DEALLOCATE(tauL)
    DEALLOCATE(dtdr)
    DEALLOCATE(rhor_g)
    !
  END SUBROUTINE atomnl_print
  !
END MODULE atomnl_module

