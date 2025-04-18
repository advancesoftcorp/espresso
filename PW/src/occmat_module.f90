!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE occmat_module
  !--------------------------------------------------------------------------
  !
  ! ... the module for Occupation Matrix
  !
  USE cell_base,        ONLY : at, alat
  USE io_files,         ONLY : tmp_dir, prefix
  USE io_global,        ONLY : ionode, stdout
  USE ions_base,        ONLY : nat, atm, ntyp => nsp, ityp, tau
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_sum
  USE mp_images,        ONLY : intra_image_comm
  USE noncollin_module, ONLY : noncolin
  USE uspp_param,       ONLY : upf, nhm
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL  :: do_occmat = .FALSE.
  INTEGER  :: iunoccmat
  !
  PUBLIC :: do_occmat
  PUBLIC :: occmat_print
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE occmat_open()
    !----------------------------------------------------------------------------
    !
    ! ... open file of Occupation Matrix
    !
    IMPLICIT NONE
    !
    INTEGER            :: ios
    CHARACTER(LEN=8)   :: str
    CHARACTER(LEN=256) :: filename
    !
    INTEGER, EXTERNAL  :: find_free_unit
    !
    IF (.NOT. do_occmat) THEN
      RETURN
    END IF
    !
    iunoccmat = find_free_unit()
    !
    filename = TRIM(tmp_dir) // TRIM(prefix) // '.occ'
    !
    IF (ionode) THEN
      !
      OPEN(unit=iunoccmat, file=TRIM(filename), &
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
    CALL errore('occmat_open', 'cannot open file: ' // TRIM(filename), ios)
    !
  END SUBROUTINE occmat_open
  !
  !----------------------------------------------------------------------------
  SUBROUTINE occmat_close()
    !----------------------------------------------------------------------------
    !
    ! ... close file of Occupation Matrix
    !
    IMPLICIT NONE
    !
    LOGICAL :: opnd
    !
    IF (.NOT. do_occmat) THEN
      RETURN
    END IF
    !
    IF (ionode) THEN
      !
      INQUIRE(unit=iunoccmat, opened=opnd)
      !
      IF (opnd) CLOSE(unit=iunoccmat)
      !
    END IF
    !
  END SUBROUTINE occmat_close
  !
  !----------------------------------------------------------------------------
  SUBROUTINE occmat_print()
    !----------------------------------------------------------------------------
    !
    ! ... print data for Occupation Matrix
    !
    IMPLICIT NONE
    !
    INTEGER  :: ia, it
    INTEGER  :: ib, jb
    INTEGER  :: iorb, jorb
    REAL(DP) :: becsum_t
    !
    INTEGER,  ALLOCATABLE :: i_beta(:,:)
    REAL(DP), ALLOCATABLE :: becsum(:,:) ! \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
    !
    IF (.NOT. do_occmat) THEN
      RETURN
    END IF
    !
    IF (ionode) THEN
      !
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Writing Occupation Matrix to file.")')
      !
    END IF
    !
    IF (noncolin) THEN
      !
      CALL errore('occmat_print', 'Occupation Matrix does not support noncolin', 1)
      !
    END IF
    !
    ! ... allocate memory
    !
    ALLOCATE(i_beta(4, ntyp))
    ALLOCATE(becsum(nhm, nhm, nat)) ! w/o spin
    !
    ! ... calculate Occupation Matrix
    !
    CALL occmat_sum_band(becsum)
    !
    ! ... print data, in Hartree/Bohr unit
    !
    CALL occmat_open()
    !
    IF (ionode) THEN
      !
      ! ... Lattice
      WRITE(iunoccmat, '("#Lattice")')
      WRITE(iunoccmat, '(3E25.16)') alat * at(1, 1), alat * at(2, 1), alat * at(3, 1)
      WRITE(iunoccmat, '(3E25.16)') alat * at(1, 2), alat * at(2, 2), alat * at(3, 2)
      WRITE(iunoccmat, '(3E25.16)') alat * at(1, 3), alat * at(2, 3), alat * at(3, 3)
      !
      ! ... Atoms
      WRITE(iunoccmat, '("#Number of Atoms")')
      WRITE(iunoccmat, '(I5)') nat
      WRITE(iunoccmat, '("#Atoms")')
      !
      DO ia = 1, nat
        !
        it = ityp(ia)
        !
        WRITE(iunoccmat, '(I5,A6,3E25.16)') ia, atm(it), &
        alat * tau(1, ia), alat * tau(2, ia), alat * tau(3, ia)
        !
      END DO
      !
      ! ... Occupation Matrix
      WRITE(iunoccmat, '("#Occupation Matrix (only s+p orbital)")')
      !
      DO it = 1, ntyp
        !
        CALL index_of_beta(it, i_beta(:, it))
        !
      END DO
      !
      DO ia = 1, nat
        !
        it = ityp(ia)
        !
        WRITE(iunoccmat, '(I5)') ia
        !
        DO iorb = 1, 4
          !
          ib = i_beta(iorb, it)
          !
          DO jorb = 1, 4
            !
            jb = i_beta(jorb, it)
            !
            IF (ib > 0 .AND. jb > 0) THEN
              becsum_t = becsum(ib, jb, na)
            ELSE
              becsum_t = 0.0_DP
            END IF
            !
            WRITE(iunoccmat, '(E25.16)', advance='no') becsum_t
            !
          END DO
          !
          WRITE(iunoccmat, '()')
          !
        END DO
        !
      END DO
      !
    END IF
    !
    CALL occmat_close()
    !
    CALL mp_barrier(intra_image_comm)
    !
    ! ... deallocate memory
    !
    DEALLOCATE(iorb_s)
    DEALLOCATE(iorb_p)
    DEALLOCATE(becsum)
    !
  END SUBROUTINE occmat_print
  !
  !----------------------------------------------------------------------------
  SUBROUTINE index_of_beta(it, i_beta)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: it
    INTEGER, INTENT(OUT) :: i_beta(4)
    !
    INTEGER :: ib
    INTEGER :: l
    INTEGER :: nht
    !
    INTEGER, PARAMETER :: i_s  = 1
    INTEGER, PARAMETER :: i_px = 2
    INTEGER, PARAMETER :: i_py = 3
    INTEGER, PARAMETER :: i_pz = 4
    !
    i_beta(1:4) = 0
    !
    IF (upf(it)%tcoulombp) RETURN
    !
    nht = 0
    !
    DO ib = 1, upf(it)%nbeta
      !
      l = upf(it)%lll(ib)
      !
      IF (l == 0) THEN
        !
        IF (i_beta(i_s) == 0) THEN
          !
          i_beta(i_s) = nht + 1 ! s
          !
        END IF
        !
      ELSE IF (l == 1) THEN
        !
        IF (i_beta(i_pz) == 0) THEN
          !
          i_beta(i_pz) = nht + 1 ! pz
          i_beta(i_px) = nht + 2 ! px
          i_beta(i_py) = nht + 3 ! py
          !
        END IF
        !
      END IF
      !
      nht = nht + 2 * l + 1
      !
    END DO
    !
  END SUBROUTINE index_of_beta
  !
END MODULE occmat_module
