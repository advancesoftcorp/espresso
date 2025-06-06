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
  USE mp,               ONLY : mp_sum, mp_barrier
  USE mp_images,        ONLY : intra_image_comm
  USE noncollin_module, ONLY : noncolin
  USE spin_orb,         ONLY : lspinorb
  USE uspp,             ONLY : dvan, okvan
  USE uspp_param,       ONLY : upf, nh, nhm
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL :: do_occmat = .FALSE.
  LOGICAL :: fhi98_mat = .FALSE.
  INTEGER :: iunoccmat
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
    INTEGER  :: ih, jh
    INTEGER  :: iorb, jorb
    REAL(DP) :: occmat(4, 4)
    !
    INTEGER,  ALLOCATABLE :: indx_h(:,:)
    REAL(DP), ALLOCATABLE :: becsum(:,:,:) ! \sum_i f(i) <psi(i)|beta_l><beta_m|psi(i)>
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
    IF (lspinorb) THEN
      !
      CALL errore('occmat_print', 'Occupation Matrix does not support lspinorb', 1)
      !
    END IF
    !
    IF (fhi98_mat .AND. okvan) THEN
      !
      CALL errore('occmat_print', 'Occupation Matrix of FHI98-type supports only NCPP, not USPP/PAW.', 1)
      !
    END IF
    !
    ! ... allocate memory
    !
    ALLOCATE(indx_h(4, ntyp))
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
      IF (fhi98_mat) THEN
        CALL fhi98_mat_print(becsum)
      ELSE
        CALL normalmat_print(becsum)
      END IF
      !
    END IF
    !
    CALL occmat_close()
    !
    CALL mp_barrier(intra_image_comm)
    !
    ! ... deallocate memory
    !
    DEALLOCATE(indx_h)
    DEALLOCATE(becsum)
    !
  END SUBROUTINE occmat_print
  !
  !----------------------------------------------------------------------------
  SUBROUTINE normalmat_print(becsum)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: becsum(nhm, nhm, nat)
    !
    INTEGER :: ia, it
    INTEGER :: ih, jh, nht
    !
    REAL(DP), ALLOCATABLE :: occmat(:,:)
    REAL(DP), ALLOCATABLE :: tmpmat(:,:)
    REAL(DP), ALLOCATABLE :: rotylm(:,:,:)
    !
    ! ... allocate memory
    !
    ALLOCATE(occmat(nhm, nhm))
    ALLOCATE(tmpmat(nhm, nhm))
    ALLOCATE(rotylm(nhm, nhm, ntyp))
    !
    ! ... print data
    !
    WRITE(iunoccmat, '("#Occupation Matrix")')
    !
    DO it = 1, ntyp
      !
      CALL rotation_ylm(it, rotylm(:, :, it))
      !
    END DO
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      !
      nht = nh(it)
      !
      occmat(1:nht, 1:nht) = becsum(1:nht, 1:nht, ia)
      !
      CALL DGEMM("N", "N", nht, nht, nht, 1.0_DP, &
                rotylm, nhm, occmat, nhm, 0.0_DP, tmpmat, nhm)
      !
      CALL DGEMM("N", "T", nht, nht, nht, 1.0_DP, &
                tmpmat, nhm, rotylm, nhm, 0.0_DP, occmat, nhm)
      !
      WRITE(iunoccmat, '("   iatom:", I5, "  nbasis:", I5)') ia, nht
      !
      DO ih = 1, nht
        !
        WRITE(iunoccmat, '(64E25.16)') (occmat(ih, jh), jh = 1, nht)
        !
      END DO
      !
    END DO
    !
    ! ... deallocate memory
    !
    DEALLOCATE(occmat)
    DEALLOCATE(tmpmat)
    DEALLOCATE(rotylm)
    !
  END SUBROUTINE normalmat_print
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fhi98_mat_print(becsum)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: becsum(nhm, nhm, nat)
    !
    INTEGER  :: ia, it
    INTEGER  :: ih, jh
    INTEGER  :: iorb, jorb
    REAL(DP) :: occmat(4, 4)
    !
    INTEGER,  ALLOCATABLE :: indx_h(:,:)
    !
    ! ... allocate memory
    !
    ALLOCATE(indx_h(4, ntyp))
    !
    ! ... print data
    !
    WRITE(iunoccmat, '("#Occupation Matrix (only s+p orbital)")')
    !
    DO it = 1, ntyp
      !
      CALL index_of_beta(it, indx_h(:, it))
      !
    END DO
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      !
      ! QE's Ylm [z, -x, -y] --> Standard Ylm [x, y, z]
      DO iorb = 1, 4
        !
        ih = indx_h(iorb, it)
        !
        DO jorb = 1, 4
          !
          jh = indx_h(jorb, it)
          !
          IF (ih > 0 .AND. jh > 0) THEN
            occmat(iorb, jorb) = becsum(ih, jh, ia) * dvan(ih, ih, it) * dvan(jh, jh, it)
          ELSE
            occmat(iorb, jorb) = 0.0_DP
          END IF
          !
        END DO
        !
      END DO
      !
      occmat(2, :) = -1.0_DP * occmat(2, :) ! -x -> x
      occmat(3, :) = -1.0_DP * occmat(3, :) ! -y -> y
      occmat(:, 2) = -1.0_DP * occmat(:, 2) ! -x -> x
      occmat(:, 3) = -1.0_DP * occmat(:, 3) ! -y -> y
      !
      WRITE(iunoccmat, '("   iatom:", I5, "  nbasis:", I5)') ia, 4
      !
      DO iorb = 1, 4
        !
        WRITE(iunoccmat, '(4E25.16)') (occmat(iorb, jorb), jorb = 1, 4)
        !
      END DO
      !
    END DO
    !
    ! ... deallocate memory
    !
    DEALLOCATE(indx_h)
    !
  END SUBROUTINE fhi98_mat_print
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rotation_ylm(it, rot)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: it
    INTEGER, INTENT(OUT) :: rot(nhm, nhm)
    !
    INTEGER :: ib
    INTEGER :: l
    INTEGER :: nht
    INTEGER :: i1, i2, i3, i4, i5,
    !
    REAL(DP), PARAMETER :: ROOT3 = DSQRT(3.0_DP)
    !
    rot(:, :) = 0.0_DP
    !
    nht = 0
    !
    DO ib = 1, upf(it)%nbeta
      !
      l = upf(it)%lll(ib)
      !
      IF (l == 0) THEN
        !
        i1 = nht + 1
        !
        rot(i1, i1) = +1.0_DP
        !
      ELSE IF (l == 1) THEN
        !
        i1 = nht + 1
        i2 = nht + 2
        i3 = nht + 3
        !
        ! { z, -x, -y } -> { x, y, z }
        !
        rot(i1, i2) = -1.0_DP ! -> x
        rot(i2, i3) = -1.0_DP ! -> y
        rot(i3, i1) =  1.0_DP ! -> z
        !
      ELSE IF (l == 2) THEN
        !
        i1 = nht + 1
        i2 = nht + 2
        i3 = nht + 3
        i4 = nht + 4
        i5 = nht + 5
        !
        ! { (3z2-r2)/(2*sqrt(3)), -xz, -yz, (x2-y2)/2, xy }
        ! -> { xz, xy, (3y2-r2)/(2*sqrt(3)), yz, (z2-x2)/2 }
        !
        rot(i1, i2) = -1.0_DP ! -> xz
        rot(i2, i5) =  1.0_DP ! -> xy
        rot(i3, i1) = -0.5_DP ! -> (3y2-r2)/(2*sqrt(3)) = (-1/2)*(3z2-r2)/(2*sqrt(3)) - (sqrt(3)/2)*(x2-y2)/2
        rot(i3, i4) = -0.5_DP * ROOT3
        rot(i4, i3) = -1.0_DP ! -> yz
        rot(i5, i1) =  0.5_DP * ROOT3
        rot(i5, i4) = -0.5_DP ! -> (z2-x2)/2 = (sqrt(3)/2)*(3z2-r2)/(2*sqrt(3)) - (1/2)*(x2-y2)/2
        !
      ELSE
        !
        CALL errore('occmat_print', 'Occupation Matrix does not support the case of l >= 3', 1)
        !
      END IF
      !
      nht = nht + 2 * l + 1
      !
    END DO
    !
  END SUBROUTINE rotation_ylm
  !
  !----------------------------------------------------------------------------
  SUBROUTINE index_of_beta(it, indx_h)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)  :: it
    INTEGER, INTENT(OUT) :: indx_h(4)
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
    indx_h(1:4) = 0
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
        IF (indx_h(i_s) == 0) THEN
          !
          indx_h(i_s) = nht + 1 ! s
          !
        END IF
        !
      ELSE IF (l == 1) THEN
        !
        IF (indx_h(i_pz) == 0) THEN
          !
          indx_h(i_pz) = nht + 1 ! pz
          indx_h(i_px) = nht + 2 ! px
          indx_h(i_py) = nht + 3 ! py
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
