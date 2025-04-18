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
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_sum
  USE noncollin_module, ONLY : noncolin
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
    ! ... calculate Occupation Matrix
    !
    CALL occmat_sum_band()
    !
    ! TODO
    ! TODO
    ! TODO
    !
  END SUBROUTINE occmat_print
  !
END MODULE occmat_module
