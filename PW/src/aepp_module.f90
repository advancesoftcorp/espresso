!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE aepp_module
  !--------------------------------------------------------------------------
  !
  ! ... the module for Adaptive Effective Pseudo-Potential (AEPP),
  ! ... while effective local potential is read from external file.
  !
  USE constants,   ONLY : e2
  USE fft_base,    ONLY : dfftp
  USE io_files,    ONLY : tmp_dir, prefix
  USE io_global,   ONLY : ionode
  USE kinds,       ONLY : DP
  USE mp,          ONLY : mp_sum, mp_barrier
  USE mp_images,   ONLY : intra_image_comm
  USE scatter_mod, ONLY : scatter_grid
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL               :: do_aepp   = .FALSE.
  LOGICAL               :: has_vaepp = .FALSE.
  CHARACTER(LEN=256)    :: filaepp   = ''
  REAL(DP), ALLOCATABLE :: vaepp(:)
  !
  PUBLIC :: do_aepp
  PUBLIC :: filaepp
  PUBLIC :: aepp_initialize
  PUBLIC :: aepp_finalize
  PUBLIC :: aepp_add_vloc
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE aepp_initialize()
    !----------------------------------------------------------------------------
    !
    ! ... read effective local potential from file
    !
    IMPLICIT NONE
    !
    IF (.NOT. do_aepp) RETURN
    !
    IF (has_vaepp) THEN
      !
      CALL aepp_finalize()
      !
    END IF
    !
    IF (LEN(TRIM(filaepp)) == 0) THEN
      !
      filaepp = TRIM(tmp_dir) // TRIM(prefix) // '.aepp'
      !
    END IF
    !
    ALLOCATE(vaepp(dfftp%nnr))
    !
    CALL read_cube_file(filaepp, vaepp)
    !
    has_vaepp = .TRUE.
    !
  END SUBROUTINE aepp_initialize
  !
  !----------------------------------------------------------------------------
  SUBROUTINE aepp_finalize()
    !----------------------------------------------------------------------------
    !
    ! ... release effective local potential
    !
    IMPLICIT NONE
    !
    IF (.NOT. do_aepp) RETURN
    !
    IF (.NOT. has_vaepp) RETURN
    !
    DEALLOCATE(vaepp)
    !
    has_vaepp = .FALSE.
    !
  END SUBROUTINE aepp_finalize
  !
  !----------------------------------------------------------------------------
  SUBROUTINE aepp_add_vloc(vloc)
    !----------------------------------------------------------------------------
    !
    ! ... add effective local potential
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: vloc(dfftp%nnr)
    !
    IF (.NOT. do_aepp) RETURN
    !
    IF (.NOT. has_vaepp) THEN
      !
      CALL aepp_initialize()
      !
    END IF
    !
    vloc = vloc + vaepp
    !
  END SUBROUTINE aepp_add_vloc
  !
  !----------------------------------------------------------------------------
  SUBROUTINE read_cube_file(filename, v)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN)  :: filename
    REAL(DP),         INTENT(OUT) :: v(dfftp%nnr)
    !
    INTEGER  :: iun
    INTEGER  :: ios
    INTEGER  :: nat, ia
    INTEGER  :: nr1, nr2, nr3
    INTEGER  :: ir1, ir2, ir3, ir, jr
    REAL(DP) :: xyz(3)
    !
    REAL(DP), ALLOCATABLE :: vaux(:)
    REAL(DP), ALLOCATABLE :: vcub(:)
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    ios = 0
    !
    ALLOCATE(vaux(dfftp%nr1x * dfftp%nr2x * dfftp%nr3x))
    ALLOCATE(vcub(dfftp%nr1  * dfftp%nr2  * dfftp%nr3 ))
    !
    vaux = 0.0_DP
    !
    IF (ionode) THEN
      !
      iun = find_free_unit()
      !
      OPEN(unit=iun, file=filename, status='old', form='formatted', action='read', iostat=ios)
      !
      IF (ios == 0) THEN ! opened
        !
        READ(iun, '()')
        READ(iun, '()')
        READ(iun, *) nat, xyz
        READ(iun, *) nr1, xyz
        READ(iun, *) nr2, xyz
        READ(iun, *) nr3, xyz
        !
        IF (nr1 /= dfftp%nr1 .OR. nr2 /= dfftp%nr2 .OR. nr3 /= dfftp%nr3) THEN
          !
          ios = 1
          !
          CALL infomsg('aepp_initialize', 'incorrect FFT-mesh at: ' // TRIM(filename))
          !
        END IF
        !
        IF (ios == 0) THEN ! correct mesh
          !
          DO ia = 1, ABS(nat)
            !
            READ(iun, '()')
            !
          END DO
          !
          READ(iun, *) vcub(:)
          !
          DO ir1 = 1, nr1
            !
            DO ir2 = 1, nr2
              !
              DO ir3 = 1, nr3
                !
                ir = ir1 + (ir2 - 1) * dfftp%nr1x + (ir3 - 1) * dfftp%nr1x * dfftp%nr2x
                jr = ir1 + (ir2 - 1) * dfftp%nr1  + (ir3 - 1) * dfftp%nr1  * dfftp%nr2
                !
                vaux(ir) = e2 * vcub(jr) ! Hartree -> Rydberg
                !
              END DO
              !
            END DO
            !
          END DO
          !
        END IF ! correct mesh
        !
        CLOSE(unit=iun)
        !
      END IF ! opened
      !
    END IF
    !
    CALL mp_barrier(intra_image_comm)
    !
    CALL mp_sum(ios, intra_image_comm)
    !
    IF (ios /= 0) THEN
      !
      DEALLOCATE(vaux)
      !
      CALL errore('aepp_initialize', 'cannot open file: ' // TRIM(filename), ios)
      !
      RETURN
      !
    END IF
    !
#if defined(__MPI)
    CALL scatter_grid(dfftp, vaux, v)
#else
    v = vaux
#endif
    !
    DEALLOCATE(vaux)
    DEALLOCATE(vcub)
    !
  END SUBROUTINE read_cube_file
  !
END MODULE aepp_module
