!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE zmp_module
  !--------------------------------------------------------------------------
  !
  ! ... the module for Zhao-Morrison-Parr (ZMP) method,
  ! ... where the spatial regions are limited inside the pseudopotential radii.
  ! ... The coulombic interaction is screened as Yukawa or Erfc/r.
  !
  USE constants,   ONLY : e2
  USE fft_base,    ONLY : dfftp
  USE io_files,    ONLY : tmp_dir, prefix
  USE io_global,   ONLY : ionode, ionode_id
  USE kinds,       ONLY : DP
  USE mp,          ONLY : mp_sum, mp_bcast, mp_barrier
  USE mp_images,   ONLY : intra_image_comm
  USE scatter_mod, ONLY : scatter_grid
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL               :: do_zmp  = .FALSE.
  CHARACTER(LEN=256)    :: filzmp  = ''
  INTEGER               :: n_group = 0
  REAL(DP)              :: lambda
  REAL(DP)              :: omega
  LOGICAL               :: yukawa
  REAL(DP), ALLOCATABLE :: wei(:, :) ! (dfftp%nnr, n_group)
  REAL(DP), ALLOCATABLE :: rho(:, :) ! (dfftp%nnr, n_group)
  !
  PUBLIC :: do_zmp
  PUBLIC :: filzmp
  PUBLIC :: zmp_initialize
  PUBLIC :: zmp_finalize
  PUBLIC :: add_vzmp
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE zmp_initialize()
    !----------------------------------------------------------------------------
    !
    ! ... read data of ZMP-potential from file
    !
    IMPLICIT NONE
    !
    IF (.NOT. do_zmp) RETURN
    !
    IF (n_group > 0) THEN
      !
      CALL zmp_finalize()
      !
    END IF
    !
    IF (LEN(TRIM(filzmp)) == 0) THEN
      !
      filzmp = TRIM(tmp_dir) // TRIM(prefix) // '.zmp'
      !
    END IF
    !
    CALL read_zmp_file(filzmp)
    !
  END SUBROUTINE zmp_initialize
  !
  !----------------------------------------------------------------------------
  SUBROUTINE zmp_finalize()
    !----------------------------------------------------------------------------
    !
    ! ... release effective local potential
    !
    IMPLICIT NONE
    !
    IF (.NOT. do_zmp) RETURN
    !
    IF (n_group < 1)  RETURN
    !
    n_group = 0
    !
    DEALLOCATE(wei)
    !
    DEALLOCATE(rho)
    !
  END SUBROUTINE zmp_finalize
  !
  !----------------------------------------------------------------------------
  SUBROUTINE aepp_add_vloc(vloc, minus)
    !----------------------------------------------------------------------------
    !
    ! ... add effective local potential
    !
    IMPLICIT NONE
    !
    REAL(DP),          INTENT(INOUT) :: vloc(dfftp%nnr)
    LOGICAL, OPTIONAL, INTENT(IN)    :: minus
    !
    REAL(DP) :: fac
    !
    IF (.NOT. do_aepp) RETURN
    !
    IF (.NOT. has_vaepp) THEN
      !
      CALL aepp_initialize()
      !
    END IF
    !
    fac = 1.0_DP
    !
    IF (PRESENT(minus)) THEN
      IF (minus) fac = -1.0_DP
    END IF
    !
    vloc = vloc + fac * vaepp
    !
  END SUBROUTINE aepp_add_vloc
  !
  !----------------------------------------------------------------------------
  SUBROUTINE read_zmp_file(filename)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: filename
    !
    INTEGER  :: iun
    INTEGER  :: ios
    INTEGER  :: ir1, ir2, ir3
    INTEGER  :: ir, jr
    INTEGER  :: i_group
    !
    CHARACTER(LEN=256) :: line
    !
    INTEGER, EXTERNAL :: find_free_unit
    !
    REAL(DP), ALLOCATABLE :: wei_t(:, :)
    REAL(DP), ALLOCATABLE :: wei_x(:, :)
    REAL(DP), ALLOCATABLE :: rho_t(:, :)
    REAL(DP), ALLOCATABLE :: rho_x(:, :)
    !
    ! ... read data from file
    ios = 0
    !
    IF (ionode) THEN
      !
      iun = find_free_unit()
      !
      OPEN(unit=iun, file=filename, status='old', form='formatted', action='read', iostat=ios)
      !
      IF (ios == 0) THEN ! opened
        !
        READ(iun, '(A)') line
        READ(line, *) lambda
        !
        READ(iun, '(A)') line
        READ(line, *) yukawa
        !
        READ(iun, '(A)') line
        READ(line, *) omega
        !
        READ(iun, '(A)') line
        READ(line, *) n_group
        !
        READ(iun, '(A)') line
        READ(line, *) ir1, ir2, ir3
        !
        IF (n_group < 1) THEN
          !
          ios = 1
          !
          CALL infomsg('zmp_initialize', 'there is no group at: ' // TRIM(filename))
          !
        END IF
        !
        IF (ir1 /= dfftp%nr1 .OR. ir2 /= dfftp%nr2 .OR. ir3 /= dfftp%nr3) THEN
          !
          ios = 2
          !
          CALL infomsg('zmp_initialize', 'incorrect FFT-mesh at: ' // TRIM(filename))
          !
        END IF
        !
        IF (ios == 0) THEN ! correct group and mesh
          !
          ALLOCATE(wei_t(dfftp%nr1 * dfftp%nr2 * dfftp%nr3, n_group))
          ALLOCATE(rho_t(dfftp%nr1 * dfftp%nr2 * dfftp%nr3, n_group))
          !
          READ(iun, *) wei_t
          READ(iun, *) rho_t
          !
        END IF ! correct group and mesh
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
      IF (ALLOCATED(wei_t)) DEALLOCATE(wei_t)
      IF (ALLOCATED(rho_t)) DEALLOCATE(rho_t)
      !
      CALL errore('zmp_initialize', 'error to open/read file: ' // TRIM(filename), ios)
      !
      RETURN
      !
    END IF
    !
    ! ... share data for all node
    CALL mp_bcast(lambda,  ionode_id, intra_image_comm)
    CALL mp_bcast(yukawa,  ionode_id, intra_image_comm)
    CALL mp_bcast(omega,   ionode_id, intra_image_comm)
    CALL mp_bcast(n_group, ionode_id, intra_image_comm)
    !
    ALLOCATE(wei_x(dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, n_group))
    ALLOCATE(rho_x(dfftp%nr1x * dfftp%nr2x * dfftp%nr3x, n_group))
    !
    wei_x = 0.0_DP
    rho_x = 0.0_DP
    !
    IF (ionode) THEN
      !
      DO i_group = 1, n_group
        !
        DO ir3 = 1, dfftp%nr3
          !
          DO ir2 = 1, dfftp%nr2
            !
            ir = (ir2 - 1) * dfftp%nr1x + (ir3 - 1) * dfftp%nr1x * dfftp%nr2x
            jr = (ir2 - 1) * dfftp%nr1  + (ir3 - 1) * dfftp%nr1  * dfftp%nr2
            !
            wei_x((ir+1):(ir+dfftp%nr1), i_group) = wei_t((jr+1):(jr+dfftp%nr1), i_group)
            rho_x((ir+1):(ir+dfftp%nr1), i_group) = rho_t((jr+1):(jr+dfftp%nr1), i_group)
            !
          END DO
          !
        END DO
        !
      END DO
      !
    END IF
    !
    ALLOCATE(wei(dfftp%nnr, n_group))
    ALLOCATE(rho(dfftp%nnr, n_group))
    !
#if defined(__MPI)
    CALL scatter_grid(dfftp, wei_x, wei)
    CALL scatter_grid(dfftp, rho_x, rho)
#else
    wei = wei_x
    rho = rho_x
#endif
    !
    DEALLOCATE(wei_x)
    DEALLOCATE(rho_x)
    !
  END SUBROUTINE read_zmp_file
  !
END MODULE zmp_module
