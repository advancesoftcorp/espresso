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
  USE constants,      ONLY : fpi, e2
  USE control_flags,  ONLY : gamma_only
  USE cell_base,      ONLY : tpiba2
  USE fft_base,       ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE force_mod,      ONLY : lforce, lstres
  USE gvect,          ONLY : ngm, gg, gstart
  USE io_files,       ONLY : tmp_dir, prefix
  USE io_global,      ONLY : ionode, ionode_id
  USE kinds,          ONLY : DP
  USE lsda_mod,       ONLY : nspin
  USE mp,             ONLY : mp_sum, mp_bcast, mp_barrier
  USE mp_images,      ONLY : intra_image_comm
  USE paw_variables,  ONLY : okpaw
  USE scatter_mod,    ONLY : scatter_grid
  USE uspp,           ONLY : okvan
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
  REAL(DP), ALLOCATABLE :: wei_group(:, :) ! (dfftp%nnr, n_group)
  REAL(DP), ALLOCATABLE :: rho_group(:, :) ! (dfftp%nnr, n_group)
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
    IF (lforce) THEN
      CALL errore('zmp_initialize', 'you cannot calculate force for ZMP', 1)
    END IF
    !
    IF (lstres) THEN
      CALL errore('zmp_initialize', 'you cannot calculate stress for ZMP', 1)
    END IF
    !
    IF (nspin /= 1) THEN
      CALL errore('zmp_initialize', 'ZMP does not support spin-polarized calculation', 1)
    END IF
    !
    IF (okvan) THEN
      CALL errore('zmp_initialize', 'ZMP does not support USPP', 1)
    END IF
    !
    IF (okpaw) THEN
      CALL errore('zmp_initialize', 'ZMP does not support PAW', 1)
    END IF
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
    DEALLOCATE(wei_group)
    DEALLOCATE(rho_group)
    !
  END SUBROUTINE zmp_finalize
  !
  !----------------------------------------------------------------------------
  SUBROUTINE add_vzmp(v, rho)
    !----------------------------------------------------------------------------
    !
    ! ... add ZMP-potential
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: v  (dfftp%nnr)
    REAL(DP), INTENT(IN)    :: rho(dfftp%nnr)
    !
    INTEGER  :: ig
    INTEGER  :: i_group
    REAL(DP) :: fac
    !
    REAL(DP),    ALLOCATABLE :: drho(:)
    COMPLEX(DP), ALLOCATABLE :: rhog(:)
    COMPLEX(DP), ALLOCATABLE :: vg  (:)
    COMPLEX(DP), ALLOCATABLE :: aux (:)
    !
    IF (.NOT. do_zmp) RETURN
    !
    ! ... read and initialize data
    IF (n_group < 1) THEN
      !
      CALL zmp_initialize()
      !
    END IF
    !
    ! ... calculate and add potentials for each group
    ALLOCATE(drho(dfftp%nnr))
    ALLOCATE(rhog(dfftp%nnr))
    ALLOCATE(vg  (dfftp%nnr))
    ALLOCATE(aux (dfftp%nnr))
    !
    fac = e2 * fpi / tpiba2
    !
    DO i_group = 1, n_group
      !
      drho(:) = wei_group(:, i_group) * (rho(:) - rho_group(:, i_group))
      !
      aux(:) = drho(:)
      !
      CALL fwfft('Rho', aux, dfftp)
      !
      rhog(1:ngm) = aux(dfftp%nl(1:ngm))
      !
      DO ig = gstart, ngm
        !
        ! >>> TODO
        ! >>> TODO
        ! >>> TODO
        !
        vg(ig) = rhog(ig) * fac / gg(ig)
        !
        ! <<< TODO
        ! <<< TODO
        ! <<< TODO
        !
      END DO
      !
      aux(:) = 0.0_DP
      !
      aux(dfftp%nl(1:ngm)) = vg(1:ngm)
      !
      IF (gamma_only) THEN
        !
        aux(dfftp%nlm(1:ngm)) = CONJG(vg(1:ngm))
        !
      END IF
      !
      CALL invfft('Rho', aux, dfftp)
      !
      v(:) = v(:) + lambda * wei_group(:, i_group) * DBLE(aux(:))
      !
    END DO
    !
    DEALLOCATE(drho)
    DEALLOCATE(rhog)
    DEALLOCATE(vg)
    DEALLOCATE(aux)
    !
  END SUBROUTINE add_vzmp
  !
  !----------------------------------------------------------------------------
  SUBROUTINE read_zmp_file(filename)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(IN) :: filename
    !
    INTEGER :: iun
    INTEGER :: ios
    INTEGER :: ir1, ir2, ir3
    INTEGER :: ir, jr
    INTEGER :: i_group
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
      DEALLOCATE(wei_t)
      DEALLOCATE(rho_t)
      !
    END IF
    !
    ALLOCATE(wei_group(dfftp%nnr, n_group))
    ALLOCATE(rho_group(dfftp%nnr, n_group))
    !
#if defined(__MPI)
    DO i_group = 1, n_group
      CALL scatter_grid(dfftp, wei_x(:, i_group), wei_group(:, i_group))
      CALL scatter_grid(dfftp, rho_x(:, i_group), rho_group(:, i_group))
    END DO
#else
    wei_group = wei_x
    rho_group = rho_x
#endif
    !
    DEALLOCATE(wei_x)
    DEALLOCATE(rho_x)
    !
  END SUBROUTINE read_zmp_file
  !
END MODULE zmp_module
