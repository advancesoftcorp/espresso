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
  ! ... The coulombic interaction is screened as Yukawa, Erfc/r or Gaussian.
  !
  USE constants,      ONLY : pi, fpi, e2
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
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE mp_images,      ONLY : intra_image_comm
  USE paw_variables,  ONLY : okpaw
  USE scatter_mod,    ONLY : scatter_grid
  USE uspp,           ONLY : okvan
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  LOGICAL               :: do_zmp     = .FALSE.
  CHARACTER(LEN=256)    :: filzmp     = ''
  REAL(DP)              :: zmp_mixing = 1.0_DP
  INTEGER               :: n_group    = 0
  INTEGER               :: kernel
  REAL(DP)              :: lambda
  REAL(DP)              :: omega ! in 1/Bohr
  REAL(DP), ALLOCATABLE :: wei_group(:, :) ! (dfftp%nnr, n_group)
  REAL(DP), ALLOCATABLE :: rho_group(:, :) ! (dfftp%nnr, n_group)
  !
  PUBLIC :: do_zmp
  PUBLIC :: filzmp
  PUBLIC :: zmp_mixing
  PUBLIC :: zmp_initialize
  PUBLIC :: zmp_finalize
  PUBLIC :: add_vzmp
  PUBLIC :: mix_rho_zmp
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
    !IF (okvan) THEN
    !  CALL errore('zmp_initialize', 'ZMP does not support USPP', 1)
    !END IF
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
    REAL(DP) :: fac1, fac2
    REAL(DP) :: ww
    REAL(DP) :: beta
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
    ALLOCATE(drho(dfftp%nnr))
    ALLOCATE(rhog(dfftp%nnr))
    ALLOCATE(vg  (dfftp%nnr))
    ALLOCATE(aux (dfftp%nnr))
    !
    ! ... calculate potentials for each group
    fac1 = e2 * fpi / tpiba2
    fac2 = e2 * (pi ** (3.0_DP / 2.0_DP)) / tpiba2
    ww   = omega * omega / tpiba2
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
      IF (kernel == 1) THEN
        !
        ! ... exp(-w*r)/r  ->  4pi/(g2+w2)
        !
        DO ig = gstart, ngm
          !
          vg(ig) = rhog(ig) * fac1 / (gg(ig) + ww)
          !
        END DO
        !
        IF (gstart > 1) THEN
          !
          vg(1) = rhog(1) * fac1 / ww
          !
        END IF
        !
      ELSE IF (kernel == 2) THEN
        !
        ! ... erfc(w*r)/r  ->  4pi*(1-exp(-g2/(4*w2)))/g2,  if g > 0
        !                      4pi/(4*w2),                  if g = 0
        DO ig = gstart, ngm
          !
          vg(ig) = rhog(ig) * fac1 * (1.0_DP - EXP(-0.25_DP * gg(ig) / ww)) / gg(ig)
          !
        END DO
        !
        IF (gstart > 1) THEN
          !
          vg(1) = rhog(1) * fac1 * 0.25_DP / ww
          !
        END IF
        !
      ELSE IF (kernel == 3) THEN
        !
        ! ... exp(-(w*r)^2)*w  ->  (pi^(3/2))*exp(-g2/(4*w2))/w2
        !
        DO ig = gstart, ngm
          !
          vg(ig) = rhog(ig) * fac2 * EXP(-0.25_DP * gg(ig) / ww) / ww
          !
        END DO
        !
        IF (gstart > 1) THEN
          !
          vg(1) = rhog(1) * fac2 / ww
          !
        END IF
        !
      ELSE
        !
        CALL errore('add_vzmp', 'incorrect type of kernel', MAX(1, kernel))
        !
      END IF
      !
      aux(:) = CMPLX(0.0_DP, 0.0_DP, KIND=DP)
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
  SUBROUTINE mix_rho_zmp(rho)
    !----------------------------------------------------------------------------
    !
    ! ... mixing charge density of ZMP
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: rho(dfftp%nnr)
    !
    INTEGER  :: ir
    INTEGER  :: i_group
    REAL(DP) :: beta
    REAL(DP) :: wei_sum
    REAL(DP) :: wei_tot
    REAL(DP) :: rho_tot
    REAL(DP) :: coef1, coef2
    REAL(DP) :: sum0, sum1, sum2
    !
    REAL(DP), ALLOCATABLE :: rho1(:)
    REAL(DP), ALLOCATABLE :: rho2(:)
    !
    IF (.NOT. do_zmp) RETURN
    !
    IF (zmp_mixing >= 1.0_DP) RETURN
    !
    ! ... read and initialize data
    IF (n_group < 1) THEN
      !
      CALL zmp_initialize()
      !
    END IF
    !
    ! ... mixing charge
    beta = MIN(MAX(0.0_DP, zmp_mixing), 1.0_DP)
    !
    ALLOCATE(rho1(dfftp%nnr))
    ALLOCATE(rho2(dfftp%nnr))
    !
    rho1 = rho
    rho2 = 0.0_DP
    !
    DO ir = 1, dfftp%nnr
      !
      wei_sum = 0.0_DP
      !
      DO i_group = 1, n_group
        !
        wei_sum = wei_sum + wei_group(ir, i_group)
        !
      END DO
      !
      IF (wei_sum < 1.0E-8_DP) CYCLE
      !
      wei_tot = 0.0_DP
      !
      rho_tot = 0.0_DP
      !
      DO i_group = 1, n_group
        !
        wei_tot = MAX(wei_tot, wei_group(ir, i_group))
        !
        rho_tot = rho_tot + rho_group(ir, i_group) * wei_group(ir, i_group) / wei_sum
        !
      END DO
      !
      coef2 = (1.0_DP - beta) * TANH(wei_tot)
      coef1 = 1.0_DP - coef2
      !
      rho1(ir) = coef1 * rho(ir)
      rho2(ir) = coef2 * rho_tot
      !
    END DO
    !
    sum0 = SUM(rho)
    sum1 = SUM(rho1)
    sum2 = SUM(rho2)
    !
    CALL mp_sum(sum0, intra_bgrp_comm)
    CALL mp_sum(sum1, intra_bgrp_comm)
    CALL mp_sum(sum2, intra_bgrp_comm)
    !
    rho(:) = ((sum0 - sum2) / sum1) * rho1(:) + rho2(:)
    !
    DEALLOCATE(rho1)
    DEALLOCATE(rho2)
    !
  END SUBROUTINE mix_rho_zmp
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
        READ(line, *) kernel
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
        READ(iun, '()')
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
    CALL mp_bcast(kernel,  ionode_id, intra_image_comm)
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
    wei_group = 0.0_DP
    rho_group = 0.0_DP
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
