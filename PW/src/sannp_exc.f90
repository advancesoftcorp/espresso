!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_exc(rho, rho_core, rhog_core, rhoe)
  !--------------------------------------------------------------------------
  !
  ! ... calculate exchange-correlation energy density
  !
  USE constants,   ONLY : e2
  USE fft_base,    ONLY : dfftp
  USE gvect,       ONLY : ngm
  USE kinds,       ONLY : DP
  USE lsda_mod,    ONLY : nspin
  USE scf,         ONLY : scf_type
  USE spin_orb,    ONLY : domag
  USE xc_lda_lsda, ONLY : xc
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(IN)  :: rho
  REAL(DP),       INTENT(IN)  :: rho_core(dfftp%nnr)
  COMPLEX(DP),    INTENT(IN)  :: rhog_core(ngm)
  REAL(DP),       INTENT(OUT) :: rhoe(dfftp%nnr)
  !
  INTEGER  :: ir
  REAL(DP) :: arho
  !
  REAL(DP), ALLOCATABLE :: rhox(:,:)
  !
  REAL(DP), ALLOCATABLE :: ex(:)
  REAL(DP), ALLOCATABLE :: ec(:)
  REAL(DP), ALLOCATABLE :: vx(:,:)
  REAL(DP), ALLOCATABLE :: vc(:,:)
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.0E-10_DP
  !
  ALLOCATE(rhox(dfftp%nnr, nspin))
  !
  ALLOCATE(ex(dfftp%nnr))
  ALLOCATE(ec(dfftp%nnr))
  ALLOCATE(vx(dfftp%nnr, nspin))
  ALLOCATE(vc(dfftp%nnr, nspin))
  !
  rhoe(1:dfftp%nnr)    = 0.0_DP
  rhox(1:dfftp%nnr, 1) = rho%of_r(1:dfftp%nnr, 1) + rho_core(1:dfftp%nnr)
  !
  IF (nspin == 1 .OR. (nspin == 4 .AND. .NOT. domag)) THEN
    !
    ! ... spin-unpolarized case
    !
    CALL xc(dfftp%nnr, 1, 1, rhox, ex, ec, vx, vc)
    !
    DO ir = 1, dfftp%nnr
      !
      rhoe(ir) = e2 * (ex(ir) + ec(ir)) * rhox(ir, 1)
      !
    END DO
    !
  ELSE IF (nspin == 2) THEN
    !
    ! ... spin-polarized case
    !
    rhox(1:dfftp%nnr, 2) = rho%of_r(1:dfftp%nnr, 2)
    !
    CALL xc(dfftp%nnr, 2, 2, rhox, ex, ec, vx, vc)
    !
    DO ir = 1, dfftp%nnr
      !
      rhoe(ir) = e2 * (ex(ir) + ec(ir)) * rhox(ir, 1)
      !
    END DO
    !
  ELSE IF (nspin == 4) THEN
    !
    ! ... noncolinear case
    !
    rhox(1:dfftp%nnr, 2) = rho%of_r(1:dfftp%nnr, 2)
    rhox(1:dfftp%nnr, 3) = rho%of_r(1:dfftp%nnr, 3)
    rhox(1:dfftp%nnr, 4) = rho%of_r(1:dfftp%nnr, 4)
    !
    CALL xc(dfftp%nnr, 4, 2, rhox, ex, ec, vx, vc)
    !
    DO ir = 1, dfftp%nnr
      !
      arho = ABS(rhox(ir, 1))
      !
      IF (arho >= vanishing_charge) THEN
        !
        rhoe(ir) = e2 * (ex(ir) + ec(ir)) * arho
        !
      END IF
      !
    END DO
    !
  END IF
  !
  DEALLOCATE(rhox)
  !
  DEALLOCATE(ex)
  DEALLOCATE(ec)
  DEALLOCATE(vx)
  DEALLOCATE(vc)
  !
  ! ... add gradient corrections (if any)
  !
  CALL sannp_gradcorr(rho%of_r, rho%of_g, rho_core, rhog_core, rhoe)
  !
END SUBROUTINE sannp_exc
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_gradcorr(rho, rhog, rho_core, rhog_core, rhoe)
  !--------------------------------------------------------------------------
  !
  ! ... calculate exchange-correlation energy density, for GGA
  !
  USE constants,      ONLY : e2
  USE fft_base,       ONLY : dfftp
  USE fft_rho,        ONLY : rho_r2g
  USE funct,          ONLY : dft_is_gradient
  USE gvect,          ONLY : ngm, g
  USE kinds,          ONLY : DP
  USE lsda_mod,       ONLY : nspin
  USE spin_orb,       ONLY : domag
  USE xc_gga,         ONLY : xc_gcx
  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN)    :: rho(dfftp%nnr, nspin)
  COMPLEX(DP), INTENT(IN)    :: rhog(ngm, nspin)
  REAL(DP),    INTENT(IN)    :: rho_core(dfftp%nnr)
  COMPLEX(DP), INTENT(IN)    :: rhog_core(ngm)
  REAL(DP),    INTENT(INOUT) :: rhoe(dfftp%nnr)
  !
  INTEGER  :: ir
  INTEGER  :: is, nspin0
  REAL(DP) :: fac
  REAL(DP) :: sgn(2)
  !
  REAL(DP),    ALLOCATABLE :: grho(:,:,:)
  REAL(DP),    ALLOCATABLE :: rhoaux(:,:)
  REAL(DP),    ALLOCATABLE :: segni(:)
  COMPLEX(DP), ALLOCATABLE :: rhogaux(:,:)
  !
  REAL(DP),    ALLOCATABLE :: sx(:)
  REAL(DP),    ALLOCATABLE :: sc(:)
  REAL(DP),    ALLOCATABLE :: v1x(:,:)
  REAL(DP),    ALLOCATABLE :: v2x(:,:)
  REAL(DP),    ALLOCATABLE :: v1c(:,:)
  REAL(DP),    ALLOCATABLE :: v2c(:,:)
  REAL(DP),    ALLOCATABLE :: v2c_ud(:)
  !
  IF (.NOT. dft_is_gradient()) RETURN
  !
  nspin0 = nspin
  IF (nspin == 4)             nspin0 = 1
  IF (nspin == 4 .AND. domag) nspin0 = 2
  !
  fac = 1.0_DP / DBLE(nspin0)
  !
  sgn(1) = +1.0_DP
  sgn(2) = -1.0_DP
  !
  ALLOCATE(grho(3, dfftp%nnr, nspin0))
  ALLOCATE(rhoaux( dfftp%nnr, nspin0))
  IF (nspin == 4 .AND. domag) ALLOCATE(segni(dfftp%nnr))
  ALLOCATE(rhogaux(ngm, nspin0))
  !
  ALLOCATE(sx(dfftp%nnr))
  ALLOCATE(sc(dfftp%nnr))
  ALLOCATE(v1x(dfftp%nnr, nspin0))
  ALLOCATE(v2x(dfftp%nnr, nspin0))
  ALLOCATE(v1c(dfftp%nnr, nspin0))
  ALLOCATE(v2c(dfftp%nnr, nspin0))
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  IF (nspin == 4 .AND. domag) THEN
    !
    !
    CALL compute_rho(rho, rhoaux, segni, dfftp%nnr)
    !
    ! ... bring starting rhoaux to G-space
    !
    CALL rho_r2g(dfftp, rhoaux(:, 1:nspin0), rhogaux(:, 1:nspin0))
    !
  ELSE
    !
    ! ... for convenience rhoaux and rhogaux are in (up,down) format, when LSDA
    !
    DO is = 1, nspin0
      !
      rhoaux (:, is) = (rho (:, 1) + sgn(is) * rho (:, nspin0)) * 0.5_DP
      rhogaux(:, is) = (rhog(:, 1) + sgn(is) * rhog(:, nspin0)) * 0.5_DP
      !
    END DO
    !
  END IF
  !
  DO is = 1, nspin0
    !
    rhoaux (:, is) = fac * rho_core (:) + rhoaux (:, is)
    rhogaux(:, is) = fac * rhog_core(:) + rhogaux(:, is)
    !
    CALL fft_gradient_g2r(dfftp, rhogaux(1, is), g, grho(1, 1, is))
    !
  END DO
  !
  IF (nspin0 == 1) THEN
    !
    ! ... This is the spin-unpolarised case
    !
    CALL xc_gcx(dfftp%nnr, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c)
    !
    DO ir = 1, dfftp%nnr
      !
      rhoe(ir) = rhoe(ir) + e2 * (sx(ir) + sc(ir))
      !
    END DO
    !
  ELSE
    !
    ! ... spin-polarised case
    !
    ALLOCATE(v2c_ud(dfftp%nnr))
    !
    CALL xc_gcx(dfftp%nnr, nspin0, rhoaux, grho, sx, sc, v1x, v2x, v1c, v2c, v2c_ud)
    !
    DO ir = 1, dfftp%nnr
      !
      rhoe(ir) = rhoe(ir) + e2 * (sx(ir) + sc(ir))
      !
    END DO
    !
    DEALLOCATE(v2c_ud)
    !
  END IF
  !
  DEALLOCATE(grho)
  DEALLOCATE(rhoaux)
  IF (nspin == 4 .AND. domag ) DEALLOCATE(segni)
  DEALLOCATE(rhogaux)
  !
  DEALLOCATE(sx)
  DEALLOCATE(sc)
  DEALLOCATE(v1x)
  DEALLOCATE(v2x)
  DEALLOCATE(v1c)
  DEALLOCATE(v2c)
  !
END SUBROUTINE sannp_gradcorr

