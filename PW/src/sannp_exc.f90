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
  USE constants, ONLY : e2
  USE fft_base,  ONLY : dfftp
  USE funct,     ONLY : xc, xc_spin
  USE gvect,     ONLY : ngm
  USE kinds,     ONLY : DP
  USE lsda_mod,  ONLY : nspin
  USE scf,       ONLY : scf_type
  USE spin_orb,  ONLY : domag
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(IN)  :: rho
  REAL(DP),       INTENT(IN)  :: rho_core(dfftp%nnr)
  COMPLEX(DP),    INTENT(IN)  :: rhog_core(ngm)
  REAL(DP),       INTENT(OUT) :: rhoe(dfftp%nnr)
  !
  INTEGER  :: ir
  REAL(DP) :: rhox, arhox, zeta
  REAL(DP) :: amag, xmag, ymag, zmag
  REAL(DP) :: ex, ec, vx(2), vc(2)
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.0E-10_DP
  REAL(DP), PARAMETER :: vanishing_mag    = 1.0E-20_DP
  !
  rhoe(:) = 0.0_DP
  !
  IF (nspin == 1 .OR. (nspin == 4 .AND. .NOT. domag)) THEN
    !
    ! ... spin-unpolarized case
    !
    DO ir = 1, dfftp%nnr
      !
      rhox = rho%of_r(ir, 1) + rho_core(ir)
      !
      arhox = ABS(rhox)
      !
      IF (arhox > vanishing_charge) THEN
        !
        CALL xc(arhox, ex, ec, vx(1), vc(1))
        !
        rhoe(ir) = e2 * (ex + ec) * rhox
        !
      END IF
      !
    END DO
    !
  ELSE IF (nspin == 2) THEN
    !
    ! ... spin-polarized case
    !
    DO ir = 1, dfftp%nnr
      !
      rhox = rho%of_r(ir, 1) + rho_core(ir)
      !
      arhox = ABS(rhox)
      !
      IF (arhox > vanishing_charge) THEN
        !
        zeta = rho%of_r(ir, 2) / arhox
        !
        IF (ABS(zeta) > 1.0_DP) zeta = SIGN(1.0_DP, zeta)
        !
        CALL xc_spin(arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2))
        !
        rhoe(ir) = e2 * (ex + ec) * rhox
        !
      END IF
      !
    END DO
    !
  ELSE IF (nspin == 4) THEN
    !
    ! ... noncolinear case
    !
    DO ir = 1, dfftp%nnr
      !
      xmag = rho%of_r(ir, 2)
      ymag = rho%of_r(ir, 3)
      zmag = rho%of_r(ir, 4)
      amag = SQRT(xmag * xmag + ymag * ymag + zmag * zmag)
      !
      rhox = rho%of_r(ir, 1) + rho_core(ir)
      !
      arhox = ABS(rhox)
      !
      IF (arhox > vanishing_charge) THEN
        !
        zeta = amag / arhox
        !
        IF (ABS(zeta) > 1.0_DP) zeta = SIGN(1.0_DP, zeta)
        !
        CALL xc_spin(arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2))
        !
        rhoe(ir) = e2 * (ex + ec) * rhox
        !
      END IF
      !
     END DO
     !
  END IF
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
  USE fft_interfaces, ONLY : fwfft
  USE funct,          ONLY : gcxc, gcx_spin, gcc_spin, gcc_spin_more, &
                           & igcc_is_lyp, dft_is_gradient
  USE gvect,          ONLY : ngm, g
  USE kinds,          ONLY : DP
  USE lsda_mod,       ONLY : nspin
  USE spin_orb,       ONLY : domag
  USE wavefunctions,  ONLY : psic
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
  REAL(DP) :: arho, segno, zeta, rh, grh2, sgn(2)
  REAL(DP) :: grho2(2), grhox(2), grhoy(2), grhoz(2)
  REAL(DP) :: grhoup, grhodw, grhoud
  REAL(DP) :: sx, sc, v1x, v2x, v1c, v2c
  REAL(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw
  REAL(DP) :: v2cup, v2cdw, v2cud, rup, rdw
  !
  REAL(DP),    ALLOCATABLE :: grho(:,:,:)
  REAL(DP),    ALLOCATABLE :: rhoaux(:,:)
  REAL(DP),    ALLOCATABLE :: segni(:)
  COMPLEX(DP), ALLOCATABLE :: rhogaux(:,:)
  !
  REAL(DP), PARAMETER :: epsr = 1.0E-6_DP
  REAL(DP), PARAMETER :: epsg = 1.0E-10_DP
  !
  IF (.NOT. dft_is_gradient()) RETURN
  !
  nspin0 = nspin
  IF (nspin == 4) nspin0 = 1
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
  ! ... calculate the gradient of rho + rho_core in real space
  !
  IF (nspin == 4 .AND. domag) THEN
    !
    CALL compute_rho(rho, rhoaux, segni, dfftp%nnr)
    !
    ! ... bring starting rhoaux to G-space
    !
    DO is = 1, nspin0
      !
      psic(:) = rhoaux(:, is)
      !
      CALL fwfft('Rho', psic, dfftp)
      !
      rhogaux(:, is) = psic(dfftp%nl(:))
      !
    END DO
    !
  ELSE
    !
    ! ... for convenience rhoaux and rhogaux are in (up,down) format, if LSDA
    !
    DO is = 1, nspin0
      rhoaux (:, is) = (rho (:, 1) + sgn(is) * rho (:, nspin0)) * 0.5_DP
      rhogaux(:, is) = (rhog(:, 1) + sgn(is) * rhog(:, nspin0)) * 0.5_DP
    END DO
    !
  END IF
  !
  DO is = 1, nspin0
    !
    rhoaux (:,is) = fac * rho_core (:) + rhoaux (:,is)
    rhogaux(:,is) = fac * rhog_core(:) + rhogaux(:,is)
    !
    CALL fft_gradient_g2r(dfftp, rhogaux(1, is), g, grho(1, 1, is))
    !
  END DO
  !
  IF (nspin0 == 1) THEN
    !
    ! ... spin-unpolarised case
    !
    DO ir = 1, dfftp%nnr
      !
      arho = ABS(rhoaux(ir, 1))
      !
      IF (arho > epsr) THEN
        !
        grhox(1) = grho(1, ir, 1)
        grhoy(1) = grho(2, ir, 1)
        grhoz(1) = grho(3, ir, 1)
        grho2(1) = grhox(1) * grhox(1) + grhoy(1) * grhoy(1) + grhoz(1) * grhoz(1)
        !
        IF (grho2(1) > epsg) THEN
          !
          segno = SIGN(1.0_DP, rhoaux(ir, 1))
          !
          CALL gcxc(arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c)
          !
          rhoe(ir) = rhoe(ir) + e2 * (sx + sc) * segno
          !
        END IF
        !
      END IF
      !
    END DO
    !
  ELSE
    !
    ! ... spin-polarised case
    !
    DO ir = 1, dfftp%nnr
      !
      rh = rhoaux(ir, 1) + rhoaux(ir, 2)
      !
      grhox(:) = grho(1, ir, :)
      grhoy(:) = grho(2, ir, :)
      grhoz(:) = grho(3, ir, :)
      grho2(:) = grhox(:) * grhox(:) + grhoy(:) * grhoy(:) + grhoz(:) * grhoz(:)
      !
      CALL gcx_spin(rhoaux(ir, 1), rhoaux(ir, 2), grho2(1), grho2(2), &
                  & sx, v1xup, v1xdw, v2xup, v2xdw)
      !
      IF (rh > epsr) THEN
        !
        IF (igcc_is_lyp()) THEN
          !
          rup = rhoaux(ir, 1)
          rdw = rhoaux(ir, 2)
          !
          grhoup = grho2(1)
          grhodw = grho2(2)
          grhoud = grhox(1) * grhox(2) + grhoy(1) * grhoy(2) + grhoz(1) * grhoz(2)
          !
          CALL gcc_spin_more(rup, rdw, grhoup, grhodw, grhoud, &
                           & sc, v1cup, v1cdw, v2cup, v2cdw, v2cud)
          !
        ELSE
          !
          zeta = (rhoaux(ir, 1) - rhoaux(ir, 2)) / rh
          IF (nspin == 4 .AND. domag) zeta = ABS(zeta) * segni(ir)
          !
          grhox(1) = grhox(1) + grhox(2)
          grhoy(1) = grhoy(1) + grhoy(2)
          grhoz(1) = grhoz(1) + grhoz(2)
          !
          grh2 = grhox(1) * grhox(1) + grhoy(1) * grhoy(1) + grhoz(1) * grhoz(1)
          !
          CALL gcc_spin(rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
          !
        END IF
        !
      ELSE
        !
        sc = 0.0_DP
        !
      END IF
      !
      rhoe(ir) = rhoe(ir) + e2 * (sx + sc)
      !
    END DO
    !
  END IF
  !
  DEALLOCATE(grho)
  DEALLOCATE(rhoaux)
  IF (nspin == 4 .AND. domag ) DEALLOCATE(segni)
  DEALLOCATE(rhogaux)
  !
END SUBROUTINE sannp_gradcorr

