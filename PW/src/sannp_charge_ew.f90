!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_charge_ew(energy_q, force_q, charge, ia_start, ia_end)
  !--------------------------------------------------------------------------
  !
  ! ... calculate Ewald energy and force from atomic charge
  !
  USE cell_base,     ONLY : at, bg, alat, tpiba, tpiba2, omega
  USE constants,     ONLY : fpi, sqrtpi, e2
  USE control_flags, ONLY : gamma_only
  USE gvect,         ONLY : ngm, gstart, gg, g, &
                          & mill, eigts1, eigts2, eigts3
  USE ions_base,     ONLY : nat, tau
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE mp_bands,      ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: energy_q(nat)
  REAL(DP), INTENT(INOUT) :: force_q (3, nat)
  REAL(DP), INTENT(IN)    :: charge  (nat)
  INTEGER,  INTENT(IN)    :: ia_start
  INTEGER,  INTENT(IN)    :: ia_end
  !
  INTEGER, PARAMETER :: mxr = 50
  !
  INTEGER  :: ia, ib
  INTEGER  :: ig
  INTEGER  :: ir, nr
  REAL(DP) :: qa, qb
  REAL(DP) :: ewald
  REAL(DP) :: ewaldg
  REAL(DP) :: ewaldr
  REAL(DP) :: fewald (3)
  REAL(DP) :: fewaldg(3)
  REAL(DP) :: fewaldr(3)
  REAL(DP) :: fewaldr0
  REAL(DP) :: alpha, alpha0
  REAL(DP) :: rmax
  REAL(DP) :: dtau(3)
  REAL(DP) :: r(3, mxr)
  REAL(DP) :: r2(mxr)
  REAL(DP) :: rr
  REAL(DP) :: fact
  REAL(DP) :: rhonr, rhoni
  REAL(DP) :: vionr, vioni
  !
  COMPLEX(DP)              :: rhon
  COMPLEX(DP), ALLOCATABLE :: vion(:)
  !
  REAL(DP), EXTERNAL :: qe_erfc
  !
  ALLOCATE(vion(ngm))
  !
  ! ... coulombic screening alpha
  !
  alpha  = 1.0_DP
  alpha0 = SQRT(alpha)
  !
  rmax = 5.0_DP / alpha0 / alat
  !
  ! ... long-range electrostatic potential from ions
  !
  IF (gstart == 2) THEN
    vion(1) = (0.0_DP, 0.0_DP)
  END IF
  !
  DO ig = gstart, ngm
    rhon = (0.0_DP, 0.0_DP)
    DO ia = 1, nat
      qa = charge(ia)
      rhon = rhon + qa * eigts1(mill(1, ig), ia) * &
                       & eigts2(mill(2, ig), ia) * &
                       & eigts3(mill(3, ig), ia)
    END DO
    !
    vion(ig) = rhon &
           & * EXP(-gg(ig) * tpiba2 / alpha / 4.0_DP) / gg(ig) / tpiba2
  END DO
  !
  ! ... Ewald for each atom
  !
  IF (gamma_only) THEN
    fact = 2.0_DP
  ELSE
    fact = 1.0_DP
  END IF
  !
  DO ia = ia_start, ia_end
    !
    qa = charge(ia)
    !
    ! ... G-space
    !
    ewaldg     = 0.0_DP
    fewaldg(:) = 0.0_DP
    !
    DO ig = gstart, ngm
      rhon = qa * eigts1(mill(1, ig), ia) * &
                & eigts2(mill(2, ig), ia) * &
                & eigts3(mill(3, ig), ia)
      !
      rhonr =  DBLE(rhon)
      rhoni = AIMAG(rhon)
      !
      vionr =  DBLE(vion(ig))
      vioni = AIMAG(vion(ig))
      !
      ewaldg = ewaldg + fact * (rhonr * vionr + rhoni * vioni)
      !
      fewaldg(:) = fewaldg(:) - fact * (rhoni * vionr - rhonr * vioni) &
                                   & * g(:, ig) * tpiba
    END DO
    !
    ewaldg     = ewaldg     * fpi / omega
    fewaldg(:) = fewaldg(:) * fpi / omega
    !
    IF (gstart == 2) THEN
      ewaldg = ewaldg - qa * qa * 2.0_DP * alpha0 / sqrtpi
    END IF
    !
    ! ... R-space
    !
    ewaldr     = 0.0_DP
    fewaldr(:) = 0.0_DP
    !
    IF (gstart == 2) THEN
      !
      DO ib = 1, nat
        qb = charge(ib)
        !
        dtau(:) = tau(:, ia) - tau(:, ib)
        !
        CALL rgen(dtau, rmax, mxr, at, bg, r, r2, nr)
        !
        DO ir = 1, nr
          rr = SQRT(r2(ir)) * alat
          ewaldr = ewaldr + qa * qb * qe_erfc(alpha0 * rr) / rr
        END DO
        !
        IF (ib == ia) CYCLE
        !
        DO ir = 1, nr
          rr = SQRT(r2(ir)) * alat
          fewaldr0 = qa * qb / rr * (qe_erfc(alpha0 * rr) / rr &
                 & + 2.0d0 * alpha0 / sqrtpi * exp(-alpha * rr * rr))
          fewaldr(:) = fewaldr(:) - fewaldr0 * r(:, ir) * alat / rr
        END DO
      END DO
      !
    END IF
    !
    ewald = 0.5_DP * e2 * (ewaldg + ewaldr)
    !
    fewald(:) = e2 * (fewaldg(:) + fewaldr(:))
    !
    CALL mp_sum(ewald,  intra_bgrp_comm)
    CALL mp_sum(fewald, intra_bgrp_comm)
    !
    energy_q(ia) = ewald
    !
    force_q(:, ia) = fewald(:)
    !
  END DO
  !
  DEALLOCATE(vion)
  !
END SUBROUTINE sannp_charge_ew
