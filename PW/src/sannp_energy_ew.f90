!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_energy_ew(energy, ia_start, ia_end)
  !--------------------------------------------------------------------------
  !
  ! ... calculate Ewald energy for single atom
  !
  USE cell_base,     ONLY : at, bg, alat, tpiba2, omega
  USE constants,     ONLY : tpi, fpi, sqrtpi, e2
  USE control_flags, ONLY : gamma_only
  USE gvect,         ONLY : ngm, gstart, gg, gcutm, &
                          & mill, eigts1, eigts2, eigts3
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, zv, tau
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE vlocal,        ONLY : strf
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: energy(nat)
  INTEGER,  INTENT(IN)    :: ia_start
  INTEGER,  INTENT(IN)    :: ia_end
  !
  INTEGER, PARAMETER :: mxr = 50
  !
  INTEGER  :: ia, ita
  INTEGER  :: ib, itb
  INTEGER  :: ig
  INTEGER  :: ir, nr
  REAL(DP) :: za, zb
  REAL(DP) :: ewald
  REAL(DP) :: ewaldg
  REAL(DP) :: ewaldr
  REAL(DP) :: charge
  REAL(DP) :: upperbound
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
  ! ... detect coulombic screening alpha
  !
  charge = 0.0_DP
  alpha  = 0.0_DP
  !
  IF (gstart == 2) THEN
    !
    DO ia = 1, nat
      ita = ityp(ia)
      za  = zv(ita)
      charge = charge + za
    END DO
    !
    alpha = 2.9_DP
10  alpha = alpha - 0.1_DP
    !
    IF (alpha <= 0.0_DP) THEN
      CALL errore ('sannp_energy_ew', 'optimal alpha not found', 1)
    END IF
    !
    upperbound = 2.0_DP * charge * charge * SQRT(2.0_DP * alpha / tpi) &
             & * qe_erfc(SQRT(tpiba2 * gcutm / 4.0_DP / alpha))
    !
    IF (upperbound > 1.0E-7_DP) GOTO 10
    !
  END IF
  !
  CALL mp_sum(charge, intra_bgrp_comm)
  CALL mp_sum(alpha,  intra_bgrp_comm)
  !
  alpha0 = SQRT(alpha)
  !
  rmax = 4.0_DP / alpha0 / alat
  !
  ! ... long-range electrostatic potential from ions
  !
  IF (gstart == 2) THEN
    vion(1) = (0.0_DP, 0.0_DP)
  END IF
  !
  DO ig = gstart, ngm
    rhon = (0.0_DP, 0.0_DP)
    DO ita = 1, ntyp
      za = zv(ita)
      rhon = rhon + za * strf(ig, ita)
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
    ita = ityp(ia)
    za  = zv(ita)
    !
    ! ... G-space
    !
    ewaldg = 0.0_DP
    !
    IF (gstart == 2) THEN
      ewaldg = -za * charge / alpha / 4.0_DP
    END IF
    !
    DO ig = gstart, ngm
      rhon = za * eigts1(mill(1, ig), ia) * &
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
    END DO
    !
    ewaldg = ewaldg * fpi / omega
    !
    IF (gstart == 2) THEN
      ewaldg = ewaldg - za * za * 2.0_DP * alpha0 / sqrtpi
    END IF
    !
    ! ... R-space
    !
    ewaldr = 0.0_DP
    !
    IF (gstart == 2) THEN
      !
      DO ib = 1, nat
        itb = ityp(ib)
        zb  = zv(itb)
        !
        dtau(:) = tau(:, ia) - tau(:, ib)
        !
        CALL rgen(dtau, rmax, mxr, at, bg, r, r2, nr)
        !
        DO ir = 1, nr
          rr = SQRT(r2(ir)) * alat
          ewaldr = ewaldr + za * zb * qe_erfc(alpha0 * rr) / rr
        END DO
      END DO
      !
    END IF
    !
    ewald = 0.5_DP * e2 * (ewaldg + ewaldr)
    !
    CALL mp_sum(ewald, intra_bgrp_comm)
    !
    energy(ia) = energy(ia) + ewald
    !
  END DO
  !
  DEALLOCATE(vion)
  !
END SUBROUTINE sannp_energy_ew
