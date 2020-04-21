!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_energy_lc(energy, charge, rhokin, ia_start, ia_end)
  !--------------------------------------------------------------------------
  !
  ! ... calculate local energy for single atom,
  ! ... using "Hirshfeld" spatial decomposition
  !
  USE cell_base, ONLY : omega
  USE fft_base,  ONLY : dfftp
  USE gvect,     ONLY : ngl
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, zv
  USE kinds,     ONLY : DP
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE scf,       ONLY : rho, rho_core, rhog_core, vltot
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: energy(nat)
  REAL(DP), INTENT(INOUT) :: charge(nat)
  REAL(DP), INTENT(IN)    :: rhokin(dfftp%nnr)
  INTEGER,  INTENT(IN)    :: ia_start
  INTEGER,  INTENT(IN)    :: ia_end
  !
  INTEGER  :: mat
  INTEGER  :: ia
  INTEGER  :: it
  INTEGER  :: ir
  REAL(DP) :: e
  REAL(DP) :: q
  REAL(DP) :: weight
  REAL(DP) :: rhotot0
  REAL(DP) :: vtot0
  REAL(DP) :: rscale
  REAL(DP) :: vscale
  !
  REAL(DP), ALLOCATABLE :: rhotyp(:,:)
  REAL(DP), ALLOCATABLE :: rhotot(:)
  REAL(DP), ALLOCATABLE :: rhoatm(:)
  REAL(DP), ALLOCATABLE :: rhoexc(:)
  REAL(DP), ALLOCATABLE :: vtot  (:)
  REAL(DP), ALLOCATABLE :: vatm  (:)
  REAL(DP), ALLOCATABLE :: vhart (:)
  !
  REAL(DP), PARAMETER :: ALPHA = 10.0_DP
  REAL(DP), PARAMETER :: RHO0  = 1.0E-8_DP
  REAL(DP), PARAMETER :: RHO1  = 1.0E-6_DP
  REAL(DP), PARAMETER :: LRHO1 = LOG(RHO1)
  REAL(DP), PARAMETER :: V0    = 1.0E-8_DP
  !
  REAL(DP), EXTERNAL  :: qe_erfc
  !
  ALLOCATE(rhotyp(ngl, ntyp))
  ALLOCATE(rhotot(dfftp%nnr))
  ALLOCATE(rhoatm(dfftp%nnr))
  ALLOCATE(rhoexc(dfftp%nnr))
  ALLOCATE(vtot  (dfftp%nnr))
  ALLOCATE(vatm  (dfftp%nnr))
  ALLOCATE(vhart (dfftp%nnr))
  !
  CALL sannp_rhoa_total(ALPHA, rhotyp, rhotot, vtot)
  !
  CALL sannp_exc(rho, rho_core, rhog_core, rhoexc)
  !
  CALL sannp_hartree(rho%of_g(:, 1), vhart)
  !
  DO ia = ia_start, ia_end
    !
    e = 0.0_DP
    q = 0.0_DP
    !
    CALL sannp_rhoa_single(ia, ALPHA, rhotyp, rhoatm, vatm)
    !
    ! ... Hirshfeld function
    !
    DO ir = 1, dfftp%nnr
      weight = 0.0_DP
      !
      rhotot0 = rhotot(ir)
      vtot0 = vtot(ir)
      !
      IF (rhotot0 >= RHO0) THEN
        rscale = 0.5_DP * qe_erfc(LRHO1 - LOG(rhotot0))
        vscale = 1.0_DP - rscale
      ELSE
        rscale = 0.0_DP
        vscale = 1.0_DP
      END IF
      !
      IF (rhotot0 >= RHO0) THEN
        weight = weight + rscale * rhoatm(ir) / rhotot0
      END IF
      !
      IF (vtot0 >= V0) THEN
        weight = weight + vscale * vatm(ir) / vtot0
      END IF
      !
      rhoatm(ir) = weight
    END DO
    !
    ! ... Kinetic energy
    !
    DO ir = 1, dfftp%nnr
      weight = rhoatm(ir)
      e = e + weight * rhokin(ir)
    END DO
    !
    ! ... Exchange Correlation energy
    !
    DO ir = 1, dfftp%nnr
      weight = rhoatm(ir)
      e = e + weight * rhoexc(ir)
    END DO
    !
    ! ... Local-potential and Hartree energy
    !
    DO ir = 1, dfftp%nnr
      weight = rhoatm(ir)
      e = e + weight * (vltot(ir) + 0.5_DP * vhart(ir)) * rho%of_r(ir, 1)
    END DO
    !
    ! ... Atomic Charge
    !
    DO ir = 1, dfftp%nnr
      weight = rhoatm(ir)
      q = q + weight * rho%of_r(ir, 1)
    END DO
    !
    e = e * omega / DBLE(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
    q = q * omega / DBLE(dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
    !
    CALL mp_sum(e, intra_bgrp_comm)
    CALL mp_sum(q, intra_bgrp_comm)
    !
    it = ityp(ia)
    q  = zv(it) - q
    !
    energy(ia) = energy(ia) + e
    charge(ia) = charge(ia) + q
    !
  END DO
  !
  DEALLOCATE(rhotyp)
  DEALLOCATE(rhotot)
  DEALLOCATE(rhoatm)
  DEALLOCATE(rhoexc)
  DEALLOCATE(vhart)
  DEALLOCATE(vtot)
  DEALLOCATE(vatm)
  !
END SUBROUTINE sannp_energy_lc
