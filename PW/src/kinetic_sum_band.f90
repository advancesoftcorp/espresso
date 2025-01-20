!
! Copyright (C) 2024 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE kinetic_sum_band(rhor, tauG, tauL, dtdr, by_veff)
  !--------------------------------------------------------------------------
  !
  ! ... calculate Kinetic Energy Density
  !
  USE buffers,          ONLY : get_buffer
  USE cell_base,        ONLY : omega, tpiba, tpiba2
  USE control_flags,    ONLY : gamma_only
  USE ener,             ONLY : ef
  USE fft_base,         ONLY : dffts
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE gvect,            ONLY : g, gg, gstart
  USE io_files,         ONLY : iunwfc, nwordwfc
  USE kinds,            ONLY : DP
  USE klist,            ONLY : nks, ngk, xk, igk_k
  USE mp,               ONLY : mp_sum
  USE mp_bands,         ONLY : inter_bgrp_comm
  USE mp_pools,         ONLY : inter_pool_comm
  USE noncollin_module, ONLY : noncolin, npol
  USE scf,              ONLY : vrs
  USE wavefunctions,    ONLY : evc, psic, psic_nc
  USE wvfct,            ONLY : nbnd, npwx, wg, et
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: rhor(dffts%nnr)
  REAL(DP), INTENT(OUT) :: tauG(dffts%nnr)
  REAL(DP), INTENT(OUT) :: tauL(dffts%nnr)
  REAL(DP), INTENT(OUT) :: dtdr(dffts%nnr)
  LOGICAL,  INTENT(IN)  :: by_veff
  !
  INTEGER  :: ig, ir
  INTEGER  :: ibnd_start, ibnd_end
  REAL(DP) :: fac
  REAL(DP) :: rho0
  !
  REAL(DP),    ALLOCATABLE :: kplusg(:)
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  !
  REAL(DP),    PARAMETER   :: rho_min = 1.0E-16_DP
  !
  CALL divide(inter_bgrp_comm, nbnd, ibnd_start, ibnd_end)
  !
  ALLOCATE(kplusg(npwx))
  ALLOCATE(aux(dffts%nnr))
  !
  rhor(:) = 0.0_DP
  tauG(:) = 0.0_DP
  tauL(:) = 0.0_DP
  dtdr(:) = 0.0_DP
  !
  ! ... calculate rhor = |psi|^2
  ! ... and       tauG = |grad psi|^2
  ! ... and       dtdr = (ef-e)*|psi|^2
  !
  IF (gamma_only) THEN
    !
    CALL sum_band_gamma()
    !
  ELSE
    !
    CALL sum_band_k()
    !
  END IF
  !
  CALL mp_sum(rhor, inter_pool_comm)
  CALL mp_sum(rhor, inter_bgrp_comm)
  !
  CALL mp_sum(tauG, inter_pool_comm)
  CALL mp_sum(tauG, inter_bgrp_comm)
  !
  CALL mp_sum(dtdr, inter_pool_comm)
  CALL mp_sum(dtdr, inter_bgrp_comm)
  !
  ! ... tauL = -psi Lap psi = |grad psi|^2 - (1/2)*Lap rho
  !
  aux(:) = rhor(:)
  !
  CALL fwfft('Rho', aux, dffts)
  !
  IF (gstart > 1) THEN
    !
    aux(dffts%nl(1)) = (0.0_DP, 0.0_DP)
    !
  END IF
  !
  fac = -0.5_DP * tpiba2
  !
  DO ig = gstart, dffts%ngm
    !
    aux(dffts%nl(ig)) = fac * gg(ig) * aux(dffts%nl(ig))
    !
  END DO
  !
  IF (gamma_only) THEN
    !
    DO ig = gstart, dffts%ngm
      !
      aux(dffts%nlm(ig)) = CONJG(aux(dffts%nl(ig)))
      !
    END DO
    !
  END IF
  !
  CALL invfft('Rho', aux, dffts)
  !
  tauL(:) = tauG(:) - DBLE(aux(:))
  !
  IF (by_veff) THEN
    !
    ! ... dtdr = ef - veff
    !
    dtdr(1:dffts%nnr) = vrs(1:dffts%nnr, 1) - ef
    !
  ELSE
    !
    ! ... dtdr = (-psi Lap psi + (ef-e)*|psi|^2) / rho
    !
    DO ir = 1, dffts%nnr
      !
      rho0 = rhor(ir)
      !
      IF (rho0 > rho_min) THEN
        !
        dtdr(ir) = (tauL(ir) + dtdr(ir)) / rho0
        !
      ELSE
        !
        dtdr(ir) = 0.0_DP
        !
      END IF
      !
    END DO
    !
  END IF
  !
  DEALLOCATE(kplusg)
  DEALLOCATE(aux)
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE sum_band_gamma()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER  :: ik
    INTEGER  :: npw
    INTEGER  :: ibnd
    INTEGER  :: ix
    INTEGER  :: ir
    REAL(DP) :: w1, w2
    REAL(DP) :: de1, de2
    REAL(DP) :: psir, psii
    REAL(DP) :: rho1, rho2
    !
    DO ik = 1, nks
      !
      npw = ngk(ik)
      !
      IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
      !
      DO ibnd = ibnd_start, ibnd_end, 2
        !
        w1 = wg(ibnd, ik) / omega
        de1 = ef - et(ibnd, ik)
        !
        IF (ibnd < ibnd_end) THEN
          !
          w2 = wg(ibnd + 1, ik) / omega
          de2 = ef - et(ibnd + 1, ik)
          !
        ELSE
          !
          w2 = w1
          de2 = de1
          !
        END IF
        !
        ! ... |grad psi|^2
        DO ix = 1, 3
          !
          psic(:) = (0.0_DP, 0.0_DP)
          !
          kplusg(1:npw) = (xk(ix, ik) + g(ix, 1:npw)) * tpiba
          !
          IF (ibnd < ibnd_end) THEN
            !
            psic(dffts%nl(1:npw))  = CMPLX(0.0_DP,  kplusg(1:npw), kind=DP) * &
                                   & (evc(1:npw, ibnd) + &
                                   & (0.0_DP, 1.0_DP) * evc(1:npw, ibnd + 1))
            !
            psic(dffts%nlm(1:npw)) = CMPLX(0.0_DP, -kplusg(1:npw), kind=DP) * &
                                   & CONJG(evc(1:npw, ibnd) - &
                                   & (0.0_DP, 1.0_DP) * evc(1:npw, ibnd + 1))
            !
          ELSE
            !
            psic(dffts%nl(1:npw))  = CMPLX(0.0_DP,  kplusg(1:npw), kind=DP) * &
                                   & evc(1:npw, ibnd)
            !
            psic(dffts%nlm(1:npw)) = CMPLX(0.0_DP, -kplusg(1:npw), kind=DP) * &
                                   & CONJG(evc(1:npw, ibnd))
            !
          END IF
          !
          CALL invfft('Wave', psic, dffts)
          !
          DO ir = 1, dffts%nnr
            psir =  DBLE(psic(ir))
            psii = AIMAG(psic(ir))
            !
            tauG(ir) = tauG(ir) + w1 * psir * psir + w2 * psii * psii
          END DO
          !
        END DO
        !
        ! ... |psi|^2 and (ef-e)*|psi|^2
        psic(:) = (0.0_DP, 0.0_DP)
        !
        IF (ibnd < ibnd_end) THEN
          !
          psic(dffts%nl(1:npw))  =       evc(1:npw, ibnd) + (0.0_DP, 1.0_DP) * evc(1:npw, ibnd + 1)
          psic(dffts%nlm(1:npw)) = CONJG(evc(1:npw, ibnd) - (0.0_DP, 1.0_DP) * evc(1:npw, ibnd + 1))
          !
        ELSE
          !
          psic(dffts%nl(1:npw))  =       evc(1:npw, ibnd)
          psic(dffts%nlm(1:npw)) = CONJG(evc(1:npw, ibnd))
          !
        END IF
        !
        CALL invfft('Wave', psic, dffts)
        !
        DO ir = 1, dffts%nnr
          psir =  DBLE(psic(ir))
          psii = AIMAG(psic(ir))
          !
          rho1 = w1 * psir * psir
          rho2 = w2 * psii * psii
          !
          rhor(ir) = rhor(ir) + rho1 + rho2
          dtdr(ir) = dtdr(ir) + de1 * rho1 + de2 * rho2
        END DO
        !
      END DO
      !
    END DO
    !
  END SUBROUTINE sum_band_gamma
  !
  !--------------------------------------------------------------------------
  SUBROUTINE sum_band_k()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER  :: ik
    INTEGER  :: npw
    INTEGER  :: ibnd
    INTEGER  :: ipol
    INTEGER  :: ix
    INTEGER  :: ir
    INTEGER  :: ig
    REAL(DP) :: w, de
    REAL(DP) :: psir, psii
    REAL(DP) :: rho0
    !
    DO ik = 1, nks
      !
      npw = ngk(ik)
      !
      IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
      !
      DO ibnd = ibnd_start, ibnd_end
        !
        w = wg(ibnd, ik) / omega
        de = ef - et(ibnd, ik)
        !
        IF (noncolin) THEN
          !
          ! ... |grad psi|^2
          DO ix = 1, 3
            !
            psic_nc(:, :) = (0.0_DP, 0.0_DP)
            !
            kplusg(1:npw) = (xk(ix, ik) + g(ix, igk_k(1:npw, ik))) * tpiba
            !
            DO ig = 1, npw
              psic_nc(dffts%nl(igk_k(ig, ik)), 1) = CMPLX(0.0_DP, kplusg(ig), kind=DP) * &
                                                  & evc(ig       , ibnd)
            END DO
            !
            DO ig = 1, npw
              psic_nc(dffts%nl(igk_k(ig, ik)), 2) = CMPLX(0.0_DP, kplusg(ig), kind=DP) * &
                                                  & evc(ig + npwx, ibnd)
            END DO
            !
            CALL invfft('Wave', psic_nc(:, 1), dffts)
            CALL invfft('Wave', psic_nc(:, 2), dffts)
            !
            DO ipol = 1, npol
              !
              DO ir = 1, dffts%nnr
                psir =  DBLE(psic_nc(ir, ipol))
                psii = AIMAG(psic_nc(ir, ipol))
                !
                tauG(ir) = tauG(ir) + w * (psir * psir + psii * psii)
              END DO
              !
            END DO
            !
          END DO
          !
          ! ... |psi|^2 and (ef-e)*|psi|^2
          psic_nc(:, :) = (0.0_DP, 0.0_DP)
          !
          DO ig = 1, npw
            psic_nc(dffts%nl(igk_k(ig, ik)), 1) = evc(ig, ibnd)
          END DO
          !
          DO ig = 1, npw
            psic_nc(dffts%nl(igk_k(ig, ik)), 2) = evc(ig + npwx, ibnd)
          END DO
          !
          CALL invfft('Wave', psic_nc(:, 1), dffts)
          CALL invfft('Wave', psic_nc(:, 2), dffts)
          !
          DO ipol = 1, npol
            !
            DO ir = 1, dffts%nnr
              psir =  DBLE(psic_nc(ir, ipol))
              psii = AIMAG(psic_nc(ir, ipol))
              !
              rho0 = w * (psir * psir + psii * psii)
              !
              rhor(ir) = rhor(ir) + rho0
              dtdr(ir) = dtdr(ir) + de * rho0
            END DO
            !
          END DO
          !
        ELSE
          !
          ! ... |grad psi|^2
          DO ix = 1, 3
            !
            psic = (0.0_DP, 0.0_DP)
            !
            kplusg(1:npw) = (xk(ix, ik) + g(ix, igk_k(1:npw, ik))) * tpiba
            !
            DO ig = 1, npw
              psic(dffts%nl(igk_k(ig, ik))) = CMPLX(0.0_DP, kplusg(ig), kind=DP) * &
                                            & evc(ig, ibnd)
            END DO
            !
            CALL invfft('Wave', psic, dffts)
            !
            DO ir = 1, dffts%nnr
              psir =  DBLE(psic(ir))
              psii = AIMAG(psic(ir))
              !
              tauG(ir) = tauG(ir) + w * (psir * psir + psii * psii)
            END DO
            !
          END DO
          !
          ! ... |psi|^2 and (ef-e)*|psi|^2
          psic = (0.0_DP, 0.0_DP)
          !
          DO ig = 1, npw
            psic(dffts%nl(igk_k(ig, ik))) = evc(ig, ibnd)
          END DO
          !
          CALL invfft('Wave', psic, dffts)
          !
          DO ir = 1, dffts%nnr
            psir =  DBLE(psic(ir))
            psii = AIMAG(psic(ir))
            !
            rho0 = w * (psir * psir + psii * psii)
            !
            rhor(ir) = rhor(ir) + rho0
            dtdr(ir) = dtdr(ir) + de * rho0
          END DO
          !
        END IF
        !
      END DO
      !
    END DO
    !
  END SUBROUTINE sum_band_k
  !
END SUBROUTINE kinetic_sum_band
