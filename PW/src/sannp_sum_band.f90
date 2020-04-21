!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_sum_band(rhok, rhob)
  !--------------------------------------------------------------------------
  !
  ! ... calculate kinetic energy density and atomic density matrix
  !
  USE becmod,           ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE buffers,          ONLY : get_buffer
  USE cell_base,        ONLY : omega, tpiba
  USE control_flags,    ONLY : gamma_only
  USE fft_base,         ONLY : dfftp, dffts
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE gvect,            ONLY : g
  USE io_files,         ONLY : iunwfc, nwordwfc
  USE ions_base,        ONLY : nat
  USE kinds,            ONLY : DP
  USE klist,            ONLY : nks, ngk, xk, igk_k
  USE lsda_mod,         ONLY : lsda, nspin, current_spin, isk
  USE mp,               ONLY : mp_sum, mp_get_comm_null
  USE mp_bands,         ONLY : inter_bgrp_comm, intra_bgrp_comm
  USE mp_pools,         ONLY : inter_pool_comm
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp_param,       ONLY : nhm
  USE uspp,             ONLY : nkb, vkb
  USE wavefunctions,    ONLY : evc, psic, psic_nc
  USE wvfct,            ONLY : nbnd, npwx, wg
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: rhok(dfftp%nnr)
  REAL(DP), INTENT(OUT) :: rhob(nhm * (nhm + 1) / 2, nat, nspin) ! noncolin is not supported
  !
  INTEGER :: ibnd_start, ibnd_end, this_bgrp_nbnd
  !
  REAL(DP),    ALLOCATABLE :: kplusg(:)
  COMPLEX(DP), ALLOCATABLE :: rhog(:)
  !
  CALL divide(inter_bgrp_comm, nbnd, ibnd_start, ibnd_end)
  !
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  ALLOCATE(kplusg(npwx))
  ALLOCATE(rhog(dfftp%ngm))
  !
  CALL allocate_bec_type(nkb, nbnd, becp, intra_bgrp_comm)
  !
  rhok(:)       = 0.0_DP
  rhob(:, :, :) = 0.0_DP
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
  CALL mp_sum(rhok, inter_pool_comm)
  CALL mp_sum(rhok, inter_bgrp_comm)
  !
  CALL mp_sum(rhob, inter_pool_comm)
  CALL mp_sum(rhob, inter_bgrp_comm)
  !
  psic(1:dffts%nnr) = rhok(1:dffts%nnr)
  psic(dffts%nnr + 1:) = 0.0_DP
  CALL fwfft('Rho', psic, dffts)
  !
  rhog(1:dffts%ngm)  = psic(dffts%nl(1:dffts%ngm))
  rhog(dffts%ngm + 1:) = (0.0_DP, 0.0_DP)
  !
  psic(:) = (0.0_DP, 0.0_DP)
  psic(dfftp%nl(:)) = rhog(:)
  IF (gamma_only) psic(dfftp%nlm(:)) = CONJG(rhog(:))
  CALL invfft('Rho', psic, dfftp)
  !
  rhok(:) = psic(:)
  !
  DEALLOCATE(kplusg)
  DEALLOCATE(rhog)
  !
  CALL deallocate_bec_type(becp)
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
    REAL(DP) :: psir, psii
    !
    DO ik = 1, nks
      !
      IF (lsda) current_spin = isk(ik)
      !
      npw = ngk(ik)
      !
      IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
      !
      IF (nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(1, ik), vkb)
      !
      DO ibnd = ibnd_start, ibnd_end, 2
        !
        w1 = wg(ibnd, ik) / omega
        !
        IF (ibnd < ibnd_end) THEN
          !
          w2 = wg(ibnd + 1, ik) / omega
          !
        ELSE
          !
          w2 = w1
          !
        END IF
        !
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
            rhok(ir) = rhok(ir) + w1 * psir * psir + w2 * psii * psii
          END DO
          !
        END DO
        !
      END DO
      !
      CALL sannp_sum_bec(rhob, ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd)
      !
    END DO
    !
    IF (becp%comm /= mp_get_comm_null()) CALL mp_sum(rhob, becp%comm)
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
    REAL(DP) :: w
    REAL(DP) :: psir, psii
    !
    DO ik = 1, nks
      !
      IF (lsda) current_spin = isk(ik)
      !
      npw = ngk(ik)
      !
      IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
      !
      IF (nkb > 0) CALL init_us_2(npw, igk_k(1, ik), xk(1, ik), vkb)
      !
      DO ibnd = ibnd_start, ibnd_end
        !
        w = wg(ibnd, ik) / omega
        !
        IF (noncolin) THEN
          !
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
                rhok(ir) = rhok(ir) + w * (psir * psir + psii * psii)
              END DO
              !
            END DO
            !
          END DO
          !
        ELSE
          !
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
              rhok(ir) = rhok(ir) + w * (psir * psir + psii * psii)
            END DO
            !
          END DO
          !
        END IF
        !
      END DO
      !
      CALL sannp_sum_bec(rhob, ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd)
      !
    END DO
    !
    IF(becp%comm /= mp_get_comm_null()) CALL mp_sum(rhob, becp%comm)
    !
  END SUBROUTINE sum_band_k
  !
END SUBROUTINE sannp_sum_band
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_sum_bec(rhob, ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd)
  !--------------------------------------------------------------------------
  !
  ! ... calculate atomic density matrix
  !
  USE becmod,           ONLY : becp, calbec
  USE control_flags,    ONLY : gamma_only, tqr
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp
  USE kinds,            ONLY : DP
  USE klist,            ONLY : ngk
  USE lsda_mod,         ONLY : nspin
  USE mp,               ONLY : mp_sum
  USE mp_bands,         ONLY : nbgrp, inter_bgrp_comm
  USE noncollin_module, ONLY : noncolin
  USE realus,           ONLY : real_space, &
                             & invfft_orbital_gamma, calbec_rs_gamma, &
                             & invfft_orbital_k,     calbec_rs_k
  USE uspp,             ONLY : nkb, vkb, indv_ijkb0
  USE uspp_param,       ONLY : nh, nhm
  USE wavefunctions,    ONLY : evc
  USE wvfct,            ONLY : wg, current_k
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: rhob(nhm * (nhm + 1) / 2, nat, nspin) ! noncolin is not supported
  INTEGER,  INTENT(IN)    :: ik, current_spin
  INTEGER,  INTENT(IN)    :: ibnd_start, ibnd_end, this_bgrp_nbnd
  !
  INTEGER :: npw
  INTEGER :: ia, it
  INTEGER :: is, js
  INTEGER :: nht, ih, jh, ijh
  INTEGER :: ikb, jkb
  INTEGER :: ibnd, ibnd_loc, nbnd_loc
  !
  COMPLEX(DP), ALLOCATABLE :: auxk1(:,:)
  COMPLEX(DP), ALLOCATABLE :: auxk2(:,:)
  REAL(DP),    ALLOCATABLE :: auxg(:,:)
  REAL(DP),    ALLOCATABLE :: aux_gk(:,:)
  !
  IF (noncolin) CALL errore('sannp_sum_bec', &
  & 'the current version of SANNP does not support non-colinear', 1)
  !
  npw = ngk(ik)
  !
  IF (.NOT. real_space) THEN
    CALL calbec(npw, vkb, evc, becp)
    !
  ELSE
    IF (gamma_only) THEN
      DO ibnd = ibnd_start, ibnd_end, 2
        CALL invfft_orbital_gamma(evc, ibnd, ibnd_end)
        CALL calbec_rs_gamma(ibnd, ibnd_end, becp%r)
      END DO
      CALL mp_sum(becp%r, inter_bgrp_comm)
    ELSE
      current_k = ik
      becp%k = (0.0_DP, 0.0_DP)
      DO ibnd = ibnd_start, ibnd_end
        CALL invfft_orbital_k(evc, ibnd, ibnd_end)
        CALL calbec_rs_k(ibnd, ibnd_end)
      END DO
      CALL mp_sum(becp%k, inter_bgrp_comm)
    END IF
  END IF
  !
  !
  DO it = 1, ntyp
    !
    nht = nh(it)
    !
    IF (gamma_only) THEN
      nbnd_loc = becp%nbnd_loc
      ALLOCATE(auxg(nbnd_loc, nht))
    ELSE
      ALLOCATE(auxk1(ibnd_start:ibnd_end, nht))
      ALLOCATE(auxk2(ibnd_start:ibnd_end, nht))
    END IF
    !
    ALLOCATE(aux_gk(nht, nht))
    !
    DO ia = 1, nat
      !
      IF (ityp(ia) /= it) CYCLE
      !
      IF (gamma_only) THEN
        !
        DO ih = 1, nht
          ikb = indv_ijkb0(ia) + ih
          DO ibnd_loc = 1, nbnd_loc
            ibnd = ibnd_loc + becp%ibnd_begin - 1
            auxg(ibnd_loc, ih) = wg(ibnd, ik) * becp%r(ikb, ibnd_loc)
          END DO
        END DO
        !
        CALL DGEMM('N', 'N', nht, nht, nbnd_loc, &
                 & 1.0_DP / DBLE(nbgrp), becp%r(indv_ijkb0(ia) + 1, 1), &
                 & nkb, auxg, nbnd_loc, 0.0_DP, aux_gk, nht)
        !
      ELSE
        !
        DO ih = 1, nht
          ikb = indv_ijkb0(ia) + ih
          DO ibnd = ibnd_start, ibnd_end
            auxk1(ibnd, ih) = becp%k(ikb, ibnd)
            auxk2(ibnd, ih) = wg(ibnd, ik) * becp%k(ikb, ibnd)
          END DO
        END DO
        !
        CALL DGEMM('C', 'N', nht, nht, 2 * this_bgrp_nbnd, &
                 & 1.0_DP, auxk1, 2 * this_bgrp_nbnd, auxk2, &
                 & 2 * this_bgrp_nbnd, 0.0_DP, aux_gk, nht)
        !
      END IF
      !
      ijh = 0
      DO ih = 1, nht
        DO jh = ih, nht
          ijh = ijh + 1
          !
          IF (jh == ih) THEN
            rhob(ijh, ia, current_spin) = &
            & rhob(ijh, ia, current_spin) + aux_gk(ih, jh)
          ELSE
            rhob(ijh, ia, current_spin) = &
            & rhob(ijh, ia, current_spin) + aux_gk(ih, jh) * 2.0_DP
          END IF
          !
        END DO
      END DO
      !
    END DO
    !
    IF (gamma_only) THEN
      DEALLOCATE(auxg)
    ELSE
      DEALLOCATE(auxk1)
      DEALLOCATE(auxk2)
    END IF
    !
    DEALLOCATE(aux_gk)
    !
  END DO
  !
END SUBROUTINE sannp_sum_bec
