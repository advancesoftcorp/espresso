!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE occmat_sum_band(becsum)
  !--------------------------------------------------------------------------
  !
  ! ... calculate Occupation Matrix
  !
  USE becmod,        ONLY : becp, allocate_bec_type, deallocate_bec_type
  USE buffers,       ONLY : get_buffer
  USE io_files,      ONLY : iunwfc, nwordwfc
  USE kinds,         ONLY : DP
  USE klist,         ONLY : nks, xk, ngk, igk_k
  USE mp,            ONLY : mp_sum
  USE mp_bands,      ONLY : inter_bgrp_comm, intra_bgrp_comm
  USE mp_pools,      ONLY : inter_pool_comm
  USE uspp,          ONLY : nkb, vkb
  USE wavefunctions, ONLY : evc
  USE wvfct,         ONLY : nbnd
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: becsum(:,:)
  !
  INTEGER :: ik, npw
  INTEGER :: ibnd_start, ibnd_end, this_bgrp_nbnd
  !
  CALL divide(inter_bgrp_comm, nbnd, ibnd_start, ibnd_end)
  !
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  becsum(:, :) = 0.0_DP
  !
  CALL allocate_bec_type(nkb, this_bgrp_nbnd, becp, intra_bgrp_comm)
  !
  DO ik = 1, nks
    !
    npw = ngk(ik)
    !
    IF (nks > 1) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
    !
    IF (nkb > 0) CALL init_us_2 (npw, igk_k(1, ik), xk(1, ik), vkb)
    !
    CALL occmat_sum_bec(ik, ibnd_start, ibnd_end, this_bgrp_nbnd, becsum)
    !
  END DO
  !
  CALL deallocate_bec_type(becp)
  !
  CALL mp_sum(becsum, inter_bgrp_comm)
  CALL mp_sum(becsum, inter_pool_comm)
  !
END SUBROUTINE occmat_sum_band
!
!--------------------------------------------------------------------------
SUBROUTINE occmat_sum_bec(ik, ibnd_start, ibnd_end, this_bgrp_nbnd, becsum)
  !----------------------------------------------------------------------------
  !
  ! ... calculate Occupation Matrix, for a k-point
  !
  USE becmod,        ONLY : becp, calbec
  USE control_flags, ONLY : gamma_only
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE kinds,         ONLY : DP
  USE klist,         ONLY : ngk
  USE uspp,          ONLY : nkb, vkb, indv_ijkb0
  USE uspp_param,    ONLY : upf, nh
  USE wavefunctions, ONLY : evc
  USE wvfct,         ONLY : wg
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: ik
  INTEGER,  INTENT(IN)  :: ibnd_start
  INTEGER,  INTENT(IN)  :: ibnd_end
  INTEGER,  INTENT(IN)  :: this_bgrp_nbnd
  REAL(DP), INTENT(OUT) :: becsum(:,:)
  !
  INTEGER :: npw, ikb
  INTEGER :: na, nt
  INTEGER :: ih, jh, ijh
  INTEGER :: ibnd, jbnd, ibnd_loc, nbnd_loc
  !
  REAL(DP),    ALLOCATABLE :: auxg  (:,:)
  REAL(DP),    ALLOCATABLE :: aux_gk(:,:)
  COMPLEX(DP), ALLOCATABLE :: auxk1 (:,:)
  COMPLEX(DP), ALLOCATABLE :: auxk2 (:,:)
  !
  npw = ngk(ik)
  !
  CALL calbec(npw, vkb, evc(:, ibnd_start:ibnd_end), becp)
  !
  DO nt = 1, ntyp
    !
    IF (.NOT. upf(nt)%tvanp) CYCLE
    !
    ! ... allocate memory
    !
    IF (gamma_only) THEN
      nbnd_loc = becp%nbnd_loc
      ALLOCATE(auxg (nbnd_loc, nh(nt)))
    ELSE
      ALLOCATE(auxk1(ibnd_start:ibnd_end, nh(nt)))
      ALLOCATE(auxk2(ibnd_start:ibnd_end, nh(nt)))
    END IF
    !
    ALLOCATE (aux_gk(nh(nt), nh(nt)))
    !
    ! ... In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
    !     run from index i=indv_ijkb0(na)+1 to i=indv_ijkb0(na)+nh(nt)
    !
    DO na = 1, nat
      !
      IF (ityp(na) /= np) CYCLE
      !
      ! ... sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
      !     copy into aux1, aux2 the needed data to perform a GEMM
      !
      IF (gamma_only) THEN
        !
!$omp parallel do default(shared), private(ih, ikb, ibnd, ibnd_loc)
        DO ih = 1, nh(nt)
          ikb = indv_ijkb0(na) + ih
          DO ibnd_loc = 1, nbnd_loc
            ibnd = (ibnd_start - 1) + ibnd_loc + becp%ibnd_begin - 1
            auxg(ibnd_loc, ih) = wg(ibnd, ik) * becp%r(ikb, ibnd_loc)
          END DO
        END DO
!$omp end parallel do
        !
        CALL DGEMM('N', 'N', nh(nt), nh(nt), nbnd_loc, &
                   1.0_DP, becp%r(indv_ijkb0(na) + 1, 1), nkb, auxg, nbnd_loc, &
                   0.0_DP, aux_gk, nh(nt))
        !
      ELSE
        !
!$omp parallel do default(shared), private(ih, ikb, ibnd, jbnd)
        DO ih = 1, nh(nt)
          ikb = indv_ijkb0(na) + ih
          DO jbnd = 1, this_bgrp_nbnd ! ibnd_start, ibnd_end
            ibnd = ibnd_start + jbnd - 1
            auxk1(ibnd, ih) = becp%k(ikb, jbnd)
            auxk2(ibnd, ih) = wg(ibnd, ik) * becp%k(ikb, jbnd)
          END DO
        END DO
!$omp end parallel do
        !
        ! only the real part is computed
        CALL DGEMM('C', 'N', nh(nt), nh(nt), 2 * this_bgrp_nbnd, &
                   1.0_DP, auxk1, 2 * this_bgrp_nbnd, auxk2, 2 * this_bgrp_nbnd, &
                   0.0_DP, aux_gk, nh(nt))
        !
      END IF
      !
      ! ... copy output from GEMM into desired format
      !
      ijh = 0
      !
      DO ih = 1, nh(nt)
        DO jh = ih, nh(nt)
          !
          ijh = ijh + 1
          !
          IF (jh == ih) THEN
            becsum(ijh, na) = becsum(ijh, na) + aux_gk(ih, jh)
          ELSE
            becsum(ijh, na) = becsum(ijh, na) + aux_gk(ih, jh) * 2.0_DP
          END IF
          !
        END DO
      END DO
      !
    END DO
    !
    ! ... deallocate memory
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
END SUBROUTINE occmat_sum_bec

