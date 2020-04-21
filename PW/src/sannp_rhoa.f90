!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_rhoa_total(alpha, rhogt, rhor, vr)
  !--------------------------------------------------------------------------
  !
  USE atom,           ONLY : rgrid, msh
  USE cell_base,      ONLY : tpiba, tpiba2, omega
  USE constants,      ONLY : eps8, e2, fpi
  USE control_flags,  ONLY : gamma_only
  USE fft_base,       ONLY : dfftp
  USE fft_interfaces, ONLY : invfft
  USE gvect,          ONLY : ngm, ngl, gstart, igtongl, gl, gg
  USE ions_base,      ONLY : ntyp => nsp
  USE kinds,          ONLY : DP
  USE uspp_param,     ONLY : upf
  USE vlocal,         ONLY : strf
  USE wavefunctions,  ONLY : psic
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: alpha
  REAL(DP), INTENT(OUT) :: rhogt(ngl, ntyp)
  REAL(DP), INTENT(OUT) :: rhor (dfftp%nnr)
  REAL(DP), INTENT(OUT) :: vr   (dfftp%nnr)
  !
  INTEGER  :: ndm
  INTEGER  :: it
  INTEGER  :: ir
  INTEGER  :: ig, igl
  REAL(DP) :: gx
  REAL(DP) :: ggt
  REAL(DP) :: aa
  REAL(DP) :: fac
  !
  REAL(DP),    ALLOCATABLE :: aux(:)
  COMPLEX(DP), ALLOCATABLE :: rhog(:)
  COMPLEX(DP), ALLOCATABLE :: vg(:)
  !
  ndm = MAXVAL(msh(1:ntyp))
  !
  ALLOCATE(aux(ndm))
  ALLOCATE(rhog(ngm))
  ALLOCATE(vg(ngm))
  !
  rhor(:) = 0.0_DP
  rhog(:) = (0.0_DP, 0.0_DP)
  !
  vr(:) = 0.0_DP
  vg(:) = (0.0_DP, 0.0_DP)
  !
  DO it = 1, ntyp
    !
    ! ... G == 0
    !
    IF (gstart == 2) then
      DO ir = 1, msh(it)
        aux (ir) = upf(it)%rho_at(ir)
      END DO
      !
      CALL simpson(msh(it), aux, rgrid(it)%rab, rhogt(1, it))
    END IF
    !
    ! ... G /= 0
    !
    DO igl = gstart, ngl
      gx = SQRT(gl(igl)) * tpiba
      DO ir = 1, msh(it)
        IF (rgrid(it)%r(ir) < eps8) then
          aux(ir) = upf(it)%rho_at(ir)
        ELSE
          aux(ir) = upf(it)%rho_at(ir) * &
                  & SIN(gx * rgrid(it)%r(ir)) / (rgrid(it)%r(ir) * gx)
        END IF
      END DO
      !
      CALL simpson(msh(it), aux, rgrid(it)%rab, rhogt(igl, it))
    END DO
    !
    ! ... map for G-space
    !
    DO ig = 1, ngm
      rhog(ig) = rhog(ig) + strf(ig, it) * rhogt(igtongl(ig), it) / omega
    END DO
    !
  END DO
  !
  fac = e2 * fpi
  aa  = alpha * alpha
  !
  IF (gstart == 2) THEN
    vg(1) = fac * rhog(1) * 0.25_DP * aa
  END IF
  !
  DO ig = gstart, ngm
    ggt = gg(ig) * tpiba2
    vg(ig) = fac * rhog(ig) * (1.0_DP - EXP(-0.25_DP * ggt * aa)) / ggt
  END DO
  !
  ! ... FFT to R-space
  !
  psic(:) = (0.0_DP, 0.0_DP)
  psic(dfftp%nl(:)) = rhog(:)
  IF (gamma_only) psic(dfftp%nlm(:)) = CONJG(rhog(:))
  !
  CALL invfft('Rho', psic, dfftp)
  !
  DO ir = 1, dfftp%nnr
    rhor(ir) = DBLE(psic(ir))
  END DO
  !
  psic(:) = (0.0_DP, 0.0_DP)
  psic(dfftp%nl(:)) = vg(:)
  IF (gamma_only) psic(dfftp%nlm(:)) = CONJG(vg(:))
  !
  CALL invfft('Rho', psic, dfftp)
  !
  DO ir = 1, dfftp%nnr
    vr(ir) = DBLE(psic(ir))
  END DO
  !
  DEALLOCATE(aux)
  DEALLOCATE(rhog)
  DEALLOCATE(vg)
  !
END SUBROUTINE sannp_rhoa_total
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_rhoa_single(ia, alpha, rhogt, rhor, vr)
  !--------------------------------------------------------------------------
  !
  USE cell_base,      ONLY : tpiba2, omega
  USE constants,      ONLY : e2, fpi
  USE control_flags,  ONLY : gamma_only
  USE fft_base,       ONLY : dfftp
  USE fft_interfaces, ONLY : invfft
  USE gvect,          ONLY : ngm, ngl, igtongl, gg, gstart, &
                           & mill, eigts1, eigts2, eigts3
  USE ions_base,      ONLY : ntyp => nsp, ityp
  USE kinds,          ONLY : DP
  USE wavefunctions,  ONLY : psic
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: ia
  REAL(DP), INTENT(IN)  :: alpha
  REAL(DP), INTENT(IN)  :: rhogt(ngl, ntyp)
  REAL(DP), INTENT(OUT) :: rhor (dfftp%nnr)
  REAL(DP), INTENT(OUT) :: vr   (dfftp%nnr)
  !
  INTEGER     :: it
  INTEGER     :: ir
  INTEGER     :: ig
  REAL(DP)    :: aa
  REAL(DP)    :: ggt
  REAL(DP)    :: fac
  COMPLEX(DP) :: sfac
  !
  COMPLEX(DP), ALLOCATABLE :: rhog(:)
  COMPLEX(DP), ALLOCATABLE :: vg(:)
  !
  ALLOCATE(rhog(ngm))
  ALLOCATE(vg(ngm))
  !
  rhor(:) = 0.0_DP
  rhog(:) = (0.0_DP, 0.0_DP)
  !
  vr(:) = 0.0_DP
  vg(:) = (0.0_DP, 0.0_DP)
  !
  ! ... map for G-space
  !
  it = ityp(ia)
  !
  DO ig = 1, ngm
    sfac = eigts1(mill(1, ig), ia) * &
         & eigts2(mill(2, ig), ia) * &
         & eigts3(mill(3, ig), ia)
    !
    rhog(ig) = rhog(ig) + sfac * rhogt(igtongl(ig), it) / omega
  END DO
  !
  fac = e2 * fpi
  aa  = alpha * alpha
  !
  IF (gstart == 2) THEN
    vg(1) = fac * rhog(1) * 0.25_DP * aa
  END IF
  !
  DO ig = gstart, ngm
    ggt = gg(ig) * tpiba2
    vg(ig) = fac * rhog(ig) * (1.0_DP - EXP(-0.25_DP * ggt * aa)) / ggt
  END DO
  !
  ! ... FFT to R-space
  !
  psic(:) = (0.0_DP, 0.0_DP)
  psic(dfftp%nl(:)) = rhog(:)
  IF (gamma_only) psic(dfftp%nlm(:)) = CONJG(rhog(:))
  !
  CALL invfft('Rho', psic, dfftp)
  !
  DO ir = 1, dfftp%nnr
    rhor(ir) = DBLE(psic(ir))
  END DO
  !
  psic(:) = (0.0_DP, 0.0_DP)
  psic(dfftp%nl(:)) = vg(:)
  IF (gamma_only) psic(dfftp%nlm(:)) = CONJG(vg(:))
  !
  CALL invfft('Rho', psic, dfftp)
  !
  DO ir = 1, dfftp%nnr
    vr(ir) = DBLE(psic(ir))
  END DO
  !
  DEALLOCATE(rhog)
  DEALLOCATE(vg)
  !
END SUBROUTINE sannp_rhoa_single
