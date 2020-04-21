!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_hartree(rhog, vh)
  !--------------------------------------------------------------------------
  !
  ! ... calculate hartree potential
  !
  USE cell_base,      ONLY : tpiba2
  USE constants,      ONLY : fpi, e2
  USE control_flags,  ONLY : gamma_only
  USE fft_base,       ONLY : dfftp
  USE fft_interfaces, ONLY : invfft
  USE gvect,          ONLY : ngm, gg, gstart
  USE kinds,          ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: rhog(ngm)
  REAL(DP),    INTENT(OUT) :: vh(dfftp%nnr)
  !
  INTEGER  :: ig
  REAL(DP) :: fac
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  !
  ALLOCATE(aux(dfftp%nnr))
  !
  aux = (0.0_DP, 0.0_DP)
  !
  fac = e2 * fpi / tpiba2
  !
  DO ig = gstart, ngm
    !
    aux(dfftp%nl(ig)) = fac * rhog(ig) / gg(ig)
    !
  END DO
  !
  IF ( gamma_only ) THEN
    !
    DO ig = gstart, ngm
      !
      aux(dfftp%nlm(ig)) = CONJG(aux(dfftp%nl(ig)))
      !
    END DO
    !
  END IF
  !
  CALL invfft('Rho', aux, dfftp)
  !
  vh(:) = DBLE(aux(:))
  !
  DEALLOCATE(aux)
  !
END SUBROUTINE sannp_hartree
