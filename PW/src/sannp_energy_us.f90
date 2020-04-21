!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_energy_us(energy, rhobec, ia_start, ia_end)
  !--------------------------------------------------------------------------
  !
  ! ... calculate non-local energy for single atom
  !
  USE ions_base,        ONLY : nat, ityp
  USE kinds,            ONLY : DP
  USE lsda_mod,         ONLY : nspin
  USE noncollin_module, ONLY : noncolin, nspin_mag
  USE uspp,             ONLY : dvan
  USE uspp_param,       ONLY : nhm, nh
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: energy(nat)
  REAL(DP), INTENT(IN)    :: rhobec(nhm * (nhm + 1) / 2, nat, nspin) ! noncolin is not supported
  INTEGER,  INTENT(IN)    :: ia_start
  INTEGER,  INTENT(IN)    :: ia_end
  !
  INTEGER  :: is
  INTEGER  :: ia, it
  INTEGER  :: nht
  INTEGER  :: ih, jh, ijh
  REAL(DP) :: e
  !
  IF (noncolin) CALL errore('sannp_energy_us', &
  & 'the current version of SANNP does not support non-colinear', 1)
  !
  DO ia = ia_start, ia_end
    !
    e = 0.0_DP
    !
    it  = ityp(ia)
    nht = nh(it)
    !
    DO is = 1, nspin_mag
      !
      ijh = 0
      DO ih = 1, nht
        DO jh = ih, nht
          ijh = ijh + 1
          !
          e = e + dvan(ih, jh, it) * rhobec(ijh, ia, is)
          !
        END DO
      END DO
      !
    END DO
    !
    energy(ia) = energy(ia) + e
    !
  END DO
  !
END SUBROUTINE sannp_energy_us

