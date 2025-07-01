!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE nonloc_energy(energy, becsum)
  !--------------------------------------------------------------------------
  !
  ! ... calculate Non-Local energy for each atom
  !
  USE ions_base,  ONLY : nat, ityp
  USE kinds,      ONLY : DP
  USE uspp,       ONLY : dvan
  USE uspp_param, ONLY : nhm, nh
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: energy(nat)
  REAL(DP), INTENT(IN)  :: becsum(nhm, nhm, nat)
  !
  INTEGER  :: ia, it
  INTEGER  :: nht
  INTEGER  :: ih, jh
  REAL(DP) :: e
  !
  energy = 0.0_DP
  !
  DO ia = 1, nat
    !
    e = 0.0_DP
    !
    it  = ityp(ia)
    nht = nh(it)
    !
    DO ih = 1, nht
      !
      DO jh = 1, nht
        !
        e = e + dvan(jh, ih, it) * becsum(jh, ih, ia)
        !
      END DO
      !
    END DO
    !
    energy(ia) = e
    !
  END DO
  !
END SUBROUTINE nonloc_energy
