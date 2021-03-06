!
! Copyright (C) 2019 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE sannp_energy_paw(energy, ia_start, ia_end)
  !--------------------------------------------------------------------------
  !
  ! ... calculate PAW energy for single atom
  !
  USE ions_base, ONLY : nat
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: energy(nat)
  INTEGER,  INTENT(IN)    :: ia_start
  INTEGER,  INTENT(IN)    :: ia_end
  !
  !CALL errore('sannp_energy_paw', 'the current version of SANNP does not support PAW', 1)
  !
END SUBROUTINE sannp_energy_paw
