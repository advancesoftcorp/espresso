!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE kin_funct
  !--------------------------------------------------------------------------
  !
  ! ... the module for Kinetic energy functionals
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  PUBLIC :: kin_tf
  PUBLIC :: kin_vw
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE kin_tf( length, fact, rho, et, vt )
    !--------------------------------------------------------------------------
    !
    ! ... Kinetic energy functional of Thomas-Fermi model.
    ! ... it supports only spin-unpolarised case (nspin=1).
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)    :: length       ! array length
    REAL(DP), INTENT(IN)    :: fact         ! scaling factor
    REAL(DP), INTENT(IN)    :: rho(length)  ! charge density
    REAL(DP), INTENT(INOUT) :: et (length)  ! energy density  (added in place)
    REAL(DP), INTENT(INOUT) :: vt (length)  ! potential       (added in place)
    !
    REAL(DP), PARAMETER :: pi   = 3.1415926535897932384626433832795_DP
    REAL(DP), PARAMETER :: C_tf = (3.0_DP / 10.0_DP) * ((3.0_DP * pi * pi) ** (2.0_DP / 3.0_DP))
    REAL(DP), PARAMETER :: rho_threshold = 1.0E-10_DP
    !
    INTEGER  :: ir
    REAL(DP) :: r, r23
    !
!$omp parallel do default(shared) private(ir, r, r23)
    DO ir = 1, length
       !
       r = rho(ir)
       !
       IF ( r > rho_threshold ) THEN
          !
          r23 = r ** (2.0_DP / 3.0_DP)
          !
          et(ir) = et(ir) + fact * C_tf * r * r23
          !
          vt(ir) = vt(ir) + fact * (5.0_DP / 3.0_DP) * C_tf * r23
          !
       END IF
       !
    END DO
!$omp end parallel do
    !
  END SUBROUTINE kin_tf
  !
  !--------------------------------------------------------------------------
  SUBROUTINE kin_vw( length, fact, rho, grho, et, v1t, v2t )
    !--------------------------------------------------------------------------
    !
    ! ... Kinetic energy functional of von Weizsacker model.
    ! ... it supports only spin-unpolarised case (nspin=1).
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN)    :: length        ! array length
    REAL(DP), INTENT(IN)    :: fact          ! scaling factor
    REAL(DP), INTENT(IN)    :: rho (length)  ! charge density
    REAL(DP), INTENT(IN)    :: grho(length)  ! gradient of rho: |\nabla rho|^2
    REAL(DP), INTENT(INOUT) :: et  (length)  ! energy density             (added in place)
    REAL(DP), INTENT(INOUT) :: v1t (length)  ! potential of density part  (added in place)
    REAL(DP), INTENT(INOUT) :: v2t (length)  ! potential of gradient part (added in place)
    !
    REAL(DP), PARAMETER :: C_vw = 1.0_DP / 8.0_DP
    REAL(DP), PARAMETER :: rho_threshold  = 1.0E-6_DP
    REAL(DP), PARAMETER :: grho_threshold = 1.0E-10_DP
    !
    INTEGER  :: ir
    REAL(DP) :: r, gr
    !
!$omp parallel do default(shared) private(ir, r, gr)
    DO ir = 1, length
       !
       r  = rho(ri)
       gr = grho(ri)
       !
       IF ( r > rho_threshold .AND. gr > grho_threshold ) THEN
          !
          et (ir) = et (ir) + fact * C_vw * gr / r
          !
          v1t(ir) = v1t(ir) - fact * C_vw * gr / (r * r)
          !
          v2t(ir) = v2t(ir) + fact * 2.0_DP * C_vw / r
          !
       END IF
       !
    END DO
!$omp end parallel do
    !
  END SUBROUTINE kin_vw
  !
END MODULE kin_funct
