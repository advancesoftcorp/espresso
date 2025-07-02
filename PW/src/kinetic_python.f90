!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE kinetic_python( rho, et, vt )
  !--------------------------------------------------------------------------
  !
  ! ... Calculate Kinetic energy functional, using external Python's code.
  ! ...
  ! ... [NOTE]
  ! ... This subroutine supports only spin-unpolarized case.
  ! ... MPI-parallelization is not available (serial only).
  ! ...
  ! ... T = \int tau(\text{rho}) dr, tau(\text{rho}) = \text{rho}\epsilon(\text{rho})\ ,
  ! ...
  ! ... where T is kinetic energy,
  ! ...       tau is kinetic energy density,
  ! ...       epsilon is kinetic energy weight.
  !
  USE fft_base,  ONLY : dfftp
  USE kinds,     ONLY : DP
  USE mp_images, ONLY : nproc_image
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)    :: rho (dfftp%nnr)  ! charge density (rho)
  REAL(DP), INTENT(INOUT) :: ekin(dfftp%nnr)  ! energy weight  (epsilon)
  REAL(DP), INTENT(INOUT) :: vkin(dfftp%nnr)  ! potential      (dT/drho)
  !
  INTEGER :: ir1, ir2, ir3, ir, jr
  INTEGER :: nr1, nr2, nr3
  INTEGER :: nr1x, nr2x, nr3x
  !
  REAL(DP), ALLOCATABLE :: rho_t (:)
  REAL(DP), ALLOCATABLE :: ekin_t(:)
  REAL(DP), ALLOCATABLE :: vkin_t(:)
  !
  REAL(DP), PARAMETER   :: rho_threshold = 1.0E-10_DP
  !
  IF (nproc_image > 1) CALL errore('kinetic_python', 'does not support MPI', 1)
  !
  nr1  = dfftp%nr1
  nr2  = dfftp%nr2
  nr3  = dfftp%nr3
  !
  nr1x = dfftp%nr1x
  nr2x = dfftp%nr2x
  nr3x = dfftp%nr3x
  !
  ALLOCATE(rho_t (nr1 * nr2 * nr3))
  ALLOCATE(ekin_t(nr1 * nr2 * nr3))
  ALLOCATE(vkin_t(nr1 * nr2 * nr3))
  !
  rho_t  = 0.0_DP
  ekin_t = 0.0_DP
  vkin_t = 0.0_DP
  !
  DO ir3 = 1, nr3
    !
    DO ir2 = 1, nr2
      !
      DO ir1 = 1, nr1
        !
        ir = ir1 + (ir2 - 1) * nr1x + (ir3 - 1) * nr1x * nr2x
        jr = ir3 + (ir2 - 1) * nr3  + (ir1 - 1) * nr3  * nr2
        !
        rho_t(jr) = rho(ir)
        !
      END DO
      !
    END DO
    !
  END DO
  !
  ! TODO
  ! TODO call C/C++ -> Python, to get ekin_t and vkin_t
  ! TODO
  !
  DO ir3 = 1, nr3
    !
    DO ir2 = 1, nr2
      !
      DO ir1 = 1, nr1
        !
        ir = ir1 + (ir2 - 1) * nr1x + (ir3 - 1) * nr1x * nr2x
        jr = ir3 + (ir2 - 1) * nr3  + (ir1 - 1) * nr3  * nr2
        !
        IF (rho(ir) > rho_threshold) THEN
          !
          ekin(ir) = ekin(ir) + ekin_t(jr) / rho(ir)
          !
        END IF
        !
        vkin(ir) = vkin(ir) + vkin_t(jr)
        !
      END DO
      !
    END DO
    !
  END DO
  !
  DEALLOCATE(rho_t)
  DEALLOCATE(ekin_t)
  DEALLOCATE(vkin_t)
  !
END SUBROUTINE kinetic_python
