!
! Copyright (C) 2025 AdvanceSoft Corp.
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE add_local_perturb(perturb, aux)
  !----------------------------------------------------------------------------
  !
  ! ... add perturbation potential to local potential, in G-space
  !
  USE fft_base,  ONLY : dfftp
  USE gvect,     ONLY : ngm, igtongl, mill, eigts1, eigts2, eigts3
  USE kinds,     ONLY : DP
  USE io_global, ONLY : ionode, stdout
  USE ions_base, ONLY : nat, ityp, zv, atm
  USE mp,        ONLY : mp_sum
  USE mp_images, ONLY : intra_image_comm
  USE vlocal,    ONLY : vloc
  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN)    :: perturb
  COMPLEX(DP), INTENT(INOUT) :: aux(dfftp%nnr)
  !
  INTEGER  :: ia, it
  INTEGER  :: ig
  REAL(DP) :: fac
  REAL(DP) :: ztot
  REAL(DP) :: qtot
  !
  REAL(DP), ALLOCATABLE :: za(:)
  REAL(DP), ALLOCATABLE :: qa(:)
  !
  INTEGER  :: seed(2)
  INTEGER  :: clock
  REAL(DP) :: rand_value
  LOGICAL, SAVE :: rand_init = .FALSE.
  !
  IF (perturb <= 0.0_DP) THEN
    RETURN
  END IF
  !
  ALLOCATE(za(nat))
  ALLOCATE(qa(nat))
  za = 0.0_DP
  qa = 0.0_DP
  !
  IF (ionode) THEN
    !
    IF (.NOT. rand_init) THEN
      rand_init = .TRUE.
      CALL system_clock(count=clock)
      seed(1) = clock
      seed(2) = clock / 2
      CALL random_seed(put=seed)
    END IF
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      !
      IF (zv(it) > 0.0_DP) THEN
        CALL random_number(rand_value)
        fac = 2.0_DP * (rand_value - 0.5_DP) * kinetic_perturb / zv(it)
        fac = MAX(fac, -1.0_DP)
      ELSE
        fac = 0.0_DP
      END IF
      !
      qa(ia) = fac * zv(it)
      za(ia) = zv(it) + qa(ia)
      !
    END DO
    !
  END IF
  !
  CALL mp_sum(za, intra_image_comm)
  CALL mp_sum(qa, intra_image_comm)
  !
  IF (ionode) THEN
    WRITE(stdout, '()')
    WRITE(stdout, '(5X,"Perturbed Atomic Charge:")')
  END IF
  !
  ztot = 0.0_DP
  qtot = 0.0_DP
  !
  DO ia = 1, nat
    !
    it = ityp(ia)
    !
    ztot = ztot + za(ia)
    qtot = qtot + qa(ia)
    !
    IF (ionode) &
    WRITE(stdout, '(5X,I3,2X,A4,2F10.4)') ia, ADJUSTL(atm(it)) // '    ', za(ia), qa(ia)
    !
    IF (zv(it) > 0.0_DP) THEN
      fac = qa(ia) / zv(it)
    ELSE
      fac = 0.0_DP
    END IF
    !
    DO ig = 1, ngm
      !
      aux(dfftp%nl(ig)) = aux(dfftp%nl(ig)) &
                        + fac * vloc(igtongl(ig), it) &
                        * eigts1(mill(1, ig), ia) &
                        * eigts2(mill(2, ig), ia) &
                        * eigts3(mill(3, ig), ia)
      !
    END DO
    !
  END DO
  !
  IF (ionode) &
  WRITE(stdout, '(10X,A4,2F10.4)') "Sum ", ztot, qtot
  !
  DEALLOCATE(za)
  DEALLOCATE(qa)
  !
END SUBROUTINE add_local_perturb

