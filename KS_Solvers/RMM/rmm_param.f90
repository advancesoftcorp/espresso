!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE rmm_param

  USE parallel_include

  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, PARAMETER :: stdout = 6    ! unit connected to standard output

  REAL(DP), PARAMETER :: eps14 = 1.0E-14_DP
  REAL(DP), PARAMETER :: eps16 = 1.0E-16_DP

END MODULE rmm_param
