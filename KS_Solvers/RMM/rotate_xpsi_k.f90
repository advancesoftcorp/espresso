!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_xpsi_k( h_psi, s_psi, overlap, &
                          npwx, npw, nstart, nbnd, npol, psi, evc, hevc, sevc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_xpsi for colinear, k-point calculations
  !
  USE rmm_param,     ONLY : DP
  USE mp_bands_util, ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE mp,            ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nstart)
  COMPLEX(DP), INTENT(OUT)   :: evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT) :: hevc(npwx*npol,nbnd), sevc(npwx*npol,nbnd)
    ! H|psi> and S|psi>
  REAL(DP), INTENT(OUT) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  INTEGER                  :: kdim, kdmx
  INTEGER                  :: n_start, n_end, my_n
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
  COMPLEX(DP), ALLOCATABLE :: tpsi(:,:), hpsi(:,:), spsi(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
  !
  EXTERNAL :: h_psi, s_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  CALL start_clock('rotxpsik')
  !
  ALLOCATE( tpsi( kdmx, nstart ) )
  ALLOCATE( hpsi( kdmx, nstart ) )
  IF ( overlap ) &
  ALLOCATE( spsi(kdmx, nstart ) )
  ALLOCATE( hc( nstart, nstart) )    
  ALLOCATE( sc( nstart, nstart) )    
  ALLOCATE( vc( nstart, nstart) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  CALL start_clock('rotxpsik:hpsi')
  !
  CALL h_psi( npwx, npw, nstart, psi, hpsi )
  !
  CALL stop_clock('rotxpsik:hpsi')
  !
  IF ( overlap ) THEN
     !
     CALL start_clock('rotxpsik:spsi')
     !
     CALL s_psi( npwx, npw, nstart, psi, spsi )
     !
     CALL stop_clock('rotxpsik:spsi')
     !
  END IF
  !
  CALL divide(inter_bgrp_comm, nstart, n_start, n_end)
  my_n = n_end - n_start + 1
  !
  CALL start_clock('rotxpsik:hc')
  !
  hc = (0.D0, 0.D0)
  !
  IF ( n_start <= n_end ) &
  CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), &
              psi, kdmx, hpsi(1,n_start), kdmx, (0.D0, 0.D0), hc(1,n_start), nstart )
  !
  CALL mp_sum( hc, inter_bgrp_comm )
  !
  CALL mp_sum( hc, intra_bgrp_comm )
  !
  CALL stop_clock('rotxpsik:hc')
  !
  CALL start_clock('rotxpsik:sc')
  !
  sc = (0.D0, 0.D0)
  !
  IF ( overlap ) THEN
     !
     IF ( n_start <= n_end ) &
     CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), &
                 psi, kdmx, spsi(1,n_start), kdmx, (0.D0, 0.D0), sc(1,n_start), nstart )
     !
  ELSE
     !
     IF ( n_start <= n_end ) &
     CALL ZGEMM( 'C', 'N', nstart, my_n, kdim, (1.D0, 0.D0), &
                 psi, kdmx, psi(1,n_start), kdmx, (0.D0, 0.D0), sc(1,n_start), nstart )
     !
  END IF
  !
  CALL mp_sum( sc, inter_bgrp_comm )
  !
  CALL mp_sum( sc, intra_bgrp_comm )
  !
  CALL stop_clock('rotxpsik:sc')
  !
  ! ... Diagonalize
  !
  CALL start_clock('rotxpsik:diag')
  !
  CALL cdiaghg( nstart, nbnd, hc, sc, nstart, en, vc )
  !
  e(:) = en(1:nbnd)
  !
  CALL stop_clock('rotxpsik:diag')
  !
  ! ...  update the basis set
  !
  CALL start_clock('rotxpsik:evc')
  !
  tpsi = psi
  !
  evc  = (0.D0, 0.D0)
  hevc = (0.D0, 0.D0)
  IF ( overlap ) &
  sevc = (0.D0, 0.D0)
  !
  IF ( n_start <= n_end ) THEN
     !
     CALL ZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 tpsi(1,n_start), kdmx, vc(n_start,1), nstart, (0.D0, 0.D0), evc,  kdmx )
     !
     CALL ZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 hpsi(1,n_start), kdmx, vc(n_start,1), nstart, (0.D0, 0.D0), hevc, kdmx )
     !
     IF ( overlap ) &
     CALL ZGEMM( 'N', 'N', kdim, nbnd, my_n, (1.D0, 0.D0), &
                 spsi(1,n_start), kdmx, vc(n_start,1), nstart, (0.D0, 0.D0), sevc, kdmx )
     !
  END IF
  !
  CALL mp_sum( evc,  inter_bgrp_comm )
  CALL mp_sum( hevc, inter_bgrp_comm )
  IF ( overlap ) &
  CALL mp_sum( sevc, inter_bgrp_comm )
  !
  CALL stop_clock('rotxpsik:evc')
  !     
  DEALLOCATE( en )
  DEALLOCATE( vc )
  DEALLOCATE( sc )
  DEALLOCATE( hc )
  IF ( overlap ) &
  DEALLOCATE( spsi )
  DEALLOCATE( hpsi )
  DEALLOCATE( tpsi )
  !
  CALL stop_clock('rotxpsik')
  !
  !CALL print_clock('rotxpsik')
  !CALL print_clock('rotxpsik:hpsi')
  !CALL print_clock('rotxpsik:spsi')
  !CALL print_clock('rotxpsik:hc')
  !CALL print_clock('rotxpsik:sc')
  !CALL print_clock('rotxpsik:diag')
  !CALL print_clock('rotxpsik:evc')
  !
  RETURN
  !
END SUBROUTINE rotate_xpsi_k
!
!
!----------------------------------------------------------------------------
SUBROUTINE protate_xpsi_k( h_psi, s_psi, overlap, &
                           npwx, npw, nstart, nbnd, npol, psi, evc, hevc, sevc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Parallel version of rotate_xpsi for colinear, k-point calculations
  ! ... Subroutine with distributed matrices, written by Carlo Cavazzoni
  !
  USE rmm_param,        ONLY : DP
  USE mp_bands_util,    ONLY : intra_bgrp_comm, inter_bgrp_comm, &
                               nbgrp, root_bgrp_id, my_bgrp_id
  USE mp_diag,          ONLY : ortho_comm, np_ortho, me_ortho, ortho_comm_id,&
                               leg_ortho, ortho_parent_comm, ortho_cntx, &
                               do_distr_diag_inside_bgrp
  USE descriptors,      ONLY : descla_init , la_descriptor
  USE parallel_toolkit, ONLY : zsqmher
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nstart)
  COMPLEX(DP), INTENT(OUT)   :: evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT) :: hevc(npwx*npol,nbnd), sevc(npwx*npol,nbnd)
    ! H|psi> and S|psi>
  REAL(DP), INTENT(OUT) :: e(nbnd)
    ! eigenvalues
  !
  ! ... local variables
  !
  INTEGER :: kdim, kdmx
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
  COMPLEX(DP), ALLOCATABLE :: tpsi(:,:), hpsi(:,:), spsi(:,:)
  REAL(DP),    ALLOCATABLE :: en(:)
  !
  TYPE(la_descriptor) :: desc
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  TYPE(la_descriptor), ALLOCATABLE :: desc_ip( :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  !
  EXTERNAL :: h_psi, s_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)

  CALL start_clock('protxpsik')
  !
  ALLOCATE( desc_ip( np_ortho(1), np_ortho(2) ) )
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
  !
  CALL desc_init( nstart, desc, desc_ip )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  ALLOCATE( tpsi( kdmx, nstart ) )
  ALLOCATE( hpsi( kdmx, nstart ) )
  IF ( overlap ) &
  ALLOCATE( spsi( kdmx, nstart ) )
  ALLOCATE( hc( nx, nx) )    
  ALLOCATE( sc( nx, nx) )    
  ALLOCATE( vc( nx, nx) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  CALL start_clock('protxpsik:hpsi')
  !
  CALL h_psi( npwx, npw, nstart, psi, hpsi )
  !
  CALL stop_clock('protxpsik:hpsi')
  !
  IF ( overlap ) THEN
     !
     CALL start_clock('protxpsik:spsi')
     !
     CALL s_psi( npwx, npw, nstart, psi, spsi )
     !
     CALL stop_clock('protxpsik:spsi')
     !
  END IF
  !
  CALL start_clock('protxpsik:hc')
  !
  CALL compute_distmat( hc, psi, hpsi )
  !
  CALL stop_clock('protxpsik:hc')
  !
  CALL start_clock('protxpsik:sc')
  !            
  IF ( overlap ) THEN
     !
     CALL compute_distmat( sc, psi, spsi )
     !
  ELSE
     !
     CALL compute_distmat( sc, psi, psi )
     !  
  END IF
  !
  CALL stop_clock('protxpsik:sc')
  !
  ! ... Diagonalize
  !
  CALL start_clock('protxpsik:diag')
  !
  IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pcdiaghg en and vc are the same across ortho_parent_comm
     ! only the first bgrp performs the diagonalization
     IF( my_bgrp_id == root_bgrp_id ) CALL pcdiaghg( nstart, hc, sc, nx, en, vc, desc )
     IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
       CALL mp_bcast( vc, root_bgrp_id, inter_bgrp_comm )
       CALL mp_bcast( en, root_bgrp_id, inter_bgrp_comm )
     ENDIF
  ELSE
     CALL pcdiaghg( nstart, hc, sc, nx, en, vc, desc )
  END IF
  !
  e(:) = en(1:nbnd)
  !
  CALL stop_clock('protxpsik:diag')
  !
  ! ...  update the basis set
  !
  CALL start_clock('protxpsik:evc')
  !
  tpsi = psi
  !
  CALL refresh_evc()
  !
  CALL stop_clock('protxpsik:evc')
  !
  DEALLOCATE( en )
  DEALLOCATE( vc )
  DEALLOCATE( sc )
  DEALLOCATE( hc )
  IF ( overlap ) &
  DEALLOCATE( spsi )
  DEALLOCATE( hpsi )
  DEALLOCATE( tpsi )
  !
  DEALLOCATE( desc_ip )
  DEALLOCATE( rank_ip )
  !
  CALL stop_clock('protxpsik')
  !
  !CALL print_clock('protxpsik')
  !CALL print_clock('protxpsik:hpsi')
  !CALL print_clock('protxpsik:spsi')
  !CALL print_clock('protxpsik:hc')
  !CALL print_clock('protxpsik:sc')
  !CALL print_clock('protxpsik:diag')
  !CALL print_clock('protxpsik:evc')
  !
  RETURN
  !
  !
CONTAINS
  !
  SUBROUTINE desc_init( nsiz, desc, desc_ip )
     !
     INTEGER, INTENT(IN)  :: nsiz
     TYPE(la_descriptor), INTENT(OUT) :: desc
     TYPE(la_descriptor), INTENT(OUT) :: desc_ip(:,:)
     INTEGER :: i, j, rank
     INTEGER :: coor_ip( 2 )
     !
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
     !
     nx = desc%nrcx
     !
     DO j = 0, desc%npc - 1
        DO i = 0, desc%npr - 1
           coor_ip( 1 ) = i
           coor_ip( 2 ) = j
           CALL descla_init( desc_ip(i+1,j+1), desc%n, desc%nx, &
                             np_ortho, coor_ip, ortho_comm, ortho_cntx, 1 )
           CALL GRID2D_RANK( 'R', desc%npr, desc%npc, i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( desc%active_node > 0 ) la_proc = .TRUE.
     !
     RETURN
  END SUBROUTINE desc_init
  !
  !
  SUBROUTINE compute_distmat( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm 
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = ( 0.0_DP, 0.0_DP )
     !
     DO ipc = 1, desc%npc !  loop on column procs 
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        DO ipr = 1, ipc ! desc%npr ! ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, kdim, ( 1.D0, 0.D0 ),  v(1,ir), kdmx, w(1,ic), kdmx, ( 0.D0, 0.D0 ), work, nx )

           ! accumulate result on dm of root proc.
           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     CALL zsqmher( nstart, dm, nx, desc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat


  SUBROUTINE refresh_evc( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           !
           beta = ( 0.D0, 0.D0 )

           DO ipr = 1, desc%npr
              !
              nr = desc_ip( ipr, ipc )%nr
              ir = desc_ip( ipr, ipc )%ir
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vc(:,1:nc), root, ortho_parent_comm )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             tpsi(1,ir), kdmx, vc, nx, beta, evc(1,ic),  kdmx )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             hpsi(1,ir), kdmx, vc, nx, beta, hevc(1,ic), kdmx )
                 !
                 IF ( overlap ) &
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             spsi(1,ir), kdmx, vc, nx, beta, sevc(1,ic), kdmx )
                 !
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             tpsi(1,ir), kdmx, vtmp, nx, beta, evc(1,ic),  kdmx )
                 !
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             hpsi(1,ir), kdmx, vtmp, nx, beta, hevc(1,ic), kdmx )
                 !
                 IF ( overlap ) &
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ( 1.D0, 0.D0 ),  &
                             spsi(1,ir), kdmx, vtmp, nx, beta, sevc(1,ic), kdmx )
                 !
              END IF
              ! 

              beta = ( 1.D0, 0.D0 )

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_evc

END SUBROUTINE protate_xpsi_k
