//    This file is part of ELPA.
//
//    The ELPA library was originally created by the ELPA consortium,
//    consisting of the following organizations:
//
//    - Max Planck Computing and Data Facility (MPCDF), formerly known as
//      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
//    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
//      Informatik,
//    - Technische Universität München, Lehrstuhl für Informatik mit
//      Schwerpunkt Wissenschaftliches Rechnen ,
//    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
//    - Max-Plack-Institut für Mathematik in den Naturwissenschaften,
//      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
//      and
//    - IBM Deutschland GmbH
//
//    This particular source code file has been developed within the ELPA-AEO //
//    project, which has been a joint effort of
//
//    - Max Planck Computing and Data Facility (MPCDF), formerly known as
//      Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
//    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
//      Informatik,
//    - Technische Universität München, Lehrstuhl für Informatik mit
//      Schwerpunkt Wissenschaftliches Rechnen ,
//    - Technische Universität München, Lehrstuhl für Theoretische Chemie,
//    - Fritz-Haber-Institut, Berlin, Abt. Theorie,

//    More information can be found here:
//    http://elpa.mpcdf.mpg.de/ and
//    http://elpa-aeo.mpcdf.mpg.de
//
//    ELPA is free software: you can redistribute it and/or modify
//    it under the terms of the version 3 of the license of the
//    GNU Lesser General Public License as published by the Free
//    Software Foundation.
//
//    ELPA is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public License
//    along with ELPA.  If not, see <http://www.gnu.org/licenses/>
//
//    ELPA reflects a substantial effort on the part of the original
//    ELPA consortium, and we ask you to respect the spirit of the
//    license that we chose: i.e., please contribute any changes you
//    may have back to the original ELPA library distribution, and keep
//    any derivatives of ELPA under the same license that we chose for
//    the original distribution, the GNU Lesser General Public License.
//
// Author: Valeriy Manin (Bergische Universität Wuppertal)
// integreated into the ELPA library Pavel Kus, Andeas Marek (MPCDF)

#include "config-f90.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#ifdef HAVE_64BIT_INTEGER_MATH_SUPPORT
#define C_INT_TYPE_PTR long int*
#define C_INT_TYPE long int
#define BLAS_KIND c_int64_t
#else
#define C_INT_TYPE_PTR int*
#define C_INT_TYPE int
#define BLAS_KIND c_int
#endif
#ifdef HAVE_64BIT_INTEGER_MPI_SUPPORT
#define C_INT_MPI_TYPE_PTR long int*
#define C_INT_MPI_TYPE long int
#define MPI_KIND c_int64_t
#else
#define C_INT_MPI_TYPE_PTR int*
#define C_INT_MPI_TYPE int
#define MPI_KIND c_int
#endif



// most of the file is not compiled if not using MPI
#ifdef WITH_MPI
#include <mpi.h>


//***********************************************************************************************************

#define REALCASE 1
#define DOUBLE_PRECISION 1
#include "../general/precision_macros.h"
#include "cannon_forw_template.c"
#include "cannon_back_template.c"
#undef DOUBLE_PRECISION
#undef REALCASE

/*
!f> interface
!f>   subroutine cannons_reduction_d(A, U, local_rowsCast, local_colsCast, a_desc, Res, toStore, row_comm, col_comm) &
!f>                             bind(C, name="cannons_reduction_c_d")
!f>     use precision
!f>     real(c_double)                :: A(local_rowsCast, local_colsCast), U(local_rowsCast, local_colsCast)
!f>     real(c_double)                :: Res(local_rowsCast, local_colsCast)
!f>     integer(kind=BLAS_KIND)       :: a_desc(9)
!f>     integer(kind=c_int),value     :: local_rowsCast, local_colsCast
!f>     integer(kind=MPI_KIND),value  :: row_comm, col_comm, ToStore
!f>   end subroutine
!f> end interface

*/
void cannons_reduction_c_d(double* A, double* U, int local_rowsCast, int local_colsCast, C_INT_TYPE_PTR a_desc,
                           double *Res, C_INT_MPI_TYPE ToStore, C_INT_MPI_TYPE row_comm, C_INT_MPI_TYPE col_comm);

/*
!f> interface
!f>   subroutine cannons_triang_rectangular_d(U, B, local_rowsCast, local_colsCast, u_desc, b_desc, Res, row_comm, col_comm) &
!f>                             bind(C, name="cannons_triang_rectangular_c_d")
!f>     use precision
!f>     real(c_double)                :: U(local_rowsCast, local_colsCast), B(local_rowsCast, local_colsCast)
!f>     real(c_double)                :: Res(local_rowsCast, local_colsCast)
!f>     integer(kind=BLAS_KIND)       :: u_desc(9), b_desc(9)
!f>     integer(kind=c_int),value     :: local_rowsCast, local_colsCast
!f>     integer(kind=MPI_KIND),value  :: row_comm, col_comm
!f>   end subroutine
!f> end interface
*/
void cannons_triang_rectangular_c_d(double* U, double* B, int local_rowsCast, int local_colsCast,
                                    C_INT_TYPE_PTR u_desc, C_INT_TYPE_PTR b_desc, double *Res, C_INT_MPI_TYPE row_comm, C_INT_MPI_TYPE col_comm);

//***********************************************************************************************************

#define REALCASE 1
#define SINGLE_PRECISION 1
#include "../general/precision_macros.h"
#include "cannon_forw_template.c"
#include "cannon_back_template.c"
#undef SINGLE_PRECISION
#undef REALCASE

/*
!f> interface
!f>   subroutine cannons_reduction_f(A, U, local_rowsCast, local_colsCast, a_desc, Res, toStore, row_comm, col_comm) &
!f>                             bind(C, name="cannons_reduction_c_f")
!f>     use precision
!f>     real(c_float)                 :: A(local_rowsCast, local_colsCast), U(local_rowsCast, local_colsCast)
!f>     real(c_float)                 :: Res(local_rowsCast, local_colsCast)
!f>     !type(c_ptr), value           :: A, U, Res
!f>     integer(kind=BLAS_KIND)       :: a_desc(9)
!f>     integer(kind=c_int),value     :: local_rowsCast, local_colsCast
!f>     integer(kind=MPI_KIND),value  :: row_comm, col_comm, ToStore
!f>   end subroutine
!f> end interface
*/
void cannons_reduction_c_f(float* A, float* U, int local_rowsCast, int local_colsCast, C_INT_TYPE_PTR a_desc,
                           float *Res, C_INT_MPI_TYPE ToStore, C_INT_MPI_TYPE row_comm, C_INT_MPI_TYPE col_comm);

/*
!f> interface
!f>   subroutine cannons_triang_rectangular_f(U, B, local_rowsCast, local_colsCast, u_desc, b_desc, Res, row_comm, col_comm) &
!f>                             bind(C, name="cannons_triang_rectangular_c_f")
!f>     use precision
!f>     real(c_float)                        :: U(local_rowsCast, local_colsCast), B(local_rowsCast, local_colsCast)
!f>     real(c_float)                        :: Res(local_rowsCast, local_colsCast)
!f>     integer(kind=BLAS_KIND)              :: u_desc(9), b_desc(9)
!f>     integer(kind=c_int),value            :: local_rowsCast, local_colsCast
!f>     integer(kind=MPI_KIND),value         :: row_comm, col_comm
!f>   end subroutine
!f> end interface
*/
void cannons_triang_rectangular_c_f(float* U, float* B, int local_rowsCast, int local_colsCast,
                                    C_INT_TYPE_PTR u_desc, C_INT_TYPE_PTR b_desc, float *Res, C_INT_MPI_TYPE row_comm, C_INT_MPI_TYPE col_comm);

//***********************************************************************************************************

#define COMPLEXCASE 1
#define DOUBLE_PRECISION 1
#include "../general/precision_macros.h"
#include "cannon_forw_template.c"
#include "cannon_back_template.c"
#undef DOUBLE_PRECISION
#undef COMPLEXCASE

/*
!f> interface
!f>   subroutine cannons_reduction_dc(A, U, local_rowsCast, local_colsCast, a_desc, Res, toStore, row_comm, col_comm) &
!f>                             bind(C, name="cannons_reduction_c_dc")
!f>     use precision
!f>     complex(c_double)             :: A(local_rowsCast, local_colsCast), U(local_rowsCast, local_colsCast)
!f>     complex(c_double)             :: Res(local_rowsCast, local_colsCast)
!f>     integer(kind=BLAS_KIND)       :: a_desc(9)
!f>     integer(kind=c_int),value     :: local_rowsCast, local_colsCast
!f>     integer(kind=MPI_KIND),value  :: row_comm, col_comm, ToStore
!f>   end subroutine
!f> end interface
*/
void cannons_reduction_c_dc(double complex* A, double complex* U, int local_rowsCast, int local_colsCasr, C_INT_TYPE_PTR a_desc,
                            double complex *Res, C_INT_MPI_TYPE ToStore, C_INT_MPI_TYPE row_comm, C_INT_MPI_TYPE col_comm);

/*
!f> interface
!f>   subroutine cannons_triang_rectangular_dc(U, B, local_rowsCast, local_colsCast, u_desc, b_desc, Res, row_comm, col_comm) &
!f>                             bind(C, name="cannons_triang_rectangular_c_dc")
!f>     use precision
!f>     complex(c_double)             :: U(local_rowsCast, local_colsCast), B(local_rowsCast, local_colsCast)
!f>     complex(c_double)             :: Res(local_rowsCast, local_colsCast)
!f>     integer(kind=BLAS_KIND)       :: u_desc(9), b_desc(9)
!f>     integer(kind=c_int),value     :: local_rowsCast, local_colsCast
!f>     integer(kind=MPI_KIND),value  :: row_comm, col_comm
!f>   end subroutine
!f> end interface
*/
void cannons_triang_rectangular_c_dc(double complex* U, double complex* B, int local_rowsCast, int local_colsCast,
                                    C_INT_TYPE_PTR u_desc, C_INT_TYPE_PTR b_desc, double complex *Res, C_INT_MPI_TYPE row_comm, C_INT_MPI_TYPE col_comm);
//***********************************************************************************************************

#define COMPLEXCASE 1
#define SINGLE_PRECISION 1
#include "../general/precision_macros.h"
#include "cannon_forw_template.c"
#include "cannon_back_template.c"
#undef SINGLE_PRECISION
#undef COMPLEXCASE

/*
!f> interface
!f>   subroutine cannons_reduction_fc(A, U, local_rowsCast, local_colsCast, a_desc, Res, toStore, row_comm, col_comm) &
!f>                             bind(C, name="cannons_reduction_c_fc")
!f>     use precision
!f>     complex(c_float)              :: A(local_rowsCast, local_colsCast), U(local_rowsCast, local_colsCast)
!f>     complex(c_float)              :: Res(local_rowsCast, local_colsCast)
!f>     !type(c_ptr), value           :: A, U, Res
!f>     integer(kind=BLAS_KIND)       :: a_desc(9)
!f>     integer(kind=c_int),value     :: local_rowsCast, local_colsCast
!f>     integer(kind=MPI_KIND),value  :: row_comm, col_comm, ToStore
!f>   end subroutine
!f> end interface
*/

void cannons_reduction_c_fc(float complex* A, float complex* U, int local_rowsCast, int local_colsCast, C_INT_TYPE_PTR a_desc,
                         float complex *Res, C_INT_MPI_TYPE ToStore, C_INT_MPI_TYPE row_comm, C_INT_MPI_TYPE col_comm);

/*
!f> interface
!f>   subroutine cannons_triang_rectangular_fc(U, B, local_rowsCast, local_colsCast, u_desc, b_desc, Res, row_comm, col_comm) &
!f>                             bind(C, name="cannons_triang_rectangular_c_fc")
!f>     use precision
!f>     complex(c_float)              :: U(local_rowsCast, local_colsCast), B(local_rowsCast, local_colsCast)
!f>     complex(c_float)              :: Res(local_rowsCast, local_colsCast)
!f>     integer(kind=BLAS_KIND)       :: u_desc(9), b_desc(9)
!f>     integer(kind=c_int),value     :: local_rowsCast, local_colsCast
!f>     integer(kind=MPI_KIND),value  :: row_comm, col_comm
!f>   end subroutine
!f> end interface
*/
void cannons_triang_rectangular_c_fc(float complex* U, float complex* B, int local_rowsCast, int local_colsCast,
                                    C_INT_TYPE_PTR u_desc, C_INT_TYPE_PTR b_desc, float complex *Res, C_INT_MPI_TYPE row_comm, C_INT_MPI_TYPE col_comm);
#endif
