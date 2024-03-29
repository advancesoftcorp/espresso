//    This file is part of ELPA.
//
//    The ELPA library was originally created by the ELPA consortium,
//    consisting of the following organizations:
//
//    - Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG),
//    - Bergische Universität Wuppertal, Lehrstuhl für angewandte
//      Informatik,
//    - Technische Universität München, Lehrstuhl für Informatik mit
//      Schwerpunkt Wissenschaftliches Rechnen ,
//    - Fritz-Haber-Institut, Berlin, Abt. Theorie,
//    - Max-Plack-Institut für Mathematik in den Naturwissenschaftrn,
//      Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition,
//      and
//    - IBM Deutschland GmbH
//
//    This particular source code file contains additions, changes and
//    enhancements authored by Intel Corporation which is not part of
//    the ELPA consortium.
//
//    More information can be found here:
//    http://elpa.rzg.mpg.de/
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
//
// --------------------------------------------------------------------------------------------------
//
// This file contains the compute intensive kernels for the Householder transformations.
// It should be compiled with the highest possible optimization level.
//
// On Intel Nehalem or Intel Westmere or AMD Magny Cours use -O3 -msse3
// On Intel Sandy Bridge use -O3 -mavx
//
// Copyright of the original code rests with the authors inside the ELPA
// consortium. The copyright of any additional modifications shall rest
// with their original authors, but shall adhere to the licensing terms
// distributed along with the original code in the file "COPYING".
//
// Author: Alexander Heinecke (alexander.heinecke@mytum.de)
// Adapted for building a shared-library by Andreas Marek, RZG (andreas.marek@rzg.mpg.de)
// --------------------------------------------------------------------------------------------------

#include <complex>
#include <x86intrin.h>

#define __forceinline __attribute__((always_inline))

#ifdef __USE_AVX128__
#undef __AVX__
#endif

#ifdef __FMA4__
#define __ELPA_USE_FMA__
#define _mm256_FMADDSUB_pd(a,b,c) _mm256_maddsub_pd(a,b,c)
#define _mm256_FMSUBADD_pd(a,b,c) _mm256_msubadd_pd(a,b,c)
#endif

#ifdef __AVX2__
#define __ELPA_USE_FMA__
#define _mm256_FMADDSUB_pd(a,b,c) _mm256_fmaddsub_pd(a,b,c)
#define _mm256_FMSUBADD_pd(a,b,c) _mm256_fmsubadd_pd(a,b,c)
#endif

extern "C" {

//Forward declaration
#ifdef __AVX__
static __forceinline void hh_trafo_complex_kernel_8_AVX_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s);
static __forceinline void hh_trafo_complex_kernel_6_AVX_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s);
static __forceinline void hh_trafo_complex_kernel_4_AVX_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s);
static __forceinline void hh_trafo_complex_kernel_2_AVX_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s);
#else
static __forceinline void hh_trafo_complex_kernel_4_SSE_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s);
static __forceinline void hh_trafo_complex_kernel_3_SSE_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s);
static __forceinline void hh_trafo_complex_kernel_2_SSE_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s);
static __forceinline void hh_trafo_complex_kernel_1_SSE_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s);
#endif

#if 0
static __forceinline void hh_trafo_complex_kernel_4_C_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s)
{
	std::complex<double> x1;
	std::complex<double> x2;
	std::complex<double> x3;
	std::complex<double> x4;
	std::complex<double> y1;
	std::complex<double> y2;
	std::complex<double> y3;
	std::complex<double> y4;
	std::complex<double> h1;
	std::complex<double> h2;
	std::complex<double> tau1;
	std::complex<double> tau2;
	int i=0;

	x1 = q[ldq+0];
	x2 = q[ldq+1];
	x3 = q[ldq+2];
	x4 = q[ldq+3];

	h2 = conj(hh[ldh+1]);

	y1 = q[0] + (x1*h2);
	y2 = q[1] + (x2*h2);
	y3 = q[2] + (x3*h2);
	y4 = q[3] + (x4*h2);

	for (i = 2; i < nb; i++)
	{
		h1 = conj(hh[i-1]);
		h2 = conj(hh[ldh+i]);

		x1 += (q[(i*ldq)+0] * h1);
		y1 += (q[(i*ldq)+0] * h2);
		x2 += (q[(i*ldq)+1] * h1);
		y2 += (q[(i*ldq)+1] * h2);
		x3 += (q[(i*ldq)+2] * h1);
		y3 += (q[(i*ldq)+2] * h2);
		x4 += (q[(i*ldq)+3] * h1);
		y4 += (q[(i*ldq)+3] * h2);
	}
	h1 = conj(hh[nb-1]);

	x1 += (q[(nb*ldq)+0] * h1);
	x2 += (q[(nb*ldq)+1] * h1);
	x3 += (q[(nb*ldq)+2] * h1);
	x4 += (q[(nb*ldq)+3] * h1);

	tau1 = hh[0];
	tau2 = hh[ldh];

	h1 = (-1.0)*tau1;

	x1 *= h1;
	x2 *= h1;
	x3 *= h1;
	x4 *= h1;

	h1 = (-1.0)*tau2;
	h2 = (-1.0)*tau2;
	h2 *= s;
	y1 = y1*h1 +x1*h2;
	y2 = y2*h1 +x2*h2;
	y3 = y3*h1 +x3*h2;
	y4 = y4*h1 +x4*h2;

	q[0] += y1;
	q[1] += y2;
	q[2] += y3;
	q[3] += y4;

	h2 = hh[ldh+1];
	q[ldq+0] += (x1 + (y1*h2));
	q[ldq+1] += (x2 + (y2*h2));
	q[ldq+2] += (x3 + (y3*h2));
	q[ldq+3] += (x4 + (y4*h2));

	for (i = 2; i < nb; i++)
	{
		h1 = hh[i-1];
		h2 = hh[ldh+i];

		q[(i*ldq)+0] += ((x1*h1) + (y1*h2));
		q[(i*ldq)+1] += ((x2*h1) + (y2*h2));
		q[(i*ldq)+2] += ((x3*h1) + (y3*h2));
		q[(i*ldq)+3] += ((x4*h1) + (y4*h2));
	}

	h1 = hh[nb-1];
	q[(nb*ldq)+0] += (x1*h1);
	q[(nb*ldq)+1] += (x2*h1);
	q[(nb*ldq)+2] += (x3*h1);
	q[(nb*ldq)+3] += (x4*h1);
}
#endif

void double_hh_trafo_complex_sse_avx_2hv_(std::complex<double>* q, std::complex<double>* hh, int* pnb, int* pnq, int* pldq, int* pldh)
{
	int i;
	int nb = *pnb;
	int nq = *pldq;
	int ldq = *pldq;
	int ldh = *pldh;

	std::complex<double> s = conj(hh[(ldh)+1])*1.0;
	for (i = 2; i < nb; i++)
	{
		s += hh[i-1] * conj(hh[(i+ldh)]);
	}

#ifdef __AVX__
#if 1
	for (i = 0; i < nq-4; i+=8)
	{
		hh_trafo_complex_kernel_8_AVX_2hv(&q[i], hh, nb, ldq, ldh, s);
	}
	if (nq-i > 0)
	{
		hh_trafo_complex_kernel_4_AVX_2hv(&q[i], hh, nb, ldq, ldh, s);
	}
#else
	for (i = 0; i < nq-4; i+=6)
	{
		hh_trafo_complex_kernel_6_AVX_2hv(&q[i], hh, nb, ldq, ldh, s);
	}
	if (nq-i > 2)
	{
		hh_trafo_complex_kernel_4_AVX_2hv(&q[i], hh, nb, ldq, ldh, s);
	}
	else if (nq-i > 0)
	{
		hh_trafo_complex_kernel_2_AVX_2hv(&q[i], hh, nb, ldq, ldh, s);
	}
#endif
#else
#if 1
	for (i = 0; i < nq; i+=4)
	{
		hh_trafo_complex_kernel_4_SSE_2hv(&q[i], hh, nb, ldq, ldh, s);
	}
#else
	for (i = 0; i < nq-2; i+=3)
	{
		hh_trafo_complex_kernel_3_SSE_2hv(&q[i], hh, nb, ldq, ldh, s);
	}
	if (nq-i > 1)
	{
		hh_trafo_complex_kernel_2_SSE_2hv(&q[i], hh, nb, ldq, ldh, s);
	}
	else if (nq-i > 0)
	{
		hh_trafo_complex_kernel_1_SSE_2hv(&q[i], hh, nb, ldq, ldh, s);
	}
#endif
#endif
}

#ifdef __AVX__
static __forceinline void hh_trafo_complex_kernel_8_AVX_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s)
{
	double* q_dbl = (double*)q;
	double* hh_dbl = (double*)hh;
	double* s_dbl = (double*)(&s);

	__m256d x1, x2, x3, x4;
	__m256d y1, y2, y3, y4;
	__m256d q1, q2, q3, q4;
	__m256d h1_real, h1_imag, h2_real, h2_imag;
	__m256d tmp1, tmp2, tmp3, tmp4;
	int i=0;

	__m256d sign = (__m256d)_mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000000, 0x8000000000000000);

	x1 = _mm256_load_pd(&q_dbl[(2*ldq)+0]);
	x2 = _mm256_load_pd(&q_dbl[(2*ldq)+4]);
	x3 = _mm256_load_pd(&q_dbl[(2*ldq)+8]);
	x4 = _mm256_load_pd(&q_dbl[(2*ldq)+12]);

	h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h2_imag = _mm256_xor_pd(h2_imag, sign);
#endif

	y1 = _mm256_load_pd(&q_dbl[0]);
	y2 = _mm256_load_pd(&q_dbl[4]);
	y3 = _mm256_load_pd(&q_dbl[8]);
	y4 = _mm256_load_pd(&q_dbl[12]);

	tmp1 = _mm256_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_add_pd(y1, _mm256_FMSUBADD_pd(h2_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm256_add_pd(y2, _mm256_FMSUBADD_pd(h2_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	y2 = _mm256_add_pd(y2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h2_imag, x3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm256_add_pd(y3, _mm256_FMSUBADD_pd(h2_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	y3 = _mm256_add_pd(y3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
	tmp4 = _mm256_mul_pd(h2_imag, x4);
#ifdef __ELPA_USE_FMA__
	y4 = _mm256_add_pd(y4, _mm256_FMSUBADD_pd(h2_real, x4, _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#else
	y4 = _mm256_add_pd(y4, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x4), _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#endif

	for (i = 2; i < nb; i++)
	{
		q1 = _mm256_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm256_load_pd(&q_dbl[(2*i*ldq)+4]);
		q3 = _mm256_load_pd(&q_dbl[(2*i*ldq)+8]);
		q4 = _mm256_load_pd(&q_dbl[(2*i*ldq)+12]);

		h1_real = _mm256_broadcast_sd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm256_broadcast_sd(&hh_dbl[((i-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h1_imag = _mm256_xor_pd(h1_imag, sign);
#endif

		tmp1 = _mm256_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
		x1 = _mm256_add_pd(x1, _mm256_FMSUBADD_pd(h1_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		x1 = _mm256_add_pd(x1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
		x2 = _mm256_add_pd(x2, _mm256_FMSUBADD_pd(h1_real, q2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		x2 = _mm256_add_pd(x2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
		tmp3 = _mm256_mul_pd(h1_imag, q3);
#ifdef __ELPA_USE_FMA__
		x3 = _mm256_add_pd(x3, _mm256_FMSUBADD_pd(h1_real, q3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
		x3 = _mm256_add_pd(x3, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
		tmp4 = _mm256_mul_pd(h1_imag, q4);
#ifdef __ELPA_USE_FMA__
		x4 = _mm256_add_pd(x4, _mm256_FMSUBADD_pd(h1_real, q4, _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#else
		x4 = _mm256_add_pd(x4, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q4), _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#endif

		h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+i)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h2_imag = _mm256_xor_pd(h2_imag, sign);
#endif

		tmp1 = _mm256_mul_pd(h2_imag, q1);
#ifdef __ELPA_USE_FMA__
		y1 = _mm256_add_pd(y1, _mm256_FMSUBADD_pd(h2_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h2_imag, q2);
#ifdef __ELPA_USE_FMA__
		y2 = _mm256_add_pd(y2, _mm256_FMSUBADD_pd(h2_real, q2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		y2 = _mm256_add_pd(y2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
		tmp3 = _mm256_mul_pd(h2_imag, q3);
#ifdef __ELPA_USE_FMA__
		y3 = _mm256_add_pd(y3, _mm256_FMSUBADD_pd(h2_real, q3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
		y3 = _mm256_add_pd(y3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
		tmp4 = _mm256_mul_pd(h2_imag, q4);
#ifdef __ELPA_USE_FMA__
		y4 = _mm256_add_pd(y4, _mm256_FMSUBADD_pd(h2_real, q4, _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#else
		y4 = _mm256_add_pd(y4, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q4), _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#endif
	}

	h1_real = _mm256_broadcast_sd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[((nb-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h1_imag = _mm256_xor_pd(h1_imag, sign);
#endif

	q1 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+4]);
	q3 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+8]);
	q4 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+12]);

	tmp1 = _mm256_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm256_add_pd(x1, _mm256_FMSUBADD_pd(h1_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	x1 = _mm256_add_pd(x1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm256_add_pd(x2, _mm256_FMSUBADD_pd(h1_real, q2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	x2 = _mm256_add_pd(x2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h1_imag, q3);
#ifdef __ELPA_USE_FMA__
	x3 = _mm256_add_pd(x3, _mm256_FMSUBADD_pd(h1_real, q3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	x3 = _mm256_add_pd(x3, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
	tmp4 = _mm256_mul_pd(h1_imag, q4);
#ifdef __ELPA_USE_FMA__
	x4 = _mm256_add_pd(x4, _mm256_FMSUBADD_pd(h1_real, q4, _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#else
	x4 = _mm256_add_pd(x4, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q4), _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#endif

	h1_real = _mm256_broadcast_sd(&hh_dbl[0]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[1]);
	h1_real = _mm256_xor_pd(h1_real, sign);
	h1_imag = _mm256_xor_pd(h1_imag, sign);

	tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	x1 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm256_FMADDSUB_pd(h1_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#else
	x2 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#endif
	tmp3 = _mm256_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
	x3 = _mm256_FMADDSUB_pd(h1_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5));
#else
	x3 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5));
#endif
	tmp4 = _mm256_mul_pd(h1_imag, x4);
#ifdef __ELPA_USE_FMA__
	x4 = _mm256_FMADDSUB_pd(h1_real, x4, _mm256_shuffle_pd(tmp4, tmp4, 0x5));
#else
	x4 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x4), _mm256_shuffle_pd(tmp4, tmp4, 0x5));
#endif

	h1_real = _mm256_broadcast_sd(&hh_dbl[ldh*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[(ldh*2)+1]);
	h2_real = _mm256_broadcast_sd(&hh_dbl[ldh*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[(ldh*2)+1]);

	h1_real = _mm256_xor_pd(h1_real, sign);
	h1_imag = _mm256_xor_pd(h1_imag, sign);
	h2_real = _mm256_xor_pd(h2_real, sign);
	h2_imag = _mm256_xor_pd(h2_imag, sign);

	__m128d tmp_s_128 = _mm_loadu_pd(s_dbl);
	tmp2 = _mm256_broadcast_pd(&tmp_s_128);
	tmp1 = _mm256_mul_pd(h2_imag, tmp2);
#ifdef __ELPA_USE_FMA__
	tmp2 = _mm256_FMADDSUB_pd(h2_real, tmp2, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	tmp2 = _mm256_addsub_pd( _mm256_mul_pd(h2_real, tmp2), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	_mm_storeu_pd(s_dbl, _mm256_castpd256_pd128(tmp2));
	h2_real = _mm256_broadcast_sd(&s_dbl[0]);
	h2_imag = _mm256_broadcast_sd(&s_dbl[1]);

	tmp1 = _mm256_mul_pd(h1_imag, y1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_FMADDSUB_pd(h1_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	y1 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, y2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm256_FMADDSUB_pd(h1_real, y2, _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#else
	y2 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y2), _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#endif
	tmp3 = _mm256_mul_pd(h1_imag, y3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm256_FMADDSUB_pd(h1_real, y3, _mm256_shuffle_pd(tmp3, tmp3, 0x5));
#else
	y3 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y3), _mm256_shuffle_pd(tmp3, tmp3, 0x5));
#endif
	tmp4 = _mm256_mul_pd(h1_imag, y4);
#ifdef __ELPA_USE_FMA__
	y4 = _mm256_FMADDSUB_pd(h1_real, y4, _mm256_shuffle_pd(tmp4, tmp4, 0x5));
#else
	y4 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y4), _mm256_shuffle_pd(tmp4, tmp4, 0x5));
#endif

	tmp1 = _mm256_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_add_pd(y1, _mm256_FMADDSUB_pd(h2_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm256_add_pd(y2, _mm256_FMADDSUB_pd(h2_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	y2 = _mm256_add_pd(y2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h2_imag, x3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm256_add_pd(y3, _mm256_FMADDSUB_pd(h2_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	y3 = _mm256_add_pd(y3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
	tmp4 = _mm256_mul_pd(h2_imag, x4);
#ifdef __ELPA_USE_FMA__
	y4 = _mm256_add_pd(y4, _mm256_FMADDSUB_pd(h2_real, x4, _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#else
	y4 = _mm256_add_pd(y4, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x4), _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#endif

	q1 = _mm256_load_pd(&q_dbl[0]);
	q2 = _mm256_load_pd(&q_dbl[4]);
	q3 = _mm256_load_pd(&q_dbl[8]);
	q4 = _mm256_load_pd(&q_dbl[12]);

	q1 = _mm256_add_pd(q1, y1);
	q2 = _mm256_add_pd(q2, y2);
	q3 = _mm256_add_pd(q3, y3);
	q4 = _mm256_add_pd(q4, y4);

	_mm256_store_pd(&q_dbl[0], q1);
	_mm256_store_pd(&q_dbl[4], q2);
	_mm256_store_pd(&q_dbl[8], q3);
	_mm256_store_pd(&q_dbl[12], q4);

	h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+1)*2)+1]);

	q1 = _mm256_load_pd(&q_dbl[(ldq*2)+0]);
	q2 = _mm256_load_pd(&q_dbl[(ldq*2)+4]);
	q3 = _mm256_load_pd(&q_dbl[(ldq*2)+8]);
	q4 = _mm256_load_pd(&q_dbl[(ldq*2)+12]);

	q1 = _mm256_add_pd(q1, x1);
	q2 = _mm256_add_pd(q2, x2);
	q3 = _mm256_add_pd(q3, x3);
	q4 = _mm256_add_pd(q4, x4);

	tmp1 = _mm256_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h2_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h2_imag, y2);
#ifdef __FMA4_
	q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h2_real, y2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h2_imag, y3);
#ifdef __ELPA_USE_FMA__
	q3 = _mm256_add_pd(q3, _mm256_FMADDSUB_pd(h2_real, y3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	q3 = _mm256_add_pd(q3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
	tmp4 = _mm256_mul_pd(h2_imag, y4);
#ifdef __ELPA_USE_FMA__
	q4 = _mm256_add_pd(q4, _mm256_FMADDSUB_pd(h2_real, y4, _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#else
	q4 = _mm256_add_pd(q4, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y4), _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#endif

	_mm256_store_pd(&q_dbl[(ldq*2)+0], q1);
	_mm256_store_pd(&q_dbl[(ldq*2)+4], q2);
	_mm256_store_pd(&q_dbl[(ldq*2)+8], q3);
	_mm256_store_pd(&q_dbl[(ldq*2)+12], q4);

	for (i = 2; i < nb; i++)
	{
		q1 = _mm256_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm256_load_pd(&q_dbl[(2*i*ldq)+4]);
		q3 = _mm256_load_pd(&q_dbl[(2*i*ldq)+8]);
		q4 = _mm256_load_pd(&q_dbl[(2*i*ldq)+12]);

		h1_real = _mm256_broadcast_sd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm256_broadcast_sd(&hh_dbl[((i-1)*2)+1]);

		tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h1_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
		tmp3 = _mm256_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
		q3 = _mm256_add_pd(q3, _mm256_FMADDSUB_pd(h1_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
		q3 = _mm256_add_pd(q3, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
		tmp4 = _mm256_mul_pd(h1_imag, x4);
#ifdef __ELPA_USE_FMA__
		q4 = _mm256_add_pd(q4, _mm256_FMADDSUB_pd(h1_real, x4, _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#else
		q4 = _mm256_add_pd(q4, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x4), _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#endif

		h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+i)*2)+1]);

		tmp1 = _mm256_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h2_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h2_imag, y2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h2_real, y2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
		tmp3 = _mm256_mul_pd(h2_imag, y3);
#ifdef __ELPA_USE_FMA__
		q3 = _mm256_add_pd(q3, _mm256_FMADDSUB_pd(h2_real, y3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
		q3 = _mm256_add_pd(q3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
		tmp4 = _mm256_mul_pd(h2_imag, y4);
#ifdef __ELPA_USE_FMA__
		q4 = _mm256_add_pd(q4, _mm256_FMADDSUB_pd(h2_real, y4, _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#else
		q4 = _mm256_add_pd(q4, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y4), _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#endif

		_mm256_store_pd(&q_dbl[(2*i*ldq)+0], q1);
		_mm256_store_pd(&q_dbl[(2*i*ldq)+4], q2);
		_mm256_store_pd(&q_dbl[(2*i*ldq)+8], q3);
		_mm256_store_pd(&q_dbl[(2*i*ldq)+12], q4);
	}
	h1_real = _mm256_broadcast_sd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[((nb-1)*2)+1]);

	q1 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+4]);
	q3 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+8]);
	q4 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+12]);

	tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h1_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
	q3 = _mm256_add_pd(q3, _mm256_FMADDSUB_pd(h1_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	q3 = _mm256_add_pd(q3, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
	tmp4 = _mm256_mul_pd(h1_imag, x4);
#ifdef __ELPA_USE_FMA__
	q4 = _mm256_add_pd(q4, _mm256_FMADDSUB_pd(h1_real, x4, _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#else
	q4 = _mm256_add_pd(q4, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x4), _mm256_shuffle_pd(tmp4, tmp4, 0x5)));
#endif

	_mm256_store_pd(&q_dbl[(2*nb*ldq)+0], q1);
	_mm256_store_pd(&q_dbl[(2*nb*ldq)+4], q2);
	_mm256_store_pd(&q_dbl[(2*nb*ldq)+8], q3);
	_mm256_store_pd(&q_dbl[(2*nb*ldq)+12], q4);
}

static __forceinline void hh_trafo_complex_kernel_6_AVX_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s)
{
	double* q_dbl = (double*)q;
	double* hh_dbl = (double*)hh;
	double* s_dbl = (double*)(&s);

	__m256d x1, x2, x3;
	__m256d y1, y2, y3;
	__m256d q1, q2, q3;
	__m256d h1_real, h1_imag, h2_real, h2_imag;
	__m256d tmp1, tmp2, tmp3;
	int i=0;

	__m256d sign = (__m256d)_mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000000, 0x8000000000000000);

	x1 = _mm256_load_pd(&q_dbl[(2*ldq)+0]);
	x2 = _mm256_load_pd(&q_dbl[(2*ldq)+4]);
	x3 = _mm256_load_pd(&q_dbl[(2*ldq)+8]);

	h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h2_imag = _mm256_xor_pd(h2_imag, sign);
#endif

	y1 = _mm256_load_pd(&q_dbl[0]);
	y2 = _mm256_load_pd(&q_dbl[4]);
	y3 = _mm256_load_pd(&q_dbl[8]);

	tmp1 = _mm256_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_add_pd(y1, _mm256_FMSUBADD_pd(h2_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm256_add_pd(y2, _mm256_FMSUBADD_pd(h2_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	y2 = _mm256_add_pd(y2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h2_imag, x3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm256_add_pd(y3, _mm256_FMSUBADD_pd(h2_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	y3 = _mm256_add_pd(y3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif

	for (i = 2; i < nb; i++)
	{
		q1 = _mm256_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm256_load_pd(&q_dbl[(2*i*ldq)+4]);
		q3 = _mm256_load_pd(&q_dbl[(2*i*ldq)+8]);

		h1_real = _mm256_broadcast_sd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm256_broadcast_sd(&hh_dbl[((i-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h1_imag = _mm256_xor_pd(h1_imag, sign);
#endif

		tmp1 = _mm256_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
		x1 = _mm256_add_pd(x1, _mm256_FMSUBADD_pd(h1_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		x1 = _mm256_add_pd(x1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
		x2 = _mm256_add_pd(x2, _mm256_FMSUBADD_pd(h1_real, q2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		x2 = _mm256_add_pd(x2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
		tmp3 = _mm256_mul_pd(h1_imag, q3);
#ifdef __ELPA_USE_FMA__
		x3 = _mm256_add_pd(x3, _mm256_FMSUBADD_pd(h1_real, q3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
		x3 = _mm256_add_pd(x3, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif

		h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+i)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h2_imag = _mm256_xor_pd(h2_imag, sign);
#endif

		tmp1 = _mm256_mul_pd(h2_imag, q1);
#ifdef __ELPA_USE_FMA__
		y1 = _mm256_add_pd(y1, _mm256_FMSUBADD_pd(h2_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h2_imag, q2);
#ifdef __ELPA_USE_FMA__
		y2 = _mm256_add_pd(y2, _mm256_FMSUBADD_pd(h2_real, q2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		y2 = _mm256_add_pd(y2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
		tmp3 = _mm256_mul_pd(h2_imag, q3);
#ifdef __ELPA_USE_FMA__
		y3 = _mm256_add_pd(y3, _mm256_FMSUBADD_pd(h2_real, q3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
		y3 = _mm256_add_pd(y3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif
	}

	h1_real = _mm256_broadcast_sd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[((nb-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h1_imag = _mm256_xor_pd(h1_imag, sign);
#endif

	q1 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+4]);
	q3 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+8]);

	tmp1 = _mm256_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm256_add_pd(x1, _mm256_FMSUBADD_pd(h1_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	x1 = _mm256_add_pd(x1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm256_add_pd(x2, _mm256_FMSUBADD_pd(h1_real, q2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	x2 = _mm256_add_pd(x2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h1_imag, q3);
#ifdef __ELPA_USE_FMA__
	x3 = _mm256_add_pd(x3, _mm256_FMSUBADD_pd(h1_real, q3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	x3 = _mm256_add_pd(x3, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif

	h1_real = _mm256_broadcast_sd(&hh_dbl[0]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[1]);
	h1_real = _mm256_xor_pd(h1_real, sign);
	h1_imag = _mm256_xor_pd(h1_imag, sign);

	tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	x1 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm256_FMADDSUB_pd(h1_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#else
	x2 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#endif
	tmp3 = _mm256_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
	x3 = _mm256_FMADDSUB_pd(h1_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5));
#else
	x3 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5));
#endif

	h1_real = _mm256_broadcast_sd(&hh_dbl[ldh*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[(ldh*2)+1]);
	h2_real = _mm256_broadcast_sd(&hh_dbl[ldh*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[(ldh*2)+1]);

	h1_real = _mm256_xor_pd(h1_real, sign);
	h1_imag = _mm256_xor_pd(h1_imag, sign);
	h2_real = _mm256_xor_pd(h2_real, sign);
	h2_imag = _mm256_xor_pd(h2_imag, sign);

	__m128d tmp_s_128 = _mm_loadu_pd(s_dbl);
	tmp2 = _mm256_broadcast_pd(&tmp_s_128);
	tmp1 = _mm256_mul_pd(h2_imag, tmp2);
#ifdef __ELPA_USE_FMA__
	tmp2 = _mm256_FMADDSUB_pd(h2_real, tmp2, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	tmp2 = _mm256_addsub_pd( _mm256_mul_pd(h2_real, tmp2), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	_mm_storeu_pd(s_dbl, _mm256_castpd256_pd128(tmp2));
	h2_real = _mm256_broadcast_sd(&s_dbl[0]);
	h2_imag = _mm256_broadcast_sd(&s_dbl[1]);

	tmp1 = _mm256_mul_pd(h1_imag, y1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_FMADDSUB_pd(h1_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	y1 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, y2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm256_FMADDSUB_pd(h1_real, y2, _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#else
	y2 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y2), _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#endif
	tmp3 = _mm256_mul_pd(h1_imag, y3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm256_FMADDSUB_pd(h1_real, y3, _mm256_shuffle_pd(tmp3, tmp3, 0x5));
#else
	y3 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y3), _mm256_shuffle_pd(tmp3, tmp3, 0x5));
#endif

	tmp1 = _mm256_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_add_pd(y1, _mm256_FMADDSUB_pd(h2_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm256_add_pd(y2, _mm256_FMADDSUB_pd(h2_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	y2 = _mm256_add_pd(y2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h2_imag, x3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm256_add_pd(y3, _mm256_FMADDSUB_pd(h2_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	y3 = _mm256_add_pd(y3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif

	q1 = _mm256_load_pd(&q_dbl[0]);
	q2 = _mm256_load_pd(&q_dbl[4]);
	q3 = _mm256_load_pd(&q_dbl[8]);

	q1 = _mm256_add_pd(q1, y1);
	q2 = _mm256_add_pd(q2, y2);
	q3 = _mm256_add_pd(q3, y3);

	_mm256_store_pd(&q_dbl[0], q1);
	_mm256_store_pd(&q_dbl[4], q2);
	_mm256_store_pd(&q_dbl[8], q3);

	h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+1)*2)+1]);

	q1 = _mm256_load_pd(&q_dbl[(ldq*2)+0]);
	q2 = _mm256_load_pd(&q_dbl[(ldq*2)+4]);
	q3 = _mm256_load_pd(&q_dbl[(ldq*2)+8]);

	q1 = _mm256_add_pd(q1, x1);
	q2 = _mm256_add_pd(q2, x2);
	q3 = _mm256_add_pd(q3, x3);

	tmp1 = _mm256_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h2_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h2_imag, y2);
#ifdef __FMA4_
	q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h2_real, y2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h2_imag, y3);
#ifdef __ELPA_USE_FMA__
	q3 = _mm256_add_pd(q3, _mm256_FMADDSUB_pd(h2_real, y3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	q3 = _mm256_add_pd(q3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif

	_mm256_store_pd(&q_dbl[(ldq*2)+0], q1);
	_mm256_store_pd(&q_dbl[(ldq*2)+4], q2);
	_mm256_store_pd(&q_dbl[(ldq*2)+8], q3);

	for (i = 2; i < nb; i++)
	{
		q1 = _mm256_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm256_load_pd(&q_dbl[(2*i*ldq)+4]);
		q3 = _mm256_load_pd(&q_dbl[(2*i*ldq)+8]);

		h1_real = _mm256_broadcast_sd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm256_broadcast_sd(&hh_dbl[((i-1)*2)+1]);

		tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h1_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
		tmp3 = _mm256_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
		q3 = _mm256_add_pd(q3, _mm256_FMADDSUB_pd(h1_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
		q3 = _mm256_add_pd(q3, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif

		h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+i)*2)+1]);

		tmp1 = _mm256_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h2_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h2_imag, y2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h2_real, y2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
		tmp3 = _mm256_mul_pd(h2_imag, y3);
#ifdef __ELPA_USE_FMA__
		q3 = _mm256_add_pd(q3, _mm256_FMADDSUB_pd(h2_real, y3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
		q3 = _mm256_add_pd(q3, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif

		_mm256_store_pd(&q_dbl[(2*i*ldq)+0], q1);
		_mm256_store_pd(&q_dbl[(2*i*ldq)+4], q2);
		_mm256_store_pd(&q_dbl[(2*i*ldq)+8], q3);
	}
	h1_real = _mm256_broadcast_sd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[((nb-1)*2)+1]);

	q1 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+4]);
	q3 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+8]);

	tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h1_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	tmp3 = _mm256_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
	q3 = _mm256_add_pd(q3, _mm256_FMADDSUB_pd(h1_real, x3, _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#else
	q3 = _mm256_add_pd(q3, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x3), _mm256_shuffle_pd(tmp3, tmp3, 0x5)));
#endif

	_mm256_store_pd(&q_dbl[(2*nb*ldq)+0], q1);
	_mm256_store_pd(&q_dbl[(2*nb*ldq)+4], q2);
	_mm256_store_pd(&q_dbl[(2*nb*ldq)+8], q3);
}

static __forceinline void hh_trafo_complex_kernel_4_AVX_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s)
{
	double* q_dbl = (double*)q;
	double* hh_dbl = (double*)hh;
	double* s_dbl = (double*)(&s);

	__m256d x1, x2;
	__m256d y1, y2;
	__m256d q1, q2;
	__m256d h1_real, h1_imag, h2_real, h2_imag;
	__m256d tmp1, tmp2;
	int i=0;

	__m256d sign = (__m256d)_mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000000, 0x8000000000000000);

	x1 = _mm256_load_pd(&q_dbl[(2*ldq)+0]);
	x2 = _mm256_load_pd(&q_dbl[(2*ldq)+4]);

	h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h2_imag = _mm256_xor_pd(h2_imag, sign);
#endif

	y1 = _mm256_load_pd(&q_dbl[0]);
	y2 = _mm256_load_pd(&q_dbl[4]);

	tmp1 = _mm256_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_add_pd(y1, _mm256_FMSUBADD_pd(h2_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm256_add_pd(y2, _mm256_FMSUBADD_pd(h2_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	y2 = _mm256_add_pd(y2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif

	for (i = 2; i < nb; i++)
	{
		q1 = _mm256_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm256_load_pd(&q_dbl[(2*i*ldq)+4]);

		h1_real = _mm256_broadcast_sd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm256_broadcast_sd(&hh_dbl[((i-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h1_imag = _mm256_xor_pd(h1_imag, sign);
#endif

		tmp1 = _mm256_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
		x1 = _mm256_add_pd(x1, _mm256_FMSUBADD_pd(h1_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		x1 = _mm256_add_pd(x1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
		x2 = _mm256_add_pd(x2, _mm256_FMSUBADD_pd(h1_real, q2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		x2 = _mm256_add_pd(x2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif

		h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+i)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h2_imag = _mm256_xor_pd(h2_imag, sign);
#endif

		tmp1 = _mm256_mul_pd(h2_imag, q1);
#ifdef __ELPA_USE_FMA__
		y1 = _mm256_add_pd(y1, _mm256_FMSUBADD_pd(h2_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h2_imag, q2);
#ifdef __ELPA_USE_FMA__
		y2 = _mm256_add_pd(y2, _mm256_FMSUBADD_pd(h2_real, q2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		y2 = _mm256_add_pd(y2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif
	}

	h1_real = _mm256_broadcast_sd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[((nb-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h1_imag = _mm256_xor_pd(h1_imag, sign);
#endif

	q1 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+4]);

	tmp1 = _mm256_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm256_add_pd(x1, _mm256_FMSUBADD_pd(h1_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	x1 = _mm256_add_pd(x1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm256_add_pd(x2, _mm256_FMSUBADD_pd(h1_real, q2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	x2 = _mm256_add_pd(x2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif

	h1_real = _mm256_broadcast_sd(&hh_dbl[0]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[1]);
	h1_real = _mm256_xor_pd(h1_real, sign);
	h1_imag = _mm256_xor_pd(h1_imag, sign);

	tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	x1 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm256_FMADDSUB_pd(h1_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#else
	x2 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#endif

	h1_real = _mm256_broadcast_sd(&hh_dbl[ldh*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[(ldh*2)+1]);
	h2_real = _mm256_broadcast_sd(&hh_dbl[ldh*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[(ldh*2)+1]);

	h1_real = _mm256_xor_pd(h1_real, sign);
	h1_imag = _mm256_xor_pd(h1_imag, sign);
	h2_real = _mm256_xor_pd(h2_real, sign);
	h2_imag = _mm256_xor_pd(h2_imag, sign);

	__m128d tmp_s_128 = _mm_loadu_pd(s_dbl);
	tmp2 = _mm256_broadcast_pd(&tmp_s_128);
	tmp1 = _mm256_mul_pd(h2_imag, tmp2);
#ifdef __ELPA_USE_FMA__
	tmp2 = _mm256_FMADDSUB_pd(h2_real, tmp2, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	tmp2 = _mm256_addsub_pd( _mm256_mul_pd(h2_real, tmp2), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	_mm_storeu_pd(s_dbl, _mm256_castpd256_pd128(tmp2));
	h2_real = _mm256_broadcast_sd(&s_dbl[0]);
	h2_imag = _mm256_broadcast_sd(&s_dbl[1]);

	tmp1 = _mm256_mul_pd(h1_imag, y1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_FMADDSUB_pd(h1_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	y1 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, y2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm256_FMADDSUB_pd(h1_real, y2, _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#else
	y2 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y2), _mm256_shuffle_pd(tmp2, tmp2, 0x5));
#endif

	tmp1 = _mm256_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_add_pd(y1, _mm256_FMADDSUB_pd(h2_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm256_add_pd(y2, _mm256_FMADDSUB_pd(h2_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	y2 = _mm256_add_pd(y2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif

	q1 = _mm256_load_pd(&q_dbl[0]);
	q2 = _mm256_load_pd(&q_dbl[4]);

	q1 = _mm256_add_pd(q1, y1);
	q2 = _mm256_add_pd(q2, y2);

	_mm256_store_pd(&q_dbl[0], q1);
	_mm256_store_pd(&q_dbl[4], q2);

	h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+1)*2)+1]);

	q1 = _mm256_load_pd(&q_dbl[(ldq*2)+0]);
	q2 = _mm256_load_pd(&q_dbl[(ldq*2)+4]);

	q1 = _mm256_add_pd(q1, x1);
	q2 = _mm256_add_pd(q2, x2);

	tmp1 = _mm256_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h2_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h2_imag, y2);
#ifdef __FMA4_
	q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h2_real, y2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif

	_mm256_store_pd(&q_dbl[(ldq*2)+0], q1);
	_mm256_store_pd(&q_dbl[(ldq*2)+4], q2);

	for (i = 2; i < nb; i++)
	{
		q1 = _mm256_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm256_load_pd(&q_dbl[(2*i*ldq)+4]);

		h1_real = _mm256_broadcast_sd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm256_broadcast_sd(&hh_dbl[((i-1)*2)+1]);

		tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h1_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif

		h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+i)*2)+1]);

		tmp1 = _mm256_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h2_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
		tmp2 = _mm256_mul_pd(h2_imag, y2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h2_real, y2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
		q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif

		_mm256_store_pd(&q_dbl[(2*i*ldq)+0], q1);
		_mm256_store_pd(&q_dbl[(2*i*ldq)+4], q2);
	}
	h1_real = _mm256_broadcast_sd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[((nb-1)*2)+1]);

	q1 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+4]);

	tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	tmp2 = _mm256_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	q2 = _mm256_add_pd(q2, _mm256_FMADDSUB_pd(h1_real, x2, _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#else
	q2 = _mm256_add_pd(q2, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x2), _mm256_shuffle_pd(tmp2, tmp2, 0x5)));
#endif

	_mm256_store_pd(&q_dbl[(2*nb*ldq)+0], q1);
	_mm256_store_pd(&q_dbl[(2*nb*ldq)+4], q2);
}

static __forceinline void hh_trafo_complex_kernel_2_AVX_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s)
{
	double* q_dbl = (double*)q;
	double* hh_dbl = (double*)hh;
	double* s_dbl = (double*)(&s);

	__m256d x1;
	__m256d y1;
	__m256d q1;
	__m256d h1_real, h1_imag, h2_real, h2_imag;
	__m256d tmp1;
	int i=0;

	__m256d sign = (__m256d)_mm256_set_epi64x(0x8000000000000000, 0x8000000000000000, 0x8000000000000000, 0x8000000000000000);

	x1 = _mm256_load_pd(&q_dbl[(2*ldq)+0]);

	h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h2_imag = _mm256_xor_pd(h2_imag, sign);
#endif

	y1 = _mm256_load_pd(&q_dbl[0]);

	tmp1 = _mm256_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_add_pd(y1, _mm256_FMSUBADD_pd(h2_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif

	for (i = 2; i < nb; i++)
	{
		q1 = _mm256_load_pd(&q_dbl[(2*i*ldq)+0]);

		h1_real = _mm256_broadcast_sd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm256_broadcast_sd(&hh_dbl[((i-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h1_imag = _mm256_xor_pd(h1_imag, sign);
#endif

		tmp1 = _mm256_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
		x1 = _mm256_add_pd(x1, _mm256_FMSUBADD_pd(h1_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		x1 = _mm256_add_pd(x1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif

		h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+i)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h2_imag = _mm256_xor_pd(h2_imag, sign);
#endif

		tmp1 = _mm256_mul_pd(h2_imag, q1);
#ifdef __ELPA_USE_FMA__
		y1 = _mm256_add_pd(y1, _mm256_FMSUBADD_pd(h2_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif
	}

	h1_real = _mm256_broadcast_sd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[((nb-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h1_imag = _mm256_xor_pd(h1_imag, sign);
#endif

	q1 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+0]);

	tmp1 = _mm256_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm256_add_pd(x1, _mm256_FMSUBADD_pd(h1_real, q1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	x1 = _mm256_add_pd(x1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, q1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif

	h1_real = _mm256_broadcast_sd(&hh_dbl[0]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[1]);
	h1_real = _mm256_xor_pd(h1_real, sign);
	h1_imag = _mm256_xor_pd(h1_imag, sign);

	tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	x1 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif

	h1_real = _mm256_broadcast_sd(&hh_dbl[ldh*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[(ldh*2)+1]);
	h2_real = _mm256_broadcast_sd(&hh_dbl[ldh*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[(ldh*2)+1]);

	h1_real = _mm256_xor_pd(h1_real, sign);
	h1_imag = _mm256_xor_pd(h1_imag, sign);
	h2_real = _mm256_xor_pd(h2_real, sign);
	h2_imag = _mm256_xor_pd(h2_imag, sign);

	__m128d tmp_s_128 = _mm_loadu_pd(s_dbl);
	__m256d tmp2 = _mm256_broadcast_pd(&tmp_s_128);
	tmp1 = _mm256_mul_pd(h2_imag, tmp2);
#ifdef __ELPA_USE_FMA__
	tmp2 = _mm256_FMADDSUB_pd(h2_real, tmp2, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	tmp2 = _mm256_addsub_pd( _mm256_mul_pd(h2_real, tmp2), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif
	_mm_storeu_pd(s_dbl, _mm256_castpd256_pd128(tmp2));
	h2_real = _mm256_broadcast_sd(&s_dbl[0]);
	h2_imag = _mm256_broadcast_sd(&s_dbl[1]);

	tmp1 = _mm256_mul_pd(h1_imag, y1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_FMADDSUB_pd(h1_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#else
	y1 = _mm256_addsub_pd( _mm256_mul_pd(h1_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5));
#endif

	tmp1 = _mm256_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm256_add_pd(y1, _mm256_FMADDSUB_pd(h2_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	y1 = _mm256_add_pd(y1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif

	q1 = _mm256_load_pd(&q_dbl[0]);

	q1 = _mm256_add_pd(q1, y1);

	_mm256_store_pd(&q_dbl[0], q1);

	h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+1)*2)+1]);

	q1 = _mm256_load_pd(&q_dbl[(ldq*2)+0]);

	q1 = _mm256_add_pd(q1, x1);

	tmp1 = _mm256_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h2_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif

	_mm256_store_pd(&q_dbl[(ldq*2)+0], q1);

	for (i = 2; i < nb; i++)
	{
		q1 = _mm256_load_pd(&q_dbl[(2*i*ldq)+0]);

		h1_real = _mm256_broadcast_sd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm256_broadcast_sd(&hh_dbl[((i-1)*2)+1]);

		tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif

		h2_real = _mm256_broadcast_sd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm256_broadcast_sd(&hh_dbl[((ldh+i)*2)+1]);

		tmp1 = _mm256_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h2_real, y1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
		q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h2_real, y1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif

		_mm256_store_pd(&q_dbl[(2*i*ldq)+0], q1);
	}
	h1_real = _mm256_broadcast_sd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm256_broadcast_sd(&hh_dbl[((nb-1)*2)+1]);

	q1 = _mm256_load_pd(&q_dbl[(2*nb*ldq)+0]);

	tmp1 = _mm256_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm256_add_pd(q1, _mm256_FMADDSUB_pd(h1_real, x1, _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#else
	q1 = _mm256_add_pd(q1, _mm256_addsub_pd( _mm256_mul_pd(h1_real, x1), _mm256_shuffle_pd(tmp1, tmp1, 0x5)));
#endif

	_mm256_store_pd(&q_dbl[(2*nb*ldq)+0], q1);
}
#else
static __forceinline void hh_trafo_complex_kernel_4_SSE_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s)
{
	double* q_dbl = (double*)q;
	double* hh_dbl = (double*)hh;
	double* s_dbl = (double*)(&s);

	__m128d x1, x2, x3, x4;
	__m128d y1, y2, y3, y4;
	__m128d q1, q2, q3, q4;
	__m128d h1_real, h1_imag, h2_real, h2_imag;
	__m128d tmp1, tmp2, tmp3, tmp4;
	int i=0;

	__m128d sign = (__m128d)_mm_set_epi64x(0x8000000000000000, 0x8000000000000000);

	x1 = _mm_load_pd(&q_dbl[(2*ldq)+0]);
	x2 = _mm_load_pd(&q_dbl[(2*ldq)+2]);
	x3 = _mm_load_pd(&q_dbl[(2*ldq)+4]);
	x4 = _mm_load_pd(&q_dbl[(2*ldq)+6]);

	h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h2_imag = _mm_xor_pd(h2_imag, sign);
#endif

	y1 = _mm_load_pd(&q_dbl[0]);
	y2 = _mm_load_pd(&q_dbl[2]);
	y3 = _mm_load_pd(&q_dbl[4]);
	y4 = _mm_load_pd(&q_dbl[6]);

	tmp1 = _mm_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_add_pd(y1, _mm_msubadd_pd(h2_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm_add_pd(y2, _mm_msubadd_pd(h2_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	y2 = _mm_add_pd(y2, _mm_addsub_pd( _mm_mul_pd(h2_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h2_imag, x3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm_add_pd(y3, _mm_msubadd_pd(h2_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	y3 = _mm_add_pd(y3, _mm_addsub_pd( _mm_mul_pd(h2_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
	tmp4 = _mm_mul_pd(h2_imag, x4);
#ifdef __ELPA_USE_FMA__
	y4 = _mm_add_pd(y4, _mm_msubadd_pd(h2_real, x4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#else
	y4 = _mm_add_pd(y4, _mm_addsub_pd( _mm_mul_pd(h2_real, x4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#endif

	for (i = 2; i < nb; i++)
	{
		q1 = _mm_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm_load_pd(&q_dbl[(2*i*ldq)+2]);
		q3 = _mm_load_pd(&q_dbl[(2*i*ldq)+4]);
		q4 = _mm_load_pd(&q_dbl[(2*i*ldq)+6]);

		h1_real = _mm_loaddup_pd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm_loaddup_pd(&hh_dbl[((i-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h1_imag = _mm_xor_pd(h1_imag, sign);
#endif

		tmp1 = _mm_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
		x1 = _mm_add_pd(x1, _mm_msubadd_pd(h1_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		x1 = _mm_add_pd(x1, _mm_addsub_pd( _mm_mul_pd(h1_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
		x2 = _mm_add_pd(x2, _mm_msubadd_pd(h1_real, q2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		x2 = _mm_add_pd(x2, _mm_addsub_pd( _mm_mul_pd(h1_real, q2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
		tmp3 = _mm_mul_pd(h1_imag, q3);
#ifdef __ELPA_USE_FMA__
		x3 = _mm_add_pd(x3, _mm_msubadd_pd(h1_real, q3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
		x3 = _mm_add_pd(x3, _mm_addsub_pd( _mm_mul_pd(h1_real, q3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
		tmp4 = _mm_mul_pd(h1_imag, q4);
#ifdef __ELPA_USE_FMA__
		x4 = _mm_add_pd(x4, _mm_msubadd_pd(h1_real, q4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#else
		x4 = _mm_add_pd(x4, _mm_addsub_pd( _mm_mul_pd(h1_real, q4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#endif

		h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+i)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h2_imag = _mm_xor_pd(h2_imag, sign);
#endif

		tmp1 = _mm_mul_pd(h2_imag, q1);
#ifdef __ELPA_USE_FMA__
		y1 = _mm_add_pd(y1, _mm_msubadd_pd(h2_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h2_imag, q2);
#ifdef __ELPA_USE_FMA__
		y2 = _mm_add_pd(y2, _mm_msubadd_pd(h2_real, q2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		y2 = _mm_add_pd(y2, _mm_addsub_pd( _mm_mul_pd(h2_real, q2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
		tmp3 = _mm_mul_pd(h2_imag, q3);
#ifdef __ELPA_USE_FMA__
		y3 = _mm_add_pd(y3, _mm_msubadd_pd(h2_real, q3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
		y3 = _mm_add_pd(y3, _mm_addsub_pd( _mm_mul_pd(h2_real, q3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
		tmp4 = _mm_mul_pd(h2_imag, q4);
#ifdef __ELPA_USE_FMA__
		y4 = _mm_add_pd(y4, _mm_msubadd_pd(h2_real, q4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#else
		y4 = _mm_add_pd(y4, _mm_addsub_pd( _mm_mul_pd(h2_real, q4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#endif
	}

	h1_real = _mm_loaddup_pd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[((nb-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h1_imag = _mm_xor_pd(h1_imag, sign);
#endif

	q1 = _mm_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm_load_pd(&q_dbl[(2*nb*ldq)+2]);
	q3 = _mm_load_pd(&q_dbl[(2*nb*ldq)+4]);
	q4 = _mm_load_pd(&q_dbl[(2*nb*ldq)+6]);

	tmp1 = _mm_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm_add_pd(x1, _mm_msubadd_pd(h1_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	x1 = _mm_add_pd(x1, _mm_addsub_pd( _mm_mul_pd(h1_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm_add_pd(x2, _mm_msubadd_pd(h1_real, q2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	x2 = _mm_add_pd(x2, _mm_addsub_pd( _mm_mul_pd(h1_real, q2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h1_imag, q3);
#ifdef __ELPA_USE_FMA__
	x3 = _mm_add_pd(x3, _mm_msubadd_pd(h1_real, q3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	x3 = _mm_add_pd(x3, _mm_addsub_pd( _mm_mul_pd(h1_real, q3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
	tmp4 = _mm_mul_pd(h1_imag, q4);
#ifdef __ELPA_USE_FMA__
	x4 = _mm_add_pd(x4, _mm_msubadd_pd(h1_real, q4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#else
	x4 = _mm_add_pd(x4, _mm_addsub_pd( _mm_mul_pd(h1_real, q4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#endif

	h1_real = _mm_loaddup_pd(&hh_dbl[0]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[1]);
	h1_real = _mm_xor_pd(h1_real, sign);
	h1_imag = _mm_xor_pd(h1_imag, sign);

	tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	x1 = _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	tmp2 = _mm_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm_maddsub_pd(h1_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#else
	x2 = _mm_addsub_pd( _mm_mul_pd(h1_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#endif
	tmp3 = _mm_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
	x3 = _mm_maddsub_pd(h1_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1)));
#else
	x3 = _mm_addsub_pd( _mm_mul_pd(h1_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1)));
#endif
	tmp4 = _mm_mul_pd(h1_imag, x4);
#ifdef __ELPA_USE_FMA__
	x4 = _mm_maddsub_pd(h1_real, x4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1)));
#else
	x4 = _mm_addsub_pd( _mm_mul_pd(h1_real, x4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1)));
#endif

	h1_real = _mm_loaddup_pd(&hh_dbl[ldh*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[(ldh*2)+1]);
	h2_real = _mm_loaddup_pd(&hh_dbl[ldh*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[(ldh*2)+1]);

	h1_real = _mm_xor_pd(h1_real, sign);
	h1_imag = _mm_xor_pd(h1_imag, sign);
	h2_real = _mm_xor_pd(h2_real, sign);
	h2_imag = _mm_xor_pd(h2_imag, sign);

	tmp2 = _mm_loadu_pd(s_dbl);
	tmp1 = _mm_mul_pd(h2_imag, tmp2);
#ifdef __ELPA_USE_FMA__
	tmp2 = _mm_maddsub_pd(h2_real, tmp2, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	tmp2 = _mm_addsub_pd( _mm_mul_pd(h2_real, tmp2), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	_mm_storeu_pd(s_dbl, tmp2);
	h2_real = _mm_loaddup_pd(&s_dbl[0]);
	h2_imag = _mm_loaddup_pd(&s_dbl[1]);

	tmp1 = _mm_mul_pd(h1_imag, y1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_maddsub_pd(h1_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	y1 = _mm_addsub_pd( _mm_mul_pd(h1_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	tmp2 = _mm_mul_pd(h1_imag, y2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm_maddsub_pd(h1_real, y2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#else
	y2 = _mm_addsub_pd( _mm_mul_pd(h1_real, y2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#endif
	tmp3 = _mm_mul_pd(h1_imag, y3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm_maddsub_pd(h1_real, y3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1)));
#else
	y3 = _mm_addsub_pd( _mm_mul_pd(h1_real, y3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1)));
#endif
	tmp4 = _mm_mul_pd(h1_imag, y4);
#ifdef __ELPA_USE_FMA__
	y4 = _mm_maddsub_pd(h1_real, y4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1)));
#else
	y4 = _mm_addsub_pd( _mm_mul_pd(h1_real, y4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1)));
#endif

	tmp1 = _mm_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_add_pd(y1, _mm_maddsub_pd(h2_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm_add_pd(y2, _mm_maddsub_pd(h2_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	y2 = _mm_add_pd(y2, _mm_addsub_pd( _mm_mul_pd(h2_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h2_imag, x3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm_add_pd(y3, _mm_maddsub_pd(h2_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	y3 = _mm_add_pd(y3, _mm_addsub_pd( _mm_mul_pd(h2_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
	tmp4 = _mm_mul_pd(h2_imag, x4);
#ifdef __ELPA_USE_FMA__
	y4 = _mm_add_pd(y4, _mm_maddsub_pd(h2_real, x4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#else
	y4 = _mm_add_pd(y4, _mm_addsub_pd( _mm_mul_pd(h2_real, x4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#endif

	q1 = _mm_load_pd(&q_dbl[0]);
	q2 = _mm_load_pd(&q_dbl[2]);
	q3 = _mm_load_pd(&q_dbl[4]);
	q4 = _mm_load_pd(&q_dbl[6]);

	q1 = _mm_add_pd(q1, y1);
	q2 = _mm_add_pd(q2, y2);
	q3 = _mm_add_pd(q3, y3);
	q4 = _mm_add_pd(q4, y4);

	_mm_store_pd(&q_dbl[0], q1);
	_mm_store_pd(&q_dbl[2], q2);
	_mm_store_pd(&q_dbl[4], q3);
	_mm_store_pd(&q_dbl[6], q4);

	h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+1)*2)+1]);

	q1 = _mm_load_pd(&q_dbl[(ldq*2)+0]);
	q2 = _mm_load_pd(&q_dbl[(ldq*2)+2]);
	q3 = _mm_load_pd(&q_dbl[(ldq*2)+4]);
	q4 = _mm_load_pd(&q_dbl[(ldq*2)+6]);

	q1 = _mm_add_pd(q1, x1);
	q2 = _mm_add_pd(q2, x2);
	q3 = _mm_add_pd(q3, x3);
	q4 = _mm_add_pd(q4, x4);

	tmp1 = _mm_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm_add_pd(q1, _mm_maddsub_pd(h2_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h2_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h2_imag, y2);
#ifdef __ELPA_USE_FMA__
	q2 = _mm_add_pd(q2, _mm_maddsub_pd(h2_real, y2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h2_real, y2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h2_imag, y3);
#ifdef __ELPA_USE_FMA__
	q3 = _mm_add_pd(q3, _mm_maddsub_pd(h2_real, y3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	q3 = _mm_add_pd(q3, _mm_addsub_pd( _mm_mul_pd(h2_real, y3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
	tmp4 = _mm_mul_pd(h2_imag, y4);
#ifdef __ELPA_USE_FMA__
	q4 = _mm_add_pd(q4, _mm_maddsub_pd(h2_real, y4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#else
	q4 = _mm_add_pd(q4, _mm_addsub_pd( _mm_mul_pd(h2_real, y4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#endif

	_mm_store_pd(&q_dbl[(ldq*2)+0], q1);
	_mm_store_pd(&q_dbl[(ldq*2)+2], q2);
	_mm_store_pd(&q_dbl[(ldq*2)+4], q3);
	_mm_store_pd(&q_dbl[(ldq*2)+6], q4);

	for (i = 2; i < nb; i++)
	{
		q1 = _mm_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm_load_pd(&q_dbl[(2*i*ldq)+2]);
		q3 = _mm_load_pd(&q_dbl[(2*i*ldq)+4]);
		q4 = _mm_load_pd(&q_dbl[(2*i*ldq)+6]);

		h1_real = _mm_loaddup_pd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm_loaddup_pd(&hh_dbl[((i-1)*2)+1]);

		tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm_add_pd(q1, _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm_add_pd(q2, _mm_maddsub_pd(h1_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h1_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
		tmp3 = _mm_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
		q3 = _mm_add_pd(q3, _mm_maddsub_pd(h1_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
		q3 = _mm_add_pd(q3, _mm_addsub_pd( _mm_mul_pd(h1_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
		tmp4 = _mm_mul_pd(h1_imag, x4);
#ifdef __ELPA_USE_FMA__
		q4 = _mm_add_pd(q4, _mm_maddsub_pd(h1_real, x4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#else
		q4 = _mm_add_pd(q4, _mm_addsub_pd( _mm_mul_pd(h1_real, x4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#endif

		h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+i)*2)+1]);

		tmp1 = _mm_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm_add_pd(q1, _mm_maddsub_pd(h2_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h2_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h2_imag, y2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm_add_pd(q2, _mm_maddsub_pd(h2_real, y2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h2_real, y2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
		tmp3 = _mm_mul_pd(h2_imag, y3);
#ifdef __ELPA_USE_FMA__
		q3 = _mm_add_pd(q3, _mm_maddsub_pd(h2_real, y3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
		q3 = _mm_add_pd(q3, _mm_addsub_pd( _mm_mul_pd(h2_real, y3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
		tmp4 = _mm_mul_pd(h2_imag, y4);
#ifdef __ELPA_USE_FMA__
		q4 = _mm_add_pd(q4, _mm_maddsub_pd(h2_real, y4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#else
		q4 = _mm_add_pd(q4, _mm_addsub_pd( _mm_mul_pd(h2_real, y4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#endif

		_mm_store_pd(&q_dbl[(2*i*ldq)+0], q1);
		_mm_store_pd(&q_dbl[(2*i*ldq)+2], q2);
		_mm_store_pd(&q_dbl[(2*i*ldq)+4], q3);
		_mm_store_pd(&q_dbl[(2*i*ldq)+6], q4);
	}

	h1_real = _mm_loaddup_pd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[((nb-1)*2)+1]);

	q1 = _mm_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm_load_pd(&q_dbl[(2*nb*ldq)+2]);
	q3 = _mm_load_pd(&q_dbl[(2*nb*ldq)+4]);
	q4 = _mm_load_pd(&q_dbl[(2*nb*ldq)+6]);

	tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm_add_pd(q1, _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	q2 = _mm_add_pd(q2, _mm_maddsub_pd(h1_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h1_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
	q3 = _mm_add_pd(q3, _mm_maddsub_pd(h1_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	q3 = _mm_add_pd(q3, _mm_addsub_pd( _mm_mul_pd(h1_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
	tmp4 = _mm_mul_pd(h1_imag, x4);
#ifdef __ELPA_USE_FMA__
	q4 = _mm_add_pd(q4, _mm_maddsub_pd(h1_real, x4, _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#else
	q4 = _mm_add_pd(q4, _mm_addsub_pd( _mm_mul_pd(h1_real, x4), _mm_shuffle_pd(tmp4, tmp4, _MM_SHUFFLE2(0,1))));
#endif

	_mm_store_pd(&q_dbl[(2*nb*ldq)+0], q1);
	_mm_store_pd(&q_dbl[(2*nb*ldq)+2], q2);
	_mm_store_pd(&q_dbl[(2*nb*ldq)+4], q3);
	_mm_store_pd(&q_dbl[(2*nb*ldq)+6], q4);
}

static __forceinline void hh_trafo_complex_kernel_3_SSE_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s)
{
	double* q_dbl = (double*)q;
	double* hh_dbl = (double*)hh;
	double* s_dbl = (double*)(&s);

	__m128d x1, x2, x3;
	__m128d y1, y2, y3;
	__m128d q1, q2, q3;
	__m128d h1_real, h1_imag, h2_real, h2_imag;
	__m128d tmp1, tmp2, tmp3;
	int i=0;

	__m128d sign = (__m128d)_mm_set_epi64x(0x8000000000000000, 0x8000000000000000);

	x1 = _mm_load_pd(&q_dbl[(2*ldq)+0]);
	x2 = _mm_load_pd(&q_dbl[(2*ldq)+2]);
	x3 = _mm_load_pd(&q_dbl[(2*ldq)+4]);

	h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h2_imag = _mm_xor_pd(h2_imag, sign);
#endif

	y1 = _mm_load_pd(&q_dbl[0]);
	y2 = _mm_load_pd(&q_dbl[2]);
	y3 = _mm_load_pd(&q_dbl[4]);

	tmp1 = _mm_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_add_pd(y1, _mm_msubadd_pd(h2_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm_add_pd(y2, _mm_msubadd_pd(h2_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	y2 = _mm_add_pd(y2, _mm_addsub_pd( _mm_mul_pd(h2_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h2_imag, x3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm_add_pd(y3, _mm_msubadd_pd(h2_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	y3 = _mm_add_pd(y3, _mm_addsub_pd( _mm_mul_pd(h2_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif

	for (i = 2; i < nb; i++)
	{
		q1 = _mm_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm_load_pd(&q_dbl[(2*i*ldq)+2]);
		q3 = _mm_load_pd(&q_dbl[(2*i*ldq)+4]);

		h1_real = _mm_loaddup_pd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm_loaddup_pd(&hh_dbl[((i-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h1_imag = _mm_xor_pd(h1_imag, sign);
#endif

		tmp1 = _mm_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
		x1 = _mm_add_pd(x1, _mm_msubadd_pd(h1_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		x1 = _mm_add_pd(x1, _mm_addsub_pd( _mm_mul_pd(h1_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
		x2 = _mm_add_pd(x2, _mm_msubadd_pd(h1_real, q2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		x2 = _mm_add_pd(x2, _mm_addsub_pd( _mm_mul_pd(h1_real, q2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
		tmp3 = _mm_mul_pd(h1_imag, q3);
#ifdef __ELPA_USE_FMA__
		x3 = _mm_add_pd(x3, _mm_msubadd_pd(h1_real, q3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
		x3 = _mm_add_pd(x3, _mm_addsub_pd( _mm_mul_pd(h1_real, q3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif

		h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+i)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h2_imag = _mm_xor_pd(h2_imag, sign);
#endif

		tmp1 = _mm_mul_pd(h2_imag, q1);
#ifdef __ELPA_USE_FMA__
		y1 = _mm_add_pd(y1, _mm_msubadd_pd(h2_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h2_imag, q2);
#ifdef __ELPA_USE_FMA__
		y2 = _mm_add_pd(y2, _mm_msubadd_pd(h2_real, q2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		y2 = _mm_add_pd(y2, _mm_addsub_pd( _mm_mul_pd(h2_real, q2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
		tmp3 = _mm_mul_pd(h2_imag, q3);
#ifdef __ELPA_USE_FMA__
		y3 = _mm_add_pd(y3, _mm_msubadd_pd(h2_real, q3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
		y3 = _mm_add_pd(y3, _mm_addsub_pd( _mm_mul_pd(h2_real, q3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif
	}

	h1_real = _mm_loaddup_pd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[((nb-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h1_imag = _mm_xor_pd(h1_imag, sign);
#endif

	q1 = _mm_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm_load_pd(&q_dbl[(2*nb*ldq)+2]);
	q3 = _mm_load_pd(&q_dbl[(2*nb*ldq)+4]);

	tmp1 = _mm_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm_add_pd(x1, _mm_msubadd_pd(h1_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	x1 = _mm_add_pd(x1, _mm_addsub_pd( _mm_mul_pd(h1_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm_add_pd(x2, _mm_msubadd_pd(h1_real, q2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	x2 = _mm_add_pd(x2, _mm_addsub_pd( _mm_mul_pd(h1_real, q2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h1_imag, q3);
#ifdef __ELPA_USE_FMA__
	x3 = _mm_add_pd(x3, _mm_msubadd_pd(h1_real, q3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	x3 = _mm_add_pd(x3, _mm_addsub_pd( _mm_mul_pd(h1_real, q3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif

	h1_real = _mm_loaddup_pd(&hh_dbl[0]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[1]);
	h1_real = _mm_xor_pd(h1_real, sign);
	h1_imag = _mm_xor_pd(h1_imag, sign);

	tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	x1 = _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	tmp2 = _mm_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm_maddsub_pd(h1_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#else
	x2 = _mm_addsub_pd( _mm_mul_pd(h1_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#endif
	tmp3 = _mm_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
	x3 = _mm_maddsub_pd(h1_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1)));
#else
	x3 = _mm_addsub_pd( _mm_mul_pd(h1_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1)));
#endif

	h1_real = _mm_loaddup_pd(&hh_dbl[ldh*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[(ldh*2)+1]);
	h2_real = _mm_loaddup_pd(&hh_dbl[ldh*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[(ldh*2)+1]);

	h1_real = _mm_xor_pd(h1_real, sign);
	h1_imag = _mm_xor_pd(h1_imag, sign);
	h2_real = _mm_xor_pd(h2_real, sign);
	h2_imag = _mm_xor_pd(h2_imag, sign);

	tmp2 = _mm_loadu_pd(s_dbl);
	tmp1 = _mm_mul_pd(h2_imag, tmp2);
#ifdef __ELPA_USE_FMA__
	tmp2 = _mm_maddsub_pd(h2_real, tmp2, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	tmp2 = _mm_addsub_pd( _mm_mul_pd(h2_real, tmp2), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	_mm_storeu_pd(s_dbl, tmp2);
	h2_real = _mm_loaddup_pd(&s_dbl[0]);
	h2_imag = _mm_loaddup_pd(&s_dbl[1]);

	tmp1 = _mm_mul_pd(h1_imag, y1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_maddsub_pd(h1_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	y1 = _mm_addsub_pd( _mm_mul_pd(h1_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	tmp2 = _mm_mul_pd(h1_imag, y2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm_maddsub_pd(h1_real, y2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#else
	y2 = _mm_addsub_pd( _mm_mul_pd(h1_real, y2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#endif
	tmp3 = _mm_mul_pd(h1_imag, y3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm_maddsub_pd(h1_real, y3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1)));
#else
	y3 = _mm_addsub_pd( _mm_mul_pd(h1_real, y3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1)));
#endif

	tmp1 = _mm_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_add_pd(y1, _mm_maddsub_pd(h2_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm_add_pd(y2, _mm_maddsub_pd(h2_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	y2 = _mm_add_pd(y2, _mm_addsub_pd( _mm_mul_pd(h2_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h2_imag, x3);
#ifdef __ELPA_USE_FMA__
	y3 = _mm_add_pd(y3, _mm_maddsub_pd(h2_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	y3 = _mm_add_pd(y3, _mm_addsub_pd( _mm_mul_pd(h2_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif

	q1 = _mm_load_pd(&q_dbl[0]);
	q2 = _mm_load_pd(&q_dbl[2]);
	q3 = _mm_load_pd(&q_dbl[4]);

	q1 = _mm_add_pd(q1, y1);
	q2 = _mm_add_pd(q2, y2);
	q3 = _mm_add_pd(q3, y3);

	_mm_store_pd(&q_dbl[0], q1);
	_mm_store_pd(&q_dbl[2], q2);
	_mm_store_pd(&q_dbl[4], q3);

	h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+1)*2)+1]);

	q1 = _mm_load_pd(&q_dbl[(ldq*2)+0]);
	q2 = _mm_load_pd(&q_dbl[(ldq*2)+2]);
	q3 = _mm_load_pd(&q_dbl[(ldq*2)+4]);

	q1 = _mm_add_pd(q1, x1);
	q2 = _mm_add_pd(q2, x2);
	q3 = _mm_add_pd(q3, x3);

	tmp1 = _mm_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm_add_pd(q1, _mm_maddsub_pd(h2_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h2_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h2_imag, y2);
#ifdef __ELPA_USE_FMA__
	q2 = _mm_add_pd(q2, _mm_maddsub_pd(h2_real, y2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h2_real, y2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h2_imag, y3);
#ifdef __ELPA_USE_FMA__
	q3 = _mm_add_pd(q3, _mm_maddsub_pd(h2_real, y3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	q3 = _mm_add_pd(q3, _mm_addsub_pd( _mm_mul_pd(h2_real, y3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif

	_mm_store_pd(&q_dbl[(ldq*2)+0], q1);
	_mm_store_pd(&q_dbl[(ldq*2)+2], q2);
	_mm_store_pd(&q_dbl[(ldq*2)+4], q3);

	for (i = 2; i < nb; i++)
	{
		q1 = _mm_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm_load_pd(&q_dbl[(2*i*ldq)+2]);
		q3 = _mm_load_pd(&q_dbl[(2*i*ldq)+4]);

		h1_real = _mm_loaddup_pd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm_loaddup_pd(&hh_dbl[((i-1)*2)+1]);

		tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm_add_pd(q1, _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm_add_pd(q2, _mm_maddsub_pd(h1_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h1_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
		tmp3 = _mm_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
		q3 = _mm_add_pd(q3, _mm_maddsub_pd(h1_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
		q3 = _mm_add_pd(q3, _mm_addsub_pd( _mm_mul_pd(h1_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif

		h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+i)*2)+1]);

		tmp1 = _mm_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm_add_pd(q1, _mm_maddsub_pd(h2_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h2_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h2_imag, y2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm_add_pd(q2, _mm_maddsub_pd(h2_real, y2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h2_real, y2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
		tmp3 = _mm_mul_pd(h2_imag, y3);
#ifdef __ELPA_USE_FMA__
		q3 = _mm_add_pd(q3, _mm_maddsub_pd(h2_real, y3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
		q3 = _mm_add_pd(q3, _mm_addsub_pd( _mm_mul_pd(h2_real, y3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif

		_mm_store_pd(&q_dbl[(2*i*ldq)+0], q1);
		_mm_store_pd(&q_dbl[(2*i*ldq)+2], q2);
		_mm_store_pd(&q_dbl[(2*i*ldq)+4], q3);
	}

	h1_real = _mm_loaddup_pd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[((nb-1)*2)+1]);

	q1 = _mm_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm_load_pd(&q_dbl[(2*nb*ldq)+2]);
	q3 = _mm_load_pd(&q_dbl[(2*nb*ldq)+4]);

	tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm_add_pd(q1, _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	q2 = _mm_add_pd(q2, _mm_maddsub_pd(h1_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h1_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	tmp3 = _mm_mul_pd(h1_imag, x3);
#ifdef __ELPA_USE_FMA__
	q3 = _mm_add_pd(q3, _mm_maddsub_pd(h1_real, x3, _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#else
	q3 = _mm_add_pd(q3, _mm_addsub_pd( _mm_mul_pd(h1_real, x3), _mm_shuffle_pd(tmp3, tmp3, _MM_SHUFFLE2(0,1))));
#endif

	_mm_store_pd(&q_dbl[(2*nb*ldq)+0], q1);
	_mm_store_pd(&q_dbl[(2*nb*ldq)+2], q2);
	_mm_store_pd(&q_dbl[(2*nb*ldq)+4], q3);
}

static __forceinline void hh_trafo_complex_kernel_2_SSE_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s)
{
	double* q_dbl = (double*)q;
	double* hh_dbl = (double*)hh;
	double* s_dbl = (double*)(&s);

	__m128d x1, x2;
	__m128d y1, y2;
	__m128d q1, q2;
	__m128d h1_real, h1_imag, h2_real, h2_imag;
	__m128d tmp1, tmp2;
	int i=0;

	__m128d sign = (__m128d)_mm_set_epi64x(0x8000000000000000, 0x8000000000000000);

	x1 = _mm_load_pd(&q_dbl[(2*ldq)+0]);
	x2 = _mm_load_pd(&q_dbl[(2*ldq)+2]);

	h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h2_imag = _mm_xor_pd(h2_imag, sign);
#endif

	y1 = _mm_load_pd(&q_dbl[0]);
	y2 = _mm_load_pd(&q_dbl[2]);

	tmp1 = _mm_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_add_pd(y1, _mm_msubadd_pd(h2_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm_add_pd(y2, _mm_msubadd_pd(h2_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	y2 = _mm_add_pd(y2, _mm_addsub_pd( _mm_mul_pd(h2_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif

	for (i = 2; i < nb; i++)
	{
		q1 = _mm_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm_load_pd(&q_dbl[(2*i*ldq)+2]);

		h1_real = _mm_loaddup_pd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm_loaddup_pd(&hh_dbl[((i-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h1_imag = _mm_xor_pd(h1_imag, sign);
#endif

		tmp1 = _mm_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
		x1 = _mm_add_pd(x1, _mm_msubadd_pd(h1_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		x1 = _mm_add_pd(x1, _mm_addsub_pd( _mm_mul_pd(h1_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
		x2 = _mm_add_pd(x2, _mm_msubadd_pd(h1_real, q2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		x2 = _mm_add_pd(x2, _mm_addsub_pd( _mm_mul_pd(h1_real, q2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif

		h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+i)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h2_imag = _mm_xor_pd(h2_imag, sign);
#endif

		tmp1 = _mm_mul_pd(h2_imag, q1);
#ifdef __ELPA_USE_FMA__
		y1 = _mm_add_pd(y1, _mm_msubadd_pd(h2_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h2_imag, q2);
#ifdef __ELPA_USE_FMA__
		y2 = _mm_add_pd(y2, _mm_msubadd_pd(h2_real, q2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		y2 = _mm_add_pd(y2, _mm_addsub_pd( _mm_mul_pd(h2_real, q2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif
	}

	h1_real = _mm_loaddup_pd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[((nb-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h1_imag = _mm_xor_pd(h1_imag, sign);
#endif

	q1 = _mm_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm_load_pd(&q_dbl[(2*nb*ldq)+2]);

	tmp1 = _mm_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm_add_pd(x1, _mm_msubadd_pd(h1_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	x1 = _mm_add_pd(x1, _mm_addsub_pd( _mm_mul_pd(h1_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h1_imag, q2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm_add_pd(x2, _mm_msubadd_pd(h1_real, q2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	x2 = _mm_add_pd(x2, _mm_addsub_pd( _mm_mul_pd(h1_real, q2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif

	h1_real = _mm_loaddup_pd(&hh_dbl[0]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[1]);
	h1_real = _mm_xor_pd(h1_real, sign);
	h1_imag = _mm_xor_pd(h1_imag, sign);

	tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	x1 = _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	tmp2 = _mm_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	x2 = _mm_maddsub_pd(h1_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#else
	x2 = _mm_addsub_pd( _mm_mul_pd(h1_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#endif

	h1_real = _mm_loaddup_pd(&hh_dbl[ldh*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[(ldh*2)+1]);
	h2_real = _mm_loaddup_pd(&hh_dbl[ldh*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[(ldh*2)+1]);

	h1_real = _mm_xor_pd(h1_real, sign);
	h1_imag = _mm_xor_pd(h1_imag, sign);
	h2_real = _mm_xor_pd(h2_real, sign);
	h2_imag = _mm_xor_pd(h2_imag, sign);

	tmp2 = _mm_loadu_pd(s_dbl);
	tmp1 = _mm_mul_pd(h2_imag, tmp2);
#ifdef __ELPA_USE_FMA__
	tmp2 = _mm_maddsub_pd(h2_real, tmp2, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	tmp2 = _mm_addsub_pd( _mm_mul_pd(h2_real, tmp2), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	_mm_storeu_pd(s_dbl, tmp2);
	h2_real = _mm_loaddup_pd(&s_dbl[0]);
	h2_imag = _mm_loaddup_pd(&s_dbl[1]);

	tmp1 = _mm_mul_pd(h1_imag, y1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_maddsub_pd(h1_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	y1 = _mm_addsub_pd( _mm_mul_pd(h1_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	tmp2 = _mm_mul_pd(h1_imag, y2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm_maddsub_pd(h1_real, y2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#else
	y2 = _mm_addsub_pd( _mm_mul_pd(h1_real, y2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1)));
#endif

	tmp1 = _mm_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_add_pd(y1, _mm_maddsub_pd(h2_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h2_imag, x2);
#ifdef __ELPA_USE_FMA__
	y2 = _mm_add_pd(y2, _mm_maddsub_pd(h2_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	y2 = _mm_add_pd(y2, _mm_addsub_pd( _mm_mul_pd(h2_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif

	q1 = _mm_load_pd(&q_dbl[0]);
	q2 = _mm_load_pd(&q_dbl[2]);

	q1 = _mm_add_pd(q1, y1);
	q2 = _mm_add_pd(q2, y2);

	_mm_store_pd(&q_dbl[0], q1);
	_mm_store_pd(&q_dbl[2], q2);

	h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+1)*2)+1]);

	q1 = _mm_load_pd(&q_dbl[(ldq*2)+0]);
	q2 = _mm_load_pd(&q_dbl[(ldq*2)+2]);

	q1 = _mm_add_pd(q1, x1);
	q2 = _mm_add_pd(q2, x2);

	tmp1 = _mm_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm_add_pd(q1, _mm_maddsub_pd(h2_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h2_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h2_imag, y2);
#ifdef __ELPA_USE_FMA__
	q2 = _mm_add_pd(q2, _mm_maddsub_pd(h2_real, y2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h2_real, y2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif

	_mm_store_pd(&q_dbl[(ldq*2)+0], q1);
	_mm_store_pd(&q_dbl[(ldq*2)+2], q2);

	for (i = 2; i < nb; i++)
	{
		q1 = _mm_load_pd(&q_dbl[(2*i*ldq)+0]);
		q2 = _mm_load_pd(&q_dbl[(2*i*ldq)+2]);

		h1_real = _mm_loaddup_pd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm_loaddup_pd(&hh_dbl[((i-1)*2)+1]);

		tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm_add_pd(q1, _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm_add_pd(q2, _mm_maddsub_pd(h1_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h1_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif

		h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+i)*2)+1]);

		tmp1 = _mm_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm_add_pd(q1, _mm_maddsub_pd(h2_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h2_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
		tmp2 = _mm_mul_pd(h2_imag, y2);
#ifdef __ELPA_USE_FMA__
		q2 = _mm_add_pd(q2, _mm_maddsub_pd(h2_real, y2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
		q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h2_real, y2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif

		_mm_store_pd(&q_dbl[(2*i*ldq)+0], q1);
		_mm_store_pd(&q_dbl[(2*i*ldq)+2], q2);
	}

	h1_real = _mm_loaddup_pd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[((nb-1)*2)+1]);

	q1 = _mm_load_pd(&q_dbl[(2*nb*ldq)+0]);
	q2 = _mm_load_pd(&q_dbl[(2*nb*ldq)+2]);

	tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm_add_pd(q1, _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	tmp2 = _mm_mul_pd(h1_imag, x2);
#ifdef __ELPA_USE_FMA__
	q2 = _mm_add_pd(q2, _mm_maddsub_pd(h1_real, x2, _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#else
	q2 = _mm_add_pd(q2, _mm_addsub_pd( _mm_mul_pd(h1_real, x2), _mm_shuffle_pd(tmp2, tmp2, _MM_SHUFFLE2(0,1))));
#endif

	_mm_store_pd(&q_dbl[(2*nb*ldq)+0], q1);
	_mm_store_pd(&q_dbl[(2*nb*ldq)+2], q2);
}

static __forceinline void hh_trafo_complex_kernel_1_SSE_2hv(std::complex<double>* q, std::complex<double>* hh, int nb, int ldq, int ldh, std::complex<double> s)
{
	double* q_dbl = (double*)q;
	double* hh_dbl = (double*)hh;
	double* s_dbl = (double*)(&s);

	__m128d x1;
	__m128d y1;
	__m128d q1;
	__m128d h1_real, h1_imag, h2_real, h2_imag;
	__m128d tmp1;
	int i=0;

	__m128d sign = (__m128d)_mm_set_epi64x(0x8000000000000000, 0x8000000000000000);

	x1 = _mm_load_pd(&q_dbl[(2*ldq)+0]);

	h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h2_imag = _mm_xor_pd(h2_imag, sign);
#endif

	y1 = _mm_load_pd(&q_dbl[0]);

	tmp1 = _mm_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_add_pd(y1, _mm_msubadd_pd(h2_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif

	for (i = 2; i < nb; i++)
	{
		q1 = _mm_load_pd(&q_dbl[(2*i*ldq)+0]);

		h1_real = _mm_loaddup_pd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm_loaddup_pd(&hh_dbl[((i-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h1_imag = _mm_xor_pd(h1_imag, sign);
#endif

		tmp1 = _mm_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
		x1 = _mm_add_pd(x1, _mm_msubadd_pd(h1_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		x1 = _mm_add_pd(x1, _mm_addsub_pd( _mm_mul_pd(h1_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif

		h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+i)*2)+1]);
#ifndef __ELPA_USE_FMA__
		// conjugate
		h2_imag = _mm_xor_pd(h2_imag, sign);
#endif

		tmp1 = _mm_mul_pd(h2_imag, q1);
#ifdef __ELPA_USE_FMA__
		y1 = _mm_add_pd(y1, _mm_msubadd_pd(h2_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif
	}

	h1_real = _mm_loaddup_pd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[((nb-1)*2)+1]);
#ifndef __ELPA_USE_FMA__
	// conjugate
	h1_imag = _mm_xor_pd(h1_imag, sign);
#endif

	q1 = _mm_load_pd(&q_dbl[(2*nb*ldq)+0]);

	tmp1 = _mm_mul_pd(h1_imag, q1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm_add_pd(x1, _mm_msubadd_pd(h1_real, q1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	x1 = _mm_add_pd(x1, _mm_addsub_pd( _mm_mul_pd(h1_real, q1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif

	h1_real = _mm_loaddup_pd(&hh_dbl[0]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[1]);
	h1_real = _mm_xor_pd(h1_real, sign);
	h1_imag = _mm_xor_pd(h1_imag, sign);

	tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	x1 = _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	x1 = _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif

	h1_real = _mm_loaddup_pd(&hh_dbl[ldh*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[(ldh*2)+1]);
	h2_real = _mm_loaddup_pd(&hh_dbl[ldh*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[(ldh*2)+1]);

	h1_real = _mm_xor_pd(h1_real, sign);
	h1_imag = _mm_xor_pd(h1_imag, sign);
	h2_real = _mm_xor_pd(h2_real, sign);
	h2_imag = _mm_xor_pd(h2_imag, sign);

	__m128d tmp2 = _mm_loadu_pd(s_dbl);
	tmp1 = _mm_mul_pd(h2_imag, tmp2);
#ifdef __ELPA_USE_FMA__
	tmp2 = _mm_maddsub_pd(h2_real, tmp2, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	tmp2 = _mm_addsub_pd( _mm_mul_pd(h2_real, tmp2), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif
	_mm_storeu_pd(s_dbl, tmp2);
	h2_real = _mm_loaddup_pd(&s_dbl[0]);
	h2_imag = _mm_loaddup_pd(&s_dbl[1]);

	tmp1 = _mm_mul_pd(h1_imag, y1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_maddsub_pd(h1_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#else
	y1 = _mm_addsub_pd( _mm_mul_pd(h1_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1)));
#endif

	tmp1 = _mm_mul_pd(h2_imag, x1);
#ifdef __ELPA_USE_FMA__
	y1 = _mm_add_pd(y1, _mm_maddsub_pd(h2_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	y1 = _mm_add_pd(y1, _mm_addsub_pd( _mm_mul_pd(h2_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif

	q1 = _mm_load_pd(&q_dbl[0]);

	q1 = _mm_add_pd(q1, y1);

	_mm_store_pd(&q_dbl[0], q1);

	h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+1)*2]);
	h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+1)*2)+1]);

	q1 = _mm_load_pd(&q_dbl[(ldq*2)+0]);

	q1 = _mm_add_pd(q1, x1);

	tmp1 = _mm_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm_add_pd(q1, _mm_maddsub_pd(h2_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h2_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif

	_mm_store_pd(&q_dbl[(ldq*2)+0], q1);

	for (i = 2; i < nb; i++)
	{
		q1 = _mm_load_pd(&q_dbl[(2*i*ldq)+0]);

		h1_real = _mm_loaddup_pd(&hh_dbl[(i-1)*2]);
		h1_imag = _mm_loaddup_pd(&hh_dbl[((i-1)*2)+1]);

		tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm_add_pd(q1, _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif

		h2_real = _mm_loaddup_pd(&hh_dbl[(ldh+i)*2]);
		h2_imag = _mm_loaddup_pd(&hh_dbl[((ldh+i)*2)+1]);

		tmp1 = _mm_mul_pd(h2_imag, y1);
#ifdef __ELPA_USE_FMA__
		q1 = _mm_add_pd(q1, _mm_maddsub_pd(h2_real, y1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
		q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h2_real, y1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif

		_mm_store_pd(&q_dbl[(2*i*ldq)+0], q1);
	}

	h1_real = _mm_loaddup_pd(&hh_dbl[(nb-1)*2]);
	h1_imag = _mm_loaddup_pd(&hh_dbl[((nb-1)*2)+1]);

	q1 = _mm_load_pd(&q_dbl[(2*nb*ldq)+0]);

	tmp1 = _mm_mul_pd(h1_imag, x1);
#ifdef __ELPA_USE_FMA__
	q1 = _mm_add_pd(q1, _mm_maddsub_pd(h1_real, x1, _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#else
	q1 = _mm_add_pd(q1, _mm_addsub_pd( _mm_mul_pd(h1_real, x1), _mm_shuffle_pd(tmp1, tmp1, _MM_SHUFFLE2(0,1))));
#endif

	_mm_store_pd(&q_dbl[(2*nb*ldq)+0], q1);
}
#endif
} // extern C
