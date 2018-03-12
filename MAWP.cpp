
#include "stdafx.h"
#include <memory.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <cstdio>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <time.h>
#include <conio.h>
#include <sys/timeb.h>
#include <windows.h>
#include <iostream>
using namespace std;

#define P_MAX 200	//maximum pattern length
#define V_MAX 20000 //maximum size of pointer and shift arrays

#define SIGMA 8		//alphabet size
#define SIGMA2 64
#define SIGMA3 512
#define SIGMA4 4096
#define SIGMA5 32768
#define SIGMA6 262144
#define SIGMA7 2097152
#define SIGMA8 16777216
#define SIGMA9 134217728
#define LOG_SIGMA 3 //logarithm of alphabet size

/*#define SIGMA 6		//alphabet size
#define SIGMA2 36
#define SIGMA3 216
#define SIGMA4 1296
#define SIGMA5 7776 
#define SIGMA6 46656
#define SIGMA7 279936
#define SIGMA8 1679616
#define SIGMA9 10077696

#define SIGMA 4		//alphabet size
#define SIGMA2 16
#define SIGMA3 64
#define SIGMA4 256
#define SIGMA5 1024
#define SIGMA6 4096
#define SIGMA7 16384
#define SIGMA8 65536
#define SIGMA9 262144
#define LOG_SIGMA 2 //logarithm of alphabet size
#define SIGMA_DOUBLE 8
#define SIGMA2_DOUBLE 32*/

const int TOTAL = 10000000 + 5 * P_MAX;
unsigned char T[TOTAL], T1[TOTAL], P[P_MAX], P1[P_MAX];
int N = TOTAL - 5 * P_MAX, ITER = 200, m = 5;

FILE * f;
LARGE_INTEGER start, _end, freq, _freq, prep_start, prep_end;
double u;
int nm2, glob = 0;

long long sum_maw22, sum_maw23, sum_maw24, sum_maw32, sum_maw33;
long long sum_prep22, sum_prep23, sum_prep24, sum_prep32, sum_prep33;
int maw22, maw23, maw24, maw32, maw33;
int maw22c, maw23c, maw24c, maw32c, maw33c;

long long sum_maw22p, sum_maw23p, sum_maw24p, sum_maw32p, sum_maw33p, sum_maw22p_mod;
long long sum_prep22p, sum_prep23p, sum_prep24p, sum_prep32p, sum_prep33p, sum_prep22p_mod;
int maw22p, maw23p, maw24p, maw32p, maw33p, maw22p_mod;
int maw22pc, maw23pc, maw24pc, maw32pc, maw33pc, maw22p_modc;

//a bytes pointed by c are repeated until c[0..b] is filled
void mem_fill(int a, int b, unsigned char* c) {
	int i;
	for (i = a; i <= (b >> 1); i <<= 1)
		memcpy(c + i, c, i);
	memcpy(c + i, c, b - i);
}

template <class T, class U>
void copy_value(T* pointer, U value, int length, int type_length = 1) {
	*pointer = value;
	mem_fill(type_length, length, (unsigned char*)pointer);
}

// The MAW22 algorithm
int MAW22(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&prep_start);

	int mp1 = m + 1, mm2 = m - 2, mm1 = m - 1, m2 = 2 * m, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, k, r, pos, count = 0;
	int D[P_MAX];
	unsigned char* M22 = (unsigned char *)malloc(SIGMA4); // MAW22 search table

	// Preprocessing
	// Build the BMH shift table D
	for (int i = 0; i < SIGMA; i++)
		D[i] = m;
	for (int i = 0; i < mm1; i++)
		D[x[i]] = mm1 - i;

	// Fill the M22 search table
	*M22 = m2;
	mem_fill(1, SIGMA, M22);

	*(M22 + x[0]) = m2m1;
	mem_fill(SIGMA, SIGMA2, M22);

	for (int k = 0; k < mm1; k++)
		*(M22 + x[k] * SIGMA + x[k + 1]) = m2m2 - k;
	mem_fill(SIGMA2, SIGMA3, M22);

	copy_value(M22 + x[0] * SIGMA2 + x[mm1] * SIGMA, mm1, SIGMA);
	mem_fill(SIGMA3, SIGMA4, M22);

	for (int k = 0; k < mm1; k++)
		copy_value(M22 + x[k] * SIGMA3 + x[k + 1] * SIGMA2, mm2 - k, SIGMA2);

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep22 += u;

	QueryPerformanceCounter(&start);

	//Search
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M22 + y[pos] * SIGMA3 + y[pos + 1] * SIGMA2 + y[pos + m] * SIGMA + y[pos + mp1]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < mm2 && y[pos + k - mm2] == x[k]; k++);
			if (k == mm2) {
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
		maw22c++;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw22 += u;

	free(M22);
	return count;
}

// The MAW23 algorithm
int MAW23(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 3) return -1;
	QueryPerformanceCounter(&prep_start);

	int mp1 = m + 1, mm2 = m - 2, mm1 = m - 1, m2 = 2 * m, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3, mm3 = m - 3, mp2 = m + 2, k, r, pos, count = 0;
	int D[P_MAX];
	unsigned char* M23 = (unsigned char*)malloc(SIGMA6); // MAW23 search table

	// Preprocessing
	// Build the BMH shift table D
	for (int i = 0; i < SIGMA; i++)
		D[i] = m;
	for (int i = 0; i < mm1; i++)
		D[x[i]] = mm1 - i;

	// Fill the M23 search table
	*M23 = m2;
	mem_fill(1, SIGMA, M23);

	*(M23 + x[0]) = m2m1;
	mem_fill(SIGMA, SIGMA2, M23);

	*(M23 + x[0] * SIGMA + x[1]) = m2m2;
	mem_fill(SIGMA2, SIGMA3, M23);

	for (int k = 0; k < mm2; k++)
		*(M23 + x[k] * SIGMA2 + x[k + 1] * SIGMA + x[k + 2]) = m2m3 - k;
	mem_fill(SIGMA3, SIGMA4, M23);

	copy_value(M23 + x[0] * SIGMA3 + x[mm2] * SIGMA2 + x[mm1] * SIGMA, mm1, SIGMA);
	mem_fill(SIGMA4, SIGMA5, M23);

	copy_value(M23 + x[0] * SIGMA4 + x[1] * SIGMA3 + x[mm1] * SIGMA2, mm2, SIGMA2);
	mem_fill(SIGMA5, SIGMA6, M23);

	for (int k = 0; k < mm2; k++)
		copy_value(M23 + x[k] * SIGMA5 + x[k + 1] * SIGMA4 + x[k + 2] * SIGMA3, mm3 - k, SIGMA3);

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep23 += u;

	QueryPerformanceCounter(&start);

	//Search
	pos = mm3;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M23 + y[pos] * SIGMA5 + y[pos + 1] * SIGMA4 + y[pos + 2] * SIGMA3 + y[pos + m] * SIGMA2 + y[pos + mp1] * SIGMA + y[pos + mp2]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < mm3 && y[pos + k - mm3] == x[k]; k++);
			if (k == mm3) {
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw23 += u;

	free(M23);
	return count;
}

// The MAW24 algorithm
int MAW24(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 4) return -1;
	QueryPerformanceCounter(&prep_start);

	int mp1 = m + 1, mp2 = m + 2, mp3 = m + 3, mm2 = m - 2, mm1 = m - 1, mm3 = m - 3, mm4 = m - 4,
		m2 = 2 * m, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3, m2m4 = 2 * m - 4,
		k, r, pos, count = 0;
	int D[P_MAX];
	unsigned char* M24 = (unsigned char *)malloc(SIGMA8); // MAW24 search table

	// Preprocessing
	// Build the BMH shift table D
	for (int i = 0; i < SIGMA; i++)
		D[i] = m;
	for (int i = 0; i < mm1; i++)
		D[x[i]] = mm1 - i;

	// Fill the M24 search table
	*M24 = m2;
	mem_fill(1, SIGMA, M24);

	*(M24 + x[0]) = m2m1;
	mem_fill(SIGMA, SIGMA2, M24);

	*(M24 + x[0] * SIGMA + x[1]) = m2m2;
	mem_fill(SIGMA2, SIGMA3, M24);

	*(M24 + x[0] * SIGMA2 + x[1] * SIGMA + x[2]) = m2m3;
	mem_fill(SIGMA3, SIGMA4, M24);

	for (int k = 0; k < mm3; k++)
		*(M24 + x[k] * SIGMA3 + x[k + 1] * SIGMA2 + x[k + 2] * SIGMA + x[k + 3]) = m2m4 - k;
	mem_fill(SIGMA4, SIGMA5, M24);

	copy_value(M24 + x[0] * SIGMA4 + x[mm3] * SIGMA3 + x[mm2] * SIGMA2 + x[mm1] * SIGMA, mm1, SIGMA);
	mem_fill(SIGMA5, SIGMA6, M24);

	copy_value(M24 + x[0] * SIGMA5 + x[1] * SIGMA4 + x[mm2] * SIGMA3 + x[mm1] * SIGMA2, mm2, SIGMA2);
	mem_fill(SIGMA6, SIGMA7, M24);

	copy_value(M24 + x[0] * SIGMA6 + x[1] * SIGMA5 + x[2] * SIGMA4 + x[mm1] * SIGMA3, mm3, SIGMA3);
	mem_fill(SIGMA7, SIGMA8, M24);

	for (int k = 0; k < mm3; k++)
		copy_value(M24 + x[k] * SIGMA7 + x[k + 1] * SIGMA6 + x[k + 2] * SIGMA5 + x[k + 3] * SIGMA4, mm4 - k, SIGMA4);

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep24 += u;

	QueryPerformanceCounter(&start);

	//Search
	pos = mm4;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M24 + y[pos] * SIGMA7 + y[pos + 1] * SIGMA6 + y[pos + 2] * SIGMA5 + y[pos + 3] * SIGMA4 +
			y[pos + m] * SIGMA3 + y[pos + mp1] * SIGMA2 + y[pos + mp2] * SIGMA + y[pos + mp3]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < mm4 && y[pos + k - mm4] == x[k]; k++);
			if (k == mm4) {
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw24 += u;

	free(M24);
	return count;
}

// The MAW32 algorithm
int MAW32(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&prep_start);

	int mp1 = m + 1, mm2 = m - 2, mm1 = m - 1, m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m3 = m * 3,
		m3m1 = 3 * m - 1, m3m2 = 3 * m - 2, k, r, pos, count = 0;
	int D[P_MAX];
	unsigned char* M32 = (unsigned char *)malloc(SIGMA6); // MAW32 search table

	// Preprocessing
	// Build the BMH shift table D
	for (int i = 0; i < SIGMA; i++)
		D[i] = m;
	for (int i = 0; i < mm1; i++)
		D[x[i]] = mm1 - i;

	// Fill the M32 search table
	*M32 = m3;
	mem_fill(1, SIGMA, M32);

	*(M32 + x[0]) = m3m1;
	mem_fill(SIGMA, SIGMA2, M32);

	for (int k = 0; k < mm1; k++)
		*(M32 + x[k] * SIGMA + x[k + 1]) = m3m2 - k;
	mem_fill(SIGMA2, SIGMA3, M32);

	copy_value(M32 + x[0] * SIGMA2 + x[mm1] * SIGMA, m2m1, SIGMA);
	mem_fill(SIGMA3, SIGMA4, M32);

	for (int k = 0; k < mm1; k++)
		copy_value(M32 + x[k] * SIGMA3 + x[k + 1] * SIGMA2, m2m2 - k, SIGMA2);
	mem_fill(SIGMA4, SIGMA5, M32);

	copy_value(M32 + x[0] * SIGMA4 + x[mm1] * SIGMA3, mm1, SIGMA3);
	mem_fill(SIGMA5, SIGMA6, M32);

	for (int k = 0; k < mm1; k++)
		copy_value(M32 + x[k] * SIGMA5 + x[k + 1] * SIGMA4, mm2 - k, SIGMA4);

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep32 += u;

	QueryPerformanceCounter(&start);

	//Search
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M32 + y[pos] * SIGMA5 + y[pos + 1] * SIGMA4 + y[pos + m] * SIGMA3 + y[pos + mp1] * SIGMA2 + y[pos + m2] * SIGMA + y[pos + m2p1]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < mm2 && y[pos + k - mm2] == x[k]; k++);
			if (k == mm2) {
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw32 += u;

	free(M32);
	return count;
}

// The MAW33 algorithm
int MAW33(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 3) return -1;
	QueryPerformanceCounter(&prep_start);

	int mp1 = m + 1, mp2 = m + 2, mm2 = m - 2, mm3 = m - 3, mm1 = m - 1,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2p2 = 2 * m + 2, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3,
		m3 = m * 3, m3m1 = 3 * m - 1, m3m2 = 3 * m - 2, m3m3 = 3 * m - 3,
		k, r, pos, count = 0;
	int D[P_MAX];
	unsigned char* M33 = (unsigned char *)malloc(SIGMA9); // MAW33 search table
	unsigned char* point;

	// Preprocessing
	// Build the BMH shift table D
	for (int i = 0; i < SIGMA; i++)
		D[i] = m;
	for (int i = 0; i < mm1; i++)
		D[x[i]] = mm1 - i;

	// Fill the M33 search table
	*M33 = m3;
	mem_fill(1, SIGMA, M33);

	*(M33 + x[0]) = m3m1;
	mem_fill(SIGMA, SIGMA2, M33);

	*(M33 + x[0] * SIGMA + x[1]) = m3m2;
	mem_fill(SIGMA2, SIGMA3, M33);

	for (int k = 0; k < mm2; k++)
		*(M33 + x[k] * SIGMA2 + x[k + 1] * SIGMA + x[k + 2]) = m3m3 - k;
	mem_fill(SIGMA3, SIGMA4, M33);

	copy_value(M33 + x[0] * SIGMA3 + x[mm2] * SIGMA2 + x[mm1] * SIGMA, m2m1, SIGMA);
	mem_fill(SIGMA4, SIGMA5, M33);

	copy_value(M33 + x[0] * SIGMA4 + x[1] * SIGMA3 + x[mm1] * SIGMA2, m2m2, SIGMA2);
	mem_fill(SIGMA5, SIGMA6, M33);

	for (int k = 0; k < mm2; k++)
		copy_value(M33 + x[k] * SIGMA5 + x[k + 1] * SIGMA4 + x[k + 2] * SIGMA3, m2m3 - k, SIGMA3);
	mem_fill(SIGMA6, SIGMA7, M33);

	copy_value(M33 + x[0] * SIGMA6 + x[mm2] * SIGMA5 + x[mm1] * SIGMA4, mm1, SIGMA4);
	mem_fill(SIGMA7, SIGMA8, M33);

	copy_value(M33 + x[0] * SIGMA7 + x[1] * SIGMA6 + x[mm1] * SIGMA5, mm2, SIGMA5);
	mem_fill(SIGMA8, SIGMA9, M33);

	for (int k = 0; k < mm2; k++)
		copy_value(M33 + x[k] * SIGMA8 + x[k + 1] * SIGMA7 + x[k + 2] * SIGMA6, mm3 - k, SIGMA6);

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep33 += u;

	QueryPerformanceCounter(&start);

	//Search
	pos = mm3;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M33 + y[pos] * SIGMA8 + y[pos + 1] * SIGMA7 + y[pos + 2] * SIGMA6
			+ y[pos + m] * SIGMA5 + y[pos + mp1] * SIGMA4 + y[pos + mp2] * SIGMA3
			+ y[pos + m2] * SIGMA2 + y[pos + m2p1] * SIGMA + y[pos + m2p2]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < mm3 && y[pos + k - mm3] == x[k]; k++);
			if (k == mm3) {
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw33 += u;

	free(M33);
	return count;
}


// The MAW22 algorithm with pointers
int MAW22P(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&prep_start);

	int ***V0[V_MAX], **V1[V_MAX], *V2[V_MAX], V3[V_MAX]; //V3 - shift array; V0, V1, V2 - pointers arrays
	int D[P_MAX], D_[P_MAX], BR_[SIGMA][SIGMA];
	int pos, r, k, count = 0, mp1 = m + 1, mp2 = m + 2, mp21 = 2 * m + 1, mp22 = 2 * m + 2, mm1 = m - 1, mm2 = m - 2, m2 = 2 * m, int_size = sizeof(int);

	//Preprocessing
	// Build BMH shift table D and modified BMH shift table D_
	for (int i = 0; i < SIGMA; i++)
		D[i] = D_[i] = m;
	D_[x[mm1]] = 0;
	for (int i = 0; i < mm1; i++) {
		D[x[i]] = mm1 - i;
		D_[x[i]] = mm1 - i - 1;
	}

	// Build the modified Berry-Ravindran shift table
	BR_[0][0] = m;
	mem_fill(int_size, SIGMA*SIGMA * int_size, (unsigned char*)BR_);
	for (int i = 0; i < SIGMA; i++)
		BR_[i][x[0]] = mm1;
	for (int i = 0; i < mm1; i++)
		BR_[x[i]][x[i + 1]] = m - i - 2;

	//Initial fill V0 with the address V1, V1 with the address V2, V2 with the address V3
	V0[0] = V1;
	mem_fill(int_size, V_MAX, (unsigned char*)V0);
	V1[0] = V2;
	mem_fill(int_size, V_MAX, (unsigned char*)V1);
	V2[0] = V3;
	mem_fill(int_size, V_MAX, (unsigned char*)V2);

	// Filling V0 with pointers to chunks of V1
	for (int i = 0; i < mm1; i++)
		V0[x[i]] = V1 + ((D_[x[i]] + 1) << LOG_SIGMA);

	// Filling V1 with pointers to chunks of V2
	for (int i = 0; i < m; i++)
		V1[(i << LOG_SIGMA) + x[0]] = V2 + (m << LOG_SIGMA);
	for (int i = 1; i < m; i++)
	for (int j = 0; j < SIGMA; j++)
	if (BR_[x[mm1 - i]][j] < m)
		V1[(i << LOG_SIGMA) + j] = V2 + ((BR_[x[mm1 - i]][j] + 1) << LOG_SIGMA);
	else
		V1[(i << LOG_SIGMA) + j] = V2;

	// Filling V2 with pointers to chunks of V3
	for (int i = 0; i < mm1; i++) {
		V2[x[i]] = V3 + ((D_[x[i]] + mp1) << LOG_SIGMA);
		V2[(m << LOG_SIGMA) + x[i]] = V3 + ((D_[x[i]] + mp1) << LOG_SIGMA);
	}
	V2[(m << LOG_SIGMA) + x[mm1]] = V3 + (m << LOG_SIGMA);
	for (int i = 1; i < m; i++) {
		V2[SIGMA*i] = V3 + SIGMA*i;
		mem_fill(int_size, SIGMA * int_size, (unsigned char*)(V2 + SIGMA*i));
	}

	// Filling V3 with shift values
	V3[0] = m2;
	mem_fill(int_size, SIGMA * int_size, (unsigned char*)V3);
	V3[x[0]] = m2 - 1;
	for (int i = 1; i <= m; i++) {
		V3[i*SIGMA] = i - 1;
		mem_fill(int_size, SIGMA * int_size, (unsigned char*)(V3 + i*SIGMA));
	}
	for (int i = m + 1; i < m2; i++)
	for (int j = 0; j < SIGMA; j++)
		V3[(i << LOG_SIGMA) + j] = BR_[x[m2 - i - 1]][j] + m;

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep22p += u;

	QueryPerformanceCounter(&start);

	//Search
	int ***p1, **p2, *p3;
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		p1 = V0[y[pos]];
		p2 = p1[y[pos + 1]];
		p3 = p2[y[pos + m]];
		r = p3[y[pos + mp1]];
		if (!r) {
			for (k = 0; k < m && y[pos - mm2 + k] == x[k]; k++);
			if (k == m) {
				if (pos >= n)
					break;
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
		maw22pc++;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw22p += u;

	return count;
}

// Build BMH shift table D and modified BMH shift table D_
void buildBMHShiftTable(int* D, int* D_, const unsigned char *x, const int& m) {
	int mm1 = m - 1, mm2 = m - 2;
	for (int i = 0; i < SIGMA; i++) {
		D[i] = m;
		D_[i] = mm1;
	}
	for (int i = 0; i < mm1; i++) {
		D[x[i]] = mm1 - i;
		D_[x[i]] = mm2 - i;
	}
}

// Build the modified Berry-Ravindran shift table
void buildBRShiftTable(int(*BR_)[SIGMA], const unsigned char *x, const int& m, const int& int_size) {
	int mm1 = m - 1, mm2 = m - 2;
	BR_[0][0] = m;
	mem_fill(int_size, int_size * SIGMA2, (unsigned char*)BR_);
	for (int i = 0; i < SIGMA; i++)
		BR_[i][x[0]] = mm1;
	for (int i = 0; i < mm1; i++)
		BR_[x[i]][x[i + 1]] = mm2 - i;
}

// Build tripple shift table
template <class T>
void buildTrippleShiftTable(T *Tripple, const unsigned char *x, const int& m, const int& int_size) {
	int mm1 = m - 1, mm2 = m - 2, mm3 = m - 3;
	Tripple[0][0][0] = m;
	mem_fill(int_size, int_size * SIGMA3, (unsigned char*)Tripple);
	for (int i = 0; i < SIGMA; i++)
	for (int j = 0; j < SIGMA; j++)
		Tripple[i][j][x[0]] = mm1;
	for (int j = 0; j < SIGMA; j++)
		Tripple[j][x[0]][x[1]] = mm2;
	for (int i = 0; i < mm2; i++)
		Tripple[x[i]][x[i + 1]][x[i + 2]] = mm3 - i;
}

template <class T, class U>
void fillBeginning(T from, U to, int amount, int length, int type_length) {
	for (int i = 0; i < amount; i++)
		copy_value(from + (i << LOG_SIGMA), to + (i << LOG_SIGMA), length, type_length);
}

template <class T>
void fillBeginningFinal(T from, int amount, int length, int type_length) {
	for (int i = 0; i < amount; i++)
		copy_value(from + (i << LOG_SIGMA), i, length, type_length);
}

template <class T, class U>
void fillFirstLetter(T from, U to, int block_length, int type_length, int letter, int length) {
	copy_value(from, to, block_length, type_length);
	*(from + letter) = to - SIGMA;
	mem_fill(block_length, length, (unsigned char*)from);
}

template <class T, class U>
void fillFirstLetterFinal(T from, U to, int block_length, int type_length, int letter, int length) {
	copy_value(from, to, block_length, type_length);
	*(from + letter) = to - 1;
	mem_fill(block_length, length, (unsigned char*)from);
}

// The MAW22 algorithm with pointers
int MAW22P_mod(unsigned char *x, const int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&prep_start);

	int ***V0[SIGMA], **V1[SIGMA * P_MAX], *V2[SIGMA * (P_MAX + 1)], V3[SIGMA * P_MAX * 2]; //V3 - shift array; V0, V1, V2 - pointers arrays
	int D[P_MAX], D_[P_MAX], BR_[SIGMA][SIGMA];
	int pos, r, k, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mp2 = m + 2, mm1 = m - 1, mm2 = m - 2,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2,
		m_sigma = m * SIGMA, mm1_sigma = mm1 * SIGMA, int_size_sigma = int_size * SIGMA, int_size_sigma_2 = int_size_sigma * 2, int_size_sigma_m = int_size_sigma * m;

	//Preprocessing
	buildBMHShiftTable(D, D_, x, m);
	buildBRShiftTable(BR_, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	copy_value(V0, V1 + (mm1 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V0[x[i]] = V1 + (D_[x[i]] << LOG_SIGMA);

	// Filling V1 with pointers to chunks of V2
	fillFirstLetter(V1, V2 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
	for (int j = 0; j < SIGMA; j++)
		V1[(i << LOG_SIGMA) + j] = V2 + ((BR_[x[mm2 - i]][j]) << LOG_SIGMA);

	// Filling V2 with pointers to chunks of V3
	fillBeginning(V2, V3, mm1, int_size_sigma, int_size);
	copy_value(V2 + mm1_sigma, V3 + SIGMA * m2m1, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V2[(mm1 << LOG_SIGMA) + x[i]] = V3 + ((D_[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V2 + mm1_sigma));
	V2[(mm1 << LOG_SIGMA) + x[mm1]] = V3 + (mm1 << LOG_SIGMA);

	// Filling V3 with shift values
	fillBeginningFinal(V3, m, int_size_sigma, int_size);
	fillFirstLetterFinal(V3 + m_sigma, m2, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = m; i < m2m1; i++)
	for (int j = 0; j < SIGMA; j++)
		V3[(i << LOG_SIGMA) + j] = BR_[x[m2m2 - i]][j] + m;

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep22p_mod += u;

	QueryPerformanceCounter(&start);

	//Search
	int ***p1, **p2, *p3;
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		p1 = V0[y[pos]];
		p2 = p1[y[pos + 1]];
		p3 = p2[y[pos + m]];
		r = p3[y[pos + mp1]];
		if (!r) {
			for (k = 0; k < mm2 && y[pos - mm2 + k] == x[k]; k++);
			if (k == mm2) {
				if (pos >= n)
					break;
				count++;
			}
			pos += D[y[pos + 1]]; // !!!
		}
		else
			pos += r;
		maw22p_modc++;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw22p_mod += u;

	return count;
}

// The MAW23 algorithm with pointers
int MAW23P(unsigned char *x, const int m, unsigned char *y, int n) {
	if (m < 3) return -1;

	QueryPerformanceCounter(&prep_start);

	int *****V0[SIGMA], ****V1[SIGMA * (P_MAX - 1)], ***V2[SIGMA * P_MAX], **V3[SIGMA * (P_MAX + 1)],
		*V4[SIGMA * (2 * P_MAX - 1)], V5[SIGMA * P_MAX * 2]; //V5 - shift array; V0, V1, V2, V3, V4 - pointers arrays
	int D[P_MAX], D_[P_MAX], Dm2[P_MAX], BR_[SIGMA][SIGMA], BRm2[SIGMA][SIGMA], Tripple[SIGMA][SIGMA][SIGMA];
	int pos, r, k, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mp2 = m + 2, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3,
		m3 = m * 3, m3m1 = m * 3 - 1,
		m_sigma = m * SIGMA, mm1_sigma = mm1 * SIGMA, m2_sigma = m * 2 * SIGMA,
		int_size_sigma = int_size * SIGMA, int_size_sigma_2 = int_size_sigma * 2, int_size_sigma_m = int_size_sigma * m, int_size_sigma_m2 = int_size_sigma_m * 2;

	for (int i = 0; i < SIGMA; i++)
		Dm2[i] = mm2;
	for (int i = 0; i < mm2; i++)
		Dm2[x[i]] = mm3 - i;

	BRm2[0][0] = mm1;
	mem_fill(int_size, int_size * SIGMA2, (unsigned char*)BRm2);
	for (int i = 0; i < SIGMA; i++)
		BRm2[i][x[0]] = mm2;
	for (int i = 0; i < mm2; i++)
		BRm2[x[i]][x[i + 1]] = mm3 - i;

	BR_[0][0] = m;
	mem_fill(int_size, int_size * SIGMA2, (unsigned char*)BR_);
	for (int i = 0; i < SIGMA; i++)
		BR_[i][x[0]] = mm1;
	for (int i = 0; i < mm2; i++)
		BR_[x[i]][x[i + 1]] = mm2 - i;

	//Preprocessing
	buildBMHShiftTable(D, D_, x, m);
	//buildBRShiftTable(BR_, x, m, int_size);
	buildTrippleShiftTable(Tripple, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	copy_value(V0, V1 + (mm2 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V0[x[i]] = V1 + (Dm2[x[i]] << LOG_SIGMA);

	// Filling V1 with pointers to chunks of V2
	fillFirstLetter(V1, V2 + mm1_sigma, int_size_sigma, int_size, x[0], int_size_sigma*mm1);
	for (int i = 0; i < mm2; i++)
	for (int j = 0; j < SIGMA; j++)
		V1[(i << LOG_SIGMA) + j] = V2 + ((BRm2[x[mm3 - i]][j]) << LOG_SIGMA);

	// Filling V2 with pointers to chunks of V3
	fillFirstLetter(V2, V3 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	V2[(mm2 << LOG_SIGMA) + x[1]] = V3 + ((BR_[x[0]][x[1]]) << LOG_SIGMA);
	for (int i = 0; i < mm2; i++)
	for (int j = 0; j < SIGMA; j++)
		V2[(i << LOG_SIGMA) + j] = V3 + ((Tripple[x[mm3 - i]][x[mm2 - i]][j]) << LOG_SIGMA);

	// Filling V3 with pointers to chunks of V4
	fillBeginning(V3, V4, mm2, int_size_sigma, int_size);
	copy_value(V3 + (mm2 << LOG_SIGMA), V4 + SIGMA * m2m2, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V3[(mm2 << LOG_SIGMA) + x[i]] = V4 + ((Dm2[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma * 3, (unsigned char*)(V3 + mm2*SIGMA));
	V3[(mm1 << LOG_SIGMA) + x[mm2]] = V4 + (mm1 << LOG_SIGMA);
	for (int i = 0; i < m; i++)
		V3[(mm2 << LOG_SIGMA) + x[i]] = V4 + ((m2m3 - i) << LOG_SIGMA);

	// Filling V4 with pointers to chunks of V5
	fillBeginning(V4, V5, mm1, int_size_sigma, int_size);
	fillFirstLetter(V4 + mm1*SIGMA, V5 + m2m1*SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm2; i++)
	for (int j = 0; j < SIGMA; j++)
		V4[((i + m) << LOG_SIGMA) + j] = V5 + ((BRm2[x[mm3 - i]][j] + m) << LOG_SIGMA);
	for (int j = 0; j < SIGMA; j++)
		V4[(mm1 << LOG_SIGMA) + j] = V5 + ((BRm2[x[mm2]][j] + m) << LOG_SIGMA);
	V4[(mm1 << LOG_SIGMA) + x[mm1]] = V5 + (mm1 << LOG_SIGMA);

	//Filling V5 with shift values
	fillBeginningFinal(V5, m, int_size_sigma, int_size);
	fillFirstLetterFinal(V5 + m_sigma, m2, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = m; i < m2m2; i++)
	for (int j = 0; j < SIGMA; j++)
		V5[(i << LOG_SIGMA) + j] = Tripple[x[m2m3 - i]][x[m2m2 - i]][j] + m;
	V5[(m2m2 << LOG_SIGMA) + x[1]] = BR_[x[0]][x[1]] + m;

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep23p += u;

	QueryPerformanceCounter(&start);

	//Search
	int *****p1, ****p2, ***p3, **p4, *p5;
	pos = mm3;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		p1 = V0[y[pos]];
		p2 = p1[y[pos + 1]];
		p3 = p2[y[pos + 2]];
		p4 = p3[y[pos + m]];
		p5 = p4[y[pos + mp1]];
		r = p5[y[pos + mp2]];
		if (!r) {
			for (k = 0; k < mm3 && y[pos - mm3 + k] == x[k]; k++);
			if (k == mm3) {
				if (pos >= n)
					break;
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw23p += u;

	return count;
}

// The MAW24 algorithm with pointers
int MAW24P(unsigned char *x, const int m, unsigned char *y, int n) {
	if (m < 4) return -1;

	QueryPerformanceCounter(&prep_start);

	int *******V0[SIGMA], ******V1[SIGMA * (P_MAX - 2)], *****V2[SIGMA * (P_MAX - 1)], ****V3[SIGMA * P_MAX], ***V4[SIGMA * (P_MAX + 1)],
		**V5[SIGMA * (2 * P_MAX - 1)], *V6[SIGMA * (2 * P_MAX - 2)], V7[SIGMA * P_MAX * 2]; //V7 - shift array; V0, V1, V2, V3, V4, V5, V6 - pointers arrays
	int D[P_MAX], D_[P_MAX], Dm2[P_MAX], Dm3[P_MAX], BR_[SIGMA][SIGMA], BRm2[SIGMA][SIGMA], BRm3[SIGMA][SIGMA], BRi[SIGMA][SIGMA], 
		Tripple[SIGMA][SIGMA][SIGMA], Quad[SIGMA][SIGMA][SIGMA][SIGMA];
	int pos, r, k, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mp2 = m + 2, mp3 = m + 3, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3, mm4 = m - 4,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3, m2m4 = 2 * m - 4,
		m3 = m * 3, m3m1 = m * 3 - 1,
		m_sigma = m * SIGMA, mm1_sigma = mm1 * SIGMA, m2_sigma = m * 2 * SIGMA,
		int_size_sigma = int_size * SIGMA, int_size_sigma_2 = int_size_sigma * 2, int_size_sigma_m = int_size_sigma * m, int_size_sigma_m2 = int_size_sigma_m * 2;

	for (int i = 0; i < SIGMA; i++)
		Dm3[i] = mm3;
	for (int i = 0; i < mm3; i++)
		Dm3[x[i]] = mm4 - i;

	for (int i = 0; i < SIGMA; i++)
		Dm2[i] = mm2;
	for (int i = 0; i < mm2; i++)
		Dm2[x[i]] = mm3 - i;

	BRm3[0][0] = mm2;
	mem_fill(int_size, int_size * SIGMA2, (unsigned char*)BRm3);
	for (int i = 0; i < SIGMA; i++)
		BRm3[i][x[0]] = mm3;
	for (int i = 0; i < mm3; i++)
		BRm3[x[i]][x[i + 1]] = mm4 - i;

	BRm2[0][0] = mm1;
	mem_fill(int_size, int_size * SIGMA2, (unsigned char*)BRm2);
	for (int i = 0; i < SIGMA; i++)
		BRm2[i][x[0]] = mm2;
	for (int i = 0; i < mm3; i++)
		BRm2[x[i]][x[i + 1]] = mm3 - i;

	BR_[0][0] = m;
	mem_fill(int_size, int_size * SIGMA2, (unsigned char*)BR_);
	for (int i = 0; i < SIGMA; i++)
		BR_[i][x[0]] = mm1;
	for (int i = 0; i < mm1; i++)
		BR_[x[i]][x[i + 1]] = mm2 - i;

	BRi[0][0] = m2m3;
	mem_fill(int_size, int_size * SIGMA2, (unsigned char*)BRi);
	for (int i = 0; i < SIGMA; i++)
	for (int j = 0; j < mm3; j++)
		BRi[i][x[j]] = m2m4 - j;
	for (int i = mm3; i < m; i++)
		BRi[x[i - mm3]][x[i]] = m2m4 - i;

	Tripple[0][0][0] = mm1;
	mem_fill(int_size, int_size * SIGMA3, (unsigned char*)Tripple);
	for (int i = 0; i < SIGMA; i++)
	for (int j = 0; j < SIGMA; j++)
		Tripple[i][j][x[0]] = mm2;
	for (int j = 0; j < SIGMA; j++)
		Tripple[j][x[0]][x[1]] = mm3;
	for (int i = 0; i < mm3; i++)
		Tripple[x[i]][x[i + 1]][x[i + 2]] = mm4 - i;

	Quad[0][0][0][0] = m;
	mem_fill(int_size, int_size * SIGMA4, (unsigned char*)Quad);
	for (int i = 0; i < SIGMA; i++)
	for (int j = 0; j < SIGMA; j++)
	for (int k = 0; k < SIGMA; k++)
		Quad[i][j][k][x[0]] = mm1;
	for (int i = 0; i < SIGMA; i++)
	for (int j = 0; j < SIGMA; j++)
		Quad[i][j][x[0]][x[1]] = mm2;
	for (int i = 0; i < SIGMA; i++)
		Quad[i][x[0]][x[1]][x[2]] = mm3;
	for (int i = 0; i < mm3; i++)
		Quad[x[i]][x[i + 1]][x[i + 2]][x[i + 3]] = mm4 - i;

	//Preprocessing
	buildBMHShiftTable(D, D_, x, m);
	//buildBRShiftTable(BR_, x, m, int_size);
	//buildTrippleShiftTable(Tripple, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	copy_value(V0, V1 + (mm3 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm3; i++)
		V0[x[i]] = V1 + (Dm3[x[i]] << LOG_SIGMA);

	// Filling V1 with pointers to chunks of V2
	fillFirstLetter(V1, V2 + mm2*SIGMA, int_size_sigma, int_size, x[0], int_size_sigma*mm2);
	for (int i = 0; i < mm3; i++)
	for (int j = 0; j < SIGMA; j++)
		V1[(i << LOG_SIGMA) + j] = V2 + ((BRm3[x[mm4 - i]][j]) << LOG_SIGMA);

	// Filling V2 with pointers to chunks of V3
	fillFirstLetter(V2, V3 + mm1_sigma, int_size_sigma, int_size, x[0], int_size_sigma*mm1);
	V2[(mm3 << LOG_SIGMA) + x[1]] = V3 + ((BRm2[x[0]][x[1]] >= mm3 ? BRm2[x[0]][x[1]] : mm3) << LOG_SIGMA);
	for (int i = 0; i < mm3; i++)
	for (int j = 0; j < SIGMA; j++)
		V2[(i << LOG_SIGMA) + j] = V3 + ((Tripple[x[mm4 - i]][x[mm3 - i]][j]) << LOG_SIGMA);

	// Filling V3 with pointers to chunks of V4
	fillFirstLetter(V3, V4 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	V3[(mm2 << LOG_SIGMA) + x[1]] = V4 + ((BR_[x[0]][x[1]] >= mm2 ? BR_[x[0]][x[1]] : mm2) << LOG_SIGMA);
	for (int j = 0; j < SIGMA; j++)
		V3[(mm3 << LOG_SIGMA) + j] = V4 + ((Tripple[x[0]][x[1]][j] + 1) << LOG_SIGMA);
	for (int i = 0; i < mm3; i++)
	for (int j = 0; j < SIGMA; j++)
		V3[(i << LOG_SIGMA) + j] = V4 + ((Quad[x[mm4 - i]][x[mm3 - i]][x[mm2 - i]][j]) << LOG_SIGMA);

	// Filling V4 with pointers to chunks of V5
	fillBeginning(V4, V5, mm3, int_size_sigma, int_size);
	copy_value(V4 + (mm3 << LOG_SIGMA), V5 + SIGMA * m2m3, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V4[(mm3 << LOG_SIGMA) + x[i]] = V5 + ((Dm3[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma * 4, (unsigned char*)(V4 + mm3*SIGMA));
	for (int i = 0; i < 2; i++)
	for (int j = 0; j < SIGMA; j++)
		V4[((i + mm3) << LOG_SIGMA) + j] = V5 + (BRi[x[2 - i]][j] << LOG_SIGMA);
	V4[(mm1 << LOG_SIGMA) + x[mm3]] = V5 + (mm1 << LOG_SIGMA);

	// Filling V5 with pointers to chunks of V6
	fillBeginning(V5, V6, mm2, int_size_sigma, int_size);
	fillFirstLetter(V5 + mm2*SIGMA, V6 + m2m2*SIGMA, int_size_sigma, int_size, x[0], int_size_sigma*m);
	for (int i = 0; i < mm3; i++)
	for (int j = 0; j < SIGMA; j++)
		V5[((m + i) << LOG_SIGMA) + j] = V6 + ((BRm3[x[mm4 - i]][j] + m) << LOG_SIGMA);
	for (int i = 0; i < 2; i++)
	for (int j = 0; j < SIGMA; j++)
		V5[((mm2 + i) << LOG_SIGMA) + j] = V6 + ((BR_[x[mm2 - i]][j] + mm2) << LOG_SIGMA);

	// Filling V6 with pointers to chunks of V7
	fillBeginning(V6, V7, mm1, int_size_sigma, int_size);
	fillFirstLetter(V6 + mm1*SIGMA, V7 + m2m1*SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
	V6[(m2m3 << LOG_SIGMA) + x[1]] = V7 + (((BRm2[x[0]][x[1]] >= mm3 ? BRm2[x[0]][x[1]] : mm3) + m) << LOG_SIGMA);
	for (int i = 0; i < mm3; i++)
	for (int j = 0; j < SIGMA; j++)
		V6[((i + m) << LOG_SIGMA) + j] = V7 + ((Tripple[x[mm4 - i]][x[mm3 - i]][j] + m) << LOG_SIGMA);
	for (int i = 0; i < mm3; i++)
	for (int j = 0; j < SIGMA; j++)
		V6[((mm1 + i) << LOG_SIGMA) + j] = V7 + ((BR_[x[mm2 - i]][j] + mm1) << LOG_SIGMA);

	//Filling V7 with shift values
	fillBeginningFinal(V7, m, int_size_sigma, int_size);
	fillFirstLetterFinal(V7 + m_sigma, m2, int_size_sigma, int_size, x[0], int_size_sigma_m);
	V7[(m2m2 << LOG_SIGMA) + x[1]] = (BR_[x[0]][x[1]] >= mm2 ? BR_[x[0]][x[1]] : mm2) + m;
	for (int j = 0; j < SIGMA; j++)
		V7[(m2m3 << LOG_SIGMA) + j] = Tripple[x[0]][x[1]][j] + mp1;
	for (int i = m; i < m2m3; i++)
	for (int j = 0; j < SIGMA; j++)
		V7[(i << LOG_SIGMA) + j] = Quad[x[m2m4 - i]][x[m2m3 - i]][x[m2m2 - i]][j] + m;
	for (int i = 0; i < mm4; i++)
	for (int j = 0; j < SIGMA; j++)
		V7[((m + i) << LOG_SIGMA) + j] = BR_[x[mm2 - i]][j] + m;

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep24p += u;

	QueryPerformanceCounter(&start);

	//Search
	int *******p1, ******p2, *****p3, ****p4, ***p5, **p6, *p7;
	pos = mm4;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		p1 = V0[y[pos]];
		p2 = p1[y[pos + 1]];
		p3 = p2[y[pos + 2]];
		p4 = p3[y[pos + 3]];
		p5 = p4[y[pos + m]];
		p6 = p5[y[pos + mp1]];
		p7 = p6[y[pos + mp2]];
		r = p7[y[pos + mp3]];
		if (!r) {
			for (k = 0; k < mm4 && y[pos - mm4 + k] == x[k]; k++);
			if (k == mm4) {
				if (pos >= n)
					break;
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw24p += u;

	return count;
}

// The MAW32 algorithm with pointers
int MAW32P(unsigned char *x, const int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&prep_start);

	int *****V0[SIGMA], ****V1[SIGMA * P_MAX], ***V2[SIGMA * (P_MAX + 1)], **V3[SIGMA * P_MAX * 2],
		*V4[SIGMA * (P_MAX * 2 + 1)], V5[SIGMA * P_MAX * 3]; //V5 - shift array; V0, V1, V2, V3, V4 - pointers arrays
	int D[P_MAX], D_[P_MAX], BR_[SIGMA][SIGMA];
	int pos, r, k, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mp2 = m + 2, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2,
		m3 = m * 3, m3m1 = m * 3 - 1,
		m_sigma = m * SIGMA, mm1_sigma = mm1 * SIGMA, m2_sigma = m * 2 * SIGMA,
		int_size_sigma = int_size * SIGMA, int_size_sigma_2 = int_size_sigma * 2, int_size_sigma_m = int_size_sigma * m, int_size_sigma_m2 = int_size_sigma * m * 2;

	//Preprocessing
	buildBMHShiftTable(D, D_, x, m);
	buildBRShiftTable(BR_, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	copy_value(V0, V1 + (mm1 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V0[x[i]] = V1 + (D_[x[i]] << LOG_SIGMA);

	// Filling V1 with pointers to chunks of V2
	fillFirstLetter(V1, V2 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
	for (int j = 0; j < SIGMA; j++)
		V1[(i << LOG_SIGMA) + j] = V2 + ((BR_[x[mm2 - i]][j]) << LOG_SIGMA);

	// Filling V2 with pointers to chunks of V3
	fillBeginning(V2, V3, mm1, int_size_sigma, int_size);
	copy_value(V2 + mm1_sigma, V3 + SIGMA * m2m1, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V2[(mm1 << LOG_SIGMA) + x[i]] = V3 + ((D_[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V2 + mm1_sigma));
	V2[(mm1 << LOG_SIGMA) + x[mm1]] = V3 + (mm1 << LOG_SIGMA);

	// Filling V3 with pointers to chunks of V4
	fillBeginning(V3, V4, m, int_size_sigma, int_size);
	fillFirstLetter(V3 + m_sigma, V4 + m2_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
	for (int j = 0; j < SIGMA; j++)
		V3[((m + i) << LOG_SIGMA) + j] = V4 + ((m + BR_[x[mm2 - i]][j]) << LOG_SIGMA);

	// Filling V4 with pointers to chunks of V5
	fillBeginning(V4, V5, m2m1, int_size_sigma, int_size);
	copy_value(V4 + (m2m1 << LOG_SIGMA), V5 + (m3m1 << LOG_SIGMA), int_size_sigma * mp1, int_size);
	for (int i = 0; i < mm1; i++)
		V4[(m2m1 << LOG_SIGMA) + x[i]] = V5 + ((D_[x[i]] + m2) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V4 + (m2m1 << LOG_SIGMA)));
	V4[(m2m1 << LOG_SIGMA) + x[mm1]] = V5 + (m2m1 << LOG_SIGMA);

	// Filling V5 with shift values
	fillBeginningFinal(V5, m2, int_size_sigma, int_size);
	fillFirstLetterFinal(V5 + m2_sigma, m3, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = m; i < m2m1; i++)
	for (int j = 0; j < SIGMA; j++)
		V5[((m + i) << LOG_SIGMA) + j] = BR_[x[m2m2 - i]][j] + m2;

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep32p += u;

	QueryPerformanceCounter(&start);

	//Search
	int *****p1, ****p2, ***p3, **p4, *p5;
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		p1 = V0[y[pos]];
		p2 = p1[y[pos + 1]];
		p3 = p2[y[pos + m]];
		p4 = p3[y[pos + mp1]];
		p5 = p4[y[pos + m2]];
		r = p5[y[pos + m2p1]];
		if (!r) {
			for (k = 0; k < mm2 && y[pos - mm2 + k] == x[k]; k++);
			if (k == mm2) {
				if (pos >= n)
					break;
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw32p += u;

	return count;
}

// The MAW33 algorithm with pointers
int MAW33P(unsigned char *x, const int m, unsigned char *y, int n) {
	if (m < 3) return -1;

	QueryPerformanceCounter(&prep_start);

	int ********V0[SIGMA], *******V1[SIGMA * (P_MAX - 1)], ******V2[SIGMA * P_MAX], *****V3[SIGMA * (P_MAX + 1)], ****V4[SIGMA * (2 * P_MAX - 1)], ***V5[SIGMA * P_MAX * 2],
		**V6[SIGMA * (P_MAX * 2 + 1)], *V7[SIGMA * (3 * P_MAX - 1)], V8[SIGMA * P_MAX * 3]; //V7 - shift array; V0, V1, V2, V3, V4, V5, V6 - pointers arrays
	int D[P_MAX], D_[P_MAX], Dm2[P_MAX], BR_[SIGMA][SIGMA], BRm2[SIGMA][SIGMA], Tripple[SIGMA][SIGMA][SIGMA];
	int pos, r, k, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mp2 = m + 2, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2p2 = 2 * m + 2, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3,
		m3 = m * 3, m3m1 = m * 3 - 1, m3m2 = m * 3 - 2, m3m3 = m * 3 - 3,
		m_sigma = m * SIGMA, mm1_sigma = mm1 * SIGMA, m2_sigma = m * 2 * SIGMA,
		int_size_sigma = int_size * SIGMA, int_size_sigma_2 = int_size_sigma * 2, int_size_sigma_m = int_size_sigma * m, int_size_sigma_m2 = int_size_sigma * m * 2;

	for (int i = 0; i < SIGMA; i++)
		Dm2[i] = mm2;
	for (int i = 0; i < mm2; i++)
		Dm2[x[i]] = mm3 - i;

	BRm2[0][0] = mm1;
	mem_fill(int_size, int_size * SIGMA2, (unsigned char*)BRm2);
	for (int i = 0; i < SIGMA; i++)
		BRm2[i][x[0]] = mm2;
	for (int i = 0; i < mm2; i++)
		BRm2[x[i]][x[i + 1]] = mm3 - i;

	BR_[0][0] = m;
	mem_fill(int_size, int_size * SIGMA2, (unsigned char*)BR_);
	for (int i = 0; i < SIGMA; i++)
		BR_[i][x[0]] = mm1;
	for (int i = 0; i < mm2; i++)
		BR_[x[i]][x[i + 1]] = mm2 - i;

	//Preprocessing
	buildBMHShiftTable(D, D_, x, m);
	//buildBRShiftTable(BR_, x, m, int_size);
	buildTrippleShiftTable(Tripple, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	copy_value(V0, V1 + (mm2 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V0[x[i]] = V1 + (Dm2[x[i]] << LOG_SIGMA);

	// Filling V1 with pointers to chunks of V2
	fillFirstLetter(V1, V2 + mm1_sigma, int_size_sigma, int_size, x[0], int_size_sigma*mm1);
	for (int i = 0; i < mm2; i++)
	for (int j = 0; j < SIGMA; j++)
		V1[(i << LOG_SIGMA) + j] = V2 + ((BRm2[x[mm3 - i]][j]) << LOG_SIGMA);

	// Filling V2 with pointers to chunks of V3
	fillFirstLetter(V2, V3 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	V2[(mm2 << LOG_SIGMA) + x[1]] = V3 + ((BR_[x[0]][x[1]]) << LOG_SIGMA);
	for (int i = 0; i < mm2; i++)
	for (int j = 0; j < SIGMA; j++)
		V2[(i << LOG_SIGMA) + j] = V3 + ((Tripple[x[mm3 - i]][x[mm2 - i]][j]) << LOG_SIGMA);

	// Filling V3 with pointers to chunks of V4
	fillBeginning(V3, V4, mm2, int_size_sigma, int_size);
	copy_value(V3 + (mm2 << LOG_SIGMA), V4 + SIGMA * m2m2, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V3[(mm2 << LOG_SIGMA) + x[i]] = V4 + ((Dm2[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma * 3, (unsigned char*)(V3 + mm2*SIGMA));
	V3[(mm1 << LOG_SIGMA) + x[mm2]] = V4 + (mm1 << LOG_SIGMA);
	for (int i = 0; i < m; i++)
		V3[(mm2 << LOG_SIGMA) + x[i]] = V4 + ((m2m3 - i) << LOG_SIGMA);

	// Filling V4 with pointers to chunks of V5
	fillBeginning(V4, V5, mm1, int_size_sigma, int_size);
	fillFirstLetter(V4 + mm1*SIGMA, V5 + m2m1*SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm2; i++)
	for (int j = 0; j < SIGMA; j++)
		V4[((i + m) << LOG_SIGMA) + j] = V5 + ((BRm2[x[mm3 - i]][j] + m) << LOG_SIGMA);
	for (int j = 0; j < SIGMA; j++)
		V4[(mm1 << LOG_SIGMA) + j] = V5 + ((BRm2[x[mm2]][j] + m) << LOG_SIGMA);
	V4[(mm1 << LOG_SIGMA) + x[mm1]] = V5 + (mm1 << LOG_SIGMA);

	//Filling V5 with shift values
	fillBeginning(V5, V6, m, int_size_sigma, int_size);
	fillFirstLetter(V5 + m_sigma, V6 + m2*SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = m; i < m2m2; i++)
	for (int j = 0; j < SIGMA; j++)
		V5[(i << LOG_SIGMA) + j] = V6 + ((Tripple[x[m2m3 - i]][x[m2m2 - i]][j] + m) << LOG_SIGMA);
	V5[(m2m2 << LOG_SIGMA) + x[1]] = V6 + ((BR_[x[0]][x[1]] + m) << LOG_SIGMA);

	// Filling V6 with pointers to chunks of V7
	fillBeginning(V6, V7, m2m2, int_size_sigma, int_size);
	copy_value(V6 + (m2m2 << LOG_SIGMA), V7 + SIGMA * m3m2, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V6[((mm2 + m) << LOG_SIGMA) + x[i]] = V7 + ((Dm2[x[i]] + m2) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma * 3, (unsigned char*)(V6 + m2m2*SIGMA));
	V6[(m2m1 << LOG_SIGMA) + x[mm2]] = V7 + (m2m1 << LOG_SIGMA);
	for (int i = 0; i < m; i++)
		V6[(m2m2 << LOG_SIGMA) + x[i]] = V7 + ((m3m3 - i) << LOG_SIGMA);

	// Filling V7 with pointers to chunks of V8
	fillBeginning(V7, V8, m2m1, int_size_sigma, int_size);
	fillFirstLetter(V7 + m2m1*SIGMA, V8 + m3m1*SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm2; i++)
	for (int j = 0; j < SIGMA; j++)
		V7[((i + m2) << LOG_SIGMA) + j] = V8 + ((BRm2[x[mm3 - i]][j] + m2) << LOG_SIGMA);
	for (int j = 0; j < SIGMA; j++)
		V7[(m2m1 << LOG_SIGMA) + j] = V8 + ((BRm2[x[mm2]][j] + m2) << LOG_SIGMA);
	V7[(m2m1 << LOG_SIGMA) + x[mm1]] = V8 + (m2m1 << LOG_SIGMA);

	//Filling V8 with shift values
	fillBeginningFinal(V8, m2, int_size_sigma, int_size);
	fillFirstLetterFinal(V8 + m2_sigma, m3, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = m2; i < m3m2; i++)
	for (int j = 0; j < SIGMA; j++)
		V8[(i << LOG_SIGMA) + j] = Tripple[x[m3m3 - i]][x[m3m2 - i]][j] + m2;
	V8[(m3m2 << LOG_SIGMA) + x[1]] = BR_[x[0]][x[1]] + m2;

	QueryPerformanceCounter(&prep_end);
	u = (prep_end.QuadPart - prep_start.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep33p += u;

	QueryPerformanceCounter(&start);

	//Search
	int ********p1, *******p2, ******p3, *****p4, ****p5, ***p6, **p7, *p8;
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		p1 = V0[y[pos]];
		p2 = p1[y[pos + 1]];
		p3 = p2[y[pos + 2]];
		p4 = p3[y[pos + m]];
		p5 = p4[y[pos + mp1]];
		p6 = p5[y[pos + mp2]];
		p7 = p6[y[pos + m2]];
		p8 = p7[y[pos + m2p1]];
		r = p8[y[pos + m2p2]];
		if (!r) {
			for (k = 0; k < mm3 && y[pos - mm3 + k] == x[k]; k++);
			if (k == mm3) {
				if (pos >= n)
					break;
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw33p += u;

	return count;
}


void generateRandom() {

	srand((unsigned)time(NULL));
	for (int i = 0; i < N; i++) {
		T[i] = (rand() + glob % 320) % SIGMA;
		glob = (glob * 11 + 30157) % 499;
	}
	for (int i = 0; i < m; i++) {
		P[i] = (rand() + glob) % SIGMA;
		glob = (glob * 123 + 3157) % 893;
	}
}

void DNA() {
	N = 4638680;
	FILE *f;
	f = fopen("ecoli.txt", "rt");
	fread(T1, 1, N, f);
	for (int i = 0; i < N; i++)
		switch (T1[i]) {
		case 'a': T[i] = 0; break;
		case 'c': T[i] = 1; break;
		case 't': T[i] = 2; break;
		case 'g': T[i] = 3;
	}
}

void testPerformance() {
	QueryPerformanceFrequency(&freq);
	QueryPerformanceFrequency(&_freq);

	f = fopen("output.csv", "wt");
	generateRandom();
	//DNA();

	fprintf(f, "b=%d N=%d ITER=%d\n", SIGMA, N, ITER);

	fprintf(f, "m,MAW22,MAW23,MAW24,MAW32,MAW33,MAW22P,MAW22p_mod,MAW23P,MAW24P,MAW32P,MAW33P,,PREP22,PREP23,PREP24,PREP32,PREP33,PREP22P,PREP22P_mod,PREP23P,PREP24P,PREP32P,PREP33P,,SUM22,SUM23,SUM24,SUM32,SUM33,SUM22P,SUM22P_mod,SUM23P,SUM24P,SUM32P,SUM33P");
	for (m = 2; m < 81; m < 10 ? m++ : m += 10) {
		sum_maw22 = sum_maw23 = sum_maw24 = sum_maw32 = sum_maw33 = 0;
		sum_prep22 = sum_prep23 = sum_prep24 = sum_prep32 = sum_prep33 = 0;

		sum_maw22p = sum_maw23p = sum_maw24p = sum_maw32p = sum_maw33p = sum_maw22p_mod = 0;
		sum_prep22p = sum_prep23p = sum_prep24p = sum_prep32p = sum_prep33p = sum_prep22p_mod = 0;

		int nm = N - m;
		memcpy(T1, T, N);
		for (int ii = 0; ii < ITER; ii++) {
			srand((unsigned)time(NULL));
			int patpos = rand() % (N - m - 2);
			for (int i = 0; i < m; i++)
				P[i] = T[patpos + i];
			for (int i = 0; i < m; i++)
				T[N - m + i] = P[i];
			
			maw22 = MAW22(P, m, T, nm);
			maw23 = MAW23(P, m, T, nm);
			maw24 = MAW24(P, m, T, nm);
			maw32 = MAW32(P, m, T, nm);
			maw33 = MAW33(P, m, T, nm);

			maw22p = MAW22P(P, m, T, nm);
			maw22p_mod = MAW22P_mod(P, m, T, nm);
			maw23p = MAW23P(P, m, T, nm);
			maw24p = MAW24P(P, m, T, nm);
			maw32p = MAW32P(P, m, T, nm);
			maw33p = MAW33P(P, m, T, nm);
		}
		printf("b=%d m=%d\n", SIGMA, m);
		printf("%d %d %d %d %d %d %d %d %d %d %d\n\n", maw22, maw23, maw24, maw32, maw33, maw22p, maw22p_mod, maw23p, maw24p, maw32p, maw33p);
		printf("%7.lld %7.lld %7.lld %7.lld %7.lld %7.lld %7.lld %7.lld %7.lld %7.lld %7.lld\n\n",
			sum_prep22, sum_prep23, sum_prep24, sum_prep32, sum_prep33, sum_prep22p, sum_prep22p_mod, sum_prep23p, sum_prep24p, sum_prep32p, sum_prep33p);
		fprintf(f, "\n%2.d,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld",
			m, sum_maw22, sum_maw23, sum_maw24, sum_maw32, sum_maw33, sum_maw22p, sum_maw22p_mod, sum_maw23p, sum_maw24p, sum_maw32p, sum_maw33p,
			sum_prep22, sum_prep23, sum_prep24, sum_prep32, sum_prep33, sum_prep22p, sum_prep22p_mod, sum_prep23p, sum_prep24p, sum_prep32p, sum_prep33p,
			sum_maw22 + sum_prep22, sum_maw23 + sum_prep23, sum_maw24 + sum_prep24, sum_maw32 + sum_prep32, sum_maw33 + sum_prep33,
			sum_maw22p + sum_prep22p, sum_maw22p_mod + sum_prep22p_mod, sum_maw23p + sum_prep23p, sum_maw24p + sum_prep24p, sum_maw32p + sum_prep32p, sum_maw33p + sum_prep33p);

	}
	fclose(f);
	system("pause");
}

int testMAW()
{
	QueryPerformanceFrequency(&freq);
	QueryPerformanceFrequency(&_freq);

	N = 10000000;

	for (int num = 1000; num < 100000000; num++)
	{
		generateRandom();

		srand((unsigned)time(NULL));
		int patpos = rand() % (N - m - 2);
		for (int i = 0; i < m; i++)
			P[i] = T[patpos + i];

		m = 0;
		int num2 = num;
		bool f = false;
		while (num2) {
			P[m] = (num2 % 10);
			if (P[m] >= SIGMA) { f = true; break; }
			m++;
			num2 /= 10;
		}
		if (f) {
			if (m > 5) cout << num << endl;
			num += 6 * pow(10, m) - 1;
			continue;
		}

		for (int i = 0; i < m; i++)
			T[N - m + i] = P[i];

		maw22 = MAW22(P, m, T, N - m);
		maw22p = MAW22P(P, m, T, N - m);
		maw22p_mod = MAW22P_mod(P, m, T, N - m);
		maw23p = MAW23P(P, m, T, N - m);
		maw24p = MAW24P(P, m, T, N - m);
		maw32p = MAW32P(P, m, T, N - m);
		maw33p = MAW33P(P, m, T, N - m);

		if (maw22 != maw22p || maw22 != maw22p_mod || maw22 != maw23p || maw22 != maw24p || maw22 != maw32p || maw22 != maw33p)
		{
			for (int i = 0; i < m; i++)
				cout << char(P[i] + '0');
			cout << " " << m << endl;

			cout << "amount " << maw22 << " " << maw22p << " " << maw22p_mod << " " << maw23p << " " << maw24p << " " << maw32p << " " << maw33p << endl;
		}
	}

	system("pause");
	return 0;
}

int simplestTest()
{
	QueryPerformanceFrequency(&freq);
	QueryPerformanceFrequency(&_freq);

	N = 10000000;
	int m = 7;

	for (m = 4; m < 9; m++)
	for (int k = 0; k < 10; k++)
	{
		generateRandom();

		srand((unsigned)time(NULL));
		int patpos = rand() % (N - m - 2);
		for (int i = 0; i < m; i++)
			P[i] = T[patpos + i];

		for (int i = 0; i < m; i++)
			T[N - m + i] = P[i];

		//P[0] = 3; P[1] = 1; P[2] = 0; P[3] = 3; 
		//P[0] = 3; P[1] = 3; P[2] = 2; P[3] = 1;
		//P[0] = 3; P[1] = 3; P[2] = 2; P[3] = 1; P[4] = 0; P[5] = 1;
		//P[0] = 3; P[1] = 1; P[2] = 2; P[3] = 2; P[4] = 0; P[5] = 0;
		//P[0] = 3; P[1] = 3; P[2] = 3; P[3] = 1; P[4] = 0; P[5] = 1;
		//P[0] = 1; P[1] = 3; P[2] = 0; P[3] = 3; P[4] = 2; P[5] = 3;
		//P[0] = 0; P[1] = 1; P[2] = 0; P[3] = 3; P[4] = 2;
		//P[0] = 2; P[1] = 0; P[2] = 2; P[3] = 0; P[4] = 0;
		//P[0] = 2; P[1] = 1; P[2] = 2; P[3] = 1; P[4] = 1;
		//P[0] = 3; P[1] = 3; P[2] = 1; P[3] = 0; P[4] = 3; P[5] = 3; P[6] = 2;
		//P[0] = 1; P[1] = 0; P[2] = 3; P[3] = 2; P[4] = 1; P[5] = 0; P[6] = 1;
		//P[0] = 2; P[1] = 1; P[2] = 3; P[3] = 0; P[4] = 2; P[5] = 1; P[6] = 1;
		//P[0] = 2; P[1] = 1; P[2] = 1; P[3] = 2; P[4] = 2; P[5] = 1; P[6] = 1;

		maw22 = MAW22(P, m, T, N - m);
		maw22p = MAW22P(P, m, T, N - m);
		maw22p_mod = MAW22P_mod(P, m, T, N - m);
		maw23p = MAW23P(P, m, T, N - m);
		maw24p = MAW24P(P, m, T, N - m);
		maw32p = MAW32P(P, m, T, N - m);
		maw33p = MAW33P(P, m, T, N - m);

		for (int i = 0; i < m; i++)
			cout << char(P[i] + '0');
		cout << " " << m << endl;

		cout << "amount " << maw22 << " " << maw22p << " " << maw22p_mod << " " << maw23p << " " << maw24p << " " << maw32p << " " << maw33p << endl;

		/*cout << "amount " << maw22 << "	" << maw22p << "	" << maw22p_mod << "	" << maw23p << "	" << maw24p << "	" << maw32p << "	" << maw33p << endl << endl;
		cout << "algo   " << sum_maw22 << "	" << sum_maw22p << "	" << sum_maw22p_mod << "	" << sum_maw23p << "	" << sum_maw24p << "	" << sum_maw32p << "	" << sum_maw33p << endl;
		cout << "prep   " << sum_prep22 << "	" << sum_prep22p << "	" << sum_prep22p_mod << "	" << sum_prep23p << "	" << sum_prep24p << "	" << sum_prep32p << "	" << sum_prep33p << endl;
		cout << "sum    " << sum_maw22 + sum_prep22 << "	" << sum_maw22p + sum_prep22p << "	" << sum_maw22p_mod + sum_prep22p_mod << "	" << sum_maw23p + sum_prep23p << "	" <<
			sum_maw24p + sum_prep24p << "	" << sum_maw32p + sum_prep32p << "	" << sum_maw33p + sum_prep33p << endl;

		cout << endl << maw22c << " " << maw22pc << " " << maw22p_modc << endl;
		*/
	}

	system("pause");
	return 0;
}

int main()
{
	testPerformance();
	//testMAW();
	//simplestTest();
	return 0;
}
