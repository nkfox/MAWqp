﻿#define _CRT_SECURE_NO_WARNINGS
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
/*
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

#define SIGMA 6		//alphabet size
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

/*#define SIGMA 32		//alphabet size
#define SIGMA2 1024
#define SIGMA3 32768
#define SIGMA4 1048576
#define LOG_SIGMA 5
*/
#define SIGMA 64		//alphabet size
#define SIGMA2 4096
#define SIGMA3 262144
#define SIGMA4 16777216
#define LOG_SIGMA 6

#define SIGMA5 7776
#define SIGMA6 46656
#define SIGMA7 279936
#define SIGMA8 1679616
#define SIGMA9 10077696


const int TOTAL = 10000000 + 5 * P_MAX;
unsigned char T[TOTAL], T1[TOTAL], P[P_MAX], P1[P_MAX];
int N = TOTAL - 5 * P_MAX, ITER = 200, m = 5;

FILE * f;
LARGE_INTEGER start, _end, freq, _freq, prep_start, prep_end;
double u;
int nm2, glob = 0;

long long sum_maw22, sum_maw23, sum_maw24, sum_maw32, sum_maw33;
int maw22, maw23, maw24, maw32, maw33;

long long sum_maw22p, sum_maw23p, sum_maw24p, sum_maw32p, sum_maw33p, sum_maw42p;
int maw22p, maw23p, maw24p, maw32p, maw33p, maw42p;

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
	QueryPerformanceCounter(&start);

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

	//Search
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M22 + y[pos] * SIGMA3 + y[pos + 1] * SIGMA2 + y[pos + m] * SIGMA + y[pos + mp1]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < m && y[pos + k - mm2] == x[k]; k++);
			if (k == m) {
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
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
	QueryPerformanceCounter(&start);

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

	//Search
	pos = mm3;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M23 + y[pos] * SIGMA5 + y[pos + 1] * SIGMA4 + y[pos + 2] * SIGMA3 + y[pos + m] * SIGMA2 + y[pos + mp1] * SIGMA + y[pos + mp2]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < m && y[pos + k - mm3] == x[k]; k++);
			if (k == m) {
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
	QueryPerformanceCounter(&start);

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

	//Search
	pos = mm4;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M24 + y[pos] * SIGMA7 + y[pos + 1] * SIGMA6 + y[pos + 2] * SIGMA5 + y[pos + 3] * SIGMA4 +
			y[pos + m] * SIGMA3 + y[pos + mp1] * SIGMA2 + y[pos + mp2] * SIGMA + y[pos + mp3]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < m && y[pos + k - mm4] == x[k]; k++);
			if (k == m) {
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
	QueryPerformanceCounter(&start);

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

	//Search
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M32 + y[pos] * SIGMA5 + y[pos + 1] * SIGMA4 + y[pos + m] * SIGMA3 + y[pos + mp1] * SIGMA2 + y[pos + m2] * SIGMA + y[pos + m2p1]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < m && y[pos + k - mm2] == x[k]; k++);
			if (k == m) {
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
	QueryPerformanceCounter(&start);

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
			for (k = 0; k < m && y[pos + k - mm3] == x[k]; k++);
			if (k == m) {
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

//-----------------------------------MAWP--------------------------------------------------

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
int MAW22P(unsigned char *x, const int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start);

	int ***V0[SIGMA], **V1[SIGMA * P_MAX], *V2[SIGMA * (P_MAX + 1)], V3[SIGMA * P_MAX * 2]; //V3 - shift array; V0, V1, V2 - pointers arrays
	int V0m[SIGMA], V1m[SIGMA * P_MAX], V2m[SIGMA * (P_MAX + 1)], V3m[SIGMA * P_MAX * 2];
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

	copy_value(V0m, (mm1 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V0m[x[i]] = (D_[x[i]] << LOG_SIGMA);

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

	//Search
	int ***p1, **p2, *p3;
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern          /// ??? memcpy or delete this line
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
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw22p += u;

	return count;
}

class Map {
public:
	Map(const int m, unsigned char* x) {
		this->x = x;
		this->m = m;
		this->mm1 = m - 1;
		this->mm2 = m - 2;
		this->mm3 = m - 3;
	}

	void init3() {
		std::string s;
		for (int i = 0; i < mm2; i++) {
			s = "";
			s.push_back(x[i]);
			s.push_back(x[i + 1]);
			s.push_back(x[i + 2]);
			mapa[s] = mm3 - i;
		}
	}

	void init3shift() {
		this->mm4 = m - 4;
		std::string s;
		for (int i = 0; i < mm3; i++) {
			s = "";
			s.push_back(x[i]);
			s.push_back(x[i + 1]);
			s.push_back(x[i + 2]);
			mapa[s] = mm4 - i;
		}
	}

	void init4() {
		this->mm4 = m - 4;
		std::string s;
		for (int i = 0; i < mm3; i++) {
			s = "";
			s.push_back(x[i]);
			s.push_back(x[i + 1]);
			s.push_back(x[i + 2]);
			s.push_back(x[i + 3]);
			mapa[s] = mm4 - i;
		}
	}

	int get3(char c1, char c2, char c3) {
		std::string s = "";
		s.push_back(c1);
		s.push_back(c2);
		s.push_back(c3);

		int res = m;
		if (c3 == x[0])
			res = mm1;
		if (c2 == x[0] && c3 == x[1])
			res = mm2;
		iter = mapa.find(s);
		if (iter != mapa.end())
			res = iter->second;
		return res;
	}

	int get3shift(char c1, char c2, char c3) {
		std::string s = "";
		s.push_back(c1);
		s.push_back(c2);
		s.push_back(c3);

		int res = mm1;
		if (c3 == x[0])
			res = mm2;
		if (c2 == x[0] && c3 == x[1])
			res = mm3;
		iter = mapa.find(s);
		if (iter != mapa.end())
			res = iter->second;
		return res;
	}

	int get4(char c1, char c2, char c3, char c4) {
		std::string s = "";
		s.push_back(c1);
		s.push_back(c2);
		s.push_back(c3);
		s.push_back(c4);

		int res = m;
		if (c4 == x[0])
			res = mm1;
		if (c3 == x[0] && c4 == x[1])
			res = mm2;
		if (c2 == x[0] && c3 == x[1] && c4 == x[2])
			res = mm3;
		iter = mapa.find(s);
		if (iter != mapa.end())
			res = iter->second;
		return res;
	}

private:
	std::map<std::string, int> mapa;
	std::map<std::string, int>::iterator iter;
	int m, mm1, mm2, mm3, mm4;
	unsigned char* x;
};

// The MAW23 algorithm with pointers
int MAW23P(unsigned char *x, const int m, unsigned char *y, int n) {
	if (m < 3) return -1;

	QueryPerformanceCounter(&start);

	int *****V0[SIGMA], ****V1[SIGMA * (P_MAX - 1)], ***V2[SIGMA * P_MAX], **V3[SIGMA * (P_MAX + 1)],
		*V4[SIGMA * (2 * P_MAX - 1)], V5[SIGMA * P_MAX * 2]; //V5 - shift array; V0, V1, V2, V3, V4 - pointers arrays
	int D[P_MAX], D_[P_MAX], Dm2[P_MAX], BR_[SIGMA][SIGMA], BRm2[SIGMA][SIGMA];
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

	Map Tripple(m, x);
	Tripple.init3();

	//Preprocessing
	buildBMHShiftTable(D, D_, x, m);
	//buildBRShiftTable(BR_, x, m, int_size);

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
			V2[(i << LOG_SIGMA) + j] = V3 + (Tripple.get3(x[mm3 - i], x[mm2 - i], j) << LOG_SIGMA);

	// Filling V3 with pointers to chunks of V4
	fillBeginning(V3, V4, mm2, int_size_sigma, int_size);
	copy_value(V3 + (mm2 << LOG_SIGMA), V4 + SIGMA * m2m2, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V3[(mm2 << LOG_SIGMA) + x[i]] = V4 + ((Dm2[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma * 3, (unsigned char*)(V3 + mm2 * SIGMA));
	V3[(mm1 << LOG_SIGMA) + x[mm2]] = V4 + (mm1 << LOG_SIGMA);
	for (int i = 0; i < m; i++)
		V3[(mm2 << LOG_SIGMA) + x[i]] = V4 + ((m2m3 - i) << LOG_SIGMA);

	// Filling V4 with pointers to chunks of V5
	fillBeginning(V4, V5, mm1, int_size_sigma, int_size);
	fillFirstLetter(V4 + mm1 * SIGMA, V5 + m2m1 * SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
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
			V5[(i << LOG_SIGMA) + j] = Tripple.get3(x[m2m3 - i], x[m2m2 - i], j) + m;
	V5[(m2m2 << LOG_SIGMA) + x[1]] = BR_[x[0]][x[1]] + m;

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
			pos += D[y[pos + 2]];
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

	QueryPerformanceCounter(&start);

	int *******V0[SIGMA], ******V1[SIGMA * (P_MAX - 2)], *****V2[SIGMA * (P_MAX - 1)], ****V3[SIGMA * P_MAX], ***V4[SIGMA * (P_MAX + 1)],
		**V5[SIGMA * (2 * P_MAX - 1)], *V6[SIGMA * (2 * P_MAX - 2)], V7[SIGMA * P_MAX * 2]; //V7 - shift array; V0, V1, V2, V3, V4, V5, V6 - pointers arrays
	int D[P_MAX], D_[P_MAX], Dm2[P_MAX], Dm3[P_MAX], BR_[SIGMA][SIGMA], BRm2[SIGMA][SIGMA], BRm3[SIGMA][SIGMA], BRi[SIGMA][SIGMA];
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

	Map Tripple(m, x);
	Tripple.init3shift();

	Map Quad(m, x);
	Quad.init4();

	//Preprocessing
	buildBMHShiftTable(D, D_, x, m);
	//buildBRShiftTable(BR_, x, m, int_size);
	//buildTrippleShiftTable(Tripple, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	copy_value(V0, V1 + (mm3 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm3; i++)
		V0[x[i]] = V1 + (Dm3[x[i]] << LOG_SIGMA);

	// Filling V1 with pointers to chunks of V2
	fillFirstLetter(V1, V2 + mm2 * SIGMA, int_size_sigma, int_size, x[0], int_size_sigma*mm2);
	for (int i = 0; i < mm3; i++)
		for (int j = 0; j < SIGMA; j++)
			V1[(i << LOG_SIGMA) + j] = V2 + ((BRm3[x[mm4 - i]][j]) << LOG_SIGMA);

	// Filling V2 with pointers to chunks of V3
	fillFirstLetter(V2, V3 + mm1_sigma, int_size_sigma, int_size, x[0], int_size_sigma*mm1);
	V2[(mm3 << LOG_SIGMA) + x[1]] = V3 + ((BRm2[x[0]][x[1]] >= mm3 ? BRm2[x[0]][x[1]] : mm3) << LOG_SIGMA);
	for (int i = 0; i < mm3; i++)
		for (int j = 0; j < SIGMA; j++)
			V2[(i << LOG_SIGMA) + j] = V3 + (Tripple.get3shift(x[mm4 - i], x[mm3 - i], j) << LOG_SIGMA);

	// Filling V3 with pointers to chunks of V4
	fillFirstLetter(V3, V4 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	V3[(mm2 << LOG_SIGMA) + x[1]] = V4 + ((BR_[x[0]][x[1]] >= mm2 ? BR_[x[0]][x[1]] : mm2) << LOG_SIGMA);
	for (int j = 0; j < SIGMA; j++)
		V3[(mm3 << LOG_SIGMA) + j] = V4 + ((Tripple.get3shift(x[0], x[1], j) + 1) << LOG_SIGMA);
	for (int i = 0; i < mm3; i++)
		for (int j = 0; j < SIGMA; j++)
			V3[(i << LOG_SIGMA) + j] = V4 + (Quad.get4(x[mm4 - i], x[mm3 - i], x[mm2 - i], j) << LOG_SIGMA);

	// Filling V4 with pointers to chunks of V5
	fillBeginning(V4, V5, mm3, int_size_sigma, int_size);
	copy_value(V4 + (mm3 << LOG_SIGMA), V5 + SIGMA * m2m3, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V4[(mm3 << LOG_SIGMA) + x[i]] = V5 + ((Dm3[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma * 4, (unsigned char*)(V4 + mm3 * SIGMA));
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < SIGMA; j++)
			V4[((i + mm3) << LOG_SIGMA) + j] = V5 + (BRi[x[2 - i]][j] << LOG_SIGMA);
	V4[(mm1 << LOG_SIGMA) + x[mm3]] = V5 + (mm1 << LOG_SIGMA);

	// Filling V5 with pointers to chunks of V6
	fillBeginning(V5, V6, mm2, int_size_sigma, int_size);
	fillFirstLetter(V5 + mm2 * SIGMA, V6 + m2m2 * SIGMA, int_size_sigma, int_size, x[0], int_size_sigma*m);
	for (int i = 0; i < mm3; i++)
		for (int j = 0; j < SIGMA; j++)
			V5[((m + i) << LOG_SIGMA) + j] = V6 + ((BRm3[x[mm4 - i]][j] + m) << LOG_SIGMA);
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < SIGMA; j++)
			V5[((mm2 + i) << LOG_SIGMA) + j] = V6 + ((BR_[x[mm2 - i]][j] + mm2) << LOG_SIGMA);

	// Filling V6 with pointers to chunks of V7
	fillBeginning(V6, V7, mm1, int_size_sigma, int_size);
	fillFirstLetter(V6 + mm1 * SIGMA, V7 + m2m1 * SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
	V6[(m2m3 << LOG_SIGMA) + x[1]] = V7 + (((BRm2[x[0]][x[1]] >= mm3 ? BRm2[x[0]][x[1]] : mm3) + m) << LOG_SIGMA);
	for (int i = 0; i < mm3; i++)
		for (int j = 0; j < SIGMA; j++)
			V6[((i + m) << LOG_SIGMA) + j] = V7 + ((Tripple.get3shift(x[mm4 - i], x[mm3 - i], j) + m) << LOG_SIGMA);
	for (int i = 0; i < mm3; i++)
		for (int j = 0; j < SIGMA; j++)
			V6[((mm1 + i) << LOG_SIGMA) + j] = V7 + ((BR_[x[mm2 - i]][j] + mm1) << LOG_SIGMA);

	//Filling V7 with shift values
	fillBeginningFinal(V7, m, int_size_sigma, int_size);
	fillFirstLetterFinal(V7 + m_sigma, m2, int_size_sigma, int_size, x[0], int_size_sigma_m);
	V7[(m2m2 << LOG_SIGMA) + x[1]] = (BR_[x[0]][x[1]] >= mm2 ? BR_[x[0]][x[1]] : mm2) + m;
	for (int j = 0; j < SIGMA; j++)
		V7[(m2m3 << LOG_SIGMA) + j] = Tripple.get3shift(x[0], x[1], j) + mp1;
	for (int i = m; i < m2m3; i++)
		for (int j = 0; j < SIGMA; j++)
			V7[(i << LOG_SIGMA) + j] = Quad.get4(x[m2m4 - i], x[m2m3 - i], x[m2m2 - i], j) + m;
	for (int i = 0; i < mm4; i++)
		for (int j = 0; j < SIGMA; j++)
			V7[((m + i) << LOG_SIGMA) + j] = BR_[x[mm2 - i]][j] + m;

	//Search
	int *******p1, ******p2, *****p3, ****p4, ***p5, **p6, *p7;
	pos = mm4;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		p1 = V0[y[pos]]; // y[pos] заменить на указатель ???
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
			pos += D[y[pos + 3]];
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
	QueryPerformanceCounter(&start);

	int *****V0[SIGMA], ****V1[SIGMA * P_MAX], ***V2[SIGMA * (P_MAX + 1)], **V3[SIGMA * P_MAX * 2],
		*V4[SIGMA * (P_MAX * 2 + 1)], V5[SIGMA * P_MAX * 3]; //V5 - shift array; V0, V1, V2, V3, V4 - pointers arrays
	int V0m[SIGMA], V1m[SIGMA * P_MAX], V2m[SIGMA * (P_MAX + 1)], V3m[SIGMA * P_MAX * 2], V4m[SIGMA * (P_MAX * 2 + 1)], V5m[SIGMA * P_MAX * 3];
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

	copy_value(V0m, (mm1 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V0m[x[i]] = (D_[x[i]] << LOG_SIGMA);

	// Filling V1 with pointers to chunks of V2
	fillFirstLetter(V1, V2 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
		for (int j = 0; j < SIGMA; j++)
			V1[(i << LOG_SIGMA) + j] = V2 + ((BR_[x[mm2 - i]][j]) << LOG_SIGMA);

	fillFirstLetter(V1m, m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
		for (int j = 0; j < SIGMA; j++)
			V1m[(i << LOG_SIGMA) + j] = ((BR_[x[mm2 - i]][j]) << LOG_SIGMA);

	// Filling V2 with pointers to chunks of V3
	fillBeginning(V2, V3, mm1, int_size_sigma, int_size);
	copy_value(V2 + mm1_sigma, V3 + SIGMA * m2m1, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V2[(mm1 << LOG_SIGMA) + x[i]] = V3 + ((D_[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V2 + mm1_sigma));
	V2[(mm1 << LOG_SIGMA) + x[mm1]] = V3 + (mm1 << LOG_SIGMA);

	fillBeginning(V2m, 0, mm1, int_size_sigma, int_size);
	copy_value(V2m + mm1_sigma, SIGMA * m2m1, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V2m[(mm1 << LOG_SIGMA) + x[i]] = ((D_[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V2m + mm1_sigma));
	V2m[(mm1 << LOG_SIGMA) + x[mm1]] = (mm1 << LOG_SIGMA);

	// Filling V3 with pointers to chunks of V4
	fillBeginning(V3, V4, m, int_size_sigma, int_size);
	fillFirstLetter(V3 + m_sigma, V4 + m2_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
		for (int j = 0; j < SIGMA; j++)
			V3[((m + i) << LOG_SIGMA) + j] = V4 + ((m + BR_[x[mm2 - i]][j]) << LOG_SIGMA);

	fillBeginning(V3m, 0, m, int_size_sigma, int_size);
	fillFirstLetter(V3m + m_sigma, m2_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
		for (int j = 0; j < SIGMA; j++)
			V3m[((m + i) << LOG_SIGMA) + j] = ((m + BR_[x[mm2 - i]][j]) << LOG_SIGMA);

	// Filling V4 with pointers to chunks of V5
	fillBeginning(V4, V5, m2m1, int_size_sigma, int_size);
	copy_value(V4 + (m2m1 << LOG_SIGMA), V5 + (m3m1 << LOG_SIGMA), int_size_sigma * mp1, int_size);
	for (int i = 0; i < mm1; i++)
		V4[(m2m1 << LOG_SIGMA) + x[i]] = V5 + ((D_[x[i]] + m2) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V4 + (m2m1 << LOG_SIGMA)));
	V4[(m2m1 << LOG_SIGMA) + x[mm1]] = V5 + (m2m1 << LOG_SIGMA);

	fillBeginning(V4m, 0, m2m1, int_size_sigma, int_size);
	copy_value(V4m + (m2m1 << LOG_SIGMA), (m3m1 << LOG_SIGMA), int_size_sigma * mp1, int_size);
	for (int i = 0; i < mm1; i++)
		V4m[(m2m1 << LOG_SIGMA) + x[i]] = ((D_[x[i]] + m2) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V4m + (m2m1 << LOG_SIGMA)));
	V4m[(m2m1 << LOG_SIGMA) + x[mm1]] = (m2m1 << LOG_SIGMA);

	// Filling V5 with shift values
	fillBeginningFinal(V5, m2, int_size_sigma, int_size);
	fillFirstLetterFinal(V5 + m2_sigma, m3, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = m; i < m2m1; i++)
		for (int j = 0; j < SIGMA; j++)
			V5[((m + i) << LOG_SIGMA) + j] = BR_[x[m2m2 - i]][j] + m2;

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
			pos += D[y[pos + 1]];
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

	QueryPerformanceCounter(&start);

	int ********V0[SIGMA], *******V1[SIGMA * (P_MAX - 1)], ******V2[SIGMA * P_MAX], *****V3[SIGMA * (P_MAX + 1)], ****V4[SIGMA * (2 * P_MAX - 1)], ***V5[SIGMA * P_MAX * 2],
		**V6[SIGMA * (P_MAX * 2 + 1)], *V7[SIGMA * (3 * P_MAX - 1)], V8[SIGMA * P_MAX * 3]; //V7 - shift array; V0, V1, V2, V3, V4, V5, V6 - pointers arrays
	int D[P_MAX], D_[P_MAX], Dm2[P_MAX], BR_[SIGMA][SIGMA], BRm2[SIGMA][SIGMA];
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

	Map Tripple(m, x);
	Tripple.init3();

	//Preprocessing
	buildBMHShiftTable(D, D_, x, m);
	//buildBRShiftTable(BR_, x, m, int_size);

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
			V2[(i << LOG_SIGMA) + j] = V3 + (Tripple.get3(x[mm3 - i], x[mm2 - i], j) << LOG_SIGMA);

	// Filling V3 with pointers to chunks of V4
	fillBeginning(V3, V4, mm2, int_size_sigma, int_size);
	copy_value(V3 + (mm2 << LOG_SIGMA), V4 + SIGMA * m2m2, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V3[(mm2 << LOG_SIGMA) + x[i]] = V4 + ((Dm2[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma * 3, (unsigned char*)(V3 + mm2 * SIGMA));
	V3[(mm1 << LOG_SIGMA) + x[mm2]] = V4 + (mm1 << LOG_SIGMA);
	for (int i = 0; i < m; i++)
		V3[(mm2 << LOG_SIGMA) + x[i]] = V4 + ((m2m3 - i) << LOG_SIGMA);

	// Filling V4 with pointers to chunks of V5
	fillBeginning(V4, V5, mm1, int_size_sigma, int_size);
	fillFirstLetter(V4 + mm1 * SIGMA, V5 + m2m1 * SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm2; i++)
		for (int j = 0; j < SIGMA; j++)
			V4[((i + m) << LOG_SIGMA) + j] = V5 + ((BRm2[x[mm3 - i]][j] + m) << LOG_SIGMA);
	for (int j = 0; j < SIGMA; j++)
		V4[(mm1 << LOG_SIGMA) + j] = V5 + ((BRm2[x[mm2]][j] + m) << LOG_SIGMA);
	V4[(mm1 << LOG_SIGMA) + x[mm1]] = V5 + (mm1 << LOG_SIGMA);

	//Filling V5 with shift values
	fillBeginning(V5, V6, m, int_size_sigma, int_size);
	fillFirstLetter(V5 + m_sigma, V6 + m2 * SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = m; i < m2m2; i++)
		for (int j = 0; j < SIGMA; j++)
			V5[(i << LOG_SIGMA) + j] = V6 + ((Tripple.get3(x[m2m3 - i], x[m2m2 - i], j) + m) << LOG_SIGMA);
	V5[(m2m2 << LOG_SIGMA) + x[1]] = V6 + ((BR_[x[0]][x[1]] + m) << LOG_SIGMA);

	// Filling V6 with pointers to chunks of V7
	fillBeginning(V6, V7, m2m2, int_size_sigma, int_size);
	copy_value(V6 + (m2m2 << LOG_SIGMA), V7 + SIGMA * m3m2, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		V6[((mm2 + m) << LOG_SIGMA) + x[i]] = V7 + ((Dm2[x[i]] + m2) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma * 3, (unsigned char*)(V6 + m2m2 * SIGMA));
	V6[(m2m1 << LOG_SIGMA) + x[mm2]] = V7 + (m2m1 << LOG_SIGMA);
	for (int i = 0; i < m; i++)
		V6[(m2m2 << LOG_SIGMA) + x[i]] = V7 + ((m3m3 - i) << LOG_SIGMA);

	// Filling V7 with pointers to chunks of V8
	fillBeginning(V7, V8, m2m1, int_size_sigma, int_size);
	fillFirstLetter(V7 + m2m1 * SIGMA, V8 + m3m1 * SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m);
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
			V8[(i << LOG_SIGMA) + j] = Tripple.get3(x[m3m3 - i], x[m3m2 - i], j) + m2;
	V8[(m3m2 << LOG_SIGMA) + x[1]] = BR_[x[0]][x[1]] + m2;

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
			pos += D[y[pos + 2]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw33p += u;

	return count;
}

// The MAW42 algorithm with pointers
int MAW42P(unsigned char *x, const int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start);

	int *******V0[SIGMA], ******V1[SIGMA * P_MAX], *****V2[SIGMA * (P_MAX + 1)], ****V3[SIGMA * P_MAX * 2], ***V4[SIGMA * (P_MAX * 2 + 1)],
		**V5[SIGMA * P_MAX * 3], *V6[SIGMA * (P_MAX * 3 + 1)], V7[SIGMA * P_MAX * 4]; //V7 - shift array; V0, V1, V2, V3, V4, V5, V6 - pointers arrays
																					  //int V0m[SIGMA], V1m[SIGMA * 4], V2m[SIGMA * (4 + 1)], V3m[SIGMA * 4 * 2], V4m[SIGMA * (4 * 2 + 1)],
																					  //	V5m[SIGMA * 4 * 3 + 3], V6m[SIGMA * (4 * 3 + 1)], V7m[SIGMA * 4 * 4]; 
	int D[P_MAX], D_[P_MAX], BR_[SIGMA][SIGMA];
	int pos, r, k, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mp2 = m + 2, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2,
		m3 = m * 3, m3m1 = m * 3 - 1, m3p1 = m * 3 + 1,
		m4 = m * 4, m4m1 = m * 4 - 1,
		m_sigma = m * SIGMA, mm1_sigma = mm1 * SIGMA, m2_sigma = m * 2 * SIGMA, m3_sigma = m * 3 * SIGMA,
		int_size_sigma = int_size * SIGMA, int_size_sigma_2 = int_size_sigma * 2, int_size_sigma_m = int_size_sigma * m, int_size_sigma_m2 = int_size_sigma * m * 2;

	//Preprocessing
	buildBMHShiftTable(D, D_, x, m);
	buildBRShiftTable(BR_, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	copy_value(V0, V1 + (mm1 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V0[x[i]] = V1 + (D_[x[i]] << LOG_SIGMA);

	/*copy_value(V0m, (mm1 << LOG_SIGMA), int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
	V0m[x[i]] = (D_[x[i]] << LOG_SIGMA);*/

	// Filling V1 with pointers to chunks of V2
	fillFirstLetter(V1, V2 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
		for (int j = 0; j < SIGMA; j++)
			V1[(i << LOG_SIGMA) + j] = V2 + ((BR_[x[mm2 - i]][j]) << LOG_SIGMA);

	/*fillFirstLetter(V1m, m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
	for (int j = 0; j < SIGMA; j++)
	V1m[(i << LOG_SIGMA) + j] = ((BR_[x[mm2 - i]][j]) << LOG_SIGMA);*/

	// Filling V2 with pointers to chunks of V3
	fillBeginning(V2, V3, mm1, int_size_sigma, int_size);
	copy_value(V2 + mm1_sigma, V3 + SIGMA * m2m1, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		V2[(mm1 << LOG_SIGMA) + x[i]] = V3 + ((D_[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V2 + mm1_sigma));
	V2[(mm1 << LOG_SIGMA) + x[mm1]] = V3 + (mm1 << LOG_SIGMA);

	/*fillBeginning(V2m, 0, mm1, int_size_sigma, int_size);
	copy_value(V2m + mm1_sigma, SIGMA * m2m1, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
	V2m[(mm1 << LOG_SIGMA) + x[i]] = ((D_[x[i]] + m) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V2m + mm1_sigma));
	V2m[(mm1 << LOG_SIGMA) + x[mm1]] = (mm1 << LOG_SIGMA);*/

	// Filling V3 with pointers to chunks of V4
	fillBeginning(V3, V4, m, int_size_sigma, int_size);
	fillFirstLetter(V3 + m_sigma, V4 + m2_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
		for (int j = 0; j < SIGMA; j++)
			V3[((m + i) << LOG_SIGMA) + j] = V4 + ((m + BR_[x[mm2 - i]][j]) << LOG_SIGMA);

	/*fillBeginning(V3m, 0, m, int_size_sigma, int_size);
	fillFirstLetter(V3m + m_sigma, m2_sigma, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = 0; i < mm1; i++)
	for (int j = 0; j < SIGMA; j++)
	V3m[((m + i) << LOG_SIGMA) + j] = ((m + BR_[x[mm2 - i]][j]) << LOG_SIGMA);*/

	// Filling V4 with pointers to chunks of V5
	fillBeginning(V4, V5, m2m1, int_size_sigma, int_size);
	copy_value(V4 + (m2m1 << LOG_SIGMA), V5 + (m3m1 << LOG_SIGMA), int_size_sigma * mp1, int_size);
	for (int i = 0; i < mm1; i++)
		V4[(m2m1 << LOG_SIGMA) + x[i]] = V5 + ((D_[x[i]] + m2) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V4 + (m2m1 << LOG_SIGMA)));
	V4[(m2m1 << LOG_SIGMA) + x[mm1]] = V5 + (m2m1 << LOG_SIGMA);

	/*fillBeginning(V4m, 0, m2m1, int_size_sigma, int_size);
	copy_value(V4m + (m2m1 << LOG_SIGMA), (m3m1 << LOG_SIGMA), int_size_sigma * mm2, int_size); //mm2
	for (int i = 0; i < mm1; i++)
	V4m[(m2m1 << LOG_SIGMA) + x[i]] = ((D_[x[i]] + m2) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V4m + (m2m1 << LOG_SIGMA)));
	V4m[(m2m1 << LOG_SIGMA) + x[mm1]] = (m2m1 << LOG_SIGMA);*/

	// Filling V5 with pointers to chunks of V6
	fillBeginning(V5, V6, m2, int_size_sigma, int_size);
	fillFirstLetter(V5 + m2_sigma, V6 + (m3 << LOG_SIGMA), int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = m; i < m2m1; i++)
		for (int j = 0; j < SIGMA; j++)
			V5[((m + i) << LOG_SIGMA) + j] = V6 + ((BR_[x[m2m2 - i]][j] + m2) << LOG_SIGMA);

	/*fillBeginning(V5m, 0, m2, int_size_sigma, int_size);
	fillFirstLetter(V5m + m2_sigma, m3<<LOG_SIGMA, int_size_sigma, int_size, x[0], int_size_sigma_m); //<<LOG_SIGMA
	for (int i = m; i < m2m1; i++)
	for (int j = 0; j < SIGMA; j++)
	V5m[((m + i) << LOG_SIGMA) + j] = ((BR_[x[m2m2 - i]][j] + m2) << LOG_SIGMA);*/

	// Filling V6 with pointers to chunks of V7
	fillBeginning(V6, V7, m3m1, int_size_sigma, int_size);
	copy_value(V6 + (m3m1 << LOG_SIGMA), V7 + (m4m1 << LOG_SIGMA), int_size_sigma * mp1, int_size); // why mp1
	for (int i = 0; i < mm1; i++)
		V6[(m3m1 << LOG_SIGMA) + x[i]] = V7 + ((D_[x[i]] + m3) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V6 + (m3m1 << LOG_SIGMA)));
	V6[(m3m1 << LOG_SIGMA) + x[mm1]] = V7 + (m3m1 << LOG_SIGMA);

	/*fillBeginning(V6m, 0, m3m1, int_size_sigma, int_size);
	copy_value(V6m + (m3m1 << LOG_SIGMA), (m4m1 << LOG_SIGMA), int_size_sigma * mm2, int_size); // why mp1
	for (int i = 0; i < mm1; i++)
	V6m[(m3m1 << LOG_SIGMA) + x[i]] = ((D_[x[i]] + m3) << LOG_SIGMA);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V6m + (m3m1 << LOG_SIGMA)));
	V6m[(m3m1 << LOG_SIGMA) + x[mm1]] = (m3m1 << LOG_SIGMA);*/

	// Filling V7 with shift values
	fillBeginningFinal(V7, m3, int_size_sigma, int_size);
	fillFirstLetterFinal(V7 + m3_sigma, m4, int_size_sigma, int_size, x[0], int_size_sigma_m);
	for (int i = m; i < m2m1; i++)
		for (int j = 0; j < SIGMA; j++)
			V7[((m2 + i) << LOG_SIGMA) + j] = BR_[x[m2m2 - i]][j] + m3;

	//Search
	int *******p1, ******p2, *****p3, ****p4, ***p5, **p6, *p7;
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		p1 = V0[y[pos]];
		p2 = p1[y[pos + 1]];
		p3 = p2[y[pos + m]];
		p4 = p3[y[pos + mp1]];
		p5 = p4[y[pos + m2]];
		p6 = p5[y[pos + m2p1]];
		p7 = p6[y[pos + m3]];
		r = p7[y[pos + m3p1]];
		if (!r) {
			for (k = 0; k < mm2 && y[pos - mm2 + k] == x[k]; k++);
			if (k == mm2) {
				if (pos >= n)
					break;
				count++;
			}
			pos += D[y[pos + 1]];
		}
		else
			pos += r;
	}

	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw42p += u;

	return count;
}

//-----------------------------------TESTING-----------------------------------------------

long long sum_qlqs, sum_hash3, sum_ebom, sum_tvsbs, sum_fsbndm, sum_sa, sum_sbndmq2, sum_sbndmq4, sum_gsbndmq2, sum_fsb31, sum_fsb41, sum_fsb51, sum_bsdm;
int qlqs, hash3, ebom, tvsbs, fsbndm, sa, sbndmq2, sbndmq4, gsbndmq2, fsb31, fsb41, fsb51, bsdm;

int TD1[1000];
void build_TD1() {
	int i;
	unsigned char *p;
	for (i = 0; i < SIGMA; i++)
		TD1[i] = m + 1;
	for (p = (unsigned char*)P, i = 0; i < m; p++, i++)
		TD1[*p] = m - i;
}

int shb[1000];
//=============== QLQS
static int qlqsSearch(const unsigned char *P, int m, const unsigned char *T, int n) {
	QueryPerformanceCounter(&start);
	int k, i, shb[SIGMA], z = 2 * m + 1, zm1 = z - 1, nm = n - m;
	build_TD1();
	for (i = 0; i < SIGMA; i++)
		shb[i] = z - TD1[i];
	i = 0;
	int count_qlqs = 0;
	while (i <= nm) {
		k = 0;
		while (k<m && P[k] == T[i + k])
			k++;
		if (k == m)
			count_qlqs++;
		i += (TD1[T[i + m]] > shb[T[i + zm1]]
			? z
			: TD1[T[i + m]]);
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;

	sum_qlqs += u;
	return count_qlqs;
}

#define WSIZE  	    256				//greater int value fitting in a computer word
#define RANK3 3 
#define OUTPUT(j)   count++

//=============== HASH3
int searchH3(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 3) return -1;
	QueryPerformanceCounter(&start);
	int count, j, i, sh, sh1, mMinus1, mMinus2, shift[WSIZE];
	unsigned char h;
	count = 0;
	mMinus1 = m - 1;
	mMinus2 = m - 2;

	/* Preprocessing */
	for (i = 0; i < WSIZE; ++i)
		shift[i] = mMinus2;

	h = x[0];
	h = ((h << 1) + x[1]);
	h = ((h << 1) + x[2]);
	shift[h] = m - RANK3;
	for (i = RANK3; i < mMinus1; ++i) {
		h = x[i - 2];
		h = ((h << 1) + x[i - 1]);
		h = ((h << 1) + x[i]);
		shift[h] = mMinus1 - i;
	}
	h = x[i - 2];
	h = ((h << 1) + x[i - 1]);
	h = ((h << 1) + x[i]);
	sh1 = shift[h];
	shift[h] = 0;
	if (sh1 == 0) sh1 = 1;


	/* Searching */
	i = mMinus1;
	memcpy(y + n, x, m);
	while (1) {
		sh = 1;
		while (sh != 0) {
			h = y[i - 2];
			h = ((h << 1) + y[i - 1]);
			h = ((h << 1) + y[i]);
			sh = shift[h];
			i += sh;
		}
		if (i < n) {
			j = 0;
			while (j < m && x[j] == y[i - mMinus1 + j]) j++;
			if (j >= m) {
				OUTPUT(i - mMinus1);
			}
			i += sh1;
		}
		else {
			QueryPerformanceCounter(&_end);
			u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
			sum_hash3 += u;
			return count;
		}
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_hash3 += u;
	return count;
}

#define XSIZE       210			//maximal length of the pattern
#define UNDEFINED       -1

//=============== EBOW
int ebomSearch(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start);

	int S[XSIZE], FT[SIGMA][SIGMA];
	int *trans[XSIZE];
	int i, j, p, q;
	int iMinus1, mMinus1, count;
	unsigned char c;
	count = 0;

	// Preprocessing
	for (i = 0; i <= m + 1; i++) trans[i] = (int *)malloc(sizeof(int)*(SIGMA));
	for (i = 0; i <= m + 1; i++) for (j = 0; j<SIGMA; j++) trans[i][j] = UNDEFINED;
	S[m] = m + 1;
	for (i = m; i > 0; --i) {
		iMinus1 = i - 1;
		c = x[iMinus1];
		trans[i][c] = iMinus1;
		p = S[i];
		while (p <= m && (q = trans[p][c]) == UNDEFINED) {
			trans[p][c] = iMinus1;
			p = S[p];
		}
		S[iMinus1] = (p == m + 1 ? m : q);
	}

	/* Construct the FirstTransition table */
	for (i = 0; i < SIGMA; i++) {
		q = trans[m][i];
		for (j = 0; j < SIGMA; j++)
			if (q >= 0) FT[i][j] = trans[q][j];
			else FT[i][j] = UNDEFINED;
	}

	/* Searching */
	for (i = 0; i < m; i++) y[n + i] = x[i];
	if (!memcmp(x, y, m)) count++;
	j = m;
	mMinus1 = m - 1;
	while (j < n) {
		while ((FT[y[j]][y[j - 1]]) == UNDEFINED) j += mMinus1;
		i = j - 2;
		p = FT[y[j]][y[j - 1]];
		while ((p = trans[p][y[i]]) != UNDEFINED) i--;
		if (i < j - mMinus1/* && j<n*/) {
			count++;
			i++;
		}
		j = i + m;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_ebom += u;
	for (i = 0; i <= m + 1; i++) free(trans[i]);
	return count;
}

void TVSBSpreBrBc(unsigned char *x, int m, int brBc[SIGMA][SIGMA]) {
	int a, b, i;
	for (a = 0; a < SIGMA; ++a)
		for (b = 0; b < SIGMA; ++b)
			brBc[a][b] = m + 2;
	for (a = 0; a < SIGMA; ++a)
		brBc[a][P1[0]] = m + 1;
	for (i = 0; i < m - 1; ++i)
		brBc[P1[i]][P1[i + 1]] = m - i;
	for (a = 0; a < SIGMA; ++a)
		brBc[P1[m - 1]][a] = 1;
}

//=============== TVSBS
int TVSBSsearch(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start);
	int count, i, j = 0;
	int BrBc[SIGMA][SIGMA];
	unsigned char firstCh, lastCh;
	int ktvsbs = 0;
	TVSBSpreBrBc(x, m, BrBc);
	firstCh = x[0];
	lastCh = x[m - 1];
	for (i = 0; i<m; i++) y[n + i] = y[n + m + i] = x[i];
	while (j <= n - m) {
		if (lastCh == y[j + m - 1] && firstCh == y[j]) {
			for (i = m - 2; i > 0 && x[i] == y[j + i]; i--);
			if (i <= 0) ktvsbs++;
		}
		j += BrBc[y[j + m]][y[j + m + 1]];
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_tvsbs += u;
	return ktvsbs;
}

//=============== FSBNDM
int FSBNDMsearch(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start);

	unsigned int B[SIGMA], D, set;
	int i, j, pos, mMinus1, count;

	//   if (m>31) return search_large(x,m,y,n); 
	/* Preprocessing */
	count = 0;
	mMinus1 = m - 1;
	set = 1;
	for (i = 0; i < SIGMA; i++) B[i] = set;
	for (i = 0; i < m; ++i) B[x[i]] |= (1 << (m - i));

	/* Searching */
	if (!memcmp(x, y, m)) OUTPUT(0);
	j = m;
	while (j < n) {
		D = (B[y[j + 1]] << 1) & B[y[j]];
		if (D != 0) {
			pos = j;
			while (D = (D << 1) & B[y[j - 1]]) --j;
			j += mMinus1;
			if (j == pos) {
				OUTPUT(j);
				++j;
			}
		}
		else j += m;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_fsbndm += u;
	return count;
}

void preSA(unsigned char *x, int m, unsigned int S[]) {
	unsigned int j, lim;
	int i;
	for (i = 0; i < SIGMA; ++i) S[i] = 0;
	for (i = 0, j = 1; i < m; ++i, j <<= 1) {
		S[x[i]] |= j;
	}
}

//=============== Shift-And
int searchSA(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start);

	unsigned int lim, D;
	unsigned int S[SIGMA], F;
	int j, count;
	//   if (m > WORD) return search_large(x,m,y,n);

	/* Preprocessing */
	preSA(x, m, S);
	F = 1 << (m - 1);

	/* Searching */
	count = 0;
	for (D = 0, j = 0; j < n; ++j) {
		D = ((D << 1) | 1) & S[y[j]];
		if (D & F) OUTPUT(j - m + 1);
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_sa += u;
	return count;
}

#define GRAM2(j) (B[y[j]]<<1)&B[y[j-1]]

//============== SBNDMq2
int searchSBNDMq2(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start);
	unsigned int B[SIGMA], D, q;
	int i, j, pos, mMinusq, mq, count, shift;
	q = 2;
	if (m < q) return -1;
	//   if (m>32) return search_large(x,m,y,n);

	/* Preprocessing */
	count = 0;
	mMinusq = m - q + 1;
	mq = m - q;
	for (i = 0; i < SIGMA; i++) B[i] = 0;
	for (i = 1; i <= m; ++i)
		B[x[m - i]] |= (1 << (i - 1));

	D = B[x[m - 2]]; j = 1; shift = 0;
	if (D & (1 << (m - 1))) shift = m - j;
	for (i = m - 3; i >= 0; i--) {
		D = (D << 1) & B[x[i]];
		j++;
		if (D & (1 << (m - 1))) shift = m - j;
	}

	/* Searching */
	if (!memcmp(x, y, m)) OUTPUT(0);
	j = m;
	while (j < n) {
		D = GRAM2(j);
		if (D != 0) {
			pos = j;
			while (D = (D << 1) & B[y[j - q]]) --j;
			j += mq;
			if (j == pos) {
				OUTPUT(j);
				j += shift;
			}
		}
		else j += mMinusq;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_sbndmq2 += u;
	return count;
}

#define GRAM4(j) (B[y[j]]<<3)&(B[y[j-1]]<<2)&(B[y[j-2]]<<1)&B[y[j-3]]

int search_large_SBNDM4(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start);
	unsigned int B[SIGMA], D, q, shift;
	int i, j, pos, mMinusq, mq, count, p_len;
	q = 4;
	if (m < q) return 0;
	p_len = m;
	m = 32;
	int diff = p_len - m;

	/* Preprocessing */
	int kfsbndm4 = 0;
	mMinusq = m - q + 1;
	mq = m - q;
	for (i = 0; i < SIGMA; i++) B[i] = 0;
	for (i = 1; i <= m; ++i)
		B[x[m - i]] |= (1 << (i - 1));

	D = B[x[m - 2]]; j = 1; shift = 0;
	if (D & (1 << (m - 1))) shift = m - j;
	for (i = m - 3; i >= 0; i--) {
		D = (D << 1) & B[x[i]];
		j++;
		if (D & (1 << (m - 1))) shift = m - j;
	}

	/* Searching */
	if (!memcmp(x, y, p_len)) kfsbndm4++;
	j = m;
	while (j + diff < n) {
		D = GRAM4(j);
		if (D != 0) {
			pos = j;
			while (D = (D << 1) & B[y[j - q]]) --j;
			j += mq;
			if (j == pos) {
				for (i = m + 1; i < p_len && x[i] == y[j - m + 1 + i]; i++);
				if (i == p_len) kfsbndm4++;
				j += shift;
			}
		}
		else j += mMinusq;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_sbndmq4 += u;
	return kfsbndm4;
}

//================= SBNDMq4
int searchSBNDMq4(unsigned char *x, int m, unsigned char *y, int n) {
	unsigned int q = 4;

	if (m<q) return -1;
	if (m>32) return search_large_SBNDM4(x, m, y, n);

	QueryPerformanceCounter(&start);

	unsigned int B[SIGMA], D;
	int i, j, pos, mMinusq, mq, count, shift;

	/* Preprocessing */
	int kfsbndm4 = 0;
	mMinusq = m - q + 1;
	mq = m - q;
	for (i = 0; i < SIGMA; i++) B[i] = 0;
	for (i = 1; i <= m; ++i)
		B[x[m - i]] |= (1 << (i - 1));

	D = B[x[m - 2]]; j = 1; shift = 0;
	if (D & (1 << (m - 1))) shift = m - j;
	for (i = m - 3; i >= 0; i--) {
		D = (D << 1) & B[x[i]];
		j++;
		if (D & (1 << (m - 1))) shift = m - j;
	}

	/* Searching */
	if (!memcmp(x, y, m)) kfsbndm4++;
	j = m;
	while (j < n) {
		D = GRAM4(j);
		if (D != 0) {
			pos = j;
			while (D = (D << 1) & B[y[j - q]]) --j;
			j += mq;
			if (j == pos) {
				kfsbndm4++;
				j += shift;
			}
		}
		else j += mMinusq;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_sbndmq4 += u;
	return kfsbndm4;
}

//============== GSBNDMq2
int GSBNDMq2(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 4) return -1;
	QueryPerformanceCounter(&start);
	unsigned int B[SIGMA], D, q, mpm2 = 2 * m - 2, mm1 = m - 1, mm2 = m - 2;
	int i, j, pos, mMinusq, mq, shift;
	q = 2;
	int kgsb2 = 0;
	//   if (m>32) return search_large(x,m,y,n);

	/* Preprocessing */
	mMinusq = m - q + 1;
	mq = m - q;
	for (i = 0; i < SIGMA; i++) B[i] = 0;
	for (i = 0; i < m; ++i)
		B[x[i]] |= (1 << (mm1 - i));

	/* Searching */
	if (!memcmp(x, y, m)) kgsb2++;
	for (int i = 0; i < m; i++)
		y[N - m * 2 + i] = P[i];
	i = m - 2;
	while (i < n) {
		while ((D = ((B[y[i + 1]] << 1)&B[y[i]])) == 0 && ((B[y[i + m]] << 1)&B[y[i + mm1]]) == 0)
			i += mpm2;
		if (D != 0) {
			j = i;
			do {
				i--;
				D = (D << 1)&B[T[i]];
			} while (D != 0);
			i += mm1;
			if (j == i) {
				kgsb2++;
				i++;
			}
		}
		else
			i += mm2;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_gsbndmq2 += u;
	return kgsb2;
}

const int Qq = 3, F = 1, W = 32;

//============== FSBNDM31
int FSBNDMsearch31(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start);
	unsigned int B[SIGMA], D;
	int i, j, pos, mm1 = m - 1, count, nq = n - Qq + 1, mqf = m - Qq + F + 1;

	//   if (m>31) return search_large(x,m,y,n); 

	/* Preprocessing */
	int kfsbqf = 0;
	for (i = 0; i < SIGMA; i++) B[i] = 1;
	for (i = 0; i < m; ++i) B[x[i]] |= (1 << (mm1 - i + F));

	/* Searching */
	if (!memcmp(x, y, m)) kfsbqf++;
	i = mm1 - Qq + F;
	while (i < nq) {
		D = B[T[i]] & (B[T[i + 1]] << 1)&(B[T[i + 2]] << 2);
		if (D != 0) {
			pos = i - mqf;
			do {
				i--;
				D = (D << 1)&B[T[i]];
			} while (D != 0);
			if (i == pos) {
				kfsbqf++;
				++i;
			}
		}
		i += mqf;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_fsb31 += u;
	return kfsbqf;
}

const int Q41 = 4, F41 = 1;

//============== FSBNDM41
int FSBNDMsearch41(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 3) return -1;
	QueryPerformanceCounter(&start);

	unsigned int B[SIGMA], D;
	int i, j, pos, mm1 = m - 1, count, nq = n - Q41 + 1, mqf = m - Q41 + F41 + 1;

	//   if (m>31) return search_large(x,m,y,n); 

	/* Preprocessing */
	int kfsbq41 = 0;
	for (i = 0; i < SIGMA; i++) B[i] = 1;
	for (i = 0; i < m; ++i) B[x[i]] |= (1 << (mm1 - i + F41));

	/* Searching */
	if (!memcmp(x, y, m)) kfsbq41++;
	i = mm1 - Q41 + F41;
	while (i < nq) {
		D = B[T[i]] & (B[T[i + 1]] << 1)&(B[T[i + 2]] << 2)&(B[T[i + 3]] << 3);
		if (D != 0) {
			pos = i - mqf;
			do {
				i--;
				D = (D << 1)&B[T[i]];
			} while (D != 0);
			if (i == pos) {
				kfsbq41++;
				++i;
			}
		}
		i += mqf;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_fsb41 += u;
	return kfsbq41;
}

const int Q51 = 4, F51 = 1;

//============== FSBNDM41
int FSBNDMsearch51(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 4) return -1;
	//   if (m>31) return search_large(x,m,y,n); 
	QueryPerformanceCounter(&start);

	unsigned int B[SIGMA], D;
	int i, j, pos, mm1 = m - 1, count, nq = n - Q51 + 1, mqf = m - Q51 + F51 + 1;

	/* Preprocessing */
	int kfsbq51 = 0;
	for (i = 0; i < SIGMA; i++) B[i] = 1;
	for (i = 0; i < m; ++i) B[x[i]] |= (1 << (mm1 - i + F51));

	/* Searching */
	if (!memcmp(x, y, m)) kfsbq51++;
	i = mm1 - Q51 + F51;
	while (i < nq) {
		D = B[T[i]] & (B[T[i + 1]] << 1)&(B[T[i + 2]] << 2)&(B[T[i + 3]] << 3);
		if (D != 0) {
			pos = i - mqf;
			do {
				i--;
				D = (D << 1)&B[T[i]];
			} while (D != 0);
			if (i == pos) {
				kfsbq51++;
				++i;
			}
		}
		i += mqf;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_fsb51 += u;
	return kfsbq51;
}

#define DSIGMA 65536
#define HS(x,i) (x[i]<<6) + (x[i+1]<<4) +  (x[i+2]<<2) + x[i+3]
#define Q 4

//============== BSDM
int searchBSDM(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < Q) return -1;
	QueryPerformanceCounter(&start);

	unsigned int B[DSIGMA];
	int i, j, k, count = 0;
	unsigned int s, d;

	/* Preprocessing */
	unsigned int occ[DSIGMA] = { 0 };
	int star = 0, len = 0;
	unsigned short c, dch;
	for (i = 0, j = 0; i < m - Q + 1; i++) {
		c = HS(x, i);
		if (occ[c]) {
			dch = HS(x, j);
			while (dch != c) {
				occ[dch] = 0;
				j++;
				dch = HS(x, j);
			}
			occ[dch] = 0;
			j++;
		}
		occ[c] = 1;
		if (len < i - j + 1) {
			len = i - j + 1;
			star = j;
		}
	}

	unsigned int pos[DSIGMA];
	for (i = 0; i < DSIGMA; i++) pos[i] = -1;
	for (i = 0; i < len; i++) {
		c = HS(x, star + i);
		pos[c] = i;
	}

	/* Searching */
	for (i = 0; i < m; i++) y[n + i] = x[i];
	unsigned char *xstart = x + star;
	int offset = len + star - 1;
	j = len - 1;
	while (j < n) {
		c = HS(y, j);
		while ((i = pos[c]) < 0) {
			j += len;
			c = HS(y, j);
		}
		k = 1;
		while (k <= i && xstart[i - k] == y[j - k]) k++;
		if (k > i) {
			if (k == len) {
				if (!memcmp(x, y + j - offset, m)) if (j - offset <= n - m) count++;
			}
			else j -= k;
		}
		j += len;
	}
	QueryPerformanceCounter(&_end);
	u = (_end.QuadPart - start.QuadPart) * 1000000 / freq.QuadPart;
	sum_bsdm += u;
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

void English() {
	N = 4912296;
	FILE *f;
	f = fopen("bible.txt", "rt");
	fread(T, 1, N, f);
}

#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds

void main() {
	QueryPerformanceFrequency(&freq);
	QueryPerformanceFrequency(&_freq);

	f = fopen("output.csv", "wt");
	generateRandom();
	//DNA();

	fprintf(f, "b=%d N=%d ITER=%d\n", SIGMA, N, ITER);

	fprintf(f, "m,MAW22,MAW22P,MAW23P,MAW24P,MAW32P,MAW33P,MAW42P,QLQS,HASH3,EBOW,TVSBS,FSBNDM,SA,SBNDMq2,SBNDMq4,GSBNDMq2,FSBNDM31,FSBNDM41,FSBNDM51,BSDM");
	for (m = 2; m < 71; m < 10 ? m++ : m += 10) {
		//for (m = 16; m < 33; m+=2) {

		//for (int ig = 0; ig < 2; ig++) {
			sum_maw22 = sum_maw23 = sum_maw24 = sum_maw32 = sum_maw33 = 0;
			sum_maw22p = sum_maw23p = sum_maw24p = sum_maw32p = sum_maw33p = sum_maw42p = 0;
			sum_qlqs = sum_hash3 = sum_ebom = sum_tvsbs = sum_fsbndm = sum_sa = sum_sbndmq2 = sum_sbndmq4 = sum_gsbndmq2 = sum_fsb31 = sum_fsb41 = sum_fsb51 = sum_bsdm = 0;
			
			nm2 = N - 2 * m;
			int nm = N - m;
			memcpy(T1, T, N);
			for (int ii = 0; ii < ITER; ii++) {
				srand((unsigned)time(NULL));
				int patpos = rand() % (N - m - 2);
				for (int i = 0; i < m; i++)
					P[i] = T[patpos + i];
				for (int i = 0; i < m; i++)
					T[N - m + i] = P[i];

				memcpy(P1, P, m);

				maw22 = MAW22(P, m, T, nm);
				maw22p = MAW22P(P, m, T, nm);
				maw23p = MAW23P(P, m, T, nm);
				maw24p = MAW24P(P, m, T, nm);
				maw32p = MAW32P(P, m, T, nm);
				maw33p = MAW33P(P, m, T, nm);
				maw42p = MAW42P(P, m, T, nm);

				qlqs = qlqsSearch(P, m, T, N);
				hash3 = searchH3(P, m, T, N);
				ebom = ebomSearch(P, m, T, N);
				tvsbs = TVSBSsearch(P, m, T, N);
				fsbndm = FSBNDMsearch(P, m, T, N);
				sa = searchSA(P, m, T, N);
				sbndmq2 = searchSBNDMq2(P, m, T, N);
				sbndmq4 = searchSBNDMq4(P, m, T, N);
				gsbndmq2 = GSBNDMq2(P, m, T, N);
				fsb31 = FSBNDMsearch31(P, m, T, N);
				fsb41 = FSBNDMsearch41(P, m, T, N);
				fsb51 = FSBNDMsearch51(P, m, T, N);
				bsdm = searchBSDM(P, m, T, N);
			}
			printf("b=%d m=%d\n", SIGMA, m);
			printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", maw22, maw22p, maw23p, maw24p, maw32p, maw33p, maw42p, qlqs, hash3, ebom, tvsbs, fsbndm, sa, sbndmq2, sbndmq4, gsbndmq2, fsb31, fsb41, fsb51, bsdm);
			fprintf(f, "\n%2.d,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld",
				m, sum_maw22, sum_maw22p, sum_maw23p, sum_maw24p, sum_maw32p, sum_maw33p, sum_maw42p, sum_qlqs, sum_hash3, sum_ebom, sum_tvsbs, sum_fsbndm, sum_sa, sum_sbndmq2, sum_sbndmq4,
				sum_gsbndmq2, sum_fsb31, sum_fsb41, sum_fsb51, sum_bsdm);
			printf("%2.d,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld\n\n",
				m, sum_maw22, sum_maw22p, sum_maw23p, sum_maw24p, sum_maw32p, sum_maw33p, sum_maw42p, sum_qlqs, sum_hash3, sum_ebom, sum_tvsbs, sum_fsbndm, sum_sa, sum_sbndmq2, sum_sbndmq4,
				sum_gsbndmq2, sum_fsb31, sum_fsb41, sum_fsb51, sum_bsdm);

		//}
	}
	fclose(f);
	system("pause");
}