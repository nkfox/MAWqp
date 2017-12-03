#define _CRT_SECURE_NO_WARNINGS
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
/*#define SIGMA 8		//alphabet size
#define SIGMA2 SIGMA*SIGMA
#define SIGMA3 SIGMA*SIGMA2
#define SIGMA4 SIGMA*SIGMA3
#define SIGMA5 SIGMA*SIGMA4
#define SIGMA6 SIGMA*SIGMA5
#define SIGMA7 SIGMA*SIGMA6
#define SIGMA8 SIGMA*SIGMA7
#define SIGMA9 SIGMA*SIGMA8*/


/**/#define SIGMA 4		//alphabet size
#define SIGMA2 16
#define SIGMA3 64
#define SIGMA4 256
#define SIGMA5 1024
#define SIGMA6 4096
#define SIGMA7 16384
#define SIGMA8 65536
#define SIGMA9 262144

/*#define SIGMA 6		//alphabet size
#define SIGMA2 36
#define SIGMA3 216
#define SIGMA4 1296
#define SIGMA5 7776
#define SIGMA6 46656
#define SIGMA7 279936
#define SIGMA8 1679616
#define SIGMA9 10077696

#define SIGMA 8		//alphabet size
#define SIGMA2 64
#define SIGMA3 512
#define SIGMA4 4096
#define SIGMA5 32768
#define SIGMA6 262144
#define SIGMA7 2097152
#define SIGMA8 16777216
#define SIGMA9 134217728*/



const int TOTAL = 10000200;
unsigned char T[TOTAL], T1[TOTAL], P[200], P1[200];
int N = TOTAL - 200, ITER = 200, m = 5;

FILE * f;
LARGE_INTEGER start, _end, freq, _freq, prep_start, prep_end;
double u;
int nm2, glob = 0;

long long sum_maw22, sum_maw23, sum_maw24, sum_maw32, sum_maw33;
int maw22, maw23, maw24, maw32, maw33;

//a bytes pointed by c are repeated until c[0..b] is filled
void mem_fill(int a, int b, unsigned char* c) {
	int i;
	for (i = a; i <= (b >> 1); i <<= 1)
		memcpy(c + i, c, i);
	memcpy(c + i, c, b - i);
}

void copy_value(unsigned char* pointer, int value, int length)
{
	*pointer = value;
	mem_fill(1, length, pointer);
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
int TVSBSsearch(unsigned char *x, int m, unsigned char *y, int n){
	QueryPerformanceCounter(&start);
	int count, i, j = 0;
	int BrBc[SIGMA][SIGMA];
	unsigned char firstCh, lastCh;
	int ktvsbs = 0;
	TVSBSpreBrBc(x, m, BrBc);
	firstCh = x[0];
	lastCh = x[m - 1];
	for (i = 0; i<m; i++) y[n + i] = y[n + m + i] = x[i];
	while (j <= n - m){
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

void main(){
	QueryPerformanceFrequency(&freq);
	QueryPerformanceFrequency(&_freq);

	f = fopen("output.csv", "wt");
	generateRandom();
	//DNA();

	fprintf(f, "b=%d N=%d ITER=%d\n", SIGMA, N, ITER);

	fprintf(f, "m,MAW22,MAW23,MAW24,MAW32,MAW33,QLQS,HASH3,EBOW,TVSBS,FSBNDM,SA,SBNDMq2,SBNDMq4,GSBNDMq2,FSBNDM31,FSBNDM41,FSBNDM51,BSDM");
	for (m = 2; m < 16; m++) {
	//for (m = 16; m < 33; m+=2) {

		for (int ig = 0; ig < 2; ig++) {
			sum_maw22 = sum_maw23 = sum_maw24 = sum_maw32 = sum_maw33 = sum_qlqs = sum_hash3 = sum_ebom = sum_tvsbs = sum_fsbndm =
				sum_sa = sum_sbndmq2 = sum_sbndmq4 = sum_gsbndmq2 = sum_fsb31 = sum_fsb41 = sum_fsb51 = sum_bsdm = 0;
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

				maw22 = MAW22(P, m, T, nm2);
				maw23 = MAW23(P, m, T, nm2);
				maw24 = MAW24(P, m, T, nm2);
				maw32 = MAW32(P, m, T, nm2);
				maw33 = MAW33(P, m, T, nm2);

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
			printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", maw22, maw23, maw24, maw32, maw33, qlqs, hash3, ebom, tvsbs, fsbndm, sa, sbndmq2, sbndmq4, gsbndmq2, fsb31, fsb41, fsb51, bsdm);
			fprintf(f, "\n%2.d,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld",
				m, sum_maw22, sum_maw23, sum_maw24, sum_maw32, sum_maw33, sum_qlqs, sum_hash3, sum_ebom, sum_tvsbs, sum_fsbndm, sum_sa, sum_sbndmq2, sum_sbndmq4,
				sum_gsbndmq2, sum_fsb31, sum_fsb41, sum_fsb51, sum_bsdm);
			printf("%2.d,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld\n\n",
				m, sum_maw22, sum_maw23, sum_maw24, sum_maw32, sum_maw33, sum_qlqs, sum_hash3, sum_ebom, sum_tvsbs, sum_fsbndm, sum_sa, sum_sbndmq2, sum_sbndmq4,
				sum_gsbndmq2, sum_fsb31, sum_fsb41, sum_fsb51, sum_bsdm);

		}
	}
	fclose(f);
	system("pause");
}

