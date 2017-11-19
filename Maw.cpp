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
#define SIGMA 4		//alphabet size
#define SIGMA2 SIGMA*SIGMA
#define SIGMA3 SIGMA*SIGMA2
#define SIGMA4 SIGMA*SIGMA3
#define SIGMA5 SIGMA*SIGMA4
#define SIGMA6 SIGMA*SIGMA5
//#define LOG_SIGMA 3 //logarithm of alphabet size

//a bytes pointed by c are repeated until c[0..b] is filled
void mem_fill(int a, int b, unsigned char* c) {
	int i;
	b >>= 1;
	for (i = a; i < b; i <<= 1)
		memcpy(c + i, c, i);
	memcpy(c + i, c, b * 2 - i);
}

// The MAW32 algorithm
int MAW32(unsigned char *x, int m, unsigned char *y, int n) {
	int mp1 = m + 1, mm2 = m - 2, mm1 = m - 1, m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m3 = m * 3, m3m1 = 3 * m - 1, k, r, pos, count = 0;
	int D[P_MAX];
	unsigned char* M32 = (unsigned char *)malloc(SIGMA6); // MAW32 search table

	// Preprocessing
	// Build the BMH shift table D
	for (int i = 0; i < SIGMA; i++)
		D[i] = m;
	for (int i = 0; i < m - 1; i++)
		D[x[i]] = m - i - 1;

	// Fill the M22 search table
	*M32 = m3;
	mem_fill(1, SIGMA6, M32);

	*(M32 + x[0]) = m3m1;
	mem_fill(SIGMA, SIGMA6, M32);

	for (int k = 0; k < mm1; k++)
		*(M32 + x[k] * SIGMA + x[k + 1]) = m3 - k - 2;
	mem_fill(SIGMA2, SIGMA6, M32);

	*(M32 + x[0] * SIGMA2 + x[mm1] * SIGMA) = m2m1;
	mem_fill(1, SIGMA, M32 + x[0] * SIGMA2 + x[mm1] * SIGMA);
	mem_fill(SIGMA3, SIGMA6, M32);

	/*for (int k = 0; k<mm1; k++)
		*(M32 + x[k] * SIGMA + x[k + 1]) = m2 - k - 2;
		mem_fill(SIGMA2, SIGMA4, M32);*/
	for (int k = 0; k < mm1; k++)
	for (int a = 0; a < SIGMA; a++)
	for (int b = 0; b < SIGMA; b++)
	for (int c = 0; c < SIGMA; c++)
	for (int d = 0; d < SIGMA; d++)
		*(M32 + a * SIGMA5 + b * SIGMA4 + x[k] * SIGMA3 + x[k + 1] * SIGMA2 + c * SIGMA + d) = m2 - k - 2;

	/**(M32 + x[0] * SIGMA4 + x[mm1] * SIGMA3) = mm1;
	mem_fill(1, SIGMA, M32 + x[0] * SIGMA4 + x[mm1] * SIGMA3);
	mem_fill(SIGMA4, SIGMA6, M32);*/
	for (int k = 0; k < mm1; k++)
	for (int a = 0; a < SIGMA; a++)
	for (int b = 0; b < SIGMA; b++)
	for (int c = 0; c < SIGMA; c++)
	for (int d = 0; d < SIGMA; d++)
		*(M32 + a * SIGMA5 + x[0] * SIGMA4 + x[mm1] * SIGMA3 + b * SIGMA2 + c * SIGMA + d) = mm1;

	for (int j = 0; j < mm1; j++) {
		*(M32 + x[j] * SIGMA5 + x[j + 1] * SIGMA4) = m - j - 2;
		mem_fill(1, SIGMA4, M32 + x[j] * SIGMA5 + x[j + 1] * SIGMA4);
	}

	//Search
	pos = mm2;
	for (int i = 0; i<m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		if (pos>n) break;
		r = *(M32 + y[pos] * SIGMA5 + y[pos + 1] * SIGMA4 + y[pos + m] * SIGMA3 + y[pos + mp1] * SIGMA2 + y[pos + m2] * SIGMA + y[pos + m2p1]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < m && y[pos + k - mm2] == x[k]; k++);
			if (k == m) {
				count++;
			}
			pos += D[y[pos + 1]];
		}
		else
			pos += r;
	}

	free(M32);
	return count;
}

// The MAW23 algorithm
int MAW23(unsigned char *x, int m, unsigned char *y, int n) {
	int mp1 = m + 1, mm2 = m - 2, mm1 = m - 1, m2 = 2 * m, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, k, r, pos, count = 0;
	int D[P_MAX];
	unsigned char* M23 = (unsigned char *)malloc(SIGMA6); // MAW23 search table

	// Preprocessing
	// Build the BMH shift table D
	for (int i = 0; i < SIGMA; i++)
		D[i] = m;
	for (int i = 0; i < m - 1; i++)
		D[x[i]] = m - i - 1;

	// Fill the M22 search table
	*M23 = m2;
	mem_fill(1, SIGMA6, M23);

	*(M23 + x[0]) = m2m1;
	mem_fill(SIGMA, SIGMA6, M23);

	*(M23 + x[0] * SIGMA + x[1]) = m2m2;
	mem_fill(SIGMA2, SIGMA6, M23);

	for (int k = 0; k < mm2; k++)
		*(M23 + x[k] * SIGMA2 + x[k + 1] * SIGMA + x[k + 2]) = m2 - k - 3;
	mem_fill(SIGMA3, SIGMA6, M23);

	*(M23 + x[0] * SIGMA3 + x[mm2] * SIGMA2 + x[mm1] * SIGMA) = mm1;
	mem_fill(1, SIGMA, M23 + x[0] * SIGMA3 + x[mm2] * SIGMA2 + x[mm1] * SIGMA);
	mem_fill(SIGMA4, SIGMA6, M23);

	/**(M23 + x[0] * SIGMA4 + x[1] * SIGMA3 + x[mm1] * SIGMA2) = mm2;
	mem_fill(1, SIGMA, M23 + x[0] * SIGMA4 + x[1] * SIGMA3 + x[mm1] * SIGMA2);
	mem_fill(SIGMA6, SIGMA6, M23);*/
	for (int a = 0; a < SIGMA; a++)
	for (int b = 0; b < SIGMA; b++)
	for (int c = 0; c < SIGMA; c++)
		*(M23 + a * SIGMA5 + x[0] * SIGMA4 + x[1] * SIGMA3 + x[mm1] * SIGMA2 + b * SIGMA + c) = mm2;


	for (int j = 0; j < mm2; j++) {
		*(M23 + x[j] * SIGMA5 + x[j + 1] * SIGMA4 + x[j + 2] * SIGMA3) = m - j - 3;
		mem_fill(1, SIGMA3, M23 + x[j] * SIGMA5 + x[j + 1] * SIGMA4 + x[j + 2] * SIGMA3);
	}

	//Search
	pos = mm2 - 1;
	for (int i = 0; i<m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		if (pos>n + m) break;
		r = *(M23 + y[pos] * SIGMA5 + y[pos + 1] * SIGMA4 + y[pos + 2] * SIGMA3 + y[pos + m] * SIGMA2 + y[pos + mp1] * SIGMA + y[pos + mp1 + 1]);
		if (!r) {
			if (pos >= n)
				break;
			for (k = 0; k < m && y[pos + k - mm2 + 1] == x[k]; k++);
			if (k == m) {
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	free(M23);
	return count;
}

// The MAW22 algorithm
int MAW22(unsigned char *x, int m, unsigned char *y, int n) {
	int mp1 = m + 1, mm2 = m - 2, mm1 = m - 1, m2 = 2 * m, k, r, pos, count = 0;
	int D[P_MAX];
	unsigned char* M22 = (unsigned char *)malloc(SIGMA4); // MAW22 search table

	// Preprocessing
	// Build the BMH shift table D
	for (int i = 0; i < SIGMA; i++)
		D[i] = m;
	for (int i = 0; i < m - 1; i++)
		D[x[i]] = m - i - 1;

	// Fill the M22 search table
	*M22 = m2;
	mem_fill(1, SIGMA4, M22);

	*(M22 + x[0]) = m2 - 1;
	mem_fill(SIGMA, SIGMA2, M22);

	for (int k = 0; k < mm1; k++)
		*(M22 + x[k] * SIGMA + x[k + 1]) = m2 - k - 2;
	mem_fill(SIGMA2, SIGMA4, M22);

	*(M22 + x[0] * SIGMA2 + x[mm1] * SIGMA) = mm1;
	mem_fill(1, SIGMA, M22 + x[0] * SIGMA2 + x[mm1] * SIGMA);
	mem_fill(SIGMA3, SIGMA4, M22);
	for (int j = 0; j < mm1; j++) {
		*(M22 + x[j] * SIGMA3 + x[j + 1] * SIGMA2) = m - j - 2;
		mem_fill(1, SIGMA2, M22 + x[j] * SIGMA3 + x[j + 1] * SIGMA2);
	}

	//Search
	pos = mm2;
	for (int i = 0; i < m; i++) y[n + i] = x[i]; //append the text with a stop pattern
	while (true) {
		r = *(M22 + y[pos] * SIGMA3 + y[pos + 1] * SIGMA2 + y[pos + m] * SIGMA + y[pos + mp1]);
		if (!r) {
			if (pos >= n)
				break;
			//for (k = 0; k < m && y[pos + k - mm2] == P[k]; k++);
			for (k = 0; k < m && y[pos + k - mm2] == x[k]; k++);
			if (k == m) {
				count++;
			}
			pos += D[y[pos]];
		}
		else
			pos += r;
	}

	free(M22);
	return count;
}

/*// The MAW22 algorithm with pointers
int MAW22P(unsigned char *x, int m, unsigned char *y, int n) {
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

//Search
int ***p1, **p2, *p3;
pos = mm2;
for (int i = 0; i<m; i++) y[n + i] = x[i]; //append the text with a stop pattern
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
}
return count;
}
*/

const int Amount = 200000;

void writeAlphabet8()
{
	FILE *f;
	f = fopen("Alphabet8.txt", "w");
	srand((unsigned)time(NULL));
	for (int i = 0; i < Amount; i++)
	{
		char c = rand() % 4;
		fprintf(f, "%c", c);
	}
	fclose(f);
}

unsigned char T[Amount * 10], P[20];

void check()
{
	int n = Amount, m = 4;
	FILE *f;
	f = fopen("Alphabet8.txt", "rt");
	fread(T, 1, n, f);
	fclose(f);
	P[0] = 3; P[1] = 1; P[2] = 1; P[3] = 0;
	T[0] = 3; T[1] = 1; T[2] = 1; T[3] = 0;
	T[5] = 3; T[6] = 1; T[7] = 1; T[8] = 0;
	int maw22 = MAW22(P, m, T, n);
	int maw32 = MAW32(P, m, T, n);
	int maw23 = MAW23(P, m, T, n);
	cout << maw22 << " " << maw32 << " " << maw23 << endl;
	system("pause");
}

int main()
{
	writeAlphabet8();
	check();
}
