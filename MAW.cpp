#include <malloc.h>
#include <windows.h>
using namespace std;

#define P_MAX 200	//maximum pattern length
#define V_MAX 20000 //maximum size of pointer and shift arrays

#define CURRENT_SIGMA 4

#if CURRENT_SIGMA == 4
#define SIGMA 4	//alphabet size
#define SIGMA2 16
#define SIGMA3 64
#define SIGMA4 256
#define SIGMA5 1024
#define SIGMA6 4096
#define SIGMA7 16384
#define SIGMA8 65536
#define SIGMA9 262144
#define LOG_SIGMA 2 //logarithm of alphabet size

#elif CURRENT_SIGMA == 6
#define SIGMA 6		//alphabet size
#define SIGMA2 36
#define SIGMA3 216
#define SIGMA4 1296
#define SIGMA5 7776
#define SIGMA6 46656
#define SIGMA7 279936
#define SIGMA8 1679616
#define SIGMA9 10077696
#define LOG_SIGMA 3 //logarithm of alphabet size

#elif CURRENT_SIGMA == 8
#define SIGMA 8	//alphabet size
#define SIGMA2 64
#define SIGMA3 512
#define SIGMA4 4096
#define SIGMA5 32768
#define SIGMA6 262144
#define SIGMA7 2097152
#define SIGMA8 16777216
#define SIGMA9 134217728
#define LOG_SIGMA 3 //logarithm of alphabet size
#endif

//-----------------Additional functions---------------------

//amount_to_copy bytes pointed by pointer are repeated until pointer[0..range] is filled
void mem_fill(int amount_to_copy, int range, unsigned char* pointer) {
	int i;
	for (i = amount_to_copy; i <= (range >> 1); i <<= 1)
		memcpy(pointer + i, pointer, i);
	memcpy(pointer + i, pointer, range - i);
}

template <class T, class U>
void set_array(T* pointer, U value, int length, int type_length = 1) {
	*pointer = value;
	mem_fill(type_length, length, (unsigned char*)pointer);
}

void build_BMH_shift_table(int* D, const unsigned char *x, const int& m, int end_index = -1) {
	int mm1 = m - 1, int_size = sizeof(int);
	if (end_index == -1) end_index = m;
	set_array(D, m, int_size * SIGMA, int_size);
	for (int i = 0; i < end_index; i++)
		D[x[i]] = mm1 - i;
}

//-----------------------MAW--------------------------------

// The MAW22 algorithm
int MAW22(unsigned char *x, int m, unsigned char *y, int n) {
	int mp1 = m + 1, mm1 = m - 1, mm2 = m - 2,
		m2 = 2 * m, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2,
		step, pos, last_pos, count = 0;
	unsigned char *y_pos;
	int D[P_MAX];
	unsigned char* M22 = (unsigned char *)malloc(SIGMA4); // MAW22 search table

	// Preprocessing
	build_BMH_shift_table(D, x, m, mm1);

	// Fill the M22 search table
	memset(M22, m2, SIGMA);

	*(M22 + x[0]) = m2m1;
	mem_fill(SIGMA, SIGMA2, M22);

	for (int k = 0; k < mm1; k++)
		*(M22 + x[k] * SIGMA + x[k + 1]) = m2m2 - k;
	mem_fill(SIGMA2, SIGMA3, M22);

	set_array(M22 + x[0] * SIGMA2 + x[mm1] * SIGMA, mm1, SIGMA);
	mem_fill(SIGMA3, SIGMA4, M22);

	for (int k = 0; k < mm1; k++)
		set_array(M22 + x[k] * SIGMA3 + x[k + 1] * SIGMA2, mm2 - k, SIGMA2);

	//Search
	pos = mm2, last_pos = n + mm2;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = *(M22 + *y_pos * SIGMA3 + *(y_pos + 1) * SIGMA2 + *(y_pos + m) * SIGMA + *(y_pos + mp1)))) {
			if (!memcmp(x, y + (pos - mm2), mm2)) {
				if (pos == last_pos)
					break;
				++count;
			}
			pos += D[*(y_pos + 1)];
		}
		else
			pos += step;
	}

	free(M22);
	return count;
}

// The MAW23 algorithm
int MAW23(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 3) return -1;
	
	int mp1 = m + 1, mp2 = m + 2, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3,
		m2 = 2 * m, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3,
		step, pos, last_pos, count = 0;
	unsigned char *y_pos;
	int D[P_MAX];
	unsigned char* M23 = (unsigned char*)malloc(SIGMA6); // MAW23 search table

	// Preprocessing
	build_BMH_shift_table(D, x, m, mm1);

	// Fill the M23 search table
	memset(M23, m2, SIGMA);

	*(M23 + x[0]) = m2m1;
	mem_fill(SIGMA, SIGMA2, M23);

	*(M23 + x[0] * SIGMA + x[1]) = m2m2;
	mem_fill(SIGMA2, SIGMA3, M23);

	for (int k = 0; k < mm2; k++)
		*(M23 + x[k] * SIGMA2 + x[k + 1] * SIGMA + x[k + 2]) = m2m3 - k;
	mem_fill(SIGMA3, SIGMA4, M23);

	set_array(M23 + x[0] * SIGMA3 + x[mm2] * SIGMA2 + x[mm1] * SIGMA, mm1, SIGMA);
	mem_fill(SIGMA4, SIGMA5, M23);

	set_array(M23 + x[0] * SIGMA4 + x[1] * SIGMA3 + x[mm1] * SIGMA2, mm2, SIGMA2);
	mem_fill(SIGMA5, SIGMA6, M23);

	for (int k = 0; k < mm2; k++)
		set_array(M23 + x[k] * SIGMA5 + x[k + 1] * SIGMA4 + x[k + 2] * SIGMA3, mm3 - k, SIGMA3);

	//Search
	pos = mm3, last_pos = n + mm3;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = *(M23 + *y_pos * SIGMA5 + *(y_pos + 1) * SIGMA4 + *(y_pos + 2) * SIGMA3 + *(y_pos + m) * SIGMA2 + *(y_pos + mp1) * SIGMA + *(y_pos + mp2)))) {
			if (!memcmp(x, y + (pos - mm3), mm3)) {
				if (pos == last_pos)
					break;
				++count;
			}
			pos += D[*(y_pos + 2)];
		}
		else
			pos += step;
	}

	free(M23);
	return count;
}

// The MAW24 algorithm
int MAW24(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 4) return -1;
	
	int mp1 = m + 1, mp2 = m + 2, mp3 = m + 3, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3, mm4 = m - 4,
		m2 = 2 * m, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3, m2m4 = 2 * m - 4,
		step, pos, last_pos, count = 0;
	unsigned char *y_pos;
	int D[P_MAX];
	unsigned char* M24 = (unsigned char *)malloc(SIGMA8); // MAW24 search table

	// Preprocessing
	build_BMH_shift_table(D, x, m, mm1);

	// Fill the M24 search table
	memset(M24, m2, SIGMA);

	*(M24 + x[0]) = m2m1;
	mem_fill(SIGMA, SIGMA2, M24);

	*(M24 + x[0] * SIGMA + x[1]) = m2m2;
	mem_fill(SIGMA2, SIGMA3, M24);

	*(M24 + x[0] * SIGMA2 + x[1] * SIGMA + x[2]) = m2m3;
	mem_fill(SIGMA3, SIGMA4, M24);

	for (int k = 0; k < mm3; k++)
		*(M24 + x[k] * SIGMA3 + x[k + 1] * SIGMA2 + x[k + 2] * SIGMA + x[k + 3]) = m2m4 - k;
	mem_fill(SIGMA4, SIGMA5, M24);

	set_array(M24 + x[0] * SIGMA4 + x[mm3] * SIGMA3 + x[mm2] * SIGMA2 + x[mm1] * SIGMA, mm1, SIGMA);
	mem_fill(SIGMA5, SIGMA6, M24);

	set_array(M24 + x[0] * SIGMA5 + x[1] * SIGMA4 + x[mm2] * SIGMA3 + x[mm1] * SIGMA2, mm2, SIGMA2);
	mem_fill(SIGMA6, SIGMA7, M24);

	set_array(M24 + x[0] * SIGMA6 + x[1] * SIGMA5 + x[2] * SIGMA4 + x[mm1] * SIGMA3, mm3, SIGMA3);
	mem_fill(SIGMA7, SIGMA8, M24);

	for (int k = 0; k < mm3; k++)
		set_array(M24 + x[k] * SIGMA7 + x[k + 1] * SIGMA6 + x[k + 2] * SIGMA5 + x[k + 3] * SIGMA4, mm4 - k, SIGMA4);

	//Search
	pos = mm4, last_pos = n + mm4;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = *(M24 + *y_pos * SIGMA7 + *(y_pos + 1) * SIGMA6 + *(y_pos + 2) * SIGMA5 + *(y_pos + 3) * SIGMA4 +
			*(y_pos + m) * SIGMA3 + *(y_pos + mp1) * SIGMA2 + *(y_pos + mp2) * SIGMA + *(y_pos + mp3)))) {
			if (!memcmp(x, y + (pos - mm4), mm4)) {
				if (pos == last_pos)
					break;
				++count;
			}
			pos += D[*(y_pos + 3)];
		}
		else
			pos += step;
	}

	free(M24);
	return count;
}

// The MAW32 algorithm
int MAW32(unsigned char *x, int m, unsigned char *y, int n) {
	
	int mp1 = m + 1, mm1 = m - 1, mm2 = m - 2,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2,
		m3 = m * 3, m3m1 = 3 * m - 1, m3m2 = 3 * m - 2,
		step, pos, last_pos, count = 0;
	unsigned char *y_pos;
	int D[P_MAX];
	unsigned char* M32 = (unsigned char *)malloc(SIGMA6); // MAW32 search table

	// Preprocessing
	build_BMH_shift_table(D, x, m, mm1);

	// Fill the M32 search table
	memset(M32, m3, SIGMA);

	*(M32 + x[0]) = m3m1;
	mem_fill(SIGMA, SIGMA2, M32);

	for (int k = 0; k < mm1; k++)
		*(M32 + x[k] * SIGMA + x[k + 1]) = m3m2 - k;
	mem_fill(SIGMA2, SIGMA3, M32);

	set_array(M32 + x[0] * SIGMA2 + x[mm1] * SIGMA, m2m1, SIGMA);
	mem_fill(SIGMA3, SIGMA4, M32);

	for (int k = 0; k < mm1; k++)
		set_array(M32 + x[k] * SIGMA3 + x[k + 1] * SIGMA2, m2m2 - k, SIGMA2);
	mem_fill(SIGMA4, SIGMA5, M32);

	set_array(M32 + x[0] * SIGMA4 + x[mm1] * SIGMA3, mm1, SIGMA3);
	mem_fill(SIGMA5, SIGMA6, M32);

	for (int k = 0; k < mm1; k++)
		set_array(M32 + x[k] * SIGMA5 + x[k + 1] * SIGMA4, mm2 - k, SIGMA4);

	//Search
	pos = mm2, last_pos = n + mm2;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = *(M32 + *y_pos * SIGMA5 + *(y_pos + 1) * SIGMA4 + *(y_pos + m) * SIGMA3 + *(y_pos + mp1) * SIGMA2 + *(y_pos + m2) * SIGMA + *(y_pos + m2p1)))) {
			if (!memcmp(x, y + (pos - mm2), mm2)) {
				if (pos == last_pos)
					break;
				++count;
			}
			pos += D[*(y_pos + 1)];
		}
		else
			pos += step;
	}

	free(M32);
	return count;
}

// The MAW33 algorithm
int MAW33(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 3) return -1;
	
	int mp1 = m + 1, mp2 = m + 2, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2p2 = 2 * m + 2, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3,
		m3 = m * 3, m3m1 = 3 * m - 1, m3m2 = 3 * m - 2, m3m3 = 3 * m - 3,
		step, pos, last_pos, count = 0;
	unsigned char *y_pos;
	int D[P_MAX];
	unsigned char* M33 = (unsigned char *)malloc(SIGMA9); // MAW33 search table

	// Preprocessing
	build_BMH_shift_table(D, x, m, mm1);

	// Fill the M33 search table
	memset(M33, m3, SIGMA);

	*(M33 + x[0]) = m3m1;
	mem_fill(SIGMA, SIGMA2, M33);

	*(M33 + x[0] * SIGMA + x[1]) = m3m2;
	mem_fill(SIGMA2, SIGMA3, M33);

	for (int k = 0; k < mm2; k++)
		*(M33 + x[k] * SIGMA2 + x[k + 1] * SIGMA + x[k + 2]) = m3m3 - k;
	mem_fill(SIGMA3, SIGMA4, M33);

	set_array(M33 + x[0] * SIGMA3 + x[mm2] * SIGMA2 + x[mm1] * SIGMA, m2m1, SIGMA);
	mem_fill(SIGMA4, SIGMA5, M33);

	set_array(M33 + x[0] * SIGMA4 + x[1] * SIGMA3 + x[mm1] * SIGMA2, m2m2, SIGMA2);
	mem_fill(SIGMA5, SIGMA6, M33);

	for (int k = 0; k < mm2; k++)
		set_array(M33 + x[k] * SIGMA5 + x[k + 1] * SIGMA4 + x[k + 2] * SIGMA3, m2m3 - k, SIGMA3);
	mem_fill(SIGMA6, SIGMA7, M33);

	set_array(M33 + x[0] * SIGMA6 + x[mm2] * SIGMA5 + x[mm1] * SIGMA4, mm1, SIGMA4);
	mem_fill(SIGMA7, SIGMA8, M33);

	set_array(M33 + x[0] * SIGMA7 + x[1] * SIGMA6 + x[mm1] * SIGMA5, mm2, SIGMA5);
	mem_fill(SIGMA8, SIGMA9, M33);

	for (int k = 0; k < mm2; k++)
		set_array(M33 + x[k] * SIGMA8 + x[k + 1] * SIGMA7 + x[k + 2] * SIGMA6, mm3 - k, SIGMA6);

	//Search
	pos = mm3, last_pos = n + mm3;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = *(M33 + *y_pos * SIGMA8 + *(y_pos + 1) * SIGMA7 + *(y_pos + 2) * SIGMA6
			+ *(y_pos + m) * SIGMA5 + *(y_pos + mp1) * SIGMA4 + *(y_pos + mp2) * SIGMA3
			+ *(y_pos + m2) * SIGMA2 + *(y_pos + m2p1) * SIGMA + *(y_pos + m2p2)))) {
			if (!memcmp(x, y + (pos - mm3), mm3)) {
				if (pos == last_pos)
					break;
				++count;
			}
			pos += D[*(y_pos + 2)];
		}
		else
			pos += step;
	}

	free(M33);
	return count;
}

// The MAW42 algorithm
int MAW42(unsigned char *x, int m, unsigned char *y, int n) {

	int mp1 = m + 1, mm1 = m - 1, mm2 = m - 2,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2,
		m3 = 3 * m, m3p1 = 3 * m + 1, m3m1 = 3 * m - 1, m3m2 = 3 * m - 2,
		m4 = 4 * m, m4m1 = 4 * m - 1, m4m2 = 4 * m - 2,
		step, pos, last_pos, count = 0;
	unsigned char *y_pos;
	int D[P_MAX];
	unsigned char* M42 = (unsigned char *)malloc(SIGMA8); // MAW42 search table

	// Preprocessing
	build_BMH_shift_table(D, x, m, mm1);

	// Fill the M42 search table
	memset(M42, m4, SIGMA);

	*(M42 + x[0]) = m4m1;
	mem_fill(SIGMA, SIGMA2, M42);

	for (int k = 0; k < mm1; k++)
		*(M42 + x[k] * SIGMA + x[k + 1]) = m4m2 - k;
	mem_fill(SIGMA2, SIGMA3, M42);

	set_array(M42 + x[0] * SIGMA2 + x[mm1] * SIGMA, m3m1, SIGMA);
	mem_fill(SIGMA3, SIGMA4, M42);

	for (int k = 0; k < mm1; k++)
		set_array(M42 + x[k] * SIGMA3 + x[k + 1] * SIGMA2, m3m2 - k, SIGMA2);
	mem_fill(SIGMA4, SIGMA5, M42);

	set_array(M42 + x[0] * SIGMA4 + x[mm1] * SIGMA3, m2m1, SIGMA3);
	mem_fill(SIGMA5, SIGMA6, M42);

	for (int k = 0; k < mm1; k++)
		set_array(M42 + x[k] * SIGMA5 + x[k + 1] * SIGMA4, m2m2 - k, SIGMA4);
	mem_fill(SIGMA6, SIGMA7, M42);

	set_array(M42 + x[0] * SIGMA6 + x[mm1] * SIGMA5, mm1, SIGMA5);
	mem_fill(SIGMA7, SIGMA8, M42);

	for (int k = 0; k < mm1; k++)
		set_array(M42 + x[k] * SIGMA7 + x[k + 1] * SIGMA6, mm2 - k, SIGMA6);

	//Search
	pos = mm2, last_pos = n + mm2;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = *(M42 + *y_pos * SIGMA7 + *(y_pos + 1) * SIGMA6 + *(y_pos + m) * SIGMA5 + *(y_pos + mp1) * SIGMA4
			+ *(y_pos + m2) * SIGMA3 + *(y_pos + m2p1) * SIGMA2 + *(y_pos + m3) * SIGMA + *(y_pos + m3p1)))) {
			if (!memcmp(x, y + (pos - mm2), mm2)) {
				if (pos == last_pos)
					break;
				++count;
			}
			pos += D[*(y_pos + 1)];
		}
		else
			pos += step;
	}

	free(M42);
	return count;
}
