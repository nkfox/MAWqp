#pragma comment(linker, "/STACK:4000000")
#include <malloc.h>
#include <time.h>
#include <windows.h>
#include <map>
#include <iostream>
using namespace std;

#define P_MAX 200	//maximum pattern length
#define V_MAX 20000 //maximum size of pointer and shift arrays

#define CURRENT_SIGMA 4
#define BIG_SIGMA 32
#define MIDDLE_SIGMA 16

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

#elif CURRENT_SIGMA == 32
#define SIGMA 32	//alphabet size
#define SIGMA2 1024
#define SIGMA3 32768
#define SIGMA4 1048576
#define LOG_SIGMA 5 //logarithm of alphabet size

#elif CURRENT_SIGMA == 64
#define SIGMA 64	//alphabet size
#define SIGMA2 4096
#define SIGMA3 262144
#define SIGMA4 16777216
#define LOG_SIGMA 6 //logarithm of alphabet size

#else // CURRENT_SIGMA == 128
#define SIGMA 128	//alphabet size
#define SIGMA2 16384
#define SIGMA3 2097152
#define SIGMA4 268435456
#define LOG_SIGMA 7 //logarithm of alphabet size

#endif
#if CURRENT_SIGMA > 16
#define SIGMA5 1
#define SIGMA6 1
#define SIGMA7 1
#define SIGMA8 1
#define SIGMA9 1
#endif

const int TOTAL = 1000000 + 5 * P_MAX;
unsigned char T[TOTAL], P[P_MAX];
int N = TOTAL - 5 * P_MAX, ITER = 200, m;

LARGE_INTEGER start_time, end_time, freq;
long long algo_time;

long long sum_maw22, sum_maw23, sum_maw24, sum_maw32, sum_maw33, sum_maw42;
long long sum_prep22, sum_prep23, sum_prep24, sum_prep32, sum_prep33, sum_prep42;
int maw22, maw23, maw24, maw32, maw33, maw42;
int maw22c, maw23c, maw24c, maw32c, maw33c, maw42c;

long long sum_maw22p, sum_maw23p, sum_maw24p, sum_maw32p, sum_maw33p, sum_maw42p;
long long sum_prep22p, sum_prep23p, sum_prep24p, sum_prep32p, sum_prep33p, sum_prep42p;
int maw22p, maw23p, maw24p, maw32p, maw33p, maw42p;
int maw22pc, maw23pc, maw24pc, maw32pc, maw33pc, maw42pc;

//-----------------Additional functions---------------------

#define mult_sigma(x) (x << LOG_SIGMA)
#define mult_sigma2(x) mult_sigma(mult_sigma(x))
#define mult_sigma3(x) mult_sigma2(mult_sigma(x))
#define mult_sigma4(x) mult_sigma3(mult_sigma(x))

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
	set_array(D, m, mult_sigma(int_size), int_size);
	for (int i = 0; i < end_index; i++)
		D[x[i]] = mm1 - i;
}

// Build the Berry-Ravindran shift table
template <class T>
void build_BR_shift_table(T *BR, const unsigned char *x, const int& m, const int& int_size, int end_index = -1) {
	int mm1 = m - 1, mm2 = m - 2;
	if (end_index == -1) end_index = mm1;
	set_array(BR[0], m, mult_sigma(int_size), int_size);
	BR[0][x[0]] = mm1;
	mem_fill(mult_sigma(int_size), mult_sigma2(int_size), (unsigned char*)BR);
	for (int i = 0; i < end_index; i++)
		BR[x[i]][x[i + 1]] = mm2 - i;
}

template <class T>
void build_modified_BR_shift_table(T *BR, const unsigned char *x, const int& m, const int& int_size) {
	int mm3 = m - 3, m2m3 = m * 2 - 3, m2m4 = m * 2 - 4;
	set_array(BR[0], m2m3, mult_sigma(int_size), int_size);
	for (int j = 0; j < mm3; j++)
		BR[0][x[j]] = m2m4 - j;
	mem_fill(mult_sigma(int_size), mult_sigma2(int_size), (unsigned char*)BR);
	for (int i = mm3; i < m; i++)
		BR[x[i - mm3]][x[i]] = m2m4 - i;
}

// Build tripple shift table
template <class T>
void build_tripple_shift_table(T *Tripple, const unsigned char *x, const int& m, const int& int_size) {
	int mm1 = m - 1, mm2 = m - 2, mm3 = m - 3;
	set_array(Tripple[0][0], m, mult_sigma(int_size), int_size);
	Tripple[0][0][x[0]] = mm1;
	mem_fill(mult_sigma(int_size), mult_sigma2(int_size), (unsigned char*)Tripple);
	Tripple[0][x[0]][x[1]] = mm2;
	mem_fill(mult_sigma2(int_size), mult_sigma3(int_size), (unsigned char*)Tripple);
	for (int i = 0; i < mm2; i++)
		Tripple[x[i]][x[i + 1]][x[i + 2]] = mm3 - i;
}

// Build quad shift table
template <class T>
void build_quad_shift_table(T *Quad, const unsigned char *x, const int& m, const int& int_size) {
	int mm1 = m - 1, mm2 = m - 2, mm3 = m - 3, mm4 = m - 4;
	set_array(Quad[0][0][0], m, mult_sigma(int_size), int_size);
	Quad[0][0][0][x[0]] = mm1;
	mem_fill(mult_sigma(int_size), mult_sigma2(int_size), (unsigned char*)Quad);
	Quad[0][0][x[0]][x[1]] = mm2;
	mem_fill(mult_sigma2(int_size), mult_sigma3(int_size), (unsigned char*)Quad);
	Quad[0][x[0]][x[1]][x[2]] = mm3;
	mem_fill(mult_sigma3(int_size), mult_sigma4(int_size), (unsigned char*)Quad);
	for (int i = 0; i < mm3; i++)
		Quad[x[i]][x[i + 1]][x[i + 2]][x[i + 3]] = mm4 - i;
}

template <class T, class U>
void fill_beginning(T from, U to, int amount, int length, int type_length) {
	for (int i = 0; i < amount; i++)
		set_array(from + mult_sigma(i), to + mult_sigma(i), length, type_length);
}

template <class T>
void fill_beginning_final(T from, int amount, int length, int type_length) {
	for (int i = 0; i < amount; i++)
		set_array(from + mult_sigma(i), i, length, type_length);
}

template <class T, class U>
void fill_first_letter(T from, U to, int block_length, int type_length, int letter, int length) {
	set_array(from, to, block_length, type_length);
	*(from + letter) = to - SIGMA;
	mem_fill(block_length, length, (unsigned char*)from);
}

template <class T, class U>
void fill_first_letter_final(T from, U value, int block_length, int type_length, int letter, int length) {
	set_array(from, value, block_length, type_length);
	*(from + letter) = value - 1;
	mem_fill(block_length, length, (unsigned char*)from);
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

template <class T>
void add_multsigma(T* destination, T source_to_add, int source_to_multsigma) {
	*destination = source_to_add + mult_sigma(source_to_multsigma);
}

template <class T>
void cyclesigma_add_multsigma(T* destination, T source_to_add, int* source_to_multsigma) {
	for (int j = 0; j < SIGMA; j++)
		add_multsigma(destination + j, source_to_add, *(source_to_multsigma + j));
}

template <class T, class U>
void add(T* destination, U* source_to_add, int value) {
	*destination = *source_to_add + value;
}

template <class T, class U>
void cyclesigma_add(T* destination, U source_to_add, int value) {
	for (int j = 0; j < SIGMA; j++)
		add(destination + j, source_to_add + j, value);
}

//-----------------------MAW--------------------------------

// The MAW22 algorithm
int MAW22(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep22 += algo_time;

	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw22 += algo_time;

	free(M22);
	return count;
}

// The MAW23 algorithm
int MAW23(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 3) return -1;
	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep23 += algo_time;

	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw23 += algo_time;

	free(M23);
	return count;
}

// The MAW24 algorithm
int MAW24(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 4) return -1;
	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep24 += algo_time;

	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw24 += algo_time;

	free(M24);
	return count;
}

// The MAW32 algorithm
int MAW32(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep32 += algo_time;

	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw32 += algo_time;

	free(M32);
	return count;
}

// The MAW33 algorithm
int MAW33(unsigned char *x, int m, unsigned char *y, int n) {
	if (m < 3) return -1;
	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep33 += algo_time;

	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw33 += algo_time;

	free(M33);
	return count;
}

// The MAW42 algorithm
int MAW42(unsigned char *x, int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep42 += algo_time;

	QueryPerformanceCounter(&start_time);

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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw42 += algo_time;

	free(M42);
	return count;
}

//-----------------------MAWP-------------------------------

// The MAW22 algorithm with pointers
int MAW22P(unsigned char *x, const int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start_time);

	int ***V0[SIGMA], **V1[SIGMA * P_MAX], *V2[SIGMA * (P_MAX + 1)], V3[SIGMA * P_MAX * 2]; //V3 - shift array; V0, V1, V2 - pointers arrays
	int D[P_MAX], Dm1[P_MAX], BRm1[SIGMA][SIGMA];
	unsigned char *y_pos;
	int pos, last_pos, step, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mm1 = m - 1, mm2 = m - 2,
		m2 = 2 * m, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2,
		m_sigma = mult_sigma(m), mm1_sigma = mult_sigma(mm1), m2m1_sigma = mult_sigma(m2m1),
		int_size_sigma = mult_sigma(int_size), int_size_sigma_2 = int_size_sigma * 2;

	//Preprocessing
	build_BMH_shift_table(D, x, m, mm1);
	build_BMH_shift_table(Dm1, x, mm1);
	build_BR_shift_table(BRm1, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	set_array(V0, V1 + mm1_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		add_multsigma(V0 + x[i], V1, Dm1[x[i]]);

	// Filling V1 with pointers to chunks of V2
	fill_first_letter(V1 + mm1_sigma, V2 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma);
	for (int i = 0; i < mm1; i++)
		cyclesigma_add_multsigma(V1 + mult_sigma(i), V2, *(BRm1 + x[mm2 - i]));

	// Filling V2 with pointers to chunks of V3
	fill_beginning(V2, V3, mm1, int_size_sigma, int_size);
	set_array(V2 + mm1_sigma, V3 + m2m1_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		add_multsigma(V2 + mm1_sigma + x[i], V3, Dm1[x[i]] + m);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V2 + mm1_sigma));
	add_multsigma(V2 + mm1_sigma + x[mm1], V3, mm1);

	// Filling V3 with shift values
	fill_beginning_final(V3, m, int_size_sigma, int_size);
	fill_first_letter_final(V3 + m2m1_sigma, m2, int_size_sigma, int_size, x[0], int_size_sigma);
	for (int i = m; i < m2m1; i++)
		cyclesigma_add(V3 + mult_sigma(i), *(BRm1 + x[m2m2 - i]), m);

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep22p += algo_time;

	QueryPerformanceCounter(&start_time);

	//Search
	pos = mm2, last_pos = n + mm2;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = V0[*y_pos][*(y_pos + 1)][*(y_pos + m)][*(y_pos + mp1)])) {
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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw22p += algo_time;

	return count;
}

// The MAW23 algorithm with pointers
int MAW23P(unsigned char *x, const int m, unsigned char *y, int n) {
	if (m < 3) return -1;

	QueryPerformanceCounter(&start_time);

	int *****V0[SIGMA], ****V1[SIGMA * (P_MAX - 1)], ***V2[SIGMA * P_MAX], **V3[SIGMA * (P_MAX + 1)],
		*V4[SIGMA * (2 * P_MAX - 1)], V5[SIGMA * P_MAX * 2]; //V5 - shift array; V0, V1, V2, V3, V4 - pointers arrays
	int D[P_MAX], Dm1[P_MAX], Dm2[P_MAX], BRm1[SIGMA][SIGMA], BRm2[SIGMA][SIGMA];
	unsigned char *y_pos;
	int pos, last_pos, step, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mp2 = m + 2, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3,
		m2 = 2 * m, m2m1 = m * 2 - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3,
		m_sigma = mult_sigma(m), mm1_sigma = mult_sigma(mm1), mm2_sigma = mult_sigma(mm2),
		m2m1_sigma = mult_sigma(m2m1), m2m2_sigma = mult_sigma(m2m2),
		int_size_sigma = mult_sigma(int_size), int_size_sigma_2 = int_size_sigma * 2, int_size_sigma_3 = int_size_sigma * 3;

#if SIGMA <= BIG_SIGMA
	int Tripple[SIGMA][SIGMA][SIGMA];
#else
	Map Tripple(m, x);
#endif

	//Preprocessing
	build_BMH_shift_table(D, x, m, mm1);
	build_BMH_shift_table(Dm1, x, mm1);
	build_BMH_shift_table(Dm2, x, mm2);

	build_BR_shift_table(BRm1, x, m, int_size, mm2);
	build_BR_shift_table(BRm2, x, mm1, int_size, mm2);

#if SIGMA <= BIG_SIGMA
	build_tripple_shift_table(Tripple, x, m, int_size);
#else
	Tripple.init3();
#endif

	// Filling V0 with pointers to chunks of V1
	set_array(V0, V1 + mm2_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		add_multsigma(V0 + x[i], V1, Dm2[x[i]]);

	// Filling V1 with pointers to chunks of V2
	for (int i = 0; i < mm2; i++)
		cyclesigma_add_multsigma(V1 + mult_sigma(i), V2, BRm2[x[mm3 - i]]);
	fill_first_letter(V1 + mm2_sigma, V2 + mm1_sigma, int_size_sigma, int_size, x[0], int_size_sigma);

	// Filling V2 with pointers to chunks of V3
	for (int i = 0; i < mm2; i++)
#if SIGMA <= BIG_SIGMA
		cyclesigma_add_multsigma(V2 + mult_sigma(i), V3, Tripple[x[mm3 - i]][x[mm2 - i]]);
#else
	for (int j = 0; j < SIGMA; j++)
		V2[mult_sigma(i) + j] = V3 + mult_sigma(Tripple.get3(x[mm3 - i], x[mm2 - i], j));
#endif
	fill_first_letter(V2 + mm2_sigma, V3 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_2);
	add_multsigma(V2 + mm2_sigma + x[1], V3, BRm1[x[0]][x[1]]);

	// Filling V3 with pointers to chunks of V4
	fill_beginning(V3, V4, mm2, int_size_sigma, int_size);
	set_array(V3 + mm2_sigma, V4 + m2m2_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		add_multsigma(V3 + mm2_sigma + x[i], V4, Dm2[x[i]] + m);
	mem_fill(int_size_sigma, int_size_sigma_3, (unsigned char*)(V3 + mm2_sigma));
	for (int i = 0; i < m; i++)
		add_multsigma(V3 + mm2_sigma + x[i], V4, m2m3 - i);
	add_multsigma(V3 + mm1_sigma + x[mm2], V4, mm1);

	// Filling V4 with pointers to chunks of V5
	fill_beginning(V4, V5, mm1, int_size_sigma, int_size);
	cyclesigma_add_multsigma(V4 + mm1_sigma, V5 + m_sigma, BRm2[x[mm2]]);
	add_multsigma(V4 + mm1_sigma + x[mm1], V5, mm1);
	for (int i = m; i < m2m2; i++)
		cyclesigma_add_multsigma(V4 + mult_sigma(i), V5 + m_sigma, BRm2[x[m2m3 - i]]);
	fill_first_letter(V4 + m2m2_sigma, V5 + m2m1_sigma, int_size_sigma, int_size, x[0], int_size_sigma);

	//Filling V5 with shift values
	fill_beginning_final(V5, m, int_size_sigma, int_size);
	for (int i = m; i < m2m2; i++)
#if SIGMA <= BIG_SIGMA
		cyclesigma_add(V5 + mult_sigma(i), Tripple[x[m2m3 - i]][x[m2m2 - i]], m);
#else
	for (int j = 0; j < SIGMA; j++)
		V5[mult_sigma(i) + j] = Tripple.get3(x[m2m3 - i], x[m2m2 - i], j) + m;
#endif
	fill_first_letter_final(V5 + m2m2_sigma, m2, int_size_sigma, int_size, x[0], int_size_sigma_2);
	add(V5 + m2m2_sigma + x[1], BRm1[x[0]] + x[1], m);

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep23p += algo_time;

	QueryPerformanceCounter(&start_time);

	//Search
	pos = mm3, last_pos = n + mm3;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = V0[*y_pos][*(y_pos + 1)][*(y_pos + 2)][*(y_pos + m)][*(y_pos + mp1)][*(y_pos + mp2)])) {
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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw23p += algo_time;

	return count;
}

// The MAW24 algorithm with pointers
int MAW24P(unsigned char *x, const int m, unsigned char *y, int n) {
	if (m < 4) return -1;

	QueryPerformanceCounter(&start_time);

	int *******V0[SIGMA], ******V1[SIGMA * (P_MAX - 2)], *****V2[SIGMA * (P_MAX - 1)], ****V3[SIGMA * P_MAX], ***V4[SIGMA * (P_MAX + 1)],
		**V5[SIGMA * (2 * P_MAX - 2)], *V6[SIGMA * (2 * P_MAX - 1)], V7[SIGMA * P_MAX * 2]; //V7 - shift array; V0, V1, V2, V3, V4, V5, V6 - pointers arrays
	int D[P_MAX], Dm1[P_MAX], Dm2[P_MAX], Dm3[P_MAX], BRm1[SIGMA][SIGMA], BRm2[SIGMA][SIGMA], BRm3[SIGMA][SIGMA], BRi[SIGMA][SIGMA];
	unsigned char *y_pos;
	int pos, last_pos, step, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mp2 = m + 2, mp3 = m + 3, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3, mm4 = m - 4,
		m2 = 2 * m, m2m1 = m * 2 - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3, m2m4 = 2 * m - 4,
		m_sigma = mult_sigma(m), mm1_sigma = mult_sigma(mm1), mm2_sigma = mult_sigma(mm2), mm3_sigma = mult_sigma(mm3),
		m2m1_sigma = mult_sigma(m2m1), m2m2_sigma = mult_sigma(m2m2), m2m3_sigma = mult_sigma(m2m3),
		int_size_sigma = mult_sigma(int_size), int_size_sigma_2 = int_size_sigma * 2;

#if SIGMA <= BIG_SIGMA
	int Tripple[SIGMA][SIGMA][SIGMA];
#else
	Map Tripple(m, x);
#endif
#if SIGMA <= MIDDLE_SIGMA
	int Quad[SIGMA][SIGMA][SIGMA][SIGMA];
#else
	Map Quad(m, x);
#endif

	//Preprocessing
	build_BMH_shift_table(D, x, m, mm1);
	build_BMH_shift_table(Dm1, x, mm1);
	build_BMH_shift_table(Dm2, x, mm2);
	build_BMH_shift_table(Dm3, x, mm3);

	build_BR_shift_table(BRm1, x, m, int_size);
	build_BR_shift_table(BRm2, x, mm1, int_size);
	build_BR_shift_table(BRm3, x, mm2, int_size);
	build_modified_BR_shift_table(BRi, x, m, int_size);

#if SIGMA <= BIG_SIGMA
	build_tripple_shift_table(Tripple, x, mm1, int_size);
#else
	Tripple.init3shift();
#endif
#if SIGMA <= MIDDLE_SIGMA
	build_quad_shift_table(Quad, x, m, int_size);
#else
	Quad.init4();
#endif

	// Filling V0 with pointers to chunks of V1
	set_array(V0, V1 + mm3_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm3; i++)
		add_multsigma(V0 + x[i], V1, Dm3[x[i]]);

	// Filling V1 with pointers to chunks of V2
	for (int i = 0; i < mm3; i++)
		cyclesigma_add_multsigma(V1 + mult_sigma(i), V2, BRm3[x[mm4 - i]]);
	fill_first_letter(V1 + mm3_sigma, V2 + mm2_sigma, int_size_sigma, int_size, x[0], int_size_sigma);

	// Filling V2 with pointers to chunks of V3
	for (int i = 0; i < mm3; i++)
#if SIGMA <= BIG_SIGMA
		cyclesigma_add_multsigma(V2 + mult_sigma(i), V3, Tripple[x[mm4 - i]][x[mm3 - i]]);
#else
	for (int j = 0; j < SIGMA; j++)
		V2[mult_sigma(i) + j] = V3 + mult_sigma(Tripple.get3shift(x[mm4 - i], x[mm3 - i], j));
#endif
	fill_first_letter(V2 + mm3_sigma, V3 + mm1_sigma, int_size_sigma, int_size, x[0], int_size_sigma_2);
	add_multsigma(V2 + mm3_sigma + x[1], V3, BRm2[x[0]][x[1]] >= mm3 ? BRm2[x[0]][x[1]] : mm3);

	// Filling V3 with pointers to chunks of V4
	for (int i = 0; i < mm3; i++)
#if SIGMA <= MIDDLE_SIGMA
		cyclesigma_add_multsigma(V3 + mult_sigma(i), V4, Quad[x[mm4 - i]][x[mm3 - i]][x[mm2 - i]]);
#else
	for (int j = 0; j < SIGMA; j++)
		V3[mult_sigma(i) + j] = V4 + mult_sigma(Quad.get4(x[mm4 - i], x[mm3 - i], x[mm2 - i], j));
#endif
#if SIGMA <= BIG_SIGMA
	cyclesigma_add_multsigma(V3 + mm3_sigma, V4 + SIGMA, Tripple[x[0]][x[1]]);
#else
	for (int j = 0; j < SIGMA; j++)
		V3[mm3_sigma + j] = V4 + SIGMA + mult_sigma(Tripple.get3shift(x[0], x[1], j));
#endif
	fill_first_letter(V3 + mm2_sigma, V4 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_2);
	add_multsigma(V3 + mm2_sigma + x[1], V4, BRm1[x[0]][x[1]] >= mm2 ? BRm1[x[0]][x[1]] : mm2);

	// Filling V4 with pointers to chunks of V5
	fill_beginning(V4, V5, mm3, int_size_sigma, int_size);
	for (int i = mm3; i < mm1; i++)
		cyclesigma_add_multsigma(V4 + mult_sigma(i), V5, BRi[x[mm1 - i]]);
	set_array(V4 + mm1_sigma, V5 + m2m3_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		add_multsigma(V4 + mm1_sigma + x[i], V5, Dm3[x[i]] + m);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V4 + mm1_sigma));
	add_multsigma(V4 + mm1_sigma + x[mm3], V5, mm1);

	// Filling V5 with pointers to chunks of V6
	fill_beginning(V5, V6, mm2, int_size_sigma, int_size);
	for (int i = mm2; i < m; i++)
		cyclesigma_add_multsigma(V5 + mult_sigma(i), V6 + mm2_sigma, BRm1[x[m2m4 - i]]);
	for (int i = m; i < m2m3; i++)
		cyclesigma_add_multsigma(V5 + mult_sigma(i), V6 + m_sigma, BRm3[x[m2m4 - i]]);
	fill_first_letter(V5 + m2m3_sigma, V6 + m2m2_sigma, int_size_sigma, int_size, x[0], int_size_sigma);

	// Filling V6 with pointers to chunks of V7
	fill_beginning(V6, V7, mm1, int_size_sigma, int_size);
	for (int i = mm1; i < m2m4; i++)
		cyclesigma_add_multsigma(V6 + mult_sigma(i), V7 + mm1_sigma, BRm1[x[m2m3 - i]]);
	for (int i = m; i < m2m3; i++)
#if SIGMA <= BIG_SIGMA
		cyclesigma_add_multsigma(V6 + mult_sigma(i), V7 + m_sigma, Tripple[x[m2m4 - i]][x[m2m3 - i]]);
#else
	for (int j = 0; j < SIGMA; j++)
		V6[mult_sigma(i) + j] = V7 + m_sigma + mult_sigma(Tripple.get3shift(x[m2m4 - i], x[m2m3 - i], j));
#endif
	fill_first_letter(V6 + m2m3_sigma, V7 + m2m1_sigma, int_size_sigma, int_size, x[0], int_size_sigma_2);
	add_multsigma(V6 + m2m3_sigma + x[1], V7 + m_sigma, (BRm2[x[0]][x[1]] >= mm3 ? BRm2[x[0]][x[1]] : mm3));

	//Filling V7 with shift values
	fill_beginning_final(V7, m, int_size_sigma, int_size);
	for (int i = m; i < m2m4; i++)
		cyclesigma_add(V7 + mult_sigma(i), BRm1[x[m2m2 - i]], m);
#if SIGMA <= BIG_SIGMA
	cyclesigma_add(V7 + m2m3_sigma, Tripple[x[0]][x[1]], mp1);
#else
	for (int j = 0; j < SIGMA; j++)
		V7[m2m3_sigma + j] = Tripple.get3shift(x[0], x[1], j) + mp1;
#endif
	for (int i = m; i < m2m3; i++)
#if SIGMA <= MIDDLE_SIGMA
		cyclesigma_add(V7 + mult_sigma(i), Quad[x[m2m4 - i]][x[m2m3 - i]][x[m2m2 - i]], m);
#else
	for (int j = 0; j < SIGMA; j++)
		V7[mult_sigma(i) + j] = Quad.get4(x[m2m4 - i], x[m2m3 - i], x[m2m2 - i], j) + m;
#endif
	fill_first_letter_final(V7 + m2m2_sigma, m2, int_size_sigma, int_size, x[0], int_size_sigma_2);
	add(V7 + m2m2_sigma + x[1], BRm1[x[0]][x[1]] >= mm2 ? BRm1[x[0]] + x[1] : &mm2, m);

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep24p += algo_time;

	QueryPerformanceCounter(&start_time);

	//Search
	pos = mm4, last_pos = n + mm4;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = V0[*y_pos][*(y_pos + 1)][*(y_pos + 2)][*(y_pos + 3)][*(y_pos + m)][*(y_pos + mp1)][*(y_pos + mp2)][*(y_pos + mp3)])) {
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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw24p += algo_time;

	return count;
}

// The MAW32 algorithm with pointers
int MAW32P(unsigned char *x, const int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start_time);

	int *****V0[SIGMA], ****V1[SIGMA * P_MAX], ***V2[SIGMA * (P_MAX + 1)], **V3[SIGMA * P_MAX * 2],
		*V4[SIGMA * (P_MAX * 2 + 1)], V5[SIGMA * P_MAX * 3]; //V5 - shift array; V0, V1, V2, V3, V4 - pointers arrays
	int D[P_MAX], Dm1[P_MAX], BRm1[SIGMA][SIGMA];
	unsigned char *y_pos;
	int pos, last_pos, step, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mm1 = m - 1, mm2 = m - 2,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2,
		m3 = m * 3, m3m1 = m * 3 - 1, m3m2 = m * 3 - 2,
		m_sigma = mult_sigma(m), mm1_sigma = mult_sigma(mm1), m2_sigma = mult_sigma(m2),
		m2m1_sigma = mult_sigma(m2m1), m3m1_sigma = mult_sigma(m3m1),
		int_size_sigma = int_size * SIGMA, int_size_sigma_2 = int_size_sigma * 2;

	//Preprocessing
	build_BMH_shift_table(D, x, m, mm1);
	build_BMH_shift_table(Dm1, x, mm1);
	build_BR_shift_table(BRm1, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	set_array(V0, V1 + mm1_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		add_multsigma(V0 + x[i], V1, Dm1[x[i]]);

	// Filling V1 with pointers to chunks of V2
	fill_first_letter(V1 + mm1_sigma, V2 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma);
	for (int i = 0; i < mm1; i++)
		cyclesigma_add_multsigma(V1 + mult_sigma(i), V2, *(BRm1 + x[mm2 - i]));

	// Filling V2 with pointers to chunks of V3
	fill_beginning(V2, V3, mm1, int_size_sigma, int_size);
	set_array(V2 + mm1_sigma, V3 + m2m1_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		add_multsigma(V2 + mm1_sigma + x[i], V3, Dm1[x[i]] + m);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V2 + mm1_sigma));
	add_multsigma(V2 + mm1_sigma + x[mm1], V3, mm1);

	// Filling V3 with pointers to chunks of V4
	fill_beginning(V3, V4, m, int_size_sigma, int_size);
	fill_first_letter(V3 + m2m1_sigma, V4 + m2_sigma, int_size_sigma, int_size, x[0], int_size_sigma);
	for (int i = m; i < m2m1; i++)
		cyclesigma_add_multsigma(V3 + mult_sigma(i), V4 + m_sigma, BRm1[x[m2m2 - i]]);

	// Filling V4 with pointers to chunks of V5
	fill_beginning(V4, V5, m2m1, int_size_sigma, int_size);
	set_array(V4 + m2m1_sigma, V5 + m3m1_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		add_multsigma(V4 + m2m1_sigma + x[i], V5 + m2_sigma, Dm1[x[i]]);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V4 + m2m1_sigma));
	add_multsigma(V4 + m2m1_sigma + x[mm1], V5, m2m1);

	// Filling V5 with shift values
	fill_beginning_final(V5, m2, int_size_sigma, int_size);
	fill_first_letter_final(V5 + m3m1_sigma, m3, int_size_sigma, int_size, x[0], int_size_sigma);
	for (int i = m2; i < m3m1; i++)
		cyclesigma_add(V5 + mult_sigma(i), BRm1[x[m3m2 - i]], m2);

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep32p += algo_time;

	QueryPerformanceCounter(&start_time);

	//Search
	pos = mm2, last_pos = n + mm2;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = V0[*y_pos][*(y_pos + 1)][*(y_pos + m)][*(y_pos + mp1)][*(y_pos + m2)][*(y_pos + m2p1)])) {
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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw32p += algo_time;

	return count;
}

// The MAW33 algorithm with pointers
int MAW33P(unsigned char *x, const int m, unsigned char *y, int n) {
	if (m < 3) return -1;

	QueryPerformanceCounter(&start_time);

	int ********V0[SIGMA], *******V1[SIGMA * (P_MAX - 1)], ******V2[SIGMA * P_MAX],
		*****V3[SIGMA * (P_MAX + 1)], ****V4[SIGMA * (2 * P_MAX - 1)], ***V5[SIGMA * P_MAX * 2],
		**V6[SIGMA * (P_MAX * 2 + 1)], *V7[SIGMA * (3 * P_MAX - 1)], V8[SIGMA * P_MAX * 3]; //V7 - shift array; V0, V1, V2, V3, V4, V5, V6 - pointers arrays
	int D[P_MAX], Dm1[P_MAX], Dm2[P_MAX], BRm1[SIGMA][SIGMA], BRm2[SIGMA][SIGMA];
	unsigned char *y_pos;
	int pos, last_pos, step, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mp2 = m + 2, mm1 = m - 1, mm2 = m - 2, mm3 = m - 3,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2p2 = 2 * m + 2, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2, m2m3 = 2 * m - 3,
		m3 = m * 3, m3m1 = m * 3 - 1, m3m2 = m * 3 - 2, m3m3 = m * 3 - 3,
		m_sigma = mult_sigma(m), mm1_sigma = mult_sigma(mm1), mm2_sigma = mult_sigma(mm2),
		m2_sigma = mult_sigma(m2), m2m1_sigma = mult_sigma(m2m1), m2m2_sigma = mult_sigma(m2m2),
		m3m1_sigma = mult_sigma(m3m1), m3m2_sigma = mult_sigma(m3m2),
		int_size_sigma = mult_sigma(int_size), int_size_sigma_2 = int_size_sigma * 2, int_size_sigma_3 = int_size_sigma * 3;

#if SIGMA <= BIG_SIGMA
	int Tripple[SIGMA][SIGMA][SIGMA];
#else
	Map Tripple(m, x);
#endif

	//Preprocessing
	build_BMH_shift_table(D, x, m, mm1);
	build_BMH_shift_table(Dm1, x, mm1);
	build_BMH_shift_table(Dm2, x, mm2);

	build_BR_shift_table(BRm1, x, m, int_size, mm2);
	build_BR_shift_table(BRm2, x, mm1, int_size, mm2);

#if SIGMA <= BIG_SIGMA
	build_tripple_shift_table(Tripple, x, m, int_size);
#else
	Tripple.init3();
#endif

	// Filling V0 with pointers to chunks of V1
	set_array(V0, V1 + mm2_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		add_multsigma(V0 + x[i], V1, Dm2[x[i]]);

	// Filling V1 with pointers to chunks of V2
	for (int i = 0; i < mm2; i++)
		cyclesigma_add_multsigma(V1 + mult_sigma(i), V2, BRm2[x[mm3 - i]]);
	fill_first_letter(V1 + mm2_sigma, V2 + mm1_sigma, int_size_sigma, int_size, x[0], int_size_sigma);

	// Filling V2 with pointers to chunks of V3
	for (int i = 0; i < mm2; i++)
#if SIGMA <= BIG_SIGMA
		cyclesigma_add_multsigma(V2 + mult_sigma(i), V3, Tripple[x[mm3 - i]][x[mm2 - i]]);
#else
	for (int j = 0; j < SIGMA; j++)
		V2[mult_sigma(i) + j] = V3 + mult_sigma(Tripple.get3(x[mm3 - i], x[mm2 - i], j));
#endif
	fill_first_letter(V2 + mm2_sigma, V3 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma_2);
	add_multsigma(V2 + mm2_sigma + x[1], V3, BRm1[x[0]][x[1]]);

	// Filling V3 with pointers to chunks of V4
	fill_beginning(V3, V4, mm2, int_size_sigma, int_size);
	set_array(V3 + mm2_sigma, V4 + m2m2_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		add_multsigma(V3 + mm2_sigma + x[i], V4, Dm2[x[i]] + m);
	mem_fill(int_size_sigma, int_size_sigma_3, (unsigned char*)(V3 + mm2_sigma));
	for (int i = 0; i < m; i++)
		add_multsigma(V3 + mm2_sigma + x[i], V4, m2m3 - i);
	add_multsigma(V3 + mm1_sigma + x[mm2], V4, mm1);

	// Filling V4 with pointers to chunks of V5
	fill_beginning(V4, V5, mm1, int_size_sigma, int_size);
	cyclesigma_add_multsigma(V4 + mm1_sigma, V5 + m_sigma, BRm2[x[mm2]]);
	add_multsigma(V4 + mm1_sigma + x[mm1], V5, mm1);
	for (int i = m; i < m2m2; i++)
		cyclesigma_add_multsigma(V4 + mult_sigma(i), V5 + m_sigma, BRm2[x[m2m3 - i]]);
	fill_first_letter(V4 + m2m2_sigma, V5 + m2m1_sigma, int_size_sigma, int_size, x[0], int_size_sigma);

	// Filling V5 with pointers to chunks of V6
	fill_beginning(V5, V6, m, int_size_sigma, int_size);
	for (int i = m; i < m2m2; i++)
#if SIGMA <= BIG_SIGMA
		cyclesigma_add_multsigma(V5 + mult_sigma(i), V6 + m_sigma, Tripple[x[m2m3 - i]][x[m2m2 - i]]);
#else
	for (int j = 0; j < SIGMA; j++)
		V5[mult_sigma(i) + j] = V6 + m_sigma + mult_sigma(Tripple.get3(x[m2m3 - i], x[m2m2 - i], j));
#endif
	fill_first_letter(V5 + m2m2_sigma, V6 + m2_sigma, int_size_sigma, int_size, x[0], int_size_sigma_2);
	add_multsigma(V5 + m2m2_sigma + x[1], V6 + m_sigma, BRm1[x[0]][x[1]]);

	// Filling V6 with pointers to chunks of V7
	fill_beginning(V6, V7, m2m2, int_size_sigma, int_size);
	set_array(V6 + m2m2_sigma, V7 + m3m2_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm2; i++)
		add_multsigma(V6 + m2m2_sigma + x[i], V7 + m2_sigma, Dm2[x[i]]);
	mem_fill(int_size_sigma, int_size_sigma_3, (unsigned char*)(V6 + m2m2_sigma));
	add_multsigma(V6 + m2m1_sigma + x[mm2], V7, m2m1);
	for (int i = 0; i < m; i++)
		add_multsigma(V6 + m2m2_sigma + x[i], V7, m3m3 - i);

	// Filling V7 with pointers to chunks of V8
	fill_beginning(V7, V8, m2m1, int_size_sigma, int_size);
	cyclesigma_add_multsigma(V7 + m2m1_sigma, V8 + m2_sigma, BRm2[x[mm2]]);
	add_multsigma(V7 + m2m1_sigma + x[mm1], V8, m2m1);
	for (int i = m2; i < m3m2; i++)
		cyclesigma_add_multsigma(V7 + mult_sigma(i), V8 + m2_sigma, BRm2[x[m3m3 - i]]);
	fill_first_letter(V7 + m3m2_sigma, V8 + m3m1_sigma, int_size_sigma, int_size, x[0], int_size_sigma);

	//Filling V8 with shift values
	fill_beginning_final(V8, m2, int_size_sigma, int_size);
	for (int i = m2; i < m3m2; i++)
#if SIGMA <= BIG_SIGMA
		cyclesigma_add(V8 + mult_sigma(i), Tripple[x[m3m3 - i]][x[m3m2 - i]], m2);
#else
	for (int j = 0; j < SIGMA; j++)
		V8[mult_sigma(i) + j] = Tripple.get3(x[m3m3 - i], x[m3m2 - i], j) + m2;
#endif
	fill_first_letter_final(V8 + m3m2_sigma, m3, int_size_sigma, int_size, x[0], int_size_sigma_2);
	add(V8 + m3m2_sigma + x[1], BRm1[x[0]] + x[1], m2);

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep33p += algo_time;

	QueryPerformanceCounter(&start_time);

	//Search
	pos = mm3, last_pos = n + mm3;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = V0[*y_pos][*(y_pos + 1)][*(y_pos + 2)][*(y_pos + m)][*(y_pos + mp1)][*(y_pos + mp2)][*(y_pos + m2)][*(y_pos + m2p1)][*(y_pos + m2p2)])) {
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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw33p += algo_time;

	return count;
}

// The MAW42 algorithm with pointers
int MAW42P(unsigned char *x, const int m, unsigned char *y, int n) {
	QueryPerformanceCounter(&start_time);

	int *******V0[SIGMA], ******V1[SIGMA * P_MAX], *****V2[SIGMA * (P_MAX + 1)], ****V3[SIGMA * P_MAX * 2], ***V4[SIGMA * (P_MAX * 2 + 1)],
		**V5[SIGMA * P_MAX * 3], *V6[SIGMA * (P_MAX * 3 + 1)], V7[SIGMA * P_MAX * 4]; //V7 - shift array; V0, V1, V2, V3, V4, V5, V6 - pointers arrays
	int D[P_MAX], Dm1[P_MAX], BRm1[SIGMA][SIGMA];
	unsigned char *y_pos;
	int pos, last_pos, step, count = 0, int_size = sizeof(int),
		mp1 = m + 1, mm1 = m - 1, mm2 = m - 2,
		m2 = 2 * m, m2p1 = 2 * m + 1, m2m1 = 2 * m - 1, m2m2 = 2 * m - 2,
		m3 = m * 3, m3p1 = m * 3 + 1, m3m1 = m * 3 - 1, m3m2 = m * 3 - 2,
		m4 = m * 4, m4m1 = m * 4 - 1, m4m2 = m * 4 - 2,
		m_sigma = mult_sigma(m), mm1_sigma = mult_sigma(mm1), m2_sigma = mult_sigma(m2), m2m1_sigma = mult_sigma(m2m1),
		m3_sigma = mult_sigma(m3), m3m1_sigma = mult_sigma(m3m1), m4m1_sigma = mult_sigma(m4m1),
		int_size_sigma = int_size * SIGMA, int_size_sigma_2 = int_size_sigma * 2;

	//Preprocessing
	build_BMH_shift_table(D, x, m, mm1);
	build_BMH_shift_table(Dm1, x, mm1);
	build_BR_shift_table(BRm1, x, m, int_size);

	// Filling V0 with pointers to chunks of V1
	set_array(V0, V1 + mm1_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		add_multsigma(V0 + x[i], V1, Dm1[x[i]]);

	// Filling V1 with pointers to chunks of V2
	fill_first_letter(V1 + mm1_sigma, V2 + m_sigma, int_size_sigma, int_size, x[0], int_size_sigma);
	for (int i = 0; i < mm1; i++)
		cyclesigma_add_multsigma(V1 + mult_sigma(i), V2, *(BRm1 + x[mm2 - i]));

	// Filling V2 with pointers to chunks of V3
	fill_beginning(V2, V3, mm1, int_size_sigma, int_size);
	set_array(V2 + mm1_sigma, V3 + m2m1_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		add_multsigma(V2 + mm1_sigma + x[i], V3, Dm1[x[i]] + m);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V2 + mm1_sigma));
	add_multsigma(V2 + mm1_sigma + x[mm1], V3, mm1);

	// Filling V3 with pointers to chunks of V4
	fill_beginning(V3, V4, m, int_size_sigma, int_size);
	fill_first_letter(V3 + m2m1_sigma, V4 + m2_sigma, int_size_sigma, int_size, x[0], int_size_sigma);
	for (int i = m; i < m2m1; i++)
		cyclesigma_add_multsigma(V3 + mult_sigma(i), V4 + m_sigma, BRm1[x[m2m2 - i]]);

	// Filling V4 with pointers to chunks of V5
	fill_beginning(V4, V5, m2m1, int_size_sigma, int_size);
	set_array(V4 + m2m1_sigma, V5 + m3m1_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		add_multsigma(V4 + m2m1_sigma + x[i], V5 + m2_sigma, Dm1[x[i]]);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V4 + m2m1_sigma));
	add_multsigma(V4 + m2m1_sigma + x[mm1], V5, m2m1);

	// Filling V5 with pointers to chunks of V6
	fill_beginning(V5, V6, m2, int_size_sigma, int_size);
	fill_first_letter(V5 + m3m1_sigma, V6 + m3_sigma, int_size_sigma, int_size, x[0], int_size_sigma);
	for (int i = m2; i < m3m1; i++)
		cyclesigma_add_multsigma(V5 + mult_sigma(i), V6 + m2_sigma, BRm1[x[m3m2 - i]]);

	// Filling V6 with pointers to chunks of V7
	fill_beginning(V6, V7, m3m1, int_size_sigma, int_size);
	set_array(V6 + m3m1_sigma, V7 + m4m1_sigma, int_size_sigma, int_size);
	for (int i = 0; i < mm1; i++)
		add_multsigma(V6 + m3m1_sigma + x[i], V7 + m3_sigma, Dm1[x[i]]);
	mem_fill(int_size_sigma, int_size_sigma_2, (unsigned char*)(V6 + m3m1_sigma));
	add_multsigma(V6 + m3m1_sigma + x[mm1], V7, m3m1);

	// Filling V7 with shift values
	fill_beginning_final(V7, m3, int_size_sigma, int_size);
	fill_first_letter_final(V7 + m4m1_sigma, m4, int_size_sigma, int_size, x[0], int_size_sigma); //!!! to the end of the block
	for (int i = m3; i < m4m1; i++)
		cyclesigma_add(V7 + mult_sigma(i), BRm1[x[m4m2 - i]], m3);

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_prep42p += algo_time;

	QueryPerformanceCounter(&start_time);

	//Search
	pos = mm2, last_pos = n + mm2;
	memcpy(y + n, x, m); //append the text with a stop pattern
	while (true) {
		y_pos = y + pos;
		if (!(step = V0[*y_pos][*(y_pos + 1)][*(y_pos + m)][*(y_pos + mp1)][*(y_pos + m2)][*(y_pos + m2p1)][*(y_pos + m3)][*(y_pos + m3p1)])) {
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

	QueryPerformanceCounter(&end_time);
	algo_time = (end_time.QuadPart - start_time.QuadPart) * 1000000 / freq.QuadPart;
	sum_maw42p += algo_time;

	return count;
}

//----------------------Testing-----------------------------

void generateRandom() {

	static int glob = 0;
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
	fread(T, 1, N, f);
	for (int i = 0; i < N; i++)
		switch (T[i]) {
		case 'a': T[i] = 0; break;
		case 'c': T[i] = 1; break;
		case 't': T[i] = 2; break;
		case 'g': T[i] = 3;
	}
}

void testPerformance() {
	QueryPerformanceFrequency(&freq);

	std::string filename = "output";
	filename += char(SIGMA % 10 + '0');
	if (SIGMA >= 10)
		filename += char(SIGMA % 10 + '0');
	if (SIGMA >= 100)
		filename += char((SIGMA / 10) % 10 + '0');
	filename += ".csv";
	FILE *f = fopen(filename.c_str(), "wt");
	generateRandom();
	//DNA();

	fprintf(f, "b=%d N=%d ITER=%d\n", SIGMA, N, ITER);

	fprintf(f, "m,MAW22,MAW23,MAW24,MAW32,MAW33,MAW42,  MAW22,MAW23P,MAW24P,MAW32P,MAW33P,MAW42P,  ,\
			   	PREP22,PREP23,PREP24,PREP32,PREP33,PREP42,  PREP22P,PREP23P,PREP24P,PREP32P,PREP33P,PREP42P,  ,\
				SUM22,SUM23,SUM24,SUM32,SUM33,SUM42,  SUM22P,SUM23P,SUM24P,SUM32P,SUM33P,SUM42P");
	for (m = 2; m < 81; m < 10 ? m++ : m += 10) {
		sum_maw22 = sum_maw23 = sum_maw24 = sum_maw32 = sum_maw33 = sum_maw42 = 0;
		sum_prep22 = sum_prep23 = sum_prep24 = sum_prep32 = sum_prep33 = sum_prep42 = 0;

		sum_maw22p = sum_maw23p = sum_maw24p = sum_maw32p = sum_maw33p = sum_maw42p = 0;
		sum_prep22p = sum_prep23p = sum_prep24p = sum_prep32p = sum_prep33p = sum_prep42p = 0;

		int nm = N - m;
		for (int ii = 0; ii < ITER; ii++) {
			srand((unsigned)time(NULL));
			int patpos = rand() % (N - m - 2);
			for (int i = 0; i < m; i++)
				P[i] = T[patpos + i];
			for (int i = 0; i < m; i++)
				T[N - m + i] = P[i];

			maw22 = MAW22(P, m, T, nm);
			if (SIGMA < 16) maw23 = MAW23(P, m, T, nm);
			if (SIGMA < 16) maw24 = MAW24(P, m, T, nm);
			if (SIGMA < 16) maw32 = MAW32(P, m, T, nm);
			if (SIGMA < 16) maw33 = MAW33(P, m, T, nm);
			if (SIGMA < 16) maw42 = MAW42(P, m, T, nm);

			maw22p = MAW22P(P, m, T, nm);
			maw23p = MAW23P(P, m, T, nm);
			maw24p = MAW24P(P, m, T, nm);
			maw32p = MAW32P(P, m, T, nm);
			maw33p = MAW33P(P, m, T, nm);
			maw42p = MAW42P(P, m, T, nm);
		}
		printf("b=%d m=%d\n", SIGMA, m);
		printf("%d %d %d %d %d %d   %d %d %d %d %d %d\n\n", maw22, maw23, maw24, maw32, maw33, maw42, maw22p, maw23p, maw24p, maw32p, maw33p, maw42p);
		printf("%7.lld %7.lld %7.lld %7.lld %7.lld %7.lld   %7.lld %7.lld %7.lld %7.lld %7.lld %7.lld\n\n",
			sum_prep22, sum_prep23, sum_prep24, sum_prep32, sum_prep33, sum_prep42, sum_prep22p, sum_prep23p, sum_prep24p, sum_prep32p, sum_prep33p, sum_prep42p);
		fprintf(f, "\n%2.d,  %7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,\
				   			,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,\
							,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld",
			m, sum_maw22, sum_maw23, sum_maw24, sum_maw32, sum_maw33, sum_maw42,
			sum_maw22p, sum_maw23p, sum_maw24p, sum_maw32p, sum_maw33p, sum_maw42p,
			sum_prep22, sum_prep23, sum_prep24, sum_prep32, sum_prep33, sum_prep42,
			sum_prep22p, sum_prep23p, sum_prep24p, sum_prep32p, sum_prep33p, sum_prep42p,
			sum_maw22 + sum_prep22, sum_maw23 + sum_prep23, sum_maw24 + sum_prep24, sum_maw32 + sum_prep32, sum_maw33 + sum_prep33, sum_maw42 + sum_prep42,
			sum_maw22p + sum_prep22p, sum_maw23p + sum_prep23p, sum_maw24p + sum_prep24p, sum_maw32p + sum_prep32p, sum_maw33p + sum_prep33p, sum_maw42p + sum_prep42p);
	}
	fclose(f);
	system("pause");
}

void testMAW42() {
	QueryPerformanceFrequency(&freq);

	FILE *f = fopen("output.csv", "wt");
	generateRandom();
	//DNA();

	fprintf(f, "b=%d N=%d ITER=%d\n", SIGMA, N, ITER);

	fprintf(f, "m,MAW1,MAW2,,,PREP1,PREP2,,,SUM1,SUM2");

	for (m = 4; m < 81; m < 10 ? m++ : m += 10)
	{
		int sum_maw1 = 0, sum_maw2 = 0;
		int sum_prep1 = 0, sum_prep2 = 0;
		int maw1, maw2;

		int nm = N - m;
		for (int ii = 0; ii < ITER; ii++)
		{
			srand((unsigned)time(NULL));
			int patpos = rand() % (N - m - 2);
			for (int i = 0; i < m; i++)
				P[i] = T[patpos + i];

			for (int i = 0; i < m; i++)
				T[N - m + i] = P[i];

			maw1 = MAW42P(P, m, T, nm);
			maw2 = MAW42P(P, m, T, nm);
		}
		printf("b=%d m=%d\n", SIGMA, m);
		printf("%d %d\n\n", maw1, maw2);
		printf("%7.lld %7.lld\n\n", sum_prep1, sum_prep2);
		fprintf(f, "\n%2.d,%7.lld,%7.lld,%7.lld,,%7.lld,%7.lld,%7.lld,,%7.lld,%7.lld",
			m, sum_maw1, sum_maw2, sum_maw2 - sum_maw1, sum_prep1, sum_prep2, sum_prep2 - sum_prep1, sum_maw1 + sum_prep1, sum_maw2 + sum_prep2);

	}
	fclose(f);
	system("pause");
}

int testMAW()
{
	QueryPerformanceFrequency(&freq);

	for (int num = 1000; num < 4000; num++)
	{
		if (num % 1000 == 0)
			cout << num << endl;

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
			num += 6 * (int)pow(10, m) - 1;
			continue;
		}

		for (int i = 0; i < m; i++)
			T[N - m + i] = P[i];

		maw22 = MAW22(P, m, T, N - m);
		maw23 = MAW23(P, m, T, N - m);
		maw24 = MAW24(P, m, T, N - m);
		maw32 = MAW32(P, m, T, N - m);
		maw33 = MAW33(P, m, T, N - m);
		maw42 = MAW42(P, m, T, N - m);

		maw22p = MAW22P(P, m, T, N - m);
		maw23p = MAW23P(P, m, T, N - m);
		maw24p = MAW24P(P, m, T, N - m);
		maw32p = MAW32P(P, m, T, N - m);
		maw33p = MAW33P(P, m, T, N - m);
		maw42p = MAW42P(P, m, T, N - m);

		if (maw22 != maw22p || maw22 != maw23 || maw22 != maw24 || maw22 != maw32 || maw22 != maw33 || maw22 != maw42
			|| maw22 != maw23p || maw22 != maw24p || maw22 != maw32p || maw22 != maw33p || maw22 != maw42p)
		{
			for (int i = 0; i < m; i++)
				cout << char(P[i] + '0');
			cout << " " << m << endl;

			printf("%d %d %d %d %d %d   %d %d %d %d %d %d\n\n", maw22, maw23, maw24, maw32, maw33, maw42, maw22p, maw23p, maw24p, maw32p, maw33p, maw42p);
		}
	}

	system("pause");
	return 0;
}

int simplestTest()
{
	QueryPerformanceFrequency(&freq);

	m = 4;

	//for (m = 4; m < 9; m++)
	//for (int k = 0; k < 10; k++)
	{
		generateRandom();

		srand((unsigned)time(NULL));
		int patpos = rand() % (N - m - 2);
		for (int i = 0; i < m; i++)
			P[i] = T[patpos + i];

		for (int i = 0; i < m; i++)
			T[N - m + i] = P[i];

		//P[0] = 0; P[1] = 1; P[2] = 0; P[3] = 3;
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
		maw23 = MAW23(P, m, T, N - m);
		maw24 = MAW24(P, m, T, N - m);
		maw32 = MAW32(P, m, T, N - m);
		maw33 = MAW33(P, m, T, N - m);
		maw42 = MAW42(P, m, T, N - m);

		maw22p = MAW22P(P, m, T, N - m);
		maw23p = MAW23P(P, m, T, N - m);
		maw24p = MAW24P(P, m, T, N - m);
		maw32p = MAW32P(P, m, T, N - m);
		maw33p = MAW33P(P, m, T, N - m);
		maw42p = MAW42P(P, m, T, N - m);

		for (int i = 0; i < m; i++)
			cout << char(P[i] + '0');
		cout << " " << m << endl;

		printf("count	%d %d %d %d %d %d   %d %d %d %d %d %d\n\n", maw22, maw23, maw24, maw32, maw33, maw42, maw22p, maw23p, maw24p, maw32p, maw33p, maw42p);
		printf("prep	%7.lld %7.lld %7.lld %7.lld %7.lld %7.lld   %7.lld %7.lld %7.lld %7.lld %7.lld %7.lld\n\n",
			sum_prep22, sum_prep23, sum_prep24, sum_prep32, sum_prep33, sum_prep42, sum_prep22p, sum_prep23p, sum_prep24p, sum_prep32p, sum_prep33p, sum_prep42p);
		printf("algo	%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld\n",
			sum_maw22, sum_maw23, sum_maw24, sum_maw32, sum_maw33, sum_maw42,
			sum_maw22p, sum_maw23p, sum_maw24p, sum_maw32p, sum_maw33p, sum_maw42p);
		printf("sum		%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld,%7.lld",
			sum_maw22 + sum_prep22, sum_maw23 + sum_prep23, sum_maw24 + sum_prep24, sum_maw32 + sum_prep32, sum_maw33 + sum_prep33, sum_maw42 + sum_prep42,
			sum_maw22p + sum_prep22p, sum_maw23p + sum_prep23p, sum_maw24p + sum_prep24p, sum_maw32p + sum_prep32p, sum_maw33p + sum_prep33p, sum_maw42p + sum_prep42p);
	}

	system("pause");
	return 0;
}

int main()
{
	testPerformance();
	//testMAW();
	//simplestTest();
	//testMAW42();
	return 0;
}