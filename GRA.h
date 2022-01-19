/*
 * This file is part of a graduate thesis to modernize the Galbraith-Rupray 
 * algorithm on equivalence classes instead of the exponential part of 
 * the Gaudry-Schost algorithm for counting points of a curve of genus 2.
 *
 * Developed for the Immanuel Kant Baltic Federal University.
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef GRA_H
#define GRA_H

#include "NTL/ZZ_p.h"
#include "NTL/ZZ_pX.h"
#include "vec_pair_long_long.h"
#include "ZZ_pJac2.h"
#include <assert.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include <time.h>
#include <math.h>
#include <cmath>
#include <thread>


NTL_CLIENT


class DivAndTrack 
{
  	public:
		ZZ_pJac2 D;
		ZZ alpha;
		ZZ beta;
		
		DivAndTrack() { }
		DivAndTrack(const DivAndTrack& dd)
		{
			D = dd.D;  alpha = dd.alpha;  beta = dd.beta;
		}
		
		~DivAndTrack() { }
		inline DivAndTrack& operator=(const DivAndTrack& dd)
		{
			this->D = dd.D;
			this->alpha = dd.alpha;
			this->beta = dd.beta;
			return *this;
		}
};

typedef struct
{
	double l1, l2;
	bool killThreads;
	int r, pD;
	int tot_cycles, sts_cycles;
	int distCnt;
	ZZ K, p;
	ZZ tot_jumps;
	ZZ B1min, B1max;
	ZZ B2min, B2max;
	ZZ_pJac2 KBaseDivisor;
	ZZ_pJac2 mBaseDivisor;
	ZZ_pJac2 pp1mBaseDivisor;
	ZZ_pJac2 BaseDivisor;
	DivAndTrack **O;
} RWData;

typedef Vec<DivAndTrack> vec_DivAndTrack;
typedef Vec<ZZ_pJac2> vec_ZZ_pJac2;

void randomize_f(ZZ_pX&);
void print_status_bar(int, unsigned int);
void clear_line();
void init_seed_time();

void equiv_class_repr(ZZ_pJac2&);
void equiv_class_repr(ZZ_pJac2&, ZZ&, ZZ&);

unsigned int hash_value(const ZZ_pJac2&);
unsigned int hash_value_2(const ZZ_pJac2&);

inline int index_from_k_kp_b_r(int, int, int, int);
inline int is_distinguished(const ZZ_pJac2&, const int);

void init_random_walk(RWData&, int);
void init_parameters(RWData&, const ZZ&, const ZZ&, const ZZ&);

void select_random_start_point(DivAndTrack&, const int, const RWData&);
void jump_one_step(DivAndTrack&, const DivAndTrack&, RWData&);

void run_thread(DivAndTrack&, int, RWData&);
ZZ run_test(RWData&, ZZ, ZZ);

#endif //GRA_H