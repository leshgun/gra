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

#include <GRA.h>

// the main parametrs that need to be selected for the considered curve
// selected as ubiversal for most curves
#define HASH_SIZE 64
#define QUEUE_SIZE 64

// default value by Galbraith-Ruprai
// can be changed if it speeds up the algorithm
#define PROBA_DIST_INC -1
#define WILD_BOUND_SCALE 1/2
#define SEARCH_SPACE_SCALE 1

// choose a type of the Galbraith-Ruprai algorithm
// which solve "side-to-side" cycle by:
// 1 : standart (restart by a new point)
// 2 : entry point (return the point of entry into the loop)
// 3 : jump out (deterministic loop jump)
// 4 : combination of 2 and 3, where walk jumps out from the loop until
//		the steps less than some value
// 5 : the same as 4, only using a hash image of the loop
#define ALGORITHM_TYPE 1

// jump out until the steps of a walk less ("search space" * JUMP_OUT_SPACE)
// 0 <= JUMP_OUT_SPACE <= 1
#define JUMP_OUT_SPACE 1/2

#define OUTPUT_SEP "========================================="


// get random curve equation: y^2 = f(x)
// f(x) = x^5 + c4*x^4 + c3*x^3 + c2*x^2 + c1*x + c0
void randomize_f(ZZ_pX& f)
{
    while (true) 
	{
        ZZ_pX h;
        random(h, 6);
        if (deg(h) != 5) continue;
        MakeMonic(h);

        ZZ_p d = resultant(h, diff(h));
        if (d != 0)
		{
            SetCoeff(f, 5, coeff(h, 5));
            SetCoeff(f, 4, coeff(h, 4));
            SetCoeff(f, 3, coeff(h, 3));
            SetCoeff(f, 2, coeff(h, 2));
            SetCoeff(f, 1, coeff(h, 1));
            SetCoeff(f, 0, coeff(h, 0));
            break;
        }
    }
}

// for more randomly computatins
void init_seed_time()
{
	pid_t pid = getpid();
	long int tm = time(NULL);
	struct timeval tv;
	struct timezone tz;
	gettimeofday(&tv, &tz);
	tm += tv.tv_usec;
	SetSeed( to_ZZ(pid + tm) );
}

// find equivalent class representative of divisor
void equiv_class_repr(ZZ_pJac2& D)
{
	// D = (u(x), v(x)); -D = (u(x), -v(x))
	// v(x) = v1*x + v0
	// first we check the highest coefficient (v1), then the lowest one (v0)
	if (rep(D.v1) > rep(-D.v1))
		D = -D;
	else if ((D.v1 == -D.v1) && (rep(D.v0) > rep(-D.v0)))
		D = -D;
}

// find equivalent class representative of divisor
// with change and its logarithms
void equiv_class_repr(ZZ_pJac2& D, ZZ& alpha, ZZ& beta)
{
	//  D = (u(x), v(x)) = (x^2 + u1x + u0, v1x + v0)
	// -D = (u(x), -v(x))
	if (rep(D.v1) > rep(-D.v1))
	{
		D = -D;
		alpha = -alpha;
		beta = -beta;
	}
	else if ((D.v1 == -D.v1) && (rep(D.v0) > rep(-D.v0)))
	{
		D = -D;
		alpha = -alpha;
		beta = -beta;
	}
} 

// get some hash-value of a divisor
unsigned int hash_value(const ZZ_pJac2& D)
{
	ZZ h = rep(D.u1) ^ rep(D.v0);
	ZZ hh = rep(D.u0) ^ rep(D.v1);
	h ^= hh << NumBits(h);
	unsigned int x;
	unsigned char *str = (unsigned char *) (&x);
	BytesFromZZ(str, h, 4);
	return x;
}

// another variant of the hash function
unsigned int hash_value_2(const ZZ_pJac2& D)
{
	ZZ h = rep(D.u0) ^ rep(D.u1);
	ZZ hh = rep(D.v0) ^ rep(D.v1);
	h ^= hh << NumBits(h);
	unsigned int x;
	unsigned char *str = (unsigned char *) (&x);
	BytesFromZZ(str, h, 4);
	return x;
}

// get index of O-table 
inline int index_from_k_kp_b_r(int k, int kp, int b, int r)
{
	return (k + r*kp + r*r*b);
}

inline int is_distinguished(const ZZ_pJac2& D, const int proba)
{
	// checking "proba"-bits of the hash values
	return (hash_value(D) & (((1UL)<<proba) - 1)) == 0;
}

// computing all possible walk steps
void init_random_walk(RWData &rwdata, int r)
{
	// memory allocation for steps
	rwdata.r = r;
	rwdata.O = (DivAndTrack **) malloc(2*r*r*sizeof(DivAndTrack *));
	ZZ ak, bkp;
	ZZ_pJac2 akmP, bkpmP;
	for (int k = 0; k < r; ++k) {
		// select a random alpha_k of average value l1
		// if l1 is too small, then alpha_k is 0, 1 or 2
		if (rwdata.l1 <= 2) {
			if (k == 0) ak = 0;
			else if (1) ak = 1;
			else ak = 2;
		} else {
			ak = RandomBnd(to_ZZ(2*rwdata.l1));
		}
		// alpha_k*(p+1)*m*BaseDivisor
		akmP = ak*rwdata.pp1mBaseDivisor;
		for (int kp = 0; kp < r; kp++) {
			// select a random vetakp of average value l2
			bkp = RandomBnd(to_ZZ(2*rwdata.l2)) + 1;
			bkpmP = bkp*rwdata.mBaseDivisor;

			// write the corresponding value in the rwdata.O table
			rwdata.O[index_from_k_kp_b_r(k, kp, 0, r)] = new DivAndTrack;
			rwdata.O[index_from_k_kp_b_r(k, kp, 0, r)]->D = akmP + bkpmP;
			rwdata.O[index_from_k_kp_b_r(k, kp, 0, r)]->alpha = -ak;
			rwdata.O[index_from_k_kp_b_r(k, kp, 0, r)]->beta = bkp;

			rwdata.O[index_from_k_kp_b_r(k, kp, 1, r)] = new DivAndTrack;
			rwdata.O[index_from_k_kp_b_r(k, kp, 1, r)]->D = bkpmP - akmP;
			rwdata.O[index_from_k_kp_b_r(k, kp, 1, r)]->alpha = ak;
			rwdata.O[index_from_k_kp_b_r(k, kp, 1, r)]->beta = bkp;

			// get a value of equivalent class representative for the table
			DivAndTrack O0 = *rwdata.O[index_from_k_kp_b_r(k, kp, 0, r)];
			DivAndTrack O1 = *rwdata.O[index_from_k_kp_b_r(k, kp, 1, r)];
			equiv_class_repr(O0.D, O0.alpha, O0.beta);
			equiv_class_repr(O1.D, O1.alpha, O1.beta);
		}
	}
}

void init_parameters(RWData &rwdata, const ZZ& m, const ZZ& s1, const ZZ& s2)
{
	ZZ p = ZZ_p::modulus();

	// select a random base divisor and get his value of equivalence class
	random(rwdata.BaseDivisor);
	equiv_class_repr(rwdata.BaseDivisor);

	// calculation of boundaries
	rwdata.B1max = (5*SqrRoot(p))/(2*m);
	rwdata.B1min = - rwdata.B1max;
	rwdata.B2min = - ((2*p)/m);
	rwdata.B2max = (3*p)/m;
	
	rwdata.K = p*p + 1 - s1*(p+1) + s2 + m*((rwdata.B2max + rwdata.B2min)/2);
	rwdata.KBaseDivisor = rwdata.K*rwdata.BaseDivisor;
	equiv_class_repr(rwdata.KBaseDivisor);
	rwdata.mBaseDivisor = m*rwdata.BaseDivisor;
	equiv_class_repr(rwdata.mBaseDivisor);
	rwdata.pp1mBaseDivisor = (p+1)*rwdata.mBaseDivisor;

	// proba of being distinguished
	// the integer pD means that the proba is 2^-pD
	ZZ aux;
	aux = 3 * SqrRoot(
		(rwdata.B1max - rwdata.B1min) * (rwdata.B2max - rwdata.B2min)
	);
	aux /= 1000;
	rwdata.pD = NumBits(aux) - 2 + PROBA_DIST_INC;

	// calculation of the average length of a walk chain
	rwdata.l1 = to_double(rwdata.B1max-rwdata.B1min)/9.0
		/ sqrt(to_double(1UL<<rwdata.pD));
	rwdata.l2 = to_double(rwdata.B2max-rwdata.B2min)/10.0
		/ to_double(1UL<<rwdata.pD);

	// in the case where l1 << 1, the random walk is different, and we have
	// to choose another value for l1
	if (rwdata.l1 <= 2)
	{
		rwdata.l1 = to_double(rwdata.B1max-rwdata.B1min)/9.0
			/ (to_double(1UL<<rwdata.pD));
		assert (rwdata.l1 <= 2);
	}
}

void select_random_start_point(DivAndTrack& D, const int b, 
	const RWData &rwdata)
{
	// choose a random point, according the bounds
	if (b)
	{
		// for a wild walk the boundaries are different
		D.alpha = (rwdata.B1min + RandomBnd(rwdata.B1max - rwdata.B1min)) 
			* WILD_BOUND_SCALE;
		D.beta = (rwdata.B2min + RandomBnd(rwdata.B2max - rwdata.B2min)) 
			* WILD_BOUND_SCALE;
		D.D = (-D.alpha)*rwdata.pp1mBaseDivisor + D.beta*rwdata.mBaseDivisor;
		D.D += rwdata.KBaseDivisor;
	}
	else
	{
		D.alpha = rwdata.B1min + RandomBnd(rwdata.B1max - rwdata.B1min);
		D.beta = rwdata.B2min + RandomBnd(rwdata.B2max - rwdata.B2min);
		D.D = (-D.alpha)*rwdata.pp1mBaseDivisor + D.beta*rwdata.mBaseDivisor;
	}
}

void jump_one_step(DivAndTrack& D1, const DivAndTrack& D, RWData &rwdata, 
	int shift)
{
	unsigned int hash = hash_value(D.D);
	int r = rwdata.r;
	int k, kp, b;

	unsigned int il1 = 0;
	if (rwdata.l1 <= 2)
		il1 = (unsigned int)(1/rwdata.l1);
	
	// Select some bits of hash: k, kp, b
	b = hash & 1; hash >>= 1;
	kp = (hash+shift) % r; hash >>= 4;
	if (rwdata.l1 <= 2) {
		if ((hash % il1) == 0) {
			hash >>= 3;
			k = (hash % (rwdata.r - 1)) + 1;
		} else k = 0;
	} else k = hash % r;

	// take a step
	D1.D = D.D + rwdata.O[index_from_k_kp_b_r(k, kp, b, r)]->D;
	D1.alpha = D.alpha + rwdata.O[index_from_k_kp_b_r(k, kp, b, r)]->alpha;
	D1.beta = D.beta + rwdata.O[index_from_k_kp_b_r(k, kp, b, r)]->beta;
}

void run_thread(DivAndTrack& distD, int b, RWData &rwdata)
{
	// buffer for temporary value
	int tot, sts;
	DivAndTrack D;
	ZZ jmp;
	tot = -1;
	sts = 0;
	jmp = 0;

	// store last QUEUE_SIZE points to detect side-to-side cycle
	unsigned int cycle[QUEUE_SIZE];
	unsigned int hashD = 0;

	// triggers for side-to-side cycles
	int shift = 0;
	bool isCycle = false;

	// max number of steps until a distignished point is reached 
	ZZ no_cycle_bound = to_ZZ((1UL)<<rwdata.pD)*10 * SEARCH_SPACE_SCALE;
	
	// counter for detection the out-of-bounds cycle
	ZZ cnt;
	do
	{
		select_random_start_point(D, b, rwdata);
		equiv_class_repr(D.D, D.alpha, D.beta);
		cycle[0] = hash_value(D.D);
		isCycle = false;
		shift = 0;
		tot++;
		
		// start walking until:
		// --- the distignished point is found
		// --- the walk in the search space
		// --- the walk not in the "side-to-side" loop
		cnt = 0;
		while (
			!is_distinguished(D.D, rwdata.pD) 
			&& (cnt < no_cycle_bound)
			&& !isCycle
		)
		{
			jump_one_step(D, D, rwdata, shift);
			equiv_class_repr(D.D, D.alpha, D.beta);
			hashD = hash_value(D.D);
			shift = 0;
			cnt++;

			// check the side-to-side walk (cycle)
			for (int i = 0; i < QUEUE_SIZE; ++i)
			{
				if (hashD == cycle[i])
				{
					jmp += cnt;
					sts++;

					switch (ALGORITHM_TYPE)
					{

						// exit from the while loop and start with a new point
						case 1:
							cnt = no_cycle_bound;
							break;

						// return the loop entry point
						case 2:
							isCycle = true;
							tot++;
							break;
						
						// jump out from the loop deterministically
						// (change 1 bit of hash of the loop entry point)
						case 3:
							shift = 1;
							tot++;
							break;

						// if the number of steps taken is less than expected
						// (divided by JUMP_OUT_SPACE)
						// then do a deterministic jump
						// else return the loop entry point
						case 4: case 5:
							if (cnt < no_cycle_bound*JUMP_OUT_SPACE) shift = 1;
							else isCycle = true;
							tot++;
							break;
					}

					break;
				}
			}
			cycle[cnt%QUEUE_SIZE] = hashD;
		}
	} while ((cnt == no_cycle_bound) && !isCycle);
	
	rwdata.tot_jumps += jmp;
	rwdata.tot_cycles += tot;
	rwdata.sts_cycles += sts;

	// return value
	distD = D;
}

int main(int argc, char **argv)
{
	assert (sizeof(int)==4);

	ZZ p, s1, s2, m, N, N1;
	ZZ nb_jump;
    ZZ_pX f;
    RWData rwdata;

	cin >> p;
	ZZ_p::init(p);
	
	cin >> s1;
	cin >> s2;
	cin >> m;
	cin >> N1;
	cin >> f;
	ZZ_pJac2::init(f);
	
	cerr << "--- Input parameters: " << endl;
	cerr << "--- p  : " << p << endl;
	cerr << "--- m  : " << m << endl;
	cerr << "--- s1 : " << s1 << " = " << s1%m << " (mod m)" << endl;
	cerr << "--- s2 : " << s2 << " = " << s2%m << " (mod m)" << endl;
	cerr << "--- N  : " << N1 << endl;
	cerr << "--- f  : x^5";
	for (int i = 4; i >= 0; --i)
	{
		cerr << " + " << coeff(f, i);
		if (i>0) cerr << 'x';
		if (i>1) cerr << '^' << i;
	}
	cerr << endl << endl;

	// initialize the parameters
	init_seed_time();
	init_parameters(rwdata, m, s1%m, s2%m);
	init_random_walk(rwdata, HASH_SIZE);

	// initialize arrays for distignished points
	vec_DivAndTrack ListWild;
	vec_DivAndTrack ListTame;
	
	// some information about tests
	nb_jump = 0;
	rwdata.tot_jumps  = 0;
	rwdata.tot_cycles = 0;
	rwdata.sts_cycles = 0;

	// buffer for temporary values
	DivAndTrack distP;
	DivAndTrack otherdistP;

	// start detection a collison
	clock_t timeGR = clock();
	int b = (RandomWord() & 1UL);
	int found = 0;
	while (!found)
	{
		// change type of the walk
		// b = {0: tame, 1: wild}
		b = !b;
		
		// recive a distinguished point
		run_thread(distP, b, rwdata);
		
		// check the tame-wild collision
		if (b)
		{
			for (int i = 0; i < ListTame.length(); ++i)
			{
				if (distP.D == ListTame[i].D)
				{
					found = 1;
					otherdistP = ListTame[i];
					break;
				}
			}
			append(ListWild, distP);
		}
		else
		{
			for (int i = 0; i < ListWild.length(); ++i)
			{
				if (distP.D == ListWild[i].D)
				{
					found = 1;
					otherdistP = ListWild[i];
					break;
				}
			}
			append(ListTame, distP);
		}
		// Out some information about the run
		// cerr << "--- We have now " << ListWild.length() << " wilds and " 
		// 	<< ListTame.length() << " tames\n";
		// if (rwdata.sts_cycles)
		// 	cerr << "--- Side-to-side cycles: " << rwdata.sts_cycles << endl;
		// if (rwdata.tot_cycles) 
		// 	cerr << "--- Out of range cycles: " 
		//		 << rwdata.tot_cycles - rwdata.sts_cycles << endl;
	}
	timeGR = clock() - timeGR;

	rwdata.distCnt = ListWild.length() + ListTame.length();

	// compute #|J(Fq)|
	if (!b) 
	{
		DivAndTrack tmp;
		tmp = distP;
		distP = otherdistP;
		otherdistP = tmp;
	}

	// accurate to minus sign
	N = rwdata.K - (distP.alpha - otherdistP.alpha)*(p+1)*m 
		 + (distP.beta - otherdistP.beta)*m;
	if (!IsZero(N*rwdata.BaseDivisor)) N = -N + 2*rwdata.K;

	cout << N << endl;

	// the expected runtime is about: (ratio * q^(3/4)/m) group operations
	double ratio = to_double(rwdata.tot_jumps) / pow(to_double(p), 0.75 ) 
		* to_double(m);

	// Output the average values
    cerr << endl << OUTPUT_SEP << OUTPUT_SEP << endl;
	cerr << "N                    : " << N
		 << " (" << std::boolalpha << (bool) (N==N1) << ')' << endl;
    cerr << "Time                 : " 
		 << (double) timeGR/1000000 << " [sec]" << endl;
    cerr << "Jumps                : " << rwdata.tot_jumps << endl;
    cerr << "Ratio                : " << ratio << endl;
    cerr << "Speed                : " 
		 << (double) rwdata.distCnt/timeGR*1000000 << " [dP/s]" << endl;
    cerr << "Cycles out of bounds : " 
		 << (rwdata.tot_cycles - rwdata.sts_cycles) << endl;
    cerr << "Side-to-side cycles  : " << rwdata.sts_cycles << endl;
    cerr << "Total cycles         : " << rwdata.tot_cycles << endl;
    cerr << OUTPUT_SEP << OUTPUT_SEP << endl << endl;

	return 0;
}