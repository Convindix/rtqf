#include <stdio.h>
//#include <stdint.h>
#include <stdlib.h>
#include <string.h>
//#include "libdivide.h"

#define MAXSAFE (uint64)0x7FFFFFFFFFFFFFFF

typedef unsigned long long int uint64;

uint64 gcd(uint64 a, uint64 b){ //Euclidean algorithm, with accelerated subtraction
	uint64 aChunk = a;
	if(a == 0){
		return b;
	}
	if(b == 0){
		return a;
	}
	if(a == 1 || b == 1){
		return 1;
	}
	while((aChunk<<1) < b){
		aChunk = aChunk << 1;
	}
	if(aChunk < b){
		return gcd(a, b - aChunk);
	}else if(a > b){
		return gcd(b, a - b); //Large numbers in second argument for speed
	}else{ //a == b
		return a;
	}
}

uint64 absminus(uint64 a, uint64 b){
	if(a < b){
		return b - a;
	}else if(a > b){
		return a - b;
	}else{
		return 0;
	}
}

/*uint64 addModN(uint64 a, uint64 b, uint64 n){ //To prevent overflows in (a + b) % n
	if(a >= n){
		a = a % n;
	}
	if(b >= n){
		b = b % n;
	}
	if(a+b does not overflow){
		return (a + b) % n;
	}else{
		reconstruct residue mod n from (a + b), which is (a + b) % 2**64
	}
}

uint64 doubleModN(uint64 a, uint64 n){
	if(a < ((uint64)2<<61)){
		return (a<<1) % n;
	}else{
		if(a % 2 == 0){
			return addModN(a, a, n);
		}
	}
}*/

uint64 multModN(uint64 a, uint64 b, uint64 n){ //To prevent overflows, Egyptian multiplication
	if(n > MAXSAFE){
		//printf("n is 2**63 or larger, this may give incorrect answers!\n");
	}
	if(b == 0){
		return 0;
	}else if(b == 1){
		return a;
	}else{
		if((b % 2) == 0){
			uint64 half = multModN(a, b/2, n) % n;
			return (half << 1) % n; //If n >= 2**63, half may be >= 2**63, then (half << 1) overflows
		}else{
			return (multModN(a, b-1, n) % n + a) % n; //TODO: can (multModN(a, b-1, n) % n + a) overflow?
		}
	}
}

uint64 sqrModN(uint64 b, uint64 n){
	return multModN(b, b, n);
}

uint64 powModN(uint64 b, uint64 e, uint64 n){ //Exponentiation by squaring
	if(e == 0){
		return 1;
	}else{
		if((e % 2) == 0){
			return sqrModN(powModN(b, e/2, n), n);
		}else{
			return multModN(powModN(b, e-1, n), b, n);
		}
	}
}

int intLog(uint64 n){
	int log = 0;
	while(n > 1){
		n /= 2;
		log++;
	}
	return log;
}

int millerRabin(uint64 n){
	//Use bases <=37, enough for correctness on inputs <2^64 (https://math.stackexchange.com/a/2481258)
	if(n > MAXSAFE){
		printf("n is 2**63 or larger, this may give incorrect answers! Happened in millerRabin with n=%llu\n", n);
	}
	const uint64 bases[12] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
	uint64 x;
	int s = 0;
	uint64 d = n - 1;
	while((d % 2) == 0){ //Extracting s, d such that 2^s*d = n - 1
		d = d >> 1;
		s++;
	}
	if((n % 2) == 0 && n > 2){
		return 0;
	}
	for(int round = 0; round < 12; round++){
		x = powModN(bases[round], d, n);
		if(x == 1 || x == n - 1){
			return 1;
		}
		int i = 0;
		while(i < s - 1){
			x = sqrModN(x, n);
			i++;
			if(x == 1){
				return 0;
			}
			if(x == n - 1){
				return 1;
			}
		}
		if(i == s - 1){
			x = sqrModN(x, n);
			i++;
			if(x != -1){
				return 0;
			}else{
				return 1;
			}
		}
	}
	return 1;
}

uint64 pRhoSingleFactor(uint64 n){ //One factor via Pollard rho
	if(millerRabin(n)){ //Seems to not halt when n is prime
		return n;
	}
	if((n % 2) == 0){ //Always fails for even n when starting with x=y=2
		return 2;
	}
	uint64 x, y, d, con;
	con = 1; //The constant term of the polynomial (https://stackoverflow.com/a/48208795)
	while(1){
		x = 2;
		y = 2;
		d = 1;
		while(d == 1){
			x = (sqrModN(x, n) + con) % n; //To ensure that x*x calculation doesn't overflow
			y = (sqrModN(y, n) + con) % n;
			y = (sqrModN(y, n) + con) % n;
			d = gcd(absminus(x, y), n);
		}
		if(d != n){
			if(millerRabin(d)){ //If d isn't prime
				return d;
			}else{
				return pRhoSingleFactor(d);
			}
		}else{ //Failure
			if(con < 1000){ //Giving up
				con++;
			}else{ //Giving up
				return 1;
			}
			continue;
		}
	}
}

uint64* pollardRho(uint64 n){
	int maxLength = intLog(n); //To avoid recalculating
	uint64* factors = (uint64*)malloc(maxLength * sizeof(uint64) + 1); //+1 used later
	int i = 0;
	while(n > 1 && i < maxLength){ //Write to factorization
		factors[i] = pRhoSingleFactor(n);
		n /= factors[i];
		i++;
	}
	while(i < maxLength){ //Pad to length floor(log_2(n)) if necessary
		factors[i] = 1;
		i++;
	}
	factors[i] = 0; //Mark end of factor list with 0
	return factors;
}

void debugDisplayArray(uint64* arr, int length){
	printf("[");
	for(int i = 0; i < length; i++){
		if(i < length-1){
			printf("%llu, ", *(arr+i));
		}else{
			printf("%llu]\n", *(arr+i));
		}
	}
}

int isSumOfTwoSqrs(uint64 n){
	/*Returns true iff n is not divisible by a prime 4m+3, with odd exponent (https://en.wikipedia.org/wiki/Sum_of_two_squares_theorem)
	This is computationally difficult (https://mathoverflow.net/q/60367) so major improvements might not be possible.
	In "Primes of the form x**2+ny**2" (p.12) there is a generalization of the two-square theorem applicable to x**2+10*y**2, however it may be better not to use it because subtracting off 10*z**2 leaves fewer values of z to check.*/
	uint64 p = 3;
	//const struct libdivide_u64_t p_ld = libdivide_u64_gen(p);
	
	uint64* factors = pollardRho(n);
	
	//Also need to check if the last block of factors is a prime 4m+3 to an odd power
	int length = 0;
	int i = 0;
	while(factors[length] > 1){
		length++;
	}
	uint64* runningLowest = factors;
	uint64 temp = 0;
	int expcounter = 0;
	//Sort factors using insertion sort (https://stackoverflow.com/q/736920)
	//Exit if while sorting, a prime 4m+3 to an odd power is found
	while(i < length){
		for(int j = i+1; j < length; j++){ //j=i can be skipped
			if(factors[j] < *runningLowest){
				runningLowest = factors + j;
			}
		}
		if(*runningLowest < factors[i]){ //Swap contents, factors[i] is what gets replaced
			temp = factors[i];
			factors[i] = *runningLowest;
			*runningLowest = temp;
		}
		if(i > 0 && factors[i] != factors[i-1]){ //Found a new block of a factor
			if((factors[i-1] % 4) == 3 && (expcounter % 2) == 1){
				free(factors);
				return 0;
			}
			expcounter = 0;
		}
		expcounter++;
		i++;
		runningLowest = factors + i; //Reset
	}
	if(i > 0 && factors[i] != factors[i-1]){ //Checking block at end
			if((factors[i-1] % 4) == 3 && (expcounter % 2) == 1){
				free(factors);
				return 0;
			}
	}
	/*expcounter = 0;
	i = 0;
	while(i < length){
		temp = factors[i];
		i++;
		if(factors[i] == temp){
			expcounter++;
		}else{ //A new block of prime factors reached, temp is the factor making up the previous block
			if((temp % 4) == 3 && (expcounter % 2) == 1){
				return 0;
			}
			expcounter = 1;
		}
	}*/
	free(factors);
	return 1;
}

int smallSqrFree(uint64 N){
	//No quick squarefree check known (https://mathoverflow.net/a/16100)
	uint64 f = 3; //Skipping even factors since n odd
	uint64 fSquared = 9;
	while(f < 30 && fSquared <= N){ //Not entire squarefree-ness, just to cut down search space
		if(N % fSquared == 0){
			return 0;
		}
		fSquared = f*f; //High-performance multipliers in modern CPUs
		f += 2;
	}
	return 1;
}

int isInTQF(uint64 N){ //N of the form x**2 + y**2 + 10*z**2
	uint64 z = 0;
	uint64 tenZSquared = 0;
	if(!smallSqrFree(N)){ //Squarefree N that are not 5 mod 10 must be in TQF (Ono, Soundararajan, "Ramanujan's Termary Quadratic Form", 1997, https://uva.theopenscholar.com/files/ken-ono/files/025_8.pdf)
		return 1;
	}
	while(tenZSquared <= N){
		if(isSumOfTwoSqrs(N - tenZSquared)){
			return 1;
		}
		tenZSquared += 20*z + 10;
		z++;
	}
	return 0;
}

int main(){
	uint64 N = 30000000001; // Total, must be initialized to something 1 mod 10
	while(N < MAXSAFE){ //Every 10k+5 is in TQF, so unrolled loop that avoids 5 mod 10
		if(!isInTQF(N)){
			break;
		}
		//printf("%llu in TQF\n", N);
		N += 2;
		if(!isInTQF(N)){
			break;
		}
		//printf("%llu in TQF\n", N);
		N += 4;
		if(!isInTQF(N)){
			break;
		}
		//printf("%llu in TQF\n", N);
		N += 2;
		if(!isInTQF(N)){
			break;
		}
		//printf("%llu in TQF\n", N);
		N += 2;
		if(N % 1024 == 1){
			printf("Everything up to %llu in TQF\n", N);
		}
	};
	printf("Counterexample! %llu. The Generalized Riemann Hypothesis is false.\n", N);
	return 0;
}
