#ifndef _primes_64_h_
#define _primes_64_h_

#ifndef ANSI_ARGS
#ifdef __STDC__
#define ANSI_ARGS(args) args
#else
#define ANSI_ARGS(args) ()
#endif
#endif

int getprime_64 ANSI_ARGS((int need, unsigned int *prime_array, int offset));
int init_prime_64(void);

#define MAXPRIMELONG 3037000501U  /* largest odd # < sqrt(2)*2^31+2 */
#define MINPRIMELONG 55108   /* sqrt(MAXPRIME) */
#define MAXPRIMEOFFSETLONG 146138719U /* Total number of available primes */

#endif
