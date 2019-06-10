/*
 * This contains some code for using extended integers.  
 */

#include <gmp.h> 

#define EINT_PLUS_INFINITY 1
#define EINT_MINUS_INFINITY 2 

typedef struct eint eint;
struct {
  char indic; // indicator of +-inf
  mpz_t val; 
} eint; 

/* add arithmetic operations */
