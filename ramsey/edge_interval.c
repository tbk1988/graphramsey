/*
 * This program computes bounds for the number of edges
 * using methods from the HW+-article.  
 */

#include <stdio.h> 
#include <stdlib.h>
//#include <limits.h>
#include <gmp.h> 
#include <stdint.h> 
#include <assert.h>

#define VERBOSE

/*
 * This is Radz. table with improvements in
 *   Bottom right 2x2 square (and R(10,13)) from trivial sum bound to HW+-like bound.
 *   R(5,7) <= 143 -> 142 (by using HW+ and new R(5,5) bound.
 */
uint64_t irub_data[90] =
  {2,3,4,5,6,7,8,9,10,11,12,13,14,15,
   6,9,14,18,23,28,36,42,50,59,68,77,87,
   18,25,41,61,84,115,149,191,238,291,349,417,
   48,87,142,216,316,442,633,848,1138,1461,1878,
   165,298,495,780,1171,1804,2566,3703,5033,6911,
   540,1031,1713,2826,4553,6954,10578,15263,22112,
   1870,3583,6090,10630,16944,27485,41525,63609,
   6588,12677,22325,38832,64864,104590,166404,
   23556,45881,81123,145567,242462,407396};

/*
 * Return the current upper bound of R(m,n).
 */
uint64_t irub(int m,int n) {
  int p1,p2;
  if(m <= n) {
    p1 = m; p2 = n;
  } else {
    p2 = m; p1 = n;
  }

  assert(p1 >= 2);
  assert(p1 <= 10);
  assert(p2 <= 15);

  return(irub_data[13*p1 + p2 - 28 - (p1-2)*(p1-3)/2]);
}

/*
 * Print the upper bound table that is used in this program.
 */
void print_ub_ramsey_table() { 
  int i,j;
  for(i = 2; i <= 10; i++) {
    printf("%-2d: ",i);
    for(j = 2; j <= 15; j++) {
      printf("%-6lu ",irub(i,j));
    }
    printf("\n");
  }
} 

int main(int argc, char* argv[]) {
#ifdef VERBOSE
  printf("Current upper bound Ramsey table:\n");
  print_ub_ramsey_table();
#endif

  if(argc != 3) {
    printf("Usage: %s m n\n",argv[0]);
    printf("   for use on e(m,n;p)-bounds where m <= n\n");
    return(1);
  }

  char* ptr;
  uint64_t m,n;
  m = strtol(argv[1],&ptr,10);
  n = strtol(argv[2],&ptr,10);
  assert(m <= n); // only m <=n

  int exact;

  mpz_t a,b,c,d,p; // alpha,beta,gamma,delta and p
  mpz_t lb1,lb2,ub1,ub2; // lower and upper bound variables
  mpz_t A,B,C; // equation solving variables
  mpz_inits(a,b,c,d,p,lb1,lb2,ub1,ub2,A,B,C,NULL); 
  mpz_set_ui(a,irub(m-2,n)-1);
  mpz_set_ui(b,irub(m,n-2)-1);
  mpz_set_ui(c,irub(m-1,n)-1);
  mpz_set_ui(d,irub(m,n-1)-1);

  mpz_add_ui(p,d,2);
  while(mpz_cmp_ui(p,irub(m,n)) < 0) {

    // setting lb1 to ceil( p(p-d-1)/2 )
    mpz_sub_ui(lb1,p,1);
    mpz_sub(lb1,lb1,d);
    mpz_mul(lb1,lb1,p);
    mpz_cdiv_q_ui(lb1,lb1,2);

    // setting ub1 to floor(pc/2)
    mpz_mul(ub1,p,c);
    mpz_fdiv_q_ui(ub1,ub1,2);

#ifdef VERBOSE
    gmp_printf("(lb1,ub1) = (%Zu,%Zu)\n",lb1,ub1);
#endif

    // solve quadratic equation in steps
    // set A to (a-b+3(p-1))p
    mpz_sub_ui(A,p,1);
    mpz_mul_ui(A,A,3);
    mpz_add(A,A,a);
    mpz_sub(A,A,b);
    mpz_mul(A,A,p); 

    // set B to 12*p*p*(p-1)*(p-2-b)
    mpz_sub_ui(C,p,1); // tmp C = p-1
    mpz_sub_ui(B,p,2);
    mpz_sub(B,B,b);
    mpz_mul(B,B,C);
    mpz_mul(B,B,p);
    mpz_mul(B,B,p);
    mpz_mul_ui(B,B,12); 

    // set C to sqrt(A^2-B) (truncated integer part)
    mpz_mul(C,A,A);
    mpz_sub(C,C,B); 
    exact = mpz_root(C,C,2); 

    // set lb2 to A-sqrt(A^2-B) rounded up 
    mpz_neg(lb2,C);
    if(exact == 0)
      mpz_add_ui(lb2,lb2,1); // add 1 if A^2-B is not perfect square
    mpz_add(lb2,lb2,A);
    // divide lb2 with 12 and round up
    mpz_cdiv_q_ui(lb2,lb2,12);

    // set ub2 to A+sqrt(A^2-B) rounded down
    mpz_add(ub2,A,C);
    // divide ub2 with 12 and round down
    mpz_fdiv_q_ui(ub2,ub2,12);

#ifdef VERBOSE
    gmp_printf("(lb2,ub2) = (%Zu,%Zu)\n",lb2,ub2);
#endif

    // set lb1 to max of lb1 and lb2.
    if(mpz_cmp(lb1,lb2) < 0)
      mpz_set(lb1,lb2);
    // set ub1 to min of ub1 and ub2.
    if(mpz_cmp(ub1,ub2) > 0)
      mpz_set(ub1,ub2);

    gmp_printf("%Zu : (%Zu,%Zu)\n",p,lb1,ub1);

    if(mpz_cmp(lb1,ub1) > 0)
      gmp_fprintf(stderr,"> p = %Zu is impossible, change irub table.\n",p);
    
    mpz_add_ui(p,p,1);
  }  

  mpz_clears(a,b,c,d,p,lb1,lb2,ub1,ub2,A,B,C,NULL);

  return(0);
}
