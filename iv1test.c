/*
 * For testing the non-negativity of the invariant
 *   iv1(-) := 19e(-) - 88n(-) + 175alpha(-) - ind(P3|-) 
 */

#include <igraph.h> 
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <gmp.h>
#include "gutils.h"
#include "graph6_utils.h" 

#define MAXSTRLEN 80 

int main(int argc, char* argv[]) { 
  char line[MAXSTRLEN+2];
  int i,k; 
  unsigned long int j,n; 

  igraph_integer_t t;

  igraph_t g; 
  igraph_vector_t ds;

  mpz_t sum; // for degree sum
  mpz_init(sum);
  mpz_t sumsq; // for square degree sum
  mpz_init(sumsq);
  mpz_t indp3; // for ind(P3|-)
  mpz_init(indp3);
  mpz_t value; // for invariant value
  mpz_init(value);
  while(scanf("%s\n",line) != EOF) {
    i = strlen(line);
    line[i] = '\n'; line[i+1] = '\0'; 
    read_graph6(&g,line);

    n = igraph_vcount(&g); // number of vertices

    k = heur_indep_lb(&g); // compute heur.ally the ind#
    igraph_vector_init(&ds,n); // initalize deg.seq. vector
    igraph_degree(&g,&ds,igraph_vss_all(),IGRAPH_ALL,IGRAPH_NO_LOOPS);
    for(i = 0; i < igraph_vector_size(&ds); i++) {
      j = *igraph_vector_e_ptr(&ds,i);
      mpz_add_ui(sum,sum,j);
      mpz_add_ui(sumsq,sumsq,j*j);
    }
    igraph_vector_destroy(&ds);
    // set indp3 = (1/2) (sumsq - sum): 
    mpz_sub(indp3,sumsq,sum);
    mpz_tdiv_q_ui(indp3,indp3,2);

    // set up heuristic invariant value
    mpz_set_ui(value,19*igraph_ecount(&g));
    mpz_sub_ui(value,value,88*n);
    t = 175*heur_indep_lb(&g);
    mpz_add_ui(value,value,t);
    mpz_sub(value,value,indp3);

    if(mpz_cmp_ui(value,0) < 0) { // heuristically negative
      mpz_sub_ui(value,value,t);
      igraph_independence_number(&g,&t);
      t = 175*t;
      mpz_add_ui(value,value,t);
      if(mpz_cmp_ui(value,0) < 0) { // actually negative
        gmp_printf("%s>> iv =  %Zd < 0\n", line,value);
      }
    } 

    igraph_destroy(&g); 
    mpz_set_ui(sum,0);
    mpz_set_ui(sumsq,0);
  }
  
  mpz_clears(sum,sumsq,indp3,value,NULL);

  return(0);
}
