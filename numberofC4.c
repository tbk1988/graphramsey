/*
 * numberofC4
 *
 * Program that computes the number of cycles of length four in
 * the non-existent (3,9;27,59)-graphs using the G_v-reconstruction
 * formula and Cauchy-Schwarz.
 *
 * Reads the output of "geng -tD5 10" from stdin.
 */

#include <igraph.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h>
#include "gutils.h"
#include "graph6_utils.h" 

#define MAXSTRLEN 161

/*
 * modified_upper_bound_NC4(d, d2,n,e)
 *
 * A specialised way to get the upper bound, just for
 * the (3,9;27,59)-case.
 */
int modified_upper_bound_NC4(int* d, int* d2, int n, int e) {
  int i,ngv,egv;
  // c = -binomial(n,4) + binomial(e,2)
  int c = -n*(n-1)*(n-2)*(n-3)/24 + e*(e-1)/2; 
  unsigned int s = 0, t = 0;
  for(i = 0; i < n; i++) {
    // c += binomial(d(v),3) - binomial(d(v),2) 
    c += d[i]*(d[i]-1)*(d[i]-2)/6 - d[i]*(d[i]-1)/2;
    // ngv = n(G_v), egv = e(G_v)
    ngv = n - d[i] - 1;
    egv = e - d2[i];
    // s+= binomial(n(G_v),3) + e(G_v)(n(G_v) - 3)
    s += ngv*(ngv-1)*(ngv-2)/6 + egv*(ngv - 3);
    // s+= upper bound
    if(d[i] == 4) {
      assert(d2[i] == 16 || d2[i] == 17);
      t += 354 + 32*(17-d2[i]);
    } else {
      assert(d2[i] >= 20 && d2[i] <= 24); 
      t += 245 + 14*(24-d2[i]);
    }
  }
  return(c + s/4 + t/8);
}

/*
 * modified_lower_bound_NC4(d, d2,n,e)
 *
 * A specialised way to get the lower bound, just for
 * the (3,9;27,59)-case.
 */
int modified_lower_bound_NC4(int* d, int* d2, int n, int e) {
  int i,ngv,egv;
  // c = -binomial(n,4) + binomial(e,2)
  int c = -n*(n-1)*(n-2)*(n-3)/24 + e*(e-1)/2; 
  unsigned int s = 0, t = 0;
  for(i = 0; i < n; i++) {
    // c += binomial(d(v),3) - binomial(d(v),2) 
    c += d[i]*(d[i]-1)*(d[i]-2)/6 - d[i]*(d[i]-1)/2;
    // ngv = n(G_v), egv = e(G_v)
    ngv = n - d[i] - 1;
    egv = e - d2[i];
    // s+= binomial(n(G_v),3) + e(G_v)(n(G_v) - 3)
    s += ngv*(ngv-1)*(ngv-2)/6 + egv*(ngv - 3);
    // s+= upper bound
    if(d[i] == 4) {
      assert(d2[i] == 16 || d2[i] == 17);
      t += 324 + 14*(17-d2[i]); 
    } else {
      assert(d2[i] >= 20 && d2[i] <= 24);
      t += 238 + 10*(24-d2[i]);
    }
  }
  return(c + s/4 + t/8);
} 

int main(void) { 
  char line[MAXSTRLEN+2];
  int i,j,n,count;
  int lbmin=INT_MAX,ubmax=INT_MIN;
  int* snd_vals_five;
  int* vals;
  int* sndvals;
  igraph_t g;
  igraph_vector_t ds; 
  while(scanf("%s\n",line) != EOF) {
    i = strlen(line);
    line[i] = '\n'; line[i+1] = '\0'; 
    read_graph6(&g,line);

    snd_vals_five = calloc(26,sizeof(int));
    n = igraph_vcount(&g);

    vals = calloc(27,sizeof(int));
    sndvals = calloc(27,sizeof(int)); 

    igraph_vector_init(&ds,n); 
    igraph_degree(&g,&ds,igraph_vss_all(),IGRAPH_ALL,IGRAPH_NO_LOOPS); 
    count = 50-igraph_vector_sum(&ds); 

    // 25 - (5 - d) gives snd val
    igraph_vector_add_constant(&ds,20); 
    for(i = 0; i < igraph_vector_size(&ds); i++) { 
      j = *igraph_vector_e_ptr(&ds,i);
      assert(0 <= j && j <= 25);
      snd_vals_five[j]++; 
    } 
    igraph_vector_destroy(&ds); 

    // at least one tetraval with sndval 17.
    // at most 17 vertices with sndval 17
    if(17-count >= 0 && count > 0) {
      for(i = 0; i < 27; i++) {
        if(i < 17) {
          vals[i] = 4;
          if(i < count)
            sndvals[i] = 17;
          else
            sndvals[i] = 16;
        } else {
          vals[i] = 5;
          j = 0;
          while(snd_vals_five[j] == 0 && j <= 25) j++;
          sndvals[i] = j;
          assert(j <= 26);
          snd_vals_five[j]--;
        }
      } 

      int lb = modified_lower_bound_NC4(vals,sndvals,27,59);
      int ub = modified_upper_bound_NC4(vals,sndvals,27,59);

      if(lb < lbmin) lbmin = lb;
      if(ub > ubmax) ubmax = ub;

      printf("%s",line);
      printf(" > lower_bound = %d\n", lb);
      printf(" > upper_bound = %d\n", ub); 
    }
    igraph_destroy(&g);
    free(snd_vals_five); 
    free(vals);
    free(sndvals); 
  } 
  printf(":: lbmin = %d, ubmax = %d\n",lbmin,ubmax);

  return(0);
}
