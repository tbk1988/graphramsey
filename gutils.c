/*
 * gutils.c
 *
 * A mix of utilitary functions concerning graphs.
 *
 * Author: Oliver Krüger
 *         okruger (at) math.su.se
 */

#include <igraph.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <assert.h> 

/*
 * lower_bound_NC4(d,d2,n,e)
 * Computes a lower bound on the number of cycles in a triangle-free
 * graph with valencies as in d, second valencies as in d2, n vertices
 * and e edges.
 *
 * Uses Cauchy-Schwarz.
 */
int lower_bound_NC4(int* d, int* d2, int n, int e) {
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
    // s+= binomial(n(G_v),3) + e(G_v)(n(G_v) - 2)
    s += ngv*(ngv-1)*(ngv-2)/6 + egv*(ngv - 3);
    // s+= ceil(4*e(G_v)^2/n(G_v))
    // using: ceil(x/y) = (x + y - 1)/y 
    t += (4*egv*egv + ngv - 1)/ngv; 
  }
  // returning c + ceil(s/4)
  return(c + (s+3)/4 + (t+7)/8); 
}

/*
 * upper_bound_NC4(d,d2,n,e)
 * Computes a lower bound on the number of cycles in a triangle-free
 * graph with valencies as in d, second valencies as in d2, n vertices
 * and e edges.
 *
 * Uses Σd^2 ≤ (Σd)^2.
 */
int upper_bound_NC4(int* d, int* d2, int n, int e) {
  int i,ngv,egv;
  // c = -binomial(n,4) + binomial(e,2)
  int c = -n*(n-1)*(n-2)*(n-3)/24 + e*(e-1)/2; 
  unsigned int s = 0;
  for(i = 0; i < n; i++) {
    // c += binomial(d(v),3) - binomial(d(v),2) 
    c += d[i]*(d[i]-1)*(d[i]-2)/6 - d[i]*(d[i]-1)/2;
    // ngv = n(G_v), egv = e(G_v)
    ngv = n - d[i] - 1;
    egv = e - d2[i];
    // s+= binomial(n(G_v),3) + e(G_v)(n(G_v) - 2)
    s += ngv*(ngv-1)*(ngv-2)/6 + egv*(ngv - 2);
    // s+= 2e(G_v)(e(G_v) - 1) 
    s += 2*egv*(egv-1);
  }
  return(c + s/4);
} 

/*
 * heur_indep_lb(g)
 * Computes a lower bound on the independence number of the graph g using
 * the "take minimum degree vertex"-heuristic.
 */
int heur_indep_lb(igraph_t* g) {
  igraph_t h;
  igraph_copy(&h,g);
//  int n = igraph_vcount(&g);
  int count = 0;
//  int d;

  igraph_vector_t vec;
  igraph_vit_t allv;
  igraph_integer_t current_mindeg_vx;

  while(igraph_vcount(&h) > 0) {
    count++;

    // find a minimum degree vertex
    int current_mindeg = INT_MAX;
    igraph_vector_init(&vec,1);
    igraph_vit_create(&h,igraph_vss_all(),&allv);
    while (!IGRAPH_VIT_END(allv)) { 
      igraph_degree(&h,&vec,igraph_vss_1(IGRAPH_VIT_GET(allv)),IGRAPH_ALL,1);
      if(current_mindeg > VECTOR(vec)[0]) {
        current_mindeg = VECTOR(vec)[0];
        current_mindeg_vx = IGRAPH_VIT_GET(allv);
      } 
      IGRAPH_VIT_NEXT(allv);
    } 
    igraph_vit_destroy(&allv); 

    // delete current_mindeg_vx and all its neighbours
    igraph_neighbors(&h,&vec,current_mindeg_vx,IGRAPH_ALL);
    igraph_vector_push_back(&vec,current_mindeg_vx);
    igraph_delete_vertices(&h, igraph_vss_vector(&vec));

    igraph_vector_destroy(&vec);
  }

  igraph_destroy(&h);
  return(count);
}
