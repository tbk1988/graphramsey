/*
 * 
 */ 


//#include "naututil.h"
#include <stdio.h> 
#include <assert.h> 
#include <limits.h>

#include "vn_graph.h"

//#define MAXN 200 
#define STDERRNUM 1
#define MAXDEG 90
#define MAXSNDVAL 300

/*
 * Write the complement of graph g to gc
 */
void graph_complement(graph_t gc, graph_t g) {
  int i,j;
  int n = nnodes(g);
  assert(nnodes(g) == nnodes(gc));
  graph_empty(gc);
  for(i = 0; i < n; i++) {
    for(j = i+1; j < n; j++) {
      if(!graph_has_edge(g,i,j))
        graph_add_edge(gc,i,j);
    }
  }
}

/*
 * Computes 2*l(g) for graph g, where
 *   l(G) = e(G) - 8n(G) + 19.5α(G)
 *
 * Uses very_nauty for independence number.
 */
int linv2(graph_t g, int n) {
  int i,j,e=0;

  /* setup very_nauty-style complement graph */
  graph_t gc = graph_new(n);
  graph_complement(gc,g);

  int k;
  k = graph_clique_number(gc); 
  graph_clear(gc);

  return(2*nedges(g) - 16*nnodes(g) + 39*k);
}

/*
 * Get the second valency of a vertex
 * g, m, n are graph stats.
 * v is the vertex to find the second valency of.
 */
int sndval(graph_t g, int n, int v) {
  int i;
  int d = graph_node_degree(g,v);
  int c = 0;
  for(i = 0; i < d; i++) {
    c += graph_node_degree(g,neighbour(g,v,i));
  }
  return(c);
}

int main(int argc, char* argv[]) {
  unsigned int num_row = 0;
  int n,i,j,l2,d,d2;
  /* NAUTY-setting up graph */
  graph_t g; 

  int table[MAXDEG][MAXSNDVAL];
  for(i = 0; i < MAXDEG; i++)
    for(j = 0; j < MAXSNDVAL; j++)
      table[i][j] = INT_MAX; 

  char* line;
  while((line = geng_getline(stdin)) != NULL) { 
    g = geng_stringtograph(line); 
    n = nnodes(g);
    num_row++;
    if(STDERRNUM)
      fprintf(stderr, "%u ;\n",num_row);

    l2 = linv2(g,n); 
    for(i = 0; i < n; i++) {
      d = graph_node_degree(g,i); 
      d2 = sndval(g,n,i); 
      if(table[d][d2] > l2)
        table[d][d2] = l2;
    }

    graph_clear(g);
    free(line); 
  } 

  for(i = 0; i < MAXDEG; i++) {
    for(j = 0; j < MAXSNDVAL; j++) {
      if(table[i][j] != INT_MAX)
        printf("(%d,%d) : %7.1f\n",i,j,table[i][j]*0.5);
    }
  }

  return(0);
}
