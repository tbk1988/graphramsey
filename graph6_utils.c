/*
 * Library of functions related to reading (and writing)
 * graph6-formatted graphs from strings.
 */

#include <igraph.h> 

#define BIAS6 63
#define SMALLN 62
#define TOPBIT6 32 

/*
 * graph6_order(str)
 * Get the order of the graph from graph6 string str.
 * (Basically as in NAUTY.)
 */
int graph6_order(char* str) { 
  char *p;
  int n; 
  p = str;
  n = *p++ - BIAS6;
  if (n > SMALLN) {
    n = *p++ - BIAS6;
    n = (n << 6) | (*p++ - BIAS6);
    n = (n << 6) | (*p++ - BIAS6);
  }
  return(n);
}

/*
 * read_graph6(g,str)
 * Makes a new graph corresponding to the graph with graph6-string
 * str and puts it in g.
 *
 * (Basically as in NAUTY.)
 * Allocates memory, remember to free the graph!
 */
int read_graph6(igraph_t* g, char* str) {
  int n = graph6_order(str);
  igraph_empty(g,n,IGRAPH_UNDIRECTED);

  char *p;
  int i, j, k, x=0; 
  p = str + 1;
  if (n > SMALLN)
    p += 3;

  k = 1;
  for (j = 1; j < n; ++j) { 
    for (i = 0; i < j; ++i) {
      if (--k == 0) {
        k = 6;
        x = *(p++) - BIAS6;
      }
      if (x & TOPBIT6) { 
        igraph_add_edge(g,i,j);
      }
      x <<= 1;
    }
  } 

  return(0);
} 
