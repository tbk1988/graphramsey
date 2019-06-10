#include <igraph.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "graph6_utils.h" 

#define MAXSTRLEN 80 

int main(void) { 
  char line[MAXSTRLEN+2];
  int i,k;
  int minind = INT_MAX; 
  igraph_t g; 
  while(scanf("%s\n",line) != EOF) {
    i = strlen(line);
    line[i] = '\n'; line[i+1] = '\0'; 
    read_graph6(&g,line); 

    igraph_independence_number(&g,&k); 
    if(k < minind) minind = k;

    igraph_destroy(&g); 
  } 
  printf(":: minimum independence = %d\n", minind);

  return(0);
}
