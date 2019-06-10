#include <igraph.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "gutils.h"
#include "graph6_utils.h" 

#define MAXSTRLEN 80 

int main(int argc, char* argv[]) { 
  char line[MAXSTRLEN+2];
  int i,k;

  if(argc > 2 || (argc == 2 && strncmp(argv[1],"-h",2) != 0))  {
    printf("usage: %s [-h]\n",argv[0]);
    printf("  -h: use heuristic to get lower bound on indep.nr\n");
  } else { 
    igraph_t g; 
    while(scanf("%s\n",line) != EOF) {
      i = strlen(line);
      line[i] = '\n'; line[i+1] = '\0'; 
      read_graph6(&g,line);

      printf("%s",line);
      if(argc == 2) {
        k = heur_indep_lb(&g);
      } else {
        igraph_independence_number(&g,&k); 
      }
      printf(" > ind = %d\n",k); 

      igraph_destroy(&g); 
    }
  }

  return(0);
}
