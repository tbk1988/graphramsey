#include <igraph.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "graph6_utils.h" 

#define MAXSTRLEN 80 

void print_vector(igraph_vector_t *v) {
  long int i, n=igraph_vector_size(v);
  for (i=0; i<n; i++) {
    printf(" %li", (long int) VECTOR(*v)[i]);
  }
  printf("\n");
}

int main(void) { 
  char line[MAXSTRLEN+2]; 
  long int n,i,m,d;
  int alphacount[5];
//  int minval = INT_MAX;
  igraph_t g; 
  igraph_vector_ptr_t res;
  igraph_vector_ptr_init(&res,0);
  igraph_vector_t nbrs;
  igraph_vector_t isect;

  igraph_vs_t allvs;
  igraph_vit_t allvit;
  while(scanf("%s\n",line) != EOF) {
    i = strlen(line);
    line[i] = '\n'; line[i+1] = '\0'; 
    read_graph6(&g,line);

   printf("%s",line);
   igraph_largest_independent_vertex_sets(&g,&res);
   n = igraph_vector_ptr_size(&res);
   printf(" > %ld maximum independent sets found\n", n);
   for(i = 0; i<n; i++) {
     alphacount[0] = 0; alphacount[1] = 0; alphacount[2] = 0;
     alphacount[3] = 0; alphacount[4] = 0;
     igraph_vector_t* v;
     v = igraph_vector_ptr_e(&res,i); 
     igraph_vs_all(&allvs);
     igraph_vit_create(&g,allvs,&allvit);
     igraph_vector_sort(v);
     while(!IGRAPH_VIT_END(allvit)) {
       igraph_vector_t ston;
       igraph_vector_init(&ston,1);
       VECTOR(ston)[0]=IGRAPH_VIT_GET(allvit);
       d = igraph_degree(&g,&ston,igraph_vss_all(),IGRAPH_ALL,IGRAPH_NO_LOOPS);
       igraph_vector_destroy(&ston); 
       igraph_vector_init(&nbrs,d);
       igraph_vector_init(&isect,d);
       igraph_neighbors(&g,&nbrs,IGRAPH_VIT_GET(allvit),IGRAPH_ALL);
       igraph_vector_sort(&nbrs);
       igraph_vector_intersect_sorted(v,&nbrs,&isect);
       m = igraph_vector_size(&isect);
       alphacount[m]++; 
       igraph_vector_destroy(&isect);
       igraph_vector_destroy(&nbrs);
       IGRAPH_VIT_NEXT(allvit); 
     }
     if(alphacount[4] >= 2) {
       printf(" >");
       print_vector((igraph_vector_t*)v);
       printf(" > a[0] = %d, a[1] = %d, a[2] = %d, a[3] = %d, a[4] = %d\n",
         alphacount[0],alphacount[1],alphacount[2],alphacount[3],alphacount[4]);
     }
     igraph_vit_destroy(&allvit);
     igraph_vs_destroy(&allvs);
     igraph_vector_destroy(v);
     free(v);
   }

    igraph_destroy(&g); 
  } 
//  printf(":: minval = %d\n", minval);

  igraph_vector_ptr_destroy(&res);

  return(0);
}
