/*
 * assuming:
 *  1) graphs have at most 32 vertices.
 *  2) graph6-strings are at most 80 characters long. 
 */
#include <igraph.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <gmp.h>
#include "gutils.h"
#include "graph6_utils.h" 

#define MAXSTRLEN 80 

mpq_t fv[32];

/*
 * set the values of fv[d] to f(d), for 3 <= d <= 31.
 */
void setf() {
  mpq_set_str(fv[3],"37/100",10);
  mpq_set_str(fv[4],"8/25",10);
  mpq_set_str(fv[5],"37/130",10);
  mpq_set_str(fv[6],"124/481",10);
  mpq_set_str(fv[7],"5689/24050",10);
  mpq_set_str(fv[8],"171317/781625",10);
  mpq_set_str(fv[9],"13116449/64093250",10);
  mpq_set_str(fv[10],"124457366/647341825",10);
  mpq_set_str(fv[11],"2867530417/15795140530",10);
  mpq_set_str(fv[12],"197154577787/1145147688425",10);
  mpq_set_str(fv[13],"2453943217169/14975008233250",10);
  mpq_set_str(fv[14],"17753564375308/113464485459625",10);
  mpq_set_str(fv[15],"768342600854861/5128594742775050",10);
  mpq_set_str(fv[16],"18953081894794169/131804884889318785",10);
  mpq_set_str(fv[17],"311002538839607809/2248436271641320450",10);
  mpq_set_str(fv[18],"48707606578280655002/365370894141714573125",10);
  mpq_set_str(fv[19],"17023372343913698583809/132264263679300675471250",10);
  mpq_set_str(fv[20],"660114575436650613731867/5303796973539957086397125",10);
  mpq_set_str(fv[21],"56510383731386642970756253/468855652460932206437505850",10);
  mpq_set_str(fv[22],"13288326468180780629463447368/113697495721776060061095168625",10);
  mpq_set_str(fv[23],"6837590688621251058569599536833/60259672732541311832380439371250",10);
  mpq_set_str(fv[24],"1917304866425735948081399691851533/17384915583338168463641756758605625",10);
  mpq_set_str(fv[25],"46710713417551189492499262874781017/435318286206787738329589589235484850",10);
  mpq_set_str(fv[26],"47380433857869324474544785319758686/453400738095377382844818695249881913",10);
  mpq_set_str(fv[27],"518684081635686817907311661380315069/5092039058609622915026425346652519946",10);
  mpq_set_str(fv[28],"39721720477518885725295404135017071211/399725066100855398829574389712222815761",10);
  mpq_set_str(fv[29],"1125991796339523814061015260253313263417/11605810539893801579810401246127296926578",10);
  mpq_set_str(fv[30],"17089977126815164134705063408043273036196/180290263731798538334640888323460250531841",10);
  mpq_set_str(fv[31],"16073968991669901183610349857803704174194121/173439233709990193877924534567168761011631042",10);
}

int main(void) { 
  char line[MAXSTRLEN+2];
  int i,j,k,n; 
  setf(); /* initialise f*/

  igraph_t g; 
  igraph_vector_t ds; 

  mpq_t sum;
  mpq_t zero;
  mpq_t alpha;
  mpq_init(sum); /* initalise sum to 0/1 */
  mpq_init(zero);
  mpq_init(alpha);
  
  while(scanf("%s\n",line) != EOF) {
    i = strlen(line);
    line[i] = '\n'; line[i+1] = '\0'; 
    read_graph6(&g,line); 

    n = igraph_vcount(&g);

    igraph_vector_init(&ds,n);
    igraph_degree(&g,&ds,igraph_vss_all(),IGRAPH_ALL,IGRAPH_NO_LOOPS); 
    for(i = 0; i < igraph_vector_size(&ds); i++) {
      j = *igraph_vector_e_ptr(&ds,i); 
      mpq_add(sum,sum,fv[j]);
    }
    igraph_vector_destroy(&ds);

    k = heur_indep_lb(&g); // get heur. lower bound on independence
    mpq_set_ui(alpha,k,1);
    if(mpq_cmp(alpha,sum) < 0) { // only if lb < sum we need to check deeper
      igraph_independence_number(&g,&k); 
      mpq_set_ui(alpha,k,1);
      if(mpq_cmp(alpha,sum) < 0) {
        printf(">> %s", line); 
        gmp_printf("  %Qd < %Qd\n", alpha,sum);
      }
    }

    igraph_destroy(&g); 
    mpq_set(sum,zero);
  } 

  mpq_clears(fv[3],fv[4],fv[5],fv[6],fv[7],fv[8],fv[9],fv[10],fv[11],fv[12],
    fv[13],fv[14],fv[15],fv[16],fv[17],fv[18],fv[19],fv[20],fv[21],fv[22],
    fv[23],fv[24],fv[25],fv[26],fv[27],fv[28],fv[29],fv[30],fv[31],NULL); 
  mpq_clears(zero,sum,NULL); 

  return(0);
}
