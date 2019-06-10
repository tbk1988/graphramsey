/*
 * Library of functions related to reading (and writing)
 * graph6-formatted graphs from strings.
 */

#include <igraph.h> 

int graph6_order(char* str);
int read_graph6(igraph_t* g, char* str); 
