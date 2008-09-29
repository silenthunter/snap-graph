#ifndef _UTILS_H
#define _UTILS_H

#include "graph_defs.h"


//utils.c
double get_seconds(void);
void prefix_sums(attr_id_t*, attr_id_t*, attr_id_t*, long);
void usage_graph_options();
void graph_ext_check(char*, char*);
void print_graph(graph_t*);




//vectorUtils.c

void printDoubleVector(double*,int,int);
void printIntVector(int*,int,int);
void print_attr_id_t_Vector(attr_id_t*,attr_id_t,attr_id_t);



//list.c
list_t* makeList();
node_t* makeNode(int);
void append(list_t*,node_t*);
node_t *getFirst(list_t*);
void deleteFirst(list_t*);
void printList(list_t*);







#endif
