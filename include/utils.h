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
list* makeList();
node* makeNode(int);
void append(list*,node*);
node *getFirst(list*);
void deleteFirst(list*);
void printList(list*);







#endif
