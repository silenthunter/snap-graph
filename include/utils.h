#ifndef _UTILS_H
#define _UTILS_H

#include "graph_defs.h"

double get_seconds(void);
void prefix_sums(attr_id_t*, attr_id_t*, attr_id_t*, long);
void usage_graph_options();
void graph_ext_check(char*, char*);

#endif
