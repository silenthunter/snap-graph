#ifndef _BINARY_HEAP_H
#define _BINARY_HEAP_H

typedef struct {
    /* value is the key */
    attr_id_t *vals;
    attr_id_t count;
    attr_id_t max_size;
} adj_bheap_t;

void adj_bheap_insert(adj_bheap_t* h, attr_id_t val);
void adj_bheap_delete(adj_bheap_t* h, attr_id_t val);
void adj_bheap_siftup(adj_bheap_t* h, attr_id_t val);

#endif
