
#include "utils.h"


list* makeList()
{
	list *l  = (list*)malloc(sizeof(list));
	l->head=NULL;
	l->tail=NULL;
	l->size=0;
	return l;
}

node* makeNode(int id)
{
	node *new = (node*)malloc(sizeof(node));
	new->id = id;
	new->next = NULL;
	return new;
}

//note append will append at the last.
void append(list *L, node *n)
{
	if(L->size==0)
	{
		L->head=n;
		L->tail=n;
	}
	else
	{
		L->tail->next = n;
		L->tail = n;
	}
	L->size++;
	
}

node* getFirst(list *L)
{
	return L->head;
}

void deleteFirst(list *L)
{
	node *n = L->head;
	L->head = n->next;
	free(n);
	L->size--;
}

void printList(list *L)
{
	node *n;
	printf("Printing list of size:%d\n",L->size);
	if(L->size==0) return;
	for(n=L->head; n!=L->tail; n=n->next)
	{
		printf("%d,",n->id);
	}
	printf("%d,",n->id);
	printf("\n\n\n\n");
}
