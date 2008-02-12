#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

void hello() {

#pragma omp parallel
{
        fprintf(stderr, "Hello world\n");

}

}

