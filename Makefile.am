ACLOCAL_AMFLAGS = -I ./m4
SUBDIRS = src lib test
EXTRA_DIST = include
lib_LTLIBRARIES = libsnap.la

libsnap_la_SOURCES = 

libsnap_la_LIBADD = \
	src/graph_generation/libgraph_gen.la \
	src/graph_kernels/libgraph_kernels.la \
	src/graph_metrics/libgraph_metrics.la \
	src/graph_partitioning/libgraph_partitioning.la \
	src/sprng2.0/libsprng.la

