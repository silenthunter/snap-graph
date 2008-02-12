AM_CPPFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_srcdir)/lib -lm

bin_PROGRAMS = eval_vertex_betweenness eval_BFS eval_biconnected_components

eval_vertex_betweenness_SOURCES = drive_vertex_betweenness.c
eval_vertex_betweenness_LDADD = $(top_srcdir)/lib/libsnap.la

eval_BFS_SOURCES = drive_BFS.c
eval_BFS_LDADD = $(top_srcdir)/lib/libsnap.la

eval_biconnected_components_SOURCES = drive_biconnected_components.c
eval_biconnected_components_LDADD   = $(top_srcdir)/lib/libsnap.la
