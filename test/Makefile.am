AM_CPPFLAGS = -I$(top_srcdir)/include -I$(top_srcdir)/include/sprng
AM_LDFLAGS = -L$(top_srcdir)/lib -lm

#bin_PROGRAMS = eval_vertex_betweenness eval_BFS eval_biconnected_components eval_hello
bin_PROGRAMS = eval_BFS eval_biconnected_components eval_modularity \
eval_vertex_cover eval_dyn_ds

#eval_vertex_betweenness_SOURCES = drive_vertex_betweenness.c
#eval_vertex_betweenness_LDADD = $(top_srcdir)/lib/libsnap.la

eval_BFS_SOURCES = drive_BFS.c
eval_BFS_LDADD = $(top_srcdir)/lib/libsnap.la

eval_biconnected_components_SOURCES = drive_biconnected_components.c
eval_biconnected_components_LDADD   = $(top_srcdir)/lib/libsnap.la

eval_modularity_SOURCES = drive_modularity.c
eval_modularity_LDADD = $(top_srcdir)/lib/libsnap.la

eval_dyn_ds_SOURCES = test_dyn_ds.c
eval_dyn_ds_LDADD = $(top_srcdir)/lib/libsnap.la

eval_vertex_cover_SOURCES = drive_vertex_cover.c
eval_vertex_cover_LDADD = $(top_srcdir)/lib/libsnap.la
