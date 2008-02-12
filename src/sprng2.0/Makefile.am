AM_CPPFLAGS = -I$(top_srcdir)/include/sprng2.0
noinst_LTLIBRARIES = libsprng.la

libsprng_la_SOURCES = primes_32.c primes_64.c fwrap_mpi.c cputime.c makeseed.c \
simple_mpi.c memory.c communicate.c checkid.c store.c lcg.c lfg.c lcg64.c \
mlfg.c cmrg.c sprng.c checkid.h cputime.h fwrap.h lcg64.h lfg.h mlfg.h \
primelist_32.h primes_32.h simple_.h sprng.h store.h cmrg.h fwrap_.h int64.h \
lcg.h memory.h multiply.h primelist_64.h primes_64.h sprng_f.h sprng_interface.h

SPRNG_CFLAGS = -DAdd__ -DLittleEndian
libsprng_la_CFLAGS = $(SPRNG_CFLAGS) 
libsprng_la_LDFLAGS = $(SPRNG_CFLAGS)
