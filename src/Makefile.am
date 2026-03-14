bin_PROGRAMS = $(top_builddir)/bin/mptp

AM_CFLAGS=-I${srcdir} -O3 -mtune=native -Wall -Wsign-compare -g

__top_builddir__bin_mptp_SOURCES = arch.c \
auto.c \
aic.c \
mptp.c \
mptp.h \
dp.c \
fasta.c \
likelihood.c \
maps.c \
multirun.c \
output.c \
random.c \
rtree.c \
svg.c \
svg_landscape.c \
treeparse.c \
util.c \
hash.c \
list.c

test: all
	cd $(top_srcdir)/test && python3 run_tests.py all --binary $(abs_top_builddir)/bin/mptp

test-smoke: all
	cd $(top_srcdir)/test && python3 run_tests.py smoke --binary $(abs_top_builddir)/bin/mptp

test-regression: all
	cd $(top_srcdir)/test && python3 run_tests.py regression --binary $(abs_top_builddir)/bin/mptp

test-parsing: all
	cd $(top_srcdir)/test && python3 run_tests.py parsing --binary $(abs_top_builddir)/bin/mptp

test-regenerate: all
	cd $(top_srcdir)/test && python3 run_tests.py regenerate --binary $(abs_top_builddir)/bin/mptp
