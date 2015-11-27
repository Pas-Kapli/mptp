bin_PROGRAMS = $(top_builddir)/bin/delimit

libparse_utree_a_SOURCES = parse_utree.y lex_utree.l
libparse_rtree_a_SOURCES = parse_rtree.y lex_rtree.l
noinst_LIBRARIES = libparse_utree.a libparse_rtree.a


AM_CFLAGS=-I${srcdir} -O3 -mtune=native -Wall -Wsign-compare -g ${LIBS}
AM_YFLAGS = -d -p `${SED} -n 's/.*_\(.*\)/\1_/p' <<<"$*"`
AM_LFLAGS = -o lex.yy.c

__top_builddir__bin_delimit_LDADD = libparse_utree.a libparse_rtree.a
__top_builddir__bin_delimit_SOURCES = arch.c \
auto.c \
bayes.c \
delimit.c \
delimit.h \
dp.c \
fasta.c \
likelihood.c \
maps.c \
multichain.c \
output.c \
priors.c \
random.c \
rtree.c \
score.c \
svg.c \
svg_landscape.c \
util.c \
utree.c
