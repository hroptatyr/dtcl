include $(top_builddir)/version.mk

## Help the Developers and yourself. Just use the C locale and settings
## for the compilation. They can still be overriden by make LANG=<whatever>
## but that is general a not very good idea
LANG = C
LC_ALL = C

AM_CFLAGS = $(EXTRA_CFLAGS)
AM_CPPFLAGS = -D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=600 -D_BSD_SOURCE
AM_CPPFLAGS += -DTEST

EXTRA_DIST = $(BUILT_SOURCES) $(TESTS)
TESTS =
TEST_EXTENSIONS =
BUILT_SOURCES =
check_PROGRAMS =
CLEANFILES = $(check_PROGRAMS)

## summon our helper
include clitoris.am
AM_CLIT_LOG_FLAGS = --builddir $(top_builddir)/src --verbose
TESTS_ENVIRONMENT = LC_ALL=C

EXTRA_DIST += molten_01.csv
EXTRA_DIST += cast_01.csv

TESTS += dtcast_01.clit
TESTS += dtcast_02.clit
TESTS += dtcast_03.clit
TESTS += dtcast_04.clit
TESTS += dtcast_05.clit
TESTS += dtcast_06.clit
TESTS += dtcast_07.clit
TESTS += dtcast_08.clit
TESTS += dtcast_09.clit
TESTS += dtcast_10.clit

## Makefile.am ends here
