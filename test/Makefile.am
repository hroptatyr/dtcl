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
EXTRA_DIST += molten_02.csv
EXTRA_DIST += molten_03.csv
EXTRA_DIST += molten_04.csv
EXTRA_DIST += cast_01.csv
EXTRA_DIST += cast_02.csv

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
TESTS += dtcast_11.clit
TESTS += dtcast_12.clit
TESTS += dtcast_13.clit
TESTS += dtcast_14.clit
TESTS += dtcast_15.clit
TESTS += dtcast_16.clit
TESTS += dtcast_17.clit
TESTS += dtcast_18.clit
TESTS += dtcast_19.clit
TESTS += dtcast_20.clit
TESTS += dtcast_21.clit
TESTS += dtcast_22.clit
TESTS += dtcast_23.clit
TESTS += dtcast_24.clit
TESTS += dtcast_25.clit
TESTS += dtcast_26.clit
TESTS += dtcast_27.clit
TESTS += dtcast_28.clit
TESTS += dtcast_29.clit

TESTS += dtmelt_01.clit
TESTS += dtmelt_02.clit
TESTS += dtmelt_03.clit
TESTS += dtmelt_04.clit
TESTS += dtmelt_05.clit
TESTS += dtmelt_06.clit
TESTS += dtmelt_07.clit
TESTS += dtmelt_08.clit
TESTS += dtmelt_09.clit
TESTS += dtmelt_10.clit
TESTS += dtmelt_11.clit
TESTS += dtmelt_12.clit

TESTS += dtrbind_01.clit
TESTS += dtrbind_02.clit

if HAVE_ASM_COROUTINES
TESTS += dtmerge_01.clit
TESTS += dtmerge_02.clit
TESTS += dtmerge_03.clit
TESTS += dtmerge_04.clit
TESTS += dtmerge_05.clit
TESTS += dtmerge_06.clit
TESTS += dtmerge_07.clit
TESTS += dtmerge_08.clit
TESTS += dtmerge_09.clit
TESTS += dtmerge_10.clit
TESTS += dtmerge_11.clit
TESTS += dtmerge_12.clit
TESTS += dtmerge_13.clit
TESTS += dtmerge_14.clit
TESTS += dtmerge_15.clit
TESTS += dtmerge_16.clit
TESTS += dtmerge_17.clit
TESTS += dtmerge_18.clit
endif  HAVE_ASM_COROUTINES
EXTRA_DIST += merge_01.csv
EXTRA_DIST += merge_02.csv
EXTRA_DIST += merge_03.csv
EXTRA_DIST += merge_04.csv
EXTRA_DIST += merge_05.csv
EXTRA_DIST += merge_06.csv

if HAVE_ASM_COROUTINES
TESTS += dtchanges_01.clit
TESTS += dtchanges_02.clit
TESTS += dtchanges_03.clit
TESTS += dtchanges_04.clit
TESTS += dtchanges_05.clit
TESTS += dtchanges_06.clit
TESTS += dtchanges_07.clit
TESTS += dtchanges_08.clit
endif  HAVE_ASM_COROUTINES
EXTRA_DIST += changes_01.csv
EXTRA_DIST += changes_02.csv
EXTRA_DIST += changes_03.csv

## Makefile.am ends here
