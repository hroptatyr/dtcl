### Makefile.am
include $(top_builddir)/version.mk

AM_CFLAGS = $(EXTRA_CFLAGS)
AM_CPPFLAGS = -D_POSIX_C_SOURCE=201001L -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
AM_LDFLAGS = $(XCCLDFLAGS)

bin_PROGRAMS =
noinst_PROGRAMS =
lib_LIBRARIES =
noinst_LIBRARIES =
pkglib_LIBRARIES =
noinst_HEADERS =
BUILT_SOURCES =
EXTRA_DIST = $(BUILT_SOURCES)
CLEANFILES = 
SUFFIXES =

CLEANFILES += version.c
EXTRA_DIST += version.c.in

noinst_HEADERS += nifty.h

noinst_LIBRARIES += libcoru.a
libcoru_a_SOURCES =
EXTRA_libcoru_a_SOURCES =
libcoru_a_SOURCES += coru/cocore.c coru/cocore.h
libcoru_a_SOURCES += coru/switch.c coru/switch.h
libcoru_a_SOURCES += coru/platform.h
libcoru_a_SOURCES += coru/valgrind.h
libcoru_a_SOURCES += coru/memcheck.h
libcoru_a_SOURCES += coru.h coruaux.h
EXTRA_libcoru_a_SOURCES += coru/switch-arm.c
EXTRA_libcoru_a_SOURCES += coru/switch-ppc_osx.c
EXTRA_libcoru_a_SOURCES += coru/switch-x86.c
EXTRA_libcoru_a_SOURCES += coru/switch-x86_64.c


bin_PROGRAMS += dtcast
BUILT_SOURCES += dtcast.yucc

bin_PROGRAMS += dtmelt
BUILT_SOURCES += dtmelt.yucc

bin_PROGRAMS += dtrbind
BUILT_SOURCES += dtrbind.yucc

bin_PROGRAMS += lines
BUILT_SOURCES += lines.yucc


## version rules
version.c: version.c.in $(top_builddir)/.version
	$(AM_V_GEN) PATH="$(top_builddir)/build-aux:$${PATH}" \
		yuck$(EXEEXT) scmver --ignore-noscm --force -o $@ \
			--use-reference --reference $(top_builddir)/.version $<

## yuck rule
SUFFIXES += .yuck
SUFFIXES += .yucc
.yuck.yucc:
	$(AM_V_GEN) PATH=$(top_builddir)/build-aux:"$${PATH}" \
		yuck$(EXEEXT) gen -o $@ $<

## Create preprocessor output (debugging purposes only)
.c.i:
	$(COMPILE) -E -o $@ $<

## Create assembler output (debugging purposes only)
.c.s:
	$(COMPILE) -S -c $(AM_CFLAGS) $<
