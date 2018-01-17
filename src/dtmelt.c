/*** dtmelt.c -- streaming implementation of data.table's ?melt
 *
 * Copyright (C) 2017-2018 Sebastian Freundt
 *
 * Author:  Sebastian Freundt <freundt@ga-group.nl>
 *
 * This file is part of dtcl.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the author nor the names of any contributors
 *    may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ***/
#if defined HAVE_CONFIG_H
# include "config.h"
#endif	/* HAVE_CONFIG_H */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <errno.h>
#include "nifty.h"

static int hdrp = 0;
/* idvars (we call them dimensions) */
static size_t ndims;
static size_t *dims;
/* measure vars */
static size_t nmeas;
static size_t *meas;


static void
__attribute__((format(printf, 1, 2)))
error(const char *fmt, ...)
{
	va_list vap;
	va_start(vap, fmt);
	vfprintf(stderr, fmt, vap);
	va_end(vap);
	if (errno) {
		fputc(':', stderr);
		fputc(' ', stderr);
		fputs(strerror(errno), stderr);
	}
	fputc('\n', stderr);
	return;
}

static int
chck(size_t ncol)
{
	for (size_t i = 0U; i < ndims; i++) {
		if (UNLIKELY(dims[i] >= ncol)) {
			return -1;
		}
	}
	for (size_t i = 0U; i < nmeas; i++) {
		if (UNLIKELY(meas[i] >= ncol)) {
			return -1;
		}
	}
	if (!nmeas) {
		/* construct the set of measure vars */
		if (UNLIKELY(!(meas = calloc(ncol, sizeof(*meas))))) {
			return -1;
		}
		for (size_t i = 0U; i < ncol; i++) {
			size_t j;
			for (j = 0U; j < ndims; j++) {
				if (i == dims[j]) {
					break;
				}
			}
			meas[nmeas] = i;
			nmeas += j >= ndims;
		}
	}
	return 0;
}


static size_t
toklng(const char *ln, size_t lz)
{
	const char *const ep = ln + lz;
	size_t ncol = 1U;

	for (const char *lp = ln, *np;
	     (np = memchr(lp, '\t', ep - lp)); lp = np + 1U) {
		ncol++;
	}
	return ncol;
}

static size_t
tokln1(size_t *restrict c, size_t nc, const char *ln, size_t lz)
{
	const char *const ep = ln + lz;
	size_t j = 0U;

	lz -= ln[lz - 1] == '\n';

	c[j++] = 0U;
	for (const char *lp = ln, *np;
	     j < nc && lp < ep && (np = memchr(lp, '\t', ep - lp));
	     lp = np + 1U, j++) {
		c[j] = np + 1U - ln;
	}
	c[j] = lz + 1U;
	return j;
}

static int
proc1(void)
{
	int rc = 0;
	char *line = NULL;
	size_t llen = 0U;
	ssize_t nrd;
	size_t ncol;
	size_t *coff;
	/* line number */
	size_t nr = 0U;
	/* constant dimension line */
	size_t ndimln = 0U;
	size_t zdimln = 0U;
	char *dimln = NULL;

	/* probe */
	if (UNLIKELY((nrd = getline(&line, &llen, stdin)) < 0)) {
		error("\
Error: cannot read lines");
		rc = -1;
		goto out;
	} else if (UNLIKELY(nrd == 0)) {
		goto out;
	} else if (UNLIKELY(!(ncol = toklng(line, nrd)))) {
		errno = 0, error("\
Error: cannot determine number of columns");
		rc = -1;
		goto out;
	} else if (UNLIKELY(chck(ncol) < 0)) {
		errno = 0, error("\
Error: fewer columns present than needed for id or measure vars");
		rc = -1;
		goto out;
	} else if (UNLIKELY(!(coff = calloc(ncol + 1U, sizeof(*coff))))) {
		error("\
Error: cannot allocate memory to hold one line");
		rc = -1;
		goto out;
	}

	if (!hdrp) {
		goto tok;
	}

	while ((nrd = getline(&line, &llen, stdin)) > 0) {
	tok:
		nr++;
		size_t nf = tokln1(coff, ncol, line, nrd);

		if (UNLIKELY(nf < ncol)) {
			errno = 0, error("\
Error: line %zu has only %zu columns, expected %zu", nr, nf, ncol);
			rc = 2;
			break;
		}

		/* construct constant dimension prefix */
		ndimln = 0U;
		for (size_t i = 0U; i < ndims; i++) {
			const size_t of = coff[dims[i] + 0U];
			const size_t eo = coff[dims[i] + 1U];

			if (UNLIKELY(ndimln + eo - of > zdimln)) {
				/* resize */
				while ((zdimln = (zdimln * 2U) ?: 256U) <=
				       ndimln + eo - of);
				dimln = realloc(dimln, zdimln * sizeof(*dimln));
			}
			memcpy(dimln + ndimln, line + of, eo - of - 1);
			dimln[(ndimln += eo - of) - 1] = '\t';
		}

		for (size_t i = 0U; i < nmeas; i++) {
			const size_t of = coff[meas[i] + 0U];
			const size_t eo = coff[meas[i] + 1U];

			fwrite(dimln, sizeof(*dimln), ndimln, stdout);
			/* header or index */
			fprintf(stdout, "%zu", meas[i] + 1U);
			fputc('\t', stdout);
			fwrite(line + of, sizeof(*line), eo - of - 1, stdout);
			fputc('\n', stdout);
		}
	}
	free(coff);
	free(dimln);
out:
	free(line);
	return rc;
}

static ssize_t
rdvars(size_t *restrict *tgt, const char *const *v, size_t n)
{
	size_t m = 0U;

	if (UNLIKELY(!n)) {
		return 0U;
	} else if (UNLIKELY((*tgt = calloc(n, sizeof(**tgt))) == NULL)) {
		error("\
Error: not enough memory.");
		return -1;
	}
	for (size_t i = 0U; i < n; i++) {
		char *on;
		if (((*tgt)[m] = strtoul(v[i], &on, 10))) {
			(*tgt)[m]--;
			m++;
		}
	}
	return m;
}


#include "dtmelt.yucc"

int
main(int argc, char *argv[])
{
	static yuck_t argi[1U];
	int rc = 0;

	if (yuck_parse(argi, argc, argv) < 0) {
		rc = 1;
		goto out;
	}

	with (ssize_t nd = rdvars(&dims, argi->idvar_args, argi->idvar_nargs)) {
		if (UNLIKELY(nd < 0)) {
			rc = 1;
			goto out;
		}
		ndims = nd;
	}

	rc = proc1() < 0;

	free(dims);
	free(meas);
out:
	yuck_free(argi);
	return rc;
}

/* dtmelt.c ends here */
