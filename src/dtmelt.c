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

/* special value for ... on RHS */
#define ELLIPSIS	((void*)-1)

static int hdrp = 0;
static int cnmp = 0;

/* idvars (our left hand side) */
static size_t nlhs;
static size_t *lhs;
/* measure vars (our right hand side) */
static size_t nrhs;
static size_t *rhs;


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

static inline const char*
memchrnul(const char *s, int c, size_t z)
{
	return memchr(s, c, z) ?: s + z;
}


static int
chck(size_t ncol)
{
	for (size_t i = 0U; i < nlhs; i++) {
		if (UNLIKELY(lhs[i] >= ncol)) {
			return -1;
		}
	}
	for (size_t i = 0U; i < nrhs; i++) {
		if (UNLIKELY(rhs[i] >= ncol)) {
			return -1;
		}
	}

	/* they can't both be ellipses */
	if (UNLIKELY(!nlhs && !nrhs)) {
		return -1;
	}

	if (rhs == ELLIPSIS) {
		/* construct the set of measure vars */
		nrhs = ncol - nlhs;
		if (UNLIKELY(!nrhs)) {
			return -1;
		}
		if (UNLIKELY((rhs = calloc(nrhs, sizeof(*rhs))) == NULL)) {
			return -1;
		} else for (size_t i = 0U, j = 0U; j < nrhs; i++) {
			for (size_t k = 0U; k < nlhs; k++) {
				if (i == lhs[k]) {
					goto found;
				}
			}
			rhs[j++] = i;
		found:
			continue;
		}
	}
	return 0;
}

static void
phdr(const char *hdrs, const size_t *hoff, size_t nxph)
{
	for (size_t i, j = 0U; j < nlhs; j++) {
		i = lhs[j];
		const size_t of = hoff[i + 0U];
		const size_t eo = hoff[i + 1U];
		fwrite(hdrs + of, sizeof(*hdrs), eo - of - 1U, stdout);
		fputc('\t' + (!nrhs && j + 1U >= nlhs), stdout);
	}
	if (!nrhs) {
		return;
	} else if (nxph <= 1U) {
		fputs("variable\t", stdout);
	} else for (size_t j = 0U; j < nxph; j++) {
		fprintf(stdout, "variable%zu", j + 1U);
		fputc('\t', stdout);
	}
	fputs("value\n", stdout);
	return;
}

static size_t
cxph(char *restrict hn, const size_t *of)
{
/* check if we're dealing with headers from a cross-product (A*B) */
	size_t n = 0U;

	if (UNLIKELY(rhs == NULL)) {
		return 1U;
	}
	with (size_t i = rhs[0U]) {
		const size_t bo = of[i + 0U];
		const size_t eo = of[i + 1U];
		const char *const ep = hn + eo - 1U;

		for (char *restrict hp = hn + bo, *np;
		     (n++, np = memchr(hp, '*', ep - hp)); hp = np + 1U) {
			*np = '\t';
		}
	}
	/* check the rest */
	for (size_t j = 1U; j < nrhs; j++) {
		const size_t bo = of[rhs[j] + 0U];
		const size_t eo = of[rhs[j] + 1U];
		const char *const ep = hn + eo - 1U;
		size_t m = 0U;
		for (char *restrict hp = hn + bo, *np;
		     (m++, np = memchr(hp, '*', ep - hp)); hp = np + 1U) {
			*np = '\t';
		}
		/* barf if they're of different length */
		if (UNLIKELY(m != n)) {
			return 0U;
		}
	}
	return n;
}

static char*
mkhdrs(size_t *restrict of, size_t nc)
{
	size_t z = 64U;
	char *r = malloc(z * sizeof(*r));
	for (size_t i = 0U, n = 0U; i < nc; i++, n++) {
		int m = snprintf(r + n, z - n, "V%zu", i + 1U);
		if (n + m >= z) {
			z *= 2U;
			r = realloc(r, z * sizeof(*r));
			/* reprint */
			snprintf(r + n, z - n, "V%zu", i + 1U);
		}
		of[i + 0U] = n;
		of[i + 1U] = (n += m) + 1U;
	}
	return r;
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

static ssize_t
find_s(const char *ss, const size_t *of, size_t nc, const char *s, size_t z)
{
/* find S of size Z in {SS + OF} */
	for (size_t i = 0U; i < nc; i++) {
		const size_t bo = of[i + 0U];
		const size_t eo = of[i + 1U];
		if (z == eo - bo - 1U &&
		    !memcmp(ss + bo, s, z)) {
			return i;
		}
	}
	/* not found */
	return -1;
}

static int
snrf(const char *formula, const char *hn, const size_t *of, size_t nc)
{
	static const char *form;
	const char *elhs, *l;
	const char *erhs, *r;
	const char *on;
	size_t nl, nr;

	if (UNLIKELY((formula = form ?: formula) == NULL)) {
		return !form - 1;
	}
	if (UNLIKELY((elhs = strchr(l = form = formula, '~')) == NULL)) {
		elhs = l + strlen(l);
		r = erhs = elhs;
	} else {
		r = elhs + 1U;
		erhs = r + strlen(r);
	}

	for (nl = 0U, on = l; (on = memchr(on, '+', elhs - on)); nl++, on++);
	for (nr = 0U, on = r; (on = memchr(on, '+', erhs - on)); nr++, on++);

	if (!nl && *l == '~') {
		lhs = lhs ?: ELLIPSIS;
		goto rhs;
	}
	if (lhs == NULL) {
		lhs = calloc(nlhs = nl + 1U, sizeof(*lhs));
	}
	/* now try and snarf the whole shebang */
	if (!nl) {
		goto one_l;
	}
	for (nl = 0U; nl < nlhs && l < elhs; nl++, l = on + 1U) {
		/* try with numbers first */
		char *tmp;
		long unsigned int x;

	one_l:
		on = memchrnul(l, '+', elhs - l);
		if ((x = strtoul(l, &tmp, 10)) && tmp == on) {
			x--;
		} else if ((x = find_s(hn, of, nc, l, on - l)) < nc) {
			;
		} else {
			/* retry next time */
			return -1;
		}

		lhs[nl] = x;
	}

rhs:
	/* snarf right hand side */
	if (!nr && !memcmp(r, "...\0", 4U)) {
		rhs = rhs ?: ELLIPSIS;
		goto fin;
	} else if (r == erhs) {
		/* no right side */
		goto fin;
	}
	if (rhs == NULL) {
		rhs = calloc(nrhs = nr + 1U, sizeof(*rhs));
	}
	if (!nr) {
		goto one_r;
	}
	for (nr = 0U; nr < nrhs && r < erhs; nr++, r = on + 1U) {
		/* try with numbers first */
		char *tmp;
		long unsigned int x;

	one_r:
		on = memchrnul(r, '+', erhs - r);
		if ((x = strtoul(r, &tmp, 10)) && tmp == on) {
			x--;
		} else if ((x = find_s(hn, of, nc, r, on - r)) < nc) {
			;
		} else {
			return -1;
		}

		rhs[nr] = x;
	}

fin:
	/* all is good, forget about the formula then */
	form = NULL;
	return 0;
}

static int
proc1(void)
{
	int rc = 0;
	char *line = NULL;
	size_t llen = 0U;
	ssize_t nrd;
	size_t ncol;
	size_t *coff = NULL;
	/* line number */
	size_t nr = 0U;
	/* offsets for header and header buffer */
	char *hn = NULL;
	size_t *hoff = NULL;
	size_t nxph;
	/* constant dimension line */
	size_t ndln = 0U;
	size_t zdln = 0U;
	char *dln = NULL;

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
	} else if (UNLIKELY(chck(ncol) < !hdrp - 1)) {
		errno = 0, error("\
Error: fewer columns present than needed for id or measure vars");
		rc = -1;
		goto out;
	} else if (UNLIKELY(!(coff = calloc(ncol + 1U, sizeof(*coff))))) {
		error("\
Error: cannot allocate memory to hold one line");
		rc = -1;
		goto out;
	} else if (UNLIKELY(!(hoff = calloc(ncol + 1U, sizeof(*hoff))))) {
		error("\
Error: cannot allocate memory to hold a copy of the header");
		rc = -1;
		goto out;
	}

	if (!hdrp) {
		if (UNLIKELY((hn = mkhdrs(hoff, ncol)) == NULL)) {
		error("\
Error: cannot allocate memory to hold a copy of the header");
		rc = -1;
		goto out;
		}
		if (cnmp) {
			/* print col names */
			phdr(hn, hoff, 0U);
		}
		goto tok;
	}
	/* otherwise snarf col names as defined in header */
	if (UNLIKELY((hn = strndup(line, nrd)) == NULL)) {
		error("\
Error: cannot allocate memory to hold a copy of the header");
		rc = -1;
		goto out;
	}
	tokln1(hoff, ncol, line, nrd);
	for (size_t i = 1U; i <= ncol; i++) {
		hn[hoff[i] - 1U] = '\0';
	}
	/* we might need to rescan the formula now */
	if (UNLIKELY(snrf(NULL, hn, hoff, ncol)) < 0) {
		error("\
Error: cannot interpret formula");
		rc = -1;
		goto out;
	} else if (UNLIKELY(chck(ncol) < 0)) {
		errno = 0, error("\
Error: fewer columns present than needed for LHS~RHS and value");
		rc = -1;
		goto out;
	}
	/* check for cross-product headers (A*B) */
	if (UNLIKELY(!(nxph = cxph(hn, hoff)))) {
		errno = 0, error("\
Error: product headers must have same number of factors");
		rc = -1;
		goto out;
	}
	if (cnmp && (nrd = getline(&line, &llen, stdin)) > 0) {
		/* print col names */
		phdr(hn, hoff, nxph);
		goto tok;
	}

	while ((nrd = getline(&line, &llen, stdin)) > 0) {
		size_t v, i;
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
		for (i = 0U, ndln = 0U; i < nlhs; i++) {
			const size_t bo = coff[lhs[i] + 0U];
			const size_t eo = coff[lhs[i] + 1U];

			if (UNLIKELY(ndln + eo - bo > zdln)) {
				/* resize */
				while ((zdln = (zdln * 2U) ?: 256U) <=
				       ndln + eo - bo);
				dln = realloc(dln, zdln * sizeof(*dln));
			}
			memcpy(dln + ndln, line + bo, eo - bo - 1);
			ndln += eo - bo - 1;
			dln[ndln++] = '\t';
		}

		if (UNLIKELY(!nrhs)) {
			/* last tab to newline */
			dln[ndln - 1U]++;
			fwrite(dln, sizeof(*dln), ndln, stdout);
			continue;
		}
		for (i = 0U; i < nrhs; i++) {
			const size_t bo = coff[rhs[i] + 0U];
			const size_t eo = coff[rhs[i] + 1U];
			v = rhs[i];

			fwrite(dln, sizeof(*dln), ndln, stdout);
			/* header or index */
			with (size_t hb = hoff[v + 0U], he = hoff[v + 1U]) {
				fwrite(hn + hb, 1, he - hb - 1, stdout);
			}
			fputc('\t', stdout);
			fwrite(line + bo, sizeof(*line), eo - bo - 1, stdout);
			fputc('\n', stdout);
		}
	}
out:
	free(coff);
	free(dln);
	free(hoff);
	free(hn);
	free(line);
	return rc;
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

	/* overread and/or expect headers? */
	hdrp = argi->header_flag;
	/* memorise that we want col names for STCC() later on */
	cnmp = argi->col_names_flag;

	/* snarf formula */
	if (UNLIKELY(!argi->nargs ||
		     snrf(*argi->args, NULL, NULL, 0U) < 0 && !hdrp)) {
		error("\
Error: cannot interpret formula");
		rc = 1;
		goto out;
	}

	rc = proc1() < 0;

	if (lhs != ELLIPSIS) {
		free(lhs);
	}
	if (rhs != ELLIPSIS) {
		free(rhs);
	}
out:
	yuck_free(argi);
	return rc;
}

/* dtmelt.c ends here */
