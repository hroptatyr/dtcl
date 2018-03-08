/*** dtrbind.c -- streaming implementation of data.table's ?rbind
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
#include <stdarg.h>
#include <errno.h>
#include "nifty.h"

/* global line buffer */
static char *line;
static size_t llen;
/* union header */
char *hdr;
size_t nhdr;
size_t zhdr;
size_t *hof;
size_t nhof;
size_t zhof;
/* file permutations
 * length + length*beef ... */
size_t *perm;
size_t nperm;
size_t zperm;


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
addhdr(const char *s, size_t n)
{
	char pivot;

	if (UNLIKELY(!n)) {
		return -1;
	}

	pivot = s[n - 1U];
	for (const char *hp = hdr, *const ep = hdr + nhdr, *x;
	     hp < ep && (x = memchr(hp, pivot, ep - hp)) != NULL; hp = x + n) {
		const char *b = x - (n - 1U);
		if (b >= hdr && !memcmp(b, s, n) && b[n] == '\t') {
			/* found him */
			const size_t o = b - hdr;
			size_t i;
			for (i = 0U; i < nhof && hof[i] < o; i++);
			return i;
		}
	}
	/* otherwise add him */
	if (UNLIKELY(nhdr + n >= zhdr)) {
		while ((zhdr *= 2U) < nhdr + n);
		hdr = realloc(hdr, zhdr * sizeof(*hdr));
	}
	if (UNLIKELY(nhof + 1U >= zhof)) {
		zhof *= 2U;
		hof = realloc(hof, zhof * sizeof(*hof));
	}
	memcpy(hdr + nhdr, s, n);
	hof[nhof++] = nhdr;
	nhdr += n;
	hdr[nhdr++] = '\t';
	hof[nhof] = nhdr;
	return nhof - 1;
}


static int
proc_hdr(FILE *fp)
{
	ssize_t nrd;
	size_t ncol;

	/* just in case we have to leave all of a sudden */
	perm[nperm++] = 0U;

	if (UNLIKELY((nrd = getline(&line, &llen, fp)) < 0)) {
		return -1;
	}
	ncol = toklng(line, nrd);

	if (UNLIKELY(nperm + ncol >= zperm)) {
		while ((zperm *= 2U) < nperm + ncol);
		perm = realloc(perm, zperm * sizeof(*perm));
	}
	with (size_t coff[ncol + 1U]) {
		size_t nf = tokln1(coff, ncol, line, nrd);
		for (size_t i = 0U; i < nf; i++) {
			const char *bo = line + coff[i + 0U];
			const char *eo = line + coff[i + 1U];
			const ssize_t k = addhdr(bo, eo - bo - 1U);

			if (UNLIKELY(k < 0)) {
				return -1;
			}
			perm[nperm + i] = k;
		}
	}
	/* record permutation */
	perm[nperm - 1U] = ncol;
	nperm += ncol;
	return 0;
}

static int
proc_res(FILE *fp)
{
	static size_t iperm;
	size_t invp[nhof];
	ssize_t nrd;
	size_t ncol;
	int rc = 0;

	ncol = perm[iperm++];
	/* construct inverse perm */
	memset(invp, 0, sizeof(invp));
	for (size_t i = 0U; i < ncol; i++) {
		invp[perm[iperm + i]] = i + 1U;
	}
	/* lest we forget */
	iperm += ncol;

	/* check for identity mapping */
	for (size_t i = 0U; i < ncol; i++) {
		if (invp[i] != i + 1U) {
			goto non_triv;
		}
	}
	/* yay, identity mapping is much simpler */
	with (size_t coff[ncol + 1U]) {
		char res[nhof - ncol + 1U];
		size_t nr = 1U;

		/* prepare rest of line */
		memset(res, '\t', nhof - ncol);
		res[nhof - ncol] = '\n';
		while ((nrd = getline(&line, &llen, fp)) > 0) {
			nr++;
			if (UNLIKELY(tokln1(coff, ncol, line, nrd) < ncol)) {
				errno = 0, error("\
Error: line %zu has fewer than %zu columns", nr, ncol);
				rc = 2;
				break;
			}
			fwrite(line, 1, coff[ncol] - 1U, stdout);
			fwrite(res, 1, nhof - ncol + 1U, stdout);
		}
	}
	return rc;

non_triv:
	with (size_t coff[ncol + 1U]) {
		size_t nr = 1U;
		while ((nrd = getline(&line, &llen, fp)) > 0) {
			nr++;
			if (UNLIKELY(tokln1(coff, ncol, line, nrd) < ncol)) {
				errno = 0, error("\
Error: line %zu has fewer than %zu columns", nr, ncol);
				rc = 2;
				break;
			}

			for (size_t i = 0U; i < nhof; i++) {
				if (invp[i]) {
					const size_t j = invp[i] - 1U;
					const char *bo = line + coff[j + 0U];
					const char *eo = line + coff[j + 1U];

					fwrite(bo, 1, eo - bo - 1U, stdout);
				}
				fputc('\t' + (i + 1U >= nhof), stdout);
			}
		}
	}
	return rc;
}


#include "dtrbind.yucc"

int
main(int argc, char *argv[])
{
	static yuck_t argi[1U];
	FILE **fps;
	int rc = 0;

	if (yuck_parse(argi, argc, argv) < 0) {
		rc = 1;
		goto out;
	}

	/* we're oper'ing on descriptors */
	if (UNLIKELY((fps = calloc(argi->nargs, sizeof(*fps))) == NULL)) {
		error("\
Error: cannot allocate space for file descriptors");
		rc = 1;
		goto out;
	}

	/* prealloc some header space */
	if (UNLIKELY((hdr = malloc(zhdr = 256U)) == NULL)) {
		error("\
Error: cannot allocate space for header");
		rc = 1;
		goto out;
	}

	if (UNLIKELY((hof = malloc((zhof = 32U) * sizeof(*hof))) == NULL)) {
		error("\
Error: cannot allocate space for header");
		rc = 1;
		goto out;
	}

	if (UNLIKELY((perm = malloc((zperm = 32U) * sizeof(*perm))) == NULL)) {
		error("\
Error: cannot allocate space for header permutation");
		rc = 1;
		goto out;
	}

	/* snarf headers from files */
	for (size_t i = 0U; i < argi->nargs; i++) {
		if (UNLIKELY((fps[i] = fopen(argi->args[i], "r")) == NULL)) {
			error("\
Error: cannot open file `%s'", argi->args[i]);
		} else if (UNLIKELY(proc_hdr(fps[i]) < 0)) {
			error("\
Warning: header unreadable in file `%s'", argi->args[i]);
		}
	}

	if (argi->col_names_flag && nhdr) {
		hdr[nhdr - 1U] = '\n';
		fwrite(hdr, 1, nhdr, stdout);		
	}

	/* snarf residuals */
	for (size_t i = 0U; i < argi->nargs; i++) {
		rc |= proc_res(fps[i]) < 0;
	}

	/* release handles (again?) */
	for (size_t i = 0U; i < argi->nargs; i++) {
		if (fps[i]) {
			fclose(fps[i]);
		}
	}
	free(fps);
	free(hdr);
	free(hof);
	free(perm);
	/* line buffer */
	free(line);

out:
	yuck_free(argi);
	return rc;
}

/* dtrbind.c ends here */
