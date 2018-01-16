/*** dtcast.c -- streaming implementation of data.table's ?dcast
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
/* number of cast columns */
static size_t ncc;
static size_t zcc;
/* cast column hash values */
static uint64_t *cc;
/* number of distinct values per cast column */
static size_t *nccv;
/* total size of allocated ccv */
static size_t *zccv;
/* base pointers per per cast column */
static char **ccv;
static size_t *zccvo;
/* offsets per cast column
 * for each molten line we determine c<-CC(id-col), the index of the
 * cast column as per id variable then store CCV[c] + CCVO[c][nccv]
 * the value */
static size_t **ccvo;
/* buffer of the dimension line (LHS) */
static size_t ndim;
static size_t zdim;
static char *dim;

static size_t nlhs;
static union {
	size_t v;
	size_t *p;
} lhs = {0U};
static size_t rhs = 1U;
/* this one is 1-based */
static size_t vhs;


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

/* murmur2 */
#define HASHSIZE	(64U / 8U)

static uint64_t
MurmurHash64A(const void *key, size_t len)
{
	const uint64_t m = 0xc6a4a7935bd1e995ULL;
	const int r = 47;

	uint64_t h = (len * m);

	const uint64_t * data = (const uint64_t *)key;
	const uint64_t * end = data + (len/8);

	while(data != end) {
		uint64_t k = *data++;

		k *= m;
		k ^= k >> r;
		k *= m;

		h ^= k;
		h *= m;
	}

	const unsigned char * data2 = (const unsigned char*)data;

	switch(len & 7) {
	case 7: h ^= (uint64_t)(data2[6]) << 48;
	case 6: h ^= (uint64_t)(data2[5]) << 40;
	case 5: h ^= (uint64_t)(data2[4]) << 32;
	case 4: h ^= (uint64_t)(data2[3]) << 24;
	case 3: h ^= (uint64_t)(data2[2]) << 16;
	case 2: h ^= (uint64_t)(data2[1]) << 8;
	case 1: h ^= (uint64_t)(data2[0]);
		h *= m;
		break;
	};

	h ^= h >> r;
	h *= m;
	h ^= h >> r;

	return h;
} 

#define hash(x, y)	MurmurHash64A((x), (y))


/* ccv operations */
static void
prnt(void)
{
	size_t m[ncc];
	size_t s[ncc];
	size_t ns = 0U;

	if (UNLIKELY(!ndim)) {
		/* no dimensions? */
		return;
	}
	if (!(rhs + 1U)) {
		/* just the dimension */
		goto more;
	}
	memset(m, 0, sizeof(m));
	for (size_t j = 0U; j < ncc; j++) {
		if (nccv[j]) {
			/* at least one value */
			goto prnt;
		}
	}
	return;

prnt:
	/* determine order of stepping */
	memset(s, 0, sizeof(s));
	for (size_t j = ncc; j > 0U; j--) {
		if (nccv[j - 1U] > 1U) {
			s[ns++] = j - 1U;
		}
	}
more:
	/* dimension line */
	fwrite(dim, 1, ndim, stdout);
	for (size_t j = 0U; j < ncc; j++) {
		fputc('\t', stdout);
		if (nccv[j]) {
			fwrite(ccv[j] + ccvo[j][m[j]], 1,
			       ccvo[j][m[j] + 1U] - ccvo[j][m[j]], stdout);
		}
	}
	if (ns) {
		/* multi-step */
		for (size_t i = 0U; i < ns; i++) {
			if (++m[s[i]] < nccv[s[i]]) {
				fputc('\n', stdout);
				goto more;
			}
			m[s[i]] = 0U;
		}
	}
	fputc('\n', stdout);
	return;
}

static void
rset(const char *line, const size_t *coff)
{
	memset(nccv, 0, ncc * sizeof(*nccv));

	for (size_t j = 0U; j < ncc; j++) {
		ccvo[j][0U] = 0U;
		ccvo[j][1U] = 0U;
	}

	/* make up dimension line */
	if (!nlhs) {
		const size_t of = coff[lhs.v + 0U];
		const size_t eo = coff[lhs.v + 1U];
		if (UNLIKELY(eo - of >= zdim)) {
			while ((zdim = (zdim * 2U) ?: 64U) < eo - of);
			dim = realloc(dim, zdim * sizeof(*dim));
		}
		memcpy(dim, line + of, ndim = eo - of);
	} else for (size_t i = 0U, n = 0U; i < nlhs; i++, ndim = n) {
		const size_t of = coff[lhs.p[i] + 0U];
		const size_t eo = coff[lhs.p[i] + 1U];
		if (UNLIKELY(n + eo - of >= zdim)) {
			while ((zdim = (zdim * 2U) ?: 64U) < n + eo - of);
			dim = realloc(dim, zdim * sizeof(*dim));
		}
		memcpy(dim + n, line + of, eo - of - 1U);
		n += eo - of - 1U;
		dim[n++] = '\t';
	}
	/* omit trailing separator */
	ndim--;
	return;
}

static void
bang(const char *str, size_t len, size_t j)
{
	size_t eo = ccvo[j][nccv[j]];

	if (UNLIKELY(eo + len >= zccv[j])) {
		/* resize */
		while ((zccv[j] = (zccv[j] * 2U) ?: 64U) < ccvo[j][nccv[j]]);
		ccv[j] = realloc(ccv[j], zccv[j] * sizeof(*ccv[j]));
	}

	memcpy(ccv[j] + eo, str, len);
	eo += len;

	if (UNLIKELY(++nccv[j] >= zccvo[j])) {
		zccvo[j] *= 2U;
		ccvo[j] = realloc(ccvo[j], zccvo[j] * sizeof(*ccvo[j]));
	}
	/* store current end */
	ccvo[j][nccv[j]] = eo;
	return;
}

static void
mtcc(void)
{
	zcc = 0U;
	return;
}

static int
stcc(char *const *args, size_t nargs)
{
	if (!(ncc = nargs)) {
		;
	} else if (UNLIKELY((cc = calloc(ncc, sizeof(*cc))) == NULL)) {
		return -1;
	} else if (UNLIKELY((nccv = calloc(ncc, sizeof(*nccv))) == NULL)) {
		return -1;
	} else if (UNLIKELY((ccv = calloc(ncc, sizeof(*ccv))) == NULL)) {
		return -1;
	} else if (UNLIKELY((ccvo = calloc(ncc, sizeof(*ccvo))) == NULL)) {
		return -1;
	} else if (UNLIKELY((zccv = calloc(ncc, sizeof(*zccv))) == NULL)) {
		return -1;
	} else if (UNLIKELY((zccvo = calloc(ncc, sizeof(*zccvo))) == NULL)) {
		return -1;
	}
	for (size_t j = 0U; j < ncc; j++) {
		ccvo[j] = calloc(8U, sizeof(*ccvo[j]));
		zccvo[j] = 8U;
	}
	for (size_t i = 0U; i < ncc; i++) {
		const char *c = args[i];
		cc[i] = hash(c, strlen(c));
	}
	return 0;
}

static int
adcc(uint64_t c)
{
	if (UNLIKELY(ncc >= zcc)) {
		const size_t nuz = (zcc * 2U) ?: 64U;

		cc = realloc(cc, nuz * sizeof(*cc));
		nccv = realloc(nccv, nuz * sizeof(*nccv));
		ccv = realloc(ccv, nuz * sizeof(*ccv));
		ccvo = realloc(ccvo, nuz * sizeof(*ccvo));
		zccv = realloc(zccv, nuz * sizeof(*zccv));
		zccvo = realloc(zccvo, nuz * sizeof(*zccvo));

		if (UNLIKELY(cc == NULL)) {
			return -1;
		} else if (UNLIKELY(nccv == NULL)) {
			return -1;
		} else if (UNLIKELY(ccv == NULL)) {
			return -1;
		} else if (UNLIKELY(ccvo == NULL)) {
			return -1;
		} else if (UNLIKELY(zccv == NULL)) {
			return -1;
		} else if (UNLIKELY(zccvo == NULL)) {
			return -1;
		}
		for (size_t j = zcc; j < nuz; j++) {
			ccvo[j] = calloc(8U, sizeof(*ccvo[j]));
			zccvo[j] = 8U;
		}

		zcc = nuz;
	}
	/* really add him now */
	cc[ncc++] = c;
	return 0;
}


static size_t
mvhs(size_t ncol)
{
/* one pass of bubble sort */
	size_t v;

	for (v = ncol; --v;) {
		if (!nlhs) {
			if (lhs.v == v) {
				goto next;
			}
		} else for (size_t i = 0U; i < nlhs; i++) {
			if (lhs.p[i] == v) {
				goto next;
			}
		}
		if (rhs == v) {
			goto next;
		}
		return v + 1U;
	next:
		continue;
	}
	return 0;
}

static int
chck(size_t ncol)
{
	if (!nlhs) {
		if (UNLIKELY(lhs.v >= ncol)) {
			return -1;
		}
	} else for (size_t i = 0U; i < nlhs; i++) {
		if (UNLIKELY(lhs.p[i] >= ncol)) {
			return -1;
		}
	}
	if (UNLIKELY(rhs >= ncol && rhs + 1U)) {
		return -1;
	}
	if (vhs && vhs > ncol) {
		return -1;
	}
	/* determine VHS as the rightmost column not used by LHS nor RHS */
	if (UNLIKELY(!vhs && !(vhs = mvhs(ncol)))) {
		return -1;
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
	/* last dimension hash */
	uint64_t last_d = 0ULL;

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
Error: less columns present than needed for LHS~RHS and value");
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

		/* hash dimension columns */
		with (uint64_t d = 0U) {
			if (!nlhs) {
				const size_t of = coff[lhs.v + 0U];
				const size_t eo = coff[lhs.v + 1U];
				d = hash(line + of, eo - of - 1U);
			} else for (size_t i = 0U; i < nlhs; i++) {
				const size_t of = coff[lhs.p[i] + 0U];
				const size_t eo = coff[lhs.p[i] + 1U];
				d ^= hash(line + of, eo - of - 1U);
			}
			if (UNLIKELY(d != last_d)) {
				prnt();
				rset(line, coff);
				last_d = d;
				/* materialise cast cols */
				if (UNLIKELY(zcc)) {
					mtcc();
				}
			}
		}

		/* store value */
		if (rhs < ncol) {
			size_t of = coff[rhs + 0U];
			size_t eo = coff[rhs + 1U];
			const uint64_t c = hash(line + of, eo - of - 1);
			size_t j;

			for (j = 0U; j < ncc; j++) {
				if (cc[j] == c) {
					/* found him */
					goto bang;
				}
			}

			if (UNLIKELY(zcc || !ncc)) {
				/* add him */
				if (UNLIKELY(adcc(c) < 0)) {
					break;
				}
			} else {
				/* column we didn't want */
				continue;
			}
		bang:
			/* bang */
			of = coff[vhs - 1U];
			eo = coff[vhs];
			bang(line + of, eo - of - 1U, j);
		}
	}
	/* print the last one */
	prnt();

	free(coff);
out:
	free(line);
	return rc;
}

static int
snrf(char *formula)
{
	size_t zlhs = 0U;
	long unsigned int x;
	char *on = formula;

redo:
	x = strtoul(on, &on, 10);
	/* skip white space */
	for (; (unsigned char)(*on - 1) < ' '; on++);
	switch (*on++) {
	case '+':
		/* more to come */
		if (UNLIKELY(nlhs >= zlhs)) {
			zlhs = (zlhs * 2U) ?: 8U;
			lhs.p = calloc(zlhs, sizeof(*lhs.p));
		}
		lhs.p[nlhs++] = x - (x > 0U);
		goto redo;
	case '~':
		/* rhs coming */
		if (!nlhs) {
			lhs.v = x - (x > 0U);
		} else {
			if (UNLIKELY(nlhs >= zlhs)) {
				zlhs = (zlhs * 2U) ?: 8U;
				lhs.p = calloc(zlhs, sizeof(*lhs.p));
			}
			lhs.p[nlhs++] = x - (x > 0U);
		}
		goto rhs;
	default:
		return -1;
	}
rhs:
	if (UNLIKELY(!(rhs = strtoul(on, &on, 10)))) {
		/* skip white space */
		for (; (unsigned char)(*on - 1) < ' '; on++);
		if (*on == '.') {
			/* bla ~ . */
			rhs--;
			return 0;
		}
		return -1;
	} else if (UNLIKELY(*on)) {
		return -1;
	}
	/* we're 0-based internally */
	rhs--;

	/* check that LHS and RHS are disjoint */
	if (!nlhs && UNLIKELY(lhs.v == rhs)) {
		return -1;
	} else for (size_t i = 0U; i < nlhs; i++) {
		if (UNLIKELY(lhs.p[i] == rhs)) {
			return -1;
		}
	}
	return 0;
}


#include "dtcast.yucc"

int
main(int argc, char *argv[])
{
	static yuck_t argi[1U];
	int rc = 0;

	if (yuck_parse(argi, argc, argv) < 0) {
		rc = 1;
		goto out;
	}

	/* snarf formula */
	if (UNLIKELY(!argi->nargs || snrf(*argi->args) < 0)) {
		error("\
Error: cannot interpret formula");
		rc = 1;
		goto out;
	}

	if (argi->value_arg) {
		char *on;
		if (!(vhs = strtoul(argi->value_arg, &on, 10)) || *on) {
			errno = 0, error("\
Error: invalide value column");
			rc = 1;
			goto out;
		}
	}

	if (!(rhs + 1U)) {
		/* don't set up cast columns regardless what they specified */
		;
	} else if (UNLIKELY(stcc(argi->cast_args, argi->cast_nargs) < 0)) {
		error("\
Error: cannot set up cast columns");
		rc = 1;
		goto out;
	}

	rc = proc1() < 0;

	/* free cast columns */
	free(cc);
	free(nccv);
	free(ccv);
	free(zccv);
	free(zccvo);
	for (size_t i = 0U; i < ncc; i++) {
		free(ccvo[i]);
	}
	free(ccvo);
	free(dim);

	if (nlhs) {
		free(lhs.p);
	}

out:
	yuck_free(argi);
	return rc;
}

/* dtcast.c ends here */
