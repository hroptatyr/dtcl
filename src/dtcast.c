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
static int cnmp = 0;
/* number of cast columns */
static size_t ncc;
static size_t zcc;
/* cast column hash values */
static uint64_t *cc;
static const char **cn;
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
} lhs;
static size_t nrhs;
static union {
	size_t v;
	size_t *p;
} rhs;
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

static inline const char*
memchrnul(const char *s, int c, size_t z)
{
	return memchr(s, c, z) ?: s + z;
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
	if (!(rhs.v + 1U)) {
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
	fwrite(dim, sizeof(*dim), ndim, stdout);
	for (size_t j = 0U; j < ncc; j++) {
		fputc('\t', stdout);
		if (nccv[j]) {
			fwrite(ccv[j] + ccvo[j][m[j]], sizeof(char),
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
phdr(const char *hdrs, const size_t *hoff)
{
	if (hdrs == NULL) {
		size_t i;
		size_t j = 0U;

		if (!nlhs) {
			i = lhs.v;
			goto one;
		}
		while (j < nlhs) {
			i = lhs.p[j];
		one:
			fputc('V', stdout);
			fprintf(stdout, "%zu", i + 1U);
			fputc('\t' + (++j >= nlhs && !ncc), stdout);
		}
	} else {
		size_t i;
		size_t j = 0U;

		if (!nlhs) {
			i = lhs.v;
			goto onh;
		}
		while (j < nlhs) {
			i = lhs.p[j];
		onh:;
			const size_t of = hoff[i + 0U];
			const size_t eo = hoff[i + 1U];
			fwrite(hdrs + of, sizeof(*hdrs), eo - of - 1U, stdout);
			fputc('\t' + (++j >= nlhs && !ncc), stdout);
		}
	}
	for (size_t i = 0U; i < ncc; i++) {
		fputs(cn[i], stdout);
		fputc('\t' + (i + 1 >= ncc), stdout);
	}
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
	for (size_t i = ncc; i < zcc; i++) {
		free(ccvo[i]);
	}
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
	} else if (cnmp && UNLIKELY((cn = calloc(ncc, sizeof(*cn))) == NULL)) {
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
	if (!cnmp) {
		return 0;
	}
	/* otherwise also store column names */
	for (size_t i = 0U; i < ncc; i++) {
		cn[i] = args[i];
	}
	return 0;
}

static int
adcc(uint64_t c, const char *n, size_t z)
{
	if (UNLIKELY(ncc >= zcc)) {
		const size_t nuz = (zcc * 2U) ?: 64U;

		cc = realloc(cc, nuz * sizeof(*cc));
		nccv = realloc(nccv, nuz * sizeof(*nccv));
		ccv = realloc(ccv, nuz * sizeof(*ccv));
		ccvo = realloc(ccvo, nuz * sizeof(*ccvo));
		zccv = realloc(zccv, nuz * sizeof(*zccv));
		zccvo = realloc(zccvo, nuz * sizeof(*zccvo));
		if (cnmp) {
			cn = realloc(cn, nuz * sizeof(*cn));
		}

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
		} else if (cnmp && UNLIKELY(cn == NULL)) {
			return -1;
		}
		memset(cc + zcc, 0, (nuz - zcc) * sizeof(*cc));
		memset(nccv + zcc, 0, (nuz - zcc) * sizeof(*nccv));
		memset(ccv + zcc, 0, (nuz - zcc) * sizeof(*ccv));
		memset(zccv + zcc, 0, (nuz - zcc) * sizeof(*zccv));
		for (size_t j = zcc; j < nuz; j++) {
			ccvo[j] = calloc(8U, sizeof(*ccvo[j]));
			zccvo[j] = 8U;
		}

		zcc = nuz;
	}
	/* really add him now */
	cc[ncc++] = c;
	if (!cnmp) {
		return 0;
	}
	/* otherwise also remember his name */
	cn[ncc - 1] = strndup(n, z);
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
		if (!nrhs) {
			if (rhs.v == v) {
				goto next;
			}
		} else for (size_t i = 0U; i < nrhs; i++) {
			if (rhs.p[i] == v) {
				goto next;
			}
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
	if (!nrhs) {
		if (UNLIKELY(rhs.v >= ncol && rhs.v + 1U)) {
			return -1;
		}
	} else for (size_t i = 0U; i < nrhs; i++) {
		if (UNLIKELY(rhs.p[i] >= ncol)) {
			return -1;
		}
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

static uint64_t
hashln(const char *ln, size_t *of)
{
	uint64_t d = 0U;

	if (!nlhs) {
		const size_t bo = of[lhs.v + 0U];
		const size_t eo = of[lhs.v + 1U];
		d = hash(ln + bo, eo - bo - 1U);
	} else for (size_t i = 0U; i < nlhs; i++) {
		const size_t bo = of[lhs.p[i] + 0U];
		const size_t eo = of[lhs.p[i] + 1U];
		d ^= hash(ln + bo, eo - bo - 1U);
	}
	return d;
}

static ssize_t
find(const char *ss, const size_t *of, size_t nc, const char *s, size_t z)
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
	if ((elhs = strchr(l = form = formula, '~')) == NULL) {
		return -1;
	}
	r = elhs + 1U;
	erhs = r + strlen(r);

	for (nl = 0U, on = l; (on = memchr(on, '+', elhs - on)); nl++, on++);
	for (nr = 0U, on = r; (on = memchr(on, '+', erhs - on)); nr++, on++);

	if (nl && !lhs.p) {
		lhs.p = calloc(nlhs = nl + 1U, sizeof(*lhs.p));
	}
	if (nr && !rhs.p) {
		rhs.p = calloc(nrhs = nr + 1U, sizeof(*rhs.p));
	}

	/* now try and snarf the whole shebang */
	if (!nl) {
		goto one_l;
	}
	for (nl = 0U; nl < nlhs; nl++, l = on + 1U) {
		/* try with numbers first */
		char *tmp;
		long unsigned int x;

	one_l:
		on = memchrnul(l, '+', elhs - l);
		if ((x = strtoul(l, &tmp, 10)) && tmp == on) {
			x--;
		} else if ((x = find(hn, of, nc, l, on - l)) < nc) {
			;
		} else {
			return -1;
		}

		if (!nlhs) {
			lhs.v = x;
		} else {
			lhs.p[nl] = x;
		}
	}

	/* snarf right hand side */
	if (!nr) {
		goto one_r;
	}
	for (nr = 0U; nr < nrhs; nr++, r = on + 1U) {
		/* try with numbers first */
		char *tmp;
		long unsigned int x;

	one_r:
		on = memchrnul(r, '+', erhs - r);
		if ((x = strtoul(r, &tmp, 10)) && tmp == on ||
		    *tmp == '.' && tmp + 1 == on && !nrhs) {
			x--;
		} else if ((x = find(hn, of, nc, r, on - r)) < nc) {
			;
		} else {
			return -1;
		}

		if (!nrhs) {
			rhs.v = x;
		} else {
			rhs.p[nr] = x;
		}
	}
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
	size_t *coff;
	/* offsets for header and header buffer */
	char *hn = NULL;
	size_t *hoff = NULL;
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
Error: fewer columns present than needed for LHS~RHS and value");
		rc = -1;
		goto out;
	} else if (UNLIKELY(!(coff = calloc(ncol + 1U, sizeof(*coff))))) {
		error("\
Error: cannot allocate memory to hold one line");
		rc = -1;
		goto out;
	}

	if (!hdrp && ncc) {
		/* go straight to tok loop */
		goto tok;
	} else if (!hdrp) {
		/* implies !ncc, we need to snarf cast cols then */
		goto scctok;
	}
	/* otherwise snarf col names as defined in header */
	if (UNLIKELY((hn = strndup(line, nrd)) == NULL ||
		     (hoff = calloc(ncol + 1U, sizeof(*hoff))) == NULL)) {
		error("\
Error: cannot allocate memory to hold a copy of the header");
		rc = -1;
		goto err;
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
		goto err;
	} else if (UNLIKELY(chck(ncol) < 0)) {
		errno = 0, error("\
Error: fewer columns present than needed for LHS~RHS and value");
		rc = -1;
		goto err;
	}

	/* depending on whether cast cols are specified explicitly */
	if (ncc) {
		/* yea we know what we want */
		goto tok;
	}

	/* snarf first complete group to obtain cast columns */
	while ((nrd = getline(&line, &llen, stdin)) > 0) {
	scctok:
		nr++;
		size_t nf = tokln1(coff, ncol, line, nrd);

		if (UNLIKELY(nf < ncol)) {
			errno = 0, error("\
Error: line %zu has only %zu columns, expected %zu", nr, nf, ncol);
			rc = 2;
			break;
		}

		/* hash dimension columns */
		with (const uint64_t d = hashln(line, coff)) {
			if (UNLIKELY(!last_d)) {
				rset(line, coff);
			} else if (UNLIKELY(d != last_d)) {
				/* materialise cast cols */
				mtcc();
				/* pretend we didn't see this line */
				nr--;
				/* and continue with main tokenisation loop */
				goto tok;
			}
			last_d = d;
		}

		/* store value */
		if (rhs.v < ncol) {
			size_t of = coff[rhs.v + 0U];
			size_t eo = coff[rhs.v + 1U];
			const uint64_t c = hash(line + of, eo - of - 1);
			size_t j;

			for (j = 0U; j < ncc; j++) {
				if (cc[j] == c) {
					/* found him */
					break;
				}
			}

			if (j >= ncc) {
				/* add him */
				const char *n = line + of;
				size_t z = eo - of - 1;
				if (UNLIKELY(adcc(c, n, z) < 0)) {
					break;
				}
			}

			/* bang */
			of = coff[vhs - 1U];
			eo = coff[vhs];
			bang(line + of, eo - of - 1U, j);
		}
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
		with (const uint64_t d = hashln(line, coff)) {
			if (UNLIKELY(!last_d)) {
				rset(line, coff);
			} else if (UNLIKELY(d != last_d)) {
				static int cprp;
				if (UNLIKELY(!cprp && cnmp)) {
					/* print col names */
					phdr(hn, hoff);
				}
				cprp = 1;
				prnt();
				rset(line, coff);
			}
			last_d = d;
		}

		/* store value */
		if (rhs.v < ncol) {
			size_t of = coff[rhs.v + 0U];
			size_t eo = coff[rhs.v + 1U];
			const uint64_t c = hash(line + of, eo - of - 1);
			size_t j;

			for (j = 0U; j < ncc; j++) {
				if (cc[j] == c) {
					/* found him */
					goto bang;
				}
			}
			/* column we didn't want */
			continue;

		bang:
			/* bang */
			of = coff[vhs - 1U];
			eo = coff[vhs];
			bang(line + of, eo - of - 1U, j);
		}
	}
	/* print the last one */
	prnt();

err:
	free(coff);
	free(hoff);
out:
	free(line);
	return rc;
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

	if (argi->value_arg) {
		char *on;
		if (!(vhs = strtoul(argi->value_arg, &on, 10)) || *on) {
			errno = 0, error("\
Error: invalide value column");
			rc = 1;
			goto out;
		}
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

	if (!(rhs.v + 1U)) {
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
	free(zccv);
	free(zccvo);
	for (size_t i = 0U; i < ncc; i++) {
		free(ccv[i]);
	}
	free(ccv);
	for (size_t i = 0U; i < ncc; i++) {
		free(ccvo[i]);
	}
	free(ccvo);
	free(dim);
	if (!argi->cast_nargs && cn) {
		for (size_t i = 0U; i < ncc; i++) {
			free(deconst(cn[i]));
		}
	}
	free(cn);

	if (nlhs) {
		free(lhs.p);
	}
	if (nrhs) {
		free(rhs.p);
	}

out:
	yuck_free(argi);
	return rc;
}

/* dtcast.c ends here */
