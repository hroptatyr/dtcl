/*** dtchanges.c -- streaming implementation of changes.R
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
#include "coru.h"
#include "nifty.h"

struct hs_s {
	size_t n;
	size_t *c;
	size_t *p;
};

static int hdrp = 0;
static int cnmp = 0;
static int brfp = 0;
static const char *form;

/* join columns in left and right file */
static struct hs_s jc[2U];
static struct hs_s xc[2U];
static struct hs_s vc[2U];
#define L	0U
#define R	1U

typedef enum {
	EQU,
	DEL,
	ADD,
	CHG,
	NCHGTYP,
} chgtyp_t;

/* union header */
char *hdr;
size_t nhdr;
size_t zhdr;
size_t *hof;
size_t nhof;
size_t zhof;


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

static inline const char*
eatws(const char *s)
{
	for (; (unsigned char)(*s - 1) < ' '; s++);
	return s;
}

static inline unsigned int
streqp(const char *x, size_t m, const char *y, size_t n)
{
	return m == n && !memcmp(x, y, n);
}


static int
chck(struct hs_s *tg, const struct hs_s *sj, const struct hs_s *sx, size_t ncol)
{
	for (size_t i = 0U; i < sj->n; i++) {
		if (UNLIKELY(sj->c[i] >= ncol)) {
			return -1;
		}
	}
	for (size_t i = 0U; i < sx->n; i++) {
		if (UNLIKELY(sx->c[i] >= ncol)) {
			return -1;
		}
	}
	if (UNLIKELY(sj->n + sx->n > ncol)) {
		return -1;
	}

	/* construct value side, quick bubble sort */
	with (const size_t nv = ncol - sj->n - sx->n) {
		if (!nv) {
			break;
		}

		tg->c = malloc(nv * sizeof(*tg->c));
		tg->n = nv;
		for (size_t i = 0U, k = 0U; i < ncol; i++) {
			for (size_t j = 0U; j < sj->n; j++) {
				if (sj->c[j] == i) {
					goto next;
				}
			}
			for (size_t j = 0U; j < sx->n; j++) {
				if (sx->c[j] == i) {
					goto next;
				}
			}
			tg->c[k++] = i;
		next:
			continue;
		}
	}
	return 0;
}

static ssize_t
addhdr(const char *s, size_t n)
{
	char pivot;

	if (UNLIKELY(!n)) {
		errno = 0, error("\
Error: zero length column name not allowed");
		return -1;
	}

	pivot = s[n - 1U];
	for (const char *hp = hdr, *const ep = hdr + nhdr, *x;
	     hp < ep && (x = memchr(hp, pivot, ep - hp)) != NULL; hp = x + 1U) {
		const char *b = x - (n - 1U);
		if (b > hdr && b[-1] == '\t' && b[n] == '\t' &&
		    !memcmp(b, s, n)) {
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
hdrs(struct hs_s *restrict x, const char *ln, const size_t *of, size_t nc)
{
	if (UNLIKELY(!x->n)) {
		return 0;
	}

	x->p = malloc(nc * sizeof(*x->p));
	memset(x->p, -1, nc * sizeof(*x->p));

	if (of) {
		for (size_t i = 0U; i < x->n; i++) {
			const size_t c = x->c[i];
			const char *s;
			size_t n;
			ssize_t k;

			s = ln + of[c];
			n = of[c + 1U] - of[c + 0U] - 1U;
			k = addhdr(s, n);

			if (UNLIKELY(k < 0)) {
				return -1;
			}
			x->p[c] = k;
		}
	} else {
		char tmp[32U];
		for (size_t i = 0U; i < x->n; i++) {
			const size_t c = x->c[i];
			ssize_t k;
			int m;

			m = snprintf(tmp, sizeof(tmp), "V%zu", c + 1U);
			if (UNLIKELY(m < 0)) {
				errno = 0, error("\
Error: cannot generate column name");
				return -1;
			}
			k = addhdr(tmp, m);

			if (UNLIKELY(k < 0)) {
				return -1;
			}
			x->p[c] = k;
		}
	}
	return 0;
}

static int
invperm(struct hs_s *restrict tg)
{
	size_t *invp = malloc(nhof * sizeof(*invp));

	memset(invp, -1, nhof * sizeof(*invp));
	for (size_t i = 0U; i < jc->n; i++) {
		invp[jc->p[jc->c[i]]] = jc->c[i];
	}
	for (size_t i = 0U; i < xc->n; i++) {
		invp[xc->p[xc->c[i]]] = xc->c[i];
	}
	for (size_t i = 0U; i < tg->n; i++) {
		invp[tg->p[tg->c[i]]] = tg->c[i];
	}

	free(tg->p);
	tg->p = invp;
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

static ssize_t
find_s(const char *ss, const size_t *of, size_t nc, const char *s, size_t z)
{
/* find S of size Z in {SS + OF} */
	/* eat whitespace */
	for (; z > 0 && (unsigned char)(*s - 1) < ' '; s++, z--);
	for (; z > 0 && (unsigned char)(s[z - 1U] - 1) < ' '; z--);
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
snrf(struct hs_s *restrict tg,
     const char *hn, const size_t *of, size_t nc, size_t i, size_t m)
{
/* snarf m-th side of FORM (global) into temporary TG based on header line
   HN and header OF over NC columns, I-th file, i.e. skip I `=' tokens
   in formula */
	const char *ej, *j;
	const char *on;
	size_t nj;

	j = form;
	ej = j + strlen(j);

	for (size_t k = 0U; j < ej && k < m;
	     k++, j = memchrnul(j, '~', ej - j) + 1U);
	if (j >= ej) {
		return !!m - 1;
	}
	/* set end of snarf string accordingly */
	ej = memchrnul(j, '~', ej - j);

	for (nj = 0U, on = j; (on = memchr(on, '+', ej - on)); nj++, on++);

	tg->c = calloc(tg->n = nj + 1U, sizeof(*tg->c));
	/* now try and snarf the whole shebang */
	for (nj = 0U; nj < tg->n; nj++, j = on + 1U) {
		const char *om;
		char *tmp;
		long unsigned int x;
		size_t k;

		on = memchrnul(j, '+', ej - j);
		k = 0U;
		do {
			om = memchrnul(j, '=', on - j);
		} while (om < on && k++ < i && (j = om + 1U, true));

		/* try with numbers first */
		if ((x = strtoul(j, &tmp, 10)) && eatws(tmp) == om) {
			x--;
		} else if ((x = find_s(hn, of, nc, j, om - j)) < nc) {
			;
		} else {
			/* retry next time */
			free(tg->c);
			tg->c = NULL;
			tg->n = 0U;
			return -1;
		}

		tg->c[nj] = x;
	}
	return 0;
}

static void
prnc(const char *base, const size_t *cols, size_t i)
{
	if (LIKELY(i < nhof)) {
		const size_t bo = cols[i + 0U];
		const size_t eo = cols[i + 1U];
		fwrite(base + bo, sizeof(*base), eo - bo - 1U, stdout);
	}
	return;
}


struct beef_s {
	ssize_t nrd;
	char *line;
	/* intra-line */
	size_t ncol;
	size_t *coff;
	/* line number */
	size_t nr;
	/* constant dimension line */
	size_t ndln;
	char *dln;
};

DEFCORU(co_proc1, {
		FILE *fp;
		size_t fibre;
	}, void *arg)
{
	FILE *const fp = CORU_CLOSUR(fp);
	size_t fibre = CORU_CLOSUR(fibre);
	size_t llen = 0U;
	struct beef_s b = {0U};
	/* constant dimension line */
	size_t zdln = 0U;
	int rc = 0;

	/* probe */
	if (UNLIKELY((b.nrd = getline(&b.line, &llen, fp)) < 0)) {
		error("\
Error: cannot read lines");
		rc = -1;
		goto out;
	} else if (UNLIKELY(b.nrd == 0)) {
		goto out;
	} else if (UNLIKELY(!(b.ncol = toklng(b.line, b.nrd)))) {
		errno = 0, error("\
Error: cannot determine number of columns");
		rc = -1;
		goto out;
	} else if (UNLIKELY(!(b.coff = calloc(b.ncol + 1U, sizeof(*b.coff))))) {
		error("\
Error: cannot allocate memory to hold one line");
		rc = -1;
		goto out;
	}
	/* tokenise once */
	tokln1(b.coff, b.ncol, b.line, b.nrd);

	/* we might need to rescan the formula now */
	if (UNLIKELY(snrf(&jc[fibre], b.line, b.coff, b.ncol, fibre, 0U)) < 0) {
		error("\
Error: cannot interpret formula");
		rc = -1;
		goto out;
	} else if (UNLIKELY(snrf(&xc[fibre], b.line, b.coff, b.ncol, fibre, 1U)) < 0) {
		error("\
Error: cannot interpret formula");
		rc = -1;
		goto out;
	} else if (UNLIKELY(chck(&vc[fibre], &jc[fibre], &xc[fibre], b.ncol) < 0)) {
		errno = 0, error("\
Error: fewer columns present than needed for formula");
		rc = -1;
		goto out;
	}

	/* record header line and rbind it */
	{
		const char *ln = hdrp ? b.line : NULL;
		size_t *of = hdrp ? b.coff : NULL;

		if (!fibre && UNLIKELY(hdrs(jc, ln, of, b.ncol) < 0)) {
			rc = -1;
			goto out;
		}
		if (!fibre && UNLIKELY(hdrs(xc, ln, of, b.ncol) < 0)) {
			rc = -1;
			goto out;
		}
		if (UNLIKELY(hdrs(&vc[fibre], ln, of, b.ncol) < 0)) {
			rc = -1;
			goto out;
		}
	}

	if (!hdrp) {
		goto tok;
	}

	while ((b.nrd = getline(&b.line, &llen, fp)) > 0) {
		size_t bo, eo;
	tok:
		b.nr++;
		size_t nf = tokln1(b.coff, b.ncol, b.line, b.nrd);

		if (UNLIKELY(nf < b.ncol)) {
			errno = 0, error("\
Error: line %zu has only %zu columns, expected %zu", b.nr, nf, b.ncol);
			rc = -1;
			break;
		}

		/* construct constant dimension prefix */
		b.ndln = 0U;
		for (size_t i = 0U; i < jc[L].n; i++) {
			bo = b.coff[jc[fibre].c[i] + 0U];
			eo = b.coff[jc[fibre].c[i] + 1U];

			if (UNLIKELY(b.ndln + eo - bo > zdln)) {
				/* resize */
				while ((zdln = (zdln * 2U) ?: 256U) <=
				       b.ndln + eo - bo);
				b.dln = realloc(b.dln, zdln * sizeof(*b.dln));
			}
			memcpy(b.dln + b.ndln, b.line + bo, eo - bo - 1);
			b.ndln += eo - bo - 1;
			b.dln[b.ndln++] = '\t';
		}
		/* terminate dln */
		b.ndln -= !!jc[L].n;
		b.dln[b.ndln] = '\0';

		/* prep yield */
		*(struct beef_s*)arg = b;
		if (YIELD(1) < 0) {
			break;
		}
	}
out:
	if (jc[fibre].n) {
		free(jc[fibre].c);
		free(jc[fibre].p);
	}
	if (xc[fibre].n) {
		free(xc[fibre].c);
		free(xc[fibre].p);
	}
	if (vc[fibre].n) {
		free(vc[fibre].c);
		free(vc[fibre].p);
	}
	free(b.coff);
	free(b.dln);
	free(b.line);
	return rc;
}

static void
prnt(const struct beef_s *x, const struct beef_s *y)
{
	if (x && y) {
		/* compare cols
		 * "" ~ "SOMETHING" -> "+SOMETHING"
		 * "SOMETHING" ~ "" -> "-SOMETHING"
		 * "SOME" ~ "THING" -> "SOME => THING" */
		uint_fast8_t z[nhof];
		memset(z, 0, sizeof(z));
		for (size_t i = jc->n; i < countof(z); i++) {
			size_t cl = vc[L].p[i];
			size_t cr = vc[R].p[i];
#define na(z, w)	(w > (z)->ncol || (z)->coff[w] + 1U == (z)->coff[w + 1])
#define eq(l, r)	streqp(x->line + x->coff[l],			\
			       x->coff[l + 1] - (x->coff[l] + 1),	\
			       y->line + y->coff[r],			\
			       y->coff[r + 1] - (y->coff[r] + 1))
			uint_fast8_t s = (uint8_t)(na(x, cl) << 1U ^ na(y, cr));
			uint_fast8_t t = (uint8_t)(!s && !eq(cl, cr));
			/* two NAs is not considered a change */
			uint_fast8_t u = (uint8_t)(s & 0b1U ^ (s >> 1U) & 0b1U);

			/* massage s */
			s &= (uint_fast8_t)(u ^ u << 1U);
			z[i] = (uint_fast8_t)(s ^ t ^ t << 1U);
		}
		for (size_t j = jc->n + xc->n; j < countof(z); j++) {
			if (z[j]) {
				goto pr;
			}
		}
		return;
	pr:
		fputc(' ', stdout);
		fwrite(y->dln, 1, y->ndln, stdout);
		for (size_t i = jc->n; i < nhof; i++) {
			fputc('\t', stdout);
			switch (z[i]) {
			case 0U:
			default:
				continue;
			case 1U:
				fputc('-', stdout);
				prnc(x->line, x->coff, vc[L].p[i]);
				continue;
			case 2U:
				fputc('+', stdout);
				prnc(y->line, y->coff, vc[R].p[i]);
				continue;
			case 3U:
				break;
			}
			prnc(x->line, x->coff, vc[L].p[i]);
			fwrite(" => ", 1, 4U, stdout);
			prnc(y->line, y->coff, vc[R].p[i]);
		}
	} else if (x) {
		fputc('-', stdout);
		fwrite(x->dln, 1, x->ndln, stdout);
		for (size_t i = jc[L].n; i < nhof; i++) {
			fputc('\t', stdout);
			prnc(x->line, x->coff, vc[L].p[i]);
		}
	} else if (y) {
		fputc('+', stdout);
		fwrite(y->dln, 1, y->ndln, stdout);
		for (size_t i = jc[R].n; i < nhof; i++) {
			fputc('\t', stdout);
			prnc(y->line, y->coff, vc[R].p[i]);
		}
	} else {
		return;
	}
	fputc('\n', stdout);
	return;
}

static int
proc(FILE *fpx, FILE *fpy)
{
/* coordinator between fpx and fpy */
	int rc = 0;
	struct cocore *self = PREP();
	struct cocore *px = START_PACK(co_proc1, .next = self,
				       .clo = {.fp = fpx, .fibre = 0U});
	struct cocore *py = START_PACK(co_proc1, .next = self,
				       .clo = {.fp = fpy, .fibre = 1U});
	struct beef_s bx;
	struct beef_s by;
	int sx = NEXT1(px, &bx);
	int sy = NEXT1(py, &by);

	if (cnmp && sx > 0 && sy > 0) {
		hdr[nhdr - 1U] = '\n';
		fwrite(hdr + 1U, sizeof(*hdr), nhdr - 1U, stdout);
	}

	if (sx > 0 && sy > 0) {
		invperm(&vc[L]);
		invperm(&vc[R]);
	}

	for (int c; sx > 0 || sy > 0;
	     sx = NEXT1(px, &bx), sy = NEXT1(py, &by)) {
		if (sx > 0 && sy > 0) {
		redo:
			c = strcmp(bx.dln, by.dln);

			if (0) {
				;
			} else if (c < 0) {
				/* bx first, then by */
				prnt(&bx, NULL);

				if ((sx = NEXT1(px, &bx)) <= 0) {
					/* short circuit to by
					 * sy was guaranteed to be > 0 */
					goto rest_y;
				}
				goto redo;
			} else if (c > 0) {
				/* bx first, then by */
				prnt(NULL, &by);

				if ((sy = NEXT1(py, &by)) <= 0) {
					/* short-circuit to bx
					 * sx was guaranteed to be > 0 */
					goto rest_x;
				}
				goto redo;
			} else {
				/* keys are equal do a col-by-col comparison */
				prnt(&bx, &by);
			}
		} else if (sx > 0) {
		rest_x:
			/* we're out of BYs */
			do {
				prnt(&bx, NULL);
			} while ((sx = NEXT1(px, &bx)) > 0);
			break;
		} else if (sy > 0) {
		rest_y:
			/* we're out of BXs */
			do {
				prnt(NULL, &by);
			} while ((sy = NEXT1(py, &by)) > 0);
			break;
		}
	}
	UNPREP();
	return rc;
}

static unsigned int
csum(size_t sum[static NCHGTYP], const struct beef_s *x, const struct beef_s *y)
{
	/* compare cols
	 * "" ~ "SOMETHING" -> "+SOMETHING"
	 * "SOMETHING" ~ "" -> "-SOMETHING"
	 * "SOME" ~ "THING" -> "SOME => THING" */
	uint_fast8_t z[nhof];
	memset(z, 0, sizeof(z));
	for (size_t i = jc->n; i < countof(z); i++) {
		size_t cl = vc[L].p[i];
		size_t cr = vc[R].p[i];
#define na(z, w)	(w > (z)->ncol || (z)->coff[w] + 1U == (z)->coff[w + 1])
#define eq(l, r)	streqp(x->line + x->coff[l],			\
			       x->coff[l + 1] - (x->coff[l] + 1),	\
			       y->line + y->coff[r],			\
			       y->coff[r + 1] - (y->coff[r] + 1))
		uint_fast8_t s = (uint8_t)(na(x, cl) << 1U ^ na(y, cr));
		uint_fast8_t t = (uint8_t)(!s && !eq(cl, cr));
		/* two NAs is not considered a change */
		uint_fast8_t u = (uint8_t)(s & 0b1U ^ (s >> 1U) & 0b1U);

		/* massage s */
		s &= (uint_fast8_t)(u ^ u << 1U);
		z[i] = (uint_fast8_t)(s ^ t ^ t << 1U);
	}
	for (size_t j = jc->n + xc->n; j < countof(z); j++) {
		if (z[j]) {
			goto pr;
		}
	}
	return EQU;
pr:
	for (size_t i = jc->n; i < nhof; i++) {
		sum[z[i]]++;
	}
	return CHG;
}

static int
summ(FILE *fpx, FILE *fpy)
{
/* coordinator between fpx and fpy */
	struct cocore *self = PREP();
	struct cocore *px = START_PACK(co_proc1, .next = self,
				       .clo = {.fp = fpx, .fibre = 0U});
	struct cocore *py = START_PACK(co_proc1, .next = self,
				       .clo = {.fp = fpy, .fibre = 1U});
	struct beef_s bx;
	struct beef_s by;
	int sx = NEXT1(px, &bx);
	int sy = NEXT1(py, &by);
	size_t nl[NCHGTYP] = {0U};
	size_t nc[NCHGTYP] = {0U};
	int rc = 0;

	if (sx > 0 && sy > 0) {
		invperm(&vc[L]);
		invperm(&vc[R]);
	}

	for (int c; sx > 0 || sy > 0;
	     sx = NEXT1(px, &bx), sy = NEXT1(py, &by)) {
		if (sx > 0 && sy > 0) {
		redo:
			c = strcmp(bx.dln, by.dln);

			if (0) {
				;
			} else if (c < 0) {
				nl[DEL]++;

				if ((sx = NEXT1(px, &bx)) <= 0) {
					/* short circuit to by
					 * sy was guaranteed to be > 0 */
					goto rest_y;
				}
				goto redo;
			} else if (c > 0) {
				/* bx first, then by */
				nl[ADD]++;

				if ((sy = NEXT1(py, &by)) <= 0) {
					/* short-circuit to bx
					 * sx was guaranteed to be > 0 */
					goto rest_x;
				}
				goto redo;
			} else {
				/* keys are equal do a col-by-col comparison */
				nl[csum(nc, &bx, &by)]++;
			}
		} else if (sx > 0) {
		rest_x:
			/* we're out of BYs */
			do {
				nl[DEL]++;
			} while ((sx = NEXT1(px, &bx)) > 0);
			break;
		} else if (sy > 0) {
		rest_y:
			/* we're out of BXs */
			do {
				nl[ADD]++;
			} while ((sy = NEXT1(py, &by)) > 0);
			break;
		}
	}
	UNPREP();

	if (!brfp) {
		printf("%zu line(s) added\n", nl[ADD]);
		printf("%zu line(s) removed\n", nl[DEL]);
		printf("%zu line(s) changed\n", nl[CHG]);
		printf("  %zu value(s) added\n", nc[ADD]);
		printf("  %zu value(s) removed\n", nc[DEL]);
		printf("  %zu value(s) changed\n", nc[CHG]);
	} else {
		for (size_t i = 0U; i < NCHGTYP; i++) {
			fprintf(stdout, "%zu", nl[i]);
			fputc('\t', stdout);
		}
		for (size_t i = 0U; i < NCHGTYP; i++) {
			fprintf(stdout, "%zu", nc[i]);
			fputc('\t' + (i == CHG), stdout);
		}
	}
	return rc;
}


#include "dtchanges.yucc"

int
main(int argc, char *argv[])
{
	static yuck_t argi[1U];
	static FILE *fpx;
	static FILE *fpy;
	int rc = 0;

	if (yuck_parse(argi, argc, argv) < 0) {
		rc = 1;
		goto out;
	}

	/* overread and/or expect headers? */
	hdrp = argi->header_flag;
	/* memorise that we want col names for STCC() later on */
	cnmp = argi->col_names_flag;

	if (argi->nargs < 3U) {
		errno = 0, error("\
Error: need two files and a formula");
		rc = 1;
		goto out;
	} else if (UNLIKELY(!(fpx = fopen(argi->args[0U], "r")))) {
		error("\
Error: cannot open `%s' for reading", argi->args[0U]);
		rc = 1;
		goto clo;
	} else if (UNLIKELY(!(fpy = fopen(argi->args[1U], "r")))) {
		error("\
Error: cannot open `%s' for reading", argi->args[1U]);
		rc = 1;
		goto clo;
	}
	/* keep track of the formula */
	form = argi->args[2U];

	/* prealloc some header space */
	if (UNLIKELY((hdr = malloc(zhdr = 256U)) == NULL)) {
		error("\
Error: cannot allocate space for header");
		rc = 1;
		goto clo;
	}
	/* start with a framing character */
	hdr[nhdr++] = '\t';
	/* offsets for the header string */
	if (UNLIKELY((hof = malloc((zhof = 32U) * sizeof(*hof))) == NULL)) {
		error("\
Error: cannot allocate space for header");
		rc = 1;
		goto out;
	}

	/* get the coroutines going */
	initialise_cocore();

	if (!argi->summary_arg) {
		rc = proc(fpx, fpy) < 0;
	} else {
		brfp = argi->summary_arg != YUCK_OPTARG_NONE;
		rc = summ(fpx, fpy) < 0;
	}

	free(hdr);
	free(hof);

clo:
	if (fpx) {
		fclose(fpx);
	}
	if (fpy) {
		fclose(fpy);
	}
out:
	yuck_free(argi);
	return rc;
}

/* dtchange.c ends here */
