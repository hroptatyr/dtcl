/*** dtmerge.c -- streaming implementation of data.table's merge
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
	union {
		size_t v;
		size_t *p;
	};
};


static int hdrp = 0;
static int cnmp = 0;
static int allx;
static int ally;
static const char *form;

/* join columns in left and right file */
static struct hs_s jc[2U] = {{0, -1}, {0, -1}};
static struct hs_s vc[2U] = {{0, -1}, {0, -1}};
#define L	0U
#define R	1U

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

static inline const char*
memchrnul(const char *s, int c, size_t z)
{
	return memchr(s, c, z) ?: s + z;
}


static int
chck(struct hs_s *tg, const struct hs_s *sr, size_t ncol)
{
	if (!sr->n) {
		if (UNLIKELY(sr->v >= ncol)) {
			return -1;
		}
	} else for (size_t i = 0U; i < sr->n; i++) {
		if (UNLIKELY(sr->p[i] >= ncol)) {
			return -1;
		}
	}

	/* construct value side, quick bubble sort */
	with (const size_t nv = ncol - sr->n - !sr->n) {
		if (!nv) {
			break;
		}

		tg->p = malloc(nv * sizeof(*tg->p));
		tg->n = nv;
		if (!sr->n) {
			size_t k = 0U;
			for (size_t i = 0U; i < sr->v; i++) {
				tg->p[k++] = i;
			}
			for (size_t i = sr->v + 1U; i < ncol; i++) {
				tg->p[k++] = i;
			}
		} else {
			for (size_t i = 0U, k = 0U; i < ncol; i++) {
				for (size_t j = 0U; j < sr->n; j++) {
					if (sr->p[j] == i) {
						goto next;
					}
				}
				tg->p[k++] = i;
			next:
				continue;
			}
		}
	}
	return 0;
}

static void
phdr(const char *hdrs, const size_t *hoff)
{
	size_t i;

	if (!jc[L].n) {
		i = jc[L].v;
		goto onl;
	}
	for (size_t j = 0U; j < jc[L].n; j++) {
		i = jc[L].p[j];
	onl:;
		const size_t of = hoff[i + 0U];
		const size_t eo = hoff[i + 1U];
		fwrite(hdrs + of, sizeof(*hdrs), eo - of - 1U, stdout);
		fputc('\t', stdout);
	}
	if (!jc[R].n) {
		i = jc[R].v;
		goto onr;
	}
	for (size_t j = 0U; j < jc[R].n; j++) {
		i = jc[R].p[j];
	onr:;
		const size_t of = hoff[i + 0U];
		const size_t eo = hoff[i + 1U];
		fwrite(hdrs + of, sizeof(*hdrs), eo - of - 1U, stdout);
		fputc('\t', stdout);
	}
	fputs("value\n", stdout);
	return;
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
snrf(struct hs_s *restrict tg, const char *hn, const size_t *of, size_t nc, size_t i)
{
/* snarf FORM (global) into temporary TG based on header line HN and header OF
   over NC columns, I-th file, i.e. skip I = tokens in formula */
	const char *ej, *j;
	const char *on;
	size_t nj;

	j = form;
	ej = j + strlen(j);

	for (nj = 0U, on = j; (on = memchr(on, '+', ej - on)); nj++, on++);

	if (nj) {
		tg->p = calloc(tg->n = nj + 1U, sizeof(*tg->p));
	} else {
		tg->v = -1;
		tg->n = 0U;
		goto one_j;
	}
	/* now try and snarf the whole shebang */
	for (nj = 0U; nj < tg->n; nj++, j = on + 1U) {
		const char *om;
		char *tmp;
		long unsigned int x;
		size_t k;

	one_j:
		on = memchrnul(j, '+', ej - j);
		k = 0U;
		do {
			om = memchrnul(j, '=', on - j);
		} while (om < on && k++ < i && (j = om + 1U, true));

		/* try with numbers first */
		if ((x = strtoul(j, &tmp, 10)) && tmp == om) {
			x--;
		} else if ((x = find_s(hn, of, nc, j, om - j)) < nc) {
			;
		} else {
			/* retry next time */
			if (tg->v + 1U) {
				free(tg->p);
			}
			tg->v = -1;
			tg->n = 0U;
			return -1;
		}

		if (!tg->n) {
			tg->v = x;
		} else {
			tg->p[nj] = x;
		}
	}
	return 0;
}

static void
prnc(const char *base, const size_t *cols, size_t i)
{
	const size_t bo = cols[i + 0U];
	const size_t eo = cols[i + 1U];
	fwrite(base + bo, sizeof(*base), eo - bo - 1U, stdout);
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
	/* offsets for header and header buffer */
	char *hn = NULL;
	size_t *hoff = NULL;
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
	} else if (UNLIKELY(!(hoff = calloc(b.ncol + 1U, sizeof(*hoff))))) {
		error("\
Error: cannot allocate memory to hold a copy of the header");
		rc = -1;
		goto out;
	}

	if (hdrp) {
		/* snarf col names as defined in header */
		if (UNLIKELY((hn = strndup(b.line, b.nrd)) == NULL)) {
			error("\
Error: cannot allocate memory to hold a copy of the header");
			rc = -1;
			goto out;
		}
		tokln1(hoff, b.ncol, b.line, b.nrd);
		for (size_t i = 1U; i <= b.ncol; i++) {
			hn[hoff[i] - 1U] = '\0';
		}
	}
	/* we might need to rescan the formula now */
	if (UNLIKELY(snrf(&jc[fibre], hn, hoff, b.ncol, fibre)) < 0) {
		error("\
Error: cannot interpret formula");
		rc = -1;
		goto out;
	} else if (UNLIKELY(chck(&vc[fibre], &jc[fibre], b.ncol) < 0)) {
		errno = 0, error("\
Error: fewer columns present than needed for formula");
		rc = -1;
		goto out;
	}

	if (!hdrp) {
		if (UNLIKELY((hn = mkhdrs(hoff, b.ncol)) == NULL)) {
			error("\
Error: cannot allocate memory to hold a copy of the header");
			rc = -1;
			goto out;
		}
		if (cnmp) {
			/* print col names */
			phdr(hn, hoff);
		}
		goto tok;
	}

	if (cnmp && (b.nrd = getline(&b.line, &llen, fp)) > 0) {
		/* print col names */
		phdr(hn, hoff);
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
		if (!jc[fibre].n) {
			bo = b.coff[jc[fibre].v + 0U];
			eo = b.coff[jc[fibre].v + 1U];
			goto one_l;
		} else for (size_t i = 0U; i < jc[L].n; i++) {
			bo = b.coff[jc[fibre].p[i] + 0U];
			eo = b.coff[jc[fibre].p[i] + 1U];
		one_l:
			if (UNLIKELY(b.ndln + eo - bo > zdln)) {
				/* resize */
				while ((zdln = (zdln * 2U) ?: 256U) <=
				       b.ndln + eo - bo);
				b.dln = realloc(b.dln, zdln * sizeof(*b.dln));
			}
			memcpy(b.dln + b.ndln, b.line + bo, eo - bo - 1);
			b.ndln += eo - bo - 1;
			b.dln[b.ndln] = '\t';
		}
		/* terminate dln */
		b.dln[b.ndln] = '\0';

		/* prep yield */
		*(struct beef_s*)arg = b;
		if (YIELD(1) < 0) {
			break;
		}
	}
out:
	if (jc[fibre].n) {
		free(jc[fibre].p);
	}
	if (vc[fibre].n) {
		free(vc[fibre].p);
	}
	free(b.coff);
	free(b.dln);
	free(hoff);
	free(hn);
	free(b.line);
	return rc;
}

static void
prnt(const struct beef_s *x, const struct beef_s *y)
{
	static const char tabs[] = "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";

	if (x && y) {
		fwrite(y->dln, 1, y->ndln, stdout);
	} else if (x && allx) {
		fwrite(x->dln, 1, x->ndln, stdout);
	} else if (y && ally) {
		fwrite(y->dln, 1, y->ndln, stdout);
	} else {
		return;
	}
	if (x) {
		for (size_t i = 0U; i < vc[L].n; i++) {
			fputc('\t', stdout);
			prnc(x->line, x->coff, vc[L].p[i]);
		}
	} else if (ally) {
		for (size_t i = 0U; i < vc[L].n / 16U; i++) {
			fwrite(tabs, sizeof(*tabs), countof(tabs), stdout);
		}
		fwrite(tabs, sizeof(*tabs), vc[L].n % 16U, stdout);
	}
	if (y) {
		for (size_t i = 0U; i < vc[R].n; i++) {
			fputc('\t', stdout);
			prnc(y->line, y->coff, vc[R].p[i]);
		}
	} else if (allx) {
		for (size_t i = 0U; i < vc[R].n / 16U; i++) {
			fwrite(tabs, sizeof(*tabs), countof(tabs), stdout);
		}
		fwrite(tabs, sizeof(*tabs), vc[R].n % 16U, stdout);
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

	for (; sx > 0 || sy > 0;
	     sx = NEXT1(px, &bx), sy = NEXT1(py, &by)) {
		if (sx > 0 && sy > 0) {
			int c = strcmp(bx.dln, by.dln);

			if (0) {
				;
			} else if (c < 0) {
				/* bx first, then by */
				do {
					prnt(&bx, NULL);
				} while ((sx = NEXT1(px, &bx)) > 0 &&
					 strcmp(bx.dln, by.dln) < 0);
				if (sx <= 0) {
					/* short circuit to by
					 * sy was guaranteed to be > 0 */
					goto rest_y;
				}
			} else if (c > 0) {
				/* bx first, then by */
				do {
					prnt(NULL, &by);
				} while ((sy = NEXT1(py, &by)) > 0 &&
					 strcmp(bx.dln, by.dln) > 0);
				if (sy <= 0) {
					/* short-circuit to bx
					 * sx was guaranteed to be > 0 */
					goto rest_x;
				}
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
				prnt(&by, NULL);
			} while ((sy = NEXT1(py, &by)) > 0);
			break;
		}
	}
	UNPREP();
	return rc;
}


#include "dtmerge.yucc"

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

	if (argi->all_arg) {
		allx = argi->all_arg == YUCK_OPTARG_NONE ||
			*argi->all_arg == 'l' ||
			*argi->all_arg == 'x';
		ally = argi->all_arg == YUCK_OPTARG_NONE ||
			*argi->all_arg == 'r' ||
			*argi->all_arg == 'y';
	}

	/* get the coroutines going */
	initialise_cocore();

	/* prealloc some header space */
	if (UNLIKELY((hdr = malloc(zhdr = 256U)) == NULL)) {
		error("\
Error: cannot allocate space for header");
		rc = 1;
		goto clo;
	}
	/* start with a framing character */
	hdr[nhdr++] = '\t';

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
		goto clo;
	}

	rc = proc(fpx, fpy) < 0;

	free(hdr);
	free(hof);
	free(perm);

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

/* dtmerge.c ends here */
