/*** lines.c -- print certain lines from a bunch of files
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
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include "nifty.h"


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


#include "lines.yucc"

int
main(int argc, char *argv[])
{
	static yuck_t argi[1U];
	FILE **fp;
	char *line = NULL;
	size_t llen = 0U;
	size_t ln = 0U;
	int rc = 0;

	if (yuck_parse(argi, argc, argv) < 0) {
		rc = 1;
		goto out;
	} else if (!argi->nargs) {
		/* nothing to do */
		goto out;
	} else if (UNLIKELY((fp = calloc(argi->nargs, sizeof(*fp))) == NULL)) {
		error("\
Error: cannot allocate memory for file descriptors");
		rc = 1;
		goto out;
	}

	for (size_t i = 0U; i < argi->nargs; i++) {
		if (UNLIKELY((fp[i] = fopen(argi->args[i], "r")) == NULL)) {
			error("\
Warning: cannot open file `%s'", argi->args[i]);
		}
	}
	for (ssize_t nrd; (nrd = getline(&line, &llen, stdin)) > 0;) {
		char *on;
		size_t nx = strtoull(line, &on, 0);

		if (!nx || *on != '\n' || --nx < ln) {
			continue;
		}
		/* otherwise ffw to that line */
		for (; ln < nx; ln++) {
			for (size_t i = 0U; i < argi->nargs; i++) {
				if (LIKELY(fp[i] != NULL)) {
					getline(&line, &llen, fp[i]);
				}
			}
		}
		/* print him */
		for (size_t i = 0U; i < argi->nargs; i++) {
			if (UNLIKELY(fp[i] == NULL)) {
				continue;
			} else if ((nrd = getline(&line, &llen, fp[i])) <= 0) {
				continue;
			}
			nrd -= line[nrd - 1] == '\n';
			line[nrd++] = '\n';
			fwrite(line, 1, nrd, stdout);
		}
		ln++;
	}
	for (size_t i = 0U; i < argi->nargs; i++) {
		if (LIKELY(fp[i] != NULL)) {
			fclose(fp[i]);
		}
	}
out:
	yuck_free(argi);
	free(fp);
	free(line);
	return rc;
}

/* lines.c ends here */
