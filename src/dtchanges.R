require("data.table");

if (!exists("changes")) {
	changes <- function(object1, object2, formula, ...) UseMethod("changes")
}

changes.data.table <- function(dt1, dt2, formula)
{
	if (is.character(formula)) {
		on = formula;
		and = NULL;
	} else if (length(formula) == 2) {
		on = as.character(formula[[2]]);
		and = NULL;
	} else {
		on = as.character(formula[[3]]);
		and = as.character(formula[[2]]);
	}
	if (is.vector(on) && on[1] == "+") {
		on = on[-1];
	}
	if (is.vector(and) && and[1] == "+") {
		and = and[-1];
	}

	setkeyv(dt1, c(on, and));
	setkeyv(dt2, c(on, and));

	m1 <- suppressWarnings(melt(dt1, c(on, and)))
	m2 <- suppressWarnings(melt(dt2, c(on, and)))

	setkeyv(m1, c(on, "variable"))
	setkeyv(m2, c(on, "variable"))
	m <- m1[m2, on=c(on, "variable"), nomatch=0L][value != i.value | xor(is.na(value), is.na(i.value))]
	m[is.na(value), value := paste0("+", i.value)]
	m[is.na(i.value), value := paste0("-", value)]
	m[!is.na(value) & !is.na(i.value), value := paste(value, i.value, sep=" => ")]

	for (j in and) {
		i.j <- paste("i", j, sep=".")
		eval(substitute(m[j != i.j, j := paste(j, i.j, sep=" => ")], list(j=as.symbol(j), i.j=as.symbol(i.j))))
	}

	## on + and ~ variable
	f <- as.formula(paste(paste0(c(on, and), collapse="+"), "variable", sep="~"))
	M <- dcast(m, f, value.var="value");

	X <- rbind(dt1[!dt2, on=on][, .m := "-"], dt2[!dt1, on=on][, .m := "+"], M[, .m := "/"], fill=T, use.names=T)
	setkeyv(X, c(on, and));
	return(X);
}
