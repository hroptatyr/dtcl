data.table on the command-line
==============================

Collection of tools that mimic some of [data.table][1]'s operations
on the command-line with input and output from and to pipes.

Example?

    $ dtcast '1+2~3' <<EOF
    2009-03-12	AX	open	717.25
    2009-03-12	AX	high	718.47
    2009-03-12	AX	low	717.25
    2009-03-12	AX	close	718.42
    2009-03-12	BZX	open	715.32
    2009-03-12	BZX	high	717.57
    2009-03-12	BZX	low	714.65
    2009-03-12	BZX	close	718.35
    2009-03-13	AX	open	721.14
    2009-03-13	AX	high	721.24
    2009-03-13	AX	low	717.02
    2009-03-13	AX	close	717.14
    2009-03-13	BZX	open	717.34
    2009-03-13	BZX	high	719.26
    2009-03-13	BZX	low	717.34
    2009-03-13	BZX	close	718.00
    EOF
    =>
    2009-03-12	AX	717.25	718.47	717.25	718.42
    2009-03-12	BZX	715.32	717.57	714.65	718.35
    2009-03-13	AX	721.14	721.24	717.02	717.14
    2009-03-13	BZX	717.34	719.26	717.34	718.00

Oh no, we need to undo that:

    $ dtmelt -I 1 -I 2 <<EOF
    2009-03-12	AX	717.25	718.47	717.25	718.42
    2009-03-12	BZX	715.32	717.57	714.65	718.35
    2009-03-13	AX	721.14	721.24	717.02	717.14
    2009-03-13	BZX	717.34	719.26	717.34	718.00
    EOF
    =>
    2009-03-12	AX	3	717.25
    2009-03-12	AX	4	718.47
    2009-03-12	AX	5	717.25
    2009-03-12	AX	6	718.42
    2009-03-12	BZX	3	715.32
    2009-03-12	BZX	4	717.57
    2009-03-12	BZX	5	714.65
    2009-03-12	BZX	6	718.35
    2009-03-13	AX	3	721.14
    2009-03-13	AX	4	721.24
    2009-03-13	AX	5	717.02
    2009-03-13	AX	6	717.14
    2009-03-13	BZX	3	717.34
    2009-03-13	BZX	4	719.26
    2009-03-13	BZX	5	717.34
    2009-03-13	BZX	6	718.00


  [1]: http://github.com/Rdatatable/data.table/wiki
