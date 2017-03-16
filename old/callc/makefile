NIMFLAGS=-d:debug 
NIMFLAGS=-d:release
NIMFLAGS+=--verbosity:2

do: n.exe
	# This test would require the t/ directory.
	make -C t # LA4Falcon | ./n.exe

foo.nim: foo.h
	c2nim foo.h --out:foo.nim
common.nim:
	c2nim --header --cdecl foo.h --out:foo.nim
run-%: %.exe
	./$*.exe
%.exe: %.nim
	nim ${NIMFLAGS} --out:$*.exe c $<

# Not needed, but possible.
LIBS:=DW_banded.o falcon.o kmer_lookup.o
lib: ${LIBS}
	${CC} -shared -o libcommon.so ${LIBS}

clean:
	rm -f n.exe
	rm -rf nimcache/
