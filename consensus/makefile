NIMFLAGS=-d:debug --debugger:native
NIMFLAGS=-d:release
NIMFLAGS+=--verbosity:2

do: main.exe
	make -C t # LA4Falcon | ./n.exe
main.exe: main.nim falcon.nim common.nim DW_banded.nim kmer_lookup_c.nim poo.nim
run-main:
convert: DW_banded.nim falcon.nim kmer_lookup.nim poo.nim common.nim
run-%: %.exe
	./$*.exe
%.exe: %.nim
	nim ${NIMFLAGS} --out:$*.exe c $<
%.nim: %.c
	c2nim -o:$@ $<
%.nim: %.h
	c2nim -o:$@ $<

ifeq (${USER},cdunn2001)
CLEAN_NIMCACHE=rm -f nimcache/*.o
else
CLEAN_NIMCACHE=rm -rf nimcache/
endif

clean:
	${CLEAN_NIMCACHE}
	rm -f main.exe
