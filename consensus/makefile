NIMFLAGS+=--threads:on
#NIMFLAGS+=--threadAnalysis:off
NIMFLAGS+=--verbosity:2
NIMFLAGS+=-d:debug
#NIMFLAGS+=-d:release
NIMFLAGS+=--tlsemulation:on

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
CLEAN_NIMCACHE=rm -f nimcache/*.o nimcache/*.sha1
else
CLEAN_NIMCACHE=rm -rf nimcache/
endif

clean:
	${CLEAN_NIMCACHE}
	rm -f main.exe
