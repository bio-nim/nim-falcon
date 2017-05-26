default: run-falcon
update:
	git submodule update --init
run-%: %.exe
	./$<
%.exe: %.nim
	nim c --out:$@ $<
clean:
	rm -rf *.exe nimcache/
#.PRECIOUS: %.exe
