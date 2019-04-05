N=${CURDIR}/nimbleDir

foo:
	cd repos/nim-DALIGNER; nimble install --verbose --nimbleDir:$N -y
install-all:
	cd repos/nim-heap; nimble install --verbose --nimbleDir:$N -y  # binaryheap
	cd repos/msgpack4nim; nimble install --verbose --nimbleDir:$N -y
	cd repos/cligen; nimble install --verbose --nimbleDir:$N -y
	#cd repos/nim-DAZZ_DB; nimble install --verbose --nimbleDir:$N -y
	cd repos/nim-DALIGNER; nimble install --verbose --nimbleDir:$N -y
	cd repos/nim-htslib; nimble install --verbose --nimbleDir:$N -y
	cd repos/hts-nim; nimble install --verbose --nimbleDir:$N -y
default: run-falcon
update:
	git submodule update --init
run-%: %.exe
	./$<
%.exe: %.nim
	nim c -p:$N/pkgs/daligner-0.0.0 --out:$@ $<
clean:
	rm -rf *.exe nimcache/
#.PRECIOUS: %.exe
