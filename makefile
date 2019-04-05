PREFIX?=${CURDIR}
NIMBLE_DIR?=${CURDIR}/nimbleDir
export NIMBLE_DIR
# or use --nimbleDir:${NIMBLE_DIR} everywhere
NIMBLE_INSTALL=nimble install --debug -y

default:
	$N/bin/example
	./link-exe.sh bin/ # from ${NIMBLE_DIR}
	./bin/example.exe
	rm -rf bin/
all:
	${MAKE} update
	${MAKE} install
	${MAKE} link
update:
	git submodule update --init --recursive
install:
	cd repos/nim-heap; ${NIMBLE_INSTALL}
	cd repos/msgpack4nim; ${NIMBLE_INSTALL}
	cd repos/cligen; ${NIMBLE_INSTALL}
	#cd repos/nim-DAZZ_DB; ${NIMBLE_INSTALL}
	cd repos/nim-DALIGNER; ${NIMBLE_INSTALL}
	cd repos/nim-htslib; ${NIMBLE_INSTALL}
	#cd repos/hts-nim; ${NIMBLE_INSTALL}
	${NIMBLE_INSTALL}
link:
	./link-exe.sh ${PREFIX}
old-default: run-falcon
run-%: %.exe
	./$<
%.exe: %.nim
	nim c --out:$@ $<
clean:
	rm -rf *.exe nimcache/
#.PRECIOUS: %.exe

remotes:
	-git -C repos/nim-heap remote add up git@github.com:bluenote10/nim-heap
	-git -C repos/msgpack4nim remote add up git@github.com:jangko/msgpack4nim
	-git -C repos/cligen remote add up git@github.com:c-blake/cligen
	-git -C repos/hts-nim remote add up git@github.com:brentp/hts-nim
