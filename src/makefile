#NIMFLAGS+=--path:/home/UNIXHOME/cdunn/repo/gh/nim-multiproc/
NIMFLAGS+=--threads:on
NIMFLAGS+=--threadAnalysis:off
NIMFLAGS+=--verbosity:2
#NIMFLAGS+=-d:debug
NIMFLAGS+=-d:release
NIMFLAGS+=--tlsemulation:on
#NIMFLAGS+=-d:debugHeapLinks
#NIMFLAGS+=--profiler:on --stackTrace:on
#NIMFLAGS+=--passC:-Dexecvpe=execve

THISDIR=${CURDIR}

R=/home/UNIXHOME/cdunn/repo/bb/FALCON-integrate/FALCON-examples/run/greg200k-sv2
#R=/lustre/hpcprod/cdunn/jira/se-809/run

go: fc_rr_hctg_track.exe
	./fc_rr_hctg_track.exe --help
	cd $R; ${THISDIR}/fc_rr_hctg_track.exe --min-len=2500 --bestn=40 --n-core=4 --output=${THISDIR}/out1n
	cd $R; fc_rr_hctg_track.py --min-len=2500 --bestn=40 --n-core=4 --output=${THISDIR}/out1p --stream
all: fc_rr_hctg_track2.exe fc_rr_hctg_track.exe falcon/rr_hctg_track.exe
	falcon/rr_hctg_track.exe
	./fc_rr_hctg_track.exe --help
	#./fc_rr_hctg_track2.exe --help
	#./fc_rr_hctg_track2.exe --test
try: fc_rr_hctg_track2.exe
	cd $R; ${THISDIR}/fc_rr_hctg_track2.exe --min-len=2500 --bestn=40 --n-core=0 --output=${THISDIR}/out2n
	cd $R; fc_rr_hctg_track2.py --min-len=2500 --bestn=40 --n-core=0 --output=${THISDIR}/out2p
fc_rr_hctg_track2.exe fc_rr_hctg_track.exe: $(wildcard falcon/*.nim)
m.exe: m.nim

run-%: %.exe
	./$*.exe
%.exe: %.nim
	nim ${NIMFLAGS} --out:$*.exe c $<
clean:
	rm -rf nimcache/
	rm -f *.exe falcon/*.exe
