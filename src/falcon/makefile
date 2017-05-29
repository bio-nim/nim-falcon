default: run-rr_hctg_track
run-%: %.exe
	./$*.exe
%.exe: %.nim
	nim c ${NIMFLAGS} --out:$*.exe $*
