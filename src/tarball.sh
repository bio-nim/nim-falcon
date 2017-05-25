rm -f *.gz *.bz2
rm -rf archive/
mkdir -p archive/

FILES=( fc_consensus.exe      fc_rr_hctg_track.exe  fc_rr_hctg_track2.exe )

# hard-link all files
set +x
for F in "${FILES[@]}"; do
	ln -f ${F} archive/
done
set -x
cd archive/
tar cvzf ../nim-falcon${PLAT}.tar.gz .
rm -rf archive/
