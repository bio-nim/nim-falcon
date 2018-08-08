rm -f *.gz *.bz2
rm -rf archive/
mkdir -p archive/bin/
mkdir -p archive/dummy/
touch archive/dummy # to prevent conda from "hoisting" the bin-dir

FILES=( fc_consensus.exe      fc_rr_hctg_track.exe  fc_rr_hctg_track2.exe )

# hard-link all files
set +x
for F in "${FILES[@]}"; do
	ln -f ${F} archive/bin/
done
set -x
PLAT=$(bash ./plat.sh)
tar cvzf nim-falcon${PLAT}.tar.gz -C archive/ .
shasum -a 256 nim-falcon${PLAT}.tar.gz
rm -rf archive/
