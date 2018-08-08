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
os=${OSTYPE//[0-9.]/}

if [[ "$os" == 'linux-gnu' ]]; then
	PLAT='-linux-64'
    SHA='sha256sum'
elif [[ "$os" == 'darwin' ]]; then
	PLAT='-osx-64'
    SHA='sha -a 256'
else
	PLAT='-unknown'
    SHA='echo'
fi
set -x
tar cvzf nim-falcon${PLAT}.tar.gz -C archive/ .
${SHA} nim-falcon${PLAT}.tar.gz
rm -rf archive/
