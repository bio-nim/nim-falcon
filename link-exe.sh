#!/bin/bash
bin=$1/bin
echo $bin

if [[ -z $bin ]]; then
    exit
fi

set -vex
mkdir -p $bin

ln -f ${NIMBLE_DIR}/bin/example $bin/fc_example.exe
ln -f ${NIMBLE_DIR}/bin/fc_consensus $bin/fc_consensus.exe
ln -f ${NIMBLE_DIR}/bin/fc_rr_hctg_track $bin/fc_rr_hctg_track.exe
ln -f ${NIMBLE_DIR}/bin/fc_rr_hctg_track2 $bin/fc_rr_hctg_track2.exe
ln -f ${NIMBLE_DIR}/bin/LA4Falcon $bin/fc_la4falcon.exe
