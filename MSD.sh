#!/bin/bash

PLANES="$(sed '2!d' $1 | awk '{print $3}')"
N_TRACKS="$(sed '2!d' $1 | awk '{print $4}')"

echo "Number of planes: ${PLANES}"
echo "Number of tracks: ${N_TRACKS}"

awk 'NR > 3' $1 | awk '{print $1,$3,$4}' > tracking_noHeader.txt

python MSD.py tracking_noHeader.txt ${PLANES} ${N_TRACKS}

Rscript MSD.R $PLANES



