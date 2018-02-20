#!/bin/bash

PLANES="$(sed '2!d' $1 | awk '{print $3}')"
N_TRACKS="$(sed '2!d' $1 | awk '{print $4}')"

echo "Prefix for output: $2"
echo "Number of frames: ${PLANES}"
echo "Number of tracks: ${N_TRACKS}"
echo "Frame interval: $3"

awk 'NR > 3' $1 | awk '{print $1,$2,$3,$4,$5}' > ${2}_tracking_noHeader.txt

Rscript /usr/local/bin/binData.R $2

/usr/local/bin/MSD.py ${2}_tracking_noHeader.txt ${PLANES} ${N_TRACKS} $2

Rscript /usr/local/bin/MSD.R $2 $3



