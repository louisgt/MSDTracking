#!/bin/bash

PLANES="$(sed '2!d' $1 | awk '{print $3}')"
N_TRACKS="$(sed '2!d' $1 | awk '{print $4}')"

echo "Prefix for output: $2"
echo "Number of frames: ${PLANES}"
echo "Number of tracks: ${N_TRACKS}"

awk 'NR > 3' $1 | awk '{print $1,$3,$4}' > tracking_noHeader.txt

/usr/local/bin/MSD.py tracking_noHeader.txt ${PLANES} ${N_TRACKS} $2

Rscript /usr/local/bin/MSD.R $2

rm tracking_noHeader.txt



