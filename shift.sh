#!/bin/bash

PLANES="$(sed '2!d' $1 | awk '{print $3}')"
N_TRACKS="$(sed '2!d' $1 | awk '{print $4}')"

echo "Prefix for output: $3"
echo "Number of frames: ${PLANES}"
echo "Number of tracks: ${N_TRACKS}"

#calculating per frame displacements for all beads
awk 'NR>3{print $1, $2, $3, $4}' $2 | awk 'NR==1{x=$3;y=$4;next} {print $1,$2,$3-x, $4-y; x=$3;y=$4}' - > ${3}_displacements.txt

/usr/local/bin/shift.py $1 ${3}_displacements.txt $3 ${PLANES} ${N_TRACKS}


