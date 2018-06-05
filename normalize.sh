#!/bin/bash

PLANES="$(sed '2!d' $1 | awk '{print $3}')"
N_TRACKS="$(sed '2!d' $1 | awk '{print $4}')"

echo "Prefix for output: $2"
echo "Number of frames: ${PLANES}"
echo "Number of tracks: ${N_TRACKS}"

## remove header and normalize intensities with R script
awk 'NR > 3' $1 | awk '{print $1,$2,$3,$4,$5}' > ${2}_tracking_noHeader.txt

Rscript /usr/local/bin/normalize.R $2

## then replace original file with header + normalized
head -n 3 $1 > tracking_${2}.txt
cat ${2}_tracking_noHeader.txt >> tracking_${2}.txt
rm ${2}_tracking_noHeader.txt