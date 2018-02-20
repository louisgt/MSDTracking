#!/usr/bin/env python

import sys
import re
import math
import numpy as np

#input = (1) GFP tracking file
GFP = sys.argv[1]
#		 (2) bead displacement file
bead = sys.argv[2]
#		 (3) outpuf prefix
prefix = sys.argv[3]
#		 (4) number of frames
n_frames = int(sys.argv[4])
#		 (5) number of tracks
n_tracks = int(sys.argv[5])

shift_dx = [0 for i in range(n_frames)]
shift_dy = [0 for i in range(n_frames)]

## READ BEAD FILE
with open(bead) as f:
	#x displacements
	list_dx = [[] for i in range(n_frames)]
	#y displacements
	list_dy = [[] for i in range(n_frames)]
	for line in f:
		line = str.split(line)
		list_dx[int(line[1])-1].append(float(line[2]))
		list_dy[int(line[1])-1].append(float(line[3]))

	for frame in range(1,n_frames):
		shift_dx[frame]=np.mean(list_dx[frame],0)
		if(math.isnan(shift_dx[frame])): shift_dx[frame]=0
		shift_dy[frame]=np.mean(list_dy[frame],0)
		if(math.isnan(shift_dy[frame])): shift_dy[frame]=0
	
filename = prefix + "_corrected.txt"
with open(filename, 'w') as o:
	o.write("TIME" + "\t" + "XSHIFT" + "\t" + "YSHIFT" + "\n")
	for tau in range(n_frames):
		o.write(str(tau) + "\t" + str(shift_dx[tau]) + "\t" + str(shift_dy[tau]) + "\n")

#x coordinates vector
x_coord = [[] for i in range(n_tracks)]
#y coordinates vector
y_coord = [[] for i in range(n_tracks)]

intensity_single = [[] for j in range(n_tracks)]

with open(GFP) as f:
	head1 = next(f)
	head2 = next(f)
	head3 = next(f)
	head3 = str.split(head3)
	del head3[4]
	del head3[4]
	head3 = '\t'.join(head3)
	for line in f:
		line = str.split(line)
		x_coord[int(line[0])-1].append(float(line[2])-shift_dx[int(line[1])-1])
		y_coord[int(line[0])-1].append(float(line[3])-shift_dy[int(line[1])-1])
		intensity_single[int(line[0])-1].append(float(line[6]))

### count number of frames per track
track_length = [0 for i in range(n_tracks)]
for i in range(n_tracks):
	track_length[i] = len(x_coord[i])

filename = prefix + "_tracking_corrected.txt"
with open(filename, 'w') as o:
	o.write(head1)
	o.write(head2)
	o.write(head3+"\n")
	for i in range(n_tracks):
		for frame in range(track_length[i]):
			o.write(str(i+1) + "\t" + str(frame+1) + "\t" + str(x_coord[i][frame]) + "\t" + str(y_coord[i][frame]) + "\t" + str(intensity_single[i][frame]) + "\n")