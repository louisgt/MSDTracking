#!/usr/bin/env python

import sys
import re
import math
import numpy as np
import scipy.stats as sp

planes = int(sys.argv[2])
n_tracks = int(sys.argv[3])
prefix = sys.argv[4].replace(' ', '')
min_window = int(sys.argv[5])
pixel_size = float(sys.argv[6])
max_d = float(sys.argv[7])

#x coordinates vector
x_coord = [[] for i in range(n_tracks)]
#y coordinates vector
y_coord = [[] for i in range(n_tracks)]

#registers corresponding bin for trajectory i
track_bin = [0 for i in range(n_tracks)]

intensity_single = [[] for j in range(n_tracks)]

with open(sys.argv[1]) as f:
	for line in f:
		line = str.split(line)
		x_coord[int(line[0])-1].append(float(line[2]))
		y_coord[int(line[0])-1].append(float(line[3]))
		intensity_single[int(line[0])-1].append(float(line[4]))
		if track_bin[int(line[0])-1] is 0:
			track_bin[int(line[0])-1] = (int(line[5])-1)

#number of bins
n_bins = len(set(track_bin))

#center of gyration
center = [0 for j in range(n_tracks)]

### count number of frames per track
track_length = [0 for i in range(n_tracks)]
for i in range(n_tracks):
	track_length[i] = len(x_coord[i])
	center[i] = reduce(lambda x, y: x + y, x_coord[i]) / float(track_length[i])

max_frame = max(track_length)

### Ensemble-averaged MSD

MSD_binned = [[0 for i in range(max_frame)] for j in range(n_bins)]
MSD_single = [[0 for i in range(max_frame)] for j in range(n_tracks)]
MSD_ens = [0 for i in range(max_frame)]
MSD_ens_sem = [0 for i in range(max_frame)]

displacements = [[] for j in range(n_bins)]

#radiusg = [0 for j in range(n_tracks)]

# looping on all time points
# tau denotes the time point used to compute MSD
for tau in range(max_frame):

	nb = [0 for i in range(n_bins)]
	sum = [0 for i in range(n_bins)]

	msd = []
	single = 0
	count = 0

	# go through all trajectories
	for track in range(n_tracks):
		# is the trajectory longer than the minimum window?
		if track_length[track] > min_window:
			# is the trajectory longer than the step?
			if track_length[track] > tau:

				#calculate square displacement for the frame
				#xd : x axis displacement
				#yd : y axis displacement
				#fd : frame displacement

				xd = (x_coord[track][tau]-x_coord[track][0])*pixel_size
				yd = (y_coord[track][tau]-y_coord[track][0])*pixel_size
				fd = math.sqrt((xd**2)+(yd**2))

				MSD_single[track][tau] = (fd**2)
				msd.append(fd**2)
				count += 1

				sum[track_bin[track]] += (fd**2)
				nb[track_bin[track]] += 1



	#MSD_ens[tau] = (msd/count)
	MSD_ens[tau] = np.mean(msd,0)
	MSD_ens_sem[tau] = sp.sem(msd)

	for b in range(n_bins):
		if nb[b] is not 0:
			MSD_binned[b][tau] = (sum[b]/nb[b])

##############################
### Time-ensemble-averaged MSD

MSD_time = []
MSD_time_sem = []

# dt indexes the current time step (or interval)
# range of dt depends on the track
for dt in range(1, max_frame):
	track_count = 0
	cumul = 0
	sem = []
	# loop through each track
	for track in range(n_tracks):
			#verify that interval is computable for this track
			if track_length[track] > dt:
				track_count += 1
				dt_count = len(x_coord[track])-dt
				msdt = 0
				# k indexes the current frame
				for k in range(dt_count):
					#calculate square displacement for the frame
					#xd : x axis displacement
					#yd : y axis displacement
					#fd : frame displacement

					xd = (x_coord[track][k+dt]-x_coord[track][k])*pixel_size
					yd = (y_coord[track][k+dt]-y_coord[track][k])*pixel_size
					fd = math.sqrt((xd**2)+(yd**2))

					## record all displacements for minimal step
					if dt == 1:
						displacements[track_bin[track]].append(fd)

					msdt += (fd**2)
					sem.append(fd**2)

				cumul += (msdt/dt_count)

			else:
				continue

	MSD_time.append(cumul/track_count)
	MSD_time_sem.append(sp.sem(sem))

filename = prefix + "_MSD_all.txt"
with open(filename, 'w') as o:
	o.write("TAU" + "\t" + "MSD" + "\t" + "SE" + "\t" + "TYPE" + "\n")
	#o.write("TAU" + "\t" + "MSD" + "\n")
	for tau in range(max_frame):
		#o.write(str(tau) + "\t" + str(MSD_ens[tau]) + "\n")
		o.write(str(tau) + "\t" + str(MSD_ens[tau]) + "\t" + str(MSD_ens_sem[tau]) + "\t" + "ENSEMBLE" + "\n")

filename = prefix + "_MSD_single.txt"
with open(filename, 'w') as o:
	o.write("TRACK" + "\t" + "TAU" + "\t" + "MSD" + "\t" + "INTENSITY" + "\t" + "BIN" + "\n")
	for track in range(n_tracks):
		for tau in range(track_length[track]):
			o.write(str(track+1) + "\t" + str(tau) + "\t" + str(MSD_single[track][tau]) +
			        "\t" + str(intensity_single[track][tau]) +
			        "\t" + str(track_bin[track]) + "\n")

filename = prefix + "_MSD_binned.txt"
with open(filename, 'w') as o:
	o.write("BIN" + "\t" + "TAU" + "\t" + "MSD" + "\n")
	for b in range(n_bins):
		for tau in range(max_frame):
			o.write(str(b+1) + "\t" + str(tau) + "\t" + str(MSD_binned[b][tau]) + "\n")

filename = prefix + "_MSD_tau.txt"
with open(filename, 'w') as o:
	#o.write("TAU" + "\t" + "VALUE" + "\n")
	o.write("TAU" + "\t" + "MSD" + "\t" + "SE" + "\t" + "TYPE" + "\n")
	# i = track
	for i in range(len(MSD_time)):
		#o.write(str(i) + "\t" + str(MSD_time[i]) + "\n")
		o.write(str(i) + "\t" + str(MSD_time[i]) + "\t" + str(MSD_time_sem[i]) + "\t" + "TIME+ENSEMBLE" + "\n")

filename = prefix + "_displacements.txt"
with open(filename, 'w') as o:
	for b in range(n_bins):
		for d in displacements[b]:
			o.write(str(b) + "\t" + str(d) + "\n")