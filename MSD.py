#!/usr/bin/env python

import sys
import re

planes = int(sys.argv[2])
n_tracks = int(sys.argv[3])
prefix = sys.argv[4].replace(' ','')

#x coordinates vector
x_coord = [[] for i in range(n_tracks)]
#y coordinates vector
y_coord = [[] for i in range(n_tracks)]

#delta r squared displacement vector
dr_sq = [[] for i in range(n_tracks)]

#delta x squared displacement vector
dx_sq = [[] for i in range(n_tracks)]

#delta y squared displacement vector
dy_sq = [[] for i in range(n_tracks)]

#registers corresponding bin for trajectory i
track_bin = [0 for i in range(n_tracks)]

with open(sys.argv[1]) as f:
	for line in f:
		line = str.split(line)
		x_coord[int(line[0])-1].append(float(line[1]))
		y_coord[int(line[0])-1].append(float(line[2]))
		if track_bin[int(line[0])-1] is 0:
			track_bin[int(line[0])-1] = (int(line[4])-1)

#number of bins
n_bins = len(set(track_bin))

### count number of frames per track
track_length = [0 for i in range(n_tracks)]
for i in range(n_tracks):
	track_length[i] = len(x_coord[i])

max_frame = max(track_length)

### Ensemble-averaged MSD

MSD = [[0 for i in range(max_frame)] for j in range(n_bins)]
MSD_all = [0 for i in range(max_frame)]
x_MSD = [[0 for i in range(max_frame)] for j in range(n_bins)]
y_MSD = [[0 for i in range(max_frame)] for j in range(n_bins)]
x_MD = [[0 for i in range(max_frame)] for j in range(n_bins)]
y_MD = [[0 for i in range(max_frame)] for j in range(n_bins)]

for tau in range(max_frame):

	nb_all = 0
	sum_all = 0
	nb = [0 for i in range(n_bins)]
	sum = [0 for i in range(n_bins)]
	dx_sum = [0 for i in range(n_bins)]
	dy_sum = [0 for i in range(n_bins)]
	sx_sum = [0 for i in range(n_bins)]
	sy_sum = [0 for i in range(n_bins)]

	for track in range(n_tracks):
		if(track_length[track]>tau):
			#calculate square displacement for the frame
			sum_all +=(((x_coord[track][tau]-x_coord[track][0])**2)+((y_coord[track][tau]-y_coord[track][0])**2))
			sum[track_bin[track]]+=(((x_coord[track][tau]-x_coord[track][0])**2)+((y_coord[track][tau]-y_coord[track][0])**2))
			dx_sum[track_bin[track]]+=(x_coord[track][tau]-x_coord[track][0])
			dy_sum[track_bin[track]]+=(y_coord[track][tau]-y_coord[track][0])
			sx_sum[track_bin[track]]+=((x_coord[track][tau]-x_coord[track][0])**2)
			sy_sum[track_bin[track]]+=((y_coord[track][tau]-y_coord[track][0])**2)
			nb[track_bin[track]]+=1
			nb_all+=1

	MSD_all[tau] = (sum_all/nb_all)
	for b in range(n_bins):
		if(nb[b] is not 0):
			MSD[b][tau] = (sum[b]/nb[b])
			x_MD[b][tau] = (dx_sum[b]/nb[b])
			y_MD[b][tau] = (dy_sum[b]/nb[b])
			x_MSD[b][tau] = (sx_sum[b]/nb[b])
			y_MSD[b][tau] = (sy_sum[b]/nb[b])

filename = prefix + "_MSD_all.txt"
with open(filename, 'w') as o:
	o.write("TAU" + "\t" + "MSD" + "\n")
	for b in range(n_bins):
		for tau in range(max_frame):
			o.write(str(tau) + "\t" + str(MSD_all[tau]) + "\n")

filename = prefix + "_MSD_binned.txt"
with open(filename, 'w') as o:
	o.write("BIN" + "\t" + "TAU" + "\t" + "MSD" + "\n")
	for b in range(n_bins):
		for tau in range(max_frame):
			o.write(str(b+1) + "\t" + str(tau) + "\t" + str(MSD[b][tau]) + "\n")

filename = prefix + "_X_MSD.txt"
with open(filename, 'w') as o:
	o.write("BIN" + "\t" + "TAU" + "\t" + "MSD" + "\n")
	for b in range(n_bins):		
		for tau in range(max_frame):
			o.write(str(b+1) + "\t" + str(tau) + "\t" + str(x_MSD[b][tau]) + "\n")

filename = prefix + "_Y_MSD.txt"
with open(filename, 'w') as o:
	o.write("BIN" + "\t" + "TAU" + "\t" + "MSD" + "\n")
	for b in range(n_bins):
		for tau in range(max_frame):
			o.write(str(b+1) + "\t" + str(tau) + "\t" + str(y_MSD[b][tau]) + "\n")

filename = prefix + "_X_MD.txt"
with open(filename, 'w') as o:
	o.write("BIN" + "\t" + "TAU" + "\t" + "MD" + "\n")
	for b in range(n_bins):
		for tau in range(max_frame):
			o.write(str(b+1) + "\t" + str(tau) + "\t" + str(x_MD[b][tau]) + "\n")

filename = prefix + "_Y_MD.txt"
with open(filename, 'w') as o:
	o.write("BIN" + "\t" + "TAU" + "\t" + "MD" + "\n")
	for b in range(n_bins):
		for tau in range(max_frame):
			o.write(str(b+1) + "\t" + str(tau) + "\t" + str(y_MD[b][tau]) + "\n")

### Time-ensemble-averaged MSD

# i indexes the current track
for i in range(n_tracks):
	# j indexes the current dt (equal to the number of frames)
	for j in range(0,len(x_coord[i])):
		# k indexes the current frame
		sum = 0
		for k in range(len(x_coord[i])-j):
			#calculate square displacement for the frame
			sum+=(((x_coord[i][k+j]-x_coord[i][k]))**2+((y_coord[i][k+j]-y_coord[i][k]))**2)
		##HERE: compute MSD
		sum=sum/(len(x_coord[i])-j)
		dr_sq[i].append(sum)

filename = prefix + "_MSD_tau.txt"
with open(filename, 'w') as o:
	o.write("TRACK" + "\t" + "FRAME" + "\t" + "VALUE" + "\n")
	# i = track
	for i in range(n_tracks):
		# j = step
		for j in range(len(dr_sq[i])):
			o.write(str(i) + "\t" + str(j) + "\t" + str(dr_sq[i][j]) + "\n")
	