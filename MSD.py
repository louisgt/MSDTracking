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

with open(sys.argv[1]) as f:
	for line in f:
		line = str.split(line)
		x_coord[int(line[0])-1].append(float(line[1])*0.16)
		y_coord[int(line[0])-1].append(float(line[2])*0.16)

### count number of frames per track
track_length = [0 for i in range(n_tracks)]
for i in range(n_tracks):
	track_length[i] = len(x_coord[i])

### Ensemble-averaged MSD - Parry definition

# find length of longest trajectory
# max_time
# create array MSD of length max_time
# then
# for tau in range(max_time)
#	n = len(trajectories with length >= tau)
#	for track in range(n)
#		sum+= (xi(tau) - xi(0))^2
#	MSD[tau] = sum/n

MSD = [0 for i in range(planes)]
for tau in range(planes):
	nb=0
	sum=0
	for track in range(n_tracks):
		if(track_length[track]>tau):
			#calculate square displacement for the frame
			sum+=((x_coord[track][tau]-x_coord[track][0])**2+(y_coord[track][tau]-y_coord[track][0])**2)
			nb+=1
	MSD[tau] = (sum/nb)

filename = prefix + "_MSD.txt"
with open(filename, 'w') as o:
	o.write("TAU" + "\t" + "MSD" + "\n")
	for tau in range(planes):
		o.write(str(tau) + "\t" + str(MSD[tau]) + "\n")

# i indexes the current track
for i in range(n_tracks):
	# j indexes the current dt (equal to the number of frames)
	for j in range(0,len(x_coord[i])):
		# k indexes the current frame
		sum = 0
		for k in range(len(x_coord[i])-j):
			#calculate square displacement for the frame
			sum+=((x_coord[i][k+j]-x_coord[i][k])**2+(y_coord[i][k+j]-y_coord[i][k])**2)
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
	