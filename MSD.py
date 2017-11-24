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

# i indexes the current track
for i in range(n_tracks):
	# j indexes the current dt (equal to the number of frames)
	for j in range(1,len(x_coord[i])):
		# k indexes the current frame
		sum = 0
		for k in range(len(x_coord[i])-j):
			#calculate square displacement for the frame
			sum+=((x_coord[i][k+j]-x_coord[i][k])**2+(y_coord[i][k+j]-y_coord[i][k])**2)
		##HERE: compute MSD
		sum=sum/(len(x_coord[i])-j)
		dr_sq[i].append(sum)

filename = prefix + "_melt.txt"
with open(filename, 'w') as o:
	o.write("TRACK" + "\t" + "FRAME" + "\t" + "VALUE" + "\n")
	# i = track
	for i in range(n_tracks):
		# j = step
		for j in range(len(dr_sq[i])):
			o.write(str(i) + "\t" + str(j) + "\t" + str(dr_sq[i][j]) + "\n")
	