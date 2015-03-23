from __future__ import division
import sys
from math import floor
infile = open("../events/herwigpp_9385_events_20150225.dat",'r')
outfile = open("../events/herwigpp_9385_events_20150225.dat_GeV", 'w')

lines = infile.readlines()

for i in range(9*(int(len(lines)/9))):
	line = lines[i]
	if i%9==0:
		outfile.write(line)
	else:
		words = line.split()
		outfile.write("%d %f %f %f %f %f\n" %(int(words[0]), float(words[1])/1000, float(words[2])/1000, float(words[3])/1000, float(words[4])/1000, float(words[5])/1000))