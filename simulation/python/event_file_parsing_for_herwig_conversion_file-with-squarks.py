from __future__ import division
import sys
from math import floor
infile = open("../events/herwigpp-3820-events-momcons-including-gammas-20150316_MeV.dat",'r')
outfile = open("../events/herwigpp-3820-events-momcons-including-gammas-20150316.dat", 'w')

lines = infile.readlines()

for i in range(11*int(len(lines)/11)):
	line = lines[i]
	if i%11==0:
		outfile.write(line)
	elif i%11 in [2,3,4,5, 7,8,9,10]:
		words = line.split()
		outfile.write("%d %f %f %f %f %f\n" %(int(words[0]), float(words[1])/1000, float(words[2])/1000, float(words[3])/1000, float(words[4])/1000, float(words[5])/1000))