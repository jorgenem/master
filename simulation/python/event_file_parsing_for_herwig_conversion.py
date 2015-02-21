from __future__ import division
import sys
from math import floor
infile = open("../herwigpp/LHC-MSSM-analysis_20150116_added_gluinos_and_turned_off_threebody_and_discarded_momentum-nonconserving_events.log",'r')
outfile = open("../events/Herwig_chain_20150116_with_gluinos_and_no_threebody_decay_and_discarded_momentum-nonconservation_GeV-corrected.dat", 'w')

lines = infile.readlines()

for i in range(9*(int(len(lines)/9))):
	line = lines[i]
	if i%9==0:
		outfile.write(line)
	else:
		words = line.split()
		outfile.write("%d %f %f %f %f %f\n" %(int(words[0]), float(words[1])/1000, float(words[2])/1000, float(words[3])/1000, float(words[4])/1000, float(words[5])/1000))