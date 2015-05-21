import os
import numpy as np


file = open("../best_fit_results/TEMP-preevent_minimization.dat", 'r');
lines = file.readlines()

bestfitvalues = np.empty((len(lines)-2,6))
for i in range(len(lines)-2):
	words = lines[i+2].split()
	print words
	bestfitvalues[i,:] = int(words[0]), int(words[1]), float(words[2]), float(words[3]), float(words[4]), float(words[5])

# Count how often true combo is smallest, for each event
truecomboissmallest = []
for i in range(25):
	truecomboissmallest_current = 0
	for j in range(24):
		currentvalues = bestfitvalues[24*i+j,2], bestfitvalues[24*i+j,3], bestfitvalues[24*i+j,4], bestfitvalues[24*i+j,5]
		# print currentvalues
		if currentvalues[0] < currentvalues[1] and currentvalues[0] < currentvalues[2] and currentvalues[0] < currentvalues[3]:
			truecomboissmallest_current += 1
	truecomboissmallest.append(truecomboissmallest_current/25.0)
print truecomboissmallest

# Check how many events pass some cut
cut = 0.65
passcut = 0
for i in range(len(truecomboissmallest)):
	if truecomboissmallest[i] >= cut:
		passcut += 1
print passcut