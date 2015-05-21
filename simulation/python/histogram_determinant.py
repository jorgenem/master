import sys
import matplotlib.pyplot as plt
import numpy as np

infile = open("../best_fit_results/TEMPSUBDET.dat",'r')
lines = infile.readlines()

detA1 = []
detA2 = []
for line in lines:
	words = line.split()
	detA1.append(float(words[0]))
	detA2.append(float(words[1]))

trim = 50
detA1_trimmed = []
for detA1current in detA1:
	if abs(detA1current) < trim:
		detA1_trimmed.append(detA1current)
detA2_trimmed = []
for detA2current in detA2:
	if abs(detA2current) < trim:
		detA2_trimmed.append(detA2current)

# plt.figure(1)
# plt.hist(detA1_trimmed, bins=100, normed=1, color='crimson')
# # plt.xlabel(r'$\mathrm{det}(A)$')
# plt.xlabel(r'$\mathrm{subdet}(A,1,3) / \mathrm{subdet}(A,5,7)$')
# plt.ylabel('Occurences/total')
# plt.show()

# plt.figure(2)
# plt.hist(detA2_trimmed, bins=100, normed=1, color='dodgerblue')
# # plt.xlabel(r'$\mathrm{det}(A)$')
# plt.xlabel(r'$\mathrm{subdet}(A,1,3) / \mathrm{subdet}(A,5,7)$')
# plt.ylabel('Occurences/total')
# plt.show()

print np.mean(detA1), np.mean(detA2), np.std(detA1), np.std(detA2)
print np.mean(detA1_trimmed), np.mean(detA2_trimmed), np.std(detA1_trimmed), np.std(detA2_trimmed)
print np.percentile(abs(np.array(detA1)), [1,2,3,5,10,20,30,50,70,90,100])
print np.percentile(abs(np.array(detA2)), [1,2,3,5,10,20,30,50,70,90,100])


# Quantification of distribution of passcuts
detAcut = 10
Nbins = 100
Nevents = 25

number_of_passcuts_per_event = []
for iBin in range(Nbins):
	Npasscut_current = 0
	for iEvent in range(Nevents):
		if np.abs(detA1[iEvent + iBin*Nevents]) > detAcut:
			Npasscut_current += 1
	number_of_passcuts_per_event.append(Npasscut_current)
print number_of_passcuts_per_event
print np.min(number_of_passcuts_per_event), np.max(number_of_passcuts_per_event)
print np.mean(number_of_passcuts_per_event), np.std(number_of_passcuts_per_event)