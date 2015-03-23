# Program to make a scan over many minimizations, and investigate which give the lowest best-fit point
import numpy as np 
from subprocess import call
import sys
import os


# Give starting points as lists
sp_sq_list = [200,400,568,600,800]
sp_chi2_list = [100,180,250,400]
sp_sl_list = [50,100,144,200,250]
sp_chi1_list = [20,50,97,150,200]
# Specify maxiter and tolerance
maxiter = 500
tol = 1e-12

summaryfile = open("summaryfile.txt",'w')
splist = []
best_fit_list = []
Nbinslist = []

# Loop over all starting point combinations, omitting unphysical ones if there are any
counter = 0
for sp_sq in sp_sq_list:

	for sp_chi2 in sp_chi2_list:
		if sp_chi2 >= sp_sq:
			continue

		for sp_sl in sp_sl_list:
			if sp_sl >= sp_chi2:
				continue

			for sp_chi1 in sp_chi1_list:
				if sp_chi1 >= sp_sl:
					continue

				counter += 1

				# RUN MINIMIZATION

				# Compile only on first run
				if counter == 1:
					os.system("rm main_ms")
					os.system("g++ -std=c++11 minimization6.1_combinatorics-count_avoid_unphysical_regions_minimization_scan.cpp -o main_ms -O3 -I/usr/local/include  -I/usr/local/include -larmadillo")

				print sp_sq, sp_chi2, sp_sl, sp_chi1

				# Delete the starting point file if it exists
				try:
				    os.remove("currentsp.tmp")
				except OSError:
				    pass
				currentspfile = open("currentsp.tmp",'w')
				currentspfile.write("%f %f %f %f" %(sp_sq, sp_chi2, sp_sl, sp_chi1))
				currentspfile.close()

				# os.system("./run_minimization_scan.sh %d %e" %(maxiter, tol))
				print "Calling minimization run number",counter
				os.system("./main_ms")


				# READ MINIMIZATION RESULTS
				file = open("../best_fit_results/TEMP_MS.dat",'r')
				lines = file.readlines()

				Nbins = len(lines) - 2
				Nbinslist.append(Nbins)

				best_fit = np.empty((Nbins, 7))
				for i in range(Nbins):
					words = lines[i+2].split()
					best_fit[i,:] = float(words[1]), float(words[2]), float(words[3]), float(words[4]), int(words[5]), float(words[6]), float(words[7])
					# best_fit[i,0:4] *= 100
					print "%3d %3.1f   %2.1f   %2.1f   %2.1f   %3d   % .6e	%f" %(i+1, best_fit[i,0], best_fit[i,1], best_fit[i,2], best_fit[i,3], best_fit[i,4], best_fit[i,5], best_fit[i,6])


				# STORE DATA IN LISTS
				best_fit_list.append(best_fit)
				splist.append([sp_sq, sp_chi2, sp_sl, sp_chi1])

# LOOP OVER ALL BINS, FIND SMALLEST FIT VALUE FOR EACH BIN
smallest_fit_list = []
for iBin in range(Nbins):
	iFit_smallest = 0
	for iFit in range(len(splist)):
		print best_fit_list[iFit][iBin,5]
		# For each bin, get the iFit index of the smallest value
		if best_fit_list[iFit][iBin,5] < best_fit_list[iFit_smallest][iBin,5]:
			iFit_smallest = iFit
	smallest_fit_list.append(iFit_smallest)
print smallest_fit_list
print "len(splist) = ", len(splist)

# Write stuff to file!
outfile = open("minimization_scan_results.txt",'w')
print>>outfile, len(splist)
print>>outfile, smallest_fit_list
print>>outfile, best_fit_list
print>>outfile, splist