from __future__ import division
import numpy as np
import sys
from matplotlib import pyplot as plt

Nevents = 25

plot_counter = 1

file1 = open("../best_fit_results/TEMP_TMP.dat",'r')
file2 = open("../best_fit_results/TEMP_400-300-200-100.dat",'r')
file3 = open("../best_fit_results/TEMP_800-500-300-50.dat",'r')
file4 = open("../best_fit_results/TEMP_1000-100-80-30.dat",'r')
lines1 = file1.readlines()
lines2 = file2.readlines()
lines3 = file3.readlines()
lines4 = file4.readlines()

# Fill each matrix of best-fit values separately
N1 = len(lines1)-2
N2 = len(lines2)-2
N3 = len(lines3)-2
N4 = len(lines4)-2
best_fit1 = np.empty((N1, 7))
best_fit2 = np.empty((N2, 7))
best_fit3 = np.empty((N3, 7))
best_fit4 = np.empty((N4, 7))
for i in range(N1):
	words1 = lines1[i+2].split()
	best_fit1[i,:] = float(words1[1]), float(words1[2]), float(words1[3]), float(words1[4]), float(words1[6]), float(words1[7]), int(words1[8])
for i in range(N2):
	words2 = lines2[i+2].split()
	best_fit2[i,:] = float(words2[1]), float(words2[2]), float(words2[3]), float(words2[4]), float(words2[6]), float(words2[7]), int(words2[8])
for i in range(N3):
	words3 = lines3[i+2].split()
	best_fit3[i,:] = float(words3[1]), float(words3[2]), float(words3[3]), float(words3[4]), float(words3[6]), float(words3[7]), int(words3[8])
for i in range(N4):
	words4 = lines4[i+2].split()
	best_fit4[i,:] = float(words4[1]), float(words4[2]), float(words4[3]), float(words4[4]), float(words4[6]), float(words4[7]), int(words4[8])

# Nbinsmax = max([len(lines1) - 2,len(lines2) - 2,len(lines3) - 2,len(lines4) - 2])
Nbinsmin = min([len(lines1) - 2,len(lines2) - 2,len(lines3) - 2,len(lines4) - 2])

best_fit_combined = np.empty((Nbinsmin, 28))
iCombined = 0
for i1 in range(N1):
	for i2 in range(N2):
		for i3 in range(N3):
			for i4 in range(N4):
				# print best_fit1[i1,6]
				if best_fit1[i1,6]==best_fit2[i2,6]==best_fit3[i3,6]==best_fit4[i4,6]:
					best_fit_combined[iCombined,:] = np.concatenate((best_fit1[i1,:], best_fit2[i2,:], best_fit3[i3,:], best_fit4[i4,:]))
					iCombined += 1
Ncombined = iCombined


msquark = np.empty((Ncombined,4))
mchi2 = np.empty((Ncombined,4))
mslepton = np.empty((Ncombined,4))
mchi1 = np.empty((Ncombined,4))

number_of_consistent = 0 # The number of bins where all four starting points agree
eps = 0.1
for i in range(Ncombined):
	# print "%3d %3.1f   %2.1f   %2.1f   %2.1f   %3d   %f	%d" %(i+1, best_fit_combined[i,0], best_fit_combined[i,1], best_fit_combined[i,2], best_fit_combined[i,3], best_fit_combined[i,4], best_fit_combined[i,5], best_fit_combined[i,6])
	# print best_fit_combined[i,0], best_fit_combined[i,7], best_fit_combined[i,14], best_fit_combined[i,21]
	msquark[i,:] = best_fit_combined[i,0], best_fit_combined[i,7], best_fit_combined[i,14], best_fit_combined[i,21]
	mchi2[i,:] = best_fit_combined[i,1], best_fit_combined[i,1+7], best_fit_combined[i,1+14], best_fit_combined[i,1+21]
	mslepton[i,:] = best_fit_combined[i,2], best_fit_combined[i,2+7], best_fit_combined[i,2+14], best_fit_combined[i,2+21]
	mchi1[i,:] = best_fit_combined[i,3], best_fit_combined[i,3+7], best_fit_combined[i,3+14], best_fit_combined[i,3+21]
	print int(best_fit_combined[i,6]), np.mean(msquark[i,:]), np.std(msquark[i,:]), np.mean(mchi2[i,:]), np.std(mchi2[i,:]), np.mean(mslepton[i,:]), np.std(mslepton[i,:]), np.mean(mchi1[i,:]), np.std(mchi1[i,:])
	if np.std(msquark[i,:]) < eps and np.std(mchi2[i,:]) < eps and np.std(mslepton[i,:]) < eps and np.std(mchi1[i,:]) < eps:
		number_of_consistent += 1

print "Ncombined =", Ncombined, "number_of_consistent =", number_of_consistent
# print "stdev of msquark =", np.std(msquark, axis=1)
print "Mean stdev of masses =", np.mean(np.std(msquark,axis=1)), np.mean(np.std(mchi2,axis=1)), np.mean(np.std(mslepton,axis=1)), np.mean(np.std(mchi1,axis=1))
print "Stdev of stdev of masses =", np.std(np.std(msquark,axis=1)), np.std(np.std(mchi2,axis=1)), np.std(np.std(mslepton,axis=1)), np.std(np.std(mchi1,axis=1))



print "Exit"
sys.exit(0)

# Get true mass values
MSuL = 565.312 # Mass of ~uL, ~cL
MSdL = 570.734 # Mass of ~dl, ~sL
Msquark = MZ =(MSuL+MSdL)/2.0 # mean squark mass, fit this
Mchi2 = MY = 180.337 # Mass of ~chi02
Mslepton = MX = 144.06 # Mass of ~eR, ~muR
Mchi1 = MN = 9.70071979E+01 # Mass of ~chi01 (dark matter!)


# Take out best-fit values for each bin as vectors (note minuscle m for the fit-vector)
msquark = best_fit[:,0]
mchi2 = best_fit[:,1]#*100
mslepton = best_fit[:,2]
mchi1 = best_fit[:,3]
msquark_passcut = []
mchi2_passcut = []
mslepton_passcut = []
mchi1_passcut = []
correct_combo_passcut = []
cut = 100 # xi^2 cut value in units of (100 GeV)^4
for i in range(len(best_fit[:,0])):
	if best_fit[i,5] < float(cut):
		msquark_passcut.append(best_fit[i,0])
		mchi2_passcut.append(best_fit[i,1])
		mslepton_passcut.append(best_fit[i,2])
		mchi1_passcut.append(best_fit[i,3])
		correct_combo_passcut.append(best_fit[i,6])
print "Number of events passing xi^2-cut = ", len(msquark_passcut)

# Calculation of mean values and rms error for the fit
def rmse_est(estimate_vector):
	# rms deviation from mean
	n = len(estimate_vector)
	mean = np.mean(estimate_vector)
	rmse = np.sqrt( np.mean( np.power( mean*np.ones(n) - estimate_vector , 2) ) )
	return rmse

mean_msquark = np.mean(msquark_passcut)
mean_mchi2 = np.mean(mchi2_passcut)
mean_mslepton = np.mean(mslepton_passcut)
mean_mchi1 = np.mean(mchi1_passcut)

rmse_est_msquark = rmse_est(msquark_passcut)
rmse_est_mchi2 = rmse_est(mchi2_passcut)
rmse_est_mslepton = rmse_est(mslepton_passcut)
rmse_est_mchi1 = rmse_est(mchi1_passcut)

# rmse_true_msquark = rmse_true(Msquark, msquark)
# rmse_true_mchi2 = rmse_true(Mchi2, mchi2)
# rmse_true_mslepton = rmse_true(Mslepton, mslepton)
# rmse_true_mchi1 = rmse_true(Mchi1, mchi1)

# print "Mean and rmse values:"
# print "squark:  mean = %3.3f, rmse_est = %3.3f" %(mean_msquark, rmse_est_msquark)
# print "chi2:    mean = %3.3f, rmse_est = %3.3f" %(mean_mchi2, rmse_est_mchi2)
# print "slepton: mean = %3.3f, rmse_est = %3.3f" %(mean_mslepton, rmse_est_mslepton)
# print "chi1:    mean = %3.3f, rmse_est = %3.3f" %(mean_mchi1, rmse_est_mchi1)
print "Thesis-friendly numbers: mean \pm rmse_est"
print "squark : %d \pm %d" %(round(mean_msquark), round(rmse_est_msquark))
print "chi2   : %d \pm %d" %(round(mean_mchi2), round(rmse_est_mchi2))
print "slepton: %d \pm %d" %(round(mean_mslepton), round(rmse_est_mslepton))
print "chi1   : %d \pm %d" %(round(mean_mchi1), round(rmse_est_mchi1))
print "Passcut-fraction:", (len(msquark_passcut)/float(Nbins)*100)
print "Correct-combo fraction:", (np.mean(correct_combo_passcut)*100)


# Make a nice plot like Webber - msquark on y axis, mslepton, mchi2  & mchi1 on x axis

# ylim = [np.min(msquark)-30, np.max(msquark)+30]
# xlim = [np.min(np.append(mslepton,np.append(mchi1,mchi2)))-30, np.max(np.append(mslepton,np.append(mchi1,mchi2)))+30]
ylim=[400,650] # MODIFIED 20150203
xlim=[0,300]
#print xlim, ylim

plt.figure(plot_counter)
plt.scatter(mchi2_passcut, msquark_passcut, s=12, c = 'r')
# plt.xticks([100],[r'$\pi$'],fontsize=32)
plt.xlim(xlim[0],xlim[1])
plt.ylim(ylim[0],ylim[1])
plt.hold('on')
plt.scatter(mslepton_passcut, msquark_passcut, s=12, c='b')
plt.scatter(mchi1_passcut, msquark_passcut, s=12, c='y')
plt.plot(Mchi2*np.ones(2), ylim, 'r--')
plt.plot(Mslepton*np.ones(2), ylim, 'b--')
plt.plot(Mchi1*np.ones(2), ylim, 'y--')
plt.plot(xlim, Msquark*np.ones(2), 'k--')
plt.xlabel(r'$m_i \mathrm{[GeV]}$',fontsize=20)
plt.ylabel(r'$m_{\tilde q} \mathrm{[GeV]}$',fontsize=20)
plt.text(236, 638, r"$N_\mathrm{events} = %d$" % Nevents, fontsize=14)
plt.text(236, 624, r"$N_\mathrm{bins(total)} = %d$" % (Nbins), fontsize=14)
if cut > 1e10:
	plt.text(236, 610, r"$\mathrm{\xi^2_\mathrm{max}} = \infty$", fontsize=14)
else:
	plt.text(236, 610, r"$\mathrm{\xi^2_\mathrm{max}} = %d$" % cut, fontsize=14)
plt.text(236, 596, r"$f_\xi = %d \%%  $" % (len(msquark_passcut)/float(Nbins)*100) , fontsize=14)
plt.text(236, 582, r"$f_\mathrm{corr} = %d \%%$" % (np.mean(correct_combo_passcut)*100), fontsize=14)
plt.text(50,MZ+5,r'$\tilde q$',fontsize=20)
plt.text(MY+1,420,r'$\tilde\chi_2^0$',fontsize=20)
plt.text(MX+1,420,r'$\tilde l$',fontsize=20)
plt.text(MN+1,420,r'$\tilde \chi_1^0$',fontsize=20)
# plt.text(5,460, "Best-fit values (%d bins):" %Nbins)
plt.text(5,440,r"$m_{\tilde q} = %d \pm %d$" %(round(mean_msquark), round(rmse_est_msquark)),fontsize=14)
plt.text(5,430,r"$m_{\tilde \chi_2^0} = %d \pm %d$" %(round(mean_mchi2), round(rmse_est_mchi2)),fontsize=14)
plt.text(5,420,r"$m_{\tilde l} =  %d \pm %d$" %(round(mean_mslepton), round(rmse_est_mslepton)),fontsize=14)
plt.text(5,410,r"$m_{\tilde \chi_1^0} =  %d \pm %d$" %(round(mean_mchi1), round(rmse_est_mchi1)),fontsize=14)
plt.savefig('/home/jorgenem/herwigpp-MD-dileptonedge-fit-comb-nosmear-cut-100_400-300-200-100.pdf', format='pdf')

# plt.hold('off')
# plt.close()

plt.show()