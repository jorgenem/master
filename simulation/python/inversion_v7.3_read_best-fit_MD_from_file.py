from __future__ import division
import numpy as np
import sys
from matplotlib import pyplot as plt

Nevents = 25

plot_counter = 1

file = open("../best_fit_results/MDTEMP.dat",'r')
lines = file.readlines()

Nbins = len(lines) - 2

best_fit = np.empty((Nbins, 6))
for i in range(Nbins):
	words = lines[i+2].split()
	best_fit[i,:] = float(words[1]), float(words[2]), float(words[3]), int(words[4]), float(words[5]), float(words[6])
	# best_fit[i,0:4] *= 100
	print "%3d  %2.1f   %2.1f   %2.1f   %3d   % .6e	%f" %(i+1, best_fit[i,0], best_fit[i,1], best_fit[i,2], best_fit[i,3], best_fit[i,4], best_fit[i,5])




# print "Exit"
# sys.exit(0)

# Get true mass values
MSuL = 565.312 # Mass of ~uL, ~cL
MSdL = 570.734 # Mass of ~dl, ~sL
Msquark = MZ =(MSuL+MSdL)/2.0 # mean squark mass, fit this
Mchi2 = MY = 180.337 # Mass of ~chi02
Mslepton = MX = 144.06 # Mass of ~eR, ~muR
Mchi1 = MN = 9.70071979E+01 # Mass of ~chi01 (dark matter!)


# Take out best-fit values for each bin as vectors (note minuscle m for the fit-vector)
m1 = best_fit[:,0]
m2 = best_fit[:,1]
m3 = best_fit[:,2]
m1_passcut = []
m2_passcut = []
m3_passcut = []
correct_combo_passcut = []
cut = 100**100 # xi^2 cut value in units of (100 GeV)^4
for i in range(len(best_fit[:,0])):
	if best_fit[i,4] < float(cut):
		m1_passcut.append(best_fit[i,0])
		m2_passcut.append(best_fit[i,1])
		m3_passcut.append(best_fit[i,2])
		correct_combo_passcut.append(best_fit[i,5])
print "Number of events passing xi^2-cut = ", len(m1_passcut)

# Calculation of mean values and rms error for the fit
def rmse_est(estimate_vector):
	# rms deviation from mean
	n = len(estimate_vector)
	mean = np.mean(estimate_vector)
	rmse = np.sqrt( np.mean( np.power( mean*np.ones(n) - estimate_vector , 2) ) )
	return rmse

mean_m1 = np.mean(m1_passcut)
mean_m2 = np.mean(m2_passcut)
mean_m3 = np.mean(m3_passcut)

rmse_est_m1 = rmse_est(m1_passcut)
rmse_est_m2 = rmse_est(m2_passcut)
rmse_est_m3 = rmse_est(m3_passcut)

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
print "M1: %d \pm %d" %(round(mean_m1), round(rmse_est_m1))
print "M2: %d \pm %d" %(round(mean_m2), round(rmse_est_m2))
print "M3: %d \pm %d" %(round(mean_m3), round(rmse_est_m3))
print "Passcut-fraction:", (len(m1_passcut)/float(Nbins)*100)
print "Correct-combo fraction:", (np.mean(correct_combo_passcut)*100)


# Plot M1 on y-axis, M2 and M3 on x-axis

# ylim = [np.min(msquark)-30, np.max(msquark)+30]
# xlim = [np.min(np.append(mslepton,np.append(mchi1,mchi2)))-30, np.max(np.append(mslepton,np.append(mchi1,mchi2)))+30]
ylim=[200000, 350000] # MODIFIED 20150330
xlim=[5000, 20000]
#print xlim, ylim

plt.figure(plot_counter)
# plt.ylabel(r'$M1 \mathrm{[GeV^2]}$',fontsize=20, horizontalalignment='right')

# First subplot
plt.subplot(1,2,1)
plt.xlabel(r'$M2 \mathrm{[GeV^2]}$',fontsize=20, verticalalignment='center')
plt.ylabel(r'$M1 \mathrm{[GeV^2]}$',fontsize=20)
plt.xticks( rotation=17)
plt.scatter(m2_passcut, m1_passcut, s=12, c = 'c')
# plt.xticks([100],[r'$\pi$'],fontsize=32)
plt.xlim(xlim[0],xlim[1])
plt.ylim(ylim[0],ylim[1])
plt.plot((Mchi2**2-Mslepton**2)*np.ones(2), ylim, 'c--')
plt.plot(xlim, (Msquark**2-Mchi2**2)*np.ones(2), 'k--')
plt.text(5500,221000,r"$M1 = %d \pm %d$" %(round(mean_m1), round(rmse_est_m1)),fontsize=14)
plt.text(5500,213000,r"$M2 = %d \pm %d$" %(round(mean_m2), round(rmse_est_m2)),fontsize=14)
plt.text(5500,205000,r"$M3 = %d \pm %d$" %(round(mean_m3), round(rmse_est_m3)),fontsize=14)

# Second subplot
plt.subplot(1,2,2)
# plt.ylabel(r'$m_{\tilde q} \mathrm{[GeV]}$',fontsize=20)
plt.xlim(xlim[0],xlim[1])
plt.ylim(ylim[0],ylim[1])
plt.yticks(rotation=17)
plt.xticks( rotation=17)
plt.scatter(m3_passcut, m1_passcut, s=12, c='m')
plt.plot((Mslepton**2-Mchi1**2)*np.ones(2), ylim, 'm--')
plt.plot(xlim, (Msquark**2-Mchi2**2)*np.ones(2), 'k--')
plt.xlabel(r'$M3 \mathrm{[GeV^2]}$',fontsize=20, verticalalignment='center')
plt.text(13000, 237000, r"$N_\mathrm{events} = %d$" % Nevents, fontsize=14)
plt.text(13000, 229000, r"$N_\mathrm{bins(total)} = %d$" % (Nbins), fontsize=14)
if cut > 1e10:
	plt.text(13000, 221000, r"$\mathrm{\xi^2_\mathrm{max}} = \infty$", fontsize=14)
else:
	plt.text(13000, 221000, r"$\mathrm{\xi^2_\mathrm{max}} = %d$" % cut, fontsize=14)
plt.text(13000, 213000, r"$f_\xi = %d \%%  $" % (len(m1_passcut)/float(Nbins)*100) , fontsize=14)
plt.text(13000, 205000, r"$f_\mathrm{corr} = %d \%%$" % (np.mean(correct_combo_passcut)*100), fontsize=14)
# plt.text(50,MZ+5,r'$\tilde q$',fontsize=20)
# plt.text(MY+1,420,r'$\tilde\chi_2^0$',fontsize=20)
# plt.text(MX+1,420,r'$\tilde l$',fontsize=20)
# plt.text(MN+1,420,r'$\tilde \chi_1^0$',fontsize=20)
# plt.text(5,460, "Best-fit values (%d bins):" %Nbins)
# plt.savefig('/home/jorgenem/herwigpp_5psmear_lowtol_nocomb_1000-100-80-30.pdf', format='pdf')

# plt.hold('off')
# plt.close()

plt.show()