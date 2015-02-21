from __future__ import division
import numpy as np
import sys
from matplotlib import pyplot as plt

Nevents = 25
Nbins = 100

plot_counter = 1

file = open("../best_fit_results/TEST.dat",'r')
lines = file.readlines()

best_fit = np.empty((Nbins, 6))
for i in range(Nbins):
	words = lines[i+2].split()
	best_fit[i,:] = float(words[1]), float(words[2]), float(words[3]), float(words[4]), int(words[5]), float(words[6])
	# best_fit[i,0:4] *= 100
	print "%3d %3.1f   %2.1f   %2.1f   %2.1f   %3d   % .6e" %(i+1, best_fit[i,0], best_fit[i,1], best_fit[i,2], best_fit[i,3], best_fit[i,4], best_fit[i,5])




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
msquark = best_fit[:,0]
mchi2 = best_fit[:,1]#*100
mslepton = best_fit[:,2]
mchi1 = best_fit[:,3]
msquark_passcut = []
mchi2_passcut = []
mslepton_passcut = []
mchi1_passcut = []
cut = 100**100 # xi^2 cut value in units of (100 GeV)^4
for i in range(len(best_fit[:,0])):
	if best_fit[i,5] < float(cut):
		msquark_passcut.append(best_fit[i,0])
		mchi2_passcut.append(best_fit[i,1])
		mslepton_passcut.append(best_fit[i,2])
		mchi1_passcut.append(best_fit[i,3])
print "Number of events passing xi^2-cut = ", len(msquark_passcut)

# Calculation of mean values and rms error for the fit
def rmse_est(estimate_vector):
	# rms deviation from mean
	n = len(estimate_vector)
	mean = np.mean(estimate_vector)
	rmse = np.sqrt( np.mean( np.power( mean*np.ones(n)-estimate_vector , 2) ) )
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


# Make a nice plot like Webber - msquark on y axis, mslepton, mchi2  & mchi1 on x axis

# ylim = [np.min(msquark)-30, np.max(msquark)+30]
# xlim = [np.min(np.append(mslepton,np.append(mchi1,mchi2)))-30, np.max(np.append(mslepton,np.append(mchi1,mchi2)))+30]
ylim=[400,650] # MODIFIED 20150203
xlim=[0,300]
#print xlim, ylim

plt.figure(plot_counter)
plt.plot(mchi2_passcut, msquark_passcut, 'ro')
# plt.xticks([100],[r'$\pi$'],fontsize=32)
plt.xlim(xlim[0],xlim[1])
plt.ylim(ylim[0],ylim[1])
plt.hold('on')
plt.plot(mslepton_passcut, msquark_passcut, 'bo')
plt.plot(mchi1_passcut, msquark_passcut, 'yo')
plt.plot(Mchi2*np.ones(2), ylim, 'r--')
plt.plot(Mslepton*np.ones(2), ylim, 'b--')
plt.plot(Mchi1*np.ones(2), ylim, 'y--')
plt.plot(xlim, Msquark*np.ones(2), 'k--')
plt.xlabel(r'$m_i \mathrm{[GeV]}$',fontsize=20)
plt.ylabel(r'$m_{\tilde q} \mathrm{[GeV]}$',fontsize=20)
plt.title("Nevents = %d, cut = %.2e, passcut-fraction = %.3f"%(Nevents,cut,len(msquark_passcut)/float(Nbins)))
plt.text(50,MZ+5,r'$\tilde q$',fontsize=20)
plt.text(MY+1,420,r'$\tilde\chi_2^0$',fontsize=20)
plt.text(MX+1,420,r'$\tilde l$',fontsize=20)
plt.text(MN+1,420,r'$\tilde \chi_1^0$',fontsize=20)
plt.text(5,460, "Best-fit values:")
plt.text(5,450,r"squark : %d $\pm$ %d" %(round(mean_msquark), round(rmse_est_msquark)),fontsize=10)
plt.text(5,440,r"chi2   : %d $\pm$ %d" %(round(mean_mchi2), round(rmse_est_mchi2)),fontsize=10)
plt.text(5,430,r"slepton: %d $\pm$ %d" %(round(mean_mslepton), round(rmse_est_mslepton)),fontsize=10)
plt.text(5,420,r"chi1   : %d $\pm$ %d" %(round(mean_mchi1), round(rmse_est_mchi1)),fontsize=10)
# plt.savefig('100_bins_25_events_pythia_events_nelder-mead_%1.2f_initial_guess_no_ISR_or_FSR_xisquared-cut.pdf'%mass_offset, format='pdf')

# plt.hold('off')
# plt.close()

plt.show()