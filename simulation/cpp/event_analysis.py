# Program for analyzing events: distributions of angles, momenta etc
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from math import floor

def minkowskidot(a,b):
	# Inner product in Minkowski space, E is *last* index
	return float(a[0,3]*b[0,3]-a[0,0]*b[0,0]-a[0,1]*b[0,1]-a[0,2]*b[0,2])
def minkowskinorm(a):
	# Inner product of a with a in Minkowski space
	return minkowskidot(a,a)
def etacal(theta):
	return -np.log( np.tan( theta/2 ) )
def getangles(p):
	# returns the transverse polar angle, azimuthal angle 
	# and pseudorapidity of a 4-vector relative to the beam axis
	theta = np.arccos(p[0,2]/np.sqrt( p[0,0]**2 + p[0,1]**2 + p[0,2]**2 ))
	phi = np.arctan( p[0,1]/p[0,0] )
	eta = etacal(theta)
	return phi, theta, eta


import sys
# file = open("../events/herwigpp-9563-events-complete-momcons-20150314_only_SFL.dat", 'r')
file = open("../events/Pythia_cascade_10000_events_everything_turned_on_20150210_only_opposite_flavour_leptons.dat")
lines = file.readlines()
N = int(floor(len(lines)/9))
# N = 1
events = []

for i in range(N):
	# Loop over events to get 4-vectors for each particle for each event. 
	# Particles are numbered according to Webber (arXiv:0907.5307v2) fig. 1
	# (the lepton/antilepton ordering is arbitrary in each chain, the lepton has been 
	# chosen as 2/6 and the antilepton as 3/7)

	# Read all particles from file
	try:
		# chain 1
		quark1 = lines[9*i + 1].split()
		p1 = np.matrix([ float(quark1[1]), float(quark1[2]), float(quark1[3]), float(quark1[4]), int(quark1[0]) ])
		lepton11 = lines[9*i + 2].split()
		p2 = np.matrix([ float(lepton11[1]), float(lepton11[2]), float(lepton11[3]), float(lepton11[4]), int(lepton11[0]) ])
		lepton12 = lines[9*i + 3].split()
		p3 = np.matrix([ float(lepton12[1]), float(lepton12[2]), float(lepton12[3]), float(lepton12[4]), int(lepton12[0]) ])
		neutralino1 = lines[9*i + 4].split()
		p4 = np.matrix([ float(neutralino1[1]), float(neutralino1[2]), float(neutralino1[3]), float(neutralino1[4]), int(	neutralino1[0]) ])
		#chain2
		quark2 = lines[9*i + 5].split()
		p5 = np.matrix([ float(quark2[1]), float(quark2[2]), float(quark2[3]), float(quark2[4]), int(quark2[0]) ])
		lepton21 = lines[9*i + 6].split()
		p6 = np.matrix([ float(lepton21[1]), float(lepton21[2]), float(lepton21[3]), float(lepton21[4]), int(lepton21[0]) ])
		lepton22 = lines[9*i + 7].split()
		p7 = np.matrix([ float(lepton22[1]), float(lepton22[2]), float(lepton22[3]), float(lepton22[4]), int(lepton22[0]) ])
		neutralino2 = lines[9*i + 8].split()
		p8 = np.matrix([ float(neutralino2[1]), float(neutralino2[2]), float(neutralino2[3]), float(neutralino2[4]), int(	neutralino2[0]) ])

		events.append([p1,p2,p3,p4,p5,p6,p7,p8])
		i += 1
	except:
		# End of file
		break


# DO ANALYSIS

# Dilepton invariant mass spectrum
minv_dilepton = []
for i in range(N):
	minv_dilepton.append(abs(np.sqrt(minkowskinorm(events[i][1]+events[i][2]))))	
	minv_dilepton.append(abs(np.sqrt(minkowskinorm(events[i][5]+events[i][6]))))

# plt.hist(minv_dilepton,bins=20)
# plt.show()

# # Distribution of opening angles
# angle_between_dileptons_theta = []
# angle_between_dileptons_phi = []
# angle_between_quark_and_first_lepton_theta = []
# angle_between_quark_and_second_lepton_theta = []
# angle_between_quark_and_first_lepton_phi = []
# angle_between_quark_and_second_lepton_phi = []
# angle_between_quark_and_opposite_first_lepton_theta = []
# angle_between_quark_and_opposite_second_lepton_theta = []
# angle_between_quark_and_opposite_first_lepton_phi = []
# angle_between_quark_and_opposite_second_lepton_phi = []
# angle_between_opposite_leptons_first_first_theta = []
# angle_between_opposite_leptons_first_first_phi = []
# angle_between_opposite_leptons_first_second_theta = []
# angle_between_opposite_leptons_first_second_phi = []

# for i in range(N):
# 	phi0, theta0, eta0 = getangles(events[i][0])
# 	phi1, theta1, eta1 = getangles(events[i][1])
# 	phi2, theta2, eta2 = getangles(events[i][2])
# 	phi4, theta4, eta4 = getangles(events[i][4])
# 	phi5, theta5, eta5 = getangles(events[i][5])
# 	phi6, theta6, eta6 = getangles(events[i][6])

# 	angle_between_dileptons_phi.append(abs(phi1-phi2))
# 	angle_between_dileptons_phi.append(abs(phi5-phi6))
# 	angle_between_dileptons_theta.append(abs(theta1-theta2))
# 	angle_between_dileptons_theta.append(abs(theta5-theta6))
# 	angle_between_quark_and_first_lepton_theta.append(abs(theta0 - theta1))
# 	angle_between_quark_and_first_lepton_theta.append(abs(theta4 - theta5))
# 	angle_between_quark_and_second_lepton_theta.append(abs(theta0 - theta2))
# 	angle_between_quark_and_second_lepton_theta.append(abs(theta4 - theta6))
# 	angle_between_quark_and_first_lepton_phi.append(abs(phi0 - phi1))
# 	angle_between_quark_and_first_lepton_phi.append(abs(phi4 - phi5))
# 	angle_between_quark_and_second_lepton_phi.append(abs(phi0 - phi2))
# 	angle_between_quark_and_second_lepton_phi.append(abs(phi4 - phi6))
# 	angle_between_quark_and_opposite_first_lepton_theta.append(abs(theta0-theta5))
# 	angle_between_quark_and_opposite_first_lepton_theta.append(abs(theta4-theta1))
# 	angle_between_quark_and_opposite_second_lepton_theta.append(abs(theta0-theta6))
# 	angle_between_quark_and_opposite_second_lepton_theta.append(abs(theta4-theta2))
# 	angle_between_quark_and_opposite_first_lepton_phi.append(abs(phi0-phi5))
# 	angle_between_quark_and_opposite_first_lepton_phi.append(abs(phi4-phi1))
# 	angle_between_quark_and_opposite_second_lepton_phi.append(abs(phi0-phi6))
# 	angle_between_quark_and_opposite_second_lepton_phi.append(abs(phi4-phi2))
# 	angle_between_opposite_leptons_first_first_theta.append(abs(theta1-theta5))
# 	angle_between_opposite_leptons_first_first_phi.append(abs(phi1-phi5))
# 	angle_between_opposite_leptons_first_second_theta.append(abs(theta1-theta6))
# 	angle_between_opposite_leptons_first_second_theta.append(abs(theta2-theta5))


# plt.hist(angle_between_dileptons_theta, bins=20)
# plt.show()
# plt.hist(angle_between_dileptons_phi, bins=20)
# plt.show()
# plt.hist(angle_between_quark_and_first_lepton_theta, bins=20)
# plt.show()
# plt.hist(angle_between_quark_and_second_lepton_theta, bins=20)
# plt.show()
# plt.hist(angle_between_quark_and_first_lepton_phi, bins=20)
# plt.show()
# plt.hist(angle_between_quark_and_second_lepton_phi, bins=20)
# plt.show()
# plt.hist(angle_between_quark_and_opposite_first_lepton_theta, bins=20)
# plt.show()
# plt.hist(angle_between_quark_and_opposite_second_lepton_theta, bins=20)
# plt.show()
# plt.hist(angle_between_quark_and_opposite_first_lepton_phi, bins=20)
# plt.show()
# plt.hist(angle_between_quark_and_opposite_second_lepton_phi, bins=20)
# plt.show()
# plt.hist(angle_between_opposite_leptons_first_first_theta, bins=20)
# plt.show()
# plt.hist(angle_between_opposite_leptons_first_first_phi, bins=20)
# plt.show()plt.hist(angle_between_opposite_leptons_first_first_theta, bins=20)
# plt.show()
# plt.hist(angle_between_opposite_leptons_first_second_theta, bins=20)
# plt.show()

# Dilepton invariant mass
mllmax = 58 # calculated at the end of the script
print "mllmax =", mllmax
mllsq1 = []
mllsq2 = []
mllsqwrongcomb1 = [] # save list of mllsq for wrong combo to see
mllsqwrongcomb2 = [] # how often they are above threshold
for i in range(N):
	print "%f %f | %f %f" %(minkowskidot(events[i][1],events[i][2]),minkowskidot(events[i][5],events[i][6]), minkowskidot(events[i][1],events[i][5]),minkowskidot(events[i][2],events[i][6]))
	print "%s %s | %s %s" %(minkowskidot(events[i][1],events[i][2])>mllmax**2, minkowskidot(events[i][5],events[i][6])>mllmax**2, minkowskidot(events[i][1],events[i][5])>mllmax**2,minkowskidot(events[i][2],events[i][6])>mllmax**2)
	mllsq1.append(minkowskidot(events[i][1],events[i][2]))
	mllsq2.append(minkowskidot(events[i][5],events[i][6]))
	if events[i][1][0,4]*events[i][5][0,4]>0: # Check lepton signs for correct matching
		mllsqwrongcomb1.append(minkowskidot(events[i][1],events[i][5]))
		mllsqwrongcomb2.append(minkowskidot(events[i][2],events[i][6]))
	else:
		mllsqwrongcomb1.append(minkowskidot(events[i][1],events[i][6]))
		mllsqwrongcomb2.append(minkowskidot(events[i][2],events[i][5]))
print "N =",len(mllsqwrongcomb1)
counter1 = 0
counter2 = 0
countereither = 0
threshold = (mllmax+20)**2
for i in range(len(mllsqwrongcomb1)):
	if mllsqwrongcomb1[i]>threshold:
		counter1 += 1
	if mllsqwrongcomb2[i]>threshold:
		counter2 += 1
	if mllsqwrongcomb1[i]>threshold or mllsqwrongcomb2[i]>threshold:
		countereither += 1
print "Fraction failing cut (1/2/either) =", counter1/len(mllsqwrongcomb1), counter2/len(mllsqwrongcomb2), countereither/len(mllsqwrongcomb1)
# plt.hist(np.sqrt(mllsqwrongcomb1), bins=30)
# plt.hold('on')
plt.hist(np.sqrt(mllsq1),bins=30)
plt.show()
print "mllsqmax =", max(mllsq1), max(mllsq2), np.sqrt(max(mllsq1)), np.sqrt(max(mllsq2)), np.mean( [np.sqrt(max(mllsq1)), np.sqrt(max(mllsq2))] )