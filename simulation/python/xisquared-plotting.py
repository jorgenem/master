#import stuff
import numpy as np
import matplotlib.pyplot as plt
from math import pi
np.random.seed(2) # set seed for reproducibility





# ====== define useful functions ========



def proj(v,u):
	# projects v onto u
	return np.dot(u,np.transpose(v))/float(np.dot(u,u.T)) * u
def minkowskidot(a,b):
	# Inner product in Minkowski space
	return float(a[0,0]*b[0,0]-a[0,1]*b[0,1]-a[0,2]*b[0,2]-a[0,3]*b[0,3])
def minkowskinorm(a):
	# Inner product of a with a in Minkowski space
	return minkowskidot(a,a)

def smear(p,resolution):
	# Smears 4-momentum according to AcerDET manual
	r = np.random.randn()
	p_smeared = p * ( 1 + r * resolution / np.sqrt(p[0,0]) )
	return p_smeared


generator=0
generators = ["Herwig++", "Pythia-full", "Pythia-noIFSR","Simple","Simple-14TeV"]

Nevents=25
bin_number=2


#import the Herwig .txt file of events
import sys
files = ["../events/herwigpp-9563-events-complete-momcons-20150314_only_OFL.dat",
		 "Pythia_cascade_events_no_ISR_or_FSR_20150120.dat",
		 "on-shell_decay_squarks_at_rest_2500_events_possibly_corrected.dat",
		 "on-shell_decay_squarks_with_pz_14TeV-CoM_2500_events_possibly_corrected.dat"]
file = open(files[generator],'r')
lines = file.readlines()




# True masses
MSuL = 565.312 # Mass of ~uL, ~cL
MSdL = 570.734 # Mass of ~dl, ~sL
Msquark = MZ =(MSuL+MSdL)/2 # mean squark mass, fit this
Mchi2 = MY = 180.337 # Mass of ~chi02
Mslepton = MX = 144.06 # Mass of ~eR, ~muR
# include left-handed sleptons? Must be done in the branchings before Herwig simulation in case
Mchi1 = MN = 9.70071979E+01 # Mass of ~chi01 (dark matter!)


resolution = 0 # smearing resolution, 0 means no smearing




N = Nevents*bin_number
Dlist = []
Elist = []
Adetlist = np.zeros(0)
A_nosmeardetlist = np.zeros(0)

# Define normalizing mass (characteristic mass scale of the problem)
Mnorm = 100 # EW scale
# print "Mnorm = ", Mnorm


# Save quark invariant masses
quark1mass = np.zeros((N,2));
quark2mass = np.zeros((N,2));

# N - How much loop?
for i in range(N):
	# Loop over events to get 4-vectors for each particle for each event. 
	# Particles are numbered according to Webber (arXiv:0907.5307v2) fig. 1
	# (the lepton/antilepton ordering is arbitrary in each chain, the lepton has been 
	# chosen as 2/6 and the antilepton as 3/7)

	# Read all particles from file
	# chain 1
	quark1 = lines[9*i + 1].split()
	p1 = np.matrix([ float(quark1[4]), float(quark1[1]), float(quark1[2]), float(quark1[3]), float(quark1[5]), int(quark1[0]) ])
	lepton11 = lines[9*i + 2].split()
	p2 = np.matrix([ float(lepton11[4]), float(lepton11[1]), float(lepton11[2]), float(lepton11[3]), float(lepton11[5]), int(lepton11[0]) ])
	lepton12 = lines[9*i + 3].split()
	p3 = np.matrix([ float(lepton12[4]), float(lepton12[1]), float(lepton12[2]), float(lepton12[3]), float(lepton12[5]), int(lepton12[0]) ])
	neutralino1 = lines[9*i + 4].split()
	p4 = np.matrix([ float(neutralino1[4]), float(neutralino1[1]), float(neutralino1[2]), float(neutralino1[3]), float(neutralino1[5]), int(neutralino1[0]) ])
	#chain2
	quark2 = lines[9*i + 5].split()
	p5 = np.matrix([ float(quark2[4]), float(quark2[1]), float(quark2[2]), float(quark2[3]), float(quark2[5]), int(quark2[0]) ])
	lepton21 = lines[9*i + 6].split()
	p6 = np.matrix([ float(lepton21[4]), float(lepton21[1]), float(lepton21[2]), float(lepton21[3]), float(lepton21[5]), int(lepton21[0]) ])
	lepton22 = lines[9*i + 7].split()
	p7 = np.matrix([ float(lepton22[4]), float(lepton22[1]), float(lepton22[2]), float(lepton22[3]), float(lepton22[5]), int(lepton22[0]) ])
	neutralino2 = lines[9*i + 8].split()
	p8 = np.matrix([ float(neutralino2[4]), float(neutralino2[1]), float(neutralino2[2]), float(neutralino2[3]), float(neutralino2[5]), int(neutralino2[0]) ])

	# Take care of units - Herwig++ likes MeV, we like GeV (avoid disturbing the pdg code entry)
	if generators[generator]=="Herwig":
		p1[0,0:5] /= 1000
		p2[0,0:5] /= 1000
		p3[0,0:5] /= 1000
		p4[0,0:5] /= 1000
		p5[0,0:5] /= 1000
		p6[0,0:5] /= 1000
		p7[0,0:5] /= 1000
		p8[0,0:5] /= 1000

	# Save stuff for plotting
	quark1mass[i,0] = p1[0,5]
	quark1mass[i,1] = p1[0,4]
	quark2mass[i,0] = p5[0,5]
	quark2mass[i,1] = p5[0,4]


	print p1


	# DETERMINANT TEST 
	pxmiss_nosmear = - p1[0,1] - p2[0,1] - p3[0,1] - p5[0,1] - p6[0,1] - p7[0,1]
	pymiss_nosmear = - p1[0,2] - p2[0,2] - p3[0,2] - p5[0,2] - p6[0,2] - p7[0,2]
	A_nosmear = 1/Mnorm*2*np.matrix([[ p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] , 0 , 0 , 0 , 0 ],
					[ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
					[ p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] , 0 , 0 , 0 , 0 ],
					[ 0.5*pxmiss_nosmear	  ,	0		, 0		  , 0		 , 0.5*pxmiss_nosmear,0 , 0 , 0 ],
					[ 0	, 0 , 0 , 0 , p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] ],
					[ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
					[ 0	, 0 , 0 , 0 , p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] ],
					[ 0 ,0.5*pymiss_nosmear, 0 , 0 , 0 	  , 0.5*pymiss_nosmear 	, 0		  , 0 		 ]])
	# / DETERMINANT TEST




	# Smear, r percent resolution
	r = resolution # percent/100 momentum smearing
	# p1 = smear(p1,r)
	# p2 = smear(p2,r)
	# p3 = smear(p3,r)

	# p5 = smear(p5,r)
	# p6 = smear(p6,r)
	# p7 = smear(p7,r)

	# Calculate invariant masses of measured particles
	m1square = minkowskinorm(p1)
	m2square = minkowskinorm(p2)
	m3square = minkowskinorm(p3)
	m5square = minkowskinorm(p5)
	m6square = minkowskinorm(p6)
	m7square = minkowskinorm(p7)

	# Check invariant mass of initial colliding partons?
	#print minkowskinorm(p1+p2+p3+p4+p5+p6+p7+p8)
	# Check that the invariant mass of particles is close to shell mass
	# print minkowskinorm(p1) - m1**2
	# print minkowskinorm(p2) - m2**2
	# print minkowskinorm(p3) - m3**2
	# print minkowskinorm(p4) - MN**2
	# print minkowskinorm(p5) - m5**2
	# print minkowskinorm(p6) - m6**2
	# print minkowskinorm(p7) - m7**2
	# print minkowskinorm(p8) - MNprim**2

	# Check if invariant mass of decays match mass of decaying
	# print "p3+p4 ", np.sqrt(abs(minkowskinorm(p3+p4) - MX**2))
	# print "p2+p3+p4 ", np.sqrt(abs(minkowskinorm(p2+p3+p4) - MY**2))
	# print "p1+p2+p3+p4 ", np.sqrt(abs(minkowskinorm(p1+p2+p3+p4) - MZ**2))
	# print "p7+p8 ", np.sqrt(abs(minkowskinorm(p7+p8) - MXprim**2))
	# print "p6+p7+p8 ", np.sqrt(abs(minkowskinorm(p2+p3+p4) - MYprim**2))
	# print "p5+p6+p7+p8 ", np.sqrt(abs(minkowskinorm(p1+p2+p3+p4) - MZprim**2))





	# ==== Define Webber's stuff ====

	# need the pxmiss and pymiss, taken from the actual neutralino transverse momenta 
	# (this is cheating, of course)
	# pxmiss = p4[0,1]+p8[0,1]
	# pymiss = p4[0,2]+p8[0,2]

	# Calculate missing transverse from (smeared) visible particles
	pxmiss = - p1[0,1] - p2[0,1] - p3[0,1] - p5[0,1] - p6[0,1] - p7[0,1]
	pymiss = - p1[0,2] - p2[0,2] - p3[0,2] - p5[0,2] - p6[0,2] - p7[0,2]

	# print "pxmiss", pxmisstrue - pxmiss
	# print "pymiss", pymisstrue - pymiss

	#A matrix
	A = 2*np.matrix([[ p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] , 0 , 0 , 0 , 0 ],
					[ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
					[ p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] , 0 , 0 , 0 , 0 ],
					[ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
					[ 0	, 0 , 0 , 0 , p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] ],
					[ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
					[ 0	, 0 , 0 , 0 , p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] ],
					[ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]])
	A = A/Mnorm # normalize A
	# print np.linalg.det(A)
	#A inverse
	# print A
	Ainv = A.I

	#B matrix
	B = np.matrix([[-1,1,0,0,0,0,0,0],
				   [0,-1,1,0,0,0,0,0],
				   [0,0,-1,1,0,0,0,0],
				   [0,0,0,0,0,0,0,0],
				   [0,0,0,0,-1,1,0,0],
				   [0,0,0,0,0,-1,1,0],
				   [0,0,0,0,0,0,-1,1],
				   [0,0,0,0,0,0,0,0]])

	#C vector
	C = np.matrix([ 2*minkowskidot(p1,p2) + 2*minkowskidot(p1,p3) + m1square,
					2*minkowskidot(p2,p3) + m2square,
					m3square,
					pxmiss**2,
					2*minkowskidot(p5,p6) + 2*minkowskidot(p5,p7) + m5square,
					2*minkowskidot(p6,p7) + m6square,
					m7square,
					pymiss**2])
	C = C/Mnorm**2 # normalize C
	# print C

	# Composite matrix & vector D and E
	D = np.dot(Ainv,B)
	E = np.dot(Ainv,C.T)

	# store D and E
	Dlist.append(D)
	Elist.append(E)

	# Store determinants of A w&w/o smearing
	Adetlist = np.append(Adetlist, np.linalg.det(A))
	A_nosmeardetlist = np.append(A_nosmeardetlist, np.linalg.det(A_nosmear))




# xi-squared function to minimize with identical chains
def xisquared_identical_chains(MZ, MY, MX, MN, Nevents, i): #, MZp, MYp, MXp, MNp):
	Nevents = int(Nevents)
	i = int(i)
	# Duplicate masses for primed chain
	MZp, MYp, MXp, MNp = MZ, MY, MX, MN# = Masses
	# print Masses
	# Set up Webber's M vector
	M = np.matrix([ MZ**2 , MY**2 , MX**2 , MN**2 , MZp**2 , MYp**2 , MXp**2 , MNp**2 ])
	M = M/Mnorm**2 #normalise M
	print "np.size(M)",np.size(M)
	# Calculate the "chi-squared" error of the hypothesis
	P = [] # store Pn
	xisquared = 0
	# offshell = [] # list to store p4nsquared - MN**2
	for n in range(i*Nevents, (i+1)*Nevents):
		Pn = np.dot(Dlist[n],M.T) + Elist[n]
		P.append(Pn) #store in case needed
	
		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2
	
		xisquared +=  (p4nsquared - M[0,3])**2 + (p8nsquared - M[0,7])**2 # p4/p8 is normalized by Mnorm.
		# xisquared +=  (p4nsquared - (MN/Mnorm)**2)**2 + (p8nsquared - (MN/Mnorm)**2)**2 # p4/p8 is normalized by MN.

		# offshell.append(abs(p4nsquared-MN**2))
		# offshell.append(abs(p8nsquared-MNprim**2))
	xisquared = xisquared
	return xisquared





# Plot xi^2 as function of some masses to see how bumpy
from mpl_toolkits.mplot3d import Axes3D

true_offset = 1

minm = 0.0
maxm = 3
Nlinspace = 300

msquark_linspace = np.linspace(Msquark*minm, Msquark*maxm, Nlinspace)
mchi2_linspace   = np.linspace(Mchi2*minm, Mchi2*maxm, Nlinspace)
mslepton_linspace = np.linspace(Mslepton*minm, Mslepton*maxm, Nlinspace)
mchi1_linspace = np.linspace(Mchi1*minm, Mchi1*maxm, Nlinspace)
msquark_mesh1, mchi2_mesh = np.meshgrid(msquark_linspace, mchi2_linspace)
msquark_mesh2, mslepton_mesh = np.meshgrid(msquark_linspace, mslepton_linspace)
msquark_mesh3, mchi1_mesh = np.meshgrid(msquark_linspace, mchi1_linspace)
# mslepton_mesh2, mchi1_mesh2 = np.meshgrid(mslepton_linspace, mchi1_linspace)

xi2_plot_squarkchi2 = np.log(xisquared_identical_chains(msquark_mesh1, mchi2_mesh, true_offset*Mslepton, true_offset*Mchi1, Nevents,bin_number-1))
xi2_plot_squarkslepton = np.log(xisquared_identical_chains(msquark_mesh2, true_offset*Mchi2, mslepton_mesh, true_offset*Mchi1, Nevents, bin_number-1))
xi2_plot_squarkchi1 = np.log(xisquared_identical_chains(msquark_mesh3, true_offset*Mchi2, true_offset*Mslepton, mchi1_mesh, Nevents, bin_number-1))
# xi2_plot_sleptonchi1 = np.log(xisquared_identical_chains(480, 150, mslepton_mesh2, mchi1_mesh2, Nevents, bin_number-1))

# Plot 1: squark-chi2
fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d', )
ax.set_zscale(u'linear')
ax.plot_wireframe(msquark_mesh1, mchi2_mesh, xi2_plot_squarkchi2, rstride=10, cstride=10, color='k')
plt.title(r'Generator = %s, $\tilde q - \tilde \chi_2^0$'%generators[generator],{'fontsize':16})

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde \chi_2^0}$', {'fontsize':20})
ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

plt.savefig('xisquared_surf_plots/20150123_xisquared-surface_%s_squark-chi2.pdf'%generator, format='pdf')


# plt.show()


# Plot 2: squark-slepton
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d', )
ax.set_zscale(u'linear')
ax.plot_wireframe(msquark_mesh2, mslepton_mesh, xi2_plot_squarkslepton, rstride=10, cstride=10, color='k')
plt.title(r'Generator = %s, $\tilde q - \tilde l$'%generators[generator],{'fontsize':16})

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde l}$', {'fontsize':20})
ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

plt.savefig('xisquared_surf_plots/20150123_xisquared-surface_%s_squark-slepton.pdf'%generator, format='pdf')
# plt.show()


# Plot 3: squark-chi1
fig = plt.figure(3)
ax = fig.add_subplot(111, projection='3d', )
ax.set_zscale(u'linear')
ax.plot_wireframe(msquark_mesh3, mchi1_mesh, xi2_plot_squarkchi1, rstride=10, cstride=10, color='k')
plt.title(r'Generator = %s, $\tilde q - \tilde \chi_1^0$'%generators[generator],{'fontsize':16})

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde \chi_1^0}$', {'fontsize':20})
ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

plt.savefig('xisquared_surf_plots/20150123_xisquared-surface_%s_squark-chi1.pdf'%generator, format='pdf')
plt.show()


# # Plot 4: slepton-chi1
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d', )
# ax.set_zscale(u'linear')
# ax.plot_wireframe(mslepton_mesh2, mchi1_mesh2, xi2_plot_sleptonchi1, rstride=10, cstride=10, color='k')
# # plt.title('test')

# ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
# ax.set_ylabel(r'$m_{\tilde \chi_1^0}$', {'fontsize':20})
# ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

# plt.show()

# Minitial = np.array([MZ, MY, MX, MN])*1.05
# print Minitial
# import scipy.optimize as sciopt
# minimum = sciopt.minimize(xisquared_identical_chains, Minitial, 
# 						  args=(25,0), method='SLSQP', 
# 						  #bounds=((MZ*0.5, MZ*2), (MY*0.5, MY*2), (MX*0.5, MX*2), (MN*0.0, MN*2))
# 						  )
# print minimum