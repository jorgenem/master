#import stuff
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from iminuit import Minuit
import scipy.optimize as sciopt
np.random.seed(2) # set seed for reproducibility





# ====== define useful functions ========




def minkowskidot(a,b):
	# Inner product in Minkowski space
	return float(a[0,0]*b[0,0]-a[0,1]*b[0,1]-a[0,2]*b[0,2]-a[0,3]*b[0,3])
def minkowskinorm(a):
	# Inner product of a with a in Minkowski space
	return minkowskidot(a,a)

def smear(p,resolution):
	# Smears 4-momentum according to AcerDET manual
	r = np.random.randn()
	p_smeared = p * ( 1 + r * resolution / np.sqrt(np.abs(p[0,0])) )
	return p_smeared
def smear2(p):
	# Smears 4-momentum according to AcerDET manual, differently for different particle types (the true way to do it)
	pdgid = abs(p[0,5])
	r = np.random.randn()
	if pdgid == 11: #electron
		resolution = 0.12
		p_smeared = p * ( 1 + r * resolution / np.sqrt(np.abs(p[0,0])) )
	elif pdgid <= 4: #quark
		threshold = 3.2
		pabs = np.linalg.norm(p[0,1:4])
		pseudorapidity = 0.5*np.log( ( pabs + p[0,3] )/( pabs - p[0,3] ) )
		if pseudorapidity <= threshold:
			resolution = 0.50
		else:
			resolution = 1.00
		p_smeared = p * ( 1 + r * resolution / np.sqrt(np.abs(p[0,0])) )
	elif pdgid == 13: #muon
		resolution = 0.0005
		p_smeared = p / ( 1 + r * resolution / np.sqrt(p[0,1]**2 + p[0,2]**2) )
	else:
		print "Something went wrong in the smearing function."
		print "pdg = ", pdgid
		print "pseudorapidity = ", pseudorapidity
	return p_smeared







# ==== MAIN PART ====







# #import the LHE file of events.
# import xml.etree.ElementTree as ET
# tree = ET.parse('pptiluskvarkL-produksjon/20140529-1.txt')
# root = tree.getroot()

#import the Herwig .txt file of events
import sys
# file = open("on-shell_decay_squarks_at_rest_10000_events.txt",'r')
file = open("Pythia_cascade_events_no_ISR_or_FSR_20150120.log", 'r')
herwig = False
lines = file.readlines()

# Set known parameters
# SM particle masses
# u-quark and electron mass set to zero
# mquark = m1 = m5 = 0;
# mlepton1 = m2 = m3 = 0;
# mlepton2 = m6 = m7 = 0;


# True masses
MSuL = 565.312 # Mass of ~uL, ~cL
MSdL = 570.734 # Mass of ~dl, ~sL
Msquark = MZ =(MSuL+MSdL)/2.0 # mean squark mass, fit this
MN2 = MY = 180.337 # Mass of ~chi02
MseR = MX = 144.06 # Mass of ~eR, ~muR
# include left-handed sleptons? Must be done in the branchings before Herwig simulation in case
MN1 = MN = 9.70071979E+01 # Mass of ~chi01 (dark matter!)
MZprim = MZ
MYprim = MY
MXprim = MX
MNprim = MN

true_values = np.array([MZ,MY,MX,MN])

# Define normalizing mass (characteristic mass scale of the problem)
Mnorm = 100
# print "Mnorm = ", Mnorm


def minimize(Nbins, Nevents,resolution,Minitial,Mlowbound):
	# Make lists for storing D matrices and E vectors
	N = Nbins*Nevents
	Dlist = []
	Elist = []
	Adetlist = np.zeros(0)
	A_nosmeardetlist = np.zeros(0)



	# Save invariant masses for making triangle
	invariant_mass_between_c1_leptons = [] 

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
		if herwig:
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


		# DETERMINANT TEST 
		# pxmiss_nosmear = - p1[0,1] - p2[0,1] - p3[0,1] - p5[0,1] - p6[0,1] - p7[0,1]
		# pymiss_nosmear = - p1[0,2] - p2[0,2] - p3[0,2] - p5[0,2] - p6[0,2] - p7[0,2]
		# A_nosmear = 1/Mnorm*2*np.matrix([[ p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] , 0 , 0 , 0 , 0 ],
		# 				[ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
		# 				[ p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] , 0 , 0 , 0 , 0 ],
		# 				[ 0.5*pxmiss_nosmear	  ,	0		, 0		  , 0		 , 0.5*pxmiss_nosmear,0 , 0 , 0 ],
		# 				[ 0	, 0 , 0 , 0 , p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] ],
		# 				[ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
		# 				[ 0	, 0 , 0 , 0 , p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] ],
		# 				[ 0 ,0.5*pymiss_nosmear, 0 , 0 , 0 	  , 0.5*pymiss_nosmear 	, 0		  , 0 		 ]])
		# / DETERMINANT TEST




		# Smear
		# r = resolution # percent/100 momentum smearing
		# p1 = smear2(p1)
		# p2 = smear2(p2)
		# p3 = smear2(p3)
		# p5 = smear2(p5)
		# p6 = smear2(p6)
		# p7 = smear2(p7)

		# Calculate invariant masses of measured particles
		# m1 = np.sign(minkowskinorm(p1))*np.sqrt(abs(minkowskinorm(p1)))
		# m2 = np.sign(minkowskinorm(p2))*np.sqrt(abs(minkowskinorm(p2)))
		# m3 = np.sign(minkowskinorm(p3))*np.sqrt(abs(minkowskinorm(p3)))
		# m5 = np.sign(minkowskinorm(p5))*np.sqrt(abs(minkowskinorm(p5)))
		# m6 = np.sign(minkowskinorm(p6))*np.sqrt(abs(minkowskinorm(p6)))
		# m7 = np.sign(minkowskinorm(p7))*np.sqrt(abs(minkowskinorm(p7)))
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

		# Calculate invariant mass between leptons for triangle plotting
		invmass_c1leptons = np.sqrt(minkowskinorm(p2+p3))   # calculate invariant 
		invariant_mass_between_c1_leptons.append(invmass_c1leptons) # mass between leptons in chain 1




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
						 [ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 ]])
		A = A/Mnorm # normalize A
		# print A
		# print np.linalg.det(A)
		#A inverse
		Ainv = A.I
		# print Ainv

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
		# print D
		E = np.dot(Ainv,C.T)
		# print E

		# store D and E
		Dlist.append(D)
		Elist.append(E)

		# Store determinants of A w&w/o smearing
		# Adetlist = np.append(Adetlist, np.linalg.det(A))
		# A_nosmeardetlist = np.append(A_nosmeardetlist, np.linalg.det(A_nosmear))


	 # ===========
	 # From here on we can forget about the event momenta, everything
	 # is stored in Dn and En for each event. Time to guess the masses.
	 # ===========




	# # Now to make a mass hypothesis (guess the correct one)
	# MZ = 5.45421001e+02 # Mass of ~uL
	# MY = 1.80337030e+02 # Mass of ~chi02
	# MX = 1.44059825e+02 # Mass of ~eR
	# MN = 9.70071979e+01 # Mass of ~chi01 (dark matter!)
	# MZprim = MZ
	# MYprim = MY
	# MXprim = MX
	# MNprim = MN
	# M = np.matrix([ MZ**2 , MY**2 , MX**2 , MN**2 , MZprim**2 , MYprim**2 , MXprim**2 , MNprim**2 ])
	# M = M/Mnorm
	# print M

	# # # Calculate the "chi-squared" error of the correct hypothesis
	# P = [] # store Pn
	# xisquared = 0
	# offshell = [] # list to store p4nsquared - MN**2
	# for n in range(N):
	# 	Pn = np.dot(Dlist[n],M.T) + Elist[n]
	# 	P.append(Pn) #store in case needed

	# 	#print "hei", Pn.shape

	# 	p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
	# 	p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

	# 	xisquared +=  (p4nsquared - MN**2)**2 + (p8nsquared - MNprim**2)**2
	# 	offshell.append(abs(p4nsquared-MN**2))
	# 	offshell.append(abs(p8nsquared-MNprim**2))

	# print "xisquared", xisquared/Mnorm**2
	# print np.mean(offshell)



	# === Plot some results ===

	# plt.hist(invariant_mass_between_c1_leptons, bins=100)
	# plt.title('Distribution of lepton pair invariant mass in %d events.' % N)
	# plt.xlabel(r'$m_{l^+l^-}$')
	# plt.ylabel('Occurences')
	# plt.show()


	# Plot det(A) scatter smear vs nosmear
	# plt.loglog(abs(A_nosmeardetlist), abs(Adetlist), 'bx')
	# plt.hold('on')
	# plt.plot([min(abs(Adetlist)),max(abs(Adetlist))], [min(abs(Adetlist)),max(abs(Adetlist))], 'r')	
	# plt.title('Per-event det(A) unsmeared vs smeared, smearing resolution %.1f' % resolution)
	# plt.xlabel(r'$\mathrm{det}(A_\mathrm{non-smear})$')
	# plt.ylabel(r'$\mathrm{det}(A_\mathrm{smear})$')
	# plt.show()

	# Plot quark invariant masses
	# quarkmasses = np.concatenate((quark1mass,quark2mass),axis=0)
	# plt.hist(quarkmasses[abs(quarkmasses[:,0])==4,1], bins=100)
	# plt.title('Distribution of c/cbar quark invariant masses before parton showering, %d events' % (N) )
	# plt.xlabel(r'$m_q^\mathrm{inv}$ [GeV]')
	# plt.show()

	# print len(quark1mass[quark1mass[:,0]>=5,1])
	# print len(quark1mass[abs(quark1mass[:,0])==1,1])+len(quark1mass[abs(quark1mass[:,0])==2,1])+len(quark1mass[abs(quark1mass[:,0])==3,1])+len(quark1mass[abs(quark1mass[:,0])==4,1])




	# ============ Minimization to best fit =================


	# Define xi-squared function to minimize
	# def xisquared(MZ, MY, MX, MN, MZp, MYp, MXp, MNp):
	# 	# Set up Webber's M vector
	# 	M = np.matrix([ MZ**2 , MY**2 , MX**2 , MN**2 , MZp**2 , MYp**2 , MXp**2 , MNp**2 ])
	# 	# Calculate the "chi-squared" error of the hypothesis
	# 	P = [] # store Pn
	# 	xisquared = 0
	# 	offshell = [] # list to store p4nsquared - MN**2
	# 	for n in range(N):
	# 		Pn = np.dot(Dlist[n],M.T) + Elist[n]
	# 		P.append(Pn) #store in case needed
		
		
	# 		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
	# 		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2
		
	# 		xisquared +=  (p4nsquared - MN**2)**2 + (p8nsquared - MNprim**2)**2
	# 		offshell.append(abs(p4nsquared-MN**2))
	# 		offshell.append(abs(p8nsquared-MNprim**2))
		
	# 	return xisquared

	# xi-squared function to minimize with identical chains
	def xisquared_identical_chains(Masses, Nevents, i): #, MZp, MYp, MXp, MNp):
		Nevents = int(Nevents)
		i = int(i)
		# Duplicate masses for primed chain
		MZp, MYp, MXp, MNp = MZ, MY, MX, MN = Masses
		# Set up Webber's M vector
		M = np.matrix([ MZ**2 , MY**2 , MX**2 , MN**2 , MZp**2 , MYp**2 , MXp**2 , MNp**2 ])
		M = M/Mnorm**2 #normalise M
		# Calculate the "chi-squared" error of the hypothesis
		P = [] # store Pn
		xisquared = 0
		# offshell = [] # list to store p4nsquared - MN**2
		for n in range(i*Nevents, (i+1)*Nevents):
			Pn = np.dot(Dlist[n],M.T) + Elist[n]
			# print Pn
			# P.append(Pn) #store in case needed
		
		
			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2
			# print p4nsquared,p8nsquared
			# print M[0,3]/Mnorm**2
		
			xisquared +=  (p4nsquared - M[0,3])**2 + (p8nsquared - M[0,7])**2 # p4/p8 and M is normalized by Mnorm.

			# offshell.append(abs(p4nsquared-MN**2))
			# offshell.append(abs(p8nsquared-MNprim**2))
		xisquared = xisquared

		return xisquared

	# print "calling test"
	# Mtest=[447.429941541876, 115.38434233639869, 94.00228334374647, 48.25627489449959]
	# print xisquared_identical_chains(Mtest,1,0)
	best_fit = np.zeros((Nbins,6))
	relative_fit_error = np.zeros((Nbins,4))

	# TEST to see if Mathematica has the exact same xisquared implementation. It does. (in minimization_v1.nb)	
	# for i in range(10):
	# 	print xisquared_identical_chains(true_values, 1, i)

	for i in range(Nbins):

		# m = Minuit(xisquared_identical_chains, 
		# 		MZ=Minitial[0], #limit_MZ=(300, 700),
		# 		MY=Minitial[1], #limit_MY=(100, 300),
		# 		MX=Minitial[2], #limit_MX=(100, 200), #err_MX=1000, 
		# 		MN=Minitial[3], # err_MN=10, limit_MN=(50, 150),
		# 		print_level=0,
		# 		#maxcalls=None,
		# 		error_MZ=Minitial[0]*0.01, 
		# 		error_MY=Minitial[1]*0.01, 
		# 		error_MX=Minitial[2]*0.01, 
		# 		error_MN=Minitial[3]*0.01,
		# 		Nevents=Nevents, fix_Nevents=True,
		# 		i=i, fix_i = True,
		# 		# strategy=2, # doesn't work? doesn't help.
		# 		# tol=1e-3,
		# 		# up=0.5
		# 		)
		# #m.printMode = 1
		# m.migrad()

		# # true_values = [MZ, MY, MX, MN]
		# best_fit[i,:] = m.values['MZ'], m.values['MY'], m.values['MX'], m.values['MN'], m.ncalls, m.fval
		# relative_fit_error[i,:] = [(MZ-m.values['MZ'])/MZ, (MY-m.values['MY'])/MY, (MX-m.values['MX'])/MX, (MN-m.values['MN'])/MN]

		# Scipy minimization
		m = sciopt.minimize(xisquared_identical_chains, Minitial, 
						  args=(Nevents, i), method='Nelder-Mead', 
						  #bounds=((Mlowbound[0], None), (Mlowbound[1], None), (Mlowbound[2], None), (Mlowbound[3], None)),
						  #tol=1e-40,
						  options={'maxfev': 2000,'xtol': 1e-8,
						  		   'ftol':1e-8, 'disp': True}
						  )
		best_fit[i,:] = m.x[0], m.x[1], m.x[2], m.x[3], m.nfev, m.fun
		relative_fit_error = 0
	
	return best_fit, relative_fit_error


# ==== Run: ======


# Initialize run
Nevents = 25
Nbins = 10
# mass_offset = 0.99
# Minitial = [5.5e2, 1.8e2, 1.5e2, 1e2, 5.5e2, 1.8e2, 1.5e2, 1e2] # Starting point for parameter scan. 
  # Make all mass guesses be equally far off, percentage-wise.
M_explowbound=[400,94,94,46]


smearing_resolution = 0

plot_counter = 1
for mass_offset in [1]:
	Minitial = true_values*np.array([2-mass_offset,mass_offset,mass_offset,mass_offset])
	# Minitial=np.array([447.429941541876, 115.38434233639869, 94.00228334374647, 48.25627489449959])*mass_offset
	print Minitial
	best_fit, relative_fit_error = minimize(Nbins, Nevents, smearing_resolution, Minitial,M_explowbound)
	for i in range(Nbins):
		# print "%3d % .6e   % .6e   % .6e   % .6e   %3d   % .6e" %(i+1, best_fit[i,0], best_fit[i,1], best_fit[i,2], best_fit[i,3], best_fit[i,4], best_fit[i,5])
		print "%3d %2.1f   %2.1f   %2.1f   %2.1f   %3d   % .6e" %(i+1, best_fit[i,0], best_fit[i,1], best_fit[i,2], best_fit[i,3], best_fit[i,4], best_fit[i,5])

	# Get true mass values
	Msquark = true_values[0]
	Mchi2 = true_values[1]
	Mslepton = true_values[2]
	Mchi1 = true_values[3]

	# Extra masses
	# MsquarkuL = 5.61119014E+02
	# MsquarkdL = 5.68441109E+02
	# MsquarksL = 5.68441109E+02
	# MsquarkcL = 5.61119014E+02
	# MsquarkuR = 3.00000000E+04
	# MsquarkdR = 3.00000000E+04
	# MsquarksR = 3.00000000E+04
	# MsquarkcR = 3.00000000E+04

	# Take out best-fit values for each bin as vectors (note minuscle m for the fit-vector)
	msquark = best_fit[:,0]
	mchi2 = best_fit[:,1]
	mslepton = best_fit[:,2]
	mchi1 = best_fit[:,3]

	msquark_passcut = []
	mchi2_passcut = []
	mslepton_passcut = []
	mchi1_passcut = []
	cut = 1e-9 # xi^2 cut value in units of (100 GeV)^4
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

	print "Mass offset =", mass_offset
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
	ylim=[300,700]
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
	plt.title("Mass offset = %.2f, Nevents = %d, cut = %.2e, passcut-fraction = %.3f"%(mass_offset,Nevents,cut,len(msquark_passcut)/float(Nbins)))
	plt.text(50,MZ+5,r'$\tilde q$',fontsize=20)
	plt.text(MY+1,320,r'$\tilde\chi_2^0$',fontsize=20)
	plt.text(MX+1,320,r'$\tilde l$',fontsize=20)
	plt.text(MN+1,320,r'$\tilde \chi_1^0$',fontsize=20)
	# plt.savefig('100_bins_25_events_pythia_events_nelder-mead_%1.2f_initial_guess_no_ISR_or_FSR_xisquared-cut.pdf'%mass_offset, format='pdf')

	# plt.hold('off')
	# plt.close()

	plot_counter += 1

plt.show()



	# print 'smearing', "%2.2f ," % smearing_resolution
	# print 'True masses', true_values
	# print 'Best-fit values', best_fit
	# print 'relative_fit_error', relative_fit_error, ', abs mean fit error', "%.2e" %np.mean(np.abs(relative_fit_error))
	# print "number of runs =", len(xisquaredlist), ", mean xi^2 =", np.mean(xisquaredlist), "final xi^2 =", xisquaredlist[-1]

# Minitial = [5.5e2, 1.8e2, 1.5e2, 1e2, 5.5e2, 1.8e2, 1.5e2, 1e2] # Starting point for parameter scan
# Nlist = range(100,1000, 50)
# relative_fit_error_list = []
# for N in Nlist:
# 	true_values, best_fit, relative_fit_error = minimize(N, 0, Minitial)
# 	relative_fit_error_list.append(relative_fit_error[0])

# plt.plot(relative_fit_error_list)
# plt.show()