#import stuff
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from iminuit import Minuit
import scipy.optimize as sciopt
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
	return np.matrix([p_smeared[0,0], p_smeared[0,1], p_smeared[0,2], p_smeared[0,3], p_smeared[0,4], p[0,5]])

# def xisquared_identical_chains_with_combinatorics(Masses, Ainv_list, Clist, Nevents, i, Mnorm): 
# 	Nevents = int(Nevents)
# 	i = int(i)
# 	# Duplicate masses for primed chain
# 	MZp, MYp, MXp, MNp = MZ, MY, MX, MN = Masses
# 	# Set up Webber's M vector
# 	M = np.matrix([ MZ**2 , MY**2 , MX**2 , MN**2 , MZp**2 , MYp**2 , MXp**2 , MNp**2 ])
# 	M = M/Mnorm**2 #normalise M

# 	#B matrix
# 	B = np.matrix([[-1,1,0,0,0,0,0,0],
# 				   [0,-1,1,0,0,0,0,0],
# 				   [0,0,-1,1,0,0,0,0],
# 				   [0,0,0,0,0,0,0,0],
# 				   [0,0,0,0,-1,1,0,0],
# 				   [0,0,0,0,0,-1,1,0],
# 				   [0,0,0,0,0,0,-1,1],
# 				   [0,0,0,0,0,0,0,0]])

# 	# Calculate the "chi-squared" error of the hypothesis
# 	xisquared = 0
# 	for n in range(i*Nevents, (i+1)*Nevents):
# 		xisquared_current_tmp = []


# 		for m in range(combinations):
# 			Ainv_current_permuted = np.dot(permute[m],Ainv_list[n])
# 			D = np.dot(Ainv_current_permuted, B)
# 			E = np.dot(Ainv_current_permuted, Clist[n].T)

# 			Pn = np.dot(D, M) + E
# 			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
# 			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2
# 			xisquared +=  (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 # p4/p8 is normalized by Mnorm.

# 	xisquared = xisquared/(float(Nevents))
# 	return xisquared


def xisquared_identical_chains_with_combinatorics(Masses, D_lists, E_lists, Nevents, j, all_leptons_equal_list, Mnorm):
	Nevents = int(Nevents)
	j = int(j)

	# Nx = 300
	# Ny = 300



	
	# Import masses
	MZ, MY, MX, MN = Masses
	# Set up the M vector of current mass guess
	M = np.matrix([MZ**2, MY**2, MX**2, MN**2, MZ**2, MY**2, MX**2, MN**2])/float(Mnorm**2)
	M = M.T

	best_fit_combination = [] # list of best-fit combination for each event
	# TODO: Implement a check of how many non-true best-fit combinations there are for different mass points
	#best_fit_combination = np.array(Nevents,Nx,Ny) # list of best-fit combination for each event
	xisquared = 0 # begin accumulation of the xisquared sum
	for n in range(j*Nevents, (j+1)*Nevents):

		xisquared_current_tmp = []
		# First 8 combinations are common to all cases

		# Case 1 - A1/C1

		# Subcase 1 - no permutations
		Pn = D_lists[0][n] * M + E_lists[0][n]

		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

		xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

		# Subcase 2 - permute 23
		Pn = D_lists[1][n] * M + E_lists[1][n]

		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

		xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

		# Subcase 3 - permute 67
		Pn = D_lists[2][n] * M + E_lists[2][n]

		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

		xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

		# Subcase 4 - permute 23 and 67
		Pn = D_lists[3][n] * M + E_lists[3][n]

		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

		xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )


		# Case 2 - A2/C2
		# Subcase 1 - no permutations
		Pn = D_lists[4][n] * M + E_lists[4][n]

		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

		xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

		# Subcase 2 - permute 23
		Pn = D_lists[5][n] * M + E_lists[5][n]

		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

		xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

		# Subcase 3 - permute 67
		Pn = D_lists[6][n] * M + E_lists[6][n]

		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

		xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

		# Subcase 4 - permute 23 and 67
		Pn = D_lists[7][n] * M + E_lists[7][n]

		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

		xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )


		# Must check all-leptons-equal on a per-event basis
		if all_leptons_equal_list[n]:
			# Case 3 - A3/C3
			# Subcase 1 - no permutations
			Pn = D_lists[8][n] * M + E_lists[8][n]

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

			# Subcase 2 - permute 23
			Pn = D_lists[9][n] * M + E_lists[9][n]

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

			# Subcase 3 - permute 67
			Pn = D_lists[10][n] * M + E_lists[10][n]

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

			# Subcase 4 - permute 23 and 67
			Pn = D_lists[11][n] * M + E_lists[11][n]

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )


			# Case 4 - A4/C4
			# Subcase 1 - no permutations
			Pn = D_lists[12][n] * M + E_lists[12][n]

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

			# Subcase 2 - permute 23
			Pn = D_lists[13][n] * M + E_lists[13][n]

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

			# Subcase 3 - permute 67
			Pn = D_lists[14][n] * M + E_lists[14][n]

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

			# Subcase 4 - permute 23 and 67
			Pn = D_lists[15][n] * M + E_lists[15][n]

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp.append( (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 )

			

		# END IF all_leptons_equal


		# Check which combination gives the smallest xisquared contribution
		xisquared_current = np.min(xisquared_current_tmp,0)
		best_fit_combination.append(xisquared_current_tmp.index(xisquared_current))
		if xisquared_current_tmp.count(xisquared_current) > 1:
			print "Warning: multiple best-fit combinations for event ", n

		# Add best-fit combination xisquared value to total
		xisquared += xisquared_current
	# END loop over events/n
	xisquared = xisquared
	return xisquared#, best_fit_combination
	# END xisquared definition


def best_fit(Nbins, Nevents, Mtrue, Minitial, Mnorm, file):
	# Make lists for storing D matrices and E vectors
	N = Nbins*Nevents
	Plist = []
	D11_list = []
	D12_list = []
	D13_list = []
	D14_list = []
	D21_list = []
	D22_list = []
	D23_list = []
	D24_list = []
	D31_list = []
	D32_list = []
	D33_list = []
	D34_list = []
	D41_list = []
	D42_list = []
	D43_list = []
	D44_list = []
	E11_list = []
	E12_list = []
	E13_list = []
	E14_list = []
	E21_list = []
	E22_list = []
	E23_list = []
	E24_list = []
	E31_list = []
	E32_list = []
	E33_list = []
	E34_list = []
	E41_list = []
	E42_list = []
	E43_list = []
	E44_list = []
	all_leptons_equal_list = []

	#import the Herwig .txt file of events
	import sys

	herwig = False
	lines = file.readlines()

	# Save invariant masses for making triangle
	# invariant_mass_between_c1_leptons = [] 

	# Save quark invariant masses
	# quark1mass = np.zeros((N,2));
	# quark2mass = np.zeros((N,2));

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

		# # Take care of units - Herwig++ likes MeV, we like GeV (avoid disturbing the pdg code entry)
		# if herwig:
		# 	p1[0,0:5] /= 1000
		# 	p2[0,0:5] /= 1000
		# 	p3[0,0:5] /= 1000
		# 	p4[0,0:5] /= 1000
		# 	p5[0,0:5] /= 1000
		# 	p6[0,0:5] /= 1000
		# 	p7[0,0:5] /= 1000
		# 	p8[0,0:5] /= 1000

		# Smear
		# p1 = smear2(p1)
		# p2 = smear2(p2)
		# p3 = smear2(p3)
		# p5 = smear2(p5)
		# p6 = smear2(p6)
		# p7 = smear2(p7)

		# Calculate invariant masses of measured particles after smearing and replace
		m1squared = p1[0,4] = minkowskinorm(p1)
		m2squared = p2[0,4] = minkowskinorm(p2)
		m3squared = p3[0,4] = minkowskinorm(p3)
		m5squared = p5[0,4] = minkowskinorm(p5)
		m6squared = p6[0,4] = minkowskinorm(p6)
		m7squared = p7[0,4] = minkowskinorm(p7)

		# Calculate missing transverse from (smeared) visible particles
		pxmiss = p4[0,1] + p8[0,1]
		pymiss = p4[0,2] + p8[0,2]

		# P = np.matrix(	[p1],
		# 				[p2],
		# 				[p3],
		# 				[p4],
		# 				[p5],
		# 				[p6],
		# 				[p7],
		# 				[p8]	)

		# Plist.append(P)


		# From here on we must take combinatorics into account and evaluate all possible combinations. We do it right away for each event and store all possibilities.


		# Define 8x8 permutation matrices
		permute23 = 		np.matrix([	[1,0,0,0,0,0,0,0],
										[0,0,1,0,0,0,0,0],
										[0,1,0,0,0,0,0,0],
										[0,0,0,1,0,0,0,0],
										[0,0,0,0,1,0,0,0],
										[0,0,0,0,0,1,0,0],
										[0,0,0,0,0,0,1,0],
										[0,0,0,0,0,0,0,1] 	])
		permute67 = 		np.matrix([	[1,0,0,0,0,0,0,0],
										[0,1,0,0,0,0,0,0],
										[0,0,1,0,0,0,0,0],
										[0,0,0,1,0,0,0,0],
										[0,0,0,0,1,0,0,0],
										[0,0,0,0,0,0,1,0],
										[0,0,0,0,0,1,0,0],
										[0,0,0,0,0,0,0,1] 	])
		permute23and67 = 	np.matrix([	[1,0,0,0,0,0,0,0],
										[0,0,1,0,0,0,0,0],
										[0,1,0,0,0,0,0,0],
										[0,0,0,1,0,0,0,0],
										[0,0,0,0,1,0,0,0],
										[0,0,0,0,0,0,1,0],
										[0,0,0,0,0,1,0,0],
										[0,0,0,0,0,0,0,1] 	])

		# Define the B matrix
		B = np.matrix([[-1,1,0,0,0,0,0,0],
					   [0,-1,1,0,0,0,0,0],
					   [0,0,-1,1,0,0,0,0],
					   [0,0, 0,0,0,0,0,0],
					   [0,0,0,0,-1,1,0,0],
					   [0,0,0,0,0,-1,1,0],
					   [0,0,0,0,0,0,-1,1],
					   [0,0,0,0,0,0,0,0]])



		# Check whether all leptons are same flavour
		all_leptons_equal = bool(abs(p2[0,5])==abs(p3[0,5])==abs(p6[0,5])==abs(p7[0,5]))
		all_leptons_equal_list.append(all_leptons_equal)

		# There are 16 possible combinations for all-leptons-equal: There are two possibilities for the quarks, and for each quark ordering there are two possible pairings of leptons (the options are limited by lepton sign). Within each lepton pairing we must evaluate both near-far combinations. Hence totally 2^4 = 16 possibilities. If the leptons pairs are unequal there are only 8 combinations.
		# A1 - the original ordering
		A11 = 2/float(Mnorm) * np.matrix([[ p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] , 0 , 0 , 0 , 0 ],
										 [ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
										 [ p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] , 0 , 0 , 0 , 0 ],
										 [ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
										 [ 0	, 0 , 0 , 0 , p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] ],
										 [ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
										 [ 0	, 0 , 0 , 0 , p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] ],
										 [ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]])
		A12 = permute23*A11
		A13 = permute67*A11
		A14 = permute23and67*A11

		D11_list.append(A11.I*B)
		D12_list.append(A12.I*B)
		D13_list.append(A13.I*B)
		D14_list.append(A14.I*B)



		# A2 - flip the quarks
		A21 = 2/float(Mnorm) * np.matrix([[ p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] , 0 , 0 , 0 , 0 ],
										 [ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
										 [ p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] , 0 , 0 , 0 , 0 ],
										 [ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
										 [ 0	, 0 , 0 , 0 , p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] ],
										 [ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
										 [ 0	, 0 , 0 , 0 , p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] ],
										 [ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]])


		A22 = permute23*A21
		A23 = permute67*A21
		A24 = permute23and67*A21

		D21_list.append(A21.I*B)
		D22_list.append(A22.I*B)
		D23_list.append(A23.I*B)
		D24_list.append(A24.I*B)

		#C vector
		C1 = 1/float(Mnorm)**2 * np.matrix([ 2*minkowskidot(p1,p2) + 2*minkowskidot(p1,p3) + m1squared,
						2*minkowskidot(p2,p3) + m2squared,
						m3squared,
						pxmiss**2,
						2*minkowskidot(p5,p6) + 2*minkowskidot(p5,p7) + m5squared,
						2*minkowskidot(p6,p7) + m6squared,
						m7squared,
						pymiss**2])
		C2 = 1/float(Mnorm)**2 * np.matrix([ 2*minkowskidot(p5,p2) + 2*minkowskidot(p5,p3) + m5squared,
						2*minkowskidot(p2,p3) + m2squared,
						m3squared,
						pxmiss**2,
						2*minkowskidot(p1,p6) + 2*minkowskidot(p1,p7) + m1squared,
						2*minkowskidot(p6,p7) + m6squared,
						m7squared,
						pymiss**2])

		C1 = np.transpose(C1) #make column vector
		C2 = np.transpose(C2)

		E11_list.append(A11.I*C1)
		E12_list.append(A12.I*C1)
		E13_list.append(A13.I*C1)
		E14_list.append(A14.I*C1)

		E21_list.append(A21.I*C2)
		E22_list.append(A22.I*C2)
		E23_list.append(A23.I*C2)
		E24_list.append(A24.I*C2)


		if (all_leptons_equal):
			# For lepton flipping between sides we must check signs. We choose to always flip 3 down and check which of 6,7 is the same sign as 3. Do the check by assuming it's 7, switching p6 with p7 if not.
			p67_was_switched = False
			if (np.sign(p7[0,5]) != np.sign(p3[0,5])):
				p6, p7 = p7, p6
				p67_was_switched = True
			#end if

			# A3 - original quark ordering, lepton flip
			A31 = 2/float(Mnorm) * np.matrix([[ p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] , 0 , 0 , 0 , 0 ],
											 [ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
											 [ p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] , 0 , 0 , 0 , 0 ],
											 [ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
											 [ 0	, 0 , 0 , 0 , p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] ],
											 [ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
											 [ 0	, 0 , 0 , 0 , p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] ],
											 [ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]])

			A32 = permute23*A31
			A33 = permute67*A31
			A34 = permute23and67*A31

			# A4 - flip the quarks AND the leptons
			A41 = 2/float(Mnorm) * np.matrix([[ p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] , 0 , 0 , 0 , 0 ],
											 [ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
											 [ p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] , 0 , 0 , 0 , 0 ],
											 [ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
											 [ 0	, 0 , 0 , 0 , p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] ],
											 [ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
											 [ 0	, 0 , 0 , 0 , p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] ],
											 [ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]])

			A42 = permute23*A41
			A43 = permute67*A41
			A44 = permute23and67*A41

			D31_list.append(A31.I*B)
			D32_list.append(A32.I*B)
			D33_list.append(A33.I*B)
			D34_list.append(A34.I*B)

			D41_list.append(A41.I*B)
			D42_list.append(A42.I*B)
			D43_list.append(A43.I*B)
			D44_list.append(A44.I*B)



			#C vector
			C3 = 1/float(Mnorm)**2 * np.matrix([ 2*minkowskidot(p1,p2) + 2*minkowskidot(p1,p7) + m1squared,
							2*minkowskidot(p2,p7) + m2squared,
							m7squared,
							pxmiss**2,
							2*minkowskidot(p5,p6) + 2*minkowskidot(p5,p3) + m5squared,
							2*minkowskidot(p6,p3) + m6squared,
							m3squared,
							pymiss**2])
			C4 = 1/float(Mnorm)**2 * np.matrix([ 2*minkowskidot(p5,p2) + 2*minkowskidot(p5,p7) + m5squared,
							2*minkowskidot(p2,p7) + m2squared,
							m7squared,
							pxmiss**2,
							2*minkowskidot(p1,p6) + 2*minkowskidot(p1,p3) + m1squared,
							2*minkowskidot(p6,p3) + m6squared,
							m3squared,
							pymiss**2])
			C3 = np.transpose(C3)
			C4 = np.transpose(C4)	

			E31_list.append(A31.I*C3)
			E32_list.append(A32.I*C3)
			E33_list.append(A33.I*C3)
			E34_list.append(A34.I*C3)

			E41_list.append(A41.I*C4)
			E42_list.append(A42.I*C4)
			E43_list.append(A43.I*C4)
			E44_list.append(A44.I*C4)
		else:
			# If not all_leptons_equal, we store 0's to keep list indices right


			D31_list.append(0)
			D32_list.append(0)
			D33_list.append(0)
			D34_list.append(0)

			D41_list.append(0)
			D42_list.append(0)
			D43_list.append(0)
			D44_list.append(0)

			E31_list.append(0)
			E32_list.append(0)
			E33_list.append(0)
			E34_list.append(0)

			E41_list.append(0)
			E42_list.append(0)
			E43_list.append(0)
			E44_list.append(0)

		# END IF all_leptons_equal

	# END LOOP OVER EVENTS/i


	# Collect all matrices and vectors in common lists
	D_lists = [D11_list, D12_list, D13_list, D14_list, D21_list, D22_list, D23_list, D24_list, D31_list, D32_list, D33_list, D34_list, D41_list, D42_list, D43_list, D44_list]
	E_lists = [E11_list, E12_list, E13_list, E14_list, E21_list, E22_list, E23_list, E24_list, E31_list, E32_list, E33_list, E34_list, E41_list, E42_list, E43_list, E44_list]

	# Loop over bins to analyze each bin. Store best-fit values in matrix best_fit.
	best_fit = np.zeros((Nbins,6))
	for j in range(Nbins):
		print "bin = ", j+1
		m = sciopt.minimize(xisquared_identical_chains_with_combinatorics, Minitial, 
						  args=(D_lists, E_lists ,Nevents, j, all_leptons_equal_list, Mnorm), method='Nelder-Mead', 
						  # bounds=((0, None), (0, None), (0, None), (0, None))
						  # tol=1,
						  options={'maxiter': 2000,'maxfev': 2000,'xtol': 1e-1,
						  		   'ftol':1e-1, 'disp': True}
						  )
		best_fit[j,:] = m.x[0], m.x[1], m.x[2], m.x[3], m.nfev, m.fun

	return best_fit, D_lists, E_lists, all_leptons_equal_list
	# END best_fit definition


# ===================
# ==== MAIN PART ====
# ===================



# True masses
MSuL = 565.312 # Mass of ~uL, ~cL
MSdL = 570.734 # Mass of ~dl, ~sL
Msquark = MZ =(MSuL+MSdL)/2.0 # mean squark mass, fit this
Mchi2 = MY = 180.337 # Mass of ~chi02
Mslepton = MX = 144.06 # Mass of ~eR, ~muR
# include left-handed sleptons? Must be done in the branchings before Herwig simulation in case
Mchi1 = MN = 9.70071979E+01 # Mass of ~chi01 (dark matter!)
MZprim = MZ
MYprim = MY
MXprim = MX
MNprim = MN
# Define normalizing mass (characteristic mass scale of the problem)
Mnorm = 100.0
true_values = np.array([MZ, MY, MX, MN])

# Choose bin size, number of bins and start value
Nevents = 25
Nbins = 100


## Choose dataset ##
# file = open("on-shell_decay_squarks_at_rest_10000_events.txt",'r')
# file = open("on-shell_decay_squarks_with_pz_14TeV-CoM_2500_events_possibly_corrected.dat",'r')
# file = open("Pythia_cascade_events_no_ISR_or_FSR_20150120.dat", 'r')
# file = open("Pythia_cascade_events_20150120.dat",'r')
file = open("../events/herwigpp_9385_events_20150225.dat",'r')


# === RUN BEST-FIT ===
mass_offset = 1.0
Minitial = true_values*mass_offset
best_fit, Ainv_lists, C_lists, all_leptons_equal_list = best_fit(Nbins, Nevents, true_values, Minitial, Mnorm, file)


for i in range(Nbins):
	# print "%3d % .6e   % .6e   % .6e   % .6e   %3d   % .6e" %(i+1, best_fit[i,0], best_fit[i,1], best_fit[i,2], best_fit[i,3], best_fit[i,4], best_fit[i,5])
	print "%3d %2.1f   %2.1f   %2.1f   %2.1f   %3d   % .6e" %(i+1, best_fit[i,0], best_fit[i,1], best_fit[i,2], best_fit[i,3], best_fit[i,4], best_fit[i,5])

# print "Exit"
# sys.exit(0)

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
cut = 100**10 # xi^2 cut value in units of (100 GeV)^4
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
ylim=[400,650] # MODIFIED 20150203
xlim=[0,300]
#print xlim, ylim

plt.figure(1)
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

# plot_counter += 1

plt.show()


# # Plot xi^2 as function of some masses to see how bumpy
# from mpl_toolkits.mplot3d import Axes3D

# N = Nevents
# Nx = Ny = 300
# j = 0

# msquark_linspace = np.linspace(Msquark*0.0, Msquark*3, Nx)
# mchi2_linspace   = np.linspace(Mchi2*0.0, Mchi2*3, Nx)
# mslepton_linspace = np.linspace(Mslepton*0.0, Mslepton*3, Nx)
# mchi1_linspace = np.linspace(Mchi1*0.0, Mchi1*3, Nx)
# msquark_mesh1, mchi2_mesh = np.meshgrid(msquark_linspace, mchi2_linspace)
# msquark_mesh2, mslepton_mesh = np.meshgrid(msquark_linspace, mslepton_linspace)
# msquark_mesh3, mchi1_mesh = np.meshgrid(msquark_linspace, mchi1_linspace)

# xi2_plot_squarkchi2 = np.log(xisquared_identical_chains_with_combinatorics( msquark_mesh1, mchi2_mesh, Mslepton, Mchi1, Ainv_lists, C_lists, N, j, all_leptons_equal_list, Nx, Ny))
# xi2_plot_squarkslepton = np.log(xisquared_identical_chains_with_combinatorics(msquark_mesh2, Mchi2, mslepton_mesh, Mchi1, Ainv_lists, C_lists, N, j, all_leptons_equal_list, Nx, Ny))
# xi2_plot_squarkchi1 = np.log(xisquared_identical_chains_with_combinatorics( msquark_mesh3, Mchi2, Mslepton, mchi1_mesh, Ainv_lists, C_lists, N, j, all_leptons_equal_list, Nx, Ny))

# # Plot 1: squark-chi2
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d', )
# ax.set_zscale(u'linear')
# ax.plot_wireframe(msquark_mesh1, mchi2_mesh, xi2_plot_squarkchi2, rstride=10, cstride=10, color='k')
# # plt.title('test')

# ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
# ax.set_ylabel(r'$m_{\tilde \chi_2^0}$', {'fontsize':20})
# ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

# # plt.savefig('3D_plot_xisquared_25_herwig_events_with_combinatorics_squark-chi2.pdf', format='pdf')

# plt.show()


# # Plot 2: squark-slepton
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d', )
# ax.set_zscale(u'linear')
# ax.plot_wireframe(msquark_mesh2, mslepton_mesh, xi2_plot_squarkslepton, rstride=10, cstride=10, color='k')
# # plt.title('test')

# ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
# ax.set_ylabel(r'$m_{\tilde l}$', {'fontsize':20})
# ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

# # plt.savefig('3D_plot_xisquared_25_herwig_events_with_combinatorics_squark-slepton.pdf', format='pdf')

# plt.show()


# # Plot 3: squark-chi1
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d', )
# ax.set_zscale(u'linear')
# ax.plot_wireframe(msquark_mesh3, mchi1_mesh, xi2_plot_squarkchi1, rstride=10, cstride=10, color='k')
# # plt.title('test')

# ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
# ax.set_ylabel(r'$m_{\tilde \chi_1^0}$', {'fontsize':20})
# ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

# # plt.savefig('3D_plot_xisquared_25_herwig_events_with_combinatorics_squark-chi1.pdf', format='pdf')

# plt.show()




























# === Hopefully just old crap below here ===
	# 	#A matrix
	# 	A = 2*np.matrix([	[ p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] , 0 , 0 , 0 , 0 ],
	# 						[ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
	# 						[ p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] , 0 , 0 , 0 , 0 ],
	# 						[ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
	# 						[ 0	, 0 , 0 , 0 , p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] ],
	# 						[ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
	# 						[ 0	, 0 , 0 , 0 , p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] ],
	# 						[ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]]	)
	# 	A = A/Mnorm # normalize A
	# 	# print np.linalg.det(A)
	# 	#A inverse
	# 	# print A
	# 	Ainv = A.I
	# 	Ainv_list.append(Ainv)



	# 	#C vector
	# 	C = np.transpose(np.matrix([ 	2*minkowskidot(p1,p2) + 2*minkowskidot(p1,p3) + m1**2,
	# 									2*minkowskidot(p2,p3) + m2**2,
	# 									m3**2,
	# 									pxmiss**2,
	# 									2*minkowskidot(p5,p6) + 2*minkowskidot(p5,p7) + m5**2,
	# 									2*minkowskidot(p6,p7) + m6**2,
	# 									m7**2,
	# 									pymiss**2]))
	# 	C = C/Mnorm**2 # normalize C

	# 	# Composite matrix & vector D and E
	# 	# D = np.dot(Ainv,B)
	# 	# E = np.dot(Ainv,C.T)

	# 	# store D and E
	# 	# Dlist.append(D)
	# 	# Elist.append(E)
	# 	# store A and C instead, because of combinatorical permutations
	# 	Ainv_list.append(Ainv)
	# 	Clist.append(C)

	# 	# ===========
	# 	# From here on we can forget about the event momenta, everything
	# 	# is stored in Dn and En for each event. Time to guess the masses.
	# 	# ===========

	# #### end loop over i

	# # Loop over bins
	# for i in range(Nbins):

	# 	# Calculate xisquared directly in loop, not in function

	# 	# Duplicate masses for primed chain
	# 	MZp, MYp, MXp, MNp = MZ, MY, MX, MN = Masses
	# 	# Set up Webber's M vector
	# 	M = np.matrix([ MZ**2 , MY**2 , MX**2 , MN**2 , MZp**2 , MYp**2 , MXp**2 , MNp**2 ])
	# 	M = M/Mnorm**2 #normalise M

	# 	# Calculate the "chi-squared" error of the hypothesis
	# 	xisquared = 0
	# 	for n in range(i*Nevents, (i+1)*Nevents):
	# 		xisquared_current_tmp = []


	# 		for m in range(combinations):
	# 			Ainv_current_permuted = np.dot(permute[m],Ainv_list[n])
	# 			D = np.dot(Ainv_current_permuted, B)
	# 			E = np.dot(Ainv_current_permuted, Clist[n].T)

	# 			Pn = np.dot(D, M) + E
	# 			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
	# 			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2
	# 			xisquared +=  (p4nsquared - M[3])**2 + (p8nsquared - M[3])**2 # p4/p8 is normalized by Mnorm.

	# 	xisquared = xisquared/(float(Nevents))