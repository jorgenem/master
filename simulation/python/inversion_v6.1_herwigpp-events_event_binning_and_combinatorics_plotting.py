#import stuff
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from iminuit import Minuit
import scipy.optimize as sciopt
from mpl_toolkits.mplot3d import Axes3D

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


def xisquared_identical_chains_with_combinatorics(MZ, MY, MX, MN, Ainv_lists, C_lists, Nevents, j, all_leptons_equal_list, Nx, Ny):
	Nevents = int(Nevents)
	j = int(j)

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

	
	# Import masses
	MZp, MYp, MXp, MNp = MZ, MY, MX, MN #= Masses
	# Set up the M vector of current mass guess
	M = 1/float(Mnorm)**2 * np.matrix([ MZ**2 , MY**2 , MX**2 , MN**2 , MZp**2 , MYp**2 , MXp**2 , MNp**2 ])
	M = M.T

	best_fit_combination = np.zeros((Nevents,Nx,Ny)) # array of matrices of best-fit combination indices for each event. 0 is index of true combination.
	xisquared = 0 # begin accumulation of the xisquared sum
	for n in range(j*Nevents, (j+1)*Nevents):

		# Must check all-leptons-equal on a per-event basis
		if all_leptons_equal_list[n]:
			xisquared_current_tmp = np.ones((16,Nx,Ny))*1e6
			# 16 combinations

			# Case 1 - A1/C1
			Ainv = Ainv_lists[0][n]
			C = C_lists[0][n]
			# Subcase 1 - no permutations
			Pn = np.dot( np.dot(Ainv, B), M ) + np.dot(Ainv, C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[0,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 2 - permute 23
			Pn = np.dot( np.dot(permute23 * Ainv, B), M ) + np.dot(permute23 * Ainv, permute23 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[1,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 3 - permute 67
			Pn = np.dot( np.dot(permute67 * Ainv, B), M ) + np.dot(permute67 * Ainv, permute67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[2,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 4 - permute 23 and 67
			Pn = np.dot( np.dot(permute23and67 * Ainv, B), M ) + np.dot(permute23and67 * Ainv, permute23and67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[3,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 


			# Case 2 - A2/C2
			Ainv = Ainv_lists[1][n]
			C = C_lists[1][n]
			# Subcase 1 - no permutations
			Pn = np.dot( np.dot(Ainv, B), M ) + np.dot(Ainv, C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[4,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 2 - permute 23
			Pn = np.dot( np.dot(permute23 * Ainv, B), M ) + np.dot(permute23 * Ainv, permute23 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[5,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 3 - permute 67
			Pn = np.dot( np.dot(permute67 * Ainv, B), M ) + np.dot(permute67 * Ainv, permute67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[6,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 4 - permute 23 and 67
			Pn = np.dot( np.dot(permute23and67 * Ainv, B), M ) + np.dot(permute23and67 * Ainv, permute23and67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[7,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 


			# Case 3 - A3/C3
			Ainv = Ainv_lists[2][n]
			C = C_lists[2][n]
			# Subcase 1 - no permutations
			Pn = np.dot( np.dot(Ainv, B), M ) + np.dot(Ainv, C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[8,:,:] =  (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 2 - permute 23
			Pn = np.dot( np.dot(permute23 * Ainv, B), M ) + np.dot(permute23 * Ainv, permute23 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[9,:,:] =  (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 3 - permute 67
			Pn = np.dot( np.dot(permute67 * Ainv, B), M ) + np.dot(permute67 * Ainv, permute67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[10,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 4 - permute 23 and 67
			Pn = np.dot( np.dot(permute23and67 * Ainv, B), M ) + np.dot(permute23and67 * Ainv, permute23and67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[11,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 


			# Case 4 - A4/C4
			Ainv = Ainv_lists[3][n]
			C = C_lists[3][n]
			# Subcase 1 - no permutations
			Pn = np.dot( np.dot(Ainv, B), M ) + np.dot(Ainv, C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[12,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 2 - permute 23
			Pn = np.dot( np.dot(permute23 * Ainv, B), M ) + np.dot(permute23 * Ainv, permute23 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[13,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 3 - permute 67
			Pn = np.dot( np.dot(permute67 * Ainv, B), M ) + np.dot(permute67 * Ainv, permute67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[14,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2  

			# Subcase 4 - permute 23 and 67
			Pn = np.dot( np.dot(permute23and67 * Ainv, B), M ) + np.dot(permute23and67 * Ainv, permute23and67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[15,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 



		else:
			# 8 combinations
			xisquared_current_tmp = np.ones((8,Nx,Ny))*1e6

			# Case 1 - A1/C1
			Ainv = Ainv_lists[0][n]
			C = C_lists[0][n]
			# Subcase 1 - no permutations
			Pn = np.dot( np.dot(Ainv, B), M ) + np.dot(Ainv, C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[0,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 2 - permute 23
			Pn = np.dot( np.dot(permute23 * Ainv, B), M ) + np.dot(permute23 * Ainv, permute23 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[1,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 3 - permute 67
			Pn = np.dot( np.dot(permute67 * Ainv, B), M ) + np.dot(permute67 * Ainv, permute67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[2,:,:] =  (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 4 - permute 23 and 67
			Pn = np.dot( np.dot(permute23and67 * Ainv, B), M ) + np.dot(permute23and67 * Ainv, permute23and67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[3,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 


			# Case 2 - A2/C2
			Ainv = Ainv_lists[1][n]
			C = C_lists[1][n]
			# Subcase 1 - no permutations
			Pn = np.dot( np.dot(Ainv, B), M ) + np.dot(Ainv, C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[4,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 2 - permute 23
			Pn = np.dot( np.dot(permute23 * Ainv, B), M ) + np.dot(permute23 * Ainv, permute23 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[5,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 3 - permute 67
			Pn = np.dot( np.dot(permute67 * Ainv, B), M ) + np.dot(permute67 * Ainv, permute67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[6,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 

			# Subcase 4 - permute 23 and 67
			Pn = np.dot( np.dot(permute23and67 * Ainv, B), M ) + np.dot(permute23and67 * Ainv, permute23and67 * C)

			p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
			p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

			xisquared_current_tmp[7,:,:] = (p4nsquared - M[3,0])**2 + (p8nsquared - M[7,0])**2 


			

		# END IF all_leptons_equal


		# Check which combination gives the smallest xisquared contribution
		xisquared_current = np.min(xisquared_current_tmp,0)

		# Store index of best-fit combination among the 8/16 possibilities
		for k in range(Nx):
			for l in range(Ny):
				best_fit_combination[n-j*Nevents,k,l] = np.where(xisquared_current_tmp[:,k,l]==xisquared_current[k,l])[0]
				if len(np.where(xisquared_current_tmp[:,k,l]==xisquared_current[k,l]))>1:
					print "!"
				# print np.where(xisquared_current_tmp[:,k,l]==xisquared_current[k,l])



		# print best_fit_combination
		# best_fit_combination.append(xisquared_current_tmp.index(xisquared_current))
		# if xisquared_current_tmp.count(xisquared_current) > 1:
		# 	print "Warning: multiple best-fit combinations for event ", n

		# Add best-fit combination xisquared value to total
		xisquared += xisquared_current
	# END loop over events/n

	number_of_correct_combinations = np.zeros((Nx, Ny))

	for k in range(Nx):
		for l in range(Ny):

			number_of_correct_combinations[k,l] = Nevents - np.count_nonzero(best_fit_combination[:,k,l])

	xisquared = xisquared/float(Nevents)
	return number_of_correct_combinations
	# END xisquared definition


def best_fit(Nbins, Nevents, Mtrue, Minitial, Mnorm):
	# Make lists for storing D matrices and E vectors
	N = Nbins*Nevents
	Plist = []
	A1inv_list = []
	A2inv_list = []
	A3inv_list = []
	A4inv_list = []
	C1_list = []
	C2_list = []
	C3_list = []
	C4_list = []
	all_leptons_equal_list = []

	#import the Herwig .txt file of events
	import sys
	file = open("../herwigpp/LHC-MSSM-analysis_20141128_softsusy_with_branching_to_all_four_correct_leptons_31000_events.log",'r')
	lines = file.readlines()

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
		p1[0,0:5] /= 1000
		p2[0,0:5] /= 1000
		p3[0,0:5] /= 1000
		p4[0,0:5] /= 1000
		p5[0,0:5] /= 1000
		p6[0,0:5] /= 1000
		p7[0,0:5] /= 1000
		p8[0,0:5] /= 1000

		# Smear
		p1 = smear2(p1)
		p2 = smear2(p2)
		p3 = smear2(p3)
		p5 = smear2(p5)
		p6 = smear2(p6)
		p7 = smear2(p7)

		# Calculate invariant masses of measured particles after smearing and replace
		m1 = p1[0,4] = np.sign(minkowskinorm(p1))*np.sqrt(abs(minkowskinorm(p1)))
		m2 = p2[0,4] = np.sign(minkowskinorm(p2))*np.sqrt(abs(minkowskinorm(p2)))
		m3 = p3[0,4] = np.sign(minkowskinorm(p3))*np.sqrt(abs(minkowskinorm(p3)))
		m5 = p5[0,4] = np.sign(minkowskinorm(p5))*np.sqrt(abs(minkowskinorm(p5)))
		m6 = p6[0,4] = np.sign(minkowskinorm(p6))*np.sqrt(abs(minkowskinorm(p6)))
		m7 = p7[0,4] = np.sign(minkowskinorm(p7))*np.sqrt(abs(minkowskinorm(p7)))
		m1squared = minkowskinorm(p1)
		m2squared = minkowskinorm(p2)
		m3squared = minkowskinorm(p3)
		m5squared = minkowskinorm(p5)
		m6squared = minkowskinorm(p6)
		m7squared = minkowskinorm(p7)

		# Calculate missing transverse from (smeared) visible particles
		pxmiss = - p1[0,1] - p2[0,1] - p3[0,1] - p5[0,1] - p6[0,1] - p7[0,1]
		pymiss = - p1[0,2] - p2[0,2] - p3[0,2] - p5[0,2] - p6[0,2] - p7[0,2]

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

		# Check whether all leptons are same flavour
		all_leptons_equal = bool(abs(p2[0,5])==abs(p3[0,5])==abs(p6[0,5])==abs(p7[0,5]))
		all_leptons_equal_list.append(all_leptons_equal)

		# There are 16 possible combinations for all-leptons-equal: There are two possibilities for the quarks, and for each quark ordering there are two possible pairings of leptons (the options are limited by lepton sign). Within each lepton pairing we must evaluate both near-far combinations. Hence totally 2^4 = 16 possibilities. If the leptons pairs are unequal there are only 8 combinations.
		# A1 - the original ordering
		A1 = 2/float(Mnorm) * np.matrix([[ p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] , 0 , 0 , 0 , 0 ],
										 [ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
										 [ p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] , 0 , 0 , 0 , 0 ],
										 [ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
										 [ 0	, 0 , 0 , 0 , p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] ],
										 [ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
										 [ 0	, 0 , 0 , 0 , p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] ],
										 [ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]])
		# A2 - flip the quarks
		A2 = 2/float(Mnorm) * np.matrix([[ p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] , 0 , 0 , 0 , 0 ],
										 [ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
										 [ p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] , 0 , 0 , 0 , 0 ],
										 [ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
										 [ 0	, 0 , 0 , 0 , p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] ],
										 [ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
										 [ 0	, 0 , 0 , 0 , p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] ],
										 [ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]])


		A1inv = A1.I
		A2inv = A2.I
		A1inv_list.append(A1inv)
		A2inv_list.append(A2inv)

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
		C1_list.append(C1)
		C2_list.append(C2)


		if (all_leptons_equal):
			# For lepton flipping between sides we must check signs. We choose to always flip 3 down and check which of 6,7 is the same sign as 3. Do the check by assuming it's 7, switching p6 with p7 if not.
			p67_was_switched = False
			if (np.sign(p7[0,5]) != np.sign(p3[0,5])):
				p6, p7 = p7, p6
				p67_was_switched = True
			#end if

			# A3 - original quark ordering, lepton flip
			A3 = 2/float(Mnorm) * np.matrix([[ p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] , 0 , 0 , 0 , 0 ],
											 [ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
											 [ p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] , 0 , 0 , 0 , 0 ],
											 [ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
											 [ 0	, 0 , 0 , 0 , p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] ],
											 [ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
											 [ 0	, 0 , 0 , 0 , p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] ],
											 [ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]])
			# A4 - flip the quarks AND the leptons
			A4 = 2/float(Mnorm) * np.matrix([[ p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] , 0 , 0 , 0 , 0 ],
											 [ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
											 [ p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] , 0 , 0 , 0 , 0 ],
											 [ 0.5*pxmiss,	0   , 0		  , 0		 , 0.5*pxmiss,0 , 0 , 0 ],
											 [ 0	, 0 , 0 , 0 , p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] ],
											 [ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
											 [ 0	, 0 , 0 , 0 , p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] ],
											 [ 0 ,0.5*pymiss, 0 , 0 , 0 	  , 0.5*pymiss	, 0		  , 0 		 ]])
			A3inv = A3.I
			A4inv = A4.I
			A3inv_list.append(A3inv)
			A4inv_list.append(A4inv)

			#C vector
			C3 = 1/float(Mnorm)**2 * np.matrix([ 2*minkowskidot(p1,p2) + 2*minkowskidot(p1,p7) + m1squared,
							2*minkowskidot(p2,p7) + m2squared,
							m7squared,
							pxmiss**2,
							2*minkowskidot(p5,p6) + 2*minkowskidot(p5,p3) + m5squared,
							2*minkowskidot(p6,p3) + m6squared,
							m3**2,
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
			C3_list.append(C3)
			C4_list.append(C4)
		else:
			# If not all_leptons_equal, we store 0's to keep list indices right
			A3inv_list.append(0)
			A4inv_list.append(0)
			C3_list.append(0)
			C4_list.append(0)
		# END IF all_leptons_equal

	# END LOOP OVER EVENTS/i


	# Collect all matrices and vectors in common lists
	Ainv_lists = [A1inv_list, A2inv_list, A3inv_list, A4inv_list]
	C_lists = [C1_list, C2_list, C3_list, C4_list]

	# Loop over bins to analyze each bin. Store best-fit values in matrix best_fit.
	best_fit = np.zeros((Nbins,6))
	# for j in range(Nbins):
	# 	m = sciopt.minimize(xisquared_identical_chains, Minitial, 
	# 					  args=(Ainv_lists, C_lists ,Nevents, i, Mnorm), method='TNC', 
	# 					  bounds=((0, None), (0, None), (0, None), (0, None))
	# 					  # tol=1,
	# 					  # options={'maxiter': 100}
	# 					  )
	# 	best_fit[i,:] = m.x[0], m.x[1], m.x[2], m.x[3], m.nfev, m.fun

	return best_fit, Ainv_lists, C_lists, all_leptons_equal_list
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
Mnorm = 1000
Mtrue = np.array([MZ, MY, MX, MN])

# Choose bin size, number of bins and start value
Nbins = 2
Nevents = 10


# === RUN BEST-FIT ===
mass_offset = 2.0
Minitial = Mtrue*mass_offset
null, Ainv_lists, C_lists, all_leptons_equal_list = best_fit(Nbins, Nevents, Mtrue, Minitial, Mnorm)





# Plot xi^2 as function of some masses to see how bumpy
N = Nevents
Nx = Ny = 300
j = 0

msquark_linspace = np.linspace(Msquark*0.0, Msquark*3, Nx)
mchi2_linspace   = np.linspace(Mchi2*0.0, Mchi2*3, Nx)
mslepton_linspace = np.linspace(Mslepton*0.0, Mslepton*3, Nx)
mchi1_linspace = np.linspace(Mchi1*0.0, Mchi1*3, Nx)
msquark_mesh1, mchi2_mesh = np.meshgrid(msquark_linspace, mchi2_linspace)
msquark_mesh2, mslepton_mesh = np.meshgrid(msquark_linspace, mslepton_linspace)
msquark_mesh3, mchi1_mesh = np.meshgrid(msquark_linspace, mchi1_linspace)
mslepton_mesh2, mchi1_mesh2 = np.meshgrid(mslepton_linspace, mchi1_linspace)

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



# Below is the code for plotting the amount of true best-fit combinations. Requires that the xisquared function returns number of correct combinations instead of the xisquared value.

xi2_plot_squarkchi2 = xisquared_identical_chains_with_combinatorics( msquark_mesh1, mchi2_mesh, Mslepton, Mchi1, Ainv_lists, C_lists, N, j, all_leptons_equal_list, Nx, Ny)
xi2_plot_squarkslepton = xisquared_identical_chains_with_combinatorics(msquark_mesh2, Mchi2, mslepton_mesh, Mchi1, Ainv_lists, C_lists, N, j, all_leptons_equal_list, Nx, Ny)
xi2_plot_squarkchi1 = xisquared_identical_chains_with_combinatorics( msquark_mesh3, Mchi2, Mslepton, mchi1_mesh, Ainv_lists, C_lists, N, j, all_leptons_equal_list, Nx, Ny)
xi2_plot_sleptonchi1 = xisquared_identical_chains_with_combinatorics( Msquark, Mchi2, mslepton_mesh2, mchi1_mesh2, Ainv_lists, C_lists, N, j, all_leptons_equal_list, Nx, Ny)

xi2_plot_squarkchi2 /= Nevents
xi2_plot_squarkslepton /= Nevents
xi2_plot_squarkchi1 /= Nevents
xi2_plot_sleptonchi1 /= Nevents

# Plot 1: squark-chi2
fig, ax = plt.subplots()
p = ax.pcolor(msquark_mesh1, mchi2_mesh, xi2_plot_squarkchi2, cmap='RdBu', vmin=0, vmax=abs(xi2_plot_squarkchi2).max())
plt.xlim(msquark_linspace[0],msquark_linspace[-1])
cb = fig.colorbar(p, ax=ax)

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde \chi_2^0}$', {'fontsize':20})
# ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})
plt.title("Fraction of events with correct best-fit combination, %d events"%Nevents)

plt.savefig('fraction_of_correct_combinatorics_squark-chi2.pdf', format='pdf')

plt.show()


# Plot 2: squark-slepton
fig, ax = plt.subplots()
p = ax.pcolor(msquark_mesh1, mslepton_mesh, xi2_plot_squarkslepton, cmap='RdBu', vmin=0, vmax=abs(xi2_plot_squarkslepton).max())
cb = fig.colorbar(p, ax=ax)

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde l}$', {'fontsize':20})
# ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})
plt.title("Fraction of events with correct best-fit combination, %d events"%Nevents)

plt.savefig('fraction_of_correct_combinatorics_squark-slepton.pdf', format='pdf')

plt.show()


# # Plot 3: squark-chi1
fig, ax = plt.subplots()
p = ax.pcolor(msquark_mesh1, mchi1_mesh, xi2_plot_squarkchi1, cmap='RdBu', vmin=0, vmax=abs(xi2_plot_squarkchi1).max())
cb = fig.colorbar(p, ax=ax)

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde \chi_1^0}$', {'fontsize':20})
# ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})
plt.title("Fraction of events with correct best-fit combination, %d events"%Nevents)

plt.savefig('fraction_of_correct_combinatorics_squark-chi1.pdf', format='pdf')

plt.show()


# # Plot 4: squark-chi1
fig, ax = plt.subplots()
p = ax.pcolor(mslepton_mesh2, mchi1_mesh2, xi2_plot_sleptonchi1, cmap='RdBu', vmin=0, vmax=abs(xi2_plot_sleptonchi1).max())
cb = fig.colorbar(p, ax=ax)

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde \chi_1^0}$', {'fontsize':20})
# ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})
plt.title("Fraction of events with correct best-fit combination, %d events"%Nevents)

plt.savefig('fraction_of_correct_combinatorics_slepton-chi1.pdf', format='pdf')

plt.show()

















