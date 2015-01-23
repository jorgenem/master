#import stuff
import numpy as np
import matplotlib.pyplot as plt
from math import pi
import scipy.optimize as sciopt
np.random.seed(2) # set seed for reproducibility





# ====== define useful functions ========



def proj(v,u):
	# projects v onto u
	if np.linalg.norm(u) > 0:
		return np.dot(u,np.transpose(v))/float(np.dot(u,u.T)) * u
	else: 
		return u
def minkowskidot(a,b):
	# Inner product in Minkowski space
	return float(a[0,0]*b[0,0]-a[0,1]*b[0,1]-a[0,2]*b[0,2]-a[0,3]*b[0,3])
def minkowskinorm(a):
	# Inner product of a with a in Minkowski space
	return minkowskidot(a,a)
def decayfun(m1,P1,m2,m3):
	# Calculating four-momenta of particle 2&3 going back-to-back from
	# decay of particle 1 in the frame where particle 1 has 4-mom P1
	#
	#
	# particle 1 = decaying particle
	# particle 2 & particle 3 = decay products
	# primed system is rest frame of particle 1, unprimed is lab frame
	# rotated system is at rest in lab system,
	# but rotated so particle one goes in +x direction
	p1 = P1[0,1:4]
	p1abs = np.sqrt( float( np.dot( p1 , np.transpose(p1) ) ) ) # 3-momentum 
																# of particle 1 in 
												      			# lab frame

	# == Kinematical decay in RF of particle 1 ==
	p2absprime = 1.0/(2*m1) * np.sqrt( (m1**2-m2**2-m3**2)**2- 4*m2**2*m3**2 ) # abs-val
	# of 3-momentum of particle 2/3 in RF of particle 1

	U, V = np.random.uniform(0,1,2) # random 
	phi = 2*pi*U 					# point picking 
	theta = np.arccos(2*V-1) 		# on a sphere

	# Calculate cartesian 3- and 4-momentum of particle 2&3
	p2prime = np.matrix([ p2absprime*np.sin(theta)*np.cos(phi) , 
						  p2absprime*np.sin(theta)*np.sin(phi) , 
						  p2absprime*np.cos(theta) ])
	p3prime = -p2prime
	E2prime = np.sqrt( p2absprime**2 + m2**2 )
	E3prime = np.sqrt( p2absprime**2 + m3**2 )
	P2prime = np.matrix([ E2prime , p2prime[0,0] , p2prime[0,1] , p2prime[0,2] ])
	P3prime = np.matrix([ E3prime , p3prime[0,0] , p3prime[0,1] , p3prime[0,2] ])

	# == Back-transform to lab frame ==

	# First check whether it is necessary to boost

	if p1abs > 1e-10:

		# Lorentz boost along x-direction to get to rotated lab frame
		# (lab frame moves in negative x direction)
	 	vlab = -p1abs/np.sqrt(p1abs**2 + m1**2) # velocity of particle 1 in lab frame
		gamma = 1/np.sqrt(1-vlab**2)

		P2rot = np.matrix([ gamma*(P2prime[0,0] - vlab*P2prime[0,1]) , 
				      gamma*(P2prime[0,1] - vlab*P2prime[0,0]) ,
				      P2prime[0,2] , P2prime[0,3] ])
		P3rot = np.matrix([ gamma*(P3prime[0,0] - vlab*P3prime[0,1]) , 
				      gamma*(P3prime[0,1] - vlab*P3prime[0,0]) ,
				      P3prime[0,2] , P3prime[0,3] ])

		# == Rotate back to lab frame ==

		# Calculate the unit vectors of the rotated system axes in terms of lab axes

		# The definition is that x axis is along p1.
		# For the other axes we must make a choice - y&z directions are undetermined,
		# only the yz plane is determined from x choice. But since we have drawn 
		# random angles and the yz plane is not boosted, the choice does not matter
		# as long as we are consistent from event to event.
		# So we pick two vectors orthogonal to p1 and do Gram-Schmidt orthogonalization:
		v1 = p1
		v2 = np.matrix([ p1[0,1] , -p1[0,0] , 0 ])
		v3 = np.matrix([ p1[0,2] , 0 , -p1[0,0] ])

		u1 = v1
		u2 = v2 - proj(v2,u1)
		u3 = v3 - proj(v3,u1) - proj(v3,u2)

		xrot = u1/np.linalg.norm(u1) if np.linalg.norm(u1) > 0 else np.matrix([0,0,1])
		yrot = u2/np.linalg.norm(u2) if np.linalg.norm(u2) > 0 else np.matrix([0,1,0])
		zrot = u3/np.linalg.norm(u3) if np.linalg.norm(u3) > 0 else np.matrix([1,0,0])

		# Form a matrix T which takes a vector in the lab basis to a vector 
		# in the rotated basis by
		T = np.concatenate( (xrot , yrot , zrot) , axis=0 )
		# What we need is to rotate from rotated basis to lab basis, so we need the inverse
		# - which is the transpose, since rotation matrices are orthogonal. 
		# Also, to ease calculation, we let T be the 3x3 submatrix of T4, setting the [0,0]
		#component of T4 to 1 to leave time component invariant under this spatial rotation
		T4 = np.matrix([[1,     0,     0,    0],
						[0,T[0,0],T[0,1],T[0,2]],
						[0,T[1,0],T[1,1],T[1,2]],
						[0,T[2,0],T[2,1],T[2,2]] ])

		P2 = T4.T*P2rot.T
		P3 = T4.T*P3rot.T
		P2 = P2.T
		P3 = P3.T

	# If it was unneccessary, i.e. decay happened in lab frame, then
	else:
		P2 = P2prime
		P3 = P3prime

	# Finished!

	return P2, P3
def smear(p,resolution):
	# Smears 4-momentum according to AcerDET manual
	r = np.random.randn()
	p_smeared = p * ( 1 + r * resolution / np.sqrt(p[0,0]) )
	return p_smeared








# ==== MAIN PART ====






# Set known parameters
# SM particle masses
# u-quark and electron mass set to zero
mquark = m1 = m5 = 0;
mlepton1 = m2 = m3 = 0;
mlepton2 = m6 = m7 = 0;


# Now to make a mass hypothesis (guess the correct one)
MSuL = 565.312 # Mass of ~uL, ~cL
MSdL = 570.734 # Mass of ~dl, ~sL
Msquark = MZ =(MSuL+MSdL)/2.0 # mean squark mass, fit this
MN2 = Mchi2 = MY = 180.337 # Mass of ~chi02
Mslepton = MseR = MX = 144.06 # Mass of ~eR, ~muR
# include left-handed sleptons? Must be done in the branchings before Herwig simulation in case
MN1 = Mchi1 = MN = 9.70071979E+01 # Mass of ~chi01 (dark matter!)
MZprim = MZ
MYprim = MY
MXprim = MX
MNprim = MN

true_values = np.array([MZ,MY,MX,MN])



N = 100
Dlist = []
Elist = []
Adetlist = np.zeros(0)
A_nosmeardetlist = np.zeros(0)

# Define normalizing mass (characteristic mass scale of the problem)
Mnorm = 100
# print "Mnorm = ", Mnorm

# Save invariant masses for making triangle
invariant_mass_between_c1_leptons = [] 

# N - How much loop?
for i in range(N):
	# Loop over events to get 4-vectors for each particle for each event. 
	# Particles are numbered according to Webber (arXiv:0907.5307v2) fig. 1
	# (the lepton/antilepton ordering is arbitrary in each chain, the lepton has been 
	# chosen as 2/6 and the antilepton as 3/7)
	# string = root[i+2].text
	# lines = string.splitlines()

	#1st chain, p1-4

	# Read squark 1 from file
	# psquark1 = str(lines[4]).split()
	# print "PDG number of particle 1: ",psquark1[0] # just to check
	#p1 = [float(p1[9]), float(p1[6]), float(p1[7]), float(p1[8])]
	# psquark1 = np.matrix([ float(psquark1[9]), float(psquark1[6]), float(psquark1[7]), float(psquark1[8])])
	psquark1 = np.matrix([0,0,0,Msquark])

	# #DEBUG
	# psquark1 = np.matrix([ Msquark, 0, 0, 0]) # overwrite CompHEP data to start squarks at rest
	# #/DEBUG

	# Decay squark to quark and neutralino2
	p1, pN21 = decayfun(Msquark,psquark1,mquark,MN2)
	# Decay neutralino2 to lepton1 and slepton
	p2, pseR1 = decayfun(MN2,pN21,mlepton1,MseR)
	# Decay slepton to (anti)lepton1 and neutralino1
	p3, p4 = decayfun(MseR,pseR1,mlepton1,MN1)




	#2nd chain, p5-8
	# psquark2 = str(lines[5]).split()
	# # print "PDG number of particle 5: ",psquark2[0] # just to check
	# psquark2 = np.matrix([ float(psquark2[9]), float(psquark2[6]), float(psquark2[7]), float(psquark2[8])])
	psquark2 = np.matrix([0,0,0,Msquark])

	# #DEBUG
	# psquark2 = np.matrix([ Msquark, 0, 0, 0]) # overwrite CompHEP data to start squarks at rest
	# #/DEBUG

	# See whether CompHEP produces squarks off-shell
	# print minkowskinorm(psquark1) - Msquark**2
	# print minkowskinorm(psquark2) - Msquark**2

	# Decay (anti)squark to (anti)quark and neutralino2
	p5, pN22 = decayfun(Msquark,psquark2,mquark,MN2)
	# Decay neutralino2 to lepton2 and slepton
	p6, pseR2 = decayfun(MN2,pN22,mlepton2,MseR)
	# Decay slepton to (anti)lepton2 and neutralino1
	p7, p8 = decayfun(MseR,pseR2,mlepton2,MN1)




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




	# Smear, r percent resolution
	# r = resolution # percent/100 momentum smearing
	# p1 = smear(p1,r)
	# p2 = smear(p2,r)
	# p3 = smear(p3,r)

	# p5 = smear(p5,r)
	# p6 = smear(p6,r)
	# p7 = smear(p7,r)



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
	# invmass_c1leptons = np.sqrt(minkowskinorm(p2+p3))   # calculate invariant 
	# invariant_mass_between_c1_leptons.append(invmass_c1leptons) # mass between leptons in chain 1




	# ==== Define Webber's stuff ====

	# need the pxmiss and pymiss, taken from the actual neutralino transverse momenta 
	# (this is cheating, of course)
	# pxmiss = p4[0,1]+p8[0,1]
	# pymiss = p4[0,2]+p8[0,2]

	# Calculate missing transverse from (smeared) visible particles
	pxmiss = - p1[0,1] - p2[0,1] - p3[0,1] - p5[0,1] - p6[0,1] - p7[0,1]
	pymiss = - p1[0,2] - p2[0,2] - p3[0,2] - p5[0,2] - p6[0,2] - p7[0,2]

	m1square = minkowskinorm(p1)
	m2square = minkowskinorm(p2)
	m3square = minkowskinorm(p3)
	m5square = minkowskinorm(p5)
	m6square = minkowskinorm(p6)
	m7square = minkowskinorm(p7)

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




# ============ Minimization to best fit =================
# import minuit


# xi-squared function to minimize with identical chains
def xisquared_identical_chains(MZ, MY, MX, MN, Nevents, i): #, MZp, MYp, MXp, MNp):
	Nevents = int(Nevents)
	i = int(i)
	# Duplicate masses for primed chain
	MZp, MYp, MXp, MNp = MZ, MY, MX, MN# = Masses
	# Set up Webber's M vector
	M = np.matrix([ MZ**2 , MY**2 , MX**2 , MN**2 , MZp**2 , MYp**2 , MXp**2 , MNp**2 ])
	M = M/Mnorm**2 #normalise M
	# Calculate the "chi-squared" error of the hypothesis
	# P = [] # store Pn
	xisquared = 0
	# offshell = [] # list to store p4nsquared - MN**2
	for n in range(i*Nevents, (i+1)*Nevents):
		Pn = np.dot(Dlist[n],M.T) + Elist[n]
		# P.append(Pn) #store in case needed
	
	
		p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
		p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2
	
		xisquared +=  (p4nsquared - M[0,3])**2 + (p8nsquared - M[0,7])**2 # p4/p8 is normalized by MN.

		# offshell.append(abs(p4nsquared-MN**2))
		# offshell.append(abs(p8nsquared-MNprim**2))
	xisquared = xisquared/float(Nevents)#/100**4
	# print xisquared
	# xisquaredlist.append(xisquared)
	return xisquared



# Plot xi^2 as function of some masses to see how bumpy
from mpl_toolkits.mplot3d import Axes3D

minm = 0.0
maxm = 3
Nlinspace = 300
Nevents=25
bin_number=1

msquark_linspace = np.linspace(Msquark*minm, Msquark*maxm, Nlinspace)
mchi2_linspace   = np.linspace(Mchi2*minm, Mchi2*maxm, Nlinspace)
mslepton_linspace = np.linspace(Mslepton*minm, Mslepton*maxm, Nlinspace)
mchi1_linspace = np.linspace(Mchi1*minm, Mchi1*maxm, Nlinspace)
msquark_mesh1, mchi2_mesh = np.meshgrid(msquark_linspace, mchi2_linspace)
msquark_mesh2, mslepton_mesh = np.meshgrid(msquark_linspace, mslepton_linspace)
msquark_mesh3, mchi1_mesh = np.meshgrid(msquark_linspace, mchi1_linspace)
# mslepton_mesh2, mchi1_mesh2 = np.meshgrid(mslepton_linspace, mchi1_linspace)

xi2_plot_squarkchi2 = np.log10(xisquared_identical_chains(msquark_mesh1, mchi2_mesh, Mslepton, Mchi1, Nevents,bin_number-1))
xi2_plot_squarkslepton = np.log10(xisquared_identical_chains(msquark_mesh2, Mchi2, mslepton_mesh, Mchi1, Nevents, bin_number-1))
xi2_plot_squarkchi1 = np.log10(xisquared_identical_chains(msquark_mesh3, Mchi2, Mslepton, mchi1_mesh, Nevents, bin_number-1))
# xi2_plot_sleptonchi1 = np.log10(xisquared_identical_chains(480, 150, mslepton_mesh2, mchi1_mesh2, Nevents, bin_number-1))

# Plot 1: squark-chi2
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', )
ax.set_zscale(u'linear')
ax.plot_wireframe(msquark_mesh1, mchi2_mesh, xi2_plot_squarkchi2, rstride=10, cstride=10, color='k')
# plt.title('test')

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde \chi_2^0}$', {'fontsize':20})
ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

plt.show()


# Plot 2: squark-slepton
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', )
ax.set_zscale(u'linear')
ax.plot_wireframe(msquark_mesh2, mslepton_mesh, xi2_plot_squarkslepton, rstride=10, cstride=10, color='k')
# plt.title('test')

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde l}$', {'fontsize':20})
ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

plt.show()


# Plot 3: squark-chi1
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d', )
ax.set_zscale(u'linear')
ax.plot_wireframe(msquark_mesh3, mchi1_mesh, xi2_plot_squarkchi1, rstride=10, cstride=10, color='k')
# plt.title('test')

ax.set_xlabel(r'$m_{\tilde q}$', {'fontsize':20})
ax.set_ylabel(r'$m_{\tilde \chi_1^0}$', {'fontsize':20})
ax.set_zlabel(r'$\log (\xi^2)$', {'fontsize':18})

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