#get numpy
import numpy as np

#define useful functions
def minkowskidot(a,b):
	# Inner product in Minkowski space
	return float(a[0,0]*b[0,0]-a[0,1]*b[0,1]-a[0,2]*b[0,2]-a[0,3]*b[0,3])
def minkowskinorm(a):
	# Inner product of a with a in Minkowski space
	return float(a[0,0]**2-a[0,1]**2-a[0,2]**2-a[0,3]**2)


#import the LHE file of events.
import xml.etree.ElementTree as ET
tree = ET.parse('kaskade-helveis_20140529-2.txt')
root = tree.getroot()

# Set known parameters
# SM particle masses
# u-quark and electron mass set to zero
m1=0;m2=0;m3=0;m5=0;m6=0;m7=0;


# Now to make a mass hypothesis (guess the correct one)
MZ = 545.421001 # Mass of ~uL
MY = 175.988171 # Mass of ~chi02
MX = 200.949431 # Mass of ~eL
MN = 98.361959 # Mass of ~chi01 (dark matter!)
MZprim = MZ
MYprim = MY
MXprim = MX
MNprim = MN



# Make lists for storing D matrices and E vectors
Dlist = []
Elist = []

N = 1 # How much loop?
for i in range(6,7):
	# Loop over events to get 4-vectors for each particle for each event. 
	# Particles are numbered according to Webber (arXiv:0907.5307v2) fig. 1
	# (the lepton/antilepton ordering is arbitrary in each chain, the lepton has been 
	# chosen as 2/6 and the antilepton as 3/7)
	string = root[i+2].text
	lines = string.splitlines()

	#1st chain, p1-4
	p1 = str(lines[6]).split()
	print "particle 1: ",p1[0]
	#p1 = [float(p1[9]), float(p1[6]), float(p1[7]), float(p1[8])]
	p1 = np.matrix([ float(p1[9]), float(p1[6]), float(p1[7]), float(p1[8])])

	p2 = str(lines[4]).split()
	print "particle 2:",p2[0]
	p2 = np.matrix([ float(p2[9]), float(p2[6]), float(p2[7]), float(p2[8])])

	p3 = str(lines[5]).split()
	print "particle 3: ",p3[0]
	p3 = np.matrix([ float(p3[9]), float(p3[6]), float(p3[7]), float(p3[8])])

	p4 = str(lines[7]).split()
	print "particle 4: ",p4[0]
	p4 = np.matrix([ float(p4[9]), float(p4[6]), float(p4[7]), float(p4[8])])

	#2nd chain, p5-8
	p5 = str(lines[10]).split()
	print "particle 5: ",p5[0]
	p5 = np.matrix([ float(p5[9]), float(p5[6]), float(p5[7]), float(p5[8])])

	p6 = str(lines[8]).split()
	print "particle 6: ",p6[0]
	p6 = np.matrix([ float(p6[9]), float(p6[6]), float(p6[7]), float(p6[8])])

	p7 = str(lines[9]).split()
	print "particle 7: ",p7[0]
	p7 = np.matrix([ float(p7[9]), float(p7[6]), float(p7[7]), float(p7[8])])

	p8 = str(lines[11]).split()
	print "particle 8: ",p8[0]
	p8 = np.matrix([ float(p8[9]), float(p8[6]), float(p8[7]), float(p8[8])])

	# Check invariant mass of initial colliding partons?
	#print minkowskinorm(p1+p2+p3+p4+p5+p6+p7+p8)
	# Check that the invariant mass of particles is close to shell mass
	print minkowskinorm(p1) - m1**2
	print minkowskinorm(p2) - m2**2
	print minkowskinorm(p3) - m3**2
	print minkowskinorm(p4) - MN**2
	print minkowskinorm(p5) - m5**2
	print minkowskinorm(p6) - m6**2
	print minkowskinorm(p7) - m7**2
	print minkowskinorm(p8) - MNprim**2

	# Check if invariant mass of decays match mass of decaying
	print "p3+p4 ", np.sqrt(abs(minkowskinorm(p3+p4) - MX**2))
	print "p2+p3+p4 ", np.sqrt(abs(minkowskinorm(p2+p3+p4) - MY**2))
	print "p1+p2+p3+p4 ", np.sqrt(abs(minkowskinorm(p1+p2+p3+p4) - MZ**2))


	# Define Webber's stuff

	#A matrix
	A = 2*np.matrix([[ p1[0,1] , p1[0,2] , p1[0,3] , -p1[0,0] , 0 , 0 , 0 , 0 ],
					[ p2[0,1] , p2[0,2] , p2[0,3] , -p2[0,0] , 0 , 0 , 0 , 0 ],
					[ p3[0,1] , p3[0,2] , p3[0,3] , -p3[0,0] , 0 , 0 , 0 , 0 ],
					[ 0.5	  ,	0		, 0		  , 0		 , 0.5,0 , 0 , 0 ],
					[ 0	, 0 , 0 , 0 , p5[0,1] , p5[0,2] , p5[0,3] , -p5[0,0] ],
					[ 0	, 0 , 0 , 0 , p6[0,1] , p6[0,2] , p6[0,3] , -p6[0,0] ],
					[ 0	, 0 , 0 , 0 , p7[0,1] , p7[0,2] , p7[0,3] , -p7[0,0] ],
					[ 0 ,0.5, 0 , 0 , 0 	  , 0.5 	, 0		  , 0 		 ]])
	#A inverse
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
	# need the pxmiss and pymiss, taken from the actual neutralino transverse momenta 
	# (this is cheating, of course)
	pxmiss = p4[0,1]+p8[0,1]
	pymiss = p4[0,2]+p8[0,2]

	C = np.matrix([ 2*minkowskidot(p1,p2) + 2*minkowskidot(p1,p3) + m1**2,
					2*minkowskidot(p2,p3) + m2**2,
					m3**2,
					pxmiss,
					2*minkowskidot(p5,p6) + 2*minkowskidot(p5,p7) + m5**2,
					2*minkowskidot(p6,p7) + m6**2,
					m7**2,
					pymiss])

	# Composite matrix & vector D and E
	D = Ainv*B
	E = Ainv*C.T

	# store D and E
	Dlist.append(D)
	Elist.append(E)


# Now to make a mass hypothesis (guess the correct one)
MZ = 545.421001 # Mass of ~uL
MY = 175.988171 # Mass of ~chi02
MX = 200.949431 # Mass of ~eL
MN = 98.361959 # Mass of ~chi01 (dark matter!)
MZprim = MZ
MYprim = MY
MXprim = MX
MNprim = MN
M = np.matrix([ MZ**2 , MY**2 , MX**2 , MN**2 , MZprim**2 , MYprim**2 , MXprim**2 , MNprim**2 ])

# Calculate the "chi-squared" error of the correct hypothesis
P = [] # store Pn
for n in range(N):
	Pn = Dlist[n]*M.T + Elist[n].T
	P.append(Pn) #store in case needed

	p4nsquared = Pn[3,0]**2 - Pn[0,0]**2 - Pn[1,0]**2 - Pn[2,0]**2
	p8nsquared = Pn[7,0]**2 - Pn[4,0]**2 - Pn[5,0]**2 - Pn[6,0]**2

	print p4nsquared - MN**2