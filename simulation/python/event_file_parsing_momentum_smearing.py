import sys
import numpy as np

def minkowskidot(a,b):
	# Inner product in Minkowski space
	return float(a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2])
def minkowskinorm(a):
	# Inner product of a with a in Minkowski space
	return minkowskidot(a,a)

def smear(p, res):
	smearfactor = -1
	while smearfactor <= 0:
		smearfactor = np.random.normal(1, res)
	return p*smearfactor
def smearW(p, res):
	# Smear like Webber does to keep invariant mass unchanged
	# print minkowskinorm(p)
	smearfactor = -1
	while smearfactor <= 0:
		smearfactor = np.random.normal(1, res)
	minvsq = p[3]**2 - p[0]**2 - p[1]**2 - p[2]**2
	p[0:3] = p[0:3]*smearfactor
	p[3] = np.sqrt(p[0]**2+p[1]**2+p[2]**2+minvsq)

	# print minkowskinorm(p)

	return p

infile = open("../events/herwigpp-9563-events-complete-momcons-20150314_only_OFL.dat",'r')
outfile = open("../events/herwigpp-9563-events-complete-momcons-20150314_only_OFL-5percent_WEBBERmomentum_smearing.dat", 'w')

lines = infile.readlines()

for i in range(len(lines)/9):
	i *= 9
	title = lines[i]
	quark1 = lines[i+1].split()
	lepton11 = lines[i+2].split()
	lepton12 = lines[i+3].split()
	neutralino1 = lines[i+4].split()
	quark2 = lines[i+5].split()
	lepton21 = lines[i+6].split()
	lepton22 = lines[i+7].split()
	neutralino2 = lines[i+8].split()

	# Smear momenta according to a gaussian
	res = 0.05 # percentage/100 smearing
	pquark1 = smearW(np.array([float(quark1[1]), float(quark1[2]), float(quark1[3]), float(quark1[4])]) , res)
	plepton11 = smearW(np.array([float(lepton11[1]), float(lepton11[2]), float(lepton11[3]), float(lepton11[4])]) , res)
	plepton12 = smearW(np.array([float(lepton12[1]), float(lepton12[2]), float(lepton12[3]), float(lepton12[4])]) , res)
	pneutralino1 = smearW(np.array([float(neutralino1[1]), float(neutralino1[2]), float(neutralino1[3]), float(neutralino1[4])]) , res)
	pquark2 = smearW(np.array([float(quark2[1]), float(quark2[2]), float(quark2[3]), float(quark2[4])]) , res)
	plepton21 = smearW(np.array([float(lepton21[1]), float(lepton21[2]), float(lepton21[3]), float(lepton21[4])]) , res)
	plepton22 = smearW(np.array([float(lepton22[1]), float(lepton22[2]), float(lepton22[3]), float(lepton22[4])]) , res)
	pneutralino2 = smearW(np.array([float(neutralino2[1]), float(neutralino2[2]), float(neutralino2[3]), float(neutralino2[4])]) , res)

	# Write results to file
	outfile.write(title)
	outfile.write("%s %f %f %f %f NULL\n" %(quark1[0], pquark1[0], pquark1[1], pquark1[2], pquark1[3]))
	outfile.write("%s %f %f %f %f NULL\n" %(lepton11[0], plepton11[0], plepton11[1], plepton11[2], plepton11[3]))
	outfile.write("%s %f %f %f %f NULL\n" %(lepton12[0], plepton12[0], plepton12[1], plepton12[2], plepton12[3]))
	outfile.write("%s %f %f %f %f NULL\n" %(neutralino1[0], pneutralino1[0], pneutralino1[1], pneutralino1[2], pneutralino1[3]))
	outfile.write("%s %f %f %f %f NULL\n" %(quark2[0], pquark2[0], pquark2[1], pquark2[2], pquark2[3]))
	outfile.write("%s %f %f %f %f NULL\n" %(lepton21[0], plepton21[0], plepton21[1], plepton21[2], plepton21[3]))
	outfile.write("%s %f %f %f %f NULL\n" %(lepton22[0], plepton22[0], plepton22[1], plepton22[2], plepton22[3]))
	outfile.write("%s %f %f %f %f NULL\n" %(neutralino2[0], pneutralino2[0], pneutralino2[1], pneutralino2[2], pneutralino2[3]))

infile.close()
outfile.close()