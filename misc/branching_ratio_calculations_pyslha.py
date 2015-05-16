import pyslha
# sf = pyslha.readSLHAFile("slha_decay_table_CMSSM_realistic_point3.txt")
sf = pyslha.readSLHAFile("slha_decay_table_CMSSM_realistic_point4.txt")
# sf = pyslha.readSLHAFile("susyhit_softsusy_slha.out")


squarklist = [1000001,1000002,1000003,1000004]
sleptonlist = [2000011,2000013,1000011,1000013]

branching_ratio_gluino = [] # List of branching ratios for gluino to different squarks
branching_ratio_squarks = [] # List of branching ratio for each of the squarks to neutralino2
for squark in squarklist:
	for decaypair in sf.decays[squark].decays:
		if decaypair.ids[0] == 1000023:
			branching_ratio_squarks.append(decaypair.br)
	for gluinodecays in sf.decays[1000021].decays:
		if gluinodecays.ids[0] == squark:
			branching_ratio_gluino.append(gluinodecays.br)
print "branching_ratio_squarks =", branching_ratio_squarks
print "branching_ratio_gluino =", branching_ratio_gluino

branching_ratio_neutralino2 = []
for neutralinodecays in sf.decays[1000023].decays:
	if abs(neutralinodecays.ids[0]) in sleptonlist:
		branching_ratio_neutralino2.append(neutralinodecays)
print "branching_ratio_neutralino2 =", branching_ratio_neutralino2

branching_ratio_sleptons = [] # List of branching ratio for each of the sleptons to neutralino1
for slepton in sleptonlist:
	for decaypair in sf.decays[slepton].decays:
		if decaypair.ids[0] == 1000022:
			branching_ratio_sleptons.append(decaypair.br)
print "branching_ratio_sleptons =", branching_ratio_sleptons

crossections = [0.125E-02, 0.682E-02, 0.553E-03, 0.521E-02] # [squark-antisquark, squark-squark, gluino-gluino, squark-gluino]


# Put all together:
# Calculate total xsec to produce neutralino2, then multiply by br(neutralino2->rest of chain)^2
crossection_for_neutralino2 = 0
# loop over all squark pairings:
# SHOULD WE MULTIPLY BY 2 TO COUNT ANTISQUARKS? Yes, I think.
for squark1 in range(len(squarklist)):
	for squark2 in range(len(squarklist)):
		#squark-squark:
		crossection_for_neutralino2 += crossections[0] * branching_ratio_squarks[squark1] * branching_ratio_squarks[squark2]
		#squark-antisquark
		crossection_for_neutralino2 += crossections[1] * branching_ratio_squarks[squark1] * branching_ratio_squarks[squark2]
		#squark-gluino (assume WLOG that the first squark comes from the gluino:
		crossection_for_neutralino2 += crossections[2] * branching_ratio_gluino[squark1] * branching_ratio_squarks[squark1] * branching_ratio_squarks[squark2]
		#gluino-gluino
		crossection_for_neutralino2 += crossections[3] * branching_ratio_gluino[squark1] * branching_ratio_gluino[squark2] * branching_ratio_squarks[squark1] * branching_ratio_squarks[squark2]
# Multiply by two to count antisquarks also
crossection_for_neutralino2 *= 2
print "crossection_for_neutralino2 [pb] =", crossection_for_neutralino2
# Calculate br(neutralino2 -> rest of chain)
br_neutralino2_to_chain = 0
iterator = 0
for slepton in sleptonlist:
	for br_ntl2_entry in branching_ratio_neutralino2:
		if slepton == br_ntl2_entry.ids[0]:
			br_neutralino2_to_chain += 2* (br_ntl2_entry.br * branching_ratio_sleptons[iterator]) # Multiply by two to count antiparticles also

	iterator += 1
print br_neutralino2_to_chain


# Calculate total chain cross-section
chain_crosssection = crossection_for_neutralino2 * br_neutralino2_to_chain
print "chain cross-section =", chain_crosssection

# Calculate expected number of particles given integrated luminosity
integrated_luminosity = 300 * 10^3 # pb^-1
N = chain_crosssection * integrated_luminosity
print "Expected particle yield =", N