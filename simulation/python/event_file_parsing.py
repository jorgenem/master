import sys
infile = open("../events/herwigpp-9563-events-complete-momcons-20150314.dat",'r')
outfile = open("../events/herwigpp-9563-events-complete-momcons-20150314_only_OFL.dat", 'w')

lines = infile.readlines()

for i in range(len(lines)/9):
	i *= 9
	title = lines[i]
	quark1 = lines[i+1]
	lepton11 = lines[i+2]
	lepton12 = lines[i+3]
	neutralino1 = lines[i+4]
	quark2 = lines[i+5]
	lepton21 = lines[i+6]
	lepton22 = lines[i+7]
	neutralino2 = lines[i+8]

	# throw away events where all leptons are equal
	if abs(int(lepton11.split()[0])) != abs(int(lepton21.split()[0])):
		outfile.write(title)
		outfile.write(quark1)
		outfile.write(lepton11)
		outfile.write(lepton12)
		outfile.write(neutralino1)
		outfile.write(quark2)
		outfile.write(lepton21) 
		outfile.write(lepton22)
		outfile.write(neutralino2)
